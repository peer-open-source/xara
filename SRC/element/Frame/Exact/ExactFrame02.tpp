
//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// ExactFrame02 implementation.
// See ExactFrame02.h for a summary.
//
//===----------------------------------------------------------------------===//
#include <cstddef>
#include "ExactFrame02.h"
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>

#include <utility/Unroll.h>
#include <CrdTransf.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Logging.h>
#include <Lagrange1D.cpp>
#include <quadrature/GaussLegendre1D.hpp>
#include <BeamIntegration.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <FrameLoad.h>

#include "ExactInterpolation.h"

namespace OpenSees {


template<std::size_t nen, int nwm>
ExactFrame02<nen, nwm>::ExactFrame02(int tag,
                                     std::array<int, nen>& nodes,
                                     FrameSection* section[nen - 1],
                                     CrdTransf& transf,
                                     Rotations::Parameters param)
  : FiniteElement<nen, ndm, ndf>(tag, 0, nodes, 1),
    xn{0},
    jxs(0),
    R0(),
    reference_rotation{},
    reference_curvature{},
    transform(&transf),
    parameterization(param),
    stencil(nullptr),
    parameterID(0)
{
  p.zero();
  K.zero();

  for (int i = 0; i < nip; i++) {
    pres[i].point    = 0.0;
    pres[i].weight   = 0.0;
    pres[i].material = section[i]->getFrameCopy();
  }
}


template <std::size_t nen, int nwm>
ExactFrame02<nen,nwm>::~ExactFrame02()
{
  for (GaussPoint& point : pres)
    if (point.material != nullptr)
      delete point.material;

  if (stencil != nullptr)
    delete stencil;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::setNodes()
{
  auto& theNodes = this->FiniteElement<nen,3,6+nwm>::theNodes;

  if (transform->initialize(theNodes[0], theNodes[nen-1]) != 0) {
    opserr << " -- Error initializing coordinate transformation\n";
    return -1;
  }

  const Vector& xi = theNodes[    0]->getCrds();
  const Vector& xj = theNodes[nen-1]->getCrds();
  const double L = (xj-xi).Norm();
  jxs = L;

  for (unsigned i = 0; i < nen; i++)
    xn[i] = i*L/double(nen-1);

  GaussLegendre<1, nip> gauss;
  for (unsigned i = 0; i < nip; i++) {
    pres[i].point  = (gauss.pts[i] + 1.0)*L/2.0;
    pres[i].weight =  gauss.wts[i]*L/2.0;
    lagrange<nen>(pres[i].point, xn, pres[i].shape);
  }

  this->revertToStart();
  return 0;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::revertToStart()
{
  if (transform->revertToStart() != 0)
    return -2;

  Vector E1(3), E2(3), E3(3);
  transform->getLocalAxes(E1, E2, E3);

  for (int i = 0; i < ndm; i++) {
    R0(i,0) = E1[i];
    R0(i,1) = E2[i];
    R0(i,2) = E3[i];
  }

  for (int i = 0; i < nip; ++i) {
    reference_rotation[i]  = R0;
    reference_curvature[i].zero();
  }

  for (GaussPoint& point : pres) {
    point.curvature.zero();
    point.rotation = R0;
    if (point.material->revertToStart() != 0)
      return -1;
  }
  past = pres;

  p.zero();
  K.zero();

  return 0;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::revertToLastCommit()
{
  pres = past;

  for (GaussPoint& point : pres) {
    FrameSection& section = *point.material;
    if (section.revertToLastCommit() != 0)
      return -1;
  }

  return 0;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::update()
{
  constexpr static Vector3D D{1.0, 0.0, 0.0};
  auto& theNodes = this->FiniteElement<nen,3,ndf>::theNodes;

  //
  // Collect nodal parameters
  //
  VectorND<ndf> uparam[nen];
  VectorND<ndm> xyz[nen];
  std::array<std::array<double,nwm>,nen> uwarp{};

  for (unsigned i = 0; i < nen; i++) {
    const Vector& xi = theNodes[i]->getCrds();
    const Vector& ui = theNodes[i]->getTrialDisp();

    const Vector* up = &ui;
    if (parameterization == Rotations::Parameters::Incr)
      up = &theNodes[i]->getIncrDisp();

    for (int j = 0; j < ndm; j++)
      xyz[i][j] = xi[j] + ui[j];

    for (int j = 0; j < ndf; ++j)
      uparam[i][j] = (*up)[j];

    for (int j = 0; j < nwm; j++)
      uwarp[i][j] = ui[6+j];
  }

  p.zero();
  K.zero();

  for (int i = 0; i < nip; i++) {
    Vector3D dx    {0.0};
    Vector3D theta {0.0};
    Vector3D dtheta{0.0};

    for (unsigned j = 0; j < nen; j++) {
      for (int l = 0; l < 3; l++)
        dx[l] += pres[i].shape[1][j]*xyz[j][l];

      for (int l = 0; l < 3; l++) {
        theta[l]  += pres[i].shape[0][j]*uparam[j][l+3];
        dtheta[l] += pres[i].shape[1][j]*uparam[j][l+3];
      }
    }

    std::array<double,nwm> warp{};
    std::array<double,nwm> dwarp{};
    for (int k = 0; k < nwm; k++) {
      for (unsigned j = 0; j < nen; j++) {
        warp[k]  += pres[i].shape[0][j]*uwarp[j][k];
        dwarp[k] += pres[i].shape[1][j]*uwarp[j][k];
      }
    }

    const Matrix3D& R0 =
      (parameterization == Rotations::Parameters::Init) ?
        reference_rotation[i] : past[i].rotation;

    const Vector3D& k0 =
      (parameterization == Rotations::Parameters::Init) ?
        reference_curvature[i] : past[i].curvature;

    const Matrix3D dR = ExpSO3(theta);
    const Matrix3D R  = dR*R0;

    pres[i].rotation = R;
    pres[i].curvature = dR*k0 + TExpSO3(theta)*dtheta;

    const Vector3D gamma = (R^dx) - D;
    const Vector3D kappa = R^pres[i].curvature;

    VectorND<nsr> e {
      gamma[0], gamma[1], gamma[2],
      kappa[0], kappa[1], kappa[2],
    };

    for (int j = 0; j < nwm; j++) {
      e[6+j]     = dwarp[j];
      e[6+nwm+j] = warp[j];
    }

    FrameSection& section = *pres[i].material;
    if (section.template setTrialState<nsr,scheme>(e) != 0)
      return -1;

    const VectorND<nsr> S = section.getResultant<nsr,scheme>();
    const MatrixND<nsr,nsr> Ks = section.getTangent<nsr,scheme>(State::Pres);

    MatrixND<nsr,nsr> A{};
    RotationOperator(A, R);

    const VectorND<nsr> s = ApplyOperator(A, S);
    const MatrixND<nsr,nsr> ks = PushTangent(A, Ks);

    MatrixND<nsr,ndf> B[nen];
    for (unsigned j = 0; j < nen; j++) {
      B_log<nen,nwm>(B[j], pres[i].shape, dx, theta, dtheta, j);

      const VectorND<ndf> pj = B[j]^s;
      for (int l = 0; l < ndf; l++)
        p[j*ndf+l] += pres[i].weight*pj[l];
    }

    for (unsigned j = 0; j < nen; j++) {
      for (unsigned k = 0; k < nen; k++) {
        MatrixND<ndf,ndf> Kjk{};
        AddMatrixTriple(Kjk, B[j], ks, B[k], pres[i].weight);
        K.assemble(Kjk, ndf*j, ndf*k, 1.0);
      }
    }

    for (unsigned j = 0; j < nen; j++) {
      for (unsigned k = 0; k < nen; k++) {
        MatrixND<ndf,ndf> G{};
        G_log<nen,nwm>(G, s, S, dx, pres[i].shape, j, k, theta, dtheta);
        K.assemble(G, ndf*j, ndf*k, pres[i].weight);
      }
    }
  }

  //
  //
  //
  for (FrameLoad* load : frame_loads) {
    for (auto [xp, wp] : load->quadrature()) {
      const double w  = wp*jxs;
      const double xc = xp;

      double shp[2][nen];
      lagrange<nen>(xp*jxs, xn, shp);

      Versor q;
      if (xp == 0.0)
        q = theNodes[0]->getTrialRotation();
      else if (xp == 1.0)
        q = theNodes[nen-1]->getTrialRotation();
      else
        q = theNodes[0]->getTrialRotation().slerp(
              theNodes[nen-1]->getTrialRotation(), xp);

      const Matrix3D R  = MatrixFromVersor(q);

#ifndef _MSC_VER
      Unroll<0,nen>([&](auto i_) constexpr {
        constexpr int ii = i_.value;
        load->addLoadAtPoint<ii,nen,ndf>(p, xc, w*shp[0][ii], jxs, R0, R);
        Unroll<0,nen>([&](auto j_) constexpr {
          constexpr int jj = j_.value;
          load->addTangAtPoint<ii,jj,nen,ndf>(
            K, xc, w*shp[0][ii]*shp[0][jj], jxs, R0, R);
        });
      });
#endif
    }
  }

  return 0;
}

template<std::size_t nen, int nwm>
const Vector &
ExactFrame02<nen,nwm>::getResistingForce()
{
  thread_local Vector wrapper(0);
  wrapper.setData(p);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame02<nen,nwm>::getTangentStiff()
{
  thread_local Matrix wrapper(0,0);
  wrapper.setData(K);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame02<nen,nwm>::getInitialStiff()
{
  static MatrixND<ndf*nen,ndf*nen> Ki{};
  static Matrix wrapper(Ki);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame02<nen,nwm>::getMass()
{
  static MatrixND<ndf*nen,ndf*nen> M{};
  static Matrix wrapper(M);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Vector &
ExactFrame02<nen,nwm>::getResistingForceSensitivity(int grad)
{
  static VectorND<ndf*nen> dp;
  static Vector wrapper(dp);

  dp.zero();

  auto& theNodes = this->FiniteElement<nen,3,ndf>::theNodes;

  VectorND<ndf> uparam[nen];
  VectorND<ndm> xyz[nen];
  std::array<std::array<double,nwm>,nen> uwarp{};

  for (unsigned i = 0; i < nen; i++) {
    const Vector& xi = theNodes[i]->getCrds();
    const Vector& ui = theNodes[i]->getTrialDisp();

    const Vector* up = &ui;
    if (parameterization == Rotations::Parameters::Incr)
      up = &theNodes[i]->getIncrDisp();

    for (int j = 0; j < ndm; j++)
      xyz[i][j] = xi[j] + ui[j];

    for (int j = 0; j < ndf; ++j)
      uparam[i][j] = (*up)[j];

    for (int j = 0; j < nwm; j++)
      uwarp[i][j] = ui[6+j];
  }

  for (int i = 0; i < nip; i++) {
    Vector3D dx    {0.0};
    Vector3D theta {0.0};
    Vector3D dtheta{0.0};

    for (unsigned j = 0; j < nen; j++) {
      for (int l = 0; l < 3; l++)
        dx[l] += pres[i].shape[1][j]*xyz[j][l];

      for (int l = 0; l < 3; l++) {
        theta[l]  += pres[i].shape[0][j]*uparam[j][l+3];
        dtheta[l] += pres[i].shape[1][j]*uparam[j][l+3];
      }
    }

    const Matrix3D& R0 =
      (parameterization == Rotations::Parameters::Init) ?
        reference_rotation[i] : past[i].rotation;

    const Matrix3D R = ExpSO3(theta)*R0;

    FrameSection& section = *pres[i].material;
    const VectorND<nsr> dS = section.getResultantGradient<nsr,scheme>(grad, true);

    MatrixND<nsr,nsr> A{};
    RotationOperator(A, R);
    const VectorND<nsr> ds = ApplyOperator(A, dS);

    MatrixND<nsr,ndf> B[nen];
    for (unsigned j = 0; j < nen; j++) {
      B_log<nen,nwm>(B[j], pres[i].shape, dx, theta, dtheta, j);

      const VectorND<ndf> dpi = B[j]^ds;
      for (int l = 0; l < ndf; l++)
        dp[j*ndf+l] += pres[i].weight*dpi[l];
    }
  }

  return wrapper;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  int type = theLoad->getClassTag();

  if (type == LOAD_TAG_FrameLoad && loadFactor == 0.0)
    frame_loads.erase((FrameLoad*)theLoad);

  else if (type == LOAD_TAG_FrameLoad && loadFactor == 1.0) {
    FrameLoad* frame_load = (FrameLoad*)theLoad;
    if (!frame_load->conservative())
      frame_loads.insert(frame_load);
  }
  else
    return -1;

  return 0;
}

template<std::size_t nen, int nwm>
Response*
ExactFrame02<nen,nwm>::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = nullptr;
  double L = 0.0;
  if (this->setState(State::Init) == 0) {
    auto& theNodes = this->FiniteElement<nen,3,6+nwm>::theNodes;
    const Vector& xi = theNodes[    0]->getCrds();
    const Vector& xj = theNodes[nen-1]->getCrds();
    L = (xi-xj).Norm();
  }

  const ID& node_tags = this->getExternalNodes();
  output.tag("ElementOutput");
  output.attr("eleType", this->getClassType());
  output.attr("eleTag", this->getTag());
  output.attr("node1",  node_tags(0));
  output.attr("node2",  node_tags(nen-1));

  if (strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 ||
      strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, Respond::GlobalForce, Vector(nen*ndf));

  } else if (strcmp(argv[0], "localForce") == 0 ||
             strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, Respond::LocalForce, Vector(nen*ndf));

  } else if (strcmp(argv[0], "RayleighForces") == 0 ||
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, Vector(12));

  } else if (strcmp(argv[0], "sections") == 0) {
    if (this->setState(State::Init) != 0)
      return nullptr;

    CompositeResponse* theCResponse = new CompositeResponse();
    int numResponse = 0;
    const int numSections = pres.size();

    for (int i = 0; i < numSections; i++) {
      output.tag("GaussPointOutput");
      output.attr("number", i + 1);
      output.attr("eta", pres[i].point);

      Response* theSectionResponse =
        pres[i].material->setResponse(&argv[1], argc - 1, output);

      if (theSectionResponse != nullptr)
        numResponse = theCResponse->addResponse(theSectionResponse);
    }

    if (numResponse == 0)
      delete theCResponse;
    else
      theResponse = theCResponse;
  }

  else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(pres.size()));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(pres.size()));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(pres.size()));

  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(pres.size(), 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(pres.size(), 3));
  }

  else if (strstr(argv[0], "section") != nullptr) {
    if (argc > 1) {
      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= (int)pres.size() && argc > 2) {
        if (this->setState(State::Init) != 0)
          return nullptr;

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", 2.0 * pres[sectionNum - 1].point - 1.0);

        if (strcmp(argv[2], "dsdh") != 0) {
          theResponse =
            pres[sectionNum - 1].material->setResponse(&argv[2], argc - 2, output);
        } else {
          int order = pres[sectionNum - 1].material->getOrder();
          theResponse = new ElementResponse(this, 76, Vector(order));
          Information& info = theResponse->getInformation();
          info.theInt = sectionNum;
        }

        output.endTag();
      }

      else if (sectionNum == 0) {
        CompositeResponse* theCResponse = new CompositeResponse();
        int numResponse = 0;
        const int numSections = pres.size();

        for (int i = 0; i < numSections; i++) {
          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", pres[i].point * L);

          Response* theSectionResponse =
            pres[i].material->setResponse(&argv[1], argc - 1, output);

          if (theSectionResponse != nullptr)
            numResponse = theCResponse->addResponse(theSectionResponse);
        }

        if (numResponse == 0)
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }

  else if (strcmp(argv[0], "energy") == 0)
    return new ElementResponse(this, 2000, 0.0);

  output.endTag();
  return theResponse;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::getResponse(int responseID, Information &info)
{
  double L = 0.0;
  if (this->setState(State::Init) == 0) {
    auto& theNodes = this->FiniteElement<nen,3,6+nwm>::theNodes;
    const Vector& xi = theNodes[    0]->getCrds();
    const Vector& xj = theNodes[nen-1]->getCrds();
    L = (xi-xj).Norm();
  }

  if (responseID == Respond::GlobalForce)
    return info.setVector(this->getResistingForce());

  else if (responseID == Respond::LocalForce) {
    thread_local VectorND<nen*ndf> q{0.0};
    Vector q_resp(q);

    auto p = this->getResistingForce();

    for (unsigned i = 0; i < nen; i++) {
      for (int j = 0; j < 3; j++) {
        q(i*ndf+j)   = 0.0;
        q(i*ndf+j+3) = 0.0;
        for (int k = 0; k < 3; k++) {
          q(i*ndf+j)   += R0(k,j)*p(i*6+k);
          q(i*ndf+j+3) += R0(k,j)*p(i*6+k+3);
        }
      }
    }

    return info.setVector(q_resp);
  }

  else if (responseID == 19)
    return info.setMatrix(K);

  else if (responseID == 10) {
    if (this->setState(State::Init) != 0)
      return -1;

    Vector locs(pres.size());
    for (int i = 0; i < nip; i++)
      locs[i] = pres[i].point;

    return info.setVector(locs);
  }

  else if (responseID == 11) {
    if (this->setState(State::Init) != 0)
      return -1;

    Vector weights(pres.size());
    for (unsigned i = 0; i < nip; i++)
      weights(i) = pres[i].weight * L;

    return info.setVector(weights);
  }

  else if (responseID == 110) {
    const int numSections = pres.size();
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = pres[i].material->getTag();
    return info.setID(tags);
  }

  else if (responseID == 12)
    return info.setVector(this->getRayleighDampingForces());

  return -1;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::sendSelf(int commitTag, Channel& theChannel)
{
  return -1;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::recvSelf(int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
  return -1;
}

template<std::size_t nen, int nwm>
void
ExactFrame02<nen,nwm>::Print(OPS_Stream& stream, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    stream << OPS_PRINT_JSON_ELEM_INDENT << "{";
    stream << "\"name\": " << this->getTag() << ", ";
    stream << "\"type\": \"" << this->getClassType() << "\", ";
    stream << "\"nodes\": [";
    for (int i = 0; i < (int)nen - 1; i++)
      stream << node_tags(i) << ", ";
    stream << node_tags(nen - 1) << "], ";

    stream << "\"sections\": [";
    for (decltype(pres.size()) i = 0; i < pres.size() - 1; i++)
      stream << pres[i].material->getTag() << ", ";
    stream << pres[pres.size() - 1].material->getTag() << "], ";

    stream << "\"transform\": " << transform->getTag() << ", ";
    stream << "\"parameterization\": \""
           << (parameterization == Rotations::Parameters::Init ? "Init" : "Incr")
           << "\"";
    stream << "}";
  }
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strstr(argv[0], "sectionX") != nullptr) {
    if (argc > 2 && stencil != nullptr) {
      double sectionLoc = atof(argv[1]);

      const int numSections = pres.size();
      double xi[nip];
      double L = jxs;
      stencil->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      double minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(pres[i].point - sectionLoc) < minDistance) {
          minDistance = fabs(pres[i].point - sectionLoc);
          sectionNum  = i;
        }
      }

      return pres[sectionNum].material->setParameter(&argv[2], argc - 2, param);
    }
  }

  if (strstr(argv[0], "section") != nullptr) {
    if (argc < 3)
      return -1;

    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= (int)pres.size())
      return pres[sectionNum - 1].material->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  if (strstr(argv[0], "allSections") != nullptr) {
    if (argc < 2)
      return -1;

    for (GaussPoint& point : pres) {
      int ok = point.material->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }
    return result;
  }

  if (strstr(argv[0], "integration") != nullptr) {
    if (argc < 2 || stencil == nullptr)
      return -1;

    return stencil->setParameter(&argv[1], argc - 1, param);
  }

  for (GaussPoint& point : pres) {
    int ok = point.material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::updateParameter(int parameterID, Information& info)
{
  return -1;
}

template<std::size_t nen, int nwm>
int
ExactFrame02<nen,nwm>::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}

} // namespace OpenSees
