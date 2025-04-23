//===----------------------------------------------------------------------===//
//
//        Please cite the following resources in any derivative works:
//                 https://doi.org/10.5281/zenodo.10456866
//                 https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//
// 
// Torsion warping is due to [4], shear warping is due to [5]
//
// [1] Simo J.C. (1985): A finite strain beam formulation. The three-dimensional
//     dynamic problem. Part I.
//     Computer Methods in Applied Mechanics and Engineering, 49(1):55–70.
//     https://doi.org/10.1016/0045-7825(85)90050-7
//
// [2] Simo J.C., Vu-Quoc L. (1986): A three-dimensional finite-strain rod model
//     Part II: Computational aspects.
//     Computer Methods in Applied Mechanics and Engineering, 58(1):79–116.
//     https://doi.org/10/b8wd4z
//
// [3] Perez C.M., and Filippou F.C. (2024):
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg. 
//     https://doi.org/10.1002/nme.7506
//
// [4] Simo J.C., Vu-Quoc L. (1991): A Geometrically-exact rod model incorporating shear
//     and torsion-warping deformation. 
//     International Journal of Solids and Structures, 27(3):371-393.
//     doi: 10.1016/0020-7683(91)90089-x
//
// [5] Perez C.M., Filippou F.C., Mosalam K.M. (2025): Untitled work.
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
#include <cstddef>
#include <ExactFrame3d.h>
#include <Flag.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>

#include "for_int.tpp"
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Logging.h>
#include <Lagrange1D.cpp>
#include <quadrature/GaussLegendre1D.hpp>
#include <BeamIntegration.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
//

namespace OpenSees {

template<std::size_t nen, int nwm> static void
G_matrix(MatrixND<6+nwm,6+nwm> &G, 
         const VectorND<6+2*nwm>& s, const Vector3D& dx, 
         double shape[2][nen], 
         int i, int j)
{
  //
  // This is the sum of Equation (B4), and the unnumbered equation between (B3) and (B4).
  //
  auto sn = Hat(&s[0]);
  auto sm = Hat(&s[3]);
  G.assemble(         sn, 0, 3, -shape[1][i]*shape[0][j]);
  G.assemble(         sn, 3, 0,  shape[1][j]*shape[0][i]);

  G.assemble(         sm, 3, 3, -shape[1][i]*shape[0][j]);
  G.assemble( Hat(dx)*sn, 3, 3,  shape[0][i]*shape[0][j]);
}

template<std::size_t nen, int nwm> static void
B_nat(MatrixND<6+2*nwm,6+nwm> &B, double shape[2][nen], const Vector3D& dx, int n)
{
  //
  // NOTE This is the transpose of B in Equation (B3) from the paper by 
  // Perez and Filippou (2024)
  //
  B.zero();
  for (int i=0; i<6+nwm; i++)
    B(i,i) = shape[1][n];

  for (int i=0; i<nwm; i++)
    B(6+nwm+i, 6+i) = shape[0][n];
  
  //
  // B(1:3, 4:end) = shape*Hat(dx);
  //
  B(0,3) =  0;
  B(0,4) = -shape[0][n]*dx[2];
  B(0,5) =  shape[0][n]*dx[1];

  B(1,3) =  shape[0][n]*dx[2];
  B(1,4) =  0;
  B(1,5) = -shape[0][n]*dx[0];

  B(2,3) = -shape[0][n]*dx[1];
  B(2,4) =  shape[0][n]*dx[0];
  B(2,5) =  0;

}



template<std::size_t nen, int nwm>
ExactFrame3d<nen,nwm>::ExactFrame3d(int tag, 
                                    std::array<int,nen>& nodes, 
                                    FrameSection *section[nen-1],
                                    FrameTransform3d& transf)
    : FiniteElement<nen,ndm,ndf>(tag, 0, nodes),
      transform(&transf),
      logarithm(Logarithm::None),
      stencil(nullptr)
{
//  double wt[nip];
//  double xi[nip];
//  beamIntegr->getSectionLocations(numSections, L, xi);
//  beamIntegr->getSectionWeights(numSections, L, wt);

    p.zero();
    K.zero();

    for (int i=0; i<nip; i++) {
      pres[i].point  = 0.0;
      pres[i].weight = 0.0;
      pres[i].material=section[i]->getFrameCopy(scheme);
    }

}


template<std::size_t nen, int nwm>
ExactFrame3d<nen,nwm>::~ExactFrame3d()
{
  for (GaussPoint& point : pres)
    if (point.material != nullptr)
      delete point.material;

  if (stencil != nullptr)
    delete stencil;
}


template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::setNodes()
{
  auto& theNodes = this->FiniteElement<nen,3,6+nwm>::theNodes;

  if (transform->initialize(theNodes[0], theNodes[nen-1]) != 0) {
      opserr << " -- Error initializing coordinate transformation\n";
      return -1;
  }
  const Vector& xi = theNodes[    0]->getCrds();
  const Vector& xj = theNodes[nen-1]->getCrds();
  double L = (xi-xj).Norm();
  jxs = L;


  // Node locations in local (scalar) coordinate
  for (unsigned i=0; i < nen; i++)
    xn[i] = i*L/double(nen-1);

  GaussLegendre<1, nip>    gauss;
  for (unsigned i=0; i < nip; i++) {
    pres[i].point  = (gauss.pts[i] + 1.0)*L/2.0;
    pres[i].weight =  gauss.wts[i]*L/2.0;
    lagrange<nen>(pres[i].point, xn, pres[i].shape);
  }

  // Zero out the state of the Gauss pres
  this->revertToStart();

  return 0;
}

template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::revertToStart()
{
  // Revert the transformation to start
  if (transform->revertToStart() != 0)
    return -2;

  Vector E1(3), E2(3), E3(3);
  transform->getLocalAxes(E1, E2, E3);

  // for (int i=0; i<ndm; i++) {
  //   R0(0,i) =  E1[i];
  //   R0(1,i) =  E2[i];
  //   R0(2,i) =  E3[i];
  // }
  for (int i=0; i<ndm; i++) {
    R0(i,0) =  E1[i];
    R0(i,1) =  E2[i];
    R0(i,2) =  E3[i];
  }

  // Revert the of the Gauss pres to start
  for (GaussPoint& point : pres) {
    point.curvature.zero();
    point.rotation = R0;
    if (point.material->revertToStart() != 0)
      return -1;
  }
  past = pres;

  // Revert the element state to start
  // NOTE: This assumes that there are zero initial stresses?
  p.zero();
  K.zero();

  return 0;
}


template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::revertToLastCommit()
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
ExactFrame3d<nen,nwm>::update()
{
  const Vector3D D {1, 0, 0};
  auto& theNodes = this->FiniteElement<nen,3,ndf>::theNodes;

  //
  // Collect nodal parameters
  //
  VectorND<ndf> ddu[nen];
  for (unsigned i=0; i < nen; i++) {
    const Vector& ddui = theNodes[i]->getIncrDeltaDisp();
    for (int j=0; j<ndf; j++)
      ddu[i][j] = ddui[j];
  }

  // Form displaced node locations xyz
  VectorND<ndm> xyz[nen];
  double uwarp[nen][nwm]{};
  for (unsigned i=0; i < nen; i++) {
    const Vector& xi = theNodes[i]->getCrds();
    const Vector& ui = theNodes[i]->getTrialDisp();
    for (int j=0; j<ndm; j++)
      xyz[i][j] = xi[j] + ui[j];
    for (int j=0; j<nwm; j++)
      uwarp[i][j] = ui[6+j];
  }

  //
  // Gauss loop
  //
  p.zero();
  K.zero();
  for (int i=0; i<nip; i++) {
    //
    // Interpolate
    //
    Vector3D dx     {0.0};
    Vector3D theta  {0.0};
    Vector3D dtheta {0.0};

    for (unsigned j=0; j < nen; j++) {
      for (int l=0; l<3; l++)
        dx[l]     += pres[i].shape[1][j]*xyz[j][l];
      for (int l=0; l<3; l++)
        theta[l]  += pres[i].shape[0][j]*ddu[j][l+3];
      for (int l=0; l<3; l++)
        dtheta[l] += pres[i].shape[1][j]*ddu[j][l+3];
    }

    double warp[nwm]{};
    double dwarp[nwm]{};
    if (false) {//(nsr == 7 && scheme[6] == FrameStress::Bishear) {
      for (unsigned j=0; j < nen; j++)
        ;
        // dwarp += pres[i].shape[0][j]*uwarp[j];
    }
    else {//if (nsr == 8 && scheme[6] == FrameStress::Bishear) {
      for (int k=0; k<nwm; k++) {
        for (unsigned j=0; j < nen; j++) {
          warp[k]  += pres[i].shape[0][j]*uwarp[j][k];
          dwarp[k] += pres[i].shape[1][j]*uwarp[j][k];
        }
      }
    }

    //
    //
    Matrix3D dR = ExpSO3(theta);
    Matrix3D R = dR*pres[i].rotation;

    pres[i].rotation = R;

    Vector3D omega = dR*pres[i].curvature;
    // TODO: choose 'R/L'
//  pres[i].curvature = omega + TanSO3(theta, 'R')*dtheta;
    pres[i].curvature = omega + dExpSO3(theta)*dtheta;

    Vector3D gamma = (R^dx) - D;
    Vector3D kappa = R^pres[i].curvature;

    VectorND<nsr> e {
      gamma[0], gamma[1], gamma[2],
      kappa[0], kappa[1], kappa[2],
    };
    for (int j=0; j<nwm; j++) {
      e[6+j] = dwarp[j];
      e[6+nwm+j] = warp[j];
    }

    FrameSection& section = *pres[i].material;

    if (section.template setTrialState<nsr,scheme>(e) != 0)
      return -1;

    VectorND<nsr> s = section.getResultant<nsr,scheme>();
    MatrixND<nsr,nsr> Ks = section.getTangent<nsr,scheme>(State::Pres);

    //
    //
    // A = diag(R, R);
    // Note that this is transposed
    MatrixND<nsr,nsr> A {{
      {R(0,0), R(1,0), R(2,0), 0, 0, 0},
      {R(0,1), R(1,1), R(2,1), 0, 0, 0},
      {R(0,2), R(1,2), R(2,2), 0, 0, 0},
      {0, 0, 0, R(0,0), R(1,0), R(2,0)},
      {0, 0, 0, R(0,1), R(1,1), R(2,1)},
      {0, 0, 0, R(0,2), R(1,2), R(2,2)},
    }};
    for (int j=0; j<2*nwm; j++)
      A(6+j,6+j) = 1.0;

    MatrixND<nsr,ndf> B[nen];
    for (unsigned j=0; j<nen; j++) {
      MatrixND<nsr,ndf> Bj;
      Bj.zero();
      B_nat<nen,nwm>(Bj,  pres[i].shape, dx, j);
      B[j] = A^Bj;

      // p += B s w
      VectorND<ndf> pj = B[j]^s;
      for (int l=0; l<ndf; l++)
        p[j*ndf+l] += pres[i].weight * pj[l];
    }

    // Material Tangent
    MatrixND<ndf,ndf> Kjk;
    for (unsigned j=0; j<nen; j++) {
      for (unsigned k=0; k<nen; k++) {
        Kjk.addMatrixTripleProduct(0.0, B[j], Ks, B[k], pres[i].weight);

        for (int ii=0; ii<ndf; ii++) {
          for (int jj=0; jj<ndf; jj++) {
            K(j*ndf+ii, k*ndf+jj) += Kjk(ii,jj);
          }
        }
      }
    }

    // Geometric Tangent
    MatrixND<ndf,ndf> G;
    for (unsigned j=0; j<nen; j++) {
      for (unsigned k=0; k<nen; k++) {
        G.zero();
        G_matrix<nen,nwm>(G, A*s, dx, pres[i].shape, j, k);
        K.assemble(G, ndf*j, ndf*k, pres[i].weight);
      }
    }
  } // Main Gauss loop


  for (FrameLoad* load : frame_loads) {
    for (auto [xp, wp] : load->quadrature()) {
      double shp[2][nen];
      lagrange<nen>(xp*jxs, xn, shp);
      Versor q;
      if (xp == 0.0)
        q = theNodes[0]->getTrialRotation();
      else if (xp == 1.0)
        q = theNodes[nen-1]->getTrialRotation();
      else {
        opserr << "ERROR\n";
        Vector3D v{};
        for (unsigned i=0; i<nen; i++)
          v += shp[0][i]*theNodes[i]->getTrialDisp();
        q = VersorFromMatrix(R0*ExpSO3(v));
      }
      Matrix3D R  = MatrixFromVersor(q);
      const double w = wp; // /jxs;
      const double xc = xp;
      for_int<nen>([&](auto i_) constexpr {
        constexpr int i = i_.value;
        load->addLoadAtPoint<i,nen,ndf>(p, xc, w*shp[0][i], R0, R);
        for_int<nen>([&](auto j_) constexpr {
          constexpr int j = j_.value;
          load->addTangAtPoint<i,j,nen,ndf>(K, xc, w*shp[0][i]*shp[0][j], R0, R);
        });
      });
    }
  }
  return OpenSees::Flag::Success;
}

template<std::size_t nen, int nwm>
const Vector &
ExactFrame3d<nen,nwm>::getResistingForce()
{
  thread_local Vector wrapper;
  wrapper.setData(p);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame3d<nen,nwm>::getTangentStiff()
{
  thread_local Matrix wrapper;
  wrapper.setData(K);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame3d<nen,nwm>::getInitialStiff()
{
  static MatrixND<ndf*nen,ndf*nen> Ki{};
  static Matrix wrapper(Ki);
  return wrapper;
}

template<std::size_t nen, int nwm>
const Matrix &
ExactFrame3d<nen,nwm>::getMass()
{
  // TODO
  static MatrixND<ndf*nen,ndf*nen> M{};
  static Matrix wrapper(M);
  return wrapper;
}

template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::addLoad(ElementalLoad* theLoad, double loadFactor)
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
  // TODO: compute conservative load on flag == -1?
#if 0
  else {
    c_loads[frame_load->getTag()] = VectorND<ndf*nn>{0.0};
  }
#endif
    return -1;

  return 0;
}

template<std::size_t nen, int nwm>
Response*
ExactFrame3d<nen,nwm>::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  
  Response* theResponse = nullptr;
  double L = 0;
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
  output.attr("node2",  node_tags(1));

  //
  // compare argv[0] for known response types
  //

  // Global forces
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

  // Local force
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
    int numResponse                 = 0;
    const int numSections = pres.size();

    for (int i = 0; i < numSections; i++) {
      output.tag("GaussPointOutput");
      output.attr("number", i + 1);
      output.attr("eta", pres[i].point * L);

      Response* theSectionResponse = pres[i].material->setResponse(&argv[1], argc - 1, output);

      if (theSectionResponse != 0)
        numResponse = theCResponse->addResponse(theSectionResponse);
    }

    if (numResponse == 0) // no valid responses found
      delete theCResponse;
    else
      theResponse = theCResponse;
  }

  // 10-11: Integration
  else if (strcmp(argv[0], "integrationpres") == 0)
    theResponse = new ElementResponse(this, 10, Vector(pres.size()));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(pres.size()));

  // 110-111: sections
  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(pres.size()));

  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(pres.size(), 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(pres.size(), 3));
  }

  else if (strstr(argv[0], "section") != 0) {

    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= (int)pres.size() && argc > 2) {
        if (this->setState(State::Init) != 0)
          return nullptr;

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", 2.0 * pres[sectionNum - 1].point - 1.0);

        if (strcmp(argv[2], "dsdh") != 0) {
          theResponse = pres[sectionNum - 1].material->setResponse(&argv[2], argc - 2, output);
        } else {
          int order         = pres[sectionNum - 1].material->getOrder();
          theResponse       = new ElementResponse(this, 76, Vector(order));
          Information& info = theResponse->getInformation();
          info.theInt       = sectionNum;
        }

        output.endTag();

      } else if (sectionNum == 0) { 
        // argv[1] was not an int, we want all sections,

        CompositeResponse* theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        const int numSections = pres.size();

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", pres[i].point * L);

          Response* theSectionResponse = pres[i].material->setResponse(&argv[1], argc - 1, output);

          if (theSectionResponse != 0) {
            numResponse = theCResponse->addResponse(theSectionResponse);
          }
        }

        if (numResponse == 0) // no valid responses found
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }

  else if (strcmp(argv[0], "energy") == 0) {
    return new ElementResponse(this, 2000, 0.0);
  }

  output.endTag();

  return theResponse;
}

template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::getResponse(int responseID, Information &info)
{

  double L = 0;
  if (this->setState(State::Init) == 0) {
    auto& theNodes = this->FiniteElement<nen,3,6+nwm>::theNodes;
    const Vector& xi = theNodes[    0]->getCrds();
    const Vector& xj = theNodes[nen-1]->getCrds();
    L = (xi-xj).Norm();
  }


  // NOTE: This is will never be called with Respond::GlobalForce;
  // it gets intercepted by Domain::getElementResponse
  if (responseID == Respond::GlobalForce)
    return info.setVector(this->getResistingForce());

  else if (responseID == Respond::LocalForce) {
    thread_local VectorND<nen*ndf> q{0.0};
    Vector q_resp(q);

    auto p = this->getResistingForce();
  
    for (unsigned i=0; i<nen; i++) {
      for (int j=0; j<3; j++) {
        q(i*ndf+j) = 0;
        q(i*ndf+j+3) = 0;
        for (int k=0; k<3; k++) {
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
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    Vector locs(pres.size());
    for (int i = 0; i < nip; i++)
      locs[i] = pres[i].point * L;

    return info.setVector(locs);
  }

  else if (responseID == 11) {
    // ensure we have L, xi[] and wt[]
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
ExactFrame3d<nen,nwm>::sendSelf(int commitTag, Channel& theChannel)
{
  // TODO
  return -1;
}

template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  // TODO
  return -1;
}

template<std::size_t nen, int nwm>
void
ExactFrame3d<nen,nwm>::Print(OPS_Stream& stream, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {

    stream << OPS_PRINT_JSON_ELEM_INDENT << "{";
    stream << "\"name\": " << this->getTag() << ", ";
    stream << "\"type\": \"" << this->getClassType() << "\", ";
    stream << "\"nodes\": [";
    for (int i = 0; i < (int)nen - 1; i++)
      stream << node_tags(i) << ", ";
    stream << node_tags(nen - 1) << "]";
    stream << ", ";


    stream << "\"sections\": [";
    for (decltype(pres.size()) i = 0; i < pres.size() - 1; i++)
      stream << pres[i].material->getTag() << ", ";
    stream << pres[pres.size() - 1].material->getTag() << "]";
    stream << ", ";

    stream << "\"transform\": " << transform->getTag()  ;
    stream << "}";
  }
}

template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  int result = -1;


  // Section response
  if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      const int numSections = pres.size();
      double xi[nip];
      double L = jxs;
      stencil->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(pres[i].point - sectionLoc) < minDistance) {
          minDistance = fabs(pres[i].point - sectionLoc);
          sectionNum  = i;
        }
      }

      return pres[sectionNum].material->setParameter(&argv[2], argc - 2, param);
    }
  }

  // If the parameter belongs to a particular section or lower
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= (int)pres.size())
      return pres[sectionNum - 1].material->setParameter(&argv[2], argc - 2, param);

    else
      return -1;
  }

  // If the parameter belongs to all sections or lower
  if (strstr(argv[0], "allSections") != 0) {

    if (argc < 2)
      return -1;

    for (GaussPoint& point : pres) {
      int ok = point.material->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }

    return result;
  }

  if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return stencil->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to everything

  for (GaussPoint& point : pres) {
    int ok = point.material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  int ok = stencil->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}


template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::updateParameter(int parameterID, Information& info)
{
  return -1;
}


template<std::size_t nen, int nwm>
int
ExactFrame3d<nen,nwm>::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}

} // namespace OpenSees
