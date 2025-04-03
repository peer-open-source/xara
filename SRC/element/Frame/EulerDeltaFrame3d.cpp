//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Higher-order Euler frame formulation with C1 displacement interpolation
//
// Adapted from DispBeamColumnAsym? or NL
//
#include <array>
#include <math.h>
#include <string.h>

#include <EulerDeltaFrame3d.h>

#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

#include <Node.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Domain.h>
#include <Channel.h>
#include <Parameter.h>
#include <Information.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>

#define ELE_TAG_EulerDeltaFrame3d 0 // TODO


static int
getStrainMatrix(double xi, double L, const Vector& v, MatrixND<8,12>&B, MatrixND<4,8>&A)
{
    double jsx = 1/L;

    double dshp_v1  = 1.0 + 3.0 * xi * xi - 4.0 * xi;
    double ddshp_v1 = 6.0 * xi * jsx - 4.0 * jsx;
    double dshp_v2  = 3.0 * xi * xi - 2.0 * xi;
    double ddshp_v2 = 6.0 * xi * jsx - 2.0 * jsx;
    double dshp_w1  = -dshp_v1;
    double ddshp_w1 = -ddshp_v1;
    double dshp_w2  = -dshp_v2;
    double ddshp_w2 = -ddshp_v2;
    double shp_alpha1   = xi;

    double dv  =  dshp_v1*v[1] +  dshp_v2*v[2]; // v'
    double ddv = ddshp_v1*v[1] + ddshp_v2*v[2]; // v"
    double dw  =  dshp_w1*v[3] +  dshp_w2*v[4]; // w'
    double ddw = ddshp_w1*v[3] + ddshp_w2*v[4]; // w"
    double alpha   = shp_alpha1*v[5];           // phi
    double dalpha  = jsx*v[5];                  // phi'

    A(0, 0)  = 1.0;
    A(0, 1)  = (4.0 * v[1] - v[2]) / 30.0;
    A(0, 2)  = (4.0 * v[3] - v[4]) / 30.0;
    A(0, 3)  = (4.0 * v[2] - v[1]) / 30.0;
    A(0, 4)  = (4.0 * v[4] - v[3]) / 30.0;
    A(1, 7)  = 1.0;

    B( 0, 0)  = jsx;
    B( 1, 1)  = 1.0;
    B( 2, 3)  = 1.0;
    B( 3, 2)  = 1.0;
    B( 4, 4)  = 1.0;
    B( 5, 1)  = dshp_v1;
    B( 5, 2)  = dshp_v2;
    B( 6, 3)  = dshp_w1;
    B( 6, 4)  = dshp_w2;
    B( 7, 1)  = ddshp_v1;
    B( 7, 2)  = ddshp_v2;
    B( 8, 3)  = ddshp_w1;
    B( 8, 4)  = ddshp_w2;
    B( 9, 5)  = shp_alpha1;
    B(10, 5)  = jsx;

    return 0;
}

EulerDeltaFrame3d::EulerDeltaFrame3d(int tag, std::array<int,2>& nodes,
                                     std::vector<FrameSection*> &secs,
                                     BeamIntegration &bi,
                                     FrameTransform3d &coordTransf,
                                     double r, int cm, bool use_mass_)

    : FiniteElement(tag, ELE_TAG_EulerDeltaFrame3d, nodes),
      numSections(secs.size()),  sections(nullptr),
      theCoordTransf(nullptr),   beamInt(nullptr), 
      density(r), mass_flag(cm), use_density(use_mass_),
      parameterID(0)

{
  // Allocate arrays of pointers to FrameSections
  sections = new FrameSection *[secs.size()];

  for (int i = 0; i < secs.size(); i++) {
    // Get copies of the material model for each integration point
    sections[i] = secs[i]->getFrameCopy(scheme);
  }

  beamInt = bi.getCopy();

  theCoordTransf = coordTransf.getCopy();

  q0.zero();
}

EulerDeltaFrame3d::EulerDeltaFrame3d()
    : FiniteElement(0, ELE_TAG_EulerDeltaFrame3d),
      numSections(0), sections(nullptr),
      beamInt(nullptr),
      density(0.0), mass_flag(0), parameterID(0)
{
  q0.zero();
}

EulerDeltaFrame3d::~EulerDeltaFrame3d()
{
  for (int i = 0; i < numSections; i++) {
    if (sections[i])
      delete sections[i];
  }

  // Delete the array of pointers to FrameSection pointer arrays
  if (sections)
    delete[] sections;

  if (theCoordTransf)
    delete theCoordTransf;

  if (beamInt != nullptr)
    delete beamInt;
}


int
EulerDeltaFrame3d::commitState()
{
  int status = this->Element::commitState();

  // Loop over the integration points and commit the material states
  for (int i = 0; i < numSections; i++)
    status += sections[i]->commitState();

  // NOTE: This was not done in original element
  status += theCoordTransf->commitState();

  return status;
}

int
EulerDeltaFrame3d::revertToLastCommit()
{
  int status = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < numSections; i++)
    status += sections[i]->revertToLastCommit();

  // NOTE: This was not done in original element
  status += theCoordTransf->revertToLastCommit();

  return status;
}

int
EulerDeltaFrame3d::revertToStart()
{
  int status = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < numSections; i++)
    status += sections[i]->revertToStart();

  // NOTE: This was not done in original element
  status += theCoordTransf->revertToStart();

  return status;
}

int
EulerDeltaFrame3d::update()
{
  int err = 0;

  // Update the transformation
  theCoordTransf->update();

  // Get basic deformations
  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L        = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    double xi1   = xi[i];
    double dNv[] = {
                   1.0 + 3.0*xi1*xi1 - 4.0*xi1,
                   3.0*xi1*xi1 - 2.0*xi1 
    };
    double ddNv[] = {
                   6.0*xi1*jsx - 4.0*jsx,
                   6.0*xi1*jsx - 2.0*jsx 
    };
    double dNw1  = -dNv[0];
    double dNw2  = -dNv[1];
    double ddNw1 = -ddNv[0];
    double ddNw2 = -ddNv[1];
    double Nf1   =  xi1;

    double dx   = jsx*v[0];                     // u'
    double dy   =  dNw1*v(3) +  dNw2*v(4);      // y'
    double dz   =  dNv[0]*v(1) +  dNv[1]*v(2);  // z'
    double phi  = Nf1*v(5);                 // phi
                                            //
    double ddy  = ddNw1*v(3) + ddNw2*v(4);  // y"
    double ddz  = ddNv[0]*v(1) + ddNv[1]*v(2);  // z"
    double dphi = jsx*v[5];            // phi'

    VectorND<4> e {
            dx,
            ddz + ddy*phi,
           -ddy + ddz*phi,
            dphi
    };

    // Update the section state
    err += sections[i]->setTrialState<nsr,scheme>(e);
  }

  if (err != 0) {
    opserr << "EulerDeltaFrame3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

int
EulerDeltaFrame3d::setNodes()
{
  double L  = theCoordTransf->getInitialLength();
  beamInt->getSectionLocations(numSections, L, xi);
  beamInt->getSectionWeights(numSections, L, wt);
  return 0;
}

int
EulerDeltaFrame3d::getIntegral(Field field, State state, double& total)
{

  if (this->setState(State::Init) != 0)
    return -1;

  total = 0.0;
  switch (field) {

    // Integrate density to compute total mass
    case Field::Density: {
      double value = 0.0;
      for (int i=0; i< numSections; i++) {
        // First try using section's internal density
        if (sections[i]->getIntegral(Field::Density, state, value) == 0) {
          total += wt[i]*value;
        }
        // if that didnt work, just multiply by our density
        else if (sections[i]->getIntegral(Field::Unit, state, value) == 0) {
          total += wt[i]*density;
        }
        else {
          ; // TODO: This should be written to a log
        }
      }
      return 0;
    }

    case Field::PolarInertia: {
      for (int i=0; i< numSections; i++) {
        double A;
        if (sections[i]->getIntegral(Field::UnitYY, state, A) != 0)
          continue;

        // Get \int \rho y^2
        double Iz;
        if (sections[i]->getIntegral(Field::DensityYY, state, Iz) != 0) {
          // Section does not allow integrating density; try
          // integrating product of inertia and multiplying by rho
          if (sections[i]->getIntegral(Field::UnitYY, state, Iz) == 0)
            Iz *= density/A;
          else
            continue;
        }
        // Get \int \rho z^2
        double Iy;
        if (sections[i]->getIntegral(Field::DensityZZ, state, Iy) != 0) {
          if (sections[i]->getIntegral(Field::UnitZZ, state, Iy) == 0)
            Iy *= density/A;
          else
            continue;
        }
        total += wt[i]*(Iy + Iz);
      }
      return 0;
    }

    default:
      return -1;

  }
}

const Vector &
EulerDeltaFrame3d::getResistingForce()
{
  static Vector wrap(p);
  return wrap;
}

const Matrix &
EulerDeltaFrame3d::getTangentStiff()
{
  static MatrixND<ndf*nen, ndf*nen> kb;
  static MatrixND<nsr,  8> A;
  static MatrixND<  8, 12> B;
  static MatrixND<  8,  8> ks;
  static MatrixND< 12, 12> Gm;

  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  // Zero for integral
  kb.zero();
  p.zero();

  A.zero();
  B.zero();
  Gm.zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    getStrainMatrix(xi[i], L, v, B, A);

    // Get the section tangent stiffness and stress resultant
    MatrixND<nsr,nsr> Ks = sections[i]->getTangent<nsr,scheme>(State::Pres);
    VectorND<nsr>     s  = sections[i]->getResultant<nsr,scheme>();


    // Add material stiffness
    // kb += B' * (A' ks A) * B
    ks.addMatrixTripleProduct(0.0, A, Ks, 1.0);
    kb.addMatrixTripleProduct(1.0, B, ks, L*wt[i]);

    // Beam geometric stiffness matrix
    Gm( 1,  1) = Gm(2,  2) = Gm(3, 3) = Gm(4, 4) =  s[0] * 4.0 / 30.0; // 4/30*N
    Gm( 1,  3) = Gm(2,  4) = Gm(3, 1) = Gm(4, 2) = -s[0] / 30.0;       // -1/30*N
    Gm( 9,  8) = Gm(8,  9) =  s[1];                                    // Mz
    Gm( 9,  7) = Gm(7,  9) =  s[2];                                    // My
    Gm(10, 10)             =  s[3];                                    // W

    // Add geometric stiffness
//  kb.addMatrixTripleProduct(B, Gm, L*wt[i]);


    // assemble internal force vector p
    p.addMatrixTransposeVector(1.0, B, A*s, L*wt[i]);
  }

// TODO: add q0
//q[0] += q0[0];
//q[1] += q0[1];
//q[2] += q0[2];
//q[3] += q0[3];
//q[4] += q0[4];


  static Matrix wrapper(kb);

  return wrapper;
}

const Matrix &
EulerDeltaFrame3d::getInitialStiff()
{
  return this->getTangentStiff();
}


void
EulerDeltaFrame3d::zeroLoad()
{
  q0.zero();
  return;
}

int
EulerDeltaFrame3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr
      << "EulerDeltaFrame3d::addLoad() -- load type unknown for element with tag: "
      << this->getTag() << "\n";
  return -1;
}



void
EulerDeltaFrame3d::Print(OPS_Stream &s, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"nodes\": [" << node_tags(0) << ", " 
                        << node_tags(1)
      << "], ";

    // Mass
    double mass;
    if (getIntegral(Field::Density, State::Init, mass) == 0)
      s << "\"mass\": " << mass << ", ";
    else
      s << "\"massperlength\": " << density << ", ";

    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << sections[i]->getTag() << "\", ";
    s << "\"" << sections[numSections - 1]->getTag() << "\"], ";

    s << "\"integration\": ";
    beamInt->Print(s, flag);

    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nEulerDeltaFrame3d, element id:  " 
      << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << node_tags;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << "\n";
    s << "\tmass density:  " << density 
      << ", mass_flag: " << mass_flag << "\n";

    beamInt->Print(s, flag);

    for (int i = 0; i < numSections; i++) {
      sections[i]->Print(s, flag);
    }
  }

}


Response *
EulerDeltaFrame3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType", "EulerDeltaFrame3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // Compare argv[0] for known response types
  //

  // global force
    if (strcmp(argv[0],"forces") == 0 || 
        strcmp(argv[0],"force") == 0  ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0) {

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

// TODO(cmp)
//  theResponse = new ElementResponse(this, 1, P);

    // Local force
    }  else if (strcmp(argv[0],"localForce") == 0 || 
                strcmp(argv[0],"localForces") == 0) {

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

    theResponse = new ElementResponse(this, 2, Vector(12));

    // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 ||
             strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));

  } else if (strcmp(argv[0], "RayleighForces") == 0 ||
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, Vector(12));

  } else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  // section response -
  else if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = theCoordTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i] - sectionLoc);
          sectionNum  = i;
        }
      }

      output.tag("GaussPointOutput");
      output.attr("number", sectionNum + 1);
      output.attr("eta", xi[sectionNum] * L);

      theResponse = sections[sectionNum]->setResponse(&argv[2], argc - 2, output);
    }
  }

  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);

        theResponse =
            sections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

        output.endTag();

      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

        CompositeResponse *theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);

          Response *theSectionResponse =
              sections[i]->setResponse(&argv[1], argc - 1, output);

          output.endTag();

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

  output.endTag();
  return theResponse;
}

int
EulerDeltaFrame3d::getResponse(int responseID, Information &info)
{
  double N, V, M1, M2, T;
  double L        = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  if (responseID == 1)
    return info.setVector(this->getResistingForce());

  else if (responseID == 12)
    return info.setVector(this->getRayleighDampingForces());

  else if (responseID == 2) {
    // TODO(cmp)
//  P.Zero();
//  return info.setVector(P);
    return -1;
  }

  // Chord rotation
  else if (responseID == 3)
    return info.setVector(theCoordTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    // TODO(cmp)
//  const Matrix &kb = this->getInitialBasicStiff();
//  kb.Solve(q, ve);
//  vp = theCoordTransf->getBasicTrialDisp();
//  vp -= ve;
    return info.setVector(vp);
  }

  else if (responseID == 10) {
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return info.setVector(locs);
  }

  else if (responseID == 11) {
    double L = theCoordTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i] * L;
    return info.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return info.setID(tags);
  }

  else
    return -1;
}

int
EulerDeltaFrame3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc < 3)
      return -1;

    float sectionLoc = atof(argv[1]);

    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamInt->getSectionLocations(numSections, L, xi);

    sectionLoc /= L;

    float minDistance = fabs(xi[0] - sectionLoc);
    int sectionNum    = 0;
    for (int i = 1; i < numSections; i++) {
      if (fabs(xi[i] - sectionLoc) < minDistance) {
        minDistance = fabs(xi[i] - sectionLoc);
        sectionNum  = i;
      }
    }
    return sections[sectionNum]->setParameter(&argv[2], argc - 2, param);
  }

  // If the parameter belongs to a section or lower
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  else if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to every object
  int ok     = 0;
  int result = 0;

  for (int i = 0; i < numSections; i++) {
    ok = sections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
EulerDeltaFrame3d::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    density = info.theDouble;
    return 0;
  } else
    return -1;
}

int
EulerDeltaFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}


int
EulerDeltaFrame3d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
EulerDeltaFrame3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}