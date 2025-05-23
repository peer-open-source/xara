//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Three Dimensional Mixed Beam Column Element
//
//===----------------------------------------------------------------------===
//
// element MixedFrame3d $tag $iNode $jNode $numIntgrPts $secTag $transfTag 
//   <-mass $massDens> <-integration $intType> <-damp_flag $rFlag> <-geomNonlinear>
//
// Required Input Parameters:
//   $tag                   integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts           number of integration points along the element length
//   $secTag                identifier for previously-defined section object
//   $transfTag             identifier for previously-defined coordinate-transformation (FrameTransform3d) object
//
// Optional Input:
//   -mass $massDens
//       $massDens          element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)
//   -integration $intType
//       $intType           numerical integration type, options are Lobotto, Legendre, Radau, NewtonCotes, Trapezoidal (optional, default= Lobotto)
//   -damp_flag $rFlag
//       $rFlag             optional, default = 1
//                              rFlag = 0 no rayleigh damping
//                              rFlag = 1 include rayleigh damping (default)
//   -geomNonlinear            perform analysis with internal geometric nonlinearity
//
//
// References:
//   1. Bulent N. Alemdar and Donald W. White, "Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis," Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, "Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns," Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.
//
//===----------------------------------------------------------------------===//
//
// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign
//
#include <MixedFrame3d.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <interpolate/cbdi.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Node.h>
#include <Message.h>

#include <FrameSection.h>
#include <BeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <RegularizedHingeIntegration.h>

Matrix MixedFrame3d::theMatrix(NEGD, NEGD);
Vector MixedFrame3d::theVector(NEGD);
Matrix MixedFrame3d::transformNaturalCoords(NDM_NATURAL_WITH_TORSION, NDM_NATURAL_WITH_TORSION);
Matrix MixedFrame3d::transformNaturalCoordsT(NDM_NATURAL_WITH_TORSION, NDM_NATURAL_WITH_TORSION);


Vector* MixedFrame3d::ei_trial = nullptr;
Vector* MixedFrame3d::si_trial = nullptr;
MatrixND<NDM_SECTION, NDM_NATURAL>* MixedFrame3d::nldhat   = nullptr;
MatrixND<NDM_SECTION, NDM_NATURAL>* MixedFrame3d::nd1      = nullptr;
MatrixND<NDM_SECTION, NDM_NATURAL>* MixedFrame3d::nd2      = nullptr;
MatrixND<NDM_NATURAL, NDM_SECTION>* MixedFrame3d::nd1T     = nullptr;
MatrixND<NDM_NATURAL, NDM_SECTION>* MixedFrame3d::nd2T     = nullptr;


// Constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
MixedFrame3d::MixedFrame3d(int tag, std::array<int, 2>& nodes, 
                           int numSec, FrameSection** sec,
                           BeamIntegration& bi, FrameTransform3d& coordTransf,
                           double density, int damp, int geom)
 : FiniteElement<2,3,6>(tag, ELE_TAG_MixedFrame3d, nodes),
   beamIntegr(0),
   numSections(0),
   sections(nullptr),
   crdTransf(nullptr),
   damp_flag(damp), geom_flag(geom),
   rho(density),
   L0(0.0),
   itr(0),
   state_flag(0),
   V(NDM_NATURAL), committedV(NDM_NATURAL),
   qe_pres(NDM_NATURAL_WITH_TORSION), qe_past(NDM_NATURAL_WITH_TORSION),
   naturalForce(NDM_NATURAL),    commitedNaturalForce(NDM_NATURAL),
   lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
   sp(0),
   commitedHinv(NDM_NATURAL, NDM_NATURAL),
   Ki(nullptr),
   sr_trial(nullptr),es_trial(nullptr), fs_trial(nullptr),
   sr_past(nullptr), es_past(nullptr),  fs_past(nullptr)
{
  // get copy of the beam integration object
  beamIntegr = bi.getCopy();
  if (beamIntegr == nullptr) {
    opserr << "Error: MixedFrame3d::MixedFrame3d: could not create copy of beam integration object"
           << "\n";
    exit(-1);
  }

  // get copy of the transformation object
  crdTransf = coordTransf.getCopy();
  if (crdTransf == 0) {
    opserr << "Error: MixedFrame3d::MixedFrame3d: could not create copy of coordinate "
              "transformation object"
           << "\n";
  }


  //this->setSectionPointers(numSec,sec);
  if (numSec > MAX_NUM_SECTIONS) {
    opserr << "Error: MixedFrame3d::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (sec == 0) {
    opserr << "Error: MixedFrame3d::setSectionPointers -- invalid section pointer";
  }

  sections = new FrameSection*[numSections];

  for (int i = 0; i < numSections; i++) {
    if (sec[i] == nullptr) {
      opserr << "Error: MixedFrame3d::setSectionPointers -- null section pointer " << i << "\n";
    }

    sections[i] = (FrameSection*)sec[i]->getCopy();
  }

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // Element vectors and matrices
  sr_trial = new Vector[numSections];
  sr_past  = new Vector[numSections];
  es_trial = new Vector[numSections];
  es_past  = new Vector[numSections];
  fs_trial = new MatrixND<NDM_SECTION, NDM_SECTION>[numSections];
  fs_past  = new MatrixND<NDM_SECTION, NDM_SECTION>[numSections];

  for (int i = 0; i < numSections; i++) {
    sr_trial[i] = Vector(NDM_SECTION);
    sr_trial[i].Zero();
    sr_past[i] = Vector(NDM_SECTION);
    sr_past[i].Zero();
    es_trial[i] = Vector(NDM_SECTION);
    es_trial[i].Zero();
    es_past[i] = Vector(NDM_SECTION);
    es_past[i].Zero();
    fs_trial[i].zero();
    fs_past[i].zero();
  }

  V.Zero();
  qe_pres.Zero();
  naturalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.zero();
  GMH.zero();
  kv.zero();

  committedV.Zero();
  qe_past.Zero();
  commitedNaturalForce.Zero();
  commitedLastNaturalDisp.Zero();
  commitedHinv.Zero();
  commitedGMH.zero();
  ke_past.zero();

  if (transformNaturalCoords(1, 1) != 1) {
    // if transformNaturalCoords hasn't been set yet then set it
    // This is needed because the way the state determination algorithm was formulated
    // and OpenSees assume a different order of natural degrees of freedom
    // Formulation --> { e, theta_zi, theta_yi, theta_zj, theta_yj, twist }
    // OpenSees    --> { e, theta_zi, theta_zj, theta_yi, theta_yj, twist }
    transformNaturalCoords.Zero();
    transformNaturalCoords(0, 0) = 1;
    transformNaturalCoords(1, 1) = 1;
    transformNaturalCoords(2, 3) = 1;
    transformNaturalCoords(3, 2) = 1;
    transformNaturalCoords(4, 4) = 1;
    transformNaturalCoords(5, 5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0, 0) = 1;
    transformNaturalCoordsT(1, 1) = 1;
    transformNaturalCoordsT(3, 2) = 1;
    transformNaturalCoordsT(2, 3) = 1;
    transformNaturalCoordsT(4, 4) = 1;
    transformNaturalCoordsT(5, 5) = 1;
  }

  if (ei_trial == nullptr)
    ei_trial = new Vector[MAX_NUM_SECTIONS];
  if (si_trial == nullptr)
    si_trial = new Vector[MAX_NUM_SECTIONS];
  if (nldhat == nullptr)
    nldhat = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd1 == nullptr)
    nd1 = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd2 == nullptr)
    nd2 = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd1T == nullptr)
    nd1T = new MatrixND<NDM_NATURAL,NDM_SECTION>[MAX_NUM_SECTIONS]{};
  if (nd2T == nullptr)
    nd2T = new MatrixND<NDM_NATURAL,NDM_SECTION>[MAX_NUM_SECTIONS]{};

//for (int i = 0; i < MAX_NUM_SECTIONS; i++) {
//  nd1T[i] = Matrix(NDM_NATURAL, NDM_SECTION);
//  nd2T[i] = Matrix(NDM_NATURAL, NDM_SECTION);
//}
}


// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING
MixedFrame3d::MixedFrame3d()
 : FiniteElement<2,3,6>(0, ELE_TAG_MixedFrame3d),
   beamIntegr(0),
   numSections(0),
   sections(nullptr),
   crdTransf(nullptr),
   damp_flag(0),
   geom_flag(true),
   rho(0.0),
   L0(0.0),
   itr(0),
   state_flag(0),
   V(NDM_NATURAL),
   committedV(NDM_NATURAL),
   qe_pres(NDM_NATURAL_WITH_TORSION),
   qe_past(NDM_NATURAL_WITH_TORSION),
   naturalForce(NDM_NATURAL),
   commitedNaturalForce(NDM_NATURAL),
   lastNaturalDisp(NDM_NATURAL),
   commitedLastNaturalDisp(NDM_NATURAL),
   sp(nullptr),
   commitedHinv(NDM_NATURAL, NDM_NATURAL),
   Ki(0),
   sr_trial(nullptr),
   sr_past(nullptr),
   es_trial(nullptr),
   es_past(nullptr),
   fs_trial(nullptr),
   fs_past(nullptr)
{

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // Element vectors and matrices
  sr_trial = new Vector[numSections];
  sr_past  = new Vector[numSections];
  es_trial = new Vector[numSections];
  es_past  = new Vector[numSections];
  fs_trial = new MatrixND<NDM_SECTION, NDM_SECTION>[numSections];
  fs_past  = new MatrixND<NDM_SECTION, NDM_SECTION>[numSections];

  for (int i = 0; i < numSections; i++) {
    sr_trial[i] = Vector(NDM_SECTION);
    sr_trial[i].Zero();
    sr_past[i] = Vector(NDM_SECTION);
    sr_past[i].Zero();
    es_trial[i] = Vector(NDM_SECTION);
    es_trial[i].Zero();
    es_past[i] = Vector(NDM_SECTION);
    es_past[i].Zero();
    fs_trial[i].zero();
    fs_past[i].zero();
  }

  V.Zero();
  qe_pres.Zero();
  naturalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.zero();
  GMH.zero();
  kv.zero();

  committedV.Zero();
  qe_past.Zero();
  commitedNaturalForce.Zero();
  commitedLastNaturalDisp.Zero();
  commitedHinv.Zero();
  commitedGMH.zero();
  ke_past.zero();


  if (transformNaturalCoords(1, 1) != 1) {
    // if transformNaturalCoords hasn't been set yet then set it
    static constexpr MatrixND<6,6> T {{
    // 0  1  2  3  4  5
      {1, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 0},
      {0, 0,-1, 0, 0, 0},
      {0, 0, 0, 0, 1, 0},
      {0, 0, 0, 0, 0, 1}
    }};
    transformNaturalCoords.Zero();
    transformNaturalCoords(0, 0) = 1;
    transformNaturalCoords(1, 1) = 1;
    transformNaturalCoords(2, 3) = -1;
    transformNaturalCoords(3, 2) = 1;
    transformNaturalCoords(4, 4) = -1;
    transformNaturalCoords(5, 5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0, 0) = 1;
    transformNaturalCoordsT(1, 1) = 1;
    transformNaturalCoordsT(3, 2) = -1;
    transformNaturalCoordsT(2, 3) = 1;
    transformNaturalCoordsT(4, 4) = -1;
    transformNaturalCoordsT(5, 5) = 1;
  }

  if (ei_trial == nullptr)
    ei_trial = new Vector[MAX_NUM_SECTIONS];
  if (si_trial == nullptr)
    si_trial = new Vector[MAX_NUM_SECTIONS];
  if (nldhat == nullptr)
    nldhat = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd1 == nullptr)
    nd1 = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd2 == nullptr)
    nd2 = new MatrixND<NDM_SECTION, NDM_NATURAL>[MAX_NUM_SECTIONS];
  if (nd1T == nullptr)
    nd1T = new MatrixND<NDM_NATURAL,NDM_SECTION>[MAX_NUM_SECTIONS]{};
  if (nd2T == nullptr)
    nd2T = new MatrixND<NDM_NATURAL,NDM_SECTION>[MAX_NUM_SECTIONS]{};

}

MixedFrame3d::~MixedFrame3d()
{

  if (sections) {
    for (int i = 0; i < numSections; i++) {
      if (sections[i]) {
        delete sections[i];
      }
    }
    delete[] sections;
  }

  if (crdTransf)
    delete crdTransf;

  if (beamIntegr != nullptr)
    delete beamIntegr;

  if (sp != nullptr)
    delete sp;

  if (Ki != nullptr)
    delete Ki;

  if (sr_trial != nullptr)
    delete[] sr_trial;

  if (sr_past != 0)
    delete[] sr_past;

  if (es_trial != nullptr)
    delete[] es_trial;

  if (es_past != 0)
    delete[] es_past;

  if (fs_trial != nullptr)
    delete[] fs_trial;

  if (fs_past != 0)
    delete[] fs_past;
}


int
MixedFrame3d::setNodes() // (Domain* theDomain)
{

  if (crdTransf->initialize(theNodes[0], theNodes[1]) != 0) {
      opserr << "BasicFrame3d::setDomain  tag: " 
            << this->getTag()
            << " -- Error initializing coordinate transformation\n";
      return -1;
  }

  // call the DomainComponent class method
//this->DomainComponent::setDomain(theDomain);

  // Get the numerical integration weights
  beamIntegr->getSectionWeights(numSections, L0, wt);
  beamIntegr->getSectionLocations(numSections, L0, xi);

  return 0;
}

int
MixedFrame3d::commitState()
{
  int err = 0; // error flag
  int i   = 0; // integer for loops

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "MixedFrame3d::commitState () - failed in base class";
    return err;
  }

  // commit the sections
  do {
    err = sections[i++]->commitState();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;

  // commit the element variables state
  committedV                    = V;
  qe_past = qe_pres;
  commitedNaturalForce          = naturalForce;
  commitedLastNaturalDisp       = lastNaturalDisp;
  commitedHinv                  = Hinv;
  commitedGMH                   = GMH;
  ke_past                       = kv;
  for (i = 0; i < numSections; i++) {
    sr_past[i] = sr_trial[i];
    es_past[i]   = es_trial[i];
    fs_past[i] = fs_trial[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int
MixedFrame3d::revertToLastCommit()
{
  int err;
  int i = 0;

  do {
    err = sections[i]->revertToLastCommit();
    i++;
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  // revert the element state to last commit
  V                     = committedV;
  qe_pres = qe_past;
  naturalForce          = commitedNaturalForce;
  lastNaturalDisp       = commitedLastNaturalDisp;
  Hinv                  = commitedHinv;
  GMH                   = commitedGMH;
  kv                    = ke_past;
  for (i = 0; i < numSections; i++) {
    sr_trial[i] = sr_past[i];
    es_trial[i]   = es_past[i];
    fs_trial[i] = fs_past[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int
MixedFrame3d::revertToStart()
{
  int err;

  // revert the sections state to start
  int i = 0;
  do {
    err = sections[i++]->revertToStart();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start

  // Set initial length
  L0 = crdTransf->getInitialLength();

  // Vector of zeros to use at initial natural displacements
  Vector myZeros(NDM_NATURAL);
  myZeros.Zero();

  // Set initial shape functions
  for (i = 0; i < numSections; i++) {
    nldhat[i] = this->getNld_hat(i, myZeros, L0, geom_flag);
    nd1[i]    = this->getNd1(i, myZeros, L0, geom_flag);
    nd2[i]    = this->getNd2(i, 0, L0);

    for (int j = 0; j < NDM_SECTION; j++) {
      for (int k = 0; k < NDM_NATURAL; k++) {
        nd1T[i](k, j) = nd1[i](j, k);
        nd2T[i](k, j) = nd2[i](j, k);
      }
    }
  }

  // Set initial and committed section flexibility and GJ
  MatrixND<NDM_SECTION, NDM_SECTION> Ks;
  double GJ;
  for (int i = 0; i < numSections; i++) {
    getSectionTangent(i, 2, Ks, GJ);
    Ks.invert(fs_trial[i]);
    fs_past[i] = fs_trial[i];
  }

  // Set initial and committed section forces and deformations
  for (int i = 0; i < numSections; i++) {
    sr_trial[i].Zero();
    sr_past[i].Zero();
    es_trial[i].Zero();
    es_past[i].Zero();
  }

  // Compute the following matrices: G, G2, H, H12, H22, Md, Kg
  MatrixND<NDM_NATURAL, NDM_NATURAL> G{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> G2{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H12{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H22{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> Md{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> Kg{};

  for (int i = 0; i < numSections; i++) {
    G   += L0 * wt[i] * nd1T[i] * nldhat[i];
    G2  += L0 * wt[i] * nd2T[i] * nldhat[i];
    H   += L0 * wt[i] * nd1T[i] * fs_trial[i] * nd1[i];
    H12 += L0 * wt[i] * nd1T[i] * fs_trial[i] * nd2[i];
    H22 += L0 * wt[i] * nd2T[i] * fs_trial[i] * nd2[i];
    // Md is zero since deformations are zero
    Kg.addMatrix(this->getKg(i, 0.0, L0),  L0 * wt[i]);
  }

  // Compute the inverse of the H matrix
  {
    // Matrix hinv(Hinv);
    H.invert(Hinv);
    commitedHinv = Hinv;
  }

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
  commitedGMH = GMH;

  // Compute the transposes of the following matrices: G2, GMH
  MatrixND<NDM_NATURAL, NDM_NATURAL> G2T;
  MatrixND<NDM_NATURAL, NDM_NATURAL> GMHT;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++) {
      G2T(i, j)  = G2(j, i);
      GMHT(i, j) = GMH(j, i);
    }
  }

  // Compute the stiffness matrix without the torsion term
  MatrixND<NDM_NATURAL, NDM_NATURAL> 
    K_temp_noT = (Kg + G2 + G2T - H22) + GMHT * Hinv * GMH;
  //K_temp_noT = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta

  // Add in the torsional stiffness term
  kv.zero();
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++) {
      kv(i, j) = K_temp_noT(i, j);
    }
  }

  kv(5, 5) = GJ / L0; // Torsional Stiffness GJ/L
  ke_past = kv;

  Matrix kvOpenSees = transformNaturalCoordsT * kv * transformNaturalCoords;
  if (Ki == nullptr)
    Ki = new Matrix(NEGD, NEGD);

  *Ki = crdTransf->getInitialGlobalStiffMatrix(kvOpenSees);

  // Vector V is zero at initial state
  V.Zero();
  committedV.Zero();

  // Internal force is zero at initial state
  qe_pres.Zero();
  qe_past.Zero();
  naturalForce.Zero();
  commitedNaturalForce.Zero();

  // Last natural displacement is zero at initial state
  lastNaturalDisp.Zero();
  commitedLastNaturalDisp.Zero();

  // Reset iteration counter
  itr = 0;

  // Set state_flag to 1 so update doesn't call again
  state_flag = 1;

  return err;
}

const Matrix&
MixedFrame3d::getInitialStiff()
{
  if (state_flag == 0)
    this->revertToStart();

  return *Ki;
}

const Matrix&
MixedFrame3d::getTangentStiff()
{

  if (state_flag == 0)
    this->revertToStart();

  Matrix ktOpenSees = transformNaturalCoordsT * kv * transformNaturalCoords;
  return crdTransf->getGlobalStiffMatrix(ktOpenSees, qe_pres);
}

const Vector&
MixedFrame3d::getResistingForce()
{
  // crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0Vec(p0, NDM_NATURAL);
  return crdTransf->getGlobalResistingForce(qe_pres, p0Vec);
}

int
MixedFrame3d::update()
{

  // If things haven't be initialized, then do so
  if (state_flag == 0)
    this->revertToStart();

  // Update iteration counter
  // says how many times update has been called since the last commit state
  itr++;

  // Update Coordinate Transformation
  crdTransf->update();

  // Current Length
  double L;
  if (geom_flag == Geometry::Linear) {
    L = L0;
  } else {
    L = crdTransf->getDeformedLength();
  }

  // Compute the natural displacements
  Vector naturalDispWithTorsion = crdTransf->getBasicTrialDisp();
  naturalDispWithTorsion        = transformNaturalCoords * naturalDispWithTorsion;

  // convert to the arrangement of natural deformations that the element likes
  Vector naturalDisp(NDM_NATURAL);
  for (int i = 0; i < NDM_NATURAL; i++)
    naturalDisp(i) = naturalDispWithTorsion(i); // all but the torsional component

  double twist = naturalDispWithTorsion(5);

  Vector naturalIncrDeltaDisp(NDM_NATURAL);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp      = naturalDisp;


  for (int i = 0; i < numSections; i++)
    si_trial[i] = Vector(NDM_SECTION);

  // Compute shape functions and their transposes
  for (int i = 0; i < numSections; i++) {
    // Shape Functions
    nldhat[i]   = this->getNld_hat(i, naturalDisp, L, geom_flag);
    ei_trial[i] = this->getd_hat(i, naturalDisp, L, geom_flag);
    nd1[i]      = this->getNd1(i, naturalDisp, L, geom_flag);

    if (geom_flag == Geometry::Linear) {
      nd2[i].zero();
    } else {
      nd2[i] = this->getNd2(i, qe_pres(0), L);
    }

    // Transpose of shape functions
    for (int j = 0; j < NDM_SECTION; j++) {
      for (int k = 0; k < NDM_NATURAL; k++) {
        nd1T[i](k, j) = nd1[i](j, k);
        nd2T[i](k, j) = nd2[i](j, k);
      }
    }
  }

  // Update natural force
  if (geom_flag == Geometry::Linear) {
    naturalForce += Hinv * (GMH * naturalIncrDeltaDisp + V);
  } else {
    naturalForce += Hinv * (GMH * naturalIncrDeltaDisp + V);
  }


  double GJ, torsionalForce;
  // Update sections
  for (int i = 0; i < numSections; i++) {
    // Compute section deformations
    si_trial[i] = nd1[i] * naturalForce;
    if (sp != 0) {
      const Matrix& s_p = *sp;
      for (int j = 0; j < NDM_SECTION; j++) {
        si_trial[i](j) += s_p(j, i);
      }
    }
    es_trial[i] += fs_trial[i] * (si_trial[i] - sr_trial[i]);

    // Send section deformation to section object
    double torsionalStrain = twist / L;
    setSectionDeformation(i, es_trial[i], torsionalStrain);

    // Get section force vector
    sr_trial[i] = sections[i]->getResultant<NDM_SECTION, scheme>();
    fs_trial[i] = sections[i]->getFlexibility<NDM_SECTION, scheme>();

    torsionalForce = sr_trial[i](3);

  }

  VectorND<NDM_NATURAL> V2{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> G{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> G2{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H12{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> H22{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> Md{};
  MatrixND<NDM_NATURAL, NDM_NATURAL> Kg{};

  V.Zero();
  V2.zero();
  G.zero();
  G2.zero();
  H.zero();
  H12.zero();
  H22.zero();
  Md.zero();
  Kg.zero();

  for (int i = 0; i < numSections; i++) {
    V = V + L0*wt[i] * nd1T[i]*(ei_trial[i] - es_trial[i] -
                 fs_trial[i] * (si_trial[i] - sr_trial[i]));
    V2  = V2  + L0 * wt[i] * nd2T[i] * (ei_trial[i] - es_trial[i]);
    G   = G   + L0 * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + L0 * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + L0 * wt[i] * nd1T[i] * fs_trial[i] * nd1[i];
    H12 = H12 + L0 * wt[i] * nd1T[i] * fs_trial[i] * nd2[i];
    H22 = H22 + L0 * wt[i] * nd2T[i] * fs_trial[i] * nd2[i];

    if (!geom_flag) {
      // sr_trial[i][0] is the axial load, P
      Kg += L0 * wt[i] * this->getKg(i, sr_trial[i][0], L);
      Md += L0 * wt[i] * this->getMd(i, ei_trial[i], es_trial[i], L);
    }
  }

  // Compute the inverse of the H matrix
  H.invert(Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta

  // Compute the transposes of the following matrices: G, G2, GMH
  MatrixND<NDM_NATURAL,NDM_NATURAL> GT;
  MatrixND<NDM_NATURAL,NDM_NATURAL> G2T;
  MatrixND<NDM_NATURAL,NDM_NATURAL> GMHT;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++) {
      GT(i, j)   = G(j, i);
      G2T(i, j)  = G2(j, i);
      GMHT(i, j) = GMH(j, i);
    }
  }

  // Compute new internal force
  Vector internalForce(NDM_NATURAL);
  internalForce.Zero();

  if (geom_flag == Geometry::Linear) {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  } else {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  }

  // Compute internal force for OpenSees ( i.e., add torsion and rearrange )
  for (int i = 0; i < NDM_NATURAL; i++)
    qe_pres(i) = internalForce[i];

  qe_pres(5) = torsionalForce; // Add in torsional force
  qe_pres    = transformNaturalCoordsT * qe_pres;

  // Compute the stiffness matrix without the torsion term
  MatrixND<NDM_NATURAL, NDM_NATURAL> K_temp;
  if (geom_flag == Geometry::Linear) {
    K_temp = (Kg + G2 + G2T - H22) + GMHT * Hinv * GMH;
  } else {
    K_temp = (Kg + G2 + G2T - H22) + GMHT * Hinv * GMH;
  }

  // Add in the torsional stiffness term
  kv.zero();
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++) {
      kv(i, j) = K_temp(i, j);
    }
  }

  kv(5, 5) = GJ / L; // Torsional Stiffness GJ/L

  return 0;
}



const Matrix&
MixedFrame3d::getDamp()
{
  theMatrix.Zero();

  // Add the damping forces
  if (damp_flag == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)) {
    theMatrix = this->Element::getDamp();
  }

  return theMatrix;
}

void
MixedFrame3d::zeroLoad()
{
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
}

int
MixedFrame3d::addLoad(ElementalLoad* theLoad, double loadFactor)
{

  int type;
  const Vector& data = theLoad->getData(type, loadFactor);

  if (sp == 0) {
    sp = new Matrix(NDM_SECTION, numSections);
    if (sp == 0) {
      opserr << "MixedFrame3d::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0) * loadFactor; // Transverse
    double wz = data(1) * loadFactor; // Transverse
    double wx = data(2) * loadFactor; // Axial

    Matrix& s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i] * L;
      // Axial
      s_p(0, i) += wx * (L - x);
      // Moment
      s_p(1, i) += wy * 0.5 * x * (x - L);
      // Moment
      s_p(2, i) += wz * 0.5 * x * (L - x);
    }

    // Accumulate reactions in basic system
    p0[0] -= wx * L;
    double V = 0.5 * wy * L;
    p0[1] -= V;
    p0[2] -= V;
    V = 0.5 * wz * L;
    p0[3] -= V;
    p0[4] -= V;


  } else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py     = data(0) * loadFactor;
    double Pz     = data(1) * loadFactor;
    double N      = data(2) * loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL * L;

    double Vy2 = Py * aOverL;
    double Vy1 = Py - Vy2;

    double Vz2 = Pz * aOverL;
    double Vz1 = Pz - Vz2;

    Matrix& s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i] * L;
      if (x <= a) {
        s_p(0, i) += N;
        s_p(1, i) -= x * Vy1;
        s_p(2, i) += x * Vz1;
      } else {
        s_p(1, i) -= (L - x) * Vy2;
        s_p(2, i) += (L - x) * Vz2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;


  } else {
    opserr << "MixedFrame3d::addLoad() -- load type unknown for element with tag: "
           << this->getTag() << "\n";

    return -1;
  }

  return 0;
}

const Vector&
MixedFrame3d::getResistingForceIncInertia()
{

  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Add the inertial forces
  if (rho != 0.0) {
    const Vector& accel1 = theNodes[0]->getTrialAccel();
    const Vector& accel2 = theNodes[1]->getTrialAccel();

    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L;

    theVector(0) += m * accel1(0);
    theVector(1) += m * accel1(1);
    theVector(2) += m * accel1(2);
    theVector(6) += m * accel2(0);
    theVector(7) += m * accel2(1);
    theVector(8) += m * accel2(2);
  }

  // Add the damping forces
  if (damp_flag == Damping::Rayleigh && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)) {
    theVector += this->getRayleighDampingForces();
  }

  return theVector;
}


Vector
MixedFrame3d::getd_hat(int sec, const Vector& v, double L, int geom_flag)
{
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  Vector D_hat(NDM_SECTION);
  D_hat.Zero();

  double x = L * xi[sec];
  double C = 1 / L;
  double E = -4 / L + 6 * x / (L * L);
  double F = -2 / L + 6 * x / (L * L);

  if (geom_flag == Geometry::Linear) {

    D_hat(0) = C * v(0);
    D_hat(1) = E * v(1) + F * v(3);
    D_hat(2) = E * v(2) + F * v(4);

  } else {

    double A, B;
    A = 1 - 4 * (x / L) + 3 * pow(x / L, 2);
    B = -2 * (x / L) + 3 * pow(x / L, 2);

    D_hat(0) = C * v(0) + 0.5 * (C * C * v(0)) * v(0) + 0.5 * (A * A * v(1) + A * B * v(3)) * v(1) +
               0.5 * (A * A * v(2) + A * B * v(4)) * v(2) +
               0.5 * (A * B * v(1) + B * B * v(3)) * v(3) +
               0.5 * (A * B * v(2) + B * B * v(4)) * v(4);
    D_hat(1) = E * v(1) + F * v(3);
    D_hat(2) = E * v(2) + F * v(4);
  }

  return D_hat;
}

MatrixND<NDM_NATURAL, NDM_NATURAL>
MixedFrame3d::getKg(int sec, double P, double L)
{
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x = L * xi[sec];

  MatrixND<NDM_NATURAL, NDM_NATURAL> kg;
  kg.zero();

  double temp_A = 1 - 4 * temp_x / L + 3 * (temp_x * temp_x) / (L * L);
  double temp_B = -2 * temp_x / L + 3 * (temp_x * temp_x) / (L * L);

  kg(0, 0) = P / (L * L);
  kg(1, 1) = P * temp_A * temp_A;
  kg(1, 3) = P * temp_A * temp_B;
  kg(2, 2) = P * temp_A * temp_A;
  kg(2, 4) = P * temp_A * temp_B;
  kg(3, 1) = P * temp_A * temp_B;
  kg(3, 3) = P * temp_B * temp_B;
  kg(4, 2) = P * temp_A * temp_B;
  kg(4, 4) = P * temp_B * temp_B;

  return kg;
}

MatrixND<NDM_NATURAL, NDM_NATURAL> 
MixedFrame3d::getMd(int sec, Vector& dShapeFcn, Vector& dFibers, double L)
{
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);


  MatrixND<NDM_NATURAL, NDM_NATURAL> md{};

  double x = L * xi[sec];
  double A = (x / L - 2 * pow(x / L, 2) + pow(x / L, 3)) * L;
  double B = (-pow(x / L, 2) + pow(x / L, 3)) * L;

  md(0, 1) = A * (dShapeFcn[1] - dFibers[1]);
  md(0, 2) = A * (dShapeFcn[2] - dFibers[2]);
  md(0, 3) = B * (dShapeFcn[1] - dFibers[1]);
  md(0, 4) = B * (dShapeFcn[2] - dFibers[2]);

  return md;
}

MatrixND<NDM_SECTION, NDM_NATURAL>
MixedFrame3d::getNld_hat(int sec, const Vector& v, double L, int geom_flag)
{
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  MatrixND<NDM_SECTION, NDM_NATURAL> Nld_hat;
  Nld_hat.zero();

  double x = L * xi[sec];

  double C = 1 / L;
  double E = -4 / L + 6 * x / (L * L);
  double F = -2 / L + 6 * x / (L * L);

  if (geom_flag == Geometry::Linear) {

    Nld_hat(0, 0) = C;
    Nld_hat(1, 1) = E;
    Nld_hat(1, 3) = F;
    Nld_hat(2, 2) = E;
    Nld_hat(2, 4) = F;

  } else {

    double A = 1 - 4 * (x / L) + 3 * pow((x / L), 2);
    double B = -2 * (x / L) + 3 * pow((x / L), 2);

    Nld_hat(0, 0) = C + C * C * v(0);
    Nld_hat(0, 1) = A * A * v(1) + A * B * v(3);
    Nld_hat(0, 2) = A * A * v(2) + A * B * v(4);
    Nld_hat(0, 3) = A * B * v(1) + B * B * v(3);
    Nld_hat(0, 4) = A * B * v(2) + B * B * v(4);
    Nld_hat(1, 1) = E;
    Nld_hat(1, 3) = F;
    Nld_hat(2, 2) = E;
    Nld_hat(2, 4) = F;
  }

  return Nld_hat;
}

MatrixND<NDM_SECTION, NDM_NATURAL>
MixedFrame3d::getNd2(int sec, double P, double L)
{

  double x = L * xi[sec];

  MatrixND<NDM_SECTION, NDM_NATURAL> Nd2;
  Nd2.zero();

  double A = L * (x / L - 2 * pow(x / L, 2) + pow(x / L, 3));
  double B = L * (-pow(x / L, 2) + pow(x / L, 3));

  Nd2(1, 1) = P * A;
  Nd2(1, 3) = P * B;
  Nd2(2, 2) = P * A;
  Nd2(2, 4) = P * B;

  return Nd2;
}

Matrix
MixedFrame3d::getNd1(int sec, const Vector& v, double L, int geom_flag)
{

  double x = L * xi[sec];

  MatrixND<NDM_SECTION, NDM_NATURAL> Nd1;
  Nd1.zero();

  if (geom_flag == Geometry::Linear) {

    Nd1(0, 0) = 1.0;
    Nd1(1, 1) = -x / L + 1.0;
    Nd1(1, 3) = x / L;
    Nd1(2, 2) = -x / L + 1.0;
    Nd1(2, 4) = x / L;

  } else {

    double A, B;

    A = L * (x / L - 2 * pow(x / L, 2) + pow(x / L, 3)) * v[1] +
        L * (-pow(x / L, 2) + pow(x / L, 3)) * v[3];

    B = L * (x / L - 2 * pow(x / L, 2) + pow(x / L, 3)) * v[2] +
        L * (-pow(x / L, 2) + pow(x / L, 3)) * v[4];

    Nd1(0, 0) = 1.0;
    Nd1(1, 0) = A;
    Nd1(1, 1) = -x / L + 1.0;
    Nd1(1, 3) = x / L;
    Nd1(2, 0) = B;
    Nd1(2, 2) = -x / L + 1.0;
    Nd1(2, 4) = x / L;
  }

  return Nd1;
}

void
MixedFrame3d::getSectionTangent(int sec, int type, MatrixND<NDM_SECTION, NDM_SECTION>& kSection, double& GJ)
{
  int order      = sections[sec]->getOrder();
  const ID& code = sections[sec]->getType();

  // Initialize formulation friendly variables
  kSection.zero();
  GJ = 0.0;

  // Get the stress resultant from section
  Matrix sectionTangent(order, order);
  if (type == 1) {
    sectionTangent = sections[sec]->getSectionTangent();
  } else if (type == 2) {
    sectionTangent = sections[sec]->getInitialTangent();
  } else {
    sectionTangent.Zero();
  }

  // Set Components of Section Tangent
  for (int i = 0; i < order; i++) {
    for (int j = 0; j < order; j++) {
      switch (code(i)) {
      case SECTION_RESPONSE_P:
        switch (code(j)) {
        case SECTION_RESPONSE_P:  kSection(0, 0) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MZ: kSection(0, 1) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MY: kSection(0, 2) = sectionTangent(i, j); break;
        default:                  break;
        }
        break;
      case SECTION_RESPONSE_MZ:
        switch (code(j)) {
        case SECTION_RESPONSE_P:  kSection(1, 0) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MZ: kSection(1, 1) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MY: kSection(1, 2) = sectionTangent(i, j); break;
        default:                  break;
        }
        break;
      case SECTION_RESPONSE_MY:
        switch (code(j)) {
        case SECTION_RESPONSE_P:  kSection(2, 0) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MZ: kSection(2, 1) = sectionTangent(i, j); break;
        case SECTION_RESPONSE_MY: kSection(2, 2) = sectionTangent(i, j); break;
        default:                  break;
        }
        break;
      case SECTION_RESPONSE_T: GJ = sectionTangent(i, i); break;
      default:                 break;
      }
    }
  }
}


void
MixedFrame3d::setSectionDeformation(int sec, Vector& defSection, double& twist)
{
  int order      = sections[sec]->getOrder();
  const ID& code = sections[sec]->getType();

  // Initialize Section Deformation Vector
  Vector sectionDeformation(order);
  sectionDeformation.Zero();

  // Set Components of Section Deformations
  for (int j = 0; j < order; j++) {
    switch (code(j)) {
    case SECTION_RESPONSE_P:  sectionDeformation(j) = defSection(0); break;
    case SECTION_RESPONSE_MZ: sectionDeformation(j) = defSection(1); break;
    case SECTION_RESPONSE_MY: sectionDeformation(j) = defSection(2); break;
    case SECTION_RESPONSE_T:  sectionDeformation(j) = twist; break;
    default:                  break;
    }
  }

  // Set the section deformations
  int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}



void
MixedFrame3d::Print(OPS_Stream& s, int flag)
{

  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: MixedFrame3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho;
    for (int i = 0; i < numSections; i++)
      s << "\nSection " << i << " :" << *sections[i];

  } else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: MixedFrame3d ";
    double xi[MAX_NUM_SECTIONS]; // location of sections or gauss points or integration points
    beamIntegr->getSectionLocations(numSections, L0, xi);
    double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
    beamIntegr->getSectionWeights(numSections, L0, wt);
    s << "\n section xi wt";
    for (int i = 0; i < numSections; i++)
      s << "\n" << i << " " << xi[i] << " " << wt[i];

  } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"mixedFrame2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " 
                        << connectedExternalNodes(1) 
      << "], ";

    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << sections[i]->getTag() << "\", ";
    s << "\"" << sections[numSections - 1]->getTag() << "\"], ";

    s << "\"integration\": ";
    beamIntegr->Print(s, flag);
    s << ", \"massperlength\": " << rho << ", ";
    s << "\"transform\": \"" << crdTransf->getTag() << "\"";
    if (!damp_flag)
      s << ", \"damp_flag\": false";
    if (geom_flag == Geometry::Linear)
      s << ", \"geom_flag\": true";
    s << "}";

  } else {
    s << "\nElement: " << this->getTag() << " Type: MixedFrame3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << "\n";
  }
}


Response*
MixedFrame3d::setResponse(const char** argv, int argc, OPS_Stream& output)
{

  Response* theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "MixedFrame3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

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

    theResponse = new ElementResponse(this, 1, theVector);

    // local force -
  } else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_ 1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Tz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 2, theVector);

    // basic or natural forces
  } else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 3, Vector(6));
  }
#if 0
  else if (strcmp(argv[0], "sectionDeformation_Force") == 0) {

    int i;
    char* q = new char[15];
    for (i = 0; i < numSections; i++) {
      sprintf(q, "axialStrain_%i", i + 1);
      output.tag("ResponseType", q);
      sprintf(q, "curvatureZ_%i", i + 1);
      output.tag("ResponseType", q);
      sprintf(q, "curvatureY_%i", i + 1);
      output.tag("ResponseType", q);
    }
    delete[] q;

    theResponse = new ElementResponse(this, 4, Vector(3 * numSections));
  } 

  else if (strcmp(argv[0], "plasticSectionDeformation_Force") == 0) {

    int i;
    char* q = new char[25];
    for (i = 0; i < numSections; i++) {
      sprintf(q, "plasticAxialStrain_%i", i + 1);
      output.tag("ResponseType", q);
      sprintf(q, "plasticCurvatureZ_%i", i + 1);
      output.tag("ResponseType", q);
      sprintf(q, "plasticCurvatureY_%i", i + 1);
      output.tag("ResponseType", q);
    }
    delete[] q;

    theResponse = new ElementResponse(this, 5, Vector(3 * numSections));
  } 
#endif
  else if (strcmp(argv[0], "integrationPoints") == 0) {
    theResponse = new ElementResponse(this, 100, Vector(numSections));

  } else if (strcmp(argv[0], "integrationWeights") == 0) {
    theResponse = new ElementResponse(this, 101, Vector(numSections));

  } else if (strcmp(argv[0], "sectionTags") == 0) {
    theResponse = new ElementResponse(this, 110, ID(numSections));

  } else if (strcmp(argv[0], "connectedNodes") == 0) {
    theResponse = new ElementResponse(this, 102, Vector(2));

  } else if (strcmp(argv[0], "numSections") == 0 || strcmp(argv[0], "numberOfSections") == 0) {
    theResponse = new ElementResponse(this, 103, Vector(1));

  }

  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(numSections, 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(numSections, 3));
  }

  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 2) {

      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

        double xi[MAX_NUM_SECTIONS];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);

        theResponse = sections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

        output.endTag();
      }
    }
  }

  if (theResponse == nullptr)
    theResponse = crdTransf->setResponse(argv, argc, output);

  output.endTag();
  return theResponse;
}


int
MixedFrame3d::getResponse(int responseID, Information& info)
{
  if (responseID == 1) { // global forces
    return info.setVector(this->getResistingForce());

  } else if (responseID == 2) { // local forces
    // Axial
    double N     = qe_pres(0);
    theVector(6) = N;
    theVector(0) = -N + p0[0];

    // Torsion
    double T     = qe_pres(5);
    theVector(9) = T;
    theVector(3) = -T;

    // Moments about z and shears along y
    double M1     = qe_pres(1);
    double M2     = qe_pres(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L      = crdTransf->getInitialLength();
    double V      = (M1 + M2) / L;
    theVector(1)  = V + p0[1];
    theVector(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1            = qe_pres(3);
    M2            = qe_pres(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V             = -(M1 + M2) / L;
    theVector(2)  = -V + p0[3];
    theVector(8)  = V + p0[4];

    return info.setVector(theVector);

  } else if (responseID == 3) { // basic forces
    return info.setVector(qe_pres);

  } else if (responseID == 4) { // section deformation (from forces)

    int i;
    Vector tempVector(3 * numSections);
    tempVector.Zero();
    for (i = 0; i < numSections; i++) {
      tempVector(3 * i)     = es_trial[i](0);
      tempVector(3 * i + 1) = es_trial[i](1);
      tempVector(3 * i + 2) = es_trial[i](2);
    }

    return info.setVector(tempVector);

  } else if (responseID == 5) { // plastic section deformation (from forces)

    int i;
    Vector tempVector(3 * numSections);
    Vector sectionForce(NDM_SECTION);
    Vector plasticSectionDef(NDM_SECTION);
    Matrix ks(3, 3);
    Matrix fs(3, 3);
    tempVector.Zero();
    double scratch = 0.0;
    for (i = 0; i < numSections; i++) {

      sectionForce = sections[i]->getResultant<nsr, scheme>();
//    getSectionTangent(i, 2, ks, scratch);
      ks = sections[i]->getTangent<nsr, scheme>(State::Pres);
      fs = sections[i]->getFlexibility<nsr, scheme>();

      plasticSectionDef = es_trial[i] - fs * sectionForce;

      tempVector(3 * i)     = plasticSectionDef(0);
      tempVector(3 * i + 1) = plasticSectionDef(1);
      tempVector(3 * i + 2) = plasticSectionDef(2);
    }

    return info.setVector(tempVector);

  } else if (responseID == 100) { // integration points

    double L = crdTransf->getInitialLength();
    double pts[MAX_NUM_SECTIONS];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs[i] = pts[i] * L;
    return info.setVector(locs);

  } else if (responseID == 101) { // integration weights
    double L = crdTransf->getInitialLength();
    double wts[MAX_NUM_SECTIONS];
    beamIntegr->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights[i] = wts[i] * L;
    return info.setVector(weights);

  } else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return info.setID(tags);

  } else if (responseID == 102) { // connected nodes
    Vector tempVector(2);
    tempVector(0) = connectedExternalNodes(0);
    tempVector(1) = connectedExternalNodes(1);
    return info.setVector(tempVector);

  } else if (responseID == 103) { // number of sections
    Vector tempVector(1);
    tempVector(0) = numSections;
    return info.setVector(tempVector);

  }

  else if (responseID == 111 || responseID == 1111) {
    double L = crdTransf->getInitialLength();
    double pts[MAX_NUM_SECTIONS];
    beamIntegr->getSectionLocations(numSections, L, pts);
    // CBDI influence matrix
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID& code  = sections[i]->getType();
      const Vector& e = sections[i]->getSectionDeformation();
      int order       = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz[i] += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay[i] += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(numSections); // along local y
    Vector dispsz(numSections); // along local z
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    beamIntegr->getSectionLocations(numSections, L, pts);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(numSections, 3);
    static Vector vp(6);
    vp = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i] * vp[0]; // linear shape function
      uxb(1) = dispsy[i];
      uxb(2) = dispsz[i];
      if (responseID == 111)
        uxg = crdTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      else
        uxg = crdTransf->getPointLocalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
    }
    return info.setMatrix(disps);
  }

  return -1;
}

int
MixedFrame3d::sendSelf(int commitTag, Channel& theChannel)
{
  // TODO(sendSelf)
  opserr << "Error: MixedFrame3d::sendSelf -- not yet implemented for MixedFrame3d element";
  return -1;
}

int
MixedFrame3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  // TODO(recvSelf)
  opserr << "Error: MixedFrame3d::sendSelf -- not yet implemented for MixedFrame3d element";
  return -1;
}
