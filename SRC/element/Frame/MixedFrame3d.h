//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign
//
#ifndef MixedFrame3d_h
#define MixedFrame3d_h
//
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <Frame/FiniteElement.h>
// #include <FrameTransform.h>
#include <FrameSection.h>

class Node;
class Channel;
class Response;
class CrdTransf;
class BeamIntegration;
class FrameSection;

#define ELE_TAG_MixedFrame3d 0 // TODO

class MixedFrame3d : public FiniteElement<2,3,6> {
public:
  // constructors
  MixedFrame3d(int tag, 
               std::array<int, 2>& nodes, 
               int numSections,
               FrameSection** sectionPtrs, 
               BeamIntegration& bi, 
               CrdTransf& coordTransf,
               double density, 
               int damp_flag, 
               int geom_flag = true);
  MixedFrame3d();

  ~MixedFrame3d();

  const char*
  getClassType() const
  {
    return "MixedFrame3d";
  }

  int setNodes(); // (Domain* theDomain);

  // public methods to set the state of the element
  int update();
  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getDamp();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  // MovableObject
  int sendSelf(int cTag, Channel&);
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&);

  // TaggedObject
  void Print(OPS_Stream& s, int flag);

  Response* setResponse(const char** argv, int argc, OPS_Stream& output);
  int getResponse(int responseID, Information& info);


private:

// Constants that define the dimensionality
#define NDM                      3  // dimension of the problem (3d)
#define NND                      6  // number of nodal dof's
#define NEGD                     12 // number of element global dof's
#define NDM_SECTION              3  // number of section dof's without torsion
#define NDM_NATURAL              5  // number of element dof's in the basic system without torsion
#define NDM_NATURAL_WITH_TORSION 6  // number of element dof's in the basic system with torsion
#define MAX_NUM_SECTIONS         10 // maximum number of sections allowed

  constexpr static int nsr = 3;

  enum Geometry { Linear = 1};
  enum Damping  { Rayleigh = 1};
  // Shape Functions
  MatrixND<NDM_SECTION, NDM_NATURAL> getNld_hat(int sec, const Vector& v, double L, int geom_flag);
  Vector getd_hat(int sec, const Vector& v, double L, int geom_flag);
  MatrixND<NDM_SECTION, NDM_NATURAL> getNd1(int sec, const Vector& v, double L, int geom_flag);
  MatrixND<NDM_SECTION, NDM_NATURAL> getNd2(int sec, double P, double L);
  MatrixND<NDM_NATURAL, NDM_NATURAL> getKg(int sec, double P, double L);
  MatrixND<NDM_NATURAL, NDM_NATURAL> getMd(int sec, Vector& dShapeFcn, Vector& dFibers, double L);

  // Interaction With The Sections
  void getSectionTangent(int sec, int type, MatrixND<NDM_SECTION, NDM_SECTION>& Ks, double& GJ);
  void getSectionStress(int sec, Vector& fSection, double& torsion);
  void setSectionDeformation(int sec, Vector& defSection, double& twist);


  BeamIntegration* beamIntegr;        //
  CrdTransf* crdTransf;        // pointer to coordinate transformation object

  int numSections;                    //
  FrameSection** sections;            // array of pointers to sections

  int damp_flag;   // flag for whether or not rayleigh damping is active for this element
  int geom_flag;   // flag for whether or not the internal geometric nonlinearity is active
  double rho;      // mass density per unit length

  int itr; // Count the number of iterations
  int state_flag;

  // Attributes that do NOT change during the analysis
  double L0;
  Matrix* Ki;

  // Element Load Variables
  Matrix* sp;
  double p0[5]; // Reactions in the basic system due to element loads

  // Attributes that change during the analysis
  Vector V;
  Vector qe_pres;
  Vector naturalForce;
  Vector lastNaturalDisp;
  
  MatrixND<NDM_NATURAL, NDM_NATURAL> Hinv;
  MatrixND<NDM_NATURAL, NDM_NATURAL> GMH;
  MatrixND<NDM_NATURAL_WITH_TORSION, NDM_NATURAL_WITH_TORSION> kv; // stiffness matrix in the basic system
  Vector* sr_trial;
  Vector* es_trial;
  MatrixND<NDM_SECTION, NDM_SECTION>* fs_trial;

  // Committed versions
  Vector  committedV;
  Vector  qe_past;
  Vector  commitedNaturalForce;
  Vector  commitedLastNaturalDisp;
  MatrixND<NDM_NATURAL, NDM_NATURAL> commitedHinv;
  MatrixND<NDM_NATURAL, NDM_NATURAL> commitedGMH;
  MatrixND<NDM_NATURAL_WITH_TORSION, NDM_NATURAL_WITH_TORSION>  ke_past;
  Vector* sr_past;
  Vector* es_past;
  MatrixND<NDM_SECTION, NDM_SECTION>* fs_past;

  //
  // Static data
  //
  static Matrix theMatrix;
  static Vector theVector;
  static Matrix transformNaturalCoords;
  // matrix to transform the natural coordinates from what the coordinate transformation uses and what the element uses

  // These variable are always recomputed, so there is no need to store them for each instance of the element
  static Vector* ei_trial;
  static Vector* si_trial;
  static MatrixND<NDM_SECTION, NDM_NATURAL>* nldhat;
  static MatrixND<NDM_SECTION, NDM_NATURAL>* nd1;
  static MatrixND<NDM_SECTION, NDM_NATURAL>* nd2;

  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  double xi[MAX_NUM_SECTIONS];

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Mz,
    FrameStress::My,
//  FrameStress::T,
  };
};

#endif
