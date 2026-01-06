//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef EulerDeltaFrame3d_h
#define EulerDeltaFrame3d_h

#include <array>
#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <Frame/FiniteElement.h>
#include <FrameSection.h>

class Node;
class CrdTransf;
class BeamIntegration;
class Response;

class EulerDeltaFrame3d : public FiniteElement<2, 3, 6>
{
  public:
    EulerDeltaFrame3d(int tag, 
                      std::array<int,2>& nodes,
                      std::vector<FrameSection*> &,
                      BeamIntegration  &, 
                      CrdTransf &,
                      double rho, int mass_flag, bool use_mass);

    EulerDeltaFrame3d();
    ~EulerDeltaFrame3d();

    const char *getClassType() const {
      return "EulerDeltaFrame3d";
    }
    static constexpr const char* class_name = "EulerDeltaFrame3d";


    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information    
    const Vector &getResistingForce();
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);


    // public methods for element output
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    // Parameters
    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int parameterID, Information &);
    int activateParameter(int parameterID);


  protected:
    int setNodes();
    
  private:
    int getIntegral(Field field, State state, double& total);

    constexpr static int
      nen = 2,
      ndf = 6,
      ndm = 3,
      nsr = 4,
      maxNumSections = 20; // TODO: remove
    constexpr static int max_nip = 20;

    int mass_flag;
    bool use_density;

    struct ShapeBasis {
      double eval,
             grad,
             hess;
    };
    ShapeBasis getBasis(int, double);

    double xi[max_nip];
    double wt[max_nip];

    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;
    };

    std::vector<GaussPoint> points;

    int numSections;
    FrameSection **sections;               // vector of Sections
    CrdTransf *theCoordTransf;      // pointer to manifold transformation

    BeamIntegration *beamInt;


    OpenSees::VectorND<6*nen> p;   //
    OpenSees::VectorND<6*nen> q0;  // Fixed end forces in basic system
    OpenSees::VectorND<5>     p0;  // Reactions in basic system

    double density;    // Mass density per unit length

    int parameterID;

    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Mz,
      FrameStress::My,
      FrameStress::T,
    };

};
#endif

