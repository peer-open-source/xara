//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp 2025
//
#ifndef PrismFrame3d02_h
#define PrismFrame3d02_h

#include <array>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Frame/BasicFrame3d.h>
#include <BasicFrameTransf.h>

class Channel;
class Information;
class Response;
class FrameSection;
namespace OpenSees {
  class FrameTransformBuilder;
}
using namespace OpenSees;

class PrismFrame3d02 : public BasicFrame3d, public FiniteElement<2, 3, 6>
{
  public:
    PrismFrame3d02(int tag, 
                 std::array<int, 2>& nodes,
                 double A, double E, double G, 
		             double Jx, double Iy, double Iz,
                 FrameTransformBuilder&,
                 double density, int mass_flag,
		             int releasez, int releasey,
                 int geom,
                 int shear_flag);

    PrismFrame3d02(int tag,
                 std::array<int,2>& nodes,
                 FrameSection &section, 
                 FrameTransformBuilder&,
                 double density, int mass_flag, bool use_mass,
		             int releasez, int releasey,
                 int geom,
                 int shear_flag);
    
    ~PrismFrame3d02() {
      if (basic_system != nullptr) {
        delete basic_system;
        basic_system = nullptr;
      }
    }

    const char *getClassType() const {
      return "PrismFrame3d02";
    }

    virtual void zeroLoad() {
      this->BasicFrame3d::zeroLoad();
      this->FiniteElement<2, 3, 6>::zeroLoad();
    }
    
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor) final {
      return this->BasicFrame3d::addLoad(theLoad, loadFactor);
    }
/*
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
*/
    
    int update();
    int commitState();
    int revertToLastCommit();        
    int revertToStart();
    virtual const Matrix &getTangentStiff() final;
    virtual const Matrix &getInitialStiff();
    virtual const Vector &getResistingForce() final;
    virtual const Matrix &getMass() final;

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
    
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
    virtual int getResponse(int responseID, Information &) final;
 
    // Parameter
    virtual int setParameter(const char **argv, int argc, Parameter &) final;
    virtual int updateParameter(int parameterID, Information &) final;
    virtual int activateParameter(int param)
    {
      parameterID = param; 
      return 0;
    }

    // Sensitivity
    virtual int getResponseSensitivity(int response, int grad, Information&);
    virtual const Vector& getResistingForceSensitivity(int gradNumber) override;

  protected:
    virtual  OpenSees::MatrixND<12,12>& getBasicTangent(State state, int rate) final;
    // For FiniteElement
    virtual int setNodes() final;

  private:
    constexpr static int NEN =  2;
    constexpr static int NBV = 12,
                         NDF =  6;
    struct Param {
      enum {E, G, A, Ay, Az, Iy, Iz, J, HingeY, HingeZ, Rho};
    };

    void formBasicStiffness(OpenSees::MatrixND<NBV,NBV>& kb) const;
    VectorND<NBV> getBasicForceGrad(int grad);

    double E,G;    // elastic properties

    double A;    // cross sectional area
    double Ay;   // shear area along local y axis
    double Az;   // shear area along local z axis

    double Iy;   // moment of inertia about local y axis
    double Iz;   // moment of inertia about local z axis
    double Iyz;  // product of inertia
    double Jx;   // torsion constant

    double rho;  // mass per unit length

    //
    double phiY; // ratio of bending to shear stiffness about local y axis
    double phiZ; // ratio of bending to shear stiffness about local z axis
    double L;    // element length

    int geom_flag = 0; 
    int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
    int releasey; // same for y-axis
    int mass_flag;
    int shear_flag = 0;
    int parameterID;

    double total_mass,
           twist_mass,
           density;

    int section_tag;

    OpenSees::VectorND<12>   q;
    BasicFrameTransf3d<NDF> *basic_system;

    Vector q_send;
};


#endif
