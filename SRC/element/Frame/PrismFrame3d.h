//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp 2024
//
#ifndef PrismFrame3d_h
#define PrismFrame3d_h

#include <array>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Frame/BasicFrame3d.h>

class Channel;
class Information;
class FrameTransform3d;
class Response;
class FrameSection;

class PrismFrame3d : public BasicFrame3d
{
  public:
    PrismFrame3d(int tag, 
                 std::array<int, 2>& nodes,
                 double A, double E, double G, 
		             double Jx, double Iy, double Iz,
                 FrameTransform3d &theTransf,
                 double density, int mass_flag,
		             int releasez, int releasey,
                 int geom);

    PrismFrame3d(int tag,
                 std::array<int,2>& nodes,
                 FrameSection &section, 
		             FrameTransform3d &theTransf,
                 double density, int mass_flag, bool use_mass,
		             int releasez, int releasey,
                 int geom,
                 int shear_flag);

    PrismFrame3d();

//  ~PrismFrame3d();

    const char *getClassType() const {
      return "PrismFrame3d";
    }

    // void zeroLoad();	
/*
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
*/
    
    int update();
    int commitState();
    int revertToLastCommit();        
    int revertToStart();
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
    // For BasicFrame3d
    virtual  OpenSees::VectorND<6>&   getBasicForce() final;
    virtual  OpenSees::MatrixND<6,6>& getBasicTangent(State state, int rate) final;
    // For FiniteElement
    virtual int setNodes() final;

  private:
    constexpr static int NEN = 2;
    constexpr static int NBV = 6;
    struct Param {
      enum {E, G, A, Ay, Az, Iy, Iz, J};
    };

    void formBasicStiffness(OpenSees::MatrixND<6,6>& kb) const;
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

    OpenSees::MatrixND<6,6> ke;  // TODO: make shared_ptr?
    OpenSees::MatrixND<6,6> km;
    OpenSees::MatrixND<6,6> kg;
    OpenSees::VectorND<6>   q;
    // Matrix K;
    // Vector P;
};

#endif
