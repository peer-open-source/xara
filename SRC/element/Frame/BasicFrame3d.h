//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma  once
#include <Field.h>
#include <State.h>
#include <Frame/FiniteElement.h>
#include <FrameTransform.h>

#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>
#include <VectorND.h>
#include <FrameLoad.h>
#include <ElementalLoad.h>
#include <vector>
#include <utility>

using namespace OpenSees;

class BasicFrame3d : public FiniteElement<2, 3, 6> {
  constexpr static int ndm = 3;

  public:
    virtual ~BasicFrame3d();

    BasicFrame3d(int tag, int clstag,
                 std::array<int, 2> &nodes, 
                 FrameTransform3d& tran)
      : FiniteElement<2, 3, 6> (tag, clstag, nodes), 
        theCoordTransf(tran.getCopy()),
        wx(0.0), wy(0.0), wz(0.0),
        p_iner(12),
        parameterID(0)
    {
        zeroLoad();

        if (theCoordTransf == nullptr) {
          opserr << "PrismFrame3d::PrismFrame3d -- failed to get copy of coordinate transformation\n";
        }
    }

    BasicFrame3d(int tag, int classtag)
      : FiniteElement<2, 3, 6> (tag, classtag), theCoordTransf(nullptr), 
        p_iner(12), 
        wx(0.0), wy(0.0), wz(0.0),
        parameterID(0)
    {
        zeroLoad();
    }

    // FrameElement
    virtual int getIntegral(Field field, State state, double& total) {
      return -1;
    }

    //
    // For FiniteElement
    //
    virtual int             setNodes();

    //
    // For Element
    //
    virtual int   update();
    virtual const Matrix &getTangentStiff();
    virtual const Matrix &getMass();

    virtual void  zeroLoad();
    virtual int   addLoad(ElementalLoad *theLoad, double loadFactor) final;

    virtual int   addInertiaLoadToUnbalance(const Vector &accel) final;
    virtual const Vector &getResistingForceIncInertia() final;
    virtual const Matrix &getInitialStiff();

    // Sensitivity
    const Matrix & getMassSensitivity(int gradNumber);
    virtual int setParameter(const char **argv, int argc, Parameter &);
    virtual int updateParameter(int param, Information &);
    virtual int activateParameter(int param);


protected:

    // Loads
    double wx, wy, wz;
    OpenSees::VectorND<6>   q0;  // Fixed end forces in basic system
    OpenSees::VectorND<6>   p0;  // Reactions in basic system
                                 // TODO(cmp): change to size 12

  // Implemented by children
  virtual VectorND<6>&   getBasicForce() = 0;
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate) = 0;


  // Supplied to children
  // Reactions of basic system due to element loads
  void addReactionGrad(double *dp0dh, int gradNumber);
  void computeReactions(double *p0);

// to be made private
   FrameTransform3d* theCoordTransf;

   int parameterID;

   std::vector<std::pair<ElementalLoad*,double>> eleLoads;
   std::vector<FrameLoad*> frame_loads;


   Vector p_iner;

private:
   int cMass;
   double rho;

  //  VectorND<12> pg;
   double total_mass,
          twist_mass;

   // TODO: Remove
    int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
    int releasey; // same for y-axis
};

