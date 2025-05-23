//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma  once
#include <Field.h>
#include <State.h>
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

class BasicFrame3d {
  constexpr static int ndm = 3;

  public:

    BasicFrame3d()
    : L(0.0),
      wx(0.0), wy(0.0), wz(0.0)
    {
      zeroLoad();
    }

    void setLength(double length) {
      L = length;
    }

    // virtual int getIntegral(Field field, State state, double& total) {
    //   return -1;
    // }

    //
    // For Element
    //
    void  zeroLoad();
    int   addLoad(ElementalLoad *theLoad, double loadFactor);

protected:

  // Loads
  double wx, wy, wz;
  OpenSees::VectorND<6> q0;  // Fixed end forces in basic system
  OpenSees::VectorND<6> p0;  // Reactions in basic system
                             // TODO(cmp): change to size 12

  // Supplied to children
  // Reactions of basic system due to element loads
  void addReactionGrad(double *dp0dh, int gradNumber, double dLdh);
  void computeReactions(double *p0);

   std::vector<std::pair<ElementalLoad*,double>> eleLoads;
   std::vector<FrameLoad*> frame_loads;

private:
  double L;

};

