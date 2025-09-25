//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
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
#include <BasicFrameTransf.h>

class Channel;
class Information;
class Response;
class FrameSection;
namespace OpenSees {
  class FrameTransformBuilder;
}
using namespace OpenSees;

class PrismFrame3d : public BasicFrame3d, public FiniteElement<2, 3, 6>
{
  public:
    PrismFrame3d(int tag, 
                 std::array<int, 2>& nodes,
                 double A, double E, double G, 
		             double Jx, double Iy, double Iz,
                 FrameTransformBuilder&,
                 double density, int mass_flag,
		             int releasez, int releasey,
                 int geom,
                 int shear_flag);

    PrismFrame3d(int tag,
                 std::array<int,2>& nodes,
                 FrameSection &section, 
                 FrameTransformBuilder&,
                 double density, int mass_flag, bool use_mass,
		             int releasez, int releasey,
                 int geom,
                 int shear_flag);

    ~PrismFrame3d() {
      if (basic_system != nullptr) {
        delete basic_system;
        basic_system = nullptr;
      }
    }

    const char *getClassType() const override {
      return "PrismFrame3d";
    }

    void zeroLoad() override {
      this->BasicFrame3d::zeroLoad();
      this->FiniteElement<2, 3, 6>::zeroLoad();
    }
    
    int addLoad(ElementalLoad *theLoad, double loadFactor) override {
      return this->BasicFrame3d::addLoad(theLoad, loadFactor);
    }
/*
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
*/

    int update() override;
    int commitState() override;
    int revertToLastCommit() override;        
    int revertToStart() override;
    const Matrix &getTangentStiff() final;
    const Matrix &getInitialStiff();
    const Vector &getResistingForce() final;
    const Matrix &getMass() final;

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

    void Print(OPS_Stream &s, int flag) override;

    Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
    int getResponse(int responseID, Information &) final;
 
    // Parameter
    int setParameter(const char **argv, int argc, Parameter &) final;
    int updateParameter(int parameterID, Information &) final;
    int activateParameter(int param)
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
    constexpr static int NBV = 6,
                         NDF = 6;
    struct Param {
      enum {E, G, A, Ay, Az, Iy, Iz, J, HingeY, HingeZ, Rho};
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

    OpenSees::MatrixND<6,6> ke;
    OpenSees::MatrixND<6,6> km;
    OpenSees::MatrixND<6,6> kg;
    OpenSees::VectorND<6>   q;
    BasicFrameTransf3d<NDF> *basic_system;
};


#endif
