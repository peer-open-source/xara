//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for 
// FrameFiberSection3d.h. FrameFiberSection3d provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: cmp
// Created: 2024
//
#ifndef FrameFiberSection3d_h
#define FrameFiberSection3d_h

#define SEC_TAG_FrameFiberSection3d 0

#include <FrameSection.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <memory>

class Response;
class UniaxialMaterial;

class FrameFiberSection3d : public FrameSection
{
  public:
    FrameFiberSection3d(); 
    FrameFiberSection3d(int tag, int numFibers, UniaxialMaterial &torsion, 
                        bool compCentroid,
                        double mass, bool use_mass);
    ~FrameFiberSection3d();

    const char *getClassType() const {
      return "FrameFiberSection3d";
    }

    int   setTrialSectionDeformation(const Vector &deforms);
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const override final;
    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();

    int   commitState();
    int   revertToLastCommit();    
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder () const; //  {return 4;};
 
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(UniaxialMaterial &theMat, const double area, const double y, const double z);
//  int setField(const char**, int, double);

    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector & getSectionDeformationSensitivity(int gradIndex);

    double getEnergy() const;


  protected:
    constexpr static int nsr = 4;
    constexpr static int nwm = 3;

  private:
    struct FiberData {
      double y;
      double z;
      double area,
             warp[nwm][3];
    };


    int numFibers, sizeFibers;         // number of fibers in the section
    UniaxialMaterial **theMaterials;   // array of pointers to materials
    std::shared_ptr<double[]> matData; // data for the materials [yloc, zloc, and area]

    OpenSees::MatrixND<nsr,nsr> ks;
    Matrix K_wrap;

    double QzBar, QyBar, Abar;
    double yBar;                       // Section centroid
    double zBar;
    bool computeCentroid;

    static ID code;

    OpenSees::VectorND<nsr> es, sr;
    Vector  e;         // trial section deformations 
    Vector  s;         // section resisting forces  (axial force, bending moment)

    UniaxialMaterial *theTorsion;
    void *pool;        // thread pool
};

#endif
