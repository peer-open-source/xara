//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef ElasticLinearFrameSection3d_h
#define ElasticLinearFrameSection3d_h

#include <array>
#include <memory>
#include <MatrixND.h>
#include <FrameSection.h>

class Matrix;
class Vector;
class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class ElasticLinearFrameSection3d : public FrameSection
{

 public:
  ElasticLinearFrameSection3d(int tag, 
              double E, 
              double G, 
              const FrameSectionConstants&,
              double mass,
              bool   use_mass
  );

public:
  ~ElasticLinearFrameSection3d();

  const char *getClassType() const final {
    return "ElasticFrameSection3d";
  }
  
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;
  
  int setTrialSectionDeformation(const Vector&) override;
  const Vector &getSectionDeformation();
  
  int getIntegral(Field field, State state, double&) const final;

  const Vector &getStressResultant() override;
  const Matrix &getSectionTangent() override;
  const Matrix &getInitialTangent() override;
  const Matrix &getSectionFlexibility() override;
  const Matrix &getInitialFlexibility() override;
  
  FrameSection *getFrameCopy() final;
  FrameSection* getFrameCopy(const FrameStressLayout&) final;
  const ID &getType() final;
  int getOrder() const final;
  
  int sendSelf(int commitTag, Channel &) final;
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) final;
  
  void Print(OPS_Stream &s, int flag) final;

  int setParameter(const char **argv, int argc, Parameter &);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
                                              bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);

  
 private:

  void getConstants(FrameSectionConstants& consts) const; 


  constexpr static int nr = 12;

  ElasticLinearFrameSection3d(std::shared_ptr<OpenSees::MatrixND<nr,nr>> Ks);

  double E, G;
  Vector  v;
  Matrix  M,         // Generic matrix for returning Ks or Fs (nr x nr)
         *Ksen;      // Tangent sensitivity (nr x nr)

  OpenSees::VectorND<nr> s;
  OpenSees::VectorND<nr> e;                      // section trial deformations
  std::shared_ptr<OpenSees::MatrixND<nr,nr>> Ks;
  Matrix* Fs = nullptr;

  int parameterID;

  static ID layout;
  std::array<double, 2> centroid;


//   struct Tangent {
//      OpenSees::MatrixND<3,3> nn,         nw, nv, 
//                                  mn, mm, mw, mv, 
//                                          ww,
//                                              vv;
//      void zero() {
//             nn.zero();            nw.zero(); nv.zero();
//             mn.zero(); mm.zero(); mw.zero(); mv.zero();
//                                   ww.zero();
//                                              vv.zero();
//      }
//   } K_pres;
};

#endif
