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

  ElasticLinearFrameSection3d();
  ~ElasticLinearFrameSection3d();

  const char *getClassType() const {
       return "ElasticFrameSection3d";
  }
  
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation();
  
  virtual int getIntegral(Field field, State state, double&) const override final;

  const Vector &getStressResultant();
  const Matrix &getSectionTangent();
  const Matrix &getInitialTangent();
  const Matrix &getSectionFlexibility();
  const Matrix &getInitialFlexibility();
  
  FrameSection *getFrameCopy();
  virtual FrameSection* getFrameCopy(const FrameStressLayout& layout);
  const ID &getType();
  int getOrder() const;
  
  int sendSelf(int commitTag, Channel &);
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
  
  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
                                              bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);

 protected:
  
 private:
  struct Tangent {
     OpenSees::MatrixND<3,3> nn,         nw, nv, 
                                 mn, mm, mw, mv, 
                                         ww,
                                             vv;
     void zero() {
            nn.zero();            nw.zero(); nv.zero();
            mn.zero(); mm.zero(); mw.zero(); mv.zero();
                                  ww.zero();
                                             vv.zero();
     }
  } K_pres;

  void getConstants(FrameSectionConstants& consts) const; 

  double E, G;

  constexpr static int nr = 12;

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
};

#endif
