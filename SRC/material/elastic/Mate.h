//
//
//
#pragma once
#include <Matrix3D.h>
#include <MatrixSD.h>
#include <MatrixND.h>
#include <TaggedObject.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::MatrixSD;

class Response;
class Information;

namespace OpenSees {

template <int ndim>
class Mate : public TaggedObject {
  public:
    using StrainType = MatrixSD<ndim,true>;
  
  static constexpr int ne = StrainType::size;

  Mate(int tag) : TaggedObject(tag) {}

  virtual Mate<ndim>* getCopy() = 0;

  virtual int getClassTag() { return 0;}

  virtual int commitState()
  {
    return 0;
  }

  virtual int revertToStart()
  {
    return 0;
  }

  virtual int revertToLastCommit()
  {
    return 0;
  }

  virtual double getDensity() 
  {
    return 0.0;
  }

  virtual MatrixSD<ne> getTangent() = 0;
  virtual MatrixSD<ne> getInitialTangent() = 0;

  virtual const MatrixSD<ndim>& getStress() = 0;

  virtual int setTrialStrain(const MatrixSD<ndim,true>&) = 0;

  virtual int setTrialStrain(const MatrixND<ndim,ndim>&) = 0;

#if 1
  virtual Response *setResponse (const char **argv, int argc, OPS_Stream &s) {return 0;}
  virtual int getResponse (int responseID, Information &matInformation) { return -1;}
#endif

};
}

