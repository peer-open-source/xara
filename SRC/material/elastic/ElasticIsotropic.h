#pragma once
#include "Mate.h"

//         MaterialConstraint
enum class PlaneType {None, Stress, Strain, Frame};

namespace OpenSees {

template <int ndim, PlaneType = PlaneType::None>
class ElasticIsotropic : public Mate<ndim> {
  public:
    using StrainType = MatrixSD<ndim,true>;
  
    static constexpr int ne = StrainType::size;

  ElasticIsotropic(int tag, double E, double nu, double rho);

  virtual Mate<ndim>* getCopy() final;
  virtual const char* getClassType() const final;

  virtual int revertToStart() final;

//virtual MatrixND<ndim,ndim> getStress_matrix() final;
  virtual const MatrixSD<ndim>& getStress() final;
  virtual int setTrialStrain(const MatrixSD<ndim,true>& E) final;
  virtual int setTrialStrain(const MatrixND<ndim,ndim>& F) final;

  virtual MatrixSD<ne> getTangent() final;
  virtual MatrixSD<ne> getInitialTangent() final;

  virtual void Print(OPS_Stream& s, int flag);

private:
  double E, nu, rho;
  MatrixSD<ne> ddsdde;
  MatrixSD<ndim> Smat;
};

}
#include "ElasticIsotropic.tpp"
