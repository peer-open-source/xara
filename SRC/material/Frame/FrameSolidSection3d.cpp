//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class implementation of FrameSolidSection3d.
// FrameSolidSection3d provides the abstraction of a 3D beam section discretized by fibers.
//
// Written: cmp
// Created: Spring 2025
//
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <array>
#include <Channel.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <classTags.h>
#include "FrameSolidSection3d.h"
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <NDMaterial.h>
#include <Parameter.h>

#define SEC_TAG_FrameSolidSection3d 0

using namespace OpenSees;

ID FrameSolidSection3d::code(nsr);

FrameSolidSection3d::FrameSolidSection3d(int tag, int num): 
    FrameSection(tag, SEC_TAG_FrameSolidSection3d),
    s(), e(),
    e_wrap(e), s_wrap(s),
    shear_align{},
    centroid{},
    nubar(0.0),
    parameterID(0), dedh(nsr),
    fibers(new std::vector<FiberData>),
    K_init(new Tangent),
    fiber_state(FiberState::Clean)
{
  code(inx) = FrameStress::N;
  code(iny) = FrameStress::Vy;
  code(inz) = FrameStress::Vz;
  code(imx) = FrameStress::T;
  code(imy) = FrameStress::My;
  code(imz) = FrameStress::Mz;
  code(iwx) = FrameStress::Bimoment;
  code(iwy) = FrameStress::By;
  code(iwz) = FrameStress::Bz;
  code(ivx) = FrameStress::Bishear;
  code(ivy) = FrameStress::Qy;
  code(ivz) = FrameStress::Qz;

  wagner = getenv("Wagner") != nullptr;
}


// for recvSelf
FrameSolidSection3d::FrameSolidSection3d():
  FrameSection(0, SEC_TAG_FrameSolidSection3d),
  s(), e(),
  e_wrap(e), s_wrap(s),
  shear_align{},
  centroid{},
  nubar(0.0),
  parameterID(0), dedh(nsr),
  fibers(new std::vector<FiberData>),
  fiber_state(FiberState::Clean)
{
  code(inx) = FrameStress::N;
  code(iny) = FrameStress::Vy;
  code(inz) = FrameStress::Vz;
  code(imx) = FrameStress::T;
  code(imy) = FrameStress::My;
  code(imz) = FrameStress::Mz;
  code(iwx) = FrameStress::Bimoment;
  code(iwy) = FrameStress::By;
  code(iwz) = FrameStress::Bz;
  code(ivx) = FrameStress::Bishear;
  code(ivy) = FrameStress::Qy;
  code(ivz) = FrameStress::Qz;
  wagner = getenv("Wagner") != nullptr;
}

// Used in getCopy to create an element instance from a reference instance
FrameSolidSection3d::FrameSolidSection3d(const FrameSolidSection3d &other)
  : FrameSection(other.getTag(), other.getClassTag()),
    K_pres(other.K_pres),
    K_init(other.K_init),
    fibers(other.fibers),
    s(),
    e(),
    s_wrap(s),
    e_wrap(e),
    shear_align(other.shear_align),
    centroid(other.centroid),
    wagner(other.wagner),
    fiber_state(FiberState::Clean),
    parameterID(0)
{
  materials.reserve(other.materials.size());
  for (int i = 0; i < other.materials.size(); i++)
    materials.push_back(other.materials[i]->getCopy("BeamFiber"));
    // materials[i] = other.materials[i]->getCopy("BeamFiber");

  this->revertToStart();
}


int
FrameSolidSection3d::getIntegral(Field field, State state, double& value) const
{
  value = 0.0;

  const int nf = fibers->size();
  switch (field) {
    case Field::Unit:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        value += A;
      }
      return 0;

    case Field::Density:
      // First check if density has been specified for the section
      if (this->FrameSection::getIntegral(field, state, value) == 0) 
        return 0;

      for (int i=0; i<nf; i++) {
        double density;
        const double A  = (*fibers)[i].area;
        if (materials[i]->getRho() != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;

    case Field::UnitY: // TODO: Centroid
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y; // - yBar;
        value += A*y;
      }
      return 0;

    case Field::UnitZ: // TODO: Centroid
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].z; // - yBar;
        value += A*z;
      }
      return 0;

    case Field::UnitYY:
    case Field::UnitCentroidYY:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y
                        ; //- yBar*(Field::UnitCentroidYY == field);
        value += A*y*y;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].z
                        ; //- zBar*(Field::UnitCentroidZZ == field);
        value += A*z*z;
      }

    default:
      return -1;
  }
  return -1;
}


int
FrameSolidSection3d::addFiber(NDMaterial& theMat, 
                              double Area, 
                              double yLoc, 
                              double zLoc)
{
  std::array<std::array<double,3>,3> warp{0};
  FiberData fiber {yLoc, zLoc, Area, warp, {0.0, yLoc, zLoc}};
  fibers->emplace_back(fiber);

  materials.emplace_back(theMat.getCopy("BeamFiber"));

  if (materials[materials.size()-1] == nullptr)
    return -1;

  fiber_state = FiberState::Dirty;
  return materials.size()-1;
}



FrameSolidSection3d::~FrameSolidSection3d()
{
  for (auto material : materials) {
    if (material != nullptr)
      delete material;
  }
}

int
FrameSolidSection3d::setTrialSectionDeformation(const Vector &e_trial)
{
  e = e_trial;
  s.zero();
  return stateDetermination(K_pres, &s, &e, CurrentTangent);
}


FrameSection*
FrameSolidSection3d::getFrameCopy()
{

  double area = 0.0;
  if (fiber_state == FiberState::Dirty) {
    nubar = 0.0;
    const int nf = fibers->size();
    centroid.zero();
    for (int i = 0; i < nf; i++) {
      NDMaterial &theMat = *materials[i];
      auto & fiber = (*fibers)[i];
      const Matrix & tangent = theMat.getInitialTangent();
      double nu = tangent(0,0)/(2.0*tangent(1,1)) - 1.0;
      area  += fiber.area;
      nubar += nu*fiber.area;
      centroid[1] += fiber.r[1]*fiber.area;
      centroid[2] += fiber.r[2]*fiber.area;
    }
    centroid /= area;
    nubar /= area;
  }
  if (!(getenv("Align"))) {
    fiber_state = FiberState::Clean;
  }

  if (fiber_state == FiberState::Dirty) {
    const VectorND<2> rc {
      centroid[1],
      centroid[2]
    };
    const int nf = fibers->size();
    double do_poisson = shear_align.norm() != 0.0;
    MatrixND<2,2> H{}, HF{};
    for (int i = 0; i < nf; i++) {
      NDMaterial &theMat = *materials[i];
      auto & fiber = (*fibers)[i];
      auto & w = fiber.warp;
      const VectorND<2> r {fiber.r[1], fiber.r[2]};
      const Matrix & tangent = theMat.getInitialTangent();
      double nu = tangent(0,0)/(2.0*tangent(1,1)) - 1.0;
      nu *= do_poisson;
      MatrixND<2,2> Bb{};
      // -(dev ror)Xi
      Bb.addTensorProduct(r, r, -nu);
      Bb.addDiagonal(0.5*r.dot(r)*nu);
      // + roc 
      Bb.addTensorProduct(r, rc, nu);
      for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          Bb(j,k) += w[1+k][1+j];
        }
      }
      HF.addMatrix(Bb, fiber.area*tangent(1,1));
      // HF.addTensorProduct(r-rc, r-rc, fiber.area*tangent(0,0));
      H.addMatrixTransposeProduct(1.0, Bb, Bb, fiber.area*tangent(1,1));
    }
    Cholesky<2,false> chol(H);
    MatrixND<2,2> Hi{};
    chol.invert(Hi);
    MatrixND<2,2> X = Hi*HF.transpose();

    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
        shear_align(i+1,j+1) = X(i,j)*nubar;

    for (int k=0; k<nf; k++) {
      FiberData & fiber = (*fibers)[k];
      auto & w = fiber.warp;
      MatrixND<2,2> U {{w[1][1], w[1][2],
                        w[2][1], w[2][2]}};
      MatrixND<2,2> UX = U*X;
      VectorND<2>   xu = X^(VectorND<2>{w[1][0], w[2][0]});
      for (int i=0; i<2; i++) {
        fiber.warp[i+1][0] = xu(i) - fiber.r[i+1];
        for (int j=0; j<2; j++)
          fiber.warp[i+1][j+1] = UX(j,i) - double(i==j);
      }
    }

    fiber_state = FiberState::Clean;
  }
  return new FrameSolidSection3d(*this);
}

int
FrameSolidSection3d::stateDetermination(Tangent& K, VectorND<nsr>* s_trial, const VectorND<nsr> * const e_trial, int tangentFlag)
{

  const bool conj_poisson = getenv("ConjPoisson") != nullptr;

  const Vector3D gamma {
    e_trial? (*e_trial)(inx) : 0.0,
    e_trial? (*e_trial)(iny) : 0.0,
    e_trial? (*e_trial)(inz) : 0.0 
  };
  const Vector3D kappa {
    e_trial? (*e_trial)(imx) : 0.0,
    e_trial? (*e_trial)(imy) : 0.0,
    e_trial? (*e_trial)(imz) : 0.0 
  };
  const Vector3D dalpha {
    e_trial? (*e_trial)(iwx) : 0.0,
    e_trial? (*e_trial)(iwy) : 0.0,
    e_trial? (*e_trial)(iwz) : 0.0 
  };
  const Vector3D alpha {
    e_trial? (*e_trial)(ivx) : 0.0,
    e_trial? (*e_trial)(ivy) : 0.0,
    e_trial? (*e_trial)(ivz) : 0.0 
  };

  if (s_trial != nullptr)
    s_trial->zero();

  K.zero();

  // integrate over fibers
  int res = 0;
  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {

    NDMaterial &theMat = *materials[i];
    auto & fiber = (*fibers)[i];
    auto & w = fiber.warp;
    const Vector3D& r = fiber.r;
    double tr2 = 0;

    Matrix3D Wb{};
    // -(dev ror)Xi
    Wb.addTensorProduct(r, shear_align^r, -1.0);
    Wb.addMatrix(shear_align, 0.5*r.dot(r));
    // + roc 
    Wb.addTensorProduct(r, shear_align^centroid, 1.0);
    for (int j=0; j<3; j++) Wb(0,j) = (Wb(j,0) = 0.0); // No axial Poisson effect

    if (e_trial != nullptr) {
      // Form material strain
      Vector3D eps = gamma + kappa.cross(r);
      for (int k=0; k<nwm; k++) {
        eps[0] += w[k][0]*dalpha[k];
        for (int j=1; j<3; j++)
          eps[j] += w[k][j]*alpha[k];
      }
      // Poisson effect
      eps += Wb*gamma;
      // if (eps.norm() != 0.0) {
      //   Vector3D eps_gamma = gamma + Wb*gamma;
      //   Vector3D eps_kappa = kappa.cross(r);

      //   for (int k=1; k<nwm; k++) {
      //     for (int j=1; j<3; j++)
      //       eps_gamma[j] += w[k][j]*alpha[k];
      //   }
      //   for (int j=1; j<3; j++)
      //     eps_kappa[j] += w[0][j]*alpha[0];
      //   opserr << ", eps_gamma = " << Vector(eps_gamma)
      //          << ", eps_kappa = " << Vector(eps_kappa) << "\n";
      //   exit(1);
      // }

      if (wagner) {
        tr2 = r.dot(r)*kappa[0];
        eps[0] += 0.5*kappa[0]*tr2;
      }
      res += theMat.setTrialStrain(eps);
      if (res < 0)
        return res;
    }


    const Matrix &tangent = tangentFlag==CurrentTangent
                          ? theMat.getTangent()
                          : theMat.getInitialTangent();

    Matrix3D C{};
    C.addMatrix(tangent, fiber.area);

    // NOTE: Matrix 3D is column major so these are transposed.
    const Matrix3D iow {{
      w[0][0],     0.0,     0.0,
      w[1][0],     0.0,     0.0,
      w[2][0],     0.0,     0.0
    }};

    const Matrix3D iodw {{
      0.0, w[0][1], w[0][2],
      0.0, w[1][1], w[1][2],
      0.0, w[2][1], w[2][2]
    }};

    const Matrix3D CWb   = C*Wb;
    const Matrix3D Ciow  = C*iow;
    const Matrix3D Ciodw = C*iodw;
    Matrix3D Crx{};
    Crx.addMatrixSpinProduct( C, r,  -1.0);

    {
      K.nn.addMatrix(C,    1.0);
      K.nn.addMatrix(CWb,  1.0);
      K.nw.addMatrix(Ciow,  1.0);
      K.nm.addMatrix(Crx,  1.0);
      K.nv.addMatrix(Ciodw, 1.0);
      if (conj_poisson) {
        K.nv.addMatrixTransposeProduct(1.0, Wb, Ciodw, 1.0);
        K.nm.addMatrixTransposeProduct(1.0, Wb, Crx, 1.0);
        K.nn.addMatrixTransposeProduct(1.0, Wb, C+CWb, 1.0);
        K.nw.addMatrixTransposeProduct(1.0, Wb, Ciow, 1.0);
      }
    }

    {
      K.mn.addSpinMatrixProduct(r, C,     1.0);
      K.mn.addSpinMatrixProduct(r, CWb,   1.0);
      K.mm.addSpinMatrixProduct(r, Crx,   1.0);
      K.mw.addSpinMatrixProduct(r, Ciow,  1.0);
      K.mv.addSpinMatrixProduct(r, Ciodw, 1.0);
    }
    {
      K.wn.addTranspose(Ciow, 1.0);
      K.wn.addMatrixTransposeProduct(1.0, iow,  CWb,  1.0);
      K.wm.addMatrixTransposeProduct(1.0, iow,  Crx,  1.0);
      K.ww.addMatrixTransposeProduct(1.0, iow,  Ciow, 1.0);
      K.wv.addMatrixTransposeProduct(1.0, iow,  Ciodw, 1.0);
    }
    {
      K.vn.addTranspose(Ciodw, 1.0);
      K.vn.addMatrixTransposeProduct(1.0, iodw,  CWb,   1.0);
      K.vm.addMatrixTransposeProduct(1.0, iodw,  Crx,   1.0);
      K.vw.addMatrixTransposeProduct(1.0, iodw,  Ciow,  1.0);
      K.vv.addMatrixTransposeProduct(1.0, iodw,  Ciodw, 1.0);
    }

    //
    //
    //
    const Vector &stress  = theMat.getStress();

    if (wagner && (e_trial != nullptr)) {
      constexpr Matrix3D ioi {{ 1, 0, 0 ,
                                0, 0, 0 ,
                                0, 0, 0 }};
      Matrix3D ioiC = ioi*C;
      K.mn.addMatrix(ioiC, tr2);
      K.nm.addTranspose(ioiC, tr2);
      //
      K.mm.addSpinMatrixProduct(r, ioiC.transpose(), tr2);
      K.mm.addMatrixSpinProduct(ioiC, r, -tr2);
      K.mm.addMatrixProduct(ioiC, ioi, tr2*tr2);

      // Geometric part,  equivalent to Kmm.addMatrix(ioi, r2*stress(0));
      if (kappa[0] != 0) [[likely]]
        K.mm(0,0) += (tr2/kappa[0])*stress(0)*fiber.area;

      K.mw.addMatrixProduct(ioiC, iow,  tr2);
      K.mv.addMatrixProduct(ioiC, iodw, tr2);
    }

    if (s_trial != nullptr) {
      const double y = fiber.r[1];
      const double z = fiber.r[2];
      const double sig0 = stress(0)*fiber.area;
      const double sig1 = stress(1)*fiber.area;
      const double sig2 = stress(2)*fiber.area;
      // n += s da
      (*s_trial)(inx) +=    sig0;
      (*s_trial)(iny) +=    sig1;
      (*s_trial)(inz) +=    sig2;
      // m += r.cross(s) da
      (*s_trial)(imx) += -z*sig1 + y*sig2;
      (*s_trial)(imy) +=  z*sig0;
      (*s_trial)(imz) += -y*sig0;
      for (int j=0; j<nwm; j++) {
        // w += w[j]*s da
        (*s_trial)(iwx+j) += w[j][0]*sig0;
        // v += dw[j]*s da
        (*s_trial)(ivx+j) += w[j][1]*sig1;
        (*s_trial)(ivx+j) += w[j][2]*sig2;
      }

      if (wagner && e_trial != nullptr)
        (*s_trial)(imx) += tr2*sig0;

      if (conj_poisson) {
        // Poisson effect
        (*s_trial)(iny) += Wb(1,1)*sig1 + Wb(2,1)*sig2;
        (*s_trial)(inz) += Wb(1,2)*sig1 + Wb(2,2)*sig2;
      }
    }
  }
  return res;
}


const Vector&
FrameSolidSection3d::getSectionDeformation()
{
  return e_wrap;
}


const Matrix&
FrameSolidSection3d::getSectionTangent()
{
#ifndef SEES_SECTION_THREADS
  static MatrixND<nsr,nsr> K;
  static Matrix K_wrap(K);
  // K_wrap.setData(K);
#endif

  K.zero();
  K.assemble(K_pres.nn, 0, 0, 1.0);
  K.assemble(K_pres.nm, 0, 3, 1.0);
  K.assemble(K_pres.nw, 0, 6, 1.0);
  K.assemble(K_pres.nv, 0, 9, 1.0);

  K.assemble(K_pres.mn, 3, 0, 1.0);
  K.assemble(K_pres.mm, 3, 3, 1.0);
  K.assemble(K_pres.mw, 3, 6, 1.0);
  K.assemble(K_pres.mv, 3, 9, 1.0);

  K.assemble(K_pres.wn, 6, 0, 1.0);
  K.assemble(K_pres.wm, 6, 3, 1.0);
  K.assemble(K_pres.ww, 6, 6, 1.0);
  K.assemble(K_pres.wv, 6, 9, 1.0);

  //
  K.assemble(K_pres.vn, 9, 0, 1.0);
  K.assemble(K_pres.vm, 9, 3, 1.0);
  K.assemble(K_pres.vw, 9, 6, 1.0);
  K.assemble(K_pres.vv, 9, 9, 1.0);

  // K.assembleTranspose(K_pres.mn, 0, 3, 1.0); // nm
  // K.assembleTranspose(K_pres.nw, 6, 0, 1.0);
  // K.assembleTranspose(K_pres.mw, 6, 3, 1.0);
  // // K.assembleTranspose(K_pres.nv, 9, 0, 1.0);
  // K.assembleTranspose(K_pres.mv, 9, 3, 1.0); // vm
  // K.assembleTranspose(K_pres.wv, 9, 6, 1.0); // vw
  return K_wrap;
}


const Matrix&
FrameSolidSection3d::getInitialTangent()
{
  static MatrixND<nsr,nsr> ksi;
  static Matrix wrap(ksi); //, nsr, nsr);

  ksi.zero();
  this->stateDetermination(K_pres, nullptr, nullptr, InitialTangent);

  ksi.zero();
  ksi.assemble(K_pres.nn, 0, 0, 1.0);
  ksi.assemble(K_pres.mn, 3, 0, 1.0);
  ksi.assemble(K_pres.nw, 0, 6, 1.0);
  ksi.assemble(K_pres.nv, 0, 9, 1.0);
  ksi.assemble(K_pres.mn, 3, 0, 1.0);
  ksi.assemble(K_pres.mm, 3, 3, 1.0);
  ksi.assemble(K_pres.mw, 3, 6, 1.0);
  ksi.assemble(K_pres.mv, 3, 9, 1.0);
  ksi.assemble(K_pres.ww, 6, 6, 1.0);
  ksi.assemble(K_pres.vv, 9, 9, 1.0);
  ksi.assembleTranspose(K_pres.nw, 6, 0, 1.0);
  ksi.assembleTranspose(K_pres.nv, 9, 0, 1.0);
  ksi.assembleTranspose(K_pres.mn, 0, 3, 1.0);
  ksi.assembleTranspose(K_pres.mw, 6, 3, 1.0);
  ksi.assembleTranspose(K_pres.mv, 9, 3, 1.0);
  return wrap;
}



const Vector&
FrameSolidSection3d::getStressResultant()
{
  return s_wrap;
}



const ID&
FrameSolidSection3d::getType()
{
  return code;
}


int
FrameSolidSection3d::getOrder() const
{
  return nsr;
}


int
FrameSolidSection3d::commitState()
{
  int err = 0;

  for (auto& material: materials)
    err += material->commitState();

  return err;
}

int
FrameSolidSection3d::revertToLastCommit()
{
  int err = 0;
  for (auto& material : materials)
    err += material->revertToLastCommit();

  // TODO: we may need to recompute e to be perfectly consistent
  // when Wagner is enabled.
  s.zero();
  err += this->stateDetermination(K_pres, &s, nullptr, CurrentTangent);

  return err;
}


int
FrameSolidSection3d::revertToStart()
{
  int err = 0;
  for (auto& material: materials)
    err += material->revertToStart();

  s.zero();
  e.zero();
  err += this->stateDetermination(K_pres, &s, nullptr, CurrentTangent);

  return err;
}


int
FrameSolidSection3d::sendSelf(int commitTag, Channel &)
{
  return -1;
}

int
FrameSolidSection3d::recvSelf(int , Channel &,
                              FEM_ObjectBroker &)
{
  return -1;
}


Response*
FrameSolidSection3d::setResponse(const char **argv, int argc,
                                 OPS_Stream &output)
{
  Response *theResponse = nullptr;

  if (argc > 2 && strcmp(argv[0], "fiber") == 0) {

    
    int key = fibers->size();
    int passarg = 2;
    
    if (argc <= 3) { // fiber number was input directly
      key = atoi(argv[1]);
    }

    else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag    = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      Vector3D r_given{{0, yCoord, zCoord}};
      double closestDist = 0;

      // Find first fiber with specified material tag
      const int nf = fibers->size();
      int j;
      for (j = 0; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          Vector3D dr = fiber.r - r_given;
          closestDist = dr.dot(dr);
          key = j;
          break;
        }
      }

      // Search the remaining fibers
      double distance;
      for ( ; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          Vector3D dr = fiber.r - r_given;
          distance = dr.dot(dr);
          if (distance < closestDist) {
            closestDist = distance;
            key = j;
          }
        }
      }
      passarg = 4;
    }
    
    else {
      // fiber near-to coordinate specified
      Vector3D r_given{
        0.0, atof(argv[1]), atof(argv[2])
      };
      Vector3D dr = (*fibers)[0].r - r_given;
      double closestDist = dr.dot(dr);
      key = 0;
      double distance;

      const int nf = fibers->size();
      for (int j = 1; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        Vector3D dr = fiber.r - r_given;
        distance = dr.dot(dr);
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }

    if (key < fibers->size() && key >= 0) {
      output.tag("FiberOutput");
      output.attr("y",    (*fibers)[key].r[1]);
      output.attr("z",    (*fibers)[key].r[2]);
      output.attr("area", (*fibers)[key].area);

      theResponse = materials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  }

  if (theResponse == nullptr)
    return FrameSection::setResponse(argv, argc, output);

  return theResponse;
}


int 
FrameSolidSection3d::getResponse(int responseID, Information &sectInfo)
{
  return FrameSection::getResponse(responseID, sectInfo);
}

int
FrameSolidSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strcmp(argv[0], "warp") == 0) {
    // ... warp $fiberID $warpField
    if (argc < 3) {
      opserr << "FrameSolidSection3d::setParameter - fiberID is required\n";
      return -1;
    }
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= fibers->size()) {
      opserr << "FrameSolidSection3d::setParameter - fiberID " << fiberID << " out of range\n";
      return -1;
    }

    int field = atoi(argv[2]);

    return param.addObject(Param::FiberFieldBase+fiberID*100+field, this);
  }

  if (strcmp(argv[0], "shift_shear") == 0) {
    // ... shift_shear i j
    if (argc < 3) {
      opserr << "FrameSolidSection3d::setParameter - i, j, value are required\n";
      return -1;
    }
    int i = atoi(argv[1]);
    int j = atoi(argv[2]);
    if ((i == 1) && (j == 1))
      return param.addObject(Param::ShearAlignYY, this);
    else if ((i == 2) && (j == 2))
      return param.addObject(Param::ShearAlignZZ, this);
    else if ((i == 1) && (j == 2))
      return param.addObject(Param::ShearAlignYZ, this);
    else if ((i == 2) && (j == 1))
      return param.addObject(Param::ShearAlignZY, this);
    else {
      opserr << "FrameSolidSection3d::setParameter - invalid i, j: " << i << ", " << j << "\n";
      return -1;
    }
  }

  // Check if the parameter belongs to the material
  if (strstr(argv[0], "material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    for (auto& material: materials)
      if (materialTag == material->getTag()) {
        int ok = material->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
      }
    return result;
  }

  int ok = 0; 
  for (auto& material: materials) {
    ok = material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  return result;
}

int
FrameSolidSection3d::updateParameter(int paramID, Information &info)
{

  if (paramID == Param::ShearAlignYY) {
    shear_align(1,1) = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShearAlignZZ) {
    shear_align(2,2) = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShearAlignYZ) {
    shear_align(1,2) = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShearAlignZY) {
    shear_align(2,1) = info.theDouble;
    return 0;
  }

  if (paramID >= Param::FiberFieldBase) {
    int fiberID = (paramID - Param::FiberFieldBase) / 100;
    int field   = (paramID - Param::FiberFieldBase) % 100;
  
    if (fiberID >= fibers->size()) 
      return -1;

    switch (field) {
      case Param::FiberArea:
        (*fibers)[fiberID].area = info.theDouble;
        break;
      case Param::FiberY:
        (*fibers)[fiberID].r[1] = info.theDouble;
        break;
      case Param::FiberZ:
        (*fibers)[fiberID].r[2] = info.theDouble;
        break;
      case Param::FiberWarpX:
        (*fibers)[fiberID].warp[0][0] = info.theDouble;
        break;
      case Param::FiberWarpXY:
        (*fibers)[fiberID].warp[0][1] = info.theDouble;
        break;
      case Param::FiberWarpXZ:
        (*fibers)[fiberID].warp[0][2] = info.theDouble;
        break;
      //
      case Param::FiberWarpY:
        (*fibers)[fiberID].warp[1][0] = info.theDouble;
        break;
      case Param::FiberWarpYY:
        (*fibers)[fiberID].warp[1][1] = info.theDouble;
        break;
      case Param::FiberWarpYZ:
        (*fibers)[fiberID].warp[1][2] = info.theDouble;
        break;

      case Param::FiberWarpZ:
        (*fibers)[fiberID].warp[2][0] = info.theDouble;
        break;
      case Param::FiberWarpZY:
        (*fibers)[fiberID].warp[2][1] = info.theDouble;
        break;
      case Param::FiberWarpZZ:
        (*fibers)[fiberID].warp[2][2] = info.theDouble;
        break;

      default:
        return -1;
    }
    return 0;
  }

  return -1;
}

int
FrameSolidSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector &
FrameSolidSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  return dedh;
}

const Vector &
FrameSolidSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(nsr);
  
  ds.Zero();
  
  static Vector stress(3);
  static Vector dsigdh(3);
  static Vector sig_dAdh(3);
  static Matrix tangent(3,3);

  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];
  const int nf = fibers->size();

  for (int i = 0; i < nf; i++) {
    const double y = (*fibers)[i].r[1]; // - yBar;
    const double z = (*fibers)[i].r[2]; // - zBar;
    const double A = (*fibers)[i].area;
    
    dsigdh = materials[i]->getStressSensitivity(gradIndex,true);

    ds[0] += dsigdh(0)*A;
    ds[1] += -y*dsigdh(0)*A;
    ds[2] +=  z*dsigdh(0)*A;
    ds[3] +=  dsigdh(1)*A;
    ds[4] +=  dsigdh(2)*A;
    ds[5] += (-z*dsigdh(1)+y*dsigdh(2))*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0 || parameterID == 1)
      stress = materials[i]->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0 || parameterID == 1)
      tangent = materials[i]->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh(0) = stress(0)*areaDeriv[i];
      sig_dAdh(1) = stress(1)*areaDeriv[i];
      sig_dAdh(2) = stress(2)*areaDeriv[i];
      
      ds[0] += sig_dAdh(0);
      ds[1] += -y*sig_dAdh(0);
      ds[2] +=  z*sig_dAdh(0);
      ds[3] +=    sig_dAdh(1);
      ds[4] +=    sig_dAdh(2);
      ds[5] += -z*sig_dAdh(1)+y*sig_dAdh(2);
    }

    if (dydh[i] != 0.0) {
      ds(1) += -dydh[i] * (stress(0)*A);
      ds(5) +=  dydh[i] * (stress(2)*A);
    }

    if (dzdh[i] != 0.0) {
      ds[2] +=  dzdh[i] * (stress(0)*A);
      ds[5] += -dzdh[i] * (stress(1)*A);
    }

    if (parameterID == 1) {
      ds[3] += (stress(1)*A);
      ds[4] += (stress(2)*A);
    }

    MatrixND<3,6> as;
    as(0,0) =  1;
    as(1,3) =  1;
    as(2,4) =  1;
    as(0,1) = -y;
    as(0,2) =  z;
    as(1,5) = -z;
    as(2,5) =  y;
    
    static Matrix dasdh(3,nsr);
    dasdh(0,1) = -dydh[i];
    dasdh(0,2) = dzdh[i];
    dasdh(1,3) = 0;
    dasdh(2,4) = 0;
    dasdh(1,5) = -dzdh[i];
    dasdh(2,5) = dydh[i];
#if 0 // removed to eliminate implicit cast from MatrixND to Matrix
    static Matrix tmpMatrix(nsr,nsr);
    tmpMatrix.addMatrixTripleProduct(0.0, as, tangent, dasdh, 1.0);
    
    ds.addMatrixVector(1.0, tmpMatrix, e, A);
#endif
  }

  return ds;
}

const Matrix &
FrameSolidSection3d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(6,6);
  
  dksdh.Zero();
  return dksdh;
}


int
FrameSolidSection3d::commitSensitivity(const Vector& defSens,
                                       int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);
  double d4 = defSens(4);
  double d5 = defSens(5);

  dedh = defSens;

  static double dydh[10000];
  static double dzdh[10000];
  const int nf = fibers->size();
  
  { // TODO
    for (int i = 0; i < nf; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
    }
  }

  static Vector depsdh(3);


  for (int i = 0; i < nf; i++) {
    auto& fiber = (*fibers)[i];
    const double y  = fiber.r[1]; // - yBar;
    const double z  = fiber.r[2]; // - zBar;

    // determine material strain sensitivity
    depsdh[0] = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);
    depsdh[1] = d3 - z*d5 + e(3) - dzdh[i]*e(5);
    depsdh[2] = d4 + y*d5 + e(4) + dydh[i]*e(5);

    materials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}


void
FrameSolidSection3d::Print(OPS_Stream &s, int flag)
{
  const int nf = fibers->size();
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    double mass;
    if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0)
      s << "\"mass\": " << mass << ", ";


    s << "\"shear_align\": [";
    for (int i = 1; i < 3; i++) {
      s << "[";
      for (int j = 1; j < 3; j++) {
        s << shear_align(i,j);
        if (j < 2)
          s << ", ";
      }
      s << "]";
      if (i < 2)
        s << ", ";
    }
    s << "], ";

    s << "\"fibers\": [\n";

    for (int i = 0; i < nf; i++) {
      s << OPS_PRINT_JSON_MATE_INDENT << "\t{\"location\": [" 
        << (*fibers)[i].r[1] << ", " 
        << (*fibers)[i].r[2] << "], ";
      s << "\"area\": " << (*fibers)[i].area << ", ";
      s << "\"warp\": [";
      for (int j = 0; j < nwm; j++) {
        s << "[";
        for (int k = 0; k < 3; k++) {
          s << (*fibers)[i].warp[j][k];
          if (k < 2)
            s << ", ";
        }
        s << "]";
        if (j < nwm-1)
          s << ", ";
      }
      s << "], ";

      s << "\"material\": " << materials[i]->getTag();
      if (i < nf - 1)
          s << "},\n";
      else
          s << "}\n";
    }
    s << OPS_PRINT_JSON_MATE_INDENT << "]}";
    return;
  }

  else if (flag == 1) {
    for (int i = 0; i < nf; i++) {
      auto & fiber = (*fibers)[i];
      s << "\nLocation (y,z) = " << fiber.r[1] << ' ' << fiber.r[2];
      s << "\nArea = " << fiber.area << endln;
      materials[i]->Print(s, flag);
    }
  }
}
