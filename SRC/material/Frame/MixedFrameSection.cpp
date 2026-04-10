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
// Description: This file contains the class implementation of MixedFrameSection.
// MixedFrameSection provides the abstraction of a 3D beam section discretized by fibers.
//
// Warp types:
// - UT: enhanced, based on uniform torsion assumption
// - U02: warping based on uniform torsion assumption, a la Simo and Gruttmann
//
// Written: cmp
// Created: Jan. 2026
//
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <atomic>
#include <array>
#include <Channel.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include "MixedFrameSection.h"
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>
#include <GroupSO3.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <NDMaterial.h>
#include <Parameter.h>
#include <domain/DomainStatus.h>


#include <threads/thread_pool.hpp>
#define N_FIBER_THREADS 8

#define SEC_TAG_MixedFrameSection 0

using namespace OpenSees;

ID MixedFrameSection::code(nsr);

MixedFrameSection::MixedFrameSection(int tag, int num, MixedType type, bool wagner)
  : FrameSection(tag, SEC_TAG_MixedFrameSection)
  , s{}, e{}
  , e_wrap(e)
  , s_wrap(s)
  , shear_align{}
  , shift_twist{}
  , shift_axial{}
  , centroid{}
  , mixed_type(type)
  , nubar(0.0)
  , parameterID(0), dedh(nsr)
  , fibers(std::make_shared<std::vector<FiberData>>())
  , K_init(new Tangent)
  , fiber_state(FiberState::Clean)
  , wagner(wagner || (getenv("Wagner") != nullptr))
#ifdef N_FIBER_THREADS
  , num_threads(N_FIBER_THREADS)
  , pool((void*)new OpenSees::thread_pool{N_FIBER_THREADS})
#endif
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

  fibers->reserve(num);
  materials.reserve(num);
}


// Used in getCopy to create an element instance from a reference instance
MixedFrameSection::MixedFrameSection(const MixedFrameSection &other)
  : FrameSection(other.getTag(), other.getClassTag()),
    K_pres(other.K_pres),
    K_init(other.K_init),
    fibers(other.fibers),
    s(),
    e(),
    s_wrap(s),
    e_wrap(e),
    shear_align(other.shear_align),
    shift_twist(other.shift_twist),
    shift_axial(other.shift_axial),
    centroid(other.centroid),
    nubar(other.nubar),
    mixed_type(other.mixed_type),
    wagner(other.wagner),
    fiber_state(FiberState::Clean),
    parameterID(0),
    num_threads(other.num_threads),
    pool(other.pool)
{
  materials.reserve(other.materials.size());
  for (int i = 0; i < other.materials.size(); i++)
    materials.push_back(other.materials[i]->getCopy("BeamFiber"));

  this->revertToStart();
}


MixedFrameSection::~MixedFrameSection()
{
  for (auto material : materials) {
    if (material != nullptr)
      delete material;
  }
}


int
MixedFrameSection::getIntegral(Field field, State state, double& value) const
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
        double density = materials[i]->getRho();
        const double A  = (*fibers)[i].area;
        if (density != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;

    case Field::UnitY: // TODO: Centroid
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].r[0]; // - yBar;
        value += A*y;
      }
      return 0;

    case Field::UnitZ: // TODO: Centroid
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].r[1]; // - yBar;
        value += A*z;
      }
      return 0;

    case Field::UnitYY:
    case Field::UnitCentroidYY:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].r[0]
                        ; //- yBar*(Field::UnitCentroidYY == field);
        value += A*y*y;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].r[1]
                        ; //- zBar*(Field::UnitCentroidZZ == field);
        value += A*z*z;
      }
      return 0;

    default:
      return -1;
  }
  return -1;
}


int
MixedFrameSection::addFiber(MaterialBuilder& theMat, 
                              double Area, 
                              double yLoc, 
                              double zLoc)
{
  std::array<std::array<double,3>,nwm> warp{0};
  FiberData fiber {Area, warp, {yLoc, zLoc}};
  fibers->emplace_back(fiber);

  materials.emplace_back(theMat.getCopy("BeamFiber"));

  // Check for material that cant create BeamFiber copies
  if (materials[materials.size()-1] == nullptr)
    return -1;

  if (!materials.back()->threadSafe()) {
    opserr << "Material " << materials.back()->getClassType()
           << " cannot be used with MixedFrameSection as it is not thread safe.\n";
    return -1;
  }

  fiber_state = FiberState::Dirty;
  return materials.size()-1;
}



int
MixedFrameSection::setTrialSectionDeformation(const Vector &e_trial)
{
  e = e_trial;
  s.zero();
  return stateDetermination(K_pres, &s, &e, CurrentTangent);
}


int 
MixedFrameSection::checkFiberState()
{
  double area = 0.0;
  if (fiber_state == FiberState::Dirty) {
    nubar = 0.0;
    centroid.zero();
    const int nf = fibers->size();
    for (int i = 0; i < nf; i++) {
      NDMaterial &theMat = *materials[i];
      auto & fiber = (*fibers)[i];
      const Matrix & tangent = theMat.getInitialTangent();
      double nu = tangent(0,0)/(2.0*tangent(1,1)) - 1.0;
      area  += fiber.area;
      nubar += nu*fiber.area;
      centroid[1] += fiber.r[0]*fiber.area;
      centroid[2] += fiber.r[1]*fiber.area;
    }
    centroid /= area;
    nubar /= area;
  }
  

  if (mixed_type == MixedType::Energetic || mixed_type == MixedType::Constant) {
    if ((mixed_shapes&MixedShapes::ShearY) == 0)
      mixed_type = MixedType::UT;
  }
  if (mixed_type == MixedType::UT) {
    if ((mixed_shapes&MixedShapes::TwistX) == 0)
      mixed_type = MixedType::None;
  }

  fiber_state = FiberState::Clean;
  return 0;
}


FrameSection*
MixedFrameSection::getFrameCopy()
{
  if (fiber_state == FiberState::Dirty)
    this->checkFiberState();
  return new MixedFrameSection(*this);
}



int 
MixedFrameSection::formMixedUniformL(Matrix3D& Lr, Matrix3D& Lw) //const
{
  constexpr static Matrix3D oneS {{
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  }};

  Lr.zero();
  Lw.zero();

  if (mixed_type == MixedType::UT) {
    Lr(2,2) = 1.0;
    return 0;
  }
  else if (mixed_type == MixedType::Equilibrium) {
    const int nf = fibers->size();
    double Ja = 0, Jw = 0;
    for (int i = 0; i < nf; i++) {
      const auto & fiber = (*fibers)[i];
      const FiberData::WarpArray& w = fiber.warp;
      const Matrix & tangent = materials[i]->getInitialTangent();
      double E = tangent(0,0);
      double G = tangent(1,1);
      Ja -= w[0][0]*w[1][0]*fiber.area;
      Jw += w[0][0]*w[0][0]*fiber.area;
    }
    Lr(2,2) = 1.0;
    Lw(2,2) = Jw/Ja;
    return 0;
  }

  Vector3D iesan_center{};

  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {
    auto & fiber = (*fibers)[i];
    const FiberData::WarpArray& w = fiber.warp;
    const Vector3D r = {0.0, fiber.r[0], fiber.r[1]};
  
    Matrix3D Aw {{
      0.0, w[0][1], w[0][2],
      0.0, w[1][1], w[1][2],
      0.0, w[2][1], w[2][2]
    }};
    const Vector3D rxi {0.0,  r[2], -r[1]};

    Aw.addTensorProduct(r, r, -nubar);
    Aw.addMatrix(oneS, 0.5*r.dot(r)*nubar);
    // + roc 
    Aw.addTensorProduct(r, centroid, nubar);

    iesan_center += (Aw^rxi)*fiber.area;
    iesan_center[0] -= r.dot(r)*fiber.area;
  }
  const double J = -iesan_center[0];
  iesan_center[1] /= J;
  iesan_center[2] /= J;
  
  const Vector3D zg = -1.0*(shear_align*shift_twist);

  Vector3D a = shear_align^iesan_center;
  const Vector3D zk = -1.0*(shear_align*shift_axial);
  Lr(0,0) = -1.0;
  Lr(1,1) = -1.0;
  Lr(2,0) =  a[1]+zk[1];
  Lr(2,1) =  a[2]+zk[2];
  Lr(2,2) = zg.dot(iesan_center - shift_axial);

  // Lw
  Lw(0,0) =  a[1]+zk[1];
  Lw(0,1) =  a[2]+zk[2];
  Lw(0,2) = 1.0 + zg.dot(iesan_center - shift_axial);
  Lw(1,0) = shear_align(1,1);
  Lw(1,1) = shear_align(1,2);
  Lw(2,0) = shear_align(2,1);
  Lw(2,1) = shear_align(2,2);
  Lw(1,2) = zg[1];
  Lw(2,2) = zg[2];

  return 0;
}



int
MixedFrameSection::solveMixed(const VectorND<nsr> & e_trial,
                              MatrixND<6,6> & Knne,
                              Tangent & Ks)
{

  static constexpr double mixed_tol = 1e-16;
  const bool iter_ut = getenv("MixedIterUT") != nullptr;

  const Vector3D 
    gamma { e_trial(inx), e_trial(iny), e_trial(inz) },
    kappa { e_trial(imx), e_trial(imy), e_trial(imz) },
    dalpha{ e_trial(iwx), e_trial(iwy), e_trial(iwz) },
    alpha { e_trial(ivx), e_trial(ivy), e_trial(ivz) };


  Matrix3D Gr{}, Gw{};
  this->formMixedUniformL(Gr, Gw);


  const int nf = fibers->size();

  std::atomic<int> res = 0;

  struct ThreadData {
    MatrixND<3,6> Kae;
    Matrix3D Knn, Kav, Kaw;
    VectorND<3>   r_mixed;
  };
  static std::array<ThreadData, MaxThreads> thread_data;

  int iter = 0;
  bool converged = false;
  auto& thread_pool = *(OpenSees::thread_pool*)pool;

  Vector3D eta_u = eta_past;
  // if (eta_u.norm() < 1e-14 && (mixed_type == MixedType::Equilibrium)) {
    eta_u.zero();
    if ((mixed_type != MixedType::None) && (mixed_type != MixedType::U02)) {
      eta_u = Vector3D {
        gamma[1], gamma[2], kappa[0]
      };
    }
    if (mixed_type == MixedType::Equilibrium)
      eta_u[2] -= alpha[0];
  // }

  do {

    // 1. Zero for integration
    for (auto& thread : thread_data) {
      thread.Knn.zero();
      thread.Kav.zero();
      thread.Kae.zero();
      thread.Kaw.zero();
      thread.r_mixed.zero();
    }
    //
    // 2. Loop over fibers to form Knn, Kne, and s_trial
    //
    thread_pool.submit_loop<unsigned int>(0, fibers->size(), [&](unsigned int i) {

      NDMaterial &theMat = *materials[i];
      auto & fiber = (*fibers)[i];
      const Vector3D r = {0.0, fiber.r[0], fiber.r[1]};
      double tr2 = wagner? r.dot(r)*kappa[0] : 0.0;

      const double aw = tr2;

      MatrixND<3,6> Aer{};
      RigidShape(fiber, aw, Aer);
  
      Matrix3D An{}, iow{}, iodw{};
      WarpShape(fiber, iow, iodw);
      MixedShape(fiber, Gr, Gw, An);
      // Form material strain
      // Vector3D eps = gamma + kappa.cross(r);
      Vector3D eps{};
      {
        // temp: reform Aer to avoid double-counting wagner term
        MatrixND<3,6> Aer{};
        RigidShape(fiber, 0.0, Aer);
        eps.addMatrixVector(Aer, VectorND<6>{gamma[0], gamma[1], gamma[2], kappa[0], kappa[1], kappa[2]}, 1.0);
      }
      eps.addMatrixVector(An,   eta_u, 1.0);
      eps.addMatrixVector(iow, dalpha, 1.0);
      eps.addMatrixVector(iodw, alpha, 1.0);

      if (wagner) {
        tr2 = r.dot(r)*kappa[0];
        eps[0] += 0.5*kappa[0]*tr2;
      }

      res += theMat.setTrialStrain(eps);
      if (res < 0)
        return -1;


      //
      //
      //
      const Vector &stress_  = theMat.getStress();
      const Vector3D stress {stress_(0), stress_(1), stress_(2)};
      const Matrix &tangent = theMat.getTangent();

      Matrix3D C{};
      C.addMatrix(tangent, fiber.area);

      ThreadData& thread = thread_data[*OpenSees::this_thread::get_index()];
      Matrix3D&      Knn = thread.Knn;
      Matrix3D&      Kav = thread.Kav; //
      Matrix3D&      Kaw = thread.Kaw; //
      MatrixND<3,6>& Kae = thread.Kae;
      VectorND<3>&   r_thread = thread.r_mixed;

      Kae.addMatrixTripleProduct(1.0, An, C, Aer, 1.0);
      Knn.addMatrixTripleProduct(1.0, An, C, 1.0);
      Kav.addMatrixTripleProduct(1.0, An, C, iodw, 1.0);
      Kaw.addMatrixTripleProduct(1.0, An, C, iow,  1.0);


      //
      //
      //
      r_thread.addMatrixTransposeVector(An, stress, fiber.area);

      return 0;
    }).wait();

    if (res < 0) {
      return int(DomainStatus::MaterialFailedToConverge);
    }


    //
    // 3. Solve for eta_u
    //
    Matrix3D Knn{};
    MatrixND<3,6> Kae{};
    Matrix3D Kav{}, Kaw{};
    VectorND<3> r_mixed{};
    for (int t = 0; t < num_threads; t++) {
      Knn.addMatrix(thread_data[t].Knn, 1.0);
      Kae.addMatrix(thread_data[t].Kae, 1.0);
      Kav.addMatrix(thread_data[t].Kav, 1.0);
      Kaw.addMatrix(thread_data[t].Kaw, 1.0);
      r_mixed += thread_data[t].r_mixed;
    }


    Matrix3D Knn_inv{};
    if ((mixed_type == MixedType::None) ||  (mixed_type == MixedType::U02)) {
      Knne.zero();
      converged = true;
      break;
    }
    else if (mixed_type == MixedType::Constant) {
      MatrixND<3,6> Knn_inv_Kne{};
      Knn_inv_Kne(0,1) =  1.0;
      Knn_inv_Kne(1,2) =  1.0;
      Knn_inv_Kne(2,3) =  1.0;
      Knne = Kae^Knn_inv_Kne;
      converged = true;
      break;
    }
    else if (mixed_type == MixedType::Energetic) {
      Knn.invert(Knn_inv);
      Knne.addMatrixTripleProduct(0.0, Kae, Knn_inv, -1.0);
    }
    else if (mixed_type == MixedType::UT) {
      Knn_inv.zero();// = Knn;
      // Knn(2,2) is EIv
      if (iter_ut) {
        Knn_inv(2,2) = 1.0/Knn(2,2);
        Knne.addMatrixTripleProduct(0.0, Kae, Knn_inv, -1.0);
      } else {
        MatrixND<3,6> Knn_inv_Kne{};
        Knn_inv_Kne(0,1) =  1.0;
        Knn_inv_Kne(1,2) =  1.0;
        Knn_inv_Kne(2,3) =  1.0;
        Knne = Kae^Knn_inv_Kne;
        converged = true;
        break;
      }
    }
    else if (mixed_type == MixedType::Equilibrium) {
      Knn_inv.zero();// = Knn;
      Knn_inv(2,2) = 1.0/Knn(2,2);

      Knne.addMatrixTripleProduct( 0.0, Kae, Knn_inv,   -1.0);
      if (r_mixed.dot(r_mixed) < mixed_tol) {
        Ks.vv.addMatrixTripleProduct(0.0, Kav, Knn_inv,   -1.0);
        Ks.sv.addMatrixTripleProduct(0.0, Kae, Knn_inv, Kav,  -1.0);
        Ks.sw.addMatrixTripleProduct(0.0, Kae, Knn_inv, Kaw,  -1.0);
        Ks.ww.addMatrixTripleProduct(0.0, Kaw, Knn_inv,   -1.0);
        Ks.wv.addMatrixTripleProduct(0.0, Kaw, Knn_inv, Kav,  -1.0);
      }
      // converged = true;
      // break;
    }

    if (r_mixed.dot(r_mixed) < mixed_tol) {
      converged = true;
      break;
    }

    eta_u -= Knn_inv*r_mixed;

  } while (++iter < 25); //&& converged == false);


  if (!converged) {
    opserr << "  Failed to converge in section\n";
    return int(DomainStatus::SectionFailedToConverge);
  }

  eta_past = eta_u;
  return res;
}


int
MixedFrameSection::stateDetermination(Tangent& Ks, 
                                      VectorND<nsr>* s_trial, 
                                      const VectorND<nsr> * const e_trial, 
                                      int tangentFlag)
{
  
  if (fiber_state == FiberState::Dirty) [[unlikely]] {
    int res = this->checkFiberState();
  }

  const bool do_aux_warp = true;
  // const bool do_aux_warp = (mixed_type != MixedType::Constant)
  //                        &&(mixed_type != MixedType::Energetic)
  //                        &&(mixed_type != MixedType::UT);

  

  Vector3D gamma{}, kappa{}, dalpha{}, alpha{};
  if (e_trial != nullptr) {
    gamma  = Vector3D { (*e_trial)(inx), (*e_trial)(iny), (*e_trial)(inz) };
    kappa  = Vector3D { (*e_trial)(imx), (*e_trial)(imy), (*e_trial)(imz) };
    dalpha = Vector3D { (*e_trial)(iwx), (*e_trial)(iwy), (*e_trial)(iwz) };
    alpha  = Vector3D { (*e_trial)(ivx), (*e_trial)(ivy), (*e_trial)(ivz) };
  }

  Matrix3D Gr{}, Gw{};
  this->formMixedUniformL(Gr, Gw);

  //
  if (s_trial != nullptr)
    s_trial->zero();


  //
  // 2) Form eta, Kn
  //
  Ks.zero();
  MatrixND<6,6> Knne{};
  int res = this->solveMixed(e_trial != nullptr ? *e_trial : VectorND<nsr>{}, Knne, Ks);
  if (res < 0)
    return res;


  //
  // 3) Form resultants s and tangent Ke
  //
  const int nf = fibers->size();

  auto& thread_pool = *(OpenSees::thread_pool*)pool;
  struct ThreadData {
    Tangent K;
    VectorND<nsr> s_trial;
  };
  static std::array<ThreadData, MaxThreads> thread_data;
  for (auto& thread : thread_data) {
    thread.K.zero();
    thread.s_trial.zero();
  }


  thread_pool.submit_loop<unsigned int>(0, fibers->size(), [&](unsigned int i) {

    NDMaterial &theMat = *materials[i];
    auto & fiber = (*fibers)[i];
    const FiberData::WarpArray& w = fiber.warp;
    const Vector3D r = {0.0, fiber.r[0], fiber.r[1]};
    double tr2 = wagner? r.dot(r)*kappa[0] : 0.0;

    const double aw = 0;//tr2;
    MatrixND<3,6> Aer{};
    RigidShape(fiber, aw, Aer);

    Matrix3D iow{}, iodw{};
    WarpShape(fiber, iow, iodw);


    const Vector &stress  = theMat.getStress();
    const Matrix &tangent = tangentFlag==CurrentTangent
                          ? theMat.getTangent()
                          : theMat.getInitialTangent();

    Matrix3D C{};
    C.addMatrix(tangent, fiber.area);


    ThreadData& thread = thread_data[*OpenSees::this_thread::get_index()];
    Tangent& K = thread.K;
    VectorND<nsr>* s_trial = &thread.s_trial;

    K.se.addMatrixTripleProduct(1.0, Aer, C, 1.0);

    if (do_aux_warp) [[unlikely]] {
      const Matrix3D Ciow  = C*iow;
      const Matrix3D Ciodw = C*iodw;
      {
        K.sw.assemble(Ciow,  0, 0, 1.0);
        K.sv.assemble(Ciodw, 0, 0, 1.0);
        //
        K.sw.assemble(Hat(r)*Ciow,  3, 0, 1.0);
        K.sv.assemble(Hat(r)*Ciodw, 3, 0, 1.0);
        //
        K.ww.addMatrixTransposeProduct(1.0, iow,  Ciow, 1.0);
        K.wv.addMatrixTransposeProduct(1.0, iow,  Ciodw, 1.0);
        //
        K.vv.addMatrixTransposeProduct(1.0, iodw,  Ciodw, 1.0);
      }
    }

    //
    //
    //
    if (wagner && (e_trial != nullptr)) {
      constexpr Matrix3D ioi {{ 1, 0, 0 ,
                                0, 0, 0 ,
                                0, 0, 0 }};
      const Matrix3D ioiC = ioi*C;
      K.se.assemble(ioiC, 3, 0, tr2);
      K.se.assembleTranspose(ioiC, 0, 3, tr2);
      //
      K.se.assemble(Hat(r)*ioiC.transpose() - ioiC*Hat(r), 3, 3, tr2);
      K.se.assemble(ioiC*ioi, 3, 3, tr2*tr2);
      // K.mm.addMatrixProduct(ioiC, ioi, tr2*tr2);

      // Geometric part;
      K.se(3,3) += r.dot(r)*stress(0)*fiber.area;

      K.sw.assemble(ioiC*iow, 3, 0, tr2); // 6
      // K.mw.addMatrixProduct(ioiC, iow,  tr2);
      K.sv.assemble(ioiC*iodw, 3, 0, tr2); // 9
      // K.mv.addMatrixProduct(ioiC, iodw, tr2);
    }

    if (s_trial != nullptr) {
      const double y = fiber.r[0];
      const double z = fiber.r[1];
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
      for (int j=0; j<1; j++) { //nwm-nem
        // w += w[j]*s da
        (*s_trial)(iwx+j) +=  iow(0,j)*sig0;// w[j][0]*sig0;
        // v += dw[j]*s da
        (*s_trial)(ivx+j) += iodw(1,j)*sig1;// w[j][1]*sig1;
        (*s_trial)(ivx+j) += iodw(2,j)*sig2;// w[j][2]*sig2;
      }

      if (wagner && e_trial != nullptr)
        (*s_trial)(imx) += tr2*sig0;
      // if (mixed_type == MixedType::U02) {
      //   (*s_trial)(imx) += fiber.warp[0][1]*sig1 + fiber.warp[0][2]*sig2;
      // }
    }

    return 0;
  }).wait();

  // Assemble final Ke an s
  for (int t = 0; t < num_threads; t++) {
    Ks.se.addMatrix(thread_data[t].K.se, 1.0);
    Ks.sw.addMatrix(thread_data[t].K.sw, 1.0);
    Ks.sv.addMatrix(thread_data[t].K.sv, 1.0);
    Ks.ww.addMatrix(thread_data[t].K.ww, 1.0);
    Ks.wv.addMatrix(thread_data[t].K.wv, 1.0);
    Ks.vv.addMatrix(thread_data[t].K.vv, 1.0);
    if (s_trial != nullptr) {
      for (int j = 0; j < nsr; j++)
        (*s_trial)(j) += thread_data[t].s_trial(j);
    }
  }

  //
  // 4) Add Kn to Ks
  //
  Ks.se.addMatrix(Knne,  1.0);
  return 0;
}


const Vector&
MixedFrameSection::getSectionDeformation()
{
  return e_wrap;
}


MatrixND<12,12>
MixedFrameSection::getFullTangent(State state) noexcept
{
  static MatrixND<12,12> K{};
  K.zero();

  if (state == State::Init)
    this->stateDetermination(K_pres, nullptr, nullptr, InitialTangent);
  else if (fiber_state == FiberState::Dirty) {
    this->stateDetermination(K_pres, nullptr, nullptr, CurrentTangent);
  }

  K.assemble(K_pres.se, 0, 0, 1.0);
  K.assemble(K_pres.sw, 0, 6, 1.0);
  K.assemble(K_pres.sv, 0, 9, 1.0);

  K.assemble(K_pres.ww, 6, 6, 1.0);
  K.assemble(K_pres.wv, 6, 9, 1.0);
  K.assemble(K_pres.vv, 9, 9, 1.0);
  
  K.assembleTranspose(K_pres.sw, 6, 0, 1.0);
  K.assembleTranspose(K_pres.sv, 9, 0, 1.0);
  K.assembleTranspose(K_pres.wv, 9, 6, 1.0);

  return K;
}


const Matrix&
MixedFrameSection::getSectionTangent()
{
  static MatrixND<nsr,nsr> K;
  static Matrix K_wrap(K);
  K_wrap.setData(K);
  K = this->getFullTangent(State::Pres);
  return K_wrap;
}


const Matrix&
MixedFrameSection::getInitialTangent()
{
  static MatrixND<nsr,nsr> K;
  static Matrix wrap(K);
  K = this->getFullTangent(State::Init);
  return wrap;
}



const Vector&
MixedFrameSection::getStressResultant()
{
  return s_wrap;
}



const ID&
MixedFrameSection::getType()
{
  return code;
}


int
MixedFrameSection::getOrder() const
{
  return nsr;
}


int
MixedFrameSection::commitState()
{
  int err = 0;

  for (auto& material: materials)
    err += material->commitState();
  return err;
}


int
MixedFrameSection::revertToLastCommit()
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
MixedFrameSection::revertToStart()
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
MixedFrameSection::sendSelf(int commitTag, Channel &)
{
  return -1;
}

int
MixedFrameSection::recvSelf(int , Channel &,  FEM_ObjectBroker &)
{
  return -1;
}


Response*
MixedFrameSection::setResponse(const char **argv, int argc,
                                 OPS_Stream &output)
{
  Response *theResponse = nullptr;

  if (argc > 2 && strcmp(argv[0], "fiber") == 0) {

    
    int key = fibers->size();
    int passarg = 2;
    
    if (argc <= 3) {
      // fiber number was input directly
      key = atoi(argv[1]);
    }

    else if (argc > 4) {
      // find fiber closest to coord. with mat tag
      
      int matTag    = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      VectorND<2> r_given{{yCoord, zCoord}};
      double closestDist = 0;

      // Find first fiber with specified material tag
      const int nf = fibers->size();
      int j;
      for (j = 0; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          const VectorND<2> dr = fiber.r - r_given;
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
          const VectorND<2> dr = fiber.r - r_given;
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
      VectorND<2> r_given{
        atof(argv[1]), atof(argv[2])
      };

      VectorND<2> dr = (*fibers)[0].r - r_given;
      double closestDist = dr.dot(dr);
      key = 0;
      double distance;

      const int nf = fibers->size();
      for (int j = 1; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        VectorND<2> dr = fiber.r - r_given;
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
      output.attr("y",    (*fibers)[key].r[0]);
      output.attr("z",    (*fibers)[key].r[1]);
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
MixedFrameSection::getResponse(int responseID, Information &sectInfo)
{
  return FrameSection::getResponse(responseID, sectInfo);
}


int
MixedFrameSection::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strcmp(argv[0], "warp") == 0) {
    // ... warp $fiberID $warpField
    if (argc < 3) {
      opserr << "MixedFrameSection::setParameter - fiberID is required\n";
      return -1;
    }
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= fibers->size()) {
      opserr << "MixedFrameSection::setParameter - fiberID " << fiberID << " out of range\n";
      return -1;
    }

    int field = atoi(argv[2]);

    return param.addObject(Param::FiberFieldBase+fiberID*100+field, this);
  }

  else if (strcmp(argv[0], "fiber") == 0) {
    // ... fiber $fiberID $field
    if (argc < 3) {
      opserr << "MixedFrameSection::setParameter - fiberID is required\n";
      return -1;
    }
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= fibers->size()) {
      opserr << "MixedFrameSection::setParameter - fiberID " << fiberID << " out of range\n";
      return -1;
    }

    int field;
    if (strcmp(argv[2], "y") == 0)
      field = Param::FiberY;
    else if (strcmp(argv[2], "z") == 0)
      field = Param::FiberZ;
    else if (strcmp(argv[2], "area") == 0)
      field = Param::FiberArea;
    else {
      opserr << "MixedFrameSection::setParameter - invalid fiber field: " << argv[2] << "\n";
      return -1;
    }

    return param.addObject(Param::FiberFieldBase+fiberID*100+field, this);
  }
  else if (strcmp(argv[0], "shift_shear") == 0) {
    // ... shift_shear i j
    if (argc < 3) {
      opserr << "MixedFrameSection::setParameter - i, j, value are required\n";
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
      opserr << "MixedFrameSection::setParameter - invalid i, j: " << i << ", " << j << "\n";
      return -1;
    }
  }

  else if (strcmp(argv[0], "shift_axial") == 0) {
    // ... shift_axial i 
    if (argc < 2) {
      opserr << "MixedFrameSection::setParameter - i is required\n";
      return -1;
    }
    int i = atoi(argv[1]);
    switch (i) {
      case 1:
        return param.addObject(Param::ShiftAxialY, this);
      case 2:
        return param.addObject(Param::ShiftAxialZ, this);
      default:
        opserr << "MixedFrameSection::setParameter - invalid i: " << i << "\n";
        return -1;
    }
  }

  else if (strcmp(argv[0], "shift_twist") == 0) {
    // ... shift_twist i 
    if (argc < 2) {
      opserr << "MixedFrameSection::setParameter - i is required\n";
      return -1;
    }
    int i = atoi(argv[1]);
    switch (i) {
      case 1:
        return param.addObject(Param::ShiftTwistY, this);
      case 2:
        return param.addObject(Param::ShiftTwistZ, this);
      default:
        opserr << "MixedFrameSection::setParameter - invalid i: " << i << "\n";
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
MixedFrameSection::updateParameter(int paramID, Information &info)
{
  fiber_state = FiberState::Dirty;

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
  else if (paramID == Param::ShiftAxialY) {
    shift_axial[1] = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShiftAxialZ) {
    shift_axial[2] = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShiftTwistY) {
    shift_twist[1] = info.theDouble;
    return 0;
  }
  else if (paramID == Param::ShiftTwistZ) {
    shift_twist[2] = info.theDouble;
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
        (*fibers)[fiberID].r[0] = info.theDouble;
        break;
      case Param::FiberZ:
        (*fibers)[fiberID].r[1] = info.theDouble;
        break;
      case Param::FiberWarpX:
        (*fibers)[fiberID].warp[0][0] = info.theDouble;
        if (info.theDouble != 0.0)
          mixed_shapes |= MixedShapes::TwistX;
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
        if (info.theDouble != 0.0 && mixed_type == MixedType::Equilibrium)
          mixed_shapes |= MixedShapes::TwistE;
        else if (info.theDouble != 0.0)
          mixed_shapes |= MixedShapes::ShearY;
        break;
      case Param::FiberWarpYY:
        (*fibers)[fiberID].warp[1][1] = info.theDouble;
        break;
      case Param::FiberWarpYZ:
        (*fibers)[fiberID].warp[1][2] = info.theDouble;
        break;

      case Param::FiberWarpZ:
        (*fibers)[fiberID].warp[2][0] = info.theDouble;
        if (info.theDouble != 0.0)
          mixed_shapes |= MixedShapes::ShearZ;
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
MixedFrameSection::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

const Vector &
MixedFrameSection::getSectionDeformationSensitivity(int gradIndex)
{
  return dedh;
}

const Vector &
MixedFrameSection::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  VectorND<6> ds{};
  const VectorND<6> e_rigid {
    e(0), e(1), e(2), e(3), e(4), e(5)
  };

  const Vector3D dalpha { e(iwx), e(iwy), e(iwz) };
  const Vector3D alpha  { e(ivx), e(ivy), e(ivz) };

  // find deta
  Vector3D deta{};
  Matrix3D Gr{}, Gw{}, dGr{}, dGw{};
  Vector3D dcentroid{};
  double dnubar = 0.0;
  {
    double dnubar = 0.0;
    this->formMixedUniformL(Gr, Gw);
    this->formMixedUniformLSensitivity(dGr, dGw, dcentroid, dnubar);

    const Vector3D eta = eta_past;

    this->solveMixedSensitivity(gradIndex,
                                e_rigid, dalpha, alpha,
                                Gr, Gw, dGr, dGw, dcentroid, dnubar,
                                eta, deta);
  }

  const int nf = fibers->size();

  for (int i = 0; i < nf; i++) {
    const auto& fiber = (*fibers)[i];
    double dA=0, dy=0, dz=0;
    FiberData::WarpArray dw{};
    FiberGrad(i, gradIndex, dA, dy, dz, dw);
    MatrixND<3,6> Ae{}, dAe{};
    RigidShape(fiber, 0.0, Ae);
    RigidShapeGrad(fiber, 0.0, dAe, i);

    //
    Matrix3D An{}, dAn{}, iow{}, iodw{}, diow{}, diodw{};
    WarpShape(fiber, iow, iodw);
    WarpShapeGrad(fiber, diow, diodw, i);
    MixedShape(fiber, Gr, Gw, An);
    MixedShapeGrad(fiber, Gr, Gw, dGr, dGw, dcentroid, dnubar, dAn, i);
    //
    
    const Vector& dsigdh_ = materials[i]->getStressSensitivity(gradIndex,true);
    const Vector&  stress_ = materials[i]->getStress();
    const Matrix& tangent = materials[i]->getTangent();
    const Vector3D dsigdh {dsigdh_(0), dsigdh_(1), dsigdh_(2)};
    const Vector3D stress {stress_(0), stress_(1), stress_(2)};

    Matrix3D C{};
    C.addMatrix(tangent, 1.0);

    ds.addMatrixTransposeVector(Ae, dsigdh, fiber.area);
    ds.addMatrixTransposeVector(Ae, stress, dA);
    ds.addMatrixTransposeVector(dAe, stress, fiber.area);


    MatrixND<6,6> tmpMatrix{};
    tmpMatrix.addMatrixTripleProduct(0.0, Ae, C, dAe, 1.0);
    
    ds.addMatrixVector(1.0, tmpMatrix, e_rigid, fiber.area);

  }

  static Vector wrap(nsr);
  wrap.Zero();
  for (int i=0; i<6; i++)
    wrap(i) = ds(i);
  return wrap;
}

const Matrix &
MixedFrameSection::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(nsr,nsr);
  
  dksdh.Zero();
  return dksdh;
}


int
MixedFrameSection::commitSensitivity(const Vector& de,
                                     int gradIndex, int numGrads)
{
  VectorND<6> de_rigid{}, e_rigid{};
  for (int i = 0; i < 6; i++) {
    de_rigid(i) = de(i);
    e_rigid(i) = e(i);
  }

  const int nf = fibers->size();

  for (int i = 0; i < nf; i++) {
    double dA, dy, dz;
    FiberData::WarpArray dw{};
    FiberGrad(i, gradIndex, dA, dy, dz, dw);
    auto& fiber = (*fibers)[i];
    MatrixND<3,6> Ae{}, dAe{};
    RigidShape(fiber, 0.0, Ae);
    RigidShapeGrad(fiber, 0.0, dAe, i);

    // determine material strain sensitivity
    VectorND<3> deps = Ae*de_rigid+ dAe*e_rigid;
    materials[i]->commitSensitivity(deps,gradIndex,numGrads);
  }

  return 0;
}


inline void
MixedFrameSection::formMixedUniformLSensitivity(Matrix3D& dGr, Matrix3D& dGw,
                             Vector3D& dcentroid, double& dnubar) const noexcept
{
  dGr.zero();
  dGw.zero();
  dcentroid.zero();
  dnubar = 0.0;
}


inline int
MixedFrameSection::applyMixedInverse(const Matrix3D& Knn,
                                     const Vector3D& rhs,
                                     Vector3D& x) const
{
  x.zero();

  if (mixed_type == MixedType::None || 
      mixed_type == MixedType::Constant)
    return 0;

  if (mixed_type == MixedType::Energetic) {
    Matrix3D KnnInv{};
    Knn.invert(KnnInv);
    x = KnnInv * rhs;
    return 0;
  }

  // UT and Equilibrium: only the 3rd mixed component is active.
  x(2) = rhs(2) / Knn(2,2);
  return 0;
}

int
MixedFrameSection::solveMixedSensitivity(int gradIndex,
                                         const VectorND<6>& e_rigid,
                                         const Vector3D& dalpha,
                                         const Vector3D& alpha,
                                         const Matrix3D& Gr,
                                         const Matrix3D& Gw,
                                         const Matrix3D& dGr,
                                         const Matrix3D& dGw,
                                         const Vector3D& dcentroid,
                                         double dnubar,
                                         const Vector3D& eta,
                                         Vector3D& deta) const
{
  deta.zero();

  if (mixed_type == MixedType::None || mixed_type == MixedType::Constant)
    return 0;

  Matrix3D Knn{};
  Vector3D rhs{};

  const int nf = fibers->size();
  for (int i = 0; i < nf; ++i) {
    const auto& fiber = (*fibers)[i];

    double dA = 0.0, dy = 0.0, dz = 0.0;
    FiberData::WarpArray dw{};
    FiberGrad(i, 0, dA, dy, dz, dw);

    MatrixND<3,6> Ae{}, dAe{};
    Matrix3D An{}, dAn{}, iow{}, iodw{}, diow{}, diodw{};

    RigidShape(fiber, 0.0, Ae);
    RigidShapeGrad(fiber, 0.0, dAe, i);
    WarpShape(fiber, iow, iodw);
    WarpShapeGrad(fiber, diow, diodw, i);
    MixedShape(fiber, Gr, Gw, An);
    MixedShapeGrad(fiber, Gr, Gw, dGr, dGw, dcentroid, dnubar, dAn, i);

    const Vector& stress_  = materials[i]->getStress();
    const Vector& dsigdh_  = materials[i]->getStressSensitivity(gradIndex, true);
    const Matrix& tangent_ = materials[i]->getTangent();

    const Vector3D stress    {stress_(0), stress_(1), stress_(2)};
    const Vector3D dsigdhMat {dsigdh_(0), dsigdh_(1), dsigdh_(2)};

    Matrix3D C{}, CA{};
    C.addMatrix(tangent_, 1.0);
    CA.addMatrix(tangent_, fiber.area);

    Vector3D depsExp{};
    depsExp.addMatrixVector(dAe,   e_rigid, 1.0);
    depsExp.addMatrixVector(dAn,   eta,     1.0);
    depsExp.addMatrixVector(diow,  dalpha,  1.0);
    depsExp.addMatrixVector(diodw, alpha,   1.0);

    Vector3D dsig0 = dsigdhMat;
    dsig0.addMatrixVector(C, depsExp, 1.0);

    Knn.addMatrixTripleProduct(1.0, An, CA, 1.0);
    rhs.addMatrixTransposeVector(An,  dsig0, fiber.area);
    rhs.addMatrixTransposeVector(dAn, stress, fiber.area);
    rhs.addMatrixTransposeVector(An,  stress, dA);
  }

  Vector3D tmp{};
  this->applyMixedInverse(Knn, rhs, tmp);
  deta = -1.0*tmp;
  return 0;
}

 int
MixedFrameSection::WarpShapeGrad(const FiberData& fiber,
              Matrix3D& diow, Matrix3D& diodw,
              int i) const noexcept
{
  diow.zero();
  diodw.zero();

  double dA = 0.0, dy = 0.0, dz = 0.0;
  FiberData::WarpArray dw{};
  FiberGrad(i,0, dA, dy, dz, dw);

  switch (mixed_type) {
    case MixedType::UT:
    case MixedType::U02:
    case MixedType::Energetic:
    case MixedType::Constant:
      return 0;

    case MixedType::None:
      diow(0,1)  = dw[1][0];
      diow(0,2)  = dw[2][0];
      diodw(1,1) = dw[1][1];
      diodw(2,1) = dw[1][2];
      diodw(1,2) = dw[2][1];
      diodw(2,2) = dw[2][2];
      [[fallthrough]];

    case MixedType::Equilibrium:
      diow(0,0)  = dw[0][0];
      diodw(1,0) = dw[0][1];
      diodw(2,0) = dw[0][2];
      return 1;
  }
  return 0;
}


inline void
MixedFrameSection::MixedShapeGrad(const FiberData& fiber,
               const Matrix3D& Gr,  const Matrix3D& Gw,
               const Matrix3D& dGr, const Matrix3D& dGw,
               const Vector3D& dcentroid, double dnubar,
               Matrix3D& dAn, int i) const noexcept
{
  dAn.zero();

  double dA = 0.0, dy = 0.0, dz = 0.0;
  FiberData::WarpArray dw{};
  FiberGrad(i, 0, dA, dy, dz, dw);

  if (mixed_type == MixedType::None)
    return;

  if (mixed_type == MixedType::UT) {
    dAn(1,2) = dw[0][1];
    dAn(2,2) = dw[0][2];
    return;
  }

  if (mixed_type == MixedType::Equilibrium) {
    const double b  = Gw(2,2);
    const double db = dGw(2,2);
    dAn(1,2) = dw[0][1] + db*fiber.warp[1][1] + b*dw[1][1];
    dAn(2,2) = dw[0][2] + db*fiber.warp[1][2] + b*dw[1][2];
    return;
  }

  constexpr static Matrix3D oneS {{
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  }};

  const Vector3D r  {0.0, fiber.r[0], fiber.r[1]};
  const Vector3D dr {0.0, dy,        dz       };

  Matrix3D Anr {{
    0.0,   1.0,  0.0,
    0.0,   0.0,  1.0,
    0.0, -r[2],  r[1]
  }};
  Matrix3D dAnr{};
  dAnr(2,1) = -dz;
  dAnr(2,2) =  dy;

  Matrix3D Anwo {{
    0.0, fiber.warp[0][1], fiber.warp[0][2],
    0.0, fiber.warp[1][1], fiber.warp[1][2],
    0.0, fiber.warp[2][1], fiber.warp[2][2]
  }};
  Anwo.addTensorProduct(r, r, -nubar);
  Anwo.addMatrix(oneS, 0.5*r.dot(r)*nubar);
  Anwo.addTensorProduct(r, centroid, nubar);

  Matrix3D dAnwo {{
    0.0, dw[0][1], dw[0][2],
    0.0, dw[1][1], dw[1][2],
    0.0, dw[2][1], dw[2][2]
  }};
  dAnwo.addTensorProduct(dr, r, -nubar);
  dAnwo.addTensorProduct(r, dr, -nubar);
  dAnwo.addMatrix(oneS, nubar*r.dot(dr));
  dAnwo.addTensorProduct(dr, centroid, nubar);

  dAnwo.addTensorProduct(r, r, -dnubar);
  dAnwo.addMatrix(oneS, 0.5*r.dot(r)*dnubar);
  dAnwo.addTensorProduct(r, centroid, dnubar);
  dAnwo.addTensorProduct(r, dcentroid, nubar);

  dAn.addMatrixProduct(dAnr, Gr, 1.0);
  dAn.addMatrixProduct(Anr, dGr, 1.0);
  dAn.addMatrixProduct(dAnwo, Gw, 1.0);
  dAn.addMatrixProduct(Anwo, dGw, 1.0);
}


void
MixedFrameSection::Print(OPS_Stream &s, int flag)
{
  const int nf = fibers->size();
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    double mass;
    if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0)
      s << "\"mass\": " << mass << ", ";

    s << "\"wagner\": " << (wagner ? "true" : "false") << ", ";
    s << "\"warp_type\": \"";
    switch (mixed_type) {
      case MixedType::None:        s << "None"; break;
      case MixedType::UT:          s << "UT"; break;
      case MixedType::U02:         s << "U02"; break;
      case MixedType::Equilibrium: s << "NR"; break;
      case MixedType::Energetic:   s << "UE"; break;
      case MixedType::Constant:    s << "UG"; break;
    }
    s << "\", ";

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
    s << "\"shift_axial\": [" << shift_axial[1] << ", " << shift_axial[2] << "], ";
    s << "\"shift_twist\": [" << shift_twist[1] << ", " << shift_twist[2] << "], ";

    s << "\"fibers\": [\n";

    for (int i = 0; i < nf; i++) {
      s << OPS_PRINT_JSON_MATE_INDENT << "\t{\"location\": [" 
        << (*fibers)[i].r[0] << ", " 
        << (*fibers)[i].r[1] << "], ";
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
      s << "\nLocation (y,z) = " << fiber.r[0] << ' ' << fiber.r[1];
      s << "\nArea = " << fiber.area << endln;
      materials[i]->Print(s, flag);
    }
  }
}
