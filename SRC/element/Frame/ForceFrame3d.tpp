//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//         Please cite the following resource in any derivative works:
//                 https://doi.org/10.5281/zenodo.10456866
//
//===----------------------------------------------------------------------===//
//
// Three dimensional force-interpolated frame element. 
//
// The general state-determination algorithm is developed by Spacone et al (1996).
// Extension for shear-deformable kinematics is developed by Lee and Filippou.
//
// Primary References
// ==================
//
//  Spacone, E., V. Ciampi, and F. C. Filippou (1996). 
//    "Mixed Formulation of Nonlinear Beam Finite Element."
//    Computers and Structures, 58(1):71-83.
//
//  Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
//    "Response Sensitivity for Nonlinear Beam-Column Elements."
//    Journal of Structural Engineering, 130(9):1281-1288.
//
// See also
// ========
//  Scott, M. H. and G. L. Fenves (2006). 
//    "Plastic Hinge Integration Methods for Force-Based Beam-Column Elements." 
//    Journal of Structural Engineering, 132(2):244-252.
//
//===----------------------------------------------------------------------===//
//
//  TODO:
//    how should a dense mass matrix be formed when shear is present?
//
//===----------------------------------------------------------------------===//
//
// fcf, es, mhs, fmk, cmp
//
#include <cmath>
#include <string.h>
#include <cfloat>
#include <array>

#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <ForceFrame3d.h>
#include <BeamIntegration.h>
#include <FrameSection.h>
#include <interpolate/cbdi.h>
#include <Node.h>
#include <Domain.h>
#include <Cholesky.tpp>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>

#include <BasicFrameTransf.h>
#include <runtime/commands/modeling/transform/FrameTransformBuilder.hpp>

#define ELE_TAG_ForceFrame3d 0 // TODO


template <int NIP, int nsr, int nwm>
ForceFrame3d<NIP,nsr,nwm>::ForceFrame3d(int tag, 
                           std::array<int, 2>& nodes, 
                           std::vector<FrameSection*>& sec,
                           BeamIntegration& bi,
                           FrameTransformBuilder& tb, 
                           double dens, int mass_flag, bool use_density,
                           int max_iter_, double tolerance
                           )
 :
   BasicFrame3d(),
   FiniteElement<2, 3, 6+nwm>(tag, ELE_TAG_ForceFrame3d, nodes),
   basic_system(new BasicFrameTransf3d<NDF>(tb.template create<2,NDF>())),
   stencil(nullptr),
   state_flag(0),
   Ki(nullptr),
   density(dens), mass_flag(mass_flag), use_density(use_density),
   mass_initialized(false), twist_mass(0), total_mass(0.0),
   max_iter(max_iter_), tol(tolerance),
   parameterID(0)
{
  K_pres.zero();
  K_save.zero();
  q_save.zero();
  q_pres.zero();

  stencil = bi.getCopy();

  this->setSectionPointers(sec);
}


template <int NIP, int nsr, int nwm>
ForceFrame3d<NIP,nsr,nwm>::~ForceFrame3d()
{

  for (GaussPoint& point : points)
    if (point.material != nullptr)
      delete point.material;

  delete stencil;

  delete Ki;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::setNodes()
{

  Node** theNodes = this->getNodePtrs();
  if (basic_system->initialize(theNodes[0], theNodes[1]) != 0) {
      return -1;
  }

  double L = basic_system->getInitialLength();
  this->BasicFrame3d::setLength(L);

  if (L == 0.0)
    return -1;

  int numSections = points.size();
  double *xi = new double[numSections];
  double *wt = new double[numSections];
  stencil->getSectionLocations(numSections, L, xi);
  stencil->getSectionWeights(numSections, L, wt);
  for (int i=0; i<numSections; i++) {
    points[i].point  = xi[i];
    points[i].weight = wt[i];
  }
  delete[] xi;
  delete[] wt;


  if (state_flag == 0)
    this->initializeSectionHistoryVariables();

//this->revertToStart();

  return 0;
}

template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::getIntegral(Field field, State state, double& total)
{

  if (this->setState(State::Init) != 0)
    return -1;

  total = 0.0;
  switch (field) {

    // Integrate density to compute total mass
    case Field::Density: {
      for (GaussPoint& sample : points) {
        double value = 0.0;
        // use element density if supplied
        if (use_density)
          total += sample.weight*density;

        // try using section's internal density
        else if (sample.material->getIntegral(Field::Density, state, value) == 0)
          total += sample.weight*value;

        else
          return -1;
      }
      return 0;
    }

    case Field::PolarInertia: {
      for (GaussPoint& sample : points) {
        double A;
        if (sample.material->getIntegral(Field::UnitYY, state, A) != 0)
          continue;

        // Get \int \rho y^2
        double Iz;
        if (sample.material->getIntegral(Field::DensityYY, state, Iz) != 0) {
          // Section does not support integrating density; try
          // integrating product of inertia and multiplying by rho
          if (use_density && sample.material->getIntegral(Field::UnitYY, state, Iz) == 0)
            Iz *= density/A;
          else
            continue;
        }
        // Get \int \rho z^2
        double Iy;
        if (sample.material->getIntegral(Field::DensityZZ, state, Iy) != 0) {
          if (use_density && sample.material->getIntegral(Field::UnitZZ, state, Iy) == 0)
            Iy *= density/A;
          else
            continue;
        }
        total += sample.weight*(Iy + Iz);
      }
      return 0;
    }

    default:
      return -1;
  }
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::commitState()
{
  int status = 0;

  // Call element commitState to do any base class stuff
  if ((status = this->Element::commitState()) != 0)
    opserr << "ForceFrame3d::commitState - failed in base class";

  for (GaussPoint& point : points) {
    point.es_save = point.es;
    if (point.material->commitState() != 0)
      return -1;
  }


  // commit the transformation between coord. systems
  if ((status = basic_system->commitState()) != 0)
    return status;

  // commit the element variables state
  K_save = K_pres;
  q_save = q_pres;

  // state_flag = 0; fmk - commented out, see what happens to Example3.1.tcl if uncommented
  //                     - i have not a clue why, ask remo if he ever gets in contact with us again!

  return status;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::revertToLastCommit()
{
  for (GaussPoint& point : points) {
    FrameSection& section = *point.material;

    point.es = point.es_save;
    if (section.revertToLastCommit() != 0)
      return -1;

    section.setTrialState<nsr,scheme>(point.es);
    point.sr = section.getResultant<nsr,scheme>();
    point.Fs = section.template getFlexibility<nsr,scheme>();

  }

  // Revert the transformation to last commit
  if (basic_system->revertToLastCommit() != 0)
    return -2;


  // Revert the element state to last commit
  q_pres = q_save;
  K_pres = K_save;

  state_flag = 0;

  return 0;
}



template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::revertToStart()
{
  // Revert the transformation to start
  if (basic_system->revertToStart() != 0)
    return -2;

  // Revert the sections state to start
  for (GaussPoint& point : points) {
    point.Fs.zero();
    point.es.zero();
    point.sr.zero();
    if (point.material->revertToStart() != 0)
      return -1;
  }

  // Revert the element state to start
  q_save.zero();
  K_save.zero();

  q_pres.zero();
  K_pres.zero();

  state_flag = 0;
  // this->update();
  return 0;
}


// template <int NIP, int nsr, int nwm>
// VectorND<NBV>&
// ForceFrame3d<NIP,nsr,nwm>::getBasicForce()
// {
//   return q_pres;
// }


// template <int NIP, int nsr, int nwm>
// MatrixND<NBV, NBV>&
// ForceFrame3d<NIP,nsr,nwm>::getBasicTangent(State state, int rate)
// {
//   return K_pres;
// }


template <int NIP, int nsr, int nwm>
const Matrix &
ForceFrame3d<NIP,nsr,nwm>::getMass()
{

  if (!mass_initialized) {
    total_mass = 0.0;
    if (this->getIntegral(Field::Density, State::Init, total_mass) != 0)
      total_mass = 0.0;

    // Twisting inertia
    twist_mass = 0.0;
    if (this->getIntegral(Field::PolarInertia, State::Init, twist_mass) != 0)
      twist_mass = 0.0;

    mass_initialized = true;
  }

  if (total_mass == 0.0) {
      ALWAYS_STATIC MatrixND<2*NDF,2*NDF> M{};
      ALWAYS_STATIC Matrix Wrapper{M};
      return Wrapper;
  }

  else if (mass_flag == 0)  {

      ALWAYS_STATIC MatrixND<2*NDF,2*NDF> M{};
      ALWAYS_STATIC Matrix Wrapper{M};
      // lumped mass matrix
      double m = 0.5*total_mass;
      M(0,0) = m;
      M(1,1) = m;
      M(2,2) = m;
      M(6,6) = m;
      M(7,7) = m;
      M(8,8) = m;
      return Wrapper;
  }

  else {
      // consistent (cubic, prismatic) mass matrix

      double L  = basic_system->getInitialLength();
      double m  = total_mass/420.0;
      double mx = twist_mass;
      ALWAYS_STATIC MatrixND<2*NDF,2*NDF> ml{};
      ALWAYS_STATIC Matrix mg{ml};

      ml(0,0) = ml(6,6) = m*140.0;
      ml(0,6) = ml(6,0) = m*70.0;

      ml(3,3) = ml(9,9) = mx/3.0; // Twisting
      ml(3,9) = ml(9,3) = mx/6.0;

      ml( 2, 2) = ml( 8, 8) =  m*156.0;
      ml( 2, 8) = ml( 8, 2) =  m*54.0;
      ml( 4, 4) = ml(10,10) =  m*4.0*L*L;
      ml( 4,10) = ml(10, 4) = -m*3.0*L*L;
      ml( 2, 4) = ml( 4, 2) = -m*22.0*L;
      ml( 8,10) = ml(10, 8) = -ml( 2, 4);
      ml( 2,10) = ml(10, 2) =  m*13.0*L;
      ml( 4, 8) = ml( 8, 4) = -ml( 2,10);

      ml( 1, 1) = ml( 7, 7) =  m*156.0;
      ml( 1, 7) = ml( 7, 1) =  m*54.0;
      ml( 5, 5) = ml(11,11) =  m*4.0*L*L;
      ml( 5,11) = ml(11, 5) = -m*3.0*L*L;
      ml( 1, 5) = ml( 5, 1) =  m*22.0*L;
      ml( 7,11) = ml(11, 7) = -ml(1,5);
      ml( 1,11) = ml(11, 1) = -m*13.0*L;
      ml( 5, 7) = ml( 7, 5) = -ml(1,11);

      // transform local mass matrix to global system
      ml = basic_system->t.pushConstant(ml);
      return mg;
  }
}



template <int NIP, int nsr, int nwm>
const Matrix &
ForceFrame3d<NIP,nsr,nwm>::getInitialStiff()
{
  // check for quick return
  if (Ki != nullptr)
    return *Ki;

  MatrixND<NBV,NBV> F_init;
  this->getInitialFlexibility(F_init);

  // calculate element stiffness matrix
  ALWAYS_STATIC MatrixND<NBV, NBV> K_init;
  if (Cholesky<NBV>(F_init).invert(K_init) < 0)
    opserr << "ForceFrame3d: Failed to invert flexibility";

  Matrix wrapper(K_init);
  Ki = new Matrix(basic_system->getInitialGlobalStiffMatrix(wrapper));

  return *Ki;
}



template <int NIP, int nsr, int nwm>
void
ForceFrame3d<NIP,nsr,nwm>::initializeSectionHistoryVariables()
{
  const int nip = points.size();
  for (int i = 0; i < nip; i++) {
    points[i].Fs.zero();
    points[i].es.zero();
    points[i].sr.zero();
    points[i].es_save.zero();
  }
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::update()
{
  constexpr static double TOL_SUBDIV = DBL_EPSILON;

  // TODO: remove hard limit on sections
  THREAD_LOCAL VectorND<nsr>     es_trial[NIP]; //  strain
  THREAD_LOCAL VectorND<nsr>     sr_trial[NIP]; //  stress resultant
  THREAD_LOCAL MatrixND<nsr,nsr> Fs_trial[NIP]; //  flexibility


  // If we have completed a recvSelf() do a revertToLastCommit()
  // to get sr, etc. set correctly
  if (state_flag == 2)
    this->revertToLastCommit();

  //
  // Localize deformations
  //
  basic_system->update();

  double L   = basic_system->getInitialLength();
  double jsx = 1.0 / L;

  THREAD_LOCAL VectorND<NBV> dv;
  dv = basic_system->getBasicIncrDeltaDisp();

  //
  //
  //
  if (state_flag != 0 && (dv.norm() <= DBL_EPSILON) && (eleLoads.size()==0))
    return 0;

  THREAD_LOCAL VectorND<NBV> Dv;
  {
    const Vector& v = basic_system->getBasicTrialDisp();
    Dv  = v;
    Dv -= dv;
  }

  for (int i=0; i<nwm; i++) {
    Node** nodes = this->getNodePtrs();
    for (int j=0; j<2; j++) {
      const Vector& uj  = nodes[j]->getTrialDisp();
      const Vector& duj = nodes[j]->getIncrDeltaDisp();
      dv[NNW+i*nwm+j] = duj[6+i];
      Dv[NNW+i*nwm+j] = uj[6+i] - duj[6+i];
    }
  }

  THREAD_LOCAL VectorND<NBV> dvToDo,
                             dv_trial;

  dvToDo  = dv;
  dv_trial = dvToDo;

  static constexpr double factor = 10.0;
  double dW;             // section strain energy (work) norm
  double dW0  = 0.0;



  //
  //   Iterate to find compatible forces and deformations
  //
  //   First try first a Newton iteration, if that fails we try an initial
  //   flexibility iteration on first iteration and then regular Newton, if
  //   that fails we use the initial flexiblity for all iterations.
  //   If they both fail we subdivide dV & try to get compatible forces
  //   and deformations. If they work, and we have subdivided we apply
  //   the remaining dV.
  //
  enum class Strategy {
    Newton, InitialIterations, InitialThenNewton
  };
  static constexpr std::array<Strategy,3> solve_strategy {
    Strategy::Newton, Strategy::InitialIterations, Strategy::InitialThenNewton
  };

  int subdivision = 1;
  bool converged = false;

  const int nip = points.size();
  while (converged == false && subdivision <= max_subdivision) {

    for (Strategy strategy : solve_strategy ) {

      // Allow 10 times more iterations for initial tangent strategy
      const int numIters = (strategy==Strategy::InitialIterations) ? 10*max_iter : max_iter;

      for (int i = 0; i < nip; i++) {
        es_trial[i]  = points[i].es;
        Fs_trial[i]  = points[i].Fs;
        sr_trial[i]  = points[i].sr;
      }

      if (state_flag == 2)
        continue;

      VectorND<NBV>      q_trial = q_pres;
      MatrixND<NBV, NBV> K_trial = K_pres;

      q_trial += K_pres*dv;

      for (int j = 0; j < numIters; j++) {

        THREAD_LOCAL VectorND<NBV> vr;       // element residual deformations
        THREAD_LOCAL MatrixND<NBV, NBV> F;   // element flexibility matrix
        F.zero();
        vr.zero();

        //
        // Gauss Loop
        //
        for (int i = 0; i < nip; i++) {
          double xL = points[i].point;
          double xL1 = xL - 1.0;
          double wtL = points[i].weight * L;
          // Retrieve section flexibility, deformations, and forces from last iteration
          auto& Fs = Fs_trial[i];
          auto& sr = sr_trial[i];
          FrameSection& section = *points[i].material;

          //
          // a. Calculate section force by interpolation of q_trial
          //
          //    si = b*q + bp*w;

          // Interpolation of q_trial
          //    b*q_trial
          //
          VectorND<nsr> si;
          si.zero();
          for (int ii=0; ii<nsr; ii++) {
            switch (scheme[ii]) {
              case FrameStress::N:
                si[ii] = q_trial[jnx];
                break;
              case FrameStress::Vy:
                si[ii] = (q_trial[imz] + q_trial[jmz])/L;
                break;
              case FrameStress::Vz:
                si[ii] = (q_trial[imy] + q_trial[jmy])/L;
                break;
              case FrameStress::T:
                si[ii] = q_trial[jmx];
                break;
              case FrameStress::My:
                si[ii] = xL1 * q_trial[imy] + xL * q_trial[jmy];
                break;
              case FrameStress::Mz:
                si[ii] = xL1 * q_trial[imz] + xL * q_trial[jmz];
                break;
              case FrameStress::Bimoment:
                si[ii] = xL1 * q_trial[iwx] + xL * q_trial[jwx];
                break;
              case FrameStress::Bishear:
                si[ii] = (q_trial[iwx] + q_trial[jwx])/L;
                break;
              default:
                break;
            }
          }

          // Add the particular solution
          // si += bp*w
          if (eleLoads.size() != 0)
            this->addLoadAtSection(si, points[i].point * L);


          //
          // b. Compute section deformations es_trial
          //
          //    es += Fs * ( si - sr(e) );
          //
          if (state_flag != 0) {

            // Form stress increment ds from last iteration
            // ds = si - si_past;
            VectorND<nsr> ds = si - sr;

            // Add strain correction
            //    es += Fs * ds;
            switch (strategy) {
              case Strategy::Newton:
                //  regular Newton
                es_trial[i].addMatrixVector(1.0, Fs, ds, 1.0);
                break;

              case Strategy::InitialThenNewton:
                //  Newton with initial tangent if first iteration
                //  otherwise regular Newton
                if (j == 0) {
                  MatrixND<nsr,nsr> Fs0 = section.template getFlexibility<nsr,scheme>(State::Init);
                  es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
                } else
                  es_trial[i].addMatrixVector(1.0, Fs, ds, 1.0);
                break;

              case Strategy::InitialIterations:
                //  Newton with initial tangent
                MatrixND<nsr,nsr> Fs0 = section.template getFlexibility<nsr,scheme>(State::Init);
                es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
                break;
            }
          }


          //
          // c. Set trial section state and get response
          //
          if (section.setTrialState<nsr,scheme>(es_trial[i]) < 0) {
            opserr << "ForceFrame3d: section " << i << " failed in setTrial\n";
            return -1;
          }

          sr_trial[i] = section.getResultant<nsr, scheme>();
          Fs_trial[i] = section.template getFlexibility<nsr, scheme>();

          //
          // d. Integrate element flexibility matrix
          //
          //    F += (B' * Fs * B) * wi * L;
          //
          {
            MatrixND<nsr,NBV> FsB;
            FsB.zero();
            for (int ii = 0; ii < nsr; ii++) {
              for (int jj = 0; jj < nsr; jj++) {
                switch (scheme[jj]) {
                case FrameStress::N:
                  FsB(ii, jnx) += Fs(ii, jj) * wtL;
                  break;

                case FrameStress::Vy:
                  FsB(ii, imz) += Fs(ii, jj) * wtL * jsx;
                  FsB(ii, jmz) += Fs(ii, jj) * wtL * jsx;
                  break;

                case FrameStress::Vz:
                  FsB(ii, imy) += Fs(ii, jj) * wtL * jsx;
                  FsB(ii, jmy) += Fs(ii, jj) * wtL * jsx;
                  break;

                case FrameStress::T:
                  FsB(ii, jmx) += Fs(ii, jj) * wtL;
                  break;

                case FrameStress::My:
                  FsB(ii, imy) += xL1 * Fs(ii, jj)*wtL;
                  FsB(ii, jmy) += xL  * Fs(ii, jj)*wtL;
                  break;

                case FrameStress::Mz:
                  FsB(ii, imz) += xL1 * Fs(ii, jj)*wtL;
                  FsB(ii, jmz) += xL  * Fs(ii, jj)*wtL;
                  break;

                case FrameStress::Bimoment:
                  FsB(ii, iwx) += xL1 * Fs(ii, jj)*wtL;
                  FsB(ii, jwx) += xL  * Fs(ii, jj)*wtL;
                  break;
                case FrameStress::Bishear:
                  FsB(ii, iwx) += Fs(ii, jj) * wtL * jsx;
                  FsB(ii, jwx) += Fs(ii, jj) * wtL * jsx;
                  break;
                }
              }
            }

            for (int jj = 0; jj < nsr; jj++) {
              for (int ii = 0; ii < NBV; ii++) {
                switch (scheme[jj]) {
                case FrameStress::N:
                  F(jnx, ii) += 1.0 * FsB(jj, ii); // Nj
                  break;
                case FrameStress::Vy:
                  F(imz, ii) += jsx * FsB(jj, ii);
                  F(jmz, ii) += jsx * FsB(jj, ii);
                  break;
                case FrameStress::Vz:
                  F(imy, ii) += jsx * FsB(jj, ii);
                  F(jmy, ii) += jsx * FsB(jj, ii);
                  break;
                case FrameStress::T:
                  F(jmx, ii) += 1.0 * FsB(jj, ii); //
                  break;
                case FrameStress::My:
                  F(imy, ii) += xL1 * FsB(jj, ii);
                  F(jmy, ii) += xL  * FsB(jj, ii);
                  break;
                case FrameStress::Mz:
                  F(imz, ii) += xL1 * FsB(jj, ii);
                  F(jmz, ii) += xL  * FsB(jj, ii);
                  break;
                case FrameStress::Bimoment:
                  F(iwx, ii) += xL1 * FsB(jj, ii);
                  F(jwx, ii) += xL  * FsB(jj, ii);
                  break;
                case FrameStress::Bishear:
                  F(iwx, ii) += jsx * FsB(jj, ii);
                  F(jwx, ii) += jsx * FsB(jj, ii);
                  break;
                }
              }
            }
          }


          //
          // e. Integrate residual deformations
          //
          //    vr += (B' * (es + des)) * wi * L;
          //
          {
            VectorND<nsr> des, ds;
            // calculate section residual deformations
            // des = Fs * ds,  with  ds = si - sr[i];
            ds = si;
            ds.addVector(1.0, sr, -1.0);

            des.addMatrixVector(0.0, Fs, ds, 1.0);
            des.addVector(1.0, es_trial[i], 1.0);

            // B' * des 
            for (int ii = 0; ii < nsr; ii++) {
              switch (scheme[ii]) {
              case FrameStress::N:
                vr[jnx] += des[ii]*wtL;
                break;
              case FrameStress::Vy:
                vr[imz] += jsx*des[ii]*wtL;
                vr[jmz] += jsx*des[ii]*wtL;
                break;
              case FrameStress::Vz:
                vr[imy] += jsx*des[ii]*wtL;
                vr[jmy] += jsx*des[ii]*wtL;
                break;
              case FrameStress::T:
                vr[jmx] +=       des[ii]*wtL;
                break;
              case FrameStress::My:
                vr[imy] += xL1 * des[ii]*wtL;
                vr[jmy] += xL  * des[ii]*wtL;
                break;
              case FrameStress::Mz:
                vr[imz] += xL1 * des[ii]*wtL;
                vr[jmz] += xL  * des[ii]*wtL;
                break;
              case FrameStress::Bimoment:
                vr[iwx] += xL1 * des[ii]*wtL;
                vr[jwx] += xL  * des[ii]*wtL;
                break;
              case FrameStress::Bishear:
                vr[iwx] += jsx*des[ii]*wtL;
                vr[jwx] += jsx*des[ii]*wtL;
                break;
              }
            }
          }
        } // Gauss loop


        // dv = Dv + dv_trial  - vr
        dv = Dv;
        dv += dv_trial;
        dv -= vr;

        //
        // Finalize trial element state 
        //
        //    K_trial  = inv(F)
        //    q_trial += K * (Dv + dv_trial - vr)
        //

        if (Cholesky<NBV>(F).invert(K_trial) < 0) {
          if constexpr (NBV < 7) {
            if (F.invert(K_trial) < 0)
              return -1;
          }
          if constexpr (NBV >= 7) {
            K_trial = F;
            if (K_trial.invert() < 0) {
              opserr << "ForceFrame3d: Failed to invert flexibility\n";
              return -1;
            }
          }
        }

        VectorND<NBV> dqe = K_trial * dv;

        dW = dqe.dot(dv);
        if (dW0 == 0.0)
          dW0 = dW;

        q_trial += dqe;


        //
        // Check for convergence of this interval
        //
        if (std::fabs(dW) < tol) {

          // Set the target displacement
          dvToDo -= dv_trial;
          Dv     += dv_trial;

          // Check if we have got to where we wanted
          if (dvToDo.norm() <= TOL_SUBDIV) {
            converged = true;

          } else {
            // We've converged but we have more to do;
            // reset variables for start of next subdivision
            dv_trial = dvToDo;
            // NOTE setting subdivide to 1 again maybe too much
            subdivision = 1; 
          }

          // set K_pres, es and q_pres values
          K_pres = K_trial;
          q_pres = q_trial;

          for (int k = 0; k < nip; k++) {
            points[k].es  = es_trial[k];
            points[k].Fs  = Fs_trial[k];
            points[k].sr  = sr_trial[k];
          }

          // break out of j & l loops
          goto iterations_completed;
        }
        else { //  if (fabs(dW) < tol) {

          // if we have failed to converge for all of our Newton schemes
          // - reduce step size by the factor specified
          if (j == (numIters - 1) && (strategy == Strategy::InitialThenNewton)) {
            dv_trial /= factor;
            subdivision++;
          }
        }
      } // for (iteration)
    }   // for (strategy)

iterations_completed:
        ;

  } // while (converged == false)


  if (converged == false) {
    opserr << "WARNING - ForceFrame3d failed internal state determination ";
    opserr << "for element " 
           << this->getTag() 
           << "; dW = " << dW 
           << ", dW0 = " << dW0
           << "\n";
    return -1;
  }

  state_flag = 1;

  return 0;
}


template <int NIP, int nsr, int nwm>
const Matrix &
ForceFrame3d<NIP,nsr,nwm>::getTangentStiff()
{
  MatrixND<NBV,NBV> kb = K_pres; // TODO!!! = this->getBasicTangent(State::Pres, 0);

  VectorND<NDF*2> pl{};
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];      // Nj
  pl[1*NDF+3]  =  q_pres[jmx];      // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];
  for (int i=0; i<nwm; i++) {
    pl[0*NDF+6+i] = -q_pres[NNW+i];
    pl[1*NDF+6+i] =  q_pres[NNW+i];
  }
  //
  pl[0*NDF+0]  = -q_pres[jnx];      // Ni
  pl[0*NDF+3]  = -q_pres[jmx];      // Ti
  

  MatrixND<2*NDF,2*NDF> kl;
  kl.zero();

  for (int i=0; i<NDF*2; i++) {
    int ii = std::abs(iq[i]);
    if (ii >= NBV)
      continue;

    for (int j=0; j<NDF*2; j++) {
      int jj = std::abs(iq[j]);
      if (jj >= NBV)
        continue;

      kl(i,j) = kb(ii, jj);
    }
  }

  for (int i = 0; i < 2*NDF; i++) {
    kl(0*NDF+0, i) = kl(i, 0*NDF+0) =  i==0? kl(NDF+0, NDF+0): (i==3? kl(NDF+0, NDF+3) : -kl( NDF+0, i));
    kl(0*NDF+3, i) = kl(i, 0*NDF+3) =  i==0? kl(NDF+3, NDF+0): (i==3? kl(NDF+3, NDF+3) : -kl( NDF+3, i));
    // for (int j=0; j<nwm; j++)
    //   kl(0*NDF+6+j, i) = kl(i, 0*NDF+6+j) =  i==0? kl(NDF+6+j, NDF+0): (i==3? kl(NDF+6+j, NDF+6+j) : -kl( NDF+6+j, i));
  }



  ALWAYS_STATIC MatrixND<2*NDF,2*NDF> Kg;
  ALWAYS_STATIC Matrix Wrapper(Kg);

  Kg = basic_system->t.pushResponse(kl, pl);

  return Wrapper;
}



template <int NIP, int nsr, int nwm>
void
ForceFrame3d<NIP,nsr,nwm>::addLoadAtSection(VectorND<nsr>& sp, double x)
{
  double L = basic_system->getInitialLength();

  for (auto[load, loadFactor] : eleLoads) {

    int type;
    const Vector& data = load->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wa = data(2) * loadFactor; // Axial
      double wy = data(0) * loadFactor; // Transverse
      double wz = data(1) * loadFactor; // Transverse

      for (int ii = 0; ii < nsr; ii++) {
        switch (scheme[ii]) {
        case FrameStress::N:  sp[ii] += wa * (L - x); break;
        case FrameStress::Vy: sp[ii] += wy * (x - 0.5 * L); break;
        case FrameStress::Vz: sp[ii] += wz * (0.5 * L - x); break;
        case FrameStress::Mz: sp[ii] += wy * 0.5 * x * (x - L); break;
        case FrameStress::My: sp[ii] += wz * 0.5 * x * (L - x); break;
        default: break;
        }
      }
    }
  
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * loadFactor;
      double Pz     = data(1) * loadFactor;
      double N      = data(2) * loadFactor;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double a = aOverL * L;

      double Vy1 = Py * (1.0 - aOverL);
      double Vy2 = Py * aOverL;

      double Vz1 = Pz * (1.0 - aOverL);
      double Vz2 = Pz * aOverL;

      for (int ii = 0; ii < nsr; ii++) {

        if (x <= a) {
          switch (scheme[ii]) {
          case FrameStress::N:  sp[ii] +=       N; break;
          case FrameStress::Vy: sp[ii] -=     Vy1; break;
          case FrameStress::Vz: sp[ii] -=     Vz1; break;
          case FrameStress::My: sp[ii] += x * Vz1; break;
          case FrameStress::Mz: sp[ii] -= x * Vy1; break;
          default:                  break;
          }
        } else {
          switch (scheme[ii]) {
          case FrameStress::Mz: sp[ii] -= (L - x) * Vy2; break;
          case FrameStress::Vy: sp[ii] += Vy2; break;
          case FrameStress::My: sp[ii] += (L - x) * Vz2; break;
          case FrameStress::Vz: sp[ii] += Vz2; break;
          default:                  break;
          }
        }
      }
    } 
  
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy  = data(0) * loadFactor;  // Transverse Y at start
      double wz  = data(1) * loadFactor;  // Transverse Z at start
      double wa  = data(2) * loadFactor;  // Axial at start
      double a   = data(3)*L;
      double b   = data(4)*L;
      double wyb = data(5) * loadFactor;  // Transverse Y at end
      double wzb = data(6) * loadFactor;  // Transverse Z at end
      double wab = data(7) * loadFactor;  // Axial at end
      double Fa = wa * (b - a) + 0.5 * (wab - wa) * (b - a); // resultant axial load
      double Fy = wy * (b - a); // resultant transverse load
      double Fz = wz * (b - a); // resultant transverse load
      double c = a + 0.5 * (b - a);
      double VyI = Fy * (1 - c / L);
      double VyJ = Fy * c / L;
      double VzI = Fz * (1 - c / L);
      double VzJ = Fz * c / L;
      Fy = 0.5 * (wyb - wy) * (b - a); // resultant transverse load
      Fz = 0.5 * (wzb - wz) * (b - a); // resultant transverse load
      c = a + 2.0 / 3.0 * (b - a);
      VyI += Fy * (1 - c / L);
      VyJ += Fy * c / L;
      VzI += Fz * (1 - c / L);
      VzJ += Fz * c / L;
     
      for (int ii = 0; ii < nsr; ii++) {
        if (x <= a) {
          switch(scheme[ii]) {
          case FrameStress::N:
            sp(ii) += Fa;
            break;
          case FrameStress::Mz:
            sp(ii) -= VyI*x;
            break;
          case FrameStress::My:
            sp(ii) += VzI*x;
            break;
          case FrameStress::Vy:
            sp(ii) -= VyI;
            break;
          case FrameStress::Vz:
            sp(ii) += VzI;
            break;            
          default:
            break;
          }
        }
        else if (x >= b) {
          switch(scheme[ii]) {
          case FrameStress::Mz:
            sp(ii) += VyJ*(x-L);
            break;
          case FrameStress::My:
            sp(ii) -= VzJ*(x-L);
            break;            
          case FrameStress::Vy:
            sp(ii) += VyJ;
            break;
          case FrameStress::Vz:
            sp(ii) -= VzJ;            
            break;
          default:
            break;
          }
        }
        else {
          double wyy = wy + (wyb - wy) / (b - a) * (x - a);
          double wzz = wz + (wzb - wz) / (b - a) * (x - a);
          switch(scheme[ii]) {
          case FrameStress::N:
            sp(ii) += Fa - wa * (x - a) - 0.5 * (wab - wa) / (b - a) * (x - a) * (x - a);
            break;
          case FrameStress::Mz:
            sp(ii) += -VyI * x + 0.5 * wy * (x - a) * (x - a) + 0.5 * (wyy - wy) * (x - a) * (x - a) / 3.0;
            break;
          case FrameStress::My:
            sp(ii) += VzI * x - 0.5 * wz * (x - a) * (x - a) - 0.5 * (wzz - wz) * (x - a) * (x - a) / 3.0;
            break;            
          case FrameStress::Vy:
            sp(ii) += -VyI + wy * (x - a) + 0.5 * (wyy - wy) * (x - a);
            break;
          case FrameStress::Vz:           
            sp(ii) -= -VzI + wz * (x - a) - 0.5 * (wzz - wz) * (x - a);
            break;
          default:
            break;
          }
        }
      }
    }
    else {
      opserr << "ForceFrame3d::addLoad -- load type unknown for element with tag: "
             << this->getTag() << "\n";
    }
  }
}


template <int NIP, int nsr, int nwm>
void
ForceFrame3d<NIP,nsr,nwm>::getStressGrad(VectorND<nsr>& dspdh, int isec, int gradNumber)
{

  double L    = basic_system->getInitialLength();
  double dLdh = basic_system->getLengthGrad();

  int numSections = points.size();

  double dxidh[NIP];
  stencil->getLocationsDeriv(numSections, L, dLdh, dxidh);

  double x    = points[isec].point * L;
  double dxdh = points[isec].point * dLdh + dxidh[isec] * L;

  for (auto[load, loadFactor] : eleLoads) {
    int type;
    const  Vector& data = load->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0) * 1.0; // Transverse
      double wz = data(1) * 1.0; // Transverse
      double wa = data(2) * 1.0; // Axial

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dwydh       = sens(0);
      double dwzdh       = sens(1);
      double dwadh       = sens(2);

      for (int ii = 0; ii < nsr; ii++) {
        switch (scheme[ii]) {
        case FrameStress::N:
          //sp(ii) += wa*(L-x);
          dspdh(ii) += dwadh * (L - x) + wa * (dLdh - dxdh);
          break;
        case FrameStress::Mz:
          //sp(ii) += wy*0.5*x*(x-L);
          //dspdh(ii) += 0.5*(dwydh*x*(x-L) + wy*dxdh*(x-L) + wy*x*(dxdh-dLdh));
          dspdh(ii) += 0.5 * (dwydh * x * (x - L) + wy * (dxdh * (2 * x - L) - x * dLdh));
          break;
        case FrameStress::Vy:
          //sp(ii) += wy*(x-0.5*L);
          dspdh(ii) += dwydh * (x - 0.5 * L) + wy * (dxdh - 0.5 * dLdh);
          break;
        case FrameStress::My:
          //sp(ii) += wz*0.5*x*(L-x);
          //dspdh(ii) += 0.5*(dwzdh*x*(L-x) + wz*dxdh*(L-x) + wz*x*(dLdh-dxdh));
          dspdh(ii) += 0.5 * (dwzdh * x * (L - x) + wz * (dxdh * (L - 2 * x) + x * dLdh));
          break;
        case FrameStress::Vz:
          //sp(ii) += wz*(x-0.5*L);
          dspdh(ii) += dwzdh * (0.5 * L - x) + wz * (0.5 * dLdh - dxdh);
          break;
        default: break;
        }
      }
    } else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * 1.0;
      double Pz     = data(1) * 1.0;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dPydh       = sens[0];
      double dPzdh       = sens[1];
      double dNdh        = sens[2];
      double daLdh       = sens[3];

      double a = aOverL * L;

      double Vy1    = Py * (1.0 - aOverL);
      double Vy2    = Py * aOverL;
      double dVy1dh = Py * (0.0 - daLdh) + dPydh * (1.0 - aOverL);
      double dVy2dh = Py * daLdh + dPydh * aOverL;

      double Vz1    = Pz * (1.0 - aOverL);
      double Vz2    = Pz * aOverL;
      double dVz1dh = Pz * (0.0 - daLdh) + dPzdh * (1.0 - aOverL);
      double dVz2dh = Pz * daLdh + dPzdh * aOverL;

      for (int ii = 0; ii < nsr; ii++) {

        if (x <= a) {
          switch (scheme[ii]) {
          case FrameStress::N:
            //sp(ii) += N;
            dspdh(ii) += dNdh;
            break;
          case FrameStress::Mz:
            //sp(ii) -= x*Vy1;
            dspdh(ii) -= (dxdh * Vy1 + x * dVy1dh);
            break;
          case FrameStress::Vy:
            //sp(ii) -= Vy1;
            dspdh(ii) -= dVy1dh;
            break;
          case FrameStress::My:
            //sp(ii) += x*Vz1;
            dspdh(ii) += (dxdh * Vz1 + x * dVz1dh);
            break;
          case FrameStress::Vz:
            //sp(ii) -= Vz1;
            dspdh(ii) -= dVz1dh;
            break;
          default: break;
          }
        } else {
          switch (scheme[ii]) {
          case FrameStress::Vy:
            //sp(ii) += Vy2;
            dspdh(ii) += dVy2dh;
            break;
          case FrameStress::Vz:
            //sp(ii) += Vz2;
            dspdh(ii) += dVz2dh;
            break;
          case FrameStress::My:
            //sp(ii) += (L-x)*Vz2;
            dspdh(ii) += (dLdh - dxdh) * Vz2 + (L - x) * dVz2dh;
            break;
          case FrameStress::Mz:
            //sp(ii) -= (L-x)*Vy2;
            dspdh(ii) -= (dLdh - dxdh) * Vy2 + (L - x) * dVy2dh;
            break;
          default: break;
          }
        }
      }
    } else {
      opserr << "ForceFrame3d::getStressGrad -- load type unknown for element "
                "with tag: "
             << this->getTag() 
             << "\n";
    }
  }
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::sendSelf(int commitTag, Channel& theChannel)
{
  return -1;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  return -1;
}


// addBFsB(F, s, i)
template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::getInitialFlexibility(MatrixND<NBV,NBV>& fe)
{
  fe.zero();

  double L        = basic_system->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  for (int i = 0; i < numSections; i++) {

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    const MatrixND<nsr,nsr> fSec = points[i].material->template getFlexibility<nsr, scheme>(State::Init);

    THREAD_LOCAL MatrixND<nsr, NBV> FB;
    FB.zero();
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case FrameStress::N:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 0) += fSec(jj, ii) * wtL;
        break;
      case FrameStress::Mz:
        for (int jj = 0; jj < nsr; jj++) {
          double tmp = fSec(jj, ii) * wtL;
          FB(jj, 1) += xL1 * tmp;
          FB(jj, 2) += xL * tmp;
        }
        break;
      case FrameStress::Vy:
        for (int jj = 0; jj < nsr; jj++) {
          double tmp = jsx * fSec(jj, ii) * wtL;
          FB(jj, 1) += tmp;
          FB(jj, 2) += tmp;
        }
        break;
      case FrameStress::My:
        for (int jj = 0; jj < nsr; jj++) {
          double tmp = fSec(jj, ii) * wtL;
          FB(jj, 3) += xL1 * tmp;
          FB(jj, 4) += xL * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int jj = 0; jj < nsr; jj++) {
          double tmp = jsx * fSec(jj, ii) * wtL;
          FB(jj, 3) += tmp;
          FB(jj, 4) += tmp;
        }
        break;
      case FrameStress::T:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 5) += fSec(jj, ii) * wtL;
        break;
      default: break;
      }
    }
    //
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case FrameStress::N:
        for (int jj = 0; jj < NBV; jj++)
          fe(0, jj) += FB(ii, jj);
        break;
      case FrameStress::Mz:
        for (int jj = 0; jj < NBV; jj++) {
          double tmp = FB(ii, jj);
          fe(1, jj) += xL1 * tmp;
          fe(2, jj) += xL * tmp;
        }
        break;
      case FrameStress::Vy:
        for (int jj = 0; jj < NBV; jj++) {
          double tmp = jsx * FB(ii, jj);
          fe(1, jj) += tmp;
          fe(2, jj) += tmp;
        }
        break;
      case FrameStress::My:
        for (int jj = 0; jj < NBV; jj++) {
          double tmp = FB(ii, jj);
          fe(3, jj) += xL1 * tmp;
          fe(4, jj) += xL * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int jj = 0; jj < NBV; jj++) {
          double tmp = jsx * FB(ii, jj);
          fe(3, jj) += tmp;
          fe(4, jj) += tmp;
        }
        break;
      case FrameStress::T:
        for (int jj = 0; jj < NBV; jj++)
          fe(5, jj) += FB(ii, jj);
        break;
      default: break;
      }
    }
  }
  return 0;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::getInitialDeformations(Vector& v0)
{
  v0.Zero();
  if (eleLoads.size() < 1 || (this->setState(State::Init) != 0))
    return 0;

  double L = basic_system->getInitialLength();
  double jsx = 1.0 / L;


  const int numSections = points.size();
  for (int i = 0; i < numSections; i++) {

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    THREAD_LOCAL VectorND<nsr> sp;
    sp.zero();

    this->addLoadAtSection(sp, points[i].point*L);

    MatrixND<nsr,nsr> fse = points[i].material->template getFlexibility<nsr,scheme>(State::Init);

    VectorND<nsr> e = fse*sp;

    for (int ii = 0; ii < nsr; ii++) {
      double tmp;
      double dei = e[ii] * wtL;
      switch (scheme[ii]) {
      case FrameStress::N:
        v0[0] += dei;
        break;
      case FrameStress::Vy:
        tmp = jsx * dei;
        v0[1] += tmp;
        v0[2] += tmp;
        break;
      case FrameStress::Vz:
        tmp = jsx * dei;
        v0[3] += tmp;
        v0[4] += tmp;
        break;
      case FrameStress::My:
        v0[3] += xL1 * dei;
        v0[4] += xL  * dei;
        break;
      case FrameStress::Mz:
        v0[1] += xL1 * dei;
        v0[2] += xL  * dei;
        break;
      default: break;
      }
    }
  }

  return 0;
}

template <int NIP, int nsr, int nwm>
void
ForceFrame3d<NIP,nsr,nwm>::Print(OPS_Stream& s, int flag)
{
  const int nip = points.size();
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() 
                        << "<" << NIP << ", " << nsr << ", " << nwm <<  ">" << "\"";
    s << ", ";

    s << "\"nodes\": [" << node_tags(0) << ", " 
                        << node_tags(1) << "]";
    s << ", ";

    // Mass
    double mass;
    if (getIntegral(Field::Density, State::Init, mass) == 0)
      s << "\"mass\": " << mass;
    else
      s << "\"massperlength\": " << density;
    s << ", ";
    s << "\"mass_flag\": "<< mass_flag;
    s << ", ";

    // Integration points
    s << "\"sections\": [";
    for (decltype(points.size()) i = 0; i < points.size() - 1; ++i)
      s << points[i].material->getTag() << ", ";
    s << points[points.size() - 1].material->getTag() 
      << "]";
    s << ", ";

    s << "\"transform\": " << basic_system->getTag();
    s << ", ";

    s << "\"integration\": ";
    stencil->Print(s, flag);
    s << "}";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: ForceFrame3d ";
    s << "\tNumber of Sections: " << nip;
    s << "\tMass density: " << density << "\n";
    stencil->Print(s, flag);
    double P     = q_save[0];
    double MZ1   = q_save[1];
    double MZ2   = q_save[2];
    double MY1   = q_save[3];
    double MY2   = q_save[4];
    double L     = basic_system->getInitialLength();
    double VY    = (MZ1 + MZ2) / L;
    double VZ    = (MY1 + MY2) / L;
    double T     = q_save[5];

    double p0[5]{};
    if (eleLoads.size() > 0)
      this->computeReactions(p0);

    s << "\tEnd 1 Forces (P MZ VY MY VZ T): " << -P + p0[0] << " " << MZ1 << " " << VY + p0[1]
      << " " << MY1 << " " << -VZ + p0[3] << " " << T << "\n";
    s << "\tEnd 2 Forces (P MZ VY MY VZ T): " << P << " " << MZ2 << " " << -VY + p0[2] << " " << MY2
      << " " << VZ + p0[4] << " " << -T << "\n";

    for (int i = 0; i < nip; i++)
      points[i].material->Print(s, flag);
  }
}


template <int NIP, int nsr, int nwm>
Response*
ForceFrame3d<NIP,nsr,nwm>::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = nullptr;
  const ID& node_tags = this->getExternalNodes();

  output.tag("ElementOutput");
  output.attr("eleType", this->getClassType());
  output.attr("eleTag", this->getTag());
  output.attr("node1",  node_tags(0));
  output.attr("node2",  node_tags(1));

  //
  // compare argv[0] for known response types
  //

  // Global forces
  if (strcmp(argv[0], "forces") == 0 || 
      strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 ||
      strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, Respond::GlobalForce, Vector(NEN*NDF));

  // Local force
  } else if (strcmp(argv[0], "localForce") == 0 || 
             strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, Respond::LocalForce, Vector(NEN*NDF));

  // Basic force
  } else if (strcmp(argv[0], "basicForce") == 0 || 
             strcmp(argv[0], "basicForces") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, Respond::BasicForce, Vector(NBV));

  // basic stiffness
  } else if (strcmp(argv[0], "basicStiffness") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, Respond::BasicStiff, Matrix(6, 6));

  // 3: Chord rotation
  } else if (strcmp(argv[0], "chordRotation") == 0 || 
             strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

  // 4: Plastic rotation
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, Respond::BasicPlasticDeformation, Vector(6));

  // 5: Point of inflection
  } else if (strcmp(argv[0], "inflectionPoint") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));

  // 6: Tangent drift
  } else if (strcmp(argv[0], "tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(4));

  } else if (strcmp(argv[0], "RayleighForces") == 0 || 
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, Vector(NEN*NDF));

  } else if (strcmp(argv[0], "sections") == 0) {
    if (this->setState(State::Init) != 0)
      return nullptr;

    CompositeResponse* theCResponse = new CompositeResponse();
    int numResponse                 = 0;
    const int numSections = points.size();
    double L = basic_system->getInitialLength();

    for (int i = 0; i < numSections; i++) {
      output.tag("GaussPointOutput");
      output.attr("number", i + 1);
      output.attr("eta", points[i].point * L);

      Response* theSectionResponse = points[i].material->setResponse(&argv[1], argc - 1, output);

      if (theSectionResponse != nullptr)
        numResponse = theCResponse->addResponse(theSectionResponse);
    }

    if (numResponse == 0) // no valid responses found
      delete theCResponse;
    else
      theResponse = theCResponse;
  }

  // 10-11: Integration
  else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(points.size()));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(points.size()));

  // 110-111: sections
  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(points.size()));

  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    opserr << "sectionDisplacements is not implemented\n";
    return nullptr;
  }

  else if (strcmp(argv[0], "cbdiDisplacements") == 0) {
    opserr << "cbdiDisplacements is not implemented\n";
    // theResponse = new ElementResponse(this, 112, Matrix(20, 3));
    return nullptr;
  }


  else if (strstr(argv[0], "section") != nullptr) {

    if (argc > 1) {

      const int nip = points.size();
      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= nip && argc > 2) {
        if (this->setState(State::Init) != 0)
          return nullptr;

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", 2.0 * points[sectionNum - 1].point - 1.0);

        if (strcmp(argv[2], "dsdh") != 0) {
          theResponse = points[sectionNum - 1].material->setResponse(&argv[2], argc - 2, output);
        } else {
          int order         = points[sectionNum - 1].material->getOrder();
          theResponse       = new ElementResponse(this, 76, Vector(order));
          Information& info = theResponse->getInformation();
          info.theInt       = sectionNum;
        }

        output.endTag();

      } else if (sectionNum == 0) { 
        // argv[1] was not an int, we want all sections,

        CompositeResponse* theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[NIP];
        double L = basic_system->getInitialLength();
        const int numSections = points.size();
        stencil->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", points[i].point * L);

          Response* theSectionResponse = points[i].material->setResponse(&argv[1], argc - 1, output);

          if (theSectionResponse != nullptr) {
            numResponse = theCResponse->addResponse(theSectionResponse);
          }
        }

        if (numResponse == 0) // no valid responses found
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }
  //by SAJalali
  else if (strcmp(argv[0], "energy") == 0) {
    return new ElementResponse(this, 2000, 0.0);
  }

  if (theResponse == nullptr) {
    theResponse = basic_system->setResponse(argv, argc, output);
  }

  output.endTag();

  return theResponse;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::getResponse(int responseID, Information& info)
{
  THREAD_LOCAL Vector vp(6);
  THREAD_LOCAL MatrixND<NBV,NBV> fe;

  if (responseID == Respond::GlobalForce)
    return info.setVector(this->getResistingForce());

  else if (responseID == Respond::LocalForce) {
    // TODO: Warping DOF
    THREAD_LOCAL VectorND<NEN*NDF> v_resp{};
    THREAD_LOCAL Vector v_wrap(v_resp);

    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);
    // Axial
    double N  = q_pres[0];
    v_resp(6) =  N;
    v_resp(0) = -N + p0[0];

    // Torsion
    double T  = q_pres[5];
    v_resp(9) =  T;
    v_resp(3) = -T;

    // Moments about z and shears along y
    double M1     = q_pres[imz];
    double M2     = q_pres[jmz];
    v_resp(5)  = M1;
    v_resp(11) = M2;
    double L      = basic_system->getInitialLength();
    double V      = (M1 + M2) / L;
    v_resp(1)  =  V + p0[1];
    v_resp(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1         = q_pres[imy];
    M2         = q_pres[jmy];
    v_resp(4)  = M1;
    v_resp(10) = M2;
    V          = (M1 + M2) / L;
    v_resp(2)  = -V + p0[3];
    v_resp(8)  =  V + p0[4];

    return info.setVector(v_wrap);
  }

  // Chord rotation
  else if (responseID == 3) {
    vp = basic_system->getBasicTrialDisp();
    return info.setVector(vp);
  }

  else if (responseID == Respond::BasicForce)
    return info.setVector(q_pres);

  else if (responseID == 19)
    return info.setMatrix(K_pres);

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = basic_system->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, q_pres, -1.0);
    Vector v0(6);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return info.setVector(vp);
  }

  else if (responseID == 10) {
    // "integrationPoints" response
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    double L = basic_system->getInitialLength();
    const int nip = points.size();
    Vector locs(nip);
    for (int i = 0; i < nip; i++)
      locs[i] = points[i].point * L;

    return info.setVector(locs);
  }

  else if (responseID == 11) {
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    double L = basic_system->getInitialLength();

    const int nip = points.size();
    Vector weights(nip);
    for (int i = 0; i < nip; i++)
      weights(i) = points[i].weight * L;

    return info.setVector(weights);
  }

  else if (responseID == 110) {
    const int numSections = points.size();
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = points[i].material->getTag();
    return info.setID(tags);
  }

  else if (responseID == 111 || responseID == 1111) {
  }

  else if (responseID == 112) {
  }

  else if (responseID == 12)
    return info.setVector(this->getRayleighDampingForces());

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = basic_system->getInitialLength();

    if (fabs(q_pres[imz] + q_pres[jmz]) > DBL_EPSILON)
      LI(0) = q_pres[imz] / (q_pres[imz] + q_pres[jmz]) * L;

    if (fabs(q_pres[imy] + q_pres[jmy]) > DBL_EPSILON)
      LI(1) = q_pres[imy] / (q_pres[imy] + q_pres[jmy]) * L;

    return info.setVector(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = basic_system->getInitialLength();

    double wts[NIP];
    const int numSections = points.size();
    stencil->getSectionWeights(numSections, L, wts);

    double pts[NIP];
    stencil->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(q_pres[imz] + q_pres[jmz]) > DBL_EPSILON)
      LIz = q_pres[imz] / (q_pres[imz] + q_pres[jmz]) * L;

    double LIy = 0.0;
    if (fabs(q_pres[imy] + q_pres[jmy]) > DBL_EPSILON)
      LIy = q_pres[imy] / (q_pres[imy] + q_pres[jmy]) * L;

    for (int i = 0; i < numSections; i++) {
      double x       = points[i].point * L;
      const ID& type = points[i].material->getType();
      double kappa   = 0.0;
      if (x < LIz) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == FrameStress::Mz)
            kappa += points[i].es[j];
        double b = -LIz + x;
        d2z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x < LIy) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == FrameStress::My)
            kappa += points[i].es[j];
        double b = -LIy + x;
        d2y += (wts[i] * L) * kappa * b;
      }
    }

    d2z += stencil->getTangentDriftI(L, LIz, q_pres[imz], q_pres[jmz]);
    d2y += stencil->getTangentDriftI(L, LIy, q_pres[imy], q_pres[jmy], true);

    for (int i = numSections - 1; i >= 0; i--) {
      double x       = pts[i] * L;
      double kappa   = 0.0;
      if (x > LIz) {
        for (int j = 0; j < nsr; j++)
          if (scheme[j] == FrameStress::Mz)
            kappa += points[i].es[j];
        double b = x - LIz;
        d3z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x > LIy) {
        for (int j = 0; j < nsr; j++)
          if (scheme[j] == FrameStress::My)
            kappa += points[i].es[j];
        double b = x - LIy;
        d3y += (wts[i] * L) * kappa * b;
      }
    }

    d3z += stencil->getTangentDriftJ(L, LIz, q_pres[imz], q_pres[jmz]);
    d3y += stencil->getTangentDriftJ(L, LIy, q_pres[imy], q_pres[jmy], true);

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;
    return info.setVector(d);
  }

  else if (responseID == 2000) {
    const int numSections = points.size();
    double xi[NIP];
    double L = basic_system->getInitialLength();
    stencil->getSectionWeights(numSections, L, xi);
    double energy = 0;
    for (int i = 0; i < numSections; i++) {
      energy += points[i].material->getEnergy() * points[i].point * L;
    }
    return info.setDouble(energy);
  }

  return -1;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::getResponseSensitivity(int responseID, int gradNumber, Information& info)
{
  // Basic deformation sensitivity
  if (responseID == 3) {
    const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);
    return info.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(6);
    dqdh.Zero();

    const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);

    dqdh.addMatrixVector(0.0, K_pres, dvdh, 1.0);

    dqdh.addVector(1.0, this->getBasicForceGrad(gradNumber), 1.0);

    return info.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = info.theInt;

    VectorND<nsr> dsdh;
    dsdh.zero();

    if (eleLoads.size() > 0)
      this->getStressGrad(dsdh, sectionNum - 1, gradNumber);

    static Vector dqdh(6);
    {
      // Response 7
      const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);

      dqdh.addMatrixVector(0.0, K_pres, dvdh, 1.0);

      dqdh.addVector(1.0, this->getBasicForceGrad(gradNumber), 1.0);
    }

    const int numSections = points.size();
    double L   = basic_system->getInitialLength();
    double jsx = 1.0 / L;
    double pts[NIP];
    stencil->getSectionLocations(numSections, L, pts);

    const ID& code = points[sectionNum - 1].material->getType();

    double xL  = pts[sectionNum - 1];
    double xL1 = xL - 1.0;

    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case FrameStress::N:  dsdh(ii) += dqdh(jnx); break;
      case FrameStress::Mz: dsdh(ii) += xL1 * dqdh(imz) + xL * dqdh(jmz); break;
      case FrameStress::Vy: dsdh(ii) += jsx * (dqdh(1) + dqdh(2)); break;
      case FrameStress::My: dsdh(ii) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case FrameStress::Vz: dsdh(ii) += jsx * (dqdh(3) + dqdh(4)); break;
      case FrameStress::T:  dsdh(ii) += dqdh(jmx); break;
      default:              dsdh(ii) += 0.0; break;
      }
    }

    double dLdh   = basic_system->getLengthGrad();
    double d1oLdh = basic_system->getd1overLdh();

    double dptsdh[NIP];
    stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);
    double dxLdh = dptsdh[sectionNum - 1]; // - xL/L*dLdh;

    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case FrameStress::Mz:
        dsdh(j) += dxLdh * (q_pres[imz] + q_pres[jmz]);
        //dsdh(j) -= dLdh*xL/L*(Se(1)+Se(2));
        break;
      case FrameStress::Vy: dsdh(j) += d1oLdh * (q_pres[imz] + q_pres[jmz]); break;
      case FrameStress::My: dsdh(j) += dxLdh  * (q_pres[imy] + q_pres[jmy]); break;
      case FrameStress::Vz: dsdh(j) += d1oLdh * (q_pres[imy] + q_pres[jmy]); break;
      default:                  break;
      }
    }
    return info.setVector(dsdh);
  }

  // Plastic deformation sensitivity
  else if (responseID == 4) {
    static Vector dvpdh(6);

    const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);

    dvpdh = dvdh;

    static MatrixND<NBV,NBV> fe;
    this->getInitialFlexibility(fe);

    const Vector& dqdh = this->getBasicForceGrad(gradNumber);

    dvpdh.addMatrixVector(1.0, fe, dqdh, -1.0);

    dvpdh.addMatrixVector(1.0, fe*K_pres, dvdh, -1.0);

    const Matrix& dfedh = this->computedfedh(gradNumber);

    dvpdh.addMatrixVector(1.0, dfedh, q_pres, -1.0);

    return info.setVector(dvpdh);
  }

  else
    return -1;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // If the parameter belongs to the element itself
  if ((strcmp(argv[0], "rho") == 0) ||
      (strcmp(argv[0], "density") == 0)) {
    param.setValue(density);
    return param.addObject(1, this);
  }

  // Section response
  if (strstr(argv[0], "sectionX") != nullptr) {
    if (argc > 2) {
      double sectionLoc = atof(argv[1]);

      const int numSections = points.size();
      double xi[NIP];
      double L = basic_system->getInitialLength();
      stencil->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(points[i].point - sectionLoc) < minDistance) {
          minDistance = fabs(points[i].point - sectionLoc);
          sectionNum  = i;
        }
      }

      return points[sectionNum].material->setParameter(&argv[2], argc - 2, param);
    }
  }

  // If the parameter belongs to a particular section or lower
  if (strstr(argv[0], "section") != nullptr) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= (int)points.size())
      return points[sectionNum - 1].material->setParameter(&argv[2], argc - 2, param);

    else
      return -1;
  }

  // If the parameter belongs to all sections or lower
  if (strstr(argv[0], "allSections") != nullptr) {

    if (argc < 2)
      return -1;

    for (GaussPoint& point : points) {
      int ok = point.material->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }

    return result;
  }

  if (strstr(argv[0], "integration") != nullptr) {

    if (argc < 2)
      return -1;

    return stencil->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to everything

  for (GaussPoint& point : points) {
    int ok = point.material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  int ok = stencil->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::updateParameter(int parameterID, Information& info)
{
  if (parameterID == 1) {
    this->density = info.theDouble;
    return 0;
  } else
    return -1;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}


template <int NIP, int nsr, int nwm>
const Matrix&
ForceFrame3d<NIP,nsr,nwm>::getKiSensitivity(int gradNumber)
{
  ALWAYS_STATIC MatrixND<NEN*NDF,NEN*NDF> Ksen{};
  ALWAYS_STATIC Matrix wrapper(Ksen);
  return wrapper;
}


template <int NIP, int nsr, int nwm>
const Matrix&
ForceFrame3d<NIP,nsr,nwm>::getMassSensitivity(int gradNumber)
{
  ALWAYS_STATIC MatrixND<NEN*NDF,NEN*NDF> Msen{};
  ALWAYS_STATIC Matrix wrapper(Msen);

  double L = basic_system->getInitialLength();
  if (parameterID == 1) {
    // TODO: handle consistent mass
    if (density != 0.0) {
      Msen(0, 0) = Msen(1, 1) = Msen(2, 2) = 
      Msen(6, 6) = Msen(7, 7) = Msen(8, 8) = 0.5 * L;
    }
  } else
    Msen.zero();


  return wrapper;
}


template <int NIP, int nsr, int nwm>
const Vector&
ForceFrame3d<NIP,nsr,nwm>::getResistingForceSensitivity(int gradNumber)
{
  static Vector P(NEN*NDF);
  P.Zero();

  VectorND<NBV> dqdh = this->getBasicForceGrad(gradNumber);

  // Transform forces
  double dp0dh[6]{};
  this->addReactionGrad(dp0dh, gradNumber, basic_system->getLengthGrad());
  Vector dp0dhVec(dp0dh, 6);

  if (basic_system->isShapeSensitivity()) {
    //
    // dqdh += K dvdh|_ug
    //
    dqdh.addMatrixVector(1.0, K_pres, 
                         basic_system->getBasicDisplFixedGrad(), 1.0);

    // dAdh^T q
    P = basic_system->getGlobalResistingForceShapeSensitivity(q_pres, dp0dhVec, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += basic_system->getGlobalResistingForce(dqdh, dp0dhVec);

  return P;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::commitSensitivity(int gradNumber, int numGrads)
{
  //
  //
  //

  int err = 0;

  double L   = basic_system->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();

  double pts[NIP];
  double wts[NIP];
  stencil->getSectionLocations(numSections, L, pts);
  stencil->getSectionWeights(numSections, L, wts);

  double dLdh = basic_system->getLengthGrad();

  double dptsdh[NIP];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double d1oLdh = basic_system->getd1overLdh();

  // Response 7
  VectorND<NBV> dqdh = this->getBasicForceGrad(gradNumber);

  // dvdh = A u` + A` u
  const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);
  dqdh.addMatrixVector(1.0, K_pres, dvdh, 1.0);


  //
  // Compute and commit de for each section
  //
  for (int i = 0; i < numSections; i++) {

    double xL  = pts[i];
    double xL1 = xL - 1.0;

    double dxLdh = dptsdh[i];

    VectorND<nsr> ds;
    ds.zero();

    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0)
      this->getStressGrad(ds, i, gradNumber);


    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::N:  ds(j) += dqdh(0); break;
      case FrameStress::Vy: ds(j) += jsx * (dqdh(1) + dqdh(2)); break;
      case FrameStress::Vz: ds(j) += jsx * (dqdh(3) + dqdh(4)); break;
      case FrameStress::T:  ds(j) += dqdh(5); break;
      case FrameStress::My: ds(j) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case FrameStress::Mz: ds(j) += xL1 * dqdh(1) + xL * dqdh(2); break;
      default:                  ds(j) += 0.0; break;
      }
    }

    const Vector& dsdh = points[i].material->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::Mz: ds(j) += dxLdh  * (q_pres[imz] + q_pres[jmz]); break;
      case FrameStress::Vy: ds(j) += d1oLdh * (q_pres[imz] + q_pres[jmz]); break;
      case FrameStress::My: ds(j) += dxLdh  * (q_pres[imy] + q_pres[jmy]); break;
      case FrameStress::Vz: ds(j) += d1oLdh * (q_pres[imy] + q_pres[jmy]); break;
      default: break;
      }
    }

    Vector de(nsr);
    const Matrix& fs = points[i].material->template getFlexibility<nsr,scheme>();
    de.addMatrixVector(0.0, fs, ds, 1.0);

    err += points[i].material->commitSensitivity(de, gradNumber, numGrads);
  }

  return err;
}


template <int NIP, int nsr, int nwm>
VectorND<6+nwm*2>
ForceFrame3d<NIP,nsr,nwm>::getBasicForceGrad(int gradNumber)
{

  double L   = basic_system->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  double wts[NIP];
  stencil->getSectionWeights(numSections, L, wts);

  double dLdh = basic_system->getLengthGrad();

  double dptsdh[NIP];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[NIP];
  stencil->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  double d1oLdh = basic_system->getd1overLdh();

  //
  // Integrate dvdh
  //
  static Vector dvdh(NBV);
  dvdh.Zero();
  for (int i = 0; i < numSections; i++) {
    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    double dxLdh  = dptsdh[i]; // - xL/L*dLdh;
    double dwtLdh = points[i].weight * dLdh + dwtsdh[i] * L;


    // Get section stress resultant gradient

    VectorND<nsr> dspdh{};
    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0)
      this->getStressGrad(dspdh, i, gradNumber);

    Vector dsdh(nsr);
    dsdh = points[i].material->getStressResultantSensitivity(gradNumber, true);

    dsdh.addVector(1.0, dspdh, -1.0);

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::Mz: dsdh(j) -= dxLdh  * (q_pres[imz] + q_pres[jmz]); break;
      case FrameStress::Vy: dsdh(j) -= d1oLdh * (q_pres[imz] + q_pres[jmz]); break;
      case FrameStress::My: dsdh(j) -= dxLdh  * (q_pres[imy] + q_pres[jmy]); break;
      case FrameStress::Vz: dsdh(j) -= d1oLdh * (q_pres[imy] + q_pres[jmy]); break;
      default: break;
      }
    }

    Vector dedh(nsr);
    const Matrix& fs = points[i].material->getSectionFlexibility();
    dedh.addMatrixVector(0.0, fs, dsdh, 1.0);

    for (int j = 0; j < nsr; j++) {
      double dei = dedh(j) * wtL;
      switch (scheme[j]) {
      case FrameStress::N:
        dvdh(0) += dei;
        break;
      case FrameStress::Mz:
        dvdh(1) += xL1 * dei;
        dvdh(2) += xL * dei;
        break;
      case FrameStress::Vy:
        dei = jsx * dei;
        dvdh(1) += dei;
        dvdh(2) += dei;
        break;
      case FrameStress::My:
        dvdh(3) += xL1 * dei;
        dvdh(4) += xL * dei;
        break;
      case FrameStress::Vz:
        dei = jsx * dei;
        dvdh(3) += dei;
        dvdh(4) += dei;
        break;
      case FrameStress::T:
        dvdh(jmx) += dei;
        break;
      default:
        break;
      }
    }

    const VectorND<nsr>& e = points[i].es;
    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::N:
        dvdh(0) -= e(j) * dwtLdh; 
        break;
      case FrameStress::Vy:
        dvdh(1) -= jsx * e(j) * dwtLdh;
        dvdh(2) -= jsx * e(j) * dwtLdh;

        dvdh(1) -= d1oLdh * e(j) * wtL;
        dvdh(2) -= d1oLdh * e(j) * wtL;
          break;
      case FrameStress::Vz:
        dvdh(3) -= jsx * e(j) * dwtLdh;
        dvdh(4) -= jsx * e(j) * dwtLdh;

        dvdh(3) -= d1oLdh * e(j) * wtL;
        dvdh(4) -= d1oLdh * e(j) * wtL;
        break;
      case FrameStress::T:
        dvdh(5) -= e(j) * dwtLdh; 
        break;
      case FrameStress::Mz:
        dvdh(1) -= xL1 * e(j) * dwtLdh;
        dvdh(2) -= xL * e(j) * dwtLdh;

        dvdh(1) -= dxLdh * e(j) * wtL;
        dvdh(2) -= dxLdh * e(j) * wtL;
        break;
      case FrameStress::My:
        dvdh(3) -= xL1 * e(j) * dwtLdh;
        dvdh(4) -= xL * e(j) * dwtLdh;

        dvdh(3) -= dxLdh * e(j) * wtL;
        dvdh(4) -= dxLdh * e(j) * wtL;
        break;

      default:
        break;
      }
    }
  }

  return K_pres*dvdh;
}


template <int NIP, int nsr, int nwm>
const Matrix&
ForceFrame3d<NIP,nsr,nwm>::computedfedh(int gradNumber)
{
  static Matrix dfedh(NBV, NBV);

  dfedh.Zero();
#if 0
  double L   = basic_system->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  double dLdh   = basic_system->getLengthGrad();
  double d1oLdh = basic_system->getd1overLdh();

  double dptsdh[NIP];
  double dwtsdh[NIP];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);
  stencil->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  for (int i = 0; i < numSections; i++) {

    Matrix fb(nsr, NBV);
    Matrix fb2(nsr, NBV);

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    double dxLdh  = dptsdh[i];
    double dwtLdh = points[i].weight * dLdh + dwtsdh[i] * L;

    const Matrix& fs    = points[i].material->getInitialFlexibility();
    const Matrix& dfsdh = points[i].material->getInitialFlexibilitySensitivity(gradNumber);
    fb.Zero();
    fb2.Zero();

    double tmp;
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case FrameStress::N:
        for (int jj = 0; jj < nsr; jj++) {
          fb(jj, 0) += dfsdh(jj, ii) * wtL; // 1

          //fb(jj,0) += fs(jj,ii)*dwtLdh; // 3

          //fb2(jj,0) += fs(jj,ii)*wtL; // 4
        }
        break;
      case FrameStress::Mz:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = dfsdh(jj, ii) * wtL; // 1
          fb(jj, 1) += xL1 * tmp;
          fb(jj, 2) += xL * tmp;

          tmp = fs(jj, ii) * wtL; // 2
          //fb(jj,1) += dxLdh*tmp;
          //fb(jj,2) += dxLdh*tmp;

          tmp = fs(jj, ii) * dwtLdh; // 3
          //fb(jj,1) += xL1*tmp;
          //fb(jj,2) += xL*tmp;

          tmp = fs(jj, ii) * wtL; // 4
          //fb2(jj,1) += xL1*tmp;
          //fb2(jj,2) += xL*tmp;
        }
        break;
      case FrameStress::Vy:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = jsx * dfsdh(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
          // TODO: Need to complete for dLdh != 0
        }
        break;
      default:
        break;
      }
    }
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case FrameStress::N:
        for (int jj = 0; jj < NBV; jj++)
          dfedh(0, jj) += fb(ii, jj);
        break;
      case FrameStress::Mz:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = fb(ii, jj); // 1,2,3
          dfedh(1, jj) += xL1 * tmp;
          dfedh(2, jj) += xL * tmp;

          tmp = fb2(ii, jj); // 4
          //dfedh(1,jj) += dxLdh*tmp;
          //dfedh(2,jj) += dxLdh*tmp;
        }
        break;
      case FrameStress::Vy:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = jsx * fb(ii, jj);
          dfedh(1, jj) += tmp;
          dfedh(2, jj) += tmp;
          // TODO: Need to complete for dLdh != 0
        }
        break;
      default:
        break;
      }
    }
  }
#endif
  return dfedh;
}


template <int NIP, int nsr, int nwm>
int
ForceFrame3d<NIP,nsr,nwm>::setSectionPointers(std::vector<FrameSection*>& new_sections)
{
  // Return value of 0 indicates success

  points.clear();

  for (FrameSection* section : new_sections) {
    assert(section != nullptr);

    points.push_back({
        .point=0,
        .weight=0,
        .material=section->getFrameCopy(scheme)
    });
  }
  return 0;
}


template <int NIP, int nsr, int nwm>
const Vector &
ForceFrame3d<NIP,nsr,nwm>::getResistingForce()
{
  double p0[5]{};
  if (eleLoads.size() > 0)
    this->computeReactions(p0);

  VectorND<NDF*2> pl{};
  pl[0*NDF+0]  = -q_pres[jnx];               // Ni
  #ifdef DO_BASIC
    double L = basic_system->getInitialLength();
    const double q1 = q_pres[imz],
                 q2 = q_pres[jmz],
                 q3 = q_pres[imy],
                 q4 = q_pres[jmy];
    pl[0*NDF+1]  =  (q1 + q2)/L;      // Viy
    pl[0*NDF+2]  = -(q3 + q4)/L;      // Viz
    pl[1*NDF+1]  = -pl[1];            // Vjy
    pl[1*NDF+2]  = -pl[2];            // Vjz
  #endif
    pl[0*NDF+3]  = -q_pres[jmx];      // Ti
    pl[0*NDF+4]  =  q_pres[imy];
    pl[0*NDF+5]  =  q_pres[imz];
    pl[1*NDF+0]  =  q_pres[jnx];      // Nj
    pl[1*NDF+3]  =  q_pres[jmx];      // Tj
    pl[1*NDF+4]  =  q_pres[jmy];
    pl[1*NDF+5]  =  q_pres[jmz];

  thread_local VectorND<NDF*2> pf;
  pf.zero();
  pf[0*NDF + 0] = p0[0];
  pf[0*NDF + 1] = p0[1];
  pf[0*NDF + 2] = p0[3];
  pf[1*NDF + 1] = p0[2];
  pf[1*NDF + 2] = p0[4];

  thread_local VectorND<NDF*2> pg;
  thread_local Vector wrapper(pg);

  pg  = basic_system->t.pushResponse(pl);
  pg += basic_system->t.pushConstant(pf);
  if (total_mass != 0.0)
    wrapper.addVector(1.0, this->FiniteElement<2,3,6+nwm>::p_iner, -1.0);

  return wrapper;
}



#if 0
void
ForceFrame3d<NIP,nsr,nwm>::getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code)
{
  bp.Zero();

  double L = basic_system->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case FrameStress::Mz: // Moment, Mz, interpolation
      bp(i, 1) = xi * (xi - 1) * L * L / 2;
      break;
    case FrameStress::N: // Axial, P, interpolation
      bp(i, 0) = (1 - xi) * L;
      break;
    case FrameStress::Vy: // Shear, Vy, interpolation
      bp(i, 1) = (xi - 0.5) * L;
      break;
    case FrameStress::My: // Moment, My, interpolation
      bp(i, 2) = xi * (1 - xi) * L * L / 2;
      break;
    case FrameStress::Vz: // Shear, Vz, interpolation
      bp(i, 2) = (0.5 - xi) * L;
      break;
    case FrameStress::T: // Torsion, T, interpolation
      break;
    default: break;
    }
  }
}

void
ForceFrame3d<NIP,nsr,nwm>::getForceInterpolatMatrix(double xi, Matrix& b, const ID& code)
{
  b.Zero();

  double L = basic_system->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case FrameStress::N: // Axial, P, interpolation
      b(i, 0) = 1.0;
      break;
    case FrameStress::Vy: // Shear, Vy, interpolation
      b(i, 1) = b(i, 2) = 1.0 / L;
      break;
    case FrameStress::Vz: // Shear, Vz, interpolation
      b(i, 3) = b(i, 4) = 1.0 / L;
      break;
    case FrameStress::T: // Torque, T, interpolation
      b(i, 5) = 1.0;
      break;
    case FrameStress::My: // Moment, My, interpolation
      b(i, 3) = xi - 1.0;
      b(i, 4) = xi;
      break;
    case FrameStress::Mz: // Moment, Mz, interpolation
      b(i, 1) = xi - 1.0;
      b(i, 2) = xi;
      break;
    default: break;
    }
  }
}
#endif

