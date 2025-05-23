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
// Claudio Perez
//
#include <math.h>
#include <string.h>
#include <float.h>
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
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>

#define ELE_TAG_ForceFrame3d 0 // TODO

#define THREAD_LOCAL  static
#define ALWAYS_STATIC static // used when we need to do things like return a Matrix reference


// Constructor invoked by FEM_ObjectBroker; 
// recvSelf() needs to be invoked on this object.
ForceFrame3d::ForceFrame3d()
 : BasicFrame3d(0, ELE_TAG_ForceFrame3d),
   stencil(nullptr),
   max_iter(0), tol(0.0),
   state_flag(0),
   density(0.0), mass_flag(0), use_density(false),
   mass_initialized(false),
   Ki(nullptr),
   parameterID(0)
{
  K_pres.zero();
  K_save.zero();
  q_save.zero();
  q_pres.zero();
}

ForceFrame3d::ForceFrame3d(int tag, std::array<int, 2>& nodes, 
                           std::vector<FrameSection*>& sec,
                           BeamIntegration& bi,
                           FrameTransform3d& coordTransf, 
                           double dens, int mass_flag_, bool use_density,
                           int max_iter_, double tolerance
                           )
 : BasicFrame3d(tag, ELE_TAG_ForceFrame3d, nodes, coordTransf),
   stencil(nullptr),
   state_flag(0),
   Ki(nullptr),
   density(dens), mass_flag(mass_flag_), use_density(use_density),
   mass_initialized(false),
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


ForceFrame3d::~ForceFrame3d()
{

  for (GaussPoint& point : points)
    if (point.material != nullptr)
      delete point.material;

  if (stencil != nullptr)
    delete stencil;

  if (Ki != nullptr)
    delete Ki;
}


int
ForceFrame3d::setNodes()
{
  this->BasicFrame3d::setNodes();

  double L = theCoordTransf->getInitialLength();

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

int
ForceFrame3d::getIntegral(Field field, State state, double& total)
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

int
ForceFrame3d::commitState()
{
  int status = 0;

  // Call element commitState to do any base class stuff
  if ((status = this->Element::commitState()) != 0)
    opserr << "ForceFrame3d::commitState () - failed in base class";

  for (GaussPoint& point : points) {
    point.es_save = point.es;
    if (point.material->commitState() != 0)
      return -1;
  }


  // commit the transformation between coord. systems
  if ((status = theCoordTransf->commitState()) != 0)
    return status;

  // commit the element variables state
  K_save = K_pres;
  q_save = q_pres;

  // state_flag = 0;  fmk - commented out, see what happens to Example3.1.tcl if uncommented
  //                       - i have not a clue why, ask remo if he ever gets in contact with us again!

  return status;
}

int
ForceFrame3d::revertToLastCommit()
{
  for (GaussPoint& point : points) {
    FrameSection& section = *point.material;

    point.es = point.es_save;
    if (section.revertToLastCommit() != 0)
      return -1;

    section.setTrialState<nsr,scheme>(point.es);
    point.sr = section.getResultant<nsr,scheme>();
    point.Fs = section.getFlexibility<nsr,scheme>();

  }

  // Revert the transformation to last commit
  if (theCoordTransf->revertToLastCommit() != 0)
    return -2;

  // Revert the element state to last commit
  q_pres = q_save;
  K_pres = K_save;

  state_flag = 0;

  return 0;
}


int
ForceFrame3d::revertToStart()
{
  // Revert the transformation to start
  if (theCoordTransf->revertToStart() != 0)
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

VectorND<6>&
ForceFrame3d::getBasicForce()
{
  return q_pres;
}

MatrixND<6, 6>&
ForceFrame3d::getBasicTangent(State state, int rate)
{
  return K_pres;
}

const Matrix &
ForceFrame3d::getMass()
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

      ALWAYS_STATIC MatrixND<12,12> M{0.0};
      ALWAYS_STATIC Matrix Wrapper{M};
      return Wrapper;

  } else if (mass_flag == 0)  {

      ALWAYS_STATIC MatrixND<12,12> M{0.0};
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

  } else {
      // consistent (cubic, prismatic) mass matrix

      double L  = theCoordTransf->getInitialLength();
      double m  = total_mass/420.0;
      double mx = twist_mass;
      ALWAYS_STATIC MatrixND<12,12> ml{0.0};

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
      return theCoordTransf->getGlobalMatrixFromLocal(ml);
  }
}


const Matrix &
ForceFrame3d::getInitialStiff()
{
  // check for quick return
  if (Ki != nullptr)
    return *Ki;

  THREAD_LOCAL MatrixND<NBV,NBV> F_init;
  this->getInitialFlexibility(F_init);

  // calculate element stiffness matrix
  THREAD_LOCAL MatrixND<NBV, NBV> K_init;
  if (F_init.invert(K_init) < 0)
    opserr << "ForceFrame3d::getInitialStiff -- could not invert flexibility";

  Matrix wrapper(K_init);
  Ki = new Matrix(theCoordTransf->getInitialGlobalStiffMatrix(wrapper));

  return *Ki;
}


void
ForceFrame3d::initializeSectionHistoryVariables()
{
  for (int i = 0; i < points.size(); i++) {
    points[i].Fs.zero();
    points[i].es.zero();
    points[i].sr.zero();
    points[i].es_save.zero();
  }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int
ForceFrame3d::update()
{
  this->BasicFrame3d::update();

  // TODO: remove hard limit on sections
  THREAD_LOCAL VectorND<nsr>     es_trial[maxNumSections]; //  strain
  THREAD_LOCAL VectorND<nsr>     sr_trial[maxNumSections]; //  stress resultant
  THREAD_LOCAL MatrixND<nsr,nsr> Fs_trial[maxNumSections]; //  flexibility


  // If we have completed a recvSelf() do a revertToLastCommit()
  // to get sr, etc. set correctly
  if (state_flag == 2)
    this->revertToLastCommit();

  // update the transformation
  theCoordTransf->update();

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  //
  const Vector& v = theCoordTransf->getBasicTrialDisp();

  THREAD_LOCAL VectorND<NBV> dv{0.0};
  dv = theCoordTransf->getBasicIncrDeltaDisp();

  if (state_flag != 0 && (dv.norm() <= DBL_EPSILON) && (eleLoads.size()==0))
    return 0;

  THREAD_LOCAL VectorND<NBV> Dv;
  Dv  = v;
  Dv -= dv;


  THREAD_LOCAL VectorND<NBV> vr;       // element residual deformations
  THREAD_LOCAL MatrixND<NBV, NBV> F;   // element flexibility matrix
  THREAD_LOCAL VectorND<NBV> dvToDo,
                             dv_trial;

  dvToDo  = dv;
  dv_trial = dvToDo;

  const double factor = 10.0;
  double dW;             // section strain energy (work) norm
  double dW0  = 0.0;


  //
  //   Iterate to find compatible forces and deformations
  //
  //   First try first a Newton iteration, if that fails we try an initial
  //   flexibility iteration on first iteration and then regular Newton, if
  //   that fails we use the initial flexiblity for all iterations.
  //   If they both fail we subdivide dV & try to get compatible forces
  //   and deformations. if they work and we have subdivided we apply
  //   the remaining dV.
  //
  enum class Strategy {
    Newton, InitialIterations, InitialThenNewton
  };
  static constexpr std::array<Strategy,3> solve_strategy {
    Strategy::Newton, Strategy::InitialIterations, Strategy::InitialThenNewton
  };

  int subdivision = 1;
  bool converged   = false;

  while (converged == false && subdivision <= maxSubdivisions) {

    for (Strategy strategy : solve_strategy ) {

      // Allow 10 times more iterations for initial tangent strategy
      const int numIters = (strategy==Strategy::InitialIterations) ? 10*max_iter : max_iter;

      for (int i = 0; i < points.size(); i++) {
        es_trial[i]  = points[i].es;
        Fs_trial[i]  = points[i].Fs;
        sr_trial[i]  = points[i].sr;
      }

      if (state_flag == 2)
        continue;

      VectorND<NBV>     q_trial = q_pres;
      MatrixND<NBV, NBV> K_trial = K_pres;

      q_trial += K_pres*dv;

      for (int j = 0; j < numIters; j++) {
        // Initialize F and vr for integration
        F.zero();
        vr.zero();

        //
        // Gauss Loop
        //
        for (int i = 0; i < points.size(); i++) {
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
          VectorND<nsr> si {
               q_trial[0],                           // N
               (q_trial[1] + q_trial[2])/L,          // VY
               (q_trial[3] + q_trial[4])/L,          // VZ
               q_trial[5],                           // T
               xL1 * q_trial[3] + xL * q_trial[4],   // MY
               xL1 * q_trial[1] + xL * q_trial[2],   // MZ
          };
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

            VectorND<nsr> ds; //, des;

            // Form stress increment ds from last iteration
            // ds = si - si_past;
            ds = si;
            ds.addVector(1.0, sr, -1.0);

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
                  MatrixND<nsr,nsr> Fs0 = section.getFlexibility<nsr,scheme>(State::Init);
                  es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
                } else
                  es_trial[i].addMatrixVector(1.0, Fs, ds, 1.0);
                break;

              case Strategy::InitialIterations:
                //  Newton with initial tangent
                MatrixND<nsr,nsr> Fs0 = section.getFlexibility<nsr,scheme>(State::Init);
                es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
                break;
            }
          }


          //
          // c. Set trial section state and get response
          //
          if (section.setTrialState<nsr,scheme>(es_trial[i]) < 0) {
            opserr << "ForceFrame3d::update - section " << i << " failed in setTrial\n";
            return -1;
          }

          sr_trial[i] = section.getResultant<nsr, scheme>();
          Fs_trial[i] = section.getFlexibility<nsr, scheme>();

          //
          // d. Integrate element flexibility matrix
          //
          //    F += (B' * Fs * B) * wi * L;
          //
          {
            MatrixND<nsr,NBV> FsB;
            FsB.zero();
            for (int jj = 0; jj < nsr; jj++) {
                // SECTION_RESPONSE_P:
                FsB(jj, 0) += Fs_trial[i](jj, 0) * wtL;

                // SECTION_RESPONSE_VY:
                FsB(jj, 1) += Fs(jj, 1) * wtL * jsx;
                FsB(jj, 2) += Fs(jj, 1) * wtL * jsx;

                // SECTION_RESPONSE_VZ:
                FsB(jj, 3) += Fs(jj, 2) * wtL * jsx;
                FsB(jj, 4) += Fs(jj, 2) * wtL * jsx;

                // SECTION_RESPONSE_T:
                FsB(jj, 5) += Fs(jj, 3) * wtL;

                // SECTION_RESPONSE_MY:
                FsB(jj, 3) += xL1 * Fs(jj, 4)*wtL;
                FsB(jj, 4) += xL  * Fs(jj, 4)*wtL;

                // SECTION_RESPONSE_MZ:
                FsB(jj, 1) += xL1 * Fs(jj, 5)*wtL;
                FsB(jj, 2) += xL  * Fs(jj, 5)*wtL;
            }

            for (int jj = 0; jj < NBV; jj++) {
                double tmp;
                // SECTION_RESPONSE_P:
                F(0, jj) += FsB( 0, jj);

                // SECTION_RESPONSE_VY:
                tmp = jsx * FsB( 1, jj);
                F(1, jj) += tmp;
                F(2, jj) += tmp;

                // SECTION_RESPONSE_VZ:
                tmp = jsx * FsB( 2, jj);
                F(3, jj) += tmp;
                F(4, jj) += tmp;

                // SECTION_RESPONSE_T:
                F(5, jj) += FsB( 3, jj);

                // SECTION_RESPONSE_MY:
                F(3, jj) += xL1 * FsB( 4, jj);
                F(4, jj) += xL  * FsB( 4, jj);

                // SECTION_RESPONSE_MZ:
                F(1, jj) += xL1 * FsB( 5, jj);
                F(2, jj) += xL  * FsB( 5, jj);
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

            // SECTION_RESPONSE_P:
            vr[0] += des[0]*wtL;
            // SECTION_RESPONSE_VY:
            vr[1] += jsx*des[1]*wtL;
            vr[2] += jsx*des[1]*wtL;
            // SECTION_RESPONSE_VZ:
            vr[3] += jsx*des[2]*wtL;
            vr[4] += jsx*des[2]*wtL;
            // SECTION_RESPONSE_T:
            vr[5] += des[3]*wtL;
            // SECTION_RESPONSE_MY:
            vr[3] += xL1 * des[4]*wtL;
            vr[4] += xL  * des[4]*wtL;
            // SECTION_RESPONSE_MZ:
            vr[1] += xL1 * des[5]*wtL;
            vr[2] += xL  * des[5]*wtL;
          }
        } // Gauss loop


        //
        // Finalize trial element state 
        //
        //    K_trial  = inv(F)
        //    q_trial += K * (Dv + dv_trial - vr)
        //
        if (F.invert(K_trial) < 0)
          opserr << "ForceFrame3d::update -- could not invert flexibility\n";


        // dv = Dv + dv_trial  - vr
        dv = Dv;
        dv += dv_trial;
        dv -= vr;

        VectorND<NBV> dqe = K_trial * dv;
      //dqe.addMatrixVector(0.0, K_trial, dv, 1.0);

        dW = dqe.dot(dv);
        if (dW0 == 0.0)
          dW0 = dW;

        q_trial += dqe;


        //
        // Check for convergence of this interval
        //
        if (fabs(dW) < tol) {

          // Set the target displacement
          dvToDo -= dv_trial;
          Dv     += dv_trial;

          // Check if we have got to where we wanted
          if (dvToDo.norm() <= DBL_EPSILON) {
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

          for (int k = 0; k < points.size(); k++) {
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
    opserr << "for element " << this->getTag() 
           << "; dW = " << dW << ", dW0 = " << dW0
           << "\n";
    return -1;
  }

  state_flag = 1;

  return 0;
}


const Matrix &
ForceFrame3d::getTangentStiff()
{
  MatrixND<6,6> kb = this->getBasicTangent(State::Pres, 0);

  double q0 = q_pres[0];
  double q1 = q_pres[1];
  double q2 = q_pres[2];
  double q3 = q_pres[3];
  double q4 = q_pres[4];
  double q5 = q_pres[5];

  double L = theCoordTransf->getInitialLength();

  THREAD_LOCAL VectorND<12> pl{0.0};
  pl[0]  = -q0;               // Ni
  pl[1]  =  (q1 + q2)/L;      // Viy
  pl[2]  = -(q3 + q4)/L;      // Viz
  pl[3]  = -q5;               // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;               // Nj
  pl[7]  = -pl[1];            // Vjy
  pl[8]  = -pl[2];            // Vjz
  pl[9]  =  q5;               // Tj
  pl[10] =  q4;
  pl[11] =  q2;
  

  // Transform basic stiffness to local system
  THREAD_LOCAL double tmp[12][12]{};  // Temporary storage
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][ 0] = -kb(i, 0);
    tmp[i][ 1] =  (kb(i, 1) + kb(i, 2))/L;
    tmp[i][ 2] = -(kb(i, 3) + kb(i, 4))/L;
    tmp[i][ 3] = -kb(i, 5);
    tmp[i][ 4] =  kb(i, 3);
    tmp[i][ 5] =  kb(i, 1);
    tmp[i][ 6] =  kb(i, 0);
    tmp[i][ 7] = -tmp[i][1];
    tmp[i][ 8] = -tmp[i][2];
    tmp[i][ 9] =  kb(i, 5);
    tmp[i][10] =  kb(i, 4);
    tmp[i][11] =  kb(i, 2);
  }

  THREAD_LOCAL MatrixND<12,12> kl{0.0};  // Local stiffness
  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) =  -tmp[0][i];
    kl( 1, i) =  (tmp[1][i] + tmp[2][i])/L;
    kl( 2, i) = -(tmp[3][i] + tmp[4][i])/L;
    kl( 3, i) =  -tmp[5][i];
    kl( 4, i) =   tmp[3][i];
    kl( 5, i) = tmp[1][i];
    kl( 6, i) = tmp[0][i];
    kl( 7, i) = -kl(1, i);
    kl( 8, i) = -kl(2, i);
    kl( 9, i) = tmp[5][i];
    kl(10, i) = tmp[4][i];
    kl(11, i) = tmp[2][i];
  }


  THREAD_LOCAL MatrixND<12,12> Kg;
  THREAD_LOCAL Matrix Wrapper(Kg);
  Kg = theCoordTransf->pushResponse(kl, pl);
  return Wrapper;
}



void
ForceFrame3d::addLoadAtSection(VectorND<nsr>& sp, double x)
{
  double L = theCoordTransf->getInitialLength();

  for (auto[load, loadFactor] : eleLoads) {

    int type;
    const Vector& data = load->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wa = data(2) * loadFactor; // Axial
      double wy = data(0) * loadFactor; // Transverse
      double wz = data(1) * loadFactor; // Transverse

      for (int ii = 0; ii < nsr; ii++) {
        switch (scheme[ii]) {
        case SECTION_RESPONSE_P:  sp[ii] += wa * (L - x); break;
        case SECTION_RESPONSE_VY: sp[ii] += wy * (x - 0.5 * L); break;
        case SECTION_RESPONSE_VZ: sp[ii] += wz * (0.5 * L - x); break;
        case SECTION_RESPONSE_MZ: sp[ii] += wy * 0.5 * x * (x - L); break;
        case SECTION_RESPONSE_MY: sp[ii] += wz * 0.5 * x * (L - x); break;
        default:                  break;
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
          case SECTION_RESPONSE_P:  sp[ii] +=       N; break;
          case SECTION_RESPONSE_VY: sp[ii] -=     Vy1; break;
          case SECTION_RESPONSE_VZ: sp[ii] -=     Vz1; break;
          case SECTION_RESPONSE_MY: sp[ii] += x * Vz1; break;
          case SECTION_RESPONSE_MZ: sp[ii] -= x * Vy1; break;
          default:                  break;
          }
        } else {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_MZ: sp[ii] -= (L - x) * Vy2; break;
          case SECTION_RESPONSE_VY: sp[ii] += Vy2; break;
          case SECTION_RESPONSE_MY: sp[ii] += (L - x) * Vz2; break;
          case SECTION_RESPONSE_VZ: sp[ii] += Vz2; break;
          default:                  break;
          }
        }
      }
    } 
  
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy = data(0) * loadFactor;  // Transverse Y at start
      double wz = data(1) * loadFactor;  // Transverse Z at start
      double wa = data(2) * loadFactor;  // Axial at start
      double a = data(3)*L;
      double b = data(4)*L;
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
          case SECTION_RESPONSE_P:
            sp(ii) += Fa;
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) -= VyI*x;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += VzI*x;
            break;            
          case SECTION_RESPONSE_VY:
            sp(ii) -= VyI;
            break;
          case SECTION_RESPONSE_VZ:
        sp(ii) += VzI;
            break;            
          default:
            break;
          }
        }
        else if (x >= b) {
          switch(scheme[ii]) {
          case SECTION_RESPONSE_MZ:
            sp(ii) += VyJ*(x-L);
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) -= VzJ*(x-L);
            break;            
          case SECTION_RESPONSE_VY:
            sp(ii) += VyJ;
            break;
          case SECTION_RESPONSE_VZ:
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
          case SECTION_RESPONSE_P:
            sp(ii) += Fa - wa * (x - a) - 0.5 * (wab - wa) / (b - a) * (x - a) * (x - a);
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) += -VyI * x + 0.5 * wy * (x - a) * (x - a) + 0.5 * (wyy - wy) * (x - a) * (x - a) / 3.0;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += VzI * x - 0.5 * wz * (x - a) * (x - a) - 0.5 * (wzz - wz) * (x - a) * (x - a) / 3.0;
            break;            
          case SECTION_RESPONSE_VY:
            sp(ii) += -VyI + wy * (x - a) + 0.5 * (wyy - wy) * (x - a);
            break;
          case SECTION_RESPONSE_VZ:           
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

void
ForceFrame3d::getStressGrad(VectorND<nsr>& dspdh, int isec, int gradNumber)
{

  double L    = theCoordTransf->getInitialLength();
  double dLdh = theCoordTransf->getLengthGrad();

  int numSections = points.size();

  double dxidh[maxNumSections];
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
        case SECTION_RESPONSE_P:
          //sp(ii) += wa*(L-x);
          dspdh(ii) += dwadh * (L - x) + wa * (dLdh - dxdh);
          break;
        case SECTION_RESPONSE_MZ:
          //sp(ii) += wy*0.5*x*(x-L);
          //dspdh(ii) += 0.5*(dwydh*x*(x-L) + wy*dxdh*(x-L) + wy*x*(dxdh-dLdh));
          dspdh(ii) += 0.5 * (dwydh * x * (x - L) + wy * (dxdh * (2 * x - L) - x * dLdh));
          break;
        case SECTION_RESPONSE_VY:
          //sp(ii) += wy*(x-0.5*L);
          dspdh(ii) += dwydh * (x - 0.5 * L) + wy * (dxdh - 0.5 * dLdh);
          break;
        case SECTION_RESPONSE_MY:
          //sp(ii) += wz*0.5*x*(L-x);
          //dspdh(ii) += 0.5*(dwzdh*x*(L-x) + wz*dxdh*(L-x) + wz*x*(dLdh-dxdh));
          dspdh(ii) += 0.5 * (dwzdh * x * (L - x) + wz * (dxdh * (L - 2 * x) + x * dLdh));
          break;
        case SECTION_RESPONSE_VZ:
          //sp(ii) += wz*(x-0.5*L);
          dspdh(ii) += dwzdh * (0.5 * L - x) + wz * (0.5 * dLdh - dxdh);
          break;
        default: break;
        }
      }
    } else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * 1.0;
      double Pz     = data(1) * 1.0;
      double N      = data(2) * 1.0;
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
          case SECTION_RESPONSE_P:
            //sp(ii) += N;
            dspdh(ii) += dNdh;
            break;
          case SECTION_RESPONSE_MZ:
            //sp(ii) -= x*Vy1;
            dspdh(ii) -= (dxdh * Vy1 + x * dVy1dh);
            break;
          case SECTION_RESPONSE_VY:
            //sp(ii) -= Vy1;
            dspdh(ii) -= dVy1dh;
            break;
          case SECTION_RESPONSE_MY:
            //sp(ii) += x*Vz1;
            dspdh(ii) += (dxdh * Vz1 + x * dVz1dh);
            break;
          case SECTION_RESPONSE_VZ:
            //sp(ii) -= Vz1;
            dspdh(ii) -= dVz1dh;
            break;
          default: break;
          }
        } else {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_VY:
            dspdh(ii) += dVy2dh;
            break;
          case SECTION_RESPONSE_VZ:
            dspdh(ii) += dVz2dh;
            break;
          case SECTION_RESPONSE_MY:
            //sp(ii) += (L-x)*Vz2;
            dspdh(ii) += (dLdh - dxdh) * Vz2 + (L - x) * dVz2dh;
            break;
          case SECTION_RESPONSE_MZ:
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
             << this->getTag() << "\n";
    }
  }
}

int
ForceFrame3d::sendSelf(int commitTag, Channel& theChannel)
{
  return -1;
}

int
ForceFrame3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  return -1;
}

// addBFsB(F, s, i)
int
ForceFrame3d::getInitialFlexibility(MatrixND<NBV,NBV>& fe)
{
  fe.zero();

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  for (int i = 0; i < numSections; i++) {

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    const MatrixND<nsr,nsr> fSec = points[i].material->getFlexibility<nsr, scheme>(State::Init);

    THREAD_LOCAL MatrixND<nsr, NBV> FB;
    FB.zero();
    double tmp;
    int ii, jj;
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 0) += fSec(jj, ii) * wtL;
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = fSec(jj, ii) * wtL;
          FB(jj, 1) += xL1 * tmp;
          FB(jj, 2) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = jsx * fSec(jj, ii) * wtL;
          FB(jj, 1) += tmp;
          FB(jj, 2) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = fSec(jj, ii) * wtL;
          FB(jj, 3) += xL1 * tmp;
          FB(jj, 4) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = jsx * fSec(jj, ii) * wtL;
          FB(jj, 3) += tmp;
          FB(jj, 4) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 5) += fSec(jj, ii) * wtL;
        break;
      default: break;
      }
    }
    //
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < NBV; jj++)
          fe(0, jj) += FB(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = FB(ii, jj);
          fe(1, jj) += xL1 * tmp;
          fe(2, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = jsx * FB(ii, jj);
          fe(1, jj) += tmp;
          fe(2, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = FB(ii, jj);
          fe(3, jj) += xL1 * tmp;
          fe(4, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = jsx * FB(ii, jj);
          fe(3, jj) += tmp;
          fe(4, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int jj = 0; jj < NBV; jj++)
          fe(5, jj) += FB(ii, jj);
        break;
      default: break;
      }
    }
  }
  return 0;
}

int
ForceFrame3d::getInitialDeformations(Vector& v0)
{
  v0.Zero();
  if (eleLoads.size() < 1 || (this->setState(State::Init) != 0))
    return 0;

  double L = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;


  const int numSections = points.size();
  for (int i = 0; i < numSections; i++) {

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    THREAD_LOCAL VectorND<nsr> sp;
    sp.zero();

    this->addLoadAtSection(sp, points[i].point*L);

    MatrixND<nsr,nsr> fse = points[i].material->getFlexibility<nsr,scheme>(State::Init);

    VectorND<nsr> e = fse*sp;

    double tmp;
    for (int ii = 0; ii < nsr; ii++) {
      double dei = e[ii] * wtL;
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        v0(0) += dei;
        break;
      case SECTION_RESPONSE_VY:
        tmp = jsx * dei;
        v0[1] += tmp;
        v0[2] += tmp;
        break;
      case SECTION_RESPONSE_VZ:
        tmp = jsx * dei;
        v0[3] += tmp;
        v0[4] += tmp;
        break;
      case SECTION_RESPONSE_MY:
        v0[3] += xL1 * dei;
        v0[4] += xL  * dei;
        break;
      case SECTION_RESPONSE_MZ:
        v0[1] += xL1 * dei;
        v0[2] += xL  * dei;
        break;
      default: break;
      }
    }
  }

  return 0;
}


void
ForceFrame3d::Print(OPS_Stream& s, int flag)
{
  const int nip = points.size();
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\"";
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
    for (decltype(points.size()) i = 0; i < points.size() - 1; i++)
      s << points[i].material->getTag() << ", ";
    s << points[points.size() - 1].material->getTag() << "]";
    s << ", ";

    s << "\"transform\": " << theCoordTransf->getTag();
    s << ", ";

    s << "\"integration\": ";
    stencil->Print(s, flag);
    s << "}";
  }

  // flags with negative values are used by GSA
  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << points[0].material->getTag() << "\t" << points[nip - 1].material->getTag();
    s << "\t" << node_tags(0) << "\t" << node_tags(1);
    s << "\t0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: ForceFrame3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << nip;
    s << "\tMass density: " << density << "\n";
    stencil->Print(s, flag);
    double P     = q_save[0];
    double MZ1   = q_save[1];
    double MZ2   = q_save[2];
    double MY1   = q_save[3];
    double MY2   = q_save[4];
    double L     = theCoordTransf->getInitialLength();
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

OPS_Stream&
operator<<(OPS_Stream& s, ForceFrame3d& E)
{
  E.Print(s);
  return s;
}


Response*
ForceFrame3d::setResponse(const char** argv, int argc, OPS_Stream& output)
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

    theResponse = new ElementResponse(this, Respond::GlobalForce, Vector(12));

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

  } else if (strcmp(argv[0], "getRemCriteria1") == 0) {
    theResponse = new ElementResponse(this, 77, Vector(2));

  } else if (strcmp(argv[0], "getRemCriteria2") == 0) {
    theResponse = new ElementResponse(this, 8, Vector(2), ID(6));

  } else if (strcmp(argv[0], "RayleighForces") == 0 || 
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, Vector(12));

  } else if (strcmp(argv[0], "sections") == 0) {
    if (this->setState(State::Init) != 0)
      return nullptr;

    CompositeResponse* theCResponse = new CompositeResponse();
    int numResponse                 = 0;
    const int numSections = points.size();
    double L = theCoordTransf->getInitialLength();

    for (int i = 0; i < numSections; i++) {
      output.tag("GaussPointOutput");
      output.attr("number", i + 1);
      output.attr("eta", points[i].point * L);

      Response* theSectionResponse = points[i].material->setResponse(&argv[1], argc - 1, output);

      if (theSectionResponse != 0)
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
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(points.size(), 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(points.size(), 3));
  }

  else if (strcmp(argv[0], "cbdiDisplacements") == 0)
    theResponse = new ElementResponse(this, 112, Matrix(20, 3));


  else if (strstr(argv[0], "section") != 0) {

    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= points.size() && argc > 2) {
        if (this->setState(State::Init) != 0)
          return nullptr;
        double L = theCoordTransf->getInitialLength();
        const int numSections = points.size();

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
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        const int numSections = points.size();
        stencil->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", points[i].point * L);

          Response* theSectionResponse = points[i].material->setResponse(&argv[1], argc - 1, output);

          if (theSectionResponse != 0) {
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
    theResponse = theCoordTransf->setResponse(argv, argc, output);
  }

  output.endTag();

  return theResponse;
}

int
ForceFrame3d::getResponse(int responseID, Information& info)
{
  THREAD_LOCAL Vector vp(6);
  THREAD_LOCAL MatrixND<NBV,NBV> fe;

  if (responseID == Respond::GlobalForce)
    return info.setVector(this->getResistingForce());

  else if (responseID == Respond::LocalForce) {
    THREAD_LOCAL VectorND<12> v_resp{0.0};
    THREAD_LOCAL Vector v_wrap(v_resp);

    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);
    // Axial
    double N     = q_pres[0];
    v_resp(6) =  N;
    v_resp(0) = -N + p0[0];

    // Torsion
    double T     = q_pres[5];
    v_resp(9) =  T;
    v_resp(3) = -T;

    // Moments about z and shears along y
    double M1     = q_pres[1];
    double M2     = q_pres[2];
    v_resp(5)  = M1;
    v_resp(11) = M2;
    double L      = theCoordTransf->getInitialLength();
    double V      = (M1 + M2) / L;
    v_resp(1)  =  V + p0[1];
    v_resp(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1            = q_pres[3];
    M2            = q_pres[4];
    v_resp(4)  = M1;
    v_resp(10) = M2;
    V             = (M1 + M2) / L;
    v_resp(2)  = -V + p0[3];
    v_resp(8)  =  V + p0[4];

    return info.setVector(v_wrap);

  }

  // Chord rotation
  else if (responseID == 3) {
    vp = theCoordTransf->getBasicTrialDisp();
    return info.setVector(vp);
  }

  else if (responseID == Respond::BasicForce)
    return info.setVector(q_pres);

  else if (responseID == 19)
    return info.setMatrix(K_pres);

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = theCoordTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, q_pres, -1.0);
    Vector v0(6);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return info.setVector(vp);
  }

  else if (responseID == 10) {
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    double L = theCoordTransf->getInitialLength();

    Vector locs(points.size());
    for (int i = 0; i < points.size(); i++)
      locs[i] = points[i].point * L;

    return info.setVector(locs);
  }

  else if (responseID == 11) {
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    double L = theCoordTransf->getInitialLength();

    Vector weights(points.size());
    for (int i = 0; i < points.size(); i++)
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
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    const int numSections = points.size();
    stencil->getSectionLocations(numSections, L, pts);
    // CBDI influence matrix
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID& code  = points[i].material->getType();
      const Vector& e = points[i].material->getSectionDeformation();
      int order       = points[i].material->getOrder();
      for (int j = 0; j < nsr; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(numSections); // along local y
    Vector dispsz(numSections); // along local z
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    stencil->getSectionLocations(numSections, L, pts);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(numSections, 3);
    vp = theCoordTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i] * vp(0); // linear shape function
      uxb(1) = dispsy(i);
      uxb(2) = dispsz(i);
      if (responseID == 111)
        uxg = theCoordTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      else
        uxg = theCoordTransf->getPointLocalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
    }
    return info.setMatrix(disps);
  }

  else if (responseID == 112) {
    double L = theCoordTransf->getInitialLength();
    double ipts[maxNumSections];
    const int numSections = points.size();
    stencil->getSectionLocations(numSections, L, ipts);
    // CBDI influence matrix
    double pts[20];
    for (int i = 0; i < 20; i++)
      pts[i] = 1.0 / (20 - 1) * i;
    Matrix ls(20, numSections);
    getCBDIinfluenceMatrix(20, pts, numSections, ipts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID& code  = points[i].material->getType();
      const Vector& e = points[i].material->getSectionDeformation();
      int order       = points[i].material->getOrder();
      for (int j = 0; j < nsr; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(20); // along local y
    Vector dispsz(20); // along local z
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(20, 3);
    vp = theCoordTransf->getBasicTrialDisp();
    for (int i = 0; i < 20; i++) {
      uxb(0)      = pts[i] * vp[0]; // linear shape function
      uxb(1)      = dispsy(i);
      uxb(2)      = dispsz(i);
      uxg         = theCoordTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
    }
    return info.setMatrix(disps);
  }

  else if (responseID == 12)
    return info.setVector(this->getRayleighDampingForces());

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = theCoordTransf->getInitialLength();

    if (fabs(q_pres[1] + q_pres[2]) > DBL_EPSILON)
      LI(0) = q_pres[1] / (q_pres[1] + q_pres[2]) * L;

    if (fabs(q_pres[3] + q_pres[4]) > DBL_EPSILON)
      LI(1) = q_pres[3] / (q_pres[3] + q_pres[4]) * L;

    return info.setVector(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = theCoordTransf->getInitialLength();

    double wts[maxNumSections];
    const int numSections = points.size();
    stencil->getSectionWeights(numSections, L, wts);

    double pts[maxNumSections];
    stencil->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(q_pres[1] + q_pres[2]) > DBL_EPSILON)
      LIz = q_pres[1] / (q_pres[1] + q_pres[2]) * L;

    double LIy = 0.0;
    if (fabs(q_pres[3] + q_pres[4]) > DBL_EPSILON)
      LIy = q_pres[3] / (q_pres[3] + q_pres[4]) * L;

    for (int i = 0; i < numSections; i++) {
      double x       = points[i].point * L;
      const ID& type = points[i].material->getType();
      int order      = points[i].material->getOrder();
      double kappa   = 0.0;
      if (x < LIz) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += points[i].es[j];
        double b = -LIz + x;
        d2z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x < LIy) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += points[i].es[j];
        double b = -LIy + x;
        d2y += (wts[i] * L) * kappa * b;
      }
    }

    d2z += stencil->getTangentDriftI(L, LIz, q_pres[1], q_pres[2]);
    d2y += stencil->getTangentDriftI(L, LIy, q_pres[3], q_pres[4], true);

    for (int i = numSections - 1; i >= 0; i--) {
      double x       = pts[i] * L;
      const ID& type = points[i].material->getType();
      int order      = points[i].material->getOrder();
      double kappa   = 0.0;
      if (x > LIz) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += points[i].es[j];
        double b = x - LIz;
        d3z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x > LIy) {
        for (int j = 0; j < nsr; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += points[i].es[j];
        double b = x - LIy;
        d3y += (wts[i] * L) * kappa * b;
      }
    }

    d3z += stencil->getTangentDriftJ(L, LIz, q_pres[1], q_pres[2]);
    d3y += stencil->getTangentDriftJ(L, LIy, q_pres[3], q_pres[4], true);

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;

    return info.setVector(d);

  } else if (responseID == 77) { // Why is this here?
    return -1;

  } else if (responseID == 8) {

    ID* infoID = info.theID;

    int compID     = (*infoID)(0);
    int critID     = (*infoID)(1);
    int nTagbotn11 = (*infoID)(2);
    int nTagmidn11 = (*infoID)(3);
    int nTagtopn11 = (*infoID)(4);
    int globgrav11 = (*infoID)(5);

    const char* filenamewall = info.theString;

    // int returns
    double value       = 0.0;
    double checkvalue1 = 0.0;

    if (critID == 7) {
      Domain* theDomain = this->getDomain();

      double oofwallresp;
      // determine the in plane horizontal deformation axis
      // and the out of plane horizontal deformation axis
      Node* theNode1a      = theDomain->getNode(nTagbotn11);
      Node* theNode3a      = theDomain->getNode(nTagtopn11);
      const Vector& crdIa1 = theNode1a->getCrds();
      const Vector& crdJa1 = theNode3a->getCrds();
      int indwdir1;
      int indwdir2;
      if (globgrav11 == 1) {
        indwdir1 = 1;
        indwdir2 = 2;
      } else if (globgrav11 == 2) {
        indwdir1 = 0;
        indwdir2 = 2;
      } else if (globgrav11 == 3) {
        indwdir1 = 0;
        indwdir2 = 1;
      }

      double dir1a1    = crdJa1(indwdir1) - crdIa1(indwdir1);
      double dir2a1    = crdJa1(indwdir2) - crdIa1(indwdir2);
      double dirsumdum = sqrt(dir1a1 * dir1a1 + dir2a1 * dir2a1);
      double dir1inp   = dir1a1 / dirsumdum;
      double dir2inp   = dir2a1 / dirsumdum;

      double dir1oop = -dir2inp;
      double dir2oop = dir1inp;

      Node* theNode1                = theDomain->getNode(nTagbotn11);
      const Vector& theResponsewall = theNode1->getTrialDisp();
      double valbotinfn = theResponsewall(indwdir1) * dir1inp + theResponsewall(indwdir2) * dir2inp;
      double valbotoutfn =
          theResponsewall(indwdir1) * dir1oop + theResponsewall(indwdir2) * dir2oop;

      Node* theNode2                 = theDomain->getNode(nTagmidn11);
      const Vector& theResponsewall2 = theNode2->getTrialDisp();
      double valmidinfn =
          theResponsewall2(indwdir1) * dir1inp + theResponsewall2(indwdir2) * dir2inp;
      double valmidoutfn =
          theResponsewall2(indwdir1) * dir1oop + theResponsewall2(indwdir2) * dir2oop;

      Node* theNode3                 = theDomain->getNode(nTagtopn11);
      const Vector& theResponsewall3 = theNode3->getTrialDisp();
      double valtopinfn =
          theResponsewall3(indwdir1) * dir1inp + theResponsewall3(indwdir2) * dir2inp;
      double valtopoutfn =
          theResponsewall3(indwdir1) * dir1oop + theResponsewall3(indwdir2) * dir2oop;

      value             = sqrt(pow((valtopinfn - valbotinfn), 2.0));
      double valoutchck = valmidoutfn - (valtopoutfn + valbotoutfn) / 2.0;
      oofwallresp       = sqrt(pow(valoutchck, 2.0));
      //
      double outplanevaldat;  // variable for input value
      double inplanevaldat;   // variable for input value
      double outplanevaldat1; // variable for input value
      double inplanevaldat1;  // variable for input value
      std::ifstream indata;

      if (filenamewall != nullptr) {
        //
        indata.open(filenamewall); // opens the file
        if (!indata) {             // file couldn't be opened
          opserr << "ForceFrame3d::getResponse"
                 << " file for infill wall (" << filenamewall << " could not be opened" << "\n";
          return -1;
        }
        checkvalue1    = 0.0;
        int counterdum = 0;
        while (!indata.eof()) { // keep reading until end-of-file
          counterdum = counterdum + 1;
          indata >> outplanevaldat >> inplanevaldat; // sets EOF flag if no value found
          if (counterdum != 1) {
            if (oofwallresp >= outplanevaldat1 && oofwallresp <= outplanevaldat) {
              checkvalue1 = inplanevaldat1 + (oofwallresp - outplanevaldat1) /
                                                 (outplanevaldat - outplanevaldat1) *
                                                 (inplanevaldat - inplanevaldat1);
              break;
            }
          }
          indata >> outplanevaldat1 >> inplanevaldat1;
          if (oofwallresp >= outplanevaldat && oofwallresp <= outplanevaldat1) {
            checkvalue1 = inplanevaldat + (oofwallresp - outplanevaldat) /
                                              (outplanevaldat1 - outplanevaldat) *
                                              (inplanevaldat1 - inplanevaldat);
            break;
          }
        }
        indata.close();
      }

      static Vector result8(2);
      result8(0) = value;
      result8(1) = checkvalue1;

      return info.setVector(result8);
    }

    return -1;
  }
  //by SAJalali
  else if (responseID == 2000) {
    const int numSections = points.size();
    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    stencil->getSectionWeights(numSections, L, xi);
    double energy = 0;
    for (int i = 0; i < numSections; i++) {
      energy += points[i].material->getEnergy() * points[i].point * L;
    }
    return info.setDouble(energy);
  }

  return -1;
}

int
ForceFrame3d::getResponseSensitivity(int responseID, int gradNumber, Information& info)
{
  // Basic deformation sensitivity
  if (responseID == 3) {
    const Vector& dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);
    return info.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(6);
    dqdh.Zero();

    const Vector& dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);

    dqdh.addMatrixVector(0.0, K_pres, dvdh, 1.0);

    dqdh.addVector(1.0, this->getBasicForceGrad(gradNumber), 1.0);

    return info.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = info.theInt;
    int order      = points[sectionNum - 1].material->getOrder();

    VectorND<nsr> dsdh;
    dsdh.zero();

    if (eleLoads.size() > 0)
      this->getStressGrad(dsdh, sectionNum - 1, gradNumber);

    static Vector dqdh(6);
    {
      // Response 7
      const Vector& dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);

      dqdh.addMatrixVector(0.0, K_pres, dvdh, 1.0);

      dqdh.addVector(1.0, this->getBasicForceGrad(gradNumber), 1.0);
    }

    const int numSections = points.size();
    double L   = theCoordTransf->getInitialLength();
    double jsx = 1.0 / L;
    double pts[maxNumSections];
    stencil->getSectionLocations(numSections, L, pts);

    const ID& code = points[sectionNum - 1].material->getType();

    double xL  = pts[sectionNum - 1];
    double xL1 = xL - 1.0;

    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:  dsdh(ii) += dqdh(0); break;
      case SECTION_RESPONSE_MZ: dsdh(ii) += xL1 * dqdh(1) + xL * dqdh(2); break;
      case SECTION_RESPONSE_VY: dsdh(ii) += jsx * (dqdh(1) + dqdh(2)); break;
      case SECTION_RESPONSE_MY: dsdh(ii) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case SECTION_RESPONSE_VZ: dsdh(ii) += jsx * (dqdh(3) + dqdh(4)); break;
      case SECTION_RESPONSE_T:  dsdh(ii) += dqdh(5); break;
      default:                  dsdh(ii) += 0.0; break;
      }
    }

    double dLdh   = theCoordTransf->getLengthGrad();
    double d1oLdh = theCoordTransf->getd1overLdh();

    double dptsdh[maxNumSections];
    stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);
    double dxLdh = dptsdh[sectionNum - 1]; // - xL/L*dLdh;

    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        dsdh(j) += dxLdh * (q_pres[1] + q_pres[2]);
        //dsdh(j) -= dLdh*xL/L*(Se(1)+Se(2));
        break;
      case SECTION_RESPONSE_VY: dsdh(j) += d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: dsdh(j) += dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: dsdh(j) += d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    return info.setVector(dsdh);
  }

  // Plastic deformation sensitivity
  else if (responseID == 4) {
    static Vector dvpdh(6);

    const Vector& dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);

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

int
ForceFrame3d::setParameter(const char** argv, int argc, Parameter& param)
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
  if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      const int numSections = points.size();
      double xi[maxNumSections];
      double L = theCoordTransf->getInitialLength();
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
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= points.size())
      return points[sectionNum - 1].material->setParameter(&argv[2], argc - 2, param);

    else
      return -1;
  }

  // If the parameter belongs to all sections or lower
  if (strstr(argv[0], "allSections") != 0) {

    if (argc < 2)
      return -1;

    for (GaussPoint& point : points) {
      int ok = point.material->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }

    return result;
  }

  if (strstr(argv[0], "integration") != 0) {

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

int
ForceFrame3d::updateParameter(int parameterID, Information& info)
{
  if (parameterID == 1) {
    this->density = info.theDouble;
    return 0;
  } else
    return -1;
}

int
ForceFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

const Matrix&
ForceFrame3d::getKiSensitivity(int gradNumber)
{
  THREAD_LOCAL MatrixND<12,12> Ksen{0.0};
  THREAD_LOCAL Matrix wrapper(Ksen);
  return wrapper;
}

const Matrix&
ForceFrame3d::getMassSensitivity(int gradNumber)
{
  THREAD_LOCAL MatrixND<12,12> Msen{0.0};
  THREAD_LOCAL Matrix wrapper(Msen);

  double L = theCoordTransf->getInitialLength();
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

const Vector&
ForceFrame3d::getResistingForceSensitivity(int gradNumber)
{
  static Vector P(12);
  P.Zero();

  VectorND<NBV> dqdh = this->getBasicForceGrad(gradNumber);

  // Transform forces
  double dp0dh[6];
  dp0dh[0] = 0.0;
  dp0dh[1] = 0.0;
  dp0dh[2] = 0.0;
  dp0dh[3] = 0.0;
  dp0dh[4] = 0.0;
  dp0dh[5] = 0.0;
  this->addReactionGrad(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 6);

  if (theCoordTransf->isShapeSensitivity()) {
    //
    // dqdh += K dvdh|_ug
    //
    dqdh.addMatrixVector(1.0, K_pres, 
                         theCoordTransf->getBasicDisplFixedGrad(), 1.0);

    // dAdh^T q
    P = theCoordTransf->getGlobalResistingForceShapeSensitivity(q_pres, dp0dhVec, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += theCoordTransf->getGlobalResistingForce(dqdh, dp0dhVec);

  return P;
}

int
ForceFrame3d::commitSensitivity(int gradNumber, int numGrads)
{
  //
  //
  //

  int err = 0;

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();

  double pts[maxNumSections];
  double wts[maxNumSections];
  stencil->getSectionLocations(numSections, L, pts);
  stencil->getSectionWeights(numSections, L, wts);

  double dLdh = theCoordTransf->getLengthGrad();

  double dptsdh[maxNumSections];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double d1oLdh = theCoordTransf->getd1overLdh();

  // Response 7
  VectorND<6> dqdh = this->getBasicForceGrad(gradNumber);

  // dvdh = A u` + A` u
  const Vector& dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);
  dqdh.addMatrixVector(1.0, K_pres, dvdh, 1.0);


  //
  // Compute and commit de for each section
  //
  for (int i = 0; i < numSections; i++) {

    int order      = points[i].material->getOrder();
    const ID& code = points[i].material->getType();

    double xL  = pts[i];
    double xL1 = xL - 1.0;

    double dxLdh = dptsdh[i];

    VectorND<nsr> ds;
    ds.zero();

    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0)
      this->getStressGrad(ds, i, gradNumber);


    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:  ds(j) += dqdh(0); break;
      case SECTION_RESPONSE_VY: ds(j) += jsx * (dqdh(1) + dqdh(2)); break;
      case SECTION_RESPONSE_VZ: ds(j) += jsx * (dqdh(3) + dqdh(4)); break;
      case SECTION_RESPONSE_T:  ds(j) += dqdh(5); break;
      case SECTION_RESPONSE_MY: ds(j) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case SECTION_RESPONSE_MZ: ds(j) += xL1 * dqdh(1) + xL * dqdh(2); break;
      default:                  ds(j) += 0.0; break;
      }
    }

    const Vector& dsdh = points[i].material->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ: ds(j) += dxLdh  * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_VY: ds(j) += d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: ds(j) += dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: ds(j) += d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    Vector de(order);
    const Matrix& fs = points[i].material->getSectionFlexibility();
    de.addMatrixVector(0.0, fs, ds, 1.0);

    err += points[i].material->commitSensitivity(de, gradNumber, numGrads);
  }

  return err;
}

VectorND<6>
ForceFrame3d::getBasicForceGrad(int gradNumber)
{

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  double wts[maxNumSections];
  stencil->getSectionWeights(numSections, L, wts);

  double dLdh = theCoordTransf->getLengthGrad();

  double dptsdh[maxNumSections];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  stencil->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  double d1oLdh = theCoordTransf->getd1overLdh();

  //
  // Integrate dvdh
  //
  static Vector dvdh(6);
  dvdh.Zero();
  for (int i = 0; i < numSections; i++) {

    int order      = points[i].material->getOrder();
    const ID& code = points[i].material->getType();

    double xL  = points[i].point;
    double xL1 = xL - 1.0;
    double wtL = points[i].weight * L;

    double dxLdh  = dptsdh[i]; // - xL/L*dLdh;
    double dwtLdh = points[i].weight * dLdh + dwtsdh[i] * L;


    // Get section stress resultant gradient

    VectorND<nsr> dspdh;
    dspdh.zero();
    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0)
      this->getStressGrad(dspdh, i, gradNumber);

    Vector dsdh(order);
    dsdh = points[i].material->getStressResultantSensitivity(gradNumber, true);

    dsdh.addVector(1.0, dspdh, -1.0);


    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ: dsdh(j) -= dxLdh  * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_VY: dsdh(j) -= d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: dsdh(j) -= dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: dsdh(j) -= d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    Vector dedh(order);
    const Matrix& fs = points[i].material->getSectionFlexibility();
    dedh.addMatrixVector(0.0, fs, dsdh, 1.0);

    for (int j = 0; j < nsr; j++) {
      double dei = dedh(j) * wtL;
      switch (code(j)) {
      case SECTION_RESPONSE_P: dvdh(0) += dei; break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) += xL1 * dei;
        dvdh(2) += xL * dei;
        break;
      case SECTION_RESPONSE_VY:
        dei = jsx * dei;
        dvdh(1) += dei;
        dvdh(2) += dei;
        break;
      case SECTION_RESPONSE_MY:
        dvdh(3) += xL1 * dei;
        dvdh(4) += xL * dei;
        break;
      case SECTION_RESPONSE_VZ:
        dei = jsx * dei;
        dvdh(3) += dei;
        dvdh(4) += dei;
        break;
      case SECTION_RESPONSE_T:
        dvdh(5) += dei;
        break;
      default:
        break;
      }
    }

    const VectorND<nsr>& e = points[i].es;
    for (int j = 0; j < nsr; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P: dvdh(0) -= e(j) * dwtLdh; break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) -= xL1 * e(j) * dwtLdh;
        dvdh(2) -= xL * e(j) * dwtLdh;

        dvdh(1) -= dxLdh * e(j) * wtL;
        dvdh(2) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VY:
        dvdh(1) -= jsx * e(j) * dwtLdh;
        dvdh(2) -= jsx * e(j) * dwtLdh;

        dvdh(1) -= d1oLdh * e(j) * wtL;
        dvdh(2) -= d1oLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_MY:
        dvdh(3) -= xL1 * e(j) * dwtLdh;
        dvdh(4) -= xL * e(j) * dwtLdh;

        dvdh(3) -= dxLdh * e(j) * wtL;
        dvdh(4) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VZ:
        dvdh(3) -= jsx * e(j) * dwtLdh;
        dvdh(4) -= jsx * e(j) * dwtLdh;

        dvdh(3) -= d1oLdh * e(j) * wtL;
        dvdh(4) -= d1oLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_T:
        dvdh(5) -= e(j) * dwtLdh; 
        break;

      default:
        break;
      }
    }
  }

  return K_pres*dvdh;
}

const Matrix&
ForceFrame3d::computedfedh(int gradNumber)
{
  static Matrix dfedh(6, 6);

  dfedh.Zero();

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  const int numSections = points.size();
  double dLdh   = theCoordTransf->getLengthGrad();
  double d1oLdh = theCoordTransf->getd1overLdh();

  double dptsdh[maxNumSections];
  double dwtsdh[maxNumSections];
  stencil->getLocationsDeriv(numSections, L, dLdh, dptsdh);
  stencil->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  for (int i = 0; i < numSections; i++) {

    int order      = points[i].material->getOrder();
    const ID& code = points[i].material->getType();

    Matrix fb(order, NBV);
    Matrix fb2(order, NBV);

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
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < nsr; jj++) {
          fb(jj, 0) += dfsdh(jj, ii) * wtL; // 1

          //fb(jj,0) += fs(jj,ii)*dwtLdh; // 3

          //fb2(jj,0) += fs(jj,ii)*wtL; // 4
        }
        break;
      case SECTION_RESPONSE_MZ:
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
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = jsx * dfsdh(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
          // TODO: Need to complete for dLdh != 0
        }
        break;
      default: break;
      }
    }
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < NBV; jj++)
          dfedh(0, jj) += fb(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < NBV; jj++) {
          tmp = fb(ii, jj); // 1,2,3
          dfedh(1, jj) += xL1 * tmp;
          dfedh(2, jj) += xL * tmp;

          tmp = fb2(ii, jj); // 4
          //dfedh(1,jj) += dxLdh*tmp;
          //dfedh(2,jj) += dxLdh*tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
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

  return dfedh;
}

int
ForceFrame3d::setSectionPointers(std::vector<FrameSection*>& new_sections)
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


    // Check sections
    int sectionKey1 = -1;
    int sectionKey2 = -1;
    const ID& code = section->getType();
    for (int j = 0; j < code.Size(); j++) {
      if (code(j) == SECTION_RESPONSE_MZ)
        sectionKey1 = j;
      if (code(j) == SECTION_RESPONSE_MY)
        sectionKey2 = j;
    }
    if (sectionKey1 == -1) {
      opserr << "FATAL ForceFrame3d::compSectionResponse - section does not provide Mz response\n";
      return -1;
    }
    if (sectionKey2 == -1) {
      opserr << "FATAL ForceFrame3d::compSectionResponse - section does not provide My response\n";
      return -1;
    }
  }
  return 0;
}

const Vector &
ForceFrame3d::getResistingForce()
{
  double p0[5]{};
  
  if (eleLoads.size() > 0)
    this->computeReactions(p0);
 
  double q0 = q_pres[0];
  double q1 = q_pres[1];
  double q2 = q_pres[2];
  double q3 = q_pres[3];
  double q4 = q_pres[4];
  double q5 = q_pres[5];

  double L = theCoordTransf->getInitialLength();

  thread_local VectorND<12> pl;
  pl[ 0]  = -q0;                    // Ni
  pl[ 1]  =  (q1 + q2)/L;           // Viy
  pl[ 2]  = -(q3 + q4)/L;           // Viz
  pl[ 3]  = -q5;                    // Ti
  pl[ 4]  =  q3;
  pl[ 5]  =  q1;
  pl[ 6]  =  q0;                    // Nj
  pl[ 7]  = -pl[1];                 // Vjy
  pl[ 8]  = -pl[2];                 // Vjz
  pl[ 9]  = q5;                     // Tj
  pl[10]  = q4;
  pl[11]  = q2;

  thread_local VectorND<12> pf{0.0};
  pf[0] = p0[0];
  pf[1] = p0[1];
  pf[7] = p0[2];
  pf[2] = p0[3];
  pf[8] = p0[4];

  thread_local VectorND<12> pg;
  thread_local Vector wrapper(pg);

  pg  = theCoordTransf->pushResponse(pl);
  pg += theCoordTransf->pushConstant(pf);
  if (total_mass != 0.0)
    wrapper.addVector(1.0, p_iner, -1.0);

  return wrapper;
}



#if 0

void
ForceFrame3d::getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code)
{
  bp.Zero();

  double L = theCoordTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ: // Moment, Mz, interpolation
      bp(i, 1) = xi * (xi - 1) * L * L / 2;
      break;
    case SECTION_RESPONSE_P: // Axial, P, interpolation
      bp(i, 0) = (1 - xi) * L;
      break;
    case SECTION_RESPONSE_VY: // Shear, Vy, interpolation
      bp(i, 1) = (xi - 0.5) * L;
      break;
    case SECTION_RESPONSE_MY: // Moment, My, interpolation
      bp(i, 2) = xi * (1 - xi) * L * L / 2;
      break;
    case SECTION_RESPONSE_VZ: // Shear, Vz, interpolation
      bp(i, 2) = (0.5 - xi) * L;
      break;
    case SECTION_RESPONSE_T: // Torsion, T, interpolation
      break;
    default: break;
    }
  }
}

void
ForceFrame3d::getForceInterpolatMatrix(double xi, Matrix& b, const ID& code)
{
  b.Zero();

  double L = theCoordTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_P: // Axial, P, interpolation
      b(i, 0) = 1.0;
      break;
    case SECTION_RESPONSE_VY: // Shear, Vy, interpolation
      b(i, 1) = b(i, 2) = 1.0 / L;
      break;
    case SECTION_RESPONSE_VZ: // Shear, Vz, interpolation
      b(i, 3) = b(i, 4) = 1.0 / L;
      break;
    case SECTION_RESPONSE_T: // Torque, T, interpolation
      b(i, 5) = 1.0;
      break;
    case SECTION_RESPONSE_MY: // Moment, My, interpolation
      b(i, 3) = xi - 1.0;
      b(i, 4) = xi;
      break;
    case SECTION_RESPONSE_MZ: // Moment, Mz, interpolation
      b(i, 1) = xi - 1.0;
      b(i, 2) = xi;
      break;
    default: break;
    }
  }
}

#endif

