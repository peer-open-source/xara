//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <Node.h>
#include <Parameter.h>
#include <FrameLoad.h>
#include <element/Frame/BasicFrame3d.h>

using OpenSees::VectorND;

BasicFrame3d::~BasicFrame3d()
{
  if (theCoordTransf)
    delete theCoordTransf;
}

int
BasicFrame3d::update()
{
  return 0;
}


// double
// BasicFrame3d::getLength(State state)
// {
//   if (state == State::Init)
//     return theCoordTransf->getInitialLength();

//   else
//     return theCoordTransf->getDeformedLength();
// }


int
BasicFrame3d::setNodes()
{
  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
      opserr << "PrismFrame3d::setDomain  tag: " 
             << this->getTag()
             << " -- Error initializing coordinate transformation\n";
      return -1;
  }

  return 0;
}

#if 0
const Vector &
BasicFrame3d::getResistingForce()
{
  VectorND<6> q  = this->getBasicForce();

  double q0 = q[0];
  double q1 = q[1];
  double q2 = q[2];
  double q3 = q[3];
  double q4 = q[4];
  double q5 = q[5];

  double oneOverL = 1.0 / theCoordTransf->getInitialLength();

  thread_local VectorND<12> pl;
  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  = q5;                     // Tj
  pl[10] = q4;
  pl[11] = q2;

  thread_local VectorND<12> pf{0.0};
  pf[0] = p0[0];
  pf[1] = p0[1];
  pf[7] = p0[2];
  pf[2] = p0[3];
  pf[8] = p0[4];


  thread_local VectorND<12> pg;
  thread_local Vector wrapper(pg);
//    const Vector p0Vec(p0);
//    P = theCoordTransf->getGlobalResistingForce(q, p0Vec);
  pg  = theCoordTransf->pushResponse(pl);
  pg += theCoordTransf->pushConstant(pf);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (total_mass != 0.0)
    wrapper.addVector(1.0, p_iner, -1.0);

  return wrapper;
}
#endif

const Matrix &
BasicFrame3d::getTangentStiff()
{

  VectorND<6>   q  = this->getBasicForce();
  MatrixND<6,6> kb = this->getBasicTangent(State::Pres, 0);

  return theCoordTransf->getGlobalStiffMatrix(Matrix(kb), Vector(q));
//    q += q0; // TODO!!! move this into PrismFrame and maybe DisplFrame

  double q0 = q[0];
  double q1 = q[1];
  double q2 = q[2];
  double q3 = q[3];
  double q4 = q[4];
  double q5 = q[5];

  double oneOverL = 1.0 / theCoordTransf->getInitialLength();

  thread_local VectorND<12> pl;
  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  =  q5;                    // Tj
  pl[10] =  q4;
  pl[11] =  q2;
  

  // Transform basic stiffness to local system
  static MatrixND<12,12> kl;  // Local stiffness
  static double tmp[6][12];  // Temporary storage
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][ 0] = -kb(i, 0);
    tmp[i][ 1] =  oneOverL * (kb(i, 1) + kb(i, 2));
    tmp[i][ 2] = -oneOverL * (kb(i, 3) + kb(i, 4));
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

  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) = -tmp[0][i];
    kl( 1, i) =  oneOverL * (tmp[1][i] + tmp[2][i]);
    kl( 2, i) = -oneOverL * (tmp[3][i] + tmp[4][i]);
    kl( 3, i) = -tmp[5][i];
    kl( 4, i) = tmp[3][i];
    kl( 5, i) = tmp[1][i];
    kl( 6, i) = tmp[0][i];
    kl( 7, i) = -kl(1, i);
    kl( 8, i) = -kl(2, i);
    kl( 9, i) = tmp[5][i];
    kl(10, i) = tmp[4][i];
    kl(11, i) = tmp[2][i];
  }


  static MatrixND<12,12> Kg;
  static Matrix Wrapper(Kg);
  Kg = theCoordTransf->pushResponse(kl, pl);
  return Wrapper;
}


int
BasicFrame3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (total_mass == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "BasicFrame3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  // add ( - fact * M R * accel ) to unbalance
  if (cMass == 0) {
    // take advantage of lumped mass matrix
    double m = 0.5*total_mass;

    p_iner[0] -= m * Raccel1(0);
    p_iner[1] -= m * Raccel1(1);
    p_iner[2] -= m * Raccel1(2);
    
    p_iner[6] -= m * Raccel2(0);
    p_iner[7] -= m * Raccel2(1);
    p_iner[8] -= m * Raccel2(2);

  } else  {
    // TODO: Move this to FiniteElement::getAcceleration() ?

    // use matrix vector multip. for consistent mass matrix
    static VectorND<12> Raccel;
    for (int i=0; i<6; i++)  {
      Raccel[i]   = Raccel1[i];
      Raccel[i+6] = Raccel2[i];
    }
    p_iner.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

const Vector &
BasicFrame3d::getResistingForceIncInertia()
{
  static VectorND<12> P_{0.0};
  static Vector P(P_);
  P = this->getResistingForce(); 
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);


  if (total_mass == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double m = 0.5*total_mass;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);

    P(6) += m * accel2(0);
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);

  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(12);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
  
  return P;
}

const Matrix &
BasicFrame3d::getInitialStiff()
{
  return theCoordTransf->getInitialGlobalStiffMatrix(this->getBasicTangent(State::Init, 0));
}

const Matrix &
BasicFrame3d::getMass()
{
    if (total_mass == 0.0) {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        return Wrapper;

    } else if (cMass == 0)  {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
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
        // consistent (cubic) mass matrix

        // get initial element length
        double L  = theCoordTransf->getInitialLength();
        double m  = total_mass/420.0;
        double mx = twist_mass;
        static MatrixND<12,12> ml{0};
        ml(0,0) = ml(6,6) = m*140.0;
        ml(0,6) = ml(6,0) = m*70.0;

        ml(3,3) = ml(9,9) = mx/3.0; // Twisting
        ml(3,9) = ml(9,3) = mx/6.0;

        ml( 2, 2) = ml( 8, 8) =  m*156.0;
        ml( 2, 8) = ml( 8, 2) =  m*54.0;
        ml( 4, 4) = ml(10,10) =  m*4.0*L*L;
        ml( 4,10) = ml(10, 4) = -m*3.0*L*L;
        ml( 2, 4) = ml( 4, 2) = -m*22.0*L;
        ml( 8,10) = ml(10, 8) = -ml(2,4);
        ml( 2,10) = ml(10, 2) =  m*13.0*L;
        ml( 4, 8) = ml( 8, 4) = -ml(2,10);

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

void 
BasicFrame3d::zeroLoad()
{
  q0.zero();
  p0.zero();

  wx = 0.0;
  wy = 0.0;
  wz = 0.0;
  // return BasicFrame3d::zeroLoad();
  p_iner.Zero();
  eleLoads.clear();
  return;
}

int 
BasicFrame3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  //
  // TODO:
  //
  // maintain map: {theLoad->getTag() : (theLoad, p0, q0, factor)}
  // if load not already in map, then integrate to get the shape;
  // otherwise just set the load factor
  //


  //
  // a. Store the load for computeReactions()
  //
  eleLoads.push_back({theLoad, loadFactor});

  //
  // b. Add to p0.
  //
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  // if (type == LOAD_TAG_FrameLoad && loadFactor == 1.0)
  //   frame_loads.push_back((FrameLoad*)theLoad);

  // else if (type == LOAD_TAG_FrameLoad && loadFactor == 1.0)
  //     frame_loads.push_back((FrameLoad*)theLoad);

  // else 
  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)

    this->wx += wx; // NOTE: This is not done in DispBeamColumn; why??
    this->wy += wy;
    this->wz += wz;
    
    double P  =     wx*L;
    double Vy = 0.5*wy*L;
    double Vz = 0.5*wz*L;

    // Reactions in basic system (projections on linear shape functions)
    p0[0] -=  P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    double Mz = Vy/6.0*L; // wy*L*L/12
    double My = Vz/6.0*L; // wz*L*L/12
    q0[0] -= 0.5*P;
    if (releasez == 0) {
      q0[1] -= Mz;
      q0[2] += Mz;
    }
    if (releasez == 1)
      q0[2] += wy/8*L*L;
      
    if (releasez == 2)
      q0[1] -= wy/8*L*L;
    
    if (releasey == 0) {
      q0[3] += My;
      q0[4] -= My;
    }
    if (releasey == 1)
      q0[4] -= wz/8*L*L;

    if (releasey == 2)
      q0[3] += wz/8*L*L;
  }

  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py     = data(0)*loadFactor;
    double Pz     = data(1)*loadFactor;
    double N      = data(2)*loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = Py*(1.0-a/L);
    double V2 = Py*a/L;
    p0[1] -= V1;
    p0[2] -= V2;
    p0[3] = Pz*(1.0-a/L); // V1
    p0[4] = Pz*a/L; // V2


    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*a/L;
    double M1 = -a * b2 * Py * L2;
    double M2 = a2 *  b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }

  else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      const double wa = data(2) * loadFactor;  // Axial
      const double wy = data(0) * loadFactor;  // Transverse
      const double wz = data(1) * loadFactor;  // Transverse
      const double a  = data(3) * L;
      const double b  = data(4) * L;
      const double c  = 0.5 * (b + a);
      double cOverL = c / L;

      double P = wa * (b - a);
      double Fy = wy * (b - a);
      double Fz = wz * (b - a);

      // Reactions in basic system
      p0[0] -= P;
      double V1, V2;
      V1 = Fy * (1.0 - c/L);
      V2 = Fy * c/L;
      p0[1] -= V1;
      p0[2] -= V2;
      V1 = Fz * (1.0 - c/L);
      V2 = Fz * cOverL;
      p0[3] -= V1;
      p0[4] -= V2;

      // Fixed end forces in basic system
      q0[0] -= P * cOverL;
      double beta2 = (1 - cOverL) * (1 - cOverL);
      double alfa2 = (cOverL) * (cOverL);
      double gamma2 = (b - a) / L;
      gamma2 *= gamma2;

      double M1 = -wy * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
      double M2 =  wy * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
      q0[1] += M1;
      q0[2] += M2;
      M1 = -wz * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
      M2 = wz * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
      q0[3] -= M1;
      q0[4] -= M2;
  }

  else {
    opserr << "BasicFrame3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << "\n";
    return -1;
  }

  return 0;
}

const Matrix &
BasicFrame3d::getMassSensitivity(int gradNumber)
{
  // From DispBeamColumn
  static MatrixND<12,12> M_;
  static Matrix M(M_);
  M.Zero();
  
  if (total_mass == 0.0 || parameterID != 1)
    return M;
  
  double L = theCoordTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    //double m = 0.5*rho*L;
    double m = 0.5*L;
    M(0,0) = M(1,1) = M(2,2) = M(6,6) = M(7,7) = M(8,8) = m;
  }

  else {
    // consistent mass matrix
    static Matrix ml(12,12);
    //double m = rho*L/420.0;
    double m = L/420.0;
    ml(0,0) = ml(6,6) = m*140.0;
    ml(0,6) = ml(6,0) = m*70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // TODO: CURRENTLY NO TORSIONAL MASS 
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // TODO: CURRENTLY NO TORSIONAL MASS
    
    ml(2, 2) = ml( 8, 8) =  m*156.0;
    ml(2, 8) = ml( 8, 2) =  m*54.0;
    ml(4, 4) = ml(10,10) =  m*4.0*L*L;
    ml(4,10) = ml(10, 4) = -m*3.0*L*L;
    ml(2, 4) = ml( 4, 2) = -m*22.0*L;
    ml(8,10) = ml(10, 8) = -ml(2,4);
    ml(2,10) = ml(10, 2) =  m*13.0*L;
    ml(4, 8) = ml( 8, 4) = -ml(2,10);
    
    ml(1, 1) = ml(7,7) = m*156.0;
    ml(1, 7) = ml(7,1) = m*54.0;
    ml(5, 5) = ml(11,11) = m*4.0*L*L;
    ml(5,11) = ml(11,5) = -m*3.0*L*L;
    ml(1, 5) = ml(5,1) = m*22.0*L;
    ml(7,11) = ml(11,7) = -ml(1,5);
    ml(1,11) = ml(11,1) = -m*13.0*L;
    ml(5, 7) = ml(7,5) = -ml(1,11);
    
    // transform local mass matrix to global system
    M = theCoordTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return M;
}

int
BasicFrame3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // don't do anything if MaterialStageParameter calls this element
  if (strcmp(argv[0],"updateMaterialStage") == 0) {
    return -1;
  }
  
  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }

  // moment release
  if (strcmp(argv[0],"releasez") == 0) {
    param.setValue(releasez);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"releasey") == 0) {
    param.setValue(releasey);
    return param.addObject(8, this);
  }  

  return -1;
}

int
BasicFrame3d::updateParameter(int paramID, Information &info)
{
    switch (paramID) {
      case -1:
        return -1;

      case 1:
        rho = info.theDouble;
        return 0;
      
      case 7:
        releasez = (int)info.theDouble;
        if (releasez < 0 || releasez > 3)
          releasez = 0;
        return 0;
      case 8:
        releasey = (int)info.theDouble;
        if (releasey < 0 || releasey > 3)
          releasey = 0;
        return 0;                          
      default:
        return -1;  
    }
}

int
BasicFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID; 
  return 0;
}


void
//BasicFrame3d::computeReactions(VectorND<6>& p0)
BasicFrame3d::computeReactions(double* p0)
{
  int type;
  double L = theCoordTransf->getInitialLength();

  for (auto[load, loadFactor] : eleLoads) {

    const  Vector& data = load->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data[0] * loadFactor; // Transverse
      double wz = data[1] * loadFactor; // Transverse
      double wa = data[2] * loadFactor; // Axial

      p0[0] -= wa * L;
      double V = 0.5 * wy * L;
      p0[1] -= V;
      p0[2] -= V;
      V = 0.5 * wz * L;
      p0[3] -= V;
      p0[4] -= V;

    }
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy = data(0) * loadFactor;  // Transverse Y at start
      double wz = data(1) * loadFactor;  // Transverse Z at start
      double wa = data(2) * loadFactor;  // Axial at start
      double a = data(3) * L;
      double b = data(4) * L;
      double wyb = data(5) * loadFactor;  // Transverse Y at end
      double wzb = data(6) * loadFactor;  // Transverse Z at end
      double wab = data(7) * loadFactor;  // Axial at end
      p0[0] -= wa * (b - a) + 0.5 * (wab - wa) * (b - a);
      double c = a + 0.5 * (b - a);
      double Fy = wy * (b - a); // resultant transverse load Y (uniform part)
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      double Fz = wz * (b - a); // resultant transverse load Z (uniform part)
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
      c = a + 2.0 / 3.0 * (b - a);
      Fy = 0.5 * (wyb - wy) * (b - a); // resultant transverse load Y (triang. part)
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      Fz = 0.5 * (wzb - wz) * (b - a); // resultant transverse load Z (triang. part)
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
    }

    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * loadFactor;
      double Pz     = data(1) * loadFactor;
      double N      = data(2) * loadFactor;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double V1 = Py * (1.0 - aOverL);
      double V2 = Py * aOverL;
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
      V1 = Pz * (1.0 - aOverL);
      V2 = Pz * aOverL;
      p0[3] -= V1;
      p0[4] -= V2;
    }
  }
}

void
BasicFrame3d::addReactionGrad(double* dp0dh, int gradNumber)
{
  int type;
  double L = theCoordTransf->getInitialLength();

  double dLdh = theCoordTransf->getLengthGrad();

//for (int i = 0; i < numEleLoads; i++) {
  for (auto[load, loadFactor] : eleLoads) {

//  ElementalLoad& load = *eleLoads[i];

    const Vector& data = load->getData(type, 1.0);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0) * 1.0; // Transverse
      double wz = data(1) * 1.0; // Transverse
      double wa = data(2) * 1.0; // Axial

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dwydh       = sens(0);
      double dwzdh       = sens(1);
      double dwadh       = sens(2);

      //p0[0] -= wa*L;
      dp0dh[0] -= wa * dLdh + dwadh * L;

      //double V = 0.5*wy*L;
      //p0[1] -= V;
      //p0[2] -= V;
      double dVdh = 0.5 * (wy * dLdh + dwydh * L);
      dp0dh[1] -= dVdh;
      dp0dh[2] -= dVdh;
      dVdh = 0.5 * (wz * L + dwzdh * L);
      dp0dh[3] -= dVdh;
      dp0dh[4] -= dVdh;
    } else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * 1.0;
      double Pz     = data(1) * 1.0;
      double N      = data(2) * 1.0;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dPydh       = sens(0);
      double dPzdh       = sens(1);
      double dNdh        = sens(2);
      double daLdh       = sens(3);

      //double a = aOverL*L;

      //double V1 = Py*(1.0-aOverL);
      //double V2 = Py*aOverL;
      double dV1dh = Py * (0.0 - daLdh) + dPydh * (1.0 - aOverL);
      double dV2dh = Py * daLdh + dPydh * aOverL;

      //p0[0] -= N;
      //p0[1] -= V1;
      //p0[2] -= V2;
      dp0dh[0] -= dNdh;
      dp0dh[1] -= dV1dh;
      dp0dh[2] -= dV2dh;

      dV1dh = Pz * (0.0 - daLdh) + dPzdh * (1.0 - aOverL);
      dV2dh = Pz * daLdh + dPzdh * aOverL;
      dp0dh[3] -= dV1dh;
      dp0dh[4] -= dV2dh;
    }
  }
}
