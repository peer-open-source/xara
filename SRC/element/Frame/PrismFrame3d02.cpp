//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// 3D prismatic frame element with:
//
// - Type 3 (uniform) warping from shear and torsion and 
// - symmetric cross section,
//
// Written: cmp 2024
//
#include <Frame/BasicFrame3d.h>
#include "PrismFrame3d02.h"
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>
#include <Exponential.h>
#include <BasicFrameTransf.h>
#include <runtime/commands/modeling/transform/FrameTransformBuilder.hpp>


PrismFrame3d02::PrismFrame3d02(int tag, std::array<int, 2>& nodes,
                           double  A, double  E, double  G, 
                           double jx, double iy, double iz,
                           FrameTransformBuilder& tb,
                           double r, int cm,
                           int rz, int ry,
                           int geom,
                           int shear_flag)
  : FiniteElement<2, 3, 6> (tag, ELE_TAG_ElasticBeam3d, nodes), 
    BasicFrame3d(),
    basic_system(new BasicFrameTransf3d<6>(tb.template create<2,6>())),
    A(A), E(E), G(G), Jx(jx), 
    Iy(iy), Iz(iz), Iyz(0),
    mass_flag(cm), density(r),
    releasez(rz), releasey(ry),
    geom_flag(geom),
    section_tag(-1),
    q_send(q)
{
  q.zero();
  ke.zero();
  kg.zero();
  km.zero();
}


PrismFrame3d02::PrismFrame3d02(int tag, 
                           std::array<int,2>& nodes,
                           FrameSection &section,  
                           FrameTransformBuilder& tb,
                           double density, int cMass, bool use_mass, 
                           int rz, int ry,
                           int geom,
                           int shear_flag
                           )
  : FiniteElement<2, 3, 6>(tag, ELE_TAG_ElasticBeam3d, nodes),
    BasicFrame3d(),
    basic_system(new BasicFrameTransf3d<6>(tb.template create<2,6>())),
    mass_flag(cMass), density(density),
    releasez(rz), releasey(ry),
    geom_flag(geom),
    q_send(q)
{
  q.zero();

  // 1) Get Area properties
  section_tag = section.getTag();
  section.getIntegral(Field::Unit,   State::Init, A);
  section.getIntegral(Field::UnitZZ, State::Init, Iy);
  section.getIntegral(Field::UnitYY, State::Init, Iz);

  {
    // 2) Get Young and Shear Modulus
    const ID& layout = section.getType();
    const Matrix& Ks = section.getInitialTangent();
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::N) {
        E = Ks(i,i)/A;
      }
      else if (layout(i) == FrameStress::Vy) {
        G = Ks(i,i)/A;
      } else if (layout(i) == FrameStress::T) {
        G = Ks(i,i)/(Iy + Iz);
      }
    }
  }

  // 3) Remaining Constants
  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
  };

  MatrixND<6,6> Kc = section.getTangent<6,scheme>(State::Init);

  Jx = Kc(3,3)/G;
  Iyz = 0.0;

  if (!shear_flag) {
    Ay = Az = 0.0;
  } else {
    Ay = Kc(1,1)/G;
    Az = Kc(2,2)/G;
  }

  // TODO
  if (!use_mass) {
    if (section.getIntegral(Field::Density, State::Init, density) == 0) {
      ;
    }
  }
}

int
PrismFrame3d02::setNodes()
{

  
  if (basic_system->initialize(theNodes[0], theNodes[1]) != 0) {
      return -1;
  }
  int status = 0;

  L = basic_system->getInitialLength();
  this->BasicFrame3d::setLength(L);

  if (L == 0.0) {
    opserr << "PrismFrame3d02::setDomain  tag: " << this->getTag() << " -- Element has zero length\n";
    total_mass = 0;
    twist_mass = 0;
    phiY = 0;
    phiZ = 0;
    return -1;
  }

  if (Ay != 0)
    phiY = 12.0 * E * Iz / (L * L * G * Ay);
  else
    phiY = 0.0;
  if (Az != 0)
    phiZ = 12.0 * E * Iy / (L * L * G * Az);
  else
    phiZ = 0.0;
  //

  total_mass = density*L; 
  twist_mass = (density/A)*Jx*L;

  return 0;
}

int
PrismFrame3d02::commitState()
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "PrismFrame3d02::commitState () - failed in base class";
  }    
  retVal += basic_system->commitState();
  return retVal;
}


int
PrismFrame3d02::revertToLastCommit()
{
  return basic_system->revertToLastCommit();
}


int
PrismFrame3d02::revertToStart()
{
  q.zero();
  return basic_system->revertToStart();
}


int
PrismFrame3d02::update()
{
  return basic_system->update();
}



const Matrix &
PrismFrame3d02::getTangentStiff()
{
  thread_local MatrixND<12,12> kb;
  thread_local Matrix wrapper(kb);

  kb = this->getBasicTangent(State::Pres, 0);

  basic_system->t.push(kb, q, Transform::Total);

  return wrapper;
}

const Matrix &
PrismFrame3d02::getInitialStiff()
{
  return basic_system->getInitialGlobalStiffMatrix(this->getBasicTangent(State::Init, 0));
}

const Vector &
PrismFrame3d02::getResistingForce()
{

  VectorND<12> du = basic_system->t.getStateVariation();
  q += this->getBasicTangent(State::Pres, 0)*du;

  basic_system->t.push(q);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (total_mass != 0.0)
    q_send.addVector(1.0, p_iner, -1.0);

  return q_send;

}

OpenSees::MatrixND<12,12>
PrismFrame3d02::getBasicTangent(State flag, int rate)
{
  const double EA = E*A;
  const double EIy = E*Iy;
  const double EIz = E*Iz;
  return MatrixND<12,12>{{
    EA/L,        0,      0,     0, 0, 0, -EA/L, 0, 0, 0, 0, 0,
       0,        0, 12*EIy/L/L, 0, 0, 6*EIy/L, 0, 0, -12*EIy/L/L, 0, 0, 6*EIy/L,
       0,        0,      0,     12*EIz/L/L, 0, -6*EIz/L, 0, 0, 0, -12*EIz/L/L, 0, -6*EIz/L,
       0, -6*EIz/L,      0,     0, 0, 0, 0, 6*EIz/L, 0, 0,
  }};
}


const Matrix &
PrismFrame3d02::getMass()
{
    if (total_mass == 0.0) {

        thread_local MatrixND<12,12> M;
        M.zero();
        thread_local Matrix Wrapper{M};
        return Wrapper;

    } else if (mass_flag == 0)  {
        // Lumped mass matrix

        thread_local MatrixND<12,12> M;
        M.zero();
        thread_local Matrix Wrapper{M};
        double m = 0.5*total_mass;
        M(0,0) = m;
        M(1,1) = m;
        M(2,2) = m;
        M(6,6) = m;
        M(7,7) = m;
        M(8,8) = m;
        return Wrapper;

    } else {
      // Consistent (cubic, prismatic) mass matrix

      if (!shear_flag) {
        double L  = basic_system->getInitialLength();
        double m  = total_mass/420.0;
        double mx = twist_mass;
        thread_local MatrixND<12,12> M{0};

        M(0,0) = M(6,6) = m*140.0;
        M(0,6) = M(6,0) = m*70.0;

        M(3,3) = M(9,9) = mx/3.0; // Twisting
        M(3,9) = M(9,3) = mx/6.0;

        M( 2, 2) = M( 8, 8) =  m*156.0;
        M( 2, 8) = M( 8, 2) =  m*54.0;
        M( 4, 4) = M(10,10) =  m*4.0*L*L;
        M( 4,10) = M(10, 4) = -m*3.0*L*L;
        M( 2, 4) = M( 4, 2) = -m*22.0*L;
        M( 8,10) = M(10, 8) = -M(2,4);
        M( 2,10) = M(10, 2) =  m*13.0*L;
        M( 4, 8) = M( 8, 4) = -M(2,10);

        M( 1, 1) = M( 7, 7) =  m*156.0;
        M( 1, 7) = M( 7, 1) =  m*54.0;
        M( 5, 5) = M(11,11) =  m*4.0*L*L;
        M( 5,11) = M(11, 5) = -m*3.0*L*L;
        M( 1, 5) = M( 5, 1) =  m*22.0*L;
        M( 7,11) = M(11, 7) = -M(1,5);
        M( 1,11) = M(11, 1) = -m*13.0*L;
        M( 5, 7) = M( 7, 5) = -M(1,11);

        // Transform local mass matrix to global system
        return basic_system->getGlobalMatrixFromLocal(M);
    }
    else {
      Matrix mlTrn(12, 12), mlRot(12, 12), ml(12, 12);
      mlTrn.Zero();
      mlRot.Zero();
      ml.Zero();
      double c1x  = density * L / 210.0;
      mlTrn(0, 0) = mlTrn(6, 6) = c1x * 70.0;
      mlTrn(0, 6) = mlTrn(6, 0) = c1x * 35.0;
      double c2x                = density / A * Jx * L / 210.0;
      mlTrn( 3, 3) = mlTrn( 9, 9) = c2x * 70.0;
      mlTrn( 3, 9) = mlTrn( 9, 3) = c2x * 35.0;
      double c1y                = c1x / pow(1.0 + phiY, 2);
      mlTrn( 2, 2) = mlTrn( 8, 8) = c1y * (70.0 * phiY * phiY + 147.0 * phiY + 78.0);
      mlTrn( 2, 8) = mlTrn( 8, 2) = c1y * (35.0 * phiY * phiY + 63.0 * phiY + 27.0);
      mlTrn( 4, 4) = mlTrn(10, 10) = c1y * L * L / 4.0 * (7.0 * phiY * phiY + 14.0 * phiY + 8.0);
      mlTrn( 4,10) = mlTrn(10, 4) = -c1y * L * L / 4.0 * (7.0 * phiY * phiY + 14.0 * phiY + 6.0);
      mlTrn( 2, 4) = mlTrn(4, 2) = -c1y * L / 4.0 * (35.0 * phiY * phiY + 77.0 * phiY + 44.0);
      mlTrn( 8,10) = mlTrn(10, 8) = -mlTrn(2, 4);
      mlTrn( 2,10) = mlTrn(10, 2) = c1y * L / 4.0 * (35.0 * phiY * phiY + 63.0 * phiY + 26.0);
      mlTrn( 4, 8) = mlTrn(8, 4) = -mlTrn(2, 10);
      double c2y                = density / A * Iy / (30.0 * L * pow(1.0 + phiY, 2));
      mlRot(2, 2) = mlRot(8, 8) = c2y * 36.0;
      mlRot(2, 8) = mlRot(8, 2) = -mlRot(2, 2);
      mlRot(4, 4) = mlRot(10, 10) = c2y * L * L * (10.0 * phiY * phiY + 5.0 * phiY + 4.0);
      mlRot(4, 10) = mlRot(10, 4) = c2y * L * L * (5.0 * phiY * phiY - 5.0 * phiY - 1.0);
      mlRot(2, 4) = mlRot(4, 2) = mlRot(2, 10) = mlRot(10, 2) = c2y * L * (15.0 * phiY - 3.0);
      mlRot(4, 8) = mlRot(8, 4) = mlRot(8, 10) = mlRot(10, 8) = -mlRot(2, 4);
      double c1z                                              = c1x / pow(1.0 + phiZ, 2);
      mlTrn(1, 1) = mlTrn(7, 7) = c1z * (70.0 * phiZ * phiZ + 147.0 * phiZ + 78.0);
      mlTrn(1, 7) = mlTrn(7, 1) = c1z * (35.0 * phiZ * phiZ + 63.0 * phiZ + 27.0);
      mlTrn(5, 5) = mlTrn(11, 11) = c1z * L * L / 4.0 * (7.0 * phiZ * phiZ + 14.0 * phiZ + 8.0);
      mlTrn(5, 11) = mlTrn(11, 5) = -c1z * L * L / 4.0 * (7.0 * phiZ * phiZ + 14.0 * phiZ + 6.0);
      mlTrn(1, 5) = mlTrn(5, 1) = c1z * L / 4.0 * (35.0 * phiZ * phiZ + 77.0 * phiZ + 44.0);
      mlTrn(7, 11) = mlTrn(11, 7) = -mlTrn(1, 5);
      mlTrn(1, 11) = mlTrn(11, 1) = -c1z * L / 4.0 * (35.0 * phiZ * phiZ + 63.0 * phiZ + 26.0);
      mlTrn(5, 7) = mlTrn(7, 5) = -mlTrn(1, 11);
      double c2z                = density / A * Iz / (30.0 * L * pow(1.0 + phiZ, 2));
      mlRot(1, 1) = mlRot(7, 7) = c2z * 36.0;
      mlRot(1, 7) = mlRot(7, 1) = -mlRot(1, 1);
      mlRot(5, 5) = mlRot(11, 11) = c2z * L * L * (10.0 * phiZ * phiZ + 5.0 * phiZ + 4.0);
      mlRot(5, 11) = mlRot(11, 5) = c2z * L * L * (5.0 * phiZ * phiZ - 5.0 * phiZ - 1.0);
      mlRot(1, 5) = mlRot(5, 1) = mlRot(1, 11) = mlRot(11, 1) = -c2z * L * (15.0 * phiZ - 3.0);
      mlRot(5, 7) = mlRot(7, 5) = mlRot(7, 11) = mlRot(11, 7) = -mlRot(1, 5);
      // add translational and rotational parts
      ml = mlTrn + mlRot;

      // Transform local mass matrix to global system
      return basic_system->getGlobalMatrixFromLocal(ml);
    }
  }
}


int
PrismFrame3d02::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
PrismFrame3d02::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
PrismFrame3d02::Print(OPS_Stream &s, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_ELEM_INDENT << "{";
      s << "\"name\": " << this->getTag();
      s << ", ";
      s << "\"type\": \"" << this->getClassType() << "\"";
      s << ", ";

      s << "\"nodes\": [" << node_tags(0) << ", " 
                          << node_tags(1) << "]";
      s << ", ";

      s << "\"massperlength\": " << total_mass/L;
      s << ", ";

      s << "\"releasez\": "<< releasez << ", ";
      s << "\"releasey\": "<< releasey << ", ";                
      s << "\"transform\": " << basic_system->getTag();
      s << ", ";

      // 
      if (section_tag > 0) {
        s << "\"section\": " << section_tag << ", ";
      }
      s << "\"E\": "  << E  << ", ";
      s << "\"G\": "  << G  << ", ";
      s << "\"A\": "  << A  << ", ";
      s << "\"Ay\": "  << Ay  << ", ";
      s << "\"Az\": "  << Az  << ", ";
      s << "\"Jx\": " << Jx << ", ";
      s << "\"Iy\": " << Iy << ", ";
      s << "\"Iz\": " << Iz;

      //
      s << "}";
  }
    
    this->getResistingForce(); 

    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_BEAM\t" << eleTag << "\t";
        s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
        s << "\t0\t0.0000000\n";
    }

    else if (flag == 2) {
        static Vector xAxis(3);
        static Vector yAxis(3);
        static Vector zAxis(3);

        basic_system->getLocalAxes(xAxis, yAxis, zAxis);

        s << "#PrismFrame3d02\n";
        s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
        s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
        s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << "\n";

        const Vector &xi = theNodes[0]->getCrds();
        const Vector &node2Crd = theNodes[1]->getCrds();
        const Vector &node1Disp = theNodes[0]->getDisp();
        const Vector &node2Disp = theNodes[1]->getDisp();

        s << "#NODE " << xi(0) << " " << xi(1) << " " << xi(2)
                << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
                << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << "\n";

        s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
                << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
                << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << "\n";
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\n  PrismFrame3d02: " << this->getTag() << "\n";
        s << "\tConnected Nodes: " << connectedExternalNodes;
        s << "\tCoordTransf: " << basic_system->getTag() << "\n";
        s << "\tmass density:  " << total_mass/L << ", mass_type: " << mass_flag << "\n";
        s << "\trelease about z:  " << releasez << "\n";
        s << "\trelease about y:  " << releasey << "\n";
        return;
    }
}


Response*
PrismFrame3d02::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = nullptr;
  const ID& node_tags = this->getExternalNodes();

  output.tag("ElementOutput");
  output.attr("eleType", this->getClassType());
  output.attr("eleTag",  this->getTag());
  output.attr("node1",  node_tags(0));
  output.attr("node2",  node_tags(1));
  
  // Global forces
  if (strcmp(argv[0],"force") == 0 || 
      strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || 
      strcmp(argv[0],"globalForces") == 0) {


    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Pz_1");
    output.tag("ResponseType","Mx_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","Mx_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, Vector(12));

  // local forces
  } else if (strcmp(argv[0],"localForce") == 0 || 
             strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","Vy_1");
    output.tag("ResponseType","Vz_1");
    output.tag("ResponseType","T_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","Vy_2");
    output.tag("ResponseType","Vz_2");
    output.tag("ResponseType","T_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 3, Vector(12));

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0 || 
             strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Mz_2");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","T");
    
    theResponse = new ElementResponse(this, 4, Vector(NBV));

  }  else if (strcmp(argv[0],"deformations") == 0 || 
              strcmp(argv[0],"basicDeformations") == 0) {
    
    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta11");
    output.tag("ResponseType","theta12");
    output.tag("ResponseType","theta21");
    output.tag("ResponseType","theta22");
    output.tag("ResponseType","phi");
    theResponse = new ElementResponse(this, 5, Vector(NBV));
  }

  else if (strcmp(argv[0],"sectionX") == 0) {
    if (argc > 2) {
      float xL = atof(argv[1]);
      if (xL < 0.0)
        xL = 0.0;
      if (xL > 1.0)
        xL = 1.0;
      if (strcmp(argv[2],"forces") == 0) {
        theResponse = new ElementResponse(this,6,Vector(6));
        Information &info = theResponse->getInformation();
        info.theDouble = xL;
      }
    }   
  }

  if (theResponse == nullptr)
    theResponse = basic_system->setResponse(argv, argc, output);

  output.endTag(); // ElementOutput

  return theResponse;
}

int
PrismFrame3d02::getResponse(int responseID, Information &info)
{
  double L = basic_system->getInitialLength();
  double oneOverL = 1.0/L;
  static Vector Res(12);
  Res = this->getResistingForce();
  static Vector s(6);

  switch (responseID) {
  case 1: // stiffness
    return info.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return info.setVector(Res);
    
  case 3: { // local forces
    static Vector P(12);
    P.Zero();
    double V, Mi, Mj, T;
    // Axial
    double N = q[0];
    P(6) =  N;
    P(0) = -N + p0[0];
    
    // Torsion
    T = q[5];
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shear along y
    Mi    = q[1];
    Mj    = q[2];
    P(5)  = Mi;
    P(11) = Mj;
    V     = (Mi + Mj)*oneOverL;
    P(1)  =  V+p0[1];
    P(7)  = -V+p0[2];
    
    // Moments about y and shear along z
    Mi    = q[3];
    Mj    = q[4];
    P(4)  = Mi;
    P(10) = Mj;
    V     = (Mi + Mj)*oneOverL;
    P(2)  = -V + p0[3];
    P(8)  =  V + p0[4];

    return info.setVector(P);
  } 
  case 4: // basic forces
    return info.setVector(q);

  case 5:
    return info.setVector(basic_system->getBasicTrialDisp());

  case 6: {
    double xL = info.theDouble;
    double x = xL*L;
    
    s(0) = q[0] + wx*(L-x);
    s(1) = q[1]*(xL-1.0) + q[2]*xL + 0.5*wy*x*(x-L);
    s(2) = (q[1] + q[2])/L + wy*(x-0.5*L);
    s(3) = q[3]*(xL-1.0) + q[4]*xL - 0.5*wz*x*(x-L);
    s(4) = (q[3] + q[4])/L - wz*(x-0.5*L);
    s(5) = q[5];

    return info.setVector(s);
  }
  default:
    break;
  }
  return -1;
}


const Vector&
PrismFrame3d02::getResistingForceSensitivity(int gradNumber)
{
  static Vector P(12);
  P.Zero();

  VectorND<NBV> dqdh = this->getBasicForceGrad(gradNumber);
  double dLdh = basic_system->getLengthGrad();
  // Transform forces
  double dp0dh[6];
  dp0dh[0] = 0.0;
  dp0dh[1] = 0.0;
  dp0dh[2] = 0.0;
  dp0dh[3] = 0.0;
  dp0dh[4] = 0.0;
  dp0dh[5] = 0.0;
  this->addReactionGrad(dp0dh, gradNumber, dLdh);
  Vector dp0dhVec(dp0dh, 6);

  if (basic_system->isShapeSensitivity()) {
    //
    // dqdh += K dvdh|_ug
    //
    dqdh.addMatrixVector(1.0, this->getBasicTangent(State::Pres, 0), 
                         basic_system->getBasicDisplFixedGrad(), 1.0);

    // dAdh^T q
    P = basic_system->getGlobalResistingForceShapeSensitivity(this->getBasicForce(), dp0dhVec, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += basic_system->getGlobalResistingForce(dqdh, dp0dhVec);

  return P;
}

int
PrismFrame3d02::getResponseSensitivity(int responseID, int gradNumber, Information& info)
{
  // Basic deformation sensitivity
  if (responseID == 3) {
    const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);
    return info.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(6);

    const Vector& dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);
#if 0
    dqdh.addMatrixVector(0.0, K_pres, dvdh, 1.0);
#endif
    dqdh.addVector(1.0, this->getBasicForceGrad(gradNumber), 1.0);

    return info.setVector(dqdh);
  }

  else
    return -1;
}


int
PrismFrame3d02::setParameter(const char **argv, int argc, Parameter &param)
{

  int status = this->FiniteElement<2,3,6>::setParameter(argv, argc, param);

  if (status != -1)
    return status;

  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(Param::E,  this);
  }

  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(Param::Rho, this);
  }

  if (strcmp(argv[0],"A") == 0) {
    param.setValue(A);
    return param.addObject(Param::A,  this);
  }

  if (strcmp(argv[0],"Ay") == 0) {
    param.setValue(Ay);
    return param.addObject(Param::Ay,  this);
  }

  if (strcmp(argv[0],"Az") == 0) {
    param.setValue(Az);
    return param.addObject(Param::Az,  this);
  }

  if (strcmp(argv[0],"Iy") == 0) {
    param.setValue(Iy);
    return param.addObject(Param::Iy, this);
  }

  if (strcmp(argv[0],"Iz") == 0) {
    param.setValue(Iz);
    return param.addObject(Param::Iz, this);
  }

  if (strcmp(argv[0],"G") == 0) {
    param.setValue(G);
    return param.addObject(Param::G, this);
  }

  if (strcmp(argv[0],"J") == 0) {
    param.setValue(Jx);
    return param.addObject(Param::J, this);
  }

  if (strcmp(argv[0],"releasez") == 0) {
    param.setValue(releasez);
    return param.addObject(Param::HingeZ, this);
  }

  if (strcmp(argv[0],"releasey") == 0) {
    param.setValue(releasey);
    return param.addObject(Param::HingeY, this);
  }
  return -1;
}

int
PrismFrame3d02::updateParameter(int param, Information &info)
{
    int status = this->FiniteElement<2,3,6>::updateParameter(param, info);
    if (status != -1)
      return status;

    switch (param) {
      case -1:
        return -1;
      case Param::E:
        E = info.theDouble;
        return 0;
      case Param::G:
        G = info.theDouble;
        return 0;
      case Param::A:
        A = info.theDouble;
        return 0;
      case Param::Iz:
        Iz = info.theDouble;
        return 0;
      case Param::Iy:
        Iy = info.theDouble;
        return 0;
      case Param::J:
        Jx = info.theDouble;
        return 0;

      case Param::Rho:
        rho = info.theDouble;
        return 0;
      
      case Param::HingeZ:
        releasez = (int)info.theDouble;
        if (releasez < 0 || releasez > 3)
          releasez = 0;
        return 0;

      case Param::HingeY:
        releasey = (int)info.theDouble;
        if (releasey < 0 || releasey > 3)
          releasey = 0;
        return 0;

      default:
        return -1;
    }

    // Update the element state
    formBasicStiffness(km);
}


VectorND<12>
PrismFrame3d02::getBasicForceGrad(int gradNumber)
{

  double L   = basic_system->getInitialLength();

  double
      dEIy = 0.0,
      dEIz = 0.0,
      dEA  = 0.0,
      dGAy = 0.0,
      dGAz = 0.0,
      dGJ  = 0.0;

  switch (parameterID) {
    case Param::E:
      dEIy = this->Iy;
      dEIz = this->Iz;
      dEA  = this->A;
      break;
    case Param::A:
      dEA  = this->E;
      break;
    case Param::Iy:
      dEIy = this->E;
      break;
    case Param::Iz:
      dEIz = this->E;
      break;
    default:
//    if (basic_system->isShapeSensitivity()) {
//      dEIy = E*Iy*dLi;
//      dEIz = E*Iz*dLi;
//      dEA  = E*A*dLi;
//    }
      break;
      ;
  }
 
  MatrixND<6,6> dK;
  dK.zero();
  dK(0,0) = dEA/L;
  dK(5,5) = dGJ/L;
  if (releasez == 0) {
    dK(1,1) = dK(2,2) = 4.0*dEIz/L;
    dK(2,1) = dK(1,2) = 2.0*dEIz/L;
  }
  if (releasez == 1)   // release I
    dK(2,2) = 3.0*dEIz/L;

  if (releasez == 2)   // release J
    dK(1,1) = 3.0*dEIz/L;


  if (releasey == 0) {
    dK(3,3) = dK(4,4) = 4.0*dEIy/L;
    dK(4,3) = dK(3,4) = 2.0*dEIy/L;
  }
  if (releasey == 1)   // release I
    dK(4,4) = 3.0*dEIy/L;

  if (releasey == 2)   // release J
    dK(3,3) = 3.0*dEIy/L;

  const Vector &v  = basic_system->getBasicTrialDisp();

  VectorND<NBV> dq = dK*v;

  double dL = basic_system->getLengthGrad();

  if (dL != 0.0)
    dq.addMatrixVector(1.0, this->getBasicTangent(State::Pres,0), v, -dL/L);

  return dq;
}

