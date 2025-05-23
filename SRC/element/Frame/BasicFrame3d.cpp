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

void
BasicFrame3d::zeroLoad()
{
  q0.zero();
  p0.zero();

  wx = 0.0;
  wy = 0.0;
  wz = 0.0;
  eleLoads.clear();
  return;
}

int 
BasicFrame3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  constexpr int releasey = 0;
  constexpr int releasez = 0;
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
  // double L = theCoordTransf->getInitialLength();

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

    // Fixed end forces in basic system
    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;
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

      double P  = wa * (b - a);
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
    opserr << "BasicFrame3d::addLoad -- load type unknown\n";
    return -1;
  }

  return 0;
}




void
//BasicFrame3d::computeReactions(VectorND<6>& p0)
BasicFrame3d::computeReactions(double* p0)
{

  for (auto[load, loadFactor] : eleLoads) {

    int type;
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
BasicFrame3d::addReactionGrad(double* dp0dh, int gradNumber, double dLdh)
{


  for (auto[load, loadFactor] : eleLoads) {


    int type;
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