#pragma once
class FrameSection;

// GaussPoint:
//   shape[2][nen]
//
//   rotary_inertia[]
//   linear_inertia
//
// constant
// variable
// consistent
// true


enum FrameMassType {
};

int addRigidMass();
//    // Lumped mass matrix
//    double m = 0.5*total_mass;
//    for (int i=0; i<3; i++)  {
//        M(i,i)     = m;
//        M(i+6,i+6) = m;
//    }


template <int nen, typename MatT, typename GaussArray,  typename ShapeArray>
int
addTaperMass(MatT Mass, double factor, 
             GaussArray& gauss)
//           GaussArray& sections, 
//           ShapeArray& rotary, 
//           ShapeArray& linear)
{
  for (auto point : gauss) {
  }
  return -1;
}

template <typename ElemType, typename MatT, int ndf>
int
addPrismMass(MatT Mass, double factor, FrameSection& section, double length)
{

  double total_mass; // , twist_mass;

  MatrixND<12,12> M;
  if (total_mass == 0.0)
    return 0;

  // McCalley-Archer consistent mass matrix
  double L, A, Iy, Iz, Jx, phiY, phiZ;


  MatrixND<12,12> mlTrn, 
                  mlRot, 
                  ml;

  mlTrn.zero(); mlRot.zero(); ml.zero();

  double c1x = total_mass/210.0;
  mlTrn( 0, 0) = mlTrn( 6, 6) = c1x*70.0;
  mlTrn( 0, 6) = mlTrn( 6, 0) = c1x*35.0;
  double c2x = total_mass/A*Jx/210.0;
  mlTrn( 3, 3) = mlTrn( 9, 9) = c2x*70.0;
  mlTrn( 3, 9) = mlTrn( 9, 3) = c2x*35.0;
  double c1y = c1x/pow(1.0 + phiY, 2);
  mlTrn( 2, 2) = mlTrn( 8, 8) =  c1y*(70.0*phiY*phiY + 147.0*phiY + 78.0);
  mlTrn( 2, 8) = mlTrn( 8, 2) =  c1y*(35.0*phiY*phiY + 63.0*phiY + 27.0);
  mlTrn( 4, 4) = mlTrn(10,10) =  c1y*L*L/4.0*(7.0*phiY*phiY + 14.0*phiY + 8.0);
  mlTrn( 4,10) = mlTrn(10, 4) = -c1y*L*L/4.0*(7.0*phiY*phiY + 14.0*phiY + 6.0);
  mlTrn( 2, 4) = mlTrn( 4, 2) = -c1y*L/4.0*(35.0*phiY*phiY + 77.0*phiY + 44.0);
  mlTrn( 8,10) = mlTrn(10, 8) = -mlTrn( 2, 4);
  mlTrn( 2,10) = mlTrn(10, 2) = c1y*L/4.0*(35.0*phiY*phiY + 63.0*phiY + 26.0);
  mlTrn( 4, 8) = mlTrn( 8, 4) = -mlTrn( 2,10);

  double c2y = total_mass/A*Iy/(30.0*L*pow(1.0 + phiY, 2));
  mlRot( 2, 2) = mlRot( 8, 8) = c2y*36.0;
  mlRot( 2, 8) = mlRot( 8, 2) = -mlRot( 2, 2);
  mlRot( 4, 4) = mlRot(10,10) = c2y*L*L*(10.0*phiY*phiY + 5.0*phiY + 4.0);
  mlRot( 4,10) = mlRot(10, 4) = c2y*L*L*(5.0*phiY*phiY - 5.0*phiY - 1.0);
  mlRot( 2, 4) = mlRot( 4, 2) = mlRot( 2,10) = mlRot(10, 2) = c2y*L*(15.0*phiY - 3.0);
  mlRot( 4, 8) = mlRot( 8, 4) = mlRot( 8,10) = mlRot(10, 8) = -mlRot( 2, 4);

  double c1z = c1x/pow(1.0 + phiZ, 2);
  mlTrn( 1, 1) = mlTrn( 7, 7) =  c1z*(70.0*phiZ*phiZ + 147.0*phiZ + 78.0);
  mlTrn( 1, 7) = mlTrn( 7, 1) =  c1z*(35.0*phiZ*phiZ + 63.0*phiZ + 27.0);
  mlTrn( 5, 5) = mlTrn(11,11) =  c1z*L*L/4.0*(7.0*phiZ*phiZ + 14.0*phiZ + 8.0);
  mlTrn( 5,11) = mlTrn(11, 5) = -c1z*L*L/4.0*(7.0*phiZ*phiZ + 14.0*phiZ + 6.0);
  mlTrn( 1, 5) = mlTrn( 5, 1) =  c1z*L/4.0*(35.0*phiZ*phiZ + 77.0*phiZ + 44.0);
  mlTrn( 7,11) = mlTrn(11, 7) = -mlTrn( 1, 5);
  mlTrn( 1,11) = mlTrn(11, 1) = -c1z*L/4.0*(35.0*phiZ*phiZ + 63.0*phiZ + 26.0);
  mlTrn( 5, 7) = mlTrn( 7, 5) = -mlTrn( 1,11);

  double c2z = total_mass/A*Iz/(30.0*L*pow(1.0 + phiZ, 2));
  mlRot( 1, 1) = mlRot( 7, 7) = c2z*36.0;
  mlRot( 1, 7) = mlRot( 7, 1) = -mlRot( 1, 1);
  mlRot( 5, 5) = mlRot(11,11) = c2z*L*L*(10.0*phiZ*phiZ + 5.0*phiZ + 4.0);
  mlRot( 5,11) = mlRot(11, 5) = c2z*L*L*(5.0*phiZ*phiZ - 5.0*phiZ - 1.0);
  mlRot( 1, 5) = mlRot( 5, 1) = mlRot( 1,11) = mlRot(11, 1) = -c2z*L*(15.0*phiZ - 3.0);
  mlRot( 5, 7) = mlRot( 7, 5) = mlRot( 7,11) = mlRot(11, 7) = -mlRot( 1, 5);

  //
  // Add translational and rotational parts
  //
  M = mlTrn + mlRot;
}
