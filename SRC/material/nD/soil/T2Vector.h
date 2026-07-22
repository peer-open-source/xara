//<<<<<<< T2Vector.h
//<<<<<<< T2Vector.h
// $Revision: 1.8 $
// $Date: 2009-07-23 23:57:27 $
//=======
// $Revision: 1.8 $
// $Date: 2009-07-23 23:57:27 $
//>>>>>>> 1.4
//=======
// $Revision: 1.8 $
// $Date: 2009-07-23 23:57:27 $
//>>>>>>> 1.6
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/T2Vector.h,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// T2Vector.h
// --------
//

#ifndef _T2Vector_H_
#define _T2Vector_H_

#include <Vector.h>
#include <Channel.h>
#include <float.h>
#include <VectorND.h>

#define UP_LIMIT    1.0e+30
#define LOW_LIMIT   20.*DBL_EPSILON
#define DYNAMIC_T2VECTOR
// global function: scalar product of two second order tensor vectors

namespace OpenSees {

double operator && (const Vector &, const Vector &);
void doubledotProduct (Vector & c, const Vector & a, const Matrix & b);
void tensorProduct(Matrix & c, const Vector & a, const Vector & b);
void doubledotMatrixProduct (Matrix & c, const Matrix & a, const Matrix & b);

// define second order tensor vector class
class T2Vector 
{

public:
  // constructors
  T2Vector();
  T2Vector(const Vector & T2Vector_init, int isEngrgStrain=0);
  T2Vector(const Vector & deviat_init, double volume_init);
  
  ~T2Vector();
  enum class Basis: int {Stress=0, Strain=1};

  void setData(const Vector &init, int isEngrgStrain =0);
  void setData(const Vector &deviat, double volume);
  void setData(const VectorND<6> & init, Basis);

  const Vector& t2Vector(int isEngrgStrain=0) const {return this->t2Vector(static_cast<Basis>(isEngrgStrain));}
  const Vector& deviator(int isEngrgStrain=0) const {return this->deviator(static_cast<Basis>(isEngrgStrain));}
  const Vector& t2Vector(Basis) const;
  const Vector& deviator(Basis) const;

  double volume() const {return theVolume; }
  const Vector& unitT2Vector() const;
  const Vector& unitDeviator() const;
  double t2VectorLength() const;
  double deviatorLength() const;
  double octahedralShear(int isEngrgStrain=0) const;

  // = -sqrt(3/2*(S:S))/(p+residualPress)
  double deviatorRatio(double residualPress=0.) const; 

  //next function return the angle between two T2Vectors in radians (-PI to PI)
  double angleBetweenT2Vector(const T2Vector &) const; 

  //next function return the angle between deviatoric components of
  //two vectors in radians (-PI to PI)
  double angleBetweenDeviator(const T2Vector &) const; 

  int operator == (const T2Vector & a) const;
  int isZero() const;
  int Zero();

private:
#ifdef DYNAMIC_T2VECTOR
  Vector theT2Vector;
  Vector theDeviator;
#else
  VectorND<6> theT2Vector;
  VectorND<6> theDeviator;
#endif
  double theVolume;
  static Vector engrgStrain;
};

} // namespace OpenSees

#endif
