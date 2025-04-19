/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for Vector.
// Vector is a concrete class implementing the vector abstraction.
//
// Written: fmk 
// Created: 11/96
//
#ifndef Vector_h
#define Vector_h 

#include <memory>
#include <assert.h>

#define VECTOR_VERY_LARGE_VALUE 1.0e200

class Matrix; 
class Message;
class SystemOfEqn;
class OPS_Stream;
class ID;
namespace OpenSees {
  template<int n, typename T> struct VectorND;
  template<int nr, int nc, typename T> struct MatrixND;
};

class Vector
{
  public:
    // constructors and destructor
    Vector();
    Vector(int);
    Vector(const Vector &);
    // Referencing constructors; these dont own the data
    Vector(const Vector &, int start, int size);
    Vector(double *data, int size);
    Vector(std::shared_ptr<double[]>, int size);
#if !defined(NO_CXX11_MOVE)
    Vector(Vector &&);    
#endif
    template <int n> Vector(OpenSees::VectorND<n,double>& v)
      : sz(v.size()), theData(v.values), fromFree(1)
    {
    }

    ~Vector();

    template <int n> operator OpenSees::VectorND<n,double>() const {
      OpenSees::VectorND<n,double> ret;
      ret = *this; 
      return ret;
    }

    // utility methods
    int setData(double *newData, int size);
#if 1
    template <int n>
    inline int setData(OpenSees::VectorND<n, double> v) {
      return setData(&v.values[0], n);
    }
#endif
    Vector view(int start, int end) const;
    int Assemble(const Vector &V, const ID &l, double fact = 1.0);
    double Norm() const;
    double pNorm(int p) const;
    inline int Size() const;
    int resize(int newSize);
    inline void Zero();
    int Normalize();
    
    int addVector(const Vector &other, double factOther);
    int addVector(double factThis, const Vector &other, double factOther);
    int addMatrixVector(double factThis, const Matrix &m, const Vector &v, double factOther); 
    int addMatrixTransposeVector(double factThis, const Matrix &m, const Vector &v, double factOther);

    // overloaded operators
    inline double operator()(int x) const;
    inline double &operator()(int x);
    inline double operator[](int x) const;  
    double &operator[](int x); // this operator does bounds checks
    Vector operator()(const ID &rows) const;
    Vector &operator=(const Vector  &V);
#if !defined(NO_CXX11_MOVE)   
    Vector &operator=(Vector  &&V);
#endif
    Vector &operator+=(double fact);
    Vector &operator-=(double fact);
    Vector &operator*=(double fact);
    Vector &operator/=(double fact); 

    Vector operator+(double fact) const;
    Vector operator-(double fact) const;
    Vector operator*(double fact) const;
    Vector operator/(double fact) const;

    Vector &operator+=(const Vector &V);
    Vector &operator-=(const Vector &V);
    
    Vector operator+(const Vector &V) const;
    Vector operator-(const Vector &V) const;
    double operator^(const Vector &V) const;
    Vector operator/(const Matrix &M) const;

    int operator==(const Vector &V) const;
    int operator==(double) const;
    int operator!=(const Vector &V) const;
    int operator!=(double) const;

    //operator added by Manish @ UB
    Matrix operator%(const Vector &V) const;

    // methods added by Remo
    int  Assemble(const Vector &V, int init_row, double fact = 1.0);
    int  Extract (const Vector &V, int init_row, double fact = 1.0); 
  
    friend OPS_Stream &operator<<(OPS_Stream &s, const Vector &V);
    // friend istream &operator>>(istream &s, Vector &V);    
    friend Vector operator*(double a, const Vector &V);
    
    friend class Message;
    friend class SystemOfEqn;
    friend class Matrix;
//  template<int n> friend struct OpenSees::VectorND;
    template<int n, typename> friend struct OpenSees::VectorND;
    template<int nr, int nc, typename> friend struct OpenSees::MatrixND;
    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;    
    friend class MPI_Channel;
    friend class MySqlDatastore;
    friend class BerkeleyDbDatastore;
    
  private:
    int sz;
    double *theData;
    int fromFree;
};


/********* INLINED VECTOR FUNCTIONS ***********/
inline int 
Vector::Size() const 
{
  return sz;
}


inline void
Vector::Zero() {
  for (int i=0; i<sz; i++)
    theData[i] = 0.0;
}


inline double 
Vector::operator()(int x) const
{
  assert(x >= 0 && x < sz);
  return theData[x];
}

inline double
Vector::operator[](int x) const
{
  // check if it is inside range [0,sz-1]
  assert(x >= 0 && x < sz);
  return theData[x];
}


inline double &
Vector::operator()(int x)
{
  assert(x >= 0 && x < sz);
  return theData[x];
}

#endif

