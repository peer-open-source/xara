// $Revision: 1.3 $
// $Date: 2009-05-14 22:45:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/EigenSOE.h,v $

// Written: Jun Peng
// Created: Sat Feb. 6, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenSOE.
// EigenSOE is a subclass of SystemOfEqn.
// It has pure virtual functions which must be implemented in it's derived
// subclasses. To solve the genreal eigen value equations means that
// by the given K and M, find the corresponding eigen value and eigen
// vectors.
//


#ifndef EigenSOE_h
#define EigenSOE_h

#include <MovableObject.h>

class EigenSolver;
class AnalysisModel;
class Graph;
class Matrix;
class Vector;
class ID;
class LinearSOE;


enum {
 EigenSOE_TAGS_BandArpackSOE 	=1,
 EigenSOE_TAGS_SymArpackSOE 	=2,
 EigenSOE_TAGS_SymBandEigenSOE  =3,
 EigenSOE_TAGS_FullGenEigenSOE  =4,
 EigenSOE_TAGS_ArpackSOE,
 EigenSOE_TAGS_GeneralArpackSOE,
 EigenSOE_TAGS_SymmGeneralizedEigenSOE,

 EigenSOLVER_TAGS_BandArpackSolver,
 EigenSOLVER_TAGS_SymArpackSolver ,
 EigenSOLVER_TAGS_SymBandEigenSolver,
 EigenSOLVER_TAGS_FullGenEigenSolver,
 EigenSOLVER_TAGS_ArpackSolver,
 EigenSOLVER_TAGS_GeneralArpackSolver,
 EigenSOLVER_TAGS_SymmGeneralizedEigenSolver
};

enum {
 EigenALGORITHM_TAGS_Frequency =1,
 EigenALGORITHM_TAGS_Standard  =2
};

class EigenSOE : public MovableObject
{
  public:
     EigenSOE(EigenSolver &, int classTag);
     EigenSOE(int classTag);
     virtual ~EigenSOE();
     
     virtual int solve(int numModes, bool generalized, bool findSmallest = true);
     virtual int setLinks(AnalysisModel &);    
     virtual int setLinearSOE(LinearSOE &) {return -1;};

     // pure virtual functions
     virtual int addA(const Matrix &, const ID &, double fact = 1.0) = 0;
     virtual int addM(const Matrix &, const ID &, double fact = 1.0) = 0;
     virtual void zeroA() = 0;
     virtual void zeroM() = 0;

     virtual int setSize(Graph &) = 0;
     int formSystem(AnalysisModel&, bool generalized);

     // methods to get the eigenvectors and eigenvalues
     virtual const Vector &getEigenvector(int mode);
     virtual double getEigenvalue(int mode);          

  protected:
     virtual int setSolver(EigenSolver &);
     EigenSolver *getSolver();
     
  private:
     EigenSolver *theSolver;
     
};

#endif



