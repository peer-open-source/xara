///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           stress tensor with all necessery functions                #
//# CLASS:             stresstensor                                              #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              July 22 '93                                               #
//# UPDATE HISTORY:    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                    13 septembar '96 added reportAnim  :)                     #
//#                                                                              #
//################################################################################
//*/

#ifndef STRESSTENSOR_H
#define STRESSTENSOR_H

#include "BJtensor.h"
class Material_Model;
class OPS_Stream;

// namespace ucsd {

class stresstensor : public BJtensor
{
  public:
    friend class Material_Model;

  public:
    // just send appropriate arguments to the base constructor
    stresstensor(int rank_of_tensor=2, double initval=0.0); // default constructor
    stresstensor( double *values );
    stresstensor( double initvalue );

    stresstensor(const stresstensor & x );
    stresstensor(const BJtensor & x); // copy-initializer
    stresstensor(const nDarray & x); // copy-initializer

    //~stresstensor( );
    

    stresstensor operator=(const stresstensor & rval);// stresstensor assignment
    stresstensor operator=(const BJtensor & rval);// tensor assignment to stresstensor
    stresstensor operator=(const nDarray & rval);// nDarray assignment to stresstensor

    stresstensor deep_copy();

    //ini  // use "from" and initialize already allocated stress tensor from "from" values
    //ini      void Initialize( const stresstensor & from );

    //___// operator() overloading for 3D Gauss points!
    //___    stresstensor & operator()(short ir, short is, short it,
    //___                              short tr, short ts, short tt  );
    

    double Iinvariant1( ) const;
    double Iinvariant2( ) const;
    double Iinvariant3( ) const;

    double Jinvariant1( ) const;
    double Jinvariant2( ) const;
    double Jinvariant3( ) const;

    stresstensor deviator( ) const;
    stresstensor principal( ) const;

    double sigma_octahedral( ) const;
    double tau_octahedral( ) const;

    double ksi( )     const;
    double xi( )      const;
    double ro( )      const;
    double rho( )      const;
    double theta()   const;
    double thetaPI( ) const;

    double p_hydrostatic( ) const;
    double q_deviatoric( ) const;

    BJtensor dpoverds(  ) const;
    BJtensor dqoverds( void ) const;
    BJtensor dthetaoverds( void ) const;
    BJtensor d2poverds2( ) const;
    BJtensor d2qoverds2(  ) const;
    BJtensor d2thetaoverds2( ) const;
	     	          


    //--    stresstensor yield_surface_cross(stresstensor & end_stress,
    //--                                     Material_Model & YC);

    stresstensor pqtheta2stress( double, double, double );

    void report(const char *) const;
    void reportshort(const char *) const;
    void reportshortpqtheta(const char *) const;
    void reportSHORTpqtheta(const char *) const;
    void reportSHORTs1s2s3(const char *) const;
    void reportKLOTpqtheta(const char *) const;
    void reportshortI1J2J3(const char *) const;
    void reportAnim() const;
    void reportTensor(const char *) const;

    //================================================================================
    // Overloaded Insertion Operator	  ZHaohui Added Aug. 13, 2000
    // prints an stresstensor's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const stresstensor & rhs);

    //  // routine used by root finder, takes an alfa and returns the
    //  // yield function value for that alfa
    //    public:
    //      double func( stresstensor & start_stress,
    //                   stresstensor & end_stress,
    //                   Material_Model & YC,
    //                   double alfa );
    //  
    //  
    //  //..// polinomial root solver friend functions definitions
    //  //..public:
    //  //..friend void laguer(complex *, int , complex *, double , int );
    //  //..friend void zroots(complex *, int , complex *, int );
    //  //..
    //  
    // zero of function
    friend double zbrentstress(stresstensor   & start_stress,
                             stresstensor   & end_stress,
                             Material_Model & YC,
                             double x1, double x2, double tol);
  
    //  friend double zbrent(double x1, double x2, double tol);
    //  
    //  
};

// } // namespace ucsd

#endif

