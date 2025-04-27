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
                                                                        
// $Revision: 1.3 $
// $Date: 2009-10-13 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStrainMaterial.h,v $

// Antonios Vytiniotis
//
// Generic Plane Strain Material
//


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>



class PlaneStrainMaterial: public NDMaterial{
  public : 
    PlaneStrainMaterial( ) ;
    PlaneStrainMaterial(int tag, 
			NDMaterial &the3DMaterial ) ;

    virtual ~PlaneStrainMaterial( ) ;

        NDMaterial *getCopy( ) ;
    NDMaterial *getCopy( const char *type ) ;

        int getOrder( ) const ;

        const char *getType( ) const ;

    //swap history variables
    int commitState( ) ; 

        int revertToLastCommit( ) ;

        int revertToStart( ) ;

    //get the strain 
    int setTrialStrain( const Vector &strainFromElement ) ;

    const Vector& getStrain( ) ;

    const Vector& getStress( ) ;

    const Matrix& getTangent( ) ;

    const Matrix& getInitialTangent( ) ;  // AV Not Sure if it works

    //density
    double getRho( ) ;

    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    int setParameter(const char **argv, int argc, Parameter &param);
    Response* setResponse(const char** argv, int argc, OPS_Stream& s);

private :
    NDMaterial *theMaterial ;  //pointer to three dimensional material

    Vector strain ;
    static Vector stress ;
    static Matrix tangent ;

} ;





