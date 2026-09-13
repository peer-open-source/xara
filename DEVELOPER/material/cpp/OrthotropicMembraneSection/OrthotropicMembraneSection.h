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
                                                                        
// $Revision: 1.9 $
// $Date: 2006-08-03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticMembranePlateSection.h,v $

// Ed "C++" Love
//
//  Elastic Plate Section with membrane
//

#ifndef OrthotropicMembranePlateSection_h
#define OrthotropicMembranePlateSection_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>


class OrthotropicMembraneSection : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    OrthotropicMembraneSection() ;

    //full constructor
    OrthotropicMembraneSection(int tag, double E, double E2, double v, double G, double h = 1.0, double rho = 0.0, double angle=0.0) ;

	//destructor
    ~OrthotropicMembraneSection() ;

    //make a clone of this material
    SectionForceDeformation *getCopy() ;


    const char *getClassType(void) const {return "ElasticMembranePlate";};

    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const ID& getType( ) ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //send back the initial tangent 
    const Matrix& getInitialTangent( ) ;

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    //density per unit area
    double getRho() ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);


  private :

    double E1 ;   // elastic modulus, material axis 1
	double E2 ;   // elastic modulus, material axis 2
	double v ;    // poisson ratio
    double G;     // shear ratio
    double h;     // MembranePlate thickness
    double rhoH ; // mass per unit 2D area
	double angle; // sine and cosine of the angle between material axes and local axes of the element 

	Matrix toMaterial, toLocal; // rotationMatrices from local to material axis (for strains) and from material to local (for stresses)
	// check the 1/2 factor on the strain transformation -> otherwise one is the transpoose of the other

    static const double five6 ; // =5/6 = shear correction factor

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;

    static ID array ;  

} ; 





#endif
