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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/EigenIntegrator.h,v $

// Description: This file contains the class definition of EigenIntegrator.
// EigenIntegrator is an algorithmic class for setting up the finite element 
// equations for a eigen problem analysis. 
//
// This class is inheritanted from the base class of Integrator which was
// created by fmk (Frank).
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//

#ifndef EigenIntegrator_h
#define EigenIntegrator_h

#include <Integrator.h>
#include <MovableObject.h>

class EigenSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;
class OPS_Stream;

class EigenIntegrator : public Integrator
{
  public:
     EigenIntegrator(AnalysisModel&, EigenSOE &);
     virtual ~EigenIntegrator();

     // methods to form the M and K matrices.
     int formK();
     int formM();
     int getLastResponse(Vector &result, const ID &id);


     // methods to instruct the FE_Element and DOF_Group objects
     // how to determining their contribution to M and K
     virtual int formEleTangK(FE_Element *);
     virtual int formEleTangM(FE_Element *);
     virtual int formNodTangM(DOF_Group *);

     int formEleTangent(FE_Element *) override;
     int formNodTangent(DOF_Group *) override;
     int formEleResidual(FE_Element *) override;
     int formNodUnbalance(DOF_Group *) override;

  
 private:
    EigenSOE *theSOE;
    AnalysisModel *theAnalysisModel;
    int flagK;
};

#endif




