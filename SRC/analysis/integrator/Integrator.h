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
// Integrator encapsulates the interface provided to FE_Elements and
// DOF_Groups.
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
#ifndef Integrator_h
#define Integrator_h

class FE_Element;
class DOF_Group;
class Vector;
class ID;


class Integrator
{
public:
    virtual ~Integrator(); 

    virtual int formEleTangent(FE_Element *)  =0;
    virtual int formNodTangent(DOF_Group *)   =0;
    virtual int formEleResidual(FE_Element *) =0;
    virtual int formNodUnbalance(DOF_Group *) =0;

    // Methods provided for Domain Decomposition. This
    // is implemented in IncrementalIntegrator, but is
    // part of the Integrator interface for FE_Elements
    virtual int getLastResponse(Vector &result, const ID &id) =0;
};

#endif






