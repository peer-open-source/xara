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
// $Date: 2009-05-14 22:46:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSOE.h,v $

// Written: fmk
// Created: 05/09
//
// Description: This file contains the class definition for ArpackSOE


#ifndef ArpackSOE_h
#define ArpackSOE_h

#include "eigenSOE/EigenSOE.h"
#include <Vector.h>
#include <set>

class AnalysisModel;
class ArpackSolver;
class LinearSOE;


class ArpackSOE : public EigenSOE
{
  public:
    ArpackSOE(AnalysisModel&, double shift);

    ~ArpackSOE();

    int setLinks(AnalysisModel &);
    int setLinearSOE(LinearSOE &);

    int getNumEqn() const;
    int setSize(Graph &);
    
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addK(const Matrix &, const ID &, double fact = 1.0);
    int addM(const Matrix &, const ID &, double fact = 1.0); 
    int setMask();
    void opM(int n, double *v, double *result);

    void zeroA();
    void zeroM();

    double getShift();

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

    friend class ArpackSolver;

	int checkSameInt(int);

  protected:
    
  private:
    double *M;
    int Msize;
    bool mDiagonal;
    double shift;
    AnalysisModel *theModel;
    LinearSOE *theSOE;
    std::set<int> mask;


    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;
    ID *sizeLocal;
};


#endif



