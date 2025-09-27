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
#include <SolutionAlgorithm.h>
#include <Recorder.h>
#include <stdlib.h>
#include <Logging.h>

// TODO!!!!
int SOLUTION_ALGORITHM_tangentFlag = 0;

SolutionAlgorithm::SolutionAlgorithm(int clasTag)
: MovableObject(clasTag), theRecorders(0), numRecorders(0)
{

}

SolutionAlgorithm::~SolutionAlgorithm()
{
    for (int i=0; i<numRecorders; i++)
        delete theRecorders[i];
    
    if (theRecorders != nullptr) {
        free((void *)theRecorders);    
    }
}

int  
SolutionAlgorithm::addRecorder(Recorder &theRecorder)
{
    Recorder **newRecorders = (Recorder **)malloc((numRecorders+1)*sizeof(Recorder *));
    
    for (int i=0; i<numRecorders; i++)
        newRecorders[i] = theRecorders[i];
    newRecorders[numRecorders] = &theRecorder;

    if (theRecorders != nullptr)
      free((void *)theRecorders);

    theRecorders = newRecorders;
    numRecorders++;
    return 0;
}


int  
SolutionAlgorithm::record(int cTag)
{
    for (int i=0; i<numRecorders; i++)
        theRecorders[i]->record(cTag, 0.0);
    return 0;
}

