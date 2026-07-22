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
// Description: This file contains global variables used in OpenSees files.
// if you change some of the variables, you must recompile ALL the code.
//
// Written: fmk
// Created: 11/99
//
#ifndef _OPS_Globals_h
#define _OPS_Globals_h

#define _USING_OpenSees_STREAMS
#include <handler/OPS_Stream.h>
#include <logging/Logging.h>
// namespace OpenSees {
// extern OPS_Stream *opserrPtr;
// }
// #ifndef opserr
// # define opserr (*OpenSees::opserrPtr)
// #endif
// #define endln "\n"

#include <string.h>


extern double   ops_Dt;                // current delta T for current domain doing an update
extern int ops_Creep;
class Domain;
class Element;
extern Domain  *ops_TheActiveDomain;   // current domain undergoing an update
extern Element *ops_TheActiveElement;  // current element undergoing an update

// global variable for initial state analysis
// added: Chris McGann, University of Washington
extern bool  ops_InitialStateAnalysis;

#endif
