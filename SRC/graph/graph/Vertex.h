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
#ifndef Vertex_h
#define Vertex_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for Vertex.
// Vertex is an element of a graph.
//
#include <TaggedObject.h>
#include <ID.h>

#define START_VERTEX_NUM 0
class Channel;
class FEM_ObjectBroker;

#define VIRTUAL

class Vertex: public TaggedObject
{
  public:
    Vertex(int tag, int ref, double weight=0, int color =0);
    Vertex(const Vertex &other);

    VIRTUAL ~Vertex();

    VIRTUAL void setWeight(double newWeight);
    VIRTUAL void setColor(int newColor);
    VIRTUAL void setTmp(int newTmp);    
    
    VIRTUAL int getRef(void) const;
    VIRTUAL double getWeight(void) const;
    VIRTUAL int getColor(void) const;
    VIRTUAL int getTmp(void) const;    

    VIRTUAL int addEdge(int otherTag);
    VIRTUAL int getDegree(void) const;
    VIRTUAL const ID &getAdjacency(void) const;
    VIRTUAL void setAdjacency(const ID& adj) {myAdjacency = adj;}

    VIRTUAL  void Print(OPS_Stream &s, int flag =0);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  private:
    int myRef;
    double myWeight;
    int myColor;
    int myDegree;
    int myTmp;
    ID  myAdjacency;
};

#endif

