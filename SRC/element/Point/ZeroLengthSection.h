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
// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the class definition for ZeroLengthSection.
// A ZeroLengthSection element is defined by two nodes with the same coordinate.
// A SectionForceDeformation object is associated with the nodes to
// provide the basic force-deformation relationship for the element.

#ifndef ZeroLengthSection_h
#define ZeroLengthSection_h

#include <Element.h>
#include <Matrix.h>
#include <Matrix3D.h>


class Node;
class Channel;
class FrameSection;
class Response;

namespace OpenSees {

class ZeroLengthSection : public Element
{
  public:
    
    ZeroLengthSection(int tag,                               
                      int dimension,
                      int Nd1, int Nd2, 
                      const Matrix3D& T,
                      FrameSection& theSection,
                      int doRayleighDamping = 0);

    ZeroLengthSection();
    ~ZeroLengthSection();

    const char *getClassType() const {return "ZeroLengthSection";}

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();

    int getNumDOF();        
    void setDomain(Domain *);

    // public methods to set the state of the element    
    int update();       // added by MSN to allow errors in setting section trial deformation
    int commitState();
    int revertToLastCommit();        
    int revertToStart();        

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff();
    const Matrix &getDamp();
    const Matrix &getInitialStiff();

    void zeroLoad();        
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);    

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();            

    // public methods for element output
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
    void Print(OPS_Stream &s, int flag);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &);
    
    int setParameter(const char **argv, int argc, Parameter &);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector &getResistingForceSensitivity(int gradIndex);
    int commitSensitivity(int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

  private:
    enum class Type: int {
      Truss   = 0,
      Euler2D = 1,
      Euler3D = 2,
      Shear2D = 3,
      Shear3D = 4,
    } type = Type::Truss;

    // private methods
    // void setUp (int Nd1, int Nd2, const Vector& x, const Vector& y);
    void setTransformation();
    void computeSectionDefs();

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 2 or 3 dimensions
    int numDOF;                                // number of dof for ZeroLengthSection
    Matrix3D transformation;                // transformation matrix for orientation
    int useRayleighDamping;
        
    Matrix *A;        // Transformation matrix ... e = A*(u2-u1)
    Vector *v;        // Section deformation vector, the element basic deformations
    
    Matrix *K;        // Pointer to element stiffness matrix
    Vector *P;        // Pointer to element force vector
    
    Node *theNodes[2];
    
    FrameSection *theSection;        // Pointer to section object
    int order;                // Order of the section model
    
    // Class wide matrices for return
    static Matrix K6;
    static Matrix K12;
    
    // Class wide vectors for return
    static Vector P6;
    static Vector P12;
};
} // namespace OpenSees
#endif




