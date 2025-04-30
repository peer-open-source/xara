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
// Description: This file contains the class definition for TaggedObject.
// A TaggedObject is an object with an integer identifier. It is used as
// a base class by DomainComponent, Graph and other classes in the framework.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#ifndef TaggedObject_h
#define TaggedObject_h
//
#define OPS_PRINT_CURRENTSTATE 0
#define OPS_PRINT_PRINTMODEL_SECTION  1
#define OPS_PRINT_PRINTMODEL_MATERIAL 2
#define OPS_PRINT_PRINTMODEL_JSON   25000

#define OPS_PRINT_JSON_ELEM_INDENT "       "
#define OPS_PRINT_JSON_NODE_INDENT "       "
#define OPS_PRINT_JSON_MATE_INDENT "       "

class OPS_Stream;

class TaggedObject 
{
  public:
    TaggedObject(int tag);
    virtual ~TaggedObject();

    inline int getTag() const;

    virtual void inline Print(OPS_Stream &, [[maybe_unused]] int flag =0) {}

    friend OPS_Stream &operator<<(OPS_Stream &s, TaggedObject &m);        

  protected:
    void setTag(int newTag);  // CAUTION: this is a dangerous method to call
    
    constexpr static char JsonGeometryIndent[] = "       ";
    constexpr static char JsonPropertyIndent[] = "       ";

  private:    
    int theTag;    
};

inline int 
TaggedObject::getTag() const
{
    return theTag;
}

#endif

