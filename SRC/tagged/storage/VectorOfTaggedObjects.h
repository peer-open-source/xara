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
// $Revision: 1.2 $
// $Date: 2003-02-14 23:02:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/VectorOfTaggedObjects.h,v $
// File: ~/tagged/storage/VectorOfTaggedObjects.h
// 
// Description: This file contains the class definition for 
// VectorOfTaggedObjects. 
// VectorOfTaggedObjects is a storage class. The class 
// is responsible for holding and providing access to objects of type 
// TaggedObject. 
// A std::map is used to store
// the pointers to these objects.
//
// Written: fmk 
// Created: 02/00
// Revision: A
//
#ifndef VectorOfTaggedObjects_h
#define VectorOfTaggedObjects_h
//
#include <vector>
#include <TaggedObjectStorage.h>
#include <TaggedObjectIter.h>

class VectorOfTaggedObjects : public TaggedObjectStorage
{
  public:
    VectorOfTaggedObjects();
    ~VectorOfTaggedObjects();    

    class Iterator: public TaggedObjectIter
    {
      public:
        Iterator(VectorOfTaggedObjects &theComponents);
        virtual ~Iterator();
        
        void reset() final;
        TaggedObject *operator()() final;
        
      private:
        std::vector<TaggedObject *> *container;
        std::vector<TaggedObject *>::iterator currentComponent;
    };

    // public methods to populate a domain
    int  setSize(int newSize);
    bool addComponent(TaggedObject *newComponent) final;
    TaggedObject *removeComponent(int tag) final;    
    int getNumComponents() const;
    
    TaggedObject     *getComponentPtr(int tag) final;
    TaggedObjectIter &getComponents() final {
      return this->getIterRef();
    }
    Iterator& getIterRef();
    Iterator getIter();
  
    TaggedObjectStorage *getEmptyCopy();
    void clearAll(bool invokeDestructor = true);
    
    void Print(OPS_Stream &s, int flag);

    
  private:
    std::vector<TaggedObject *> container; // the map container for storing the pointers
    Iterator  myIter;  // the iter for this object
};

#endif

