#pragma once

#include <cassert>
#include <TaggedObject.h>
#include <MovableObject.h>


namespace OpenSees {
class FrameMaterial;
}
class NDMaterial;
class UniaxialMaterial;

using namespace OpenSees;

class MaterialBuilder: public TaggedObject
{
  public:
    MaterialBuilder(int tag): TaggedObject(tag) {};
    virtual ~MaterialBuilder() {};

    virtual NDMaterial *getCopy(const char *type) =0;
    virtual FrameMaterial* getFrameFiber() {return nullptr;}

    virtual UniaxialMaterial* getAxialCopy();
};
