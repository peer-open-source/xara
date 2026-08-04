#pragma once
#include <VectorND.h>
#include <MatrixND.h>

namespace OpenSees {

class SolidMaterial
{
  public:
    SolidMaterial(int tag) : tag(tag) {}
    virtual ~SolidMaterial() {}

    virtual int setTrialStrain(const VectorND<6> &e, VectorND<6> &s, MatrixND<6,6> &dsde) = 0;

    virtual int commit() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;

    int getTag() const { return tag; }

private:
    int tag;
};

}
