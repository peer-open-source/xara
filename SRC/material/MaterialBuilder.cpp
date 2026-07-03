#include <MaterialBuilder.h>
#include <ContinuumUniaxial.h>

UniaxialMaterial* 
MaterialBuilder::getAxialCopy()
{
  return new ContinuumUniaxial(this->getTag(), *this);
}
