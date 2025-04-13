#pragma once
#include <TaggedObject.h>


template<typename MatT>
class PlaneSection: public TaggedObject {
  public:

    PlaneSection(int tag, MatT& material, double thickness) 
    : TaggedObject(tag), material(&material), thickness(thickness) 
    {
    }

    double getThickness() const
    {
      return thickness;
    }

    MatT* getMaterial() const
    {
      return material;
    }

  private:
    const double thickness;
    MatT* material;
};

