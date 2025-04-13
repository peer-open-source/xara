#pragma once
#include <map>
#include <VectorND.h>
#include <TaggedObject.h>
#include <Vector3D.h>

#include <FrameTransform.h>
#include <LinearFrameTransf.hpp>
#include <SouzaFrameTransf.hpp>


class FrameTransformBuilder : public TaggedObject {
public:
    FrameTransformBuilder(int ndm, int t, const char *n) 
    : ndm(ndm), TaggedObject(t), vz{{0, 0, 0}}, offsets{}, offset_flags(0) {
      strncpy(name, n, 128);
    }
    virtual ~FrameTransformBuilder() {}
  
    template<int nn, int ndf>
    FrameTransform<nn, ndf> *
    create()
    {
      std::array<Vector3D, nn> *offset_array = nullptr;
      std::array<Vector3D, nn>  offset_data{};

      if (offsets.size() > 0) {
        offset_array = &offset_data;
        for (int i = 0; i < nn; ++i) {
          auto it = offsets.find(i+1);
          if (it != offsets.end())
            (*offset_array)[i] = it->second;
          
        }
      }
      int tag = this->getTag();
      if (strstr(name, "Linear") != nullptr)
        return new LinearFrameTransf<nn, ndf> (tag, vz, offset_array, offset_flags);
      else if (strstr(name, "Corot") != nullptr)
        return new SouzaFrameTransf<nn, ndf> (tag, vz, offset_array, offset_flags);
      else
        return nullptr;
    }

    virtual void 
    Print(OPS_Stream&s, int flag) {
        if (flag == OPS_PRINT_PRINTMODEL_JSON) {
            s << OPS_PRINT_JSON_MATE_INDENT << "{";
            s << "\"name\": " << this->getTag() << ", ";
            s << "\"type\": \"" << name << "\"";
            s << ", ";
            s << "\"vz\": [" << vz[0] << ", " << vz[1] << ", " << vz[2] << "]";
            if (offsets.size() > 0) {
                s << ", \"offsets\": {";
                for (auto it = offsets.begin(); it != offsets.end(); ++it) {
                    s << it->first << ": [" << it->second[0] << ", " << it->second[1] << ", " << it->second[2] << "]";
                    if (std::next(it) != offsets.end())
                        s << ", ";
                }
            }
            s << "}";
        }
    }
    int ndm;
    int tag;
    char name[128];
    Vector3D vz;
    std::map<int, Vector3D> offsets;
    int offset_flags;
};
