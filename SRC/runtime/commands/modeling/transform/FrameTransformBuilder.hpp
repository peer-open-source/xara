//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: Claudio M. Perez
//
#pragma once
#include <map>
#include <TaggedObject.h>
#include <VectorND.h>
#include <Vector3D.h>

#include <FrameTransform.h>
#include <LinearFrameTransf.h>
#include <SouzaFrameTransf.h>
#include <PDeltaFrameTransf3d.h>
#include <EuclidFrameTransf.h>
#include <Isometry/RankinIsometry.h>
#include <Isometry/CrisfieldIsometry.h>
#include <Isometry/BattiniIsometry.h>
#include <Isometry/LinearIsometry.h>
#if 0
#include <Isometry/SphericalIsometry.h>
#include <Isometry/IdentityIsometry.h>
#include <Isometry/Crisfield06.h>
#endif 

namespace OpenSees {

class FrameTransformBuilder : public TaggedObject {
public:
    FrameTransformBuilder(int ndm, int t, const char *n) 
    : ndm(ndm), 
      TaggedObject(t), 
      vz{{0, 0, 0}}, offsets{}, offset_flags(0) {
      snprintf(name, sizeof(name), "%s", n);
      // offset_flags |= LogIter;
    }

    virtual ~FrameTransformBuilder() {}
  
    template<int nn, int ndf>
    FrameTransform<nn, ndf> *
    create(int flags = 0)
    {
      int c_flags = offset_flags;
      if (flags)
        c_flags |= flags;

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
      if (strstr(name, "Linear") != nullptr) {
        if (!getenv("Linear02"))
          return new LinearFrameTransf<nn, ndf> (tag, vz, offset_array, c_flags);
        else
          return new EuclidFrameTransf<nn, ndf, LinearIsometry<nn>> (tag, vz, offset_array, c_flags);
      }

      else if (strstr(name, "LinearIsometric") != nullptr)
        return new EuclidFrameTransf<nn, ndf, LinearIsometry<nn>> (tag, vz, offset_array, c_flags);

      else if (strcmp(name, "Corotational") == 0) {
        if constexpr (ndf == 6 && nn==2)
          return new SouzaFrameTransf<nn, ndf> (tag, vz, offset_array, c_flags);
        else 
          return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,false>> (tag, vz, offset_array, c_flags);
      }

      else if (strstr(name, "PDelta") != nullptr) {
        bool ctan = false;
        if (strcmp(name, "PDelta02") == 0)
          ctan = true;
        if constexpr (nn == 2)
          return new PDeltaFrameTransf<nn, ndf> (tag, vz, offset_array, c_flags, ctan);
        else {
          opserr << OpenSees::PromptValueError 
                 << "PDelta formulation only implemented for 2-node elements."
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
      }

      else if (strcmp(name, "Corotational02") == 0 || strcmp(name, "Isometric") == 0 || strstr(name, "Rigid") != nullptr)
      {
        if (getenv("Battini"))
          return new EuclidFrameTransf<nn, ndf, BattiniIsometry<nn>> (tag, vz, offset_array, c_flags);
        else if (getenv("Crisfield"))
          return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,true>> (tag, vz, offset_array, c_flags);
        else if (getenv("Crisfield02"))
          return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,false>> (tag, vz, offset_array, c_flags);
        else
          return new EuclidFrameTransf<nn, ndf, RankinIsometry<nn>> (tag, vz, offset_array, c_flags);
      }
      else if (strcmp(name, "Corotational03") == 0)
        return new EuclidFrameTransf<nn, ndf, BattiniIsometry<nn>> (tag, vz, offset_array, c_flags);

      else if (strcmp(name, "Corotational04") == 0)
        return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,true>> (tag, vz, offset_array, c_flags);

      else if (strcmp(name, "Corotational05") == 0)
        return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,false>> (tag, vz, offset_array, c_flags);
#if 0
      else if (strcmp(name, "Corotational06") == 0)
        return new EuclidFrameTransf<nn, ndf, Crisfield06<nn>> (tag, vz, offset_array, c_flags);

      else if (strcmp(name, "Spherical") == 0)
        return new EuclidFrameTransf<nn, ndf, SphericalIsometry<nn>> (tag, vz, offset_array, c_flags);

      else if (strcmp(name, "Identity") == 0)
        return new EuclidFrameTransf<nn, ndf, IdentityIsometry<nn>> (tag, vz, offset_array, c_flags);
#endif
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

} // namespace OpenSees