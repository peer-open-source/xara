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
#pragma once
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <SectionForceDeformation.h>

namespace OpenSees {
enum ShellStress : int {
    // membrane
    Fxx = 11,
    Fyy = 12,
    Fxy = 13,
    // bending
    Mxx = 14,
    Myy = 15,
    Mxy = 16,
    // shear
    Vxz = 17,
    Vyz = 18,

    Max = 19
};

typedef int ShellStressLayout[Max-ShellStress::Fxx];

class ShellSection  {
public:
  ShellSection(SectionForceDeformation& s): section(s) {}
  ~ShellSection() 
  {
    delete &section;
  }

  ShellSection* getCopy() {
    return new ShellSection(*section.getCopy());
  }

  int getTag() const {
    return section.getTag();
  }
  int commitState() {
    return section.commitState();
  }
  int revertToLastCommit() {
    return section.revertToLastCommit();
  }
  int revertToStart() {
    return section.revertToStart();
  }

  template<int n, const ShellStressLayout& layout>
  int setTrialState(const VectorND<n>& deformation) {
    double e_data[Max-ShellStress::Fxx]{};
    Vector trialDeformation(e_data, Max-ShellStress::Fxx);
    const ID& s_layout = section.getType();
    for (int i=0; i<n; i++) {
      for (int j=0; j<(ShellStress::Max-ShellStress::Fxx); j++) {
        if (s_layout(j) == layout[i])
          e_data[j] = deformation(i);
      }
    }
    return section.setTrialSectionDeformation(trialDeformation);
  }

  template<int n, const ShellStressLayout& layout>
  VectorND<n> getResultant() {
    VectorND<n> S{};
    const ID& s_layout = section.getType();
    int m = section.getOrder();
    const Vector& s = section.getStressResultant();
    for (int i=0; i<n; i++) {
      for (int j=0; j<m; j++) {
        if (s_layout(j) == layout[i])
          S(i) = s(j);
      }
    }
    return S;
  }

  template<int n, const ShellStressLayout& layout>
  MatrixND<n,n> getTangent() {
      MatrixND<n,n> Ks;
      const Matrix& ks = section.getSectionTangent();
      const ID& s_layout = section.getType();
      int m = section.getOrder();

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Ks(i,j) = 0.0;
          
          // Find where layout[i] is in the section's layout (Row)
          int row_idx = -1;
          for(int r = 0; r < m; r++) {
            if (s_layout(r) == layout[i]) { row_idx = r; break; }
          }

          // Find where layout[j] is in the section's layout (Column)
          int col_idx = -1;
          for(int c = 0; c < m; c++) {
            if (s_layout(c) == layout[j]) { col_idx = c; break; }
          }

          // Map the stiffness if both components exist in the underlying section
          if (row_idx != -1 && col_idx != -1) {
            Ks(i,j) = ks(row_idx, col_idx);
          }
        }
      }
      return Ks;
  }

  private:
    SectionForceDeformation& section;
};

} // namespace OpenSees