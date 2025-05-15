#pragma once 
#include <VectorND.h>
#include <MatrixND.h>
#include <NDMaterial.h>

namespace OpenSees {

template <int NEN, int NIP, int NST>
void StressExtrapolation(NDMaterial**materials, const MatrixND<NEN,NIP>& We,
                         VectorND<NST*NEN>& stressAtNodes)
{
    // extrapolate stress from Gauss points to element nodes
    VectorND<NST*NIP> stressGP;

    // first get stress components at Gauss points
    int cnt = 0;
    for (int i = 0; i < NIP; i++) {
      // Get material stress response
      const Vector &sigma = materials[i]->getStress();
      for (int j=0; j<NST; j++) {
        stressGP(cnt+j) = sigma(j);
      }
      cnt += NST;
    }

    for (int i = 0; i < NEN; i++) {
      for (int k = 0; k < NST; k++) { // number of stress components
        int p = NST*i + k;
        for (int j = 0; j < NIP; j++) { // nip
            stressAtNodes(p) += We(i,j) * stressGP(NST*j + k);
        }
      }
    }
}

} // namespace OpenSees