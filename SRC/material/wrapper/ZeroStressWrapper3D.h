#ifndef ZeroStressWrapper3D_h
#define ZeroStressWrapper3D_h

#include <NDMaterial.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Stream.h>

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace OpenSees {

template <const std::array<bool, 6>& IsStressZero, const std::array<int, 6>& order>>
class ZeroStressWrapper3D : public NDMaterial {
public:
    ZeroStressWrapper3D(int tag, const NDMaterial& to_wrap);
    ~ZeroStressWrapper3D();

    int setTrialStrain(const Vector &strainFromElement) override;

    const Matrix& getTangent() override;
    const Matrix& getInitialTangent() override;
    const Vector& getStress() override;
    const Vector& getStrain() override;

    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;

    NDMaterial* getCopy() override;
    NDMaterial* getCopy(const char *type) override;

    const char* getType() const override;
    int getOrder() const override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    void Print(OPS_Stream &s, int flag = 0) override;

private:
    //
    void initializeIndices() {
        idx_c.clear();
        idx_r.clear();
        for (int i = 0; i < N_TOTAL_; ++i) {
            if (IsStressZero[i]) {
            idx_c.push_back(i);
            } else {
            idx_r.push_back(i);
            }
        }
    }

    int enforceZeroStressConditions(const VectorND<6>& eps);
    int computeWrapperTangentAndStress(); // From wrapped material's current state
    int computeCondensedInitialTangent(); // Computes and stores condensed initial tangent

    // --- Member Variables ---
    NDMaterial* w_material = nullptr;

    // Internal state using fixed-size VectorND and MatrixND
    VectorND<6> eps_; // Current total strain [eps_c, eps_r] satisfying sigma_c=0
    VectorND<6> sig_; // Current stress [sigma_c=0, sigma_r]
    MatrixND<6,6> C_; // Current condensed tangent, expanded to 6x6

    struct {
        VectorND<6> eps; // strain
        VectorND<6> sig; // stress
        MatrixND<6,6> C; // tangent
    } init, past, pres;

    struct {
        VectorND<6> eps; // strain
        VectorND<6> sig; // stress
        MatrixND<6,6> C; // tangent
    } wrap_past, wrap_pres;

    // For initial tangent
    MatrixND<6,6> initialTangentFromWrappedMaterial_nd_;
    MatrixND<6,6> initialTangentCondensed_nd_; // Stored condensed initial tangent
    bool initialTangentComputed_ = false;

    // Return buffers for NDMaterial interface
    Vector retStrain;
    Vector retStress;
    Matrix retTangent;
    Matrix retInitialTangent;

    // Indices mapping
    std::vector<int> idx_c; // Global indices (0-5) for sigma_c (condensed)
    std::vector<int> idx_r; // Global indices (0-5) for sigma_r (retained)

    // Compile-time constants for dimensions
    static constexpr int N_TOTAL_ = 6;
    static constexpr int countBooleanArray(const std::array<bool, N_TOTAL_>& arr) {
        int count = 0;
        for (std::size_t i = 0; i < arr.size(); ++i) {
            if (arr[i])
              count++;
        }
        return count;
    }

    // NC = number of condensed components (stress is zero)
    // NR = number of retained components (stress is non-zero)
    static constexpr int NC = countBooleanArray(IsStressZero); 
    static constexpr int NR = N_TOTAL_ - NC;                  

    // Iteration control for internal Newton loop
    static constexpr int MAX_ITERATIONS_ = 30;
    static constexpr double CONVERGENCE_TOL_STRESS_ = 1.0e-9; // Absolute stress norm tolerance
};
} // namespace OpenSees

#include "ZeroStressWrapper3D.cpp" 

#endif // ZeroStressWrapper3D_h