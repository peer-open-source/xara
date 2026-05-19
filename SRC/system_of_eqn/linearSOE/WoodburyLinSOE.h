/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#ifndef WoodburyLinSOE_h
#define WoodburyLinSOE_h

#include <LinearSOE.h>

class IncrementalIntegrator;
class AnalysisModel;
class WoodburySolver;
class Domain;
class DirectIntegrationAnalysis;

// Wraps a structural LinearSOE for modalDampingW: A_s in inner matrix; A_modal = Q*diag(dbar)*Q^T
// applied in solve (Woodbury) and formAp. See WoodburyLinSOE.cpp for formulas.
class WoodburyLinSOE : public LinearSOE
{
  public:
    explicit WoodburyLinSOE(LinearSOE &innerSOE, bool takeOwnership = false);
    ~WoodburyLinSOE() override;

    LinearSOE &getInnerSOE(void) const { return *innerSOE; }

    // Relinquish ownership of the inner SOE (e.g. before destroying the wrapper).
    LinearSOE *releaseInner(void);

    // Build Z = A_s^{-1}*Q and G = Q^T*Z + diag(1/dbar) for current A_s.
    int formWoodburyBasis(IncrementalIntegrator &integrator,
                          const Vector *modalFactors);

    // x = x_s - Z * G^{-1} * Q^T * x_s  (inner must hold x_s = A_s^{-1} b)
    int applyWoodburyCorrection(void);

    int setSize(Graph &theGraph) override;
    int getNumEqn(void) const override;

    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addA(const Matrix &) override;
    int addB(const Vector &, const ID &, double fact = 1.0) override;
    int setB(const Vector &, double fact = 1.0) override;

    void zeroA(void) override;
    void zeroB(void) override;

    int formAp(const Vector &p, Vector &Ap) override;

    const Vector &getX(void) override;
    const Vector &getB(void) override;
    const Matrix *getA(void) override;
    double normRHS(void) override;

    void setX(int loc, double value) override;
    void setX(const Vector &x) override;

    int setLinks(AnalysisModel &theModel) override;

    int solve(void) override;
    LinearSOESolver *getSolver(void);

    int addColA(const Vector &col, int colIndex, double fact = 1.0) override;

    int saveSparseA(OPS_Stream &output, int baseIndex = 0) override;
    int getSparseA(ID &rowIndices, ID &colIndices, Vector &values,
                   int baseIndex = 0) override;
    int getSparseA(std::vector<int> &rowIndices, std::vector<int> &colIndices,
                   std::vector<double> &values, int baseIndex = 0) override;

    double getDeterminant(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) override;

    void hookWoodburySolver(void);
    void unhookWoodburySolver(void);

  private:
    void clearWoodburyBasis(void);
    int applyModalMatvec(const Vector &p, Vector &Ap);
    void resizeX(int size);

    LinearSOE *innerSOE;
    Vector *vectX;
    bool ownsInner;

    int numDOF;
    int numModes;

    double *eigenVectors;  // Q  (numDOF x numModes, column-major)
    double *dBar;          // diag(dbar)

    Matrix *woodburyZ;     // Z = A_s^{-1} * Q
    Matrix *woodburyG;     // G = Q^T * Z + diag(1/dbar)

    Vector *workV1;        // scratch (modal-sized)
    Vector *workV2;        // y = G^{-1} * Q^T * x

    WoodburySolver *woodburySolver;

    bool basisValid;
};

// Decorator: wrap *soePtr in WoodburyLinSOE and attach to transient when modalDampingW is active.
// No-op if theDomain is set and not MODAL_DAMPING_WOODBURY, or if transient is not defined yet.
int wrapLinearSOEWithWoodbury(LinearSOE *&soePtr,
                              Domain *theDomain = nullptr,
                              DirectIntegrationAnalysis *transient = nullptr);

#endif
