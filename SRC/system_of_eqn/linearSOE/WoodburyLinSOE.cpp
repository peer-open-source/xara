/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#include <WoodburyLinSOE.h>
#include <WoodburySolver.h>
#include <IncrementalIntegrator.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DirectIntegrationAnalysis.h>
#include <Matrix.h>
#include <Vector.h>
#include <Graph.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <classTags.h>
#include <cstring>
#include <math.h>

// modalDampingW: A x = b with  A = A_s + Q D Q^T  (A_s assembled in inner SOE only).
//
//   Q = M Phi   (mass-weighted mode shapes; columns with xi_i = 0 omitted)
//   D = diag(dbar),   dbar_i = c_2 * 2*xi_i*omega_i,   omega_i = sqrt(lambda_i)
//
//   x_s = A_s^{-1} b
//   Z = A_s^{-1} Q
//   G = D^{-1} + Q^T Z  
//   x = x_s - Z * G^{-1} * Q^T * x_s
//
// formAp: A*p = A_s*p + Q*D*Q^T*p.  Z,G rebuilt when A_s changes (zeroA/addA).

WoodburyLinSOE::WoodburyLinSOE(LinearSOE &inner, bool takeOwnership)
    : LinearSOE(LinSOE_TAGS_WoodburyLinSOE),
      innerSOE(&inner),
      ownsInner(takeOwnership),
      numDOF(0),
      numModes(0),
      eigenVectors(nullptr),
      dBar(nullptr),
      woodburyZ(nullptr),
      woodburyG(nullptr),
      workV1(nullptr),
      workV2(nullptr),
      woodburySolver(nullptr),
      basisValid(false),
      vectX(nullptr)
{
    hookWoodburySolver();
}

WoodburyLinSOE::~WoodburyLinSOE()
{
    unhookWoodburySolver();
    clearWoodburyBasis();
    if (vectX != nullptr) {
        delete vectX;
        vectX = nullptr;
    }
    if (ownsInner)
        delete innerSOE;
}

LinearSOE *
WoodburyLinSOE::releaseInner(void)
{
    LinearSOE *inner = innerSOE;
    innerSOE = nullptr;
    ownsInner = false;
    return inner;
}

void
WoodburyLinSOE::clearWoodburyBasis(void)
{
    if (eigenVectors != nullptr) {
        delete[] eigenVectors;
        eigenVectors = nullptr;
    }
    if (dBar != nullptr) {
        delete[] dBar;
        dBar = nullptr;
    }
    if (woodburyZ != nullptr) {
        delete woodburyZ;
        woodburyZ = nullptr;
    }
    if (woodburyG != nullptr) {
        delete woodburyG;
        woodburyG = nullptr;
    }
    if (workV1 != nullptr) {
        delete workV1;
        workV1 = nullptr;
    }
    if (workV2 != nullptr) {
        delete workV2;
        workV2 = nullptr;
    }
    numDOF = 0;
    numModes = 0;
    basisValid = false;
}

void
WoodburyLinSOE::resizeX(int size)
{
    if (size <= 0) {
        if (vectX != nullptr) {
            delete vectX;
            vectX = nullptr;
        }
        return;
    }

    if (vectX == nullptr || vectX->Size() != size) {
        if (vectX != nullptr)
            delete vectX;
        vectX = new Vector(size);
    }
}

// Form Z and G for the current A_s (call after tangent assembly; cFactor from integrator).
int
WoodburyLinSOE::formWoodburyBasis(IncrementalIntegrator &integrator,
                                  const Vector *modalFactors)
{
    basisValid = false;
    clearWoodburyBasis();

    if (modalFactors == nullptr)
        return 0;

    if (integrator.setupModal(modalFactors) < 0)
        return -1;

    if (integrator.eigenValues == nullptr || integrator.eigenVectors == nullptr) {
        opserr << "WARNING WoodburyLinSOE::formWoodburyBasis() - modal eigen data not available\n";
        return -2;
    }

    if (innerSOE == nullptr)
        return -1;

    numDOF = innerSOE->getNumEqn();
    resizeX(numDOF);

    const int numEigenModes = modalFactors->Size();
    int maxModes = numEigenModes;
    if (maxModes > integrator.eigenValues->Size())
        maxModes = integrator.eigenValues->Size();

    if (numDOF <= 0 || maxModes <= 0)
        return 0;

    // No modal tangent contribution this assembly (same as addModalDampingMatrix).
    const double cFactor = integrator.getCFactor();
    if (cFactor == 0.0)
        return 0;

    // Q stored column-wise in eigenVectors; Z in woodburyZ (numDOF x numModes).
    eigenVectors = new double[static_cast<size_t>(numDOF) * static_cast<size_t>(maxModes)];
    dBar = new double[maxModes];
    woodburyZ = new Matrix(numDOF, maxModes);

    Vector bSave(innerSOE->getB());
    int activeCol = 0;
    for (int i = 0; i < numEigenModes; ++i) {
        double eigenvalue = (*integrator.eigenValues)(i);
        double xi = (*modalFactors)(i);
        if (eigenvalue <= 0.0 || xi == 0.0)
            continue;

        // dbar_k = cFactor * 2 * xi_i * omega_i
        dBar[activeCol] = cFactor * 2.0 * xi * sqrt(eigenvalue);
        memcpy(&eigenVectors[activeCol * numDOF],
               &integrator.eigenVectors[i * numDOF],
               static_cast<size_t>(numDOF) * sizeof(double));

        Vector qcol(&eigenVectors[activeCol * numDOF], numDOF);
        // Z(:,activeCol) = A_s^{-1} * q_col
        if (innerSOE->setB(qcol) < 0) {
            innerSOE->setB(bSave);
            return -2;
        }
        if (innerSOE->solve() < 0) {
            innerSOE->setB(bSave);
            return -3;
        }
        Vector &xcol = const_cast<Vector &>(innerSOE->getX());
        memcpy(&(*woodburyZ)(0, activeCol), &xcol(0),
               static_cast<size_t>(numDOF) * sizeof(double));

        ++activeCol;
    }
    innerSOE->setB(bSave);

    if (activeCol <= 0) {
        clearWoodburyBasis();
        return 0;
    }

    numModes = activeCol;
    woodburyG = new Matrix(numModes, numModes);
    workV1 = new Vector(numModes);
    workV2 = new Vector(numModes);

    // G = Q^T * Z + diag(1/dbar)
    woodburyG->Zero();
    Matrix Qmat(eigenVectors, numDOF, numModes);
    if (woodburyG->addMatrixTransposeProduct(0.0, Qmat, *woodburyZ, 1.0) < 0)
        return -4;

    for (int i = 0; i < numModes; ++i)
        (*woodburyG)(i, i) += 1.0 / dBar[i];

    hookWoodburySolver();
    basisValid = true;
    return 0;
}

// Given inner x_s = A_s^{-1} b in innerSOE, set vectX = A_eff^{-1} b:
//   r = Q^T * x_s,   y = G^{-1} * r,   x = x_s - Z * y
int
WoodburyLinSOE::applyWoodburyCorrection(void)
{
    if (vectX == nullptr) {
        opserr << "WARNING WoodburyLinSOE::applyWoodburyCorrection() - vectX not allocated\n";
        return -1;
    }

    const Vector &innerX = innerSOE->getX();
    if (innerX.Size() != vectX->Size()) {
        opserr << "WARNING WoodburyLinSOE::applyWoodburyCorrection() - inner/wrapper X size mismatch\n";
        return -1;
    }

    *vectX = innerX;  // x_s

    // No positive eigenmodes: solve A_s x = b only (no Woodbury correction).
    if (!basisValid || numModes <= 0 || eigenVectors == nullptr)
        return 0;

    if (workV1 == nullptr || workV2 == nullptr)
        return -1;

    Matrix Qmat(eigenVectors, numDOF, numModes);
    // workV1 = Q^T * x_s
    if (workV1->addMatrixTransposeVector(0.0, Qmat, *vectX, 1.0) < 0)
        return -1;

    // workV2 = G^{-1} * workV1 = G^{-1} * Q^T * x_s
    if (woodburyG->Solve(*workV1, *workV2) < 0)
        return -1;

    // x = x_s - Z * workV2 = x_s - Z * G^{-1} * Q^T * x_s
    if (vectX->addMatrixVector(1.0, *woodburyZ, *workV2, -1.0) < 0)
        return -1;

    return 0;
}

// A*p += Q * diag(dbar) * Q^T * p  (caller must have set A*p = A_s*p via inner formAp first)
int
WoodburyLinSOE::applyModalMatvec(const Vector &p, Vector &Ap)
{
    if (!basisValid || numModes <= 0)
        return 0;

    if (workV1 == nullptr)
        return -1;

    Matrix Qmat(eigenVectors, numDOF, numModes);
    // workV1 = Q^T * p
    if (workV1->addMatrixTransposeVector(0.0, Qmat, p, 1.0) < 0)
        return -1;

    // workV1 = diag(dbar) * Q^T * p
    for (int i = 0; i < numModes; ++i)
        (*workV1)(i) *= dBar[i];

    // A*p += Q * workV1 = A*p + Q * diag(dbar) * Q^T * p
    if (Ap.addMatrixVector(1.0, Qmat, *workV1, 1.0) < 0)
        return -1;

    return 0;
}

int
WoodburyLinSOE::solve(void)
{
    if (woodburySolver == nullptr)
        hookWoodburySolver();
    if (woodburySolver == nullptr) {
        opserr << "WARNING WoodburyLinSOE::solve() - WoodburySolver not installed\n";
        return -1;
    }
    return woodburySolver->solve();
}

LinearSOESolver *
WoodburyLinSOE::getSolver(void)
{
    return woodburySolver;
}

void
WoodburyLinSOE::hookWoodburySolver(void)
{
    if (woodburySolver != nullptr)
        return;

    LinearSOESolver *innerSolver = innerSOE->getSolver();
    if (innerSolver == nullptr)
        return;

    woodburySolver = new WoodburySolver(*innerSolver, *this);
}

void
WoodburyLinSOE::unhookWoodburySolver(void)
{
    if (woodburySolver == nullptr)
        return;

    delete woodburySolver;
    woodburySolver = nullptr;
}

int
WoodburyLinSOE::addColA(const Vector &col, int colIndex, double fact)
{
    basisValid = false;
    return innerSOE->addColA(col, colIndex, fact);
}

int
WoodburyLinSOE::setSize(Graph &theGraph)
{
    int res = innerSOE->setSize(theGraph);
    clearWoodburyBasis();
    resizeX(innerSOE->getNumEqn());
    return res;
}

int
WoodburyLinSOE::getNumEqn(void) const
{
    return innerSOE->getNumEqn();
}

int
WoodburyLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // A_s changed; Z,G must be rebuilt in formWoodburyBasis (no A_modal in inner A)
    basisValid = false;
    return innerSOE->addA(m, id, fact);
}

int
WoodburyLinSOE::addA(const Matrix &m)
{
    basisValid = false;
    return innerSOE->addA(m);
}

int
WoodburyLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    return innerSOE->addB(v, id, fact);
}

int
WoodburyLinSOE::setB(const Vector &v, double fact)
{
    return innerSOE->setB(v, fact);
}

void
WoodburyLinSOE::zeroA(void)
{
    basisValid = false;
    innerSOE->zeroA();
}

void
WoodburyLinSOE::zeroB(void)
{
    innerSOE->zeroB();
}

// Ap = A_eff * p = A_s*p + Q*diag(dbar)*Q^T*p
int
WoodburyLinSOE::formAp(const Vector &p, Vector &Ap)
{
    int res = innerSOE->formAp(p, Ap);
    if (res < 0)
        return res;
    return applyModalMatvec(p, Ap);
}

const Vector &
WoodburyLinSOE::getX(void)
{
    if (vectX == nullptr) {
        opserr << "FATAL WoodburyLinSOE::getX() - vectX not allocated\n";
        exit(-1);
    }
    return *vectX;
}

const Vector &
WoodburyLinSOE::getB(void)
{
    return innerSOE->getB();
}

const Matrix *
WoodburyLinSOE::getA(void)
{
    if (innerSOE == nullptr)
        return nullptr;
    return innerSOE->getA();
}

double
WoodburyLinSOE::normRHS(void)
{
    return innerSOE->normRHS();
}

void
WoodburyLinSOE::setX(int loc, double value)
{
    if (vectX != nullptr && loc >= 0 && loc < vectX->Size())
        (*vectX)(loc) = value;
}

void
WoodburyLinSOE::setX(const Vector &x)
{
    if (vectX != nullptr && x.Size() == vectX->Size())
        *vectX = x;
}

int
WoodburyLinSOE::setLinks(AnalysisModel &theModel)
{
    LinearSOE::setLinks(theModel);
    return innerSOE->setLinks(theModel);
}

int
WoodburyLinSOE::saveSparseA(OPS_Stream &output, int baseIndex)
{
    if (innerSOE == nullptr)
        return -1;
    return innerSOE->saveSparseA(output, baseIndex);
}

int
WoodburyLinSOE::getSparseA(ID &rowIndices, ID &colIndices, Vector &values,
                           int baseIndex)
{
    if (innerSOE == nullptr)
        return -1;
    return innerSOE->getSparseA(rowIndices, colIndices, values, baseIndex);
}

int
WoodburyLinSOE::getSparseA(std::vector<int> &rowIndices,
                           std::vector<int> &colIndices,
                           std::vector<double> &values, int baseIndex)
{
    if (innerSOE == nullptr)
        return -1;
    return innerSOE->getSparseA(rowIndices, colIndices, values, baseIndex);
}

double
WoodburyLinSOE::getDeterminant(void)
{
    if (innerSOE == nullptr)
        return 0.0;
    return innerSOE->getDeterminant();
}

int
WoodburyLinSOE::sendSelf(int, Channel &)
{
    return 0;
}

int
WoodburyLinSOE::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
wrapLinearSOEWithWoodbury(LinearSOE *&soePtr, Domain *theDomain,
                          DirectIntegrationAnalysis *transient)
{
    if (soePtr == nullptr)
        return -1;

    if (theDomain != nullptr &&
        theDomain->getModalDampingOption() != MODAL_DAMPING_WOODBURY)
        return 0;

    if (transient == nullptr)
        return 0;

    if (dynamic_cast<WoodburyLinSOE *>(soePtr) == nullptr)
        soePtr = new WoodburyLinSOE(*soePtr, true);

    transient->setLinearSOE(*soePtr, false);
    return 0;
}
