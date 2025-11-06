

class LinearSOE;

enum class FactorKind { 
    SPD_Cholesky,
    General_QR, 
    Symmetric_LDLT 
};

LinearSOE* CreateAppleSparse(FactorKind kind);
