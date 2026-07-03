class AnalysisModel;
class LinearSOE;
class ConstraintHandler;
class DOF_Numberer;

class EigenSolve {
public:
  struct Options {
    int    solver     = 0;
    enum ProblemType {
      None,
      Standard,
      Generalized
    } problem         = None;
    int    numMode    = 0;
    double shift      = 0.0;
  };
  EigenSolve(int numMode, bool generalized, bool findSmallest, int typeSolver)
    : numMode(numMode)
    , generalized(generalized)
    , findSmallest(findSmallest)
    , typeSolver(typeSolver) {};

  int solve(AnalysisModel&,
            LinearSOE&,
            ConstraintHandler&,
            const Options&);
private:
  int setup(AnalysisModel&,
            LinearSOE*,
            ConstraintHandler&,
            DOF_Numberer&,
            const Options&);

  int eigen(AnalysisModel& theAnalysisModel,
                    LinearSOE& theSOE,
                    EigenSOE& theEigenSOE,
                    ConstraintHandler& theHandler);
  int clean();
  int numMode;
  bool generalized; 
  bool findSmallest;
  int typeSolver;
  double shift = 0.0;
};