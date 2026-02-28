class AnalysisModel;
class LinearSOE;
class ConstraintHandler;

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

  int solve(AnalysisModel&,
            LinearSOE&,
            ConstraintHandler&,
            const Options&);
private:
  int setup(AnalysisModel&,
            LinearSOE&,
            const Options&);
  int clean();
};