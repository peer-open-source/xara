/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#ifndef WoodburySolver_h
#define WoodburySolver_h

#include <LinearSOESolver.h>

class WoodburyLinSOE;
class LinearSOESolver;

class WoodburySolver : public LinearSOESolver
{
  public:
    WoodburySolver(LinearSOESolver &innerSolver, WoodburyLinSOE &wrap);
    ~WoodburySolver() override;

    int solve(void) override;
    int setSize(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) override;

    LinearSOESolver *releaseInner(void);

  private:
    LinearSOESolver *innerSolver;
    WoodburyLinSOE *theWrapperSOE;
};

#endif
