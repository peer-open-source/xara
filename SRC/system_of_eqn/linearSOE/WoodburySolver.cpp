/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#include <WoodburySolver.h>
#include <WoodburyLinSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

WoodburySolver::WoodburySolver(LinearSOESolver &innerSolver, WoodburyLinSOE &wrap)
    : LinearSOESolver(SOLVER_TAGS_WoodburySolver),
      innerSolver(&innerSolver),
      theWrapperSOE(&wrap)
{
}

WoodburySolver::~WoodburySolver()
{
    innerSolver = nullptr;
}

// A_eff^{-1} b:  x_s = A_s^{-1} b  (inner),  then  x = x_s - Z * G^{-1} * Q^T * x_s
int
WoodburySolver::solve(void)
{
    int r = theWrapperSOE->getInnerSOE().solve();
    if (r < 0)
        return r;

    return theWrapperSOE->applyWoodburyCorrection();
}

int
WoodburySolver::setSize(void)
{
    return innerSolver->setSize();
}

LinearSOESolver *
WoodburySolver::releaseInner(void)
{
    LinearSOESolver *s = innerSolver;
    innerSolver = nullptr;
    return s;
}

int
WoodburySolver::sendSelf(int, Channel &)
{
    return 0;
}

int
WoodburySolver::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

