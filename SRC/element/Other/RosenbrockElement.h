//
// An element implementing the following nonlinear problem from [1]:
//
// g(x)=\frac{\partial \Phi}{\partial x}=x^3+3 x^2+3 x-2.241=0 \text { or }(x+1)^3=3 \\
// g(y)=\frac{\partial \Phi}{\partial y}=y^3-3 y^2+3 y-3=0 \text { or }(y-1)^3=2
// 
// References:
//
// [1] M. A. Crisfield, “Accelerating and damping the modified Newton-Raphson method,” 
//     Computers & Structures, vol. 18, no. 3, pp. 395–407, Jan. 1984, 
//     doi: 10.1016/0045-7949(84)90059-2.
//
#include <Element.h>
#include <vector>
#include <array>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Node.h>
#include <Domain.h>

namespace OpenSees {

class RosenbrockElement : public Element
{
public:
  RosenbrockElement(int tag, const std::array<int,1>& nodeTags) :
    Element(tag, 0),
    nodeTags(1),
    node(nullptr),
    K{}, R{},
    K_wrap(K),
    R_wrap(R)
  {
    this->nodeTags(0) = nodeTags[0];
  }

  virtual ~RosenbrockElement() {}

  void setDomain(Domain *theDomain) override {
    node = theDomain->getNode(nodeTags(0));
    if (node == nullptr) {
      opserr << "No node with tag " << nodeTags(0) << " exists in the domain\n";
      return;
    }
  }

  int getNumExternalNodes() const override { return 1; }

  const ID &getExternalNodes() override { return nodeTags; }

  Node **getNodePtrs() override { return &node; }

  int getNumDOF() override { return 2; }

  int commitState() override { return 0; }
  int revertToLastCommit() override { return 0; }
  int update() override { return 0; }

  const Vector& getResistingForce() final {
    const Vector& u = node->getTrialDisp();
    double x = u(0);
    double y = u(1);
    opserr << "x = " << u;
    R[0] = x*x*x + 3.0*x*x + 3.0*x;// - 2.241;
    R[1] = y*y*y - 3.0*y*y + 3.0*y;// - 3;
    opserr << "r = " << R_wrap;
    return R_wrap;
  }

  const Matrix& getTangentStiff() final {
    const Vector& u = node->getTrialDisp();
    double x = u(0);
    double y = u(1);
    tangent(x, y, K);
    return K_wrap;
  }

  const Matrix& getInitialStiff() final {
    tangent(0.0, 0.0, K);
    return K_wrap;
  }

  void Print(OPS_Stream &s, int flag) final {
    s << "RosenbrockElement, tag: " << this->getTag() << "\n";
    s << "Node: " << nodeTags(0) << "\n";
  }


private:
  void tangent(double x, double y, MatrixND<2, 2>& K) {
      K(0,0) = 3.0*x*x + 6.0*x + 3.0;
      K(0,1) = 0;
      K(1,0) = 0;
      K(1,1) = 3.0*y*y - 6.0*y + 3.0;
  }

  Node* node;
  ID nodeTags;
  MatrixND<2, 2> K; // Stiffness matrix
  VectorND<2> R;    // Residual vector
  Vector R_wrap;
  Matrix K_wrap;
};

} // namespace OpenSees
