/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of a linear coordinate transformation 4-node shells
//
#ifndef ASDShellQ4Transformation_h
#define ASDShellQ4Transformation_h

#include <Vector3D.h>
#include <ASDShellQ4LocalCoordinateSystem.h>
#include <Node.h>
#include <ID.h>
#include <Domain.h>
#include <Logging.h>

/** \brief ASDShellQ4Transformation
*
* This class represents a basic (linear) coordinate transformation that can be used
* by any element whose geometry is a QUAD 4 in 3D space, with 6 D.O.F.s per node.
* Its main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
* 3) Transform the outgoing matrices and vectors in global coordinate system
*/
class ASDShellQ4Transformation
{

public:

    typedef Vector3D Vector3Type;
    typedef ASDQuaternion<double> QuaternionType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef std::array<Node*, 4> NodeContainerType;

public:

    ASDShellQ4Transformation()
    {
    }

    virtual ~ASDShellQ4Transformation()
    {
    }

private:

    ASDShellQ4Transformation(const ASDShellQ4Transformation& other) = delete;

    ASDShellQ4Transformation& operator = (const ASDShellQ4Transformation& other) = delete;

public:

    // virtual ASDShellQ4Transformation* create()const
    // {
    //     return new ASDShellQ4Transformation();
    // }

    virtual bool 
    isLinear() const
    {
        return true;
    }

    virtual void 
    revertToStart()
    {
    }

    virtual void 
    setDomain(Domain* domain, const ID& node_ids, bool initialized)
    {
        // if domain is null
        if (domain == nullptr) {
            for (size_t i = 0; i < 4; i++) {
                m_nodes[i] = nullptr;
            }
            return;
        }
        // get nodes and save initial displacements and rotations
        for (size_t i = 0; i < 4; i++) {
            m_nodes[i] = domain->getNode(node_ids(i));
            if (m_nodes[i] == nullptr) {
                opserr << "ASDShellQ4Transformation::setDomain - no node " << node_ids(i)
                    << " exists in the model\n";
                exit(-1);
            }
            if (!initialized) {
              const Vector& iU = m_nodes[i]->getTrialDisp();
              if (iU.Size() != 6) {
                  opserr << "ASDShellQ4Transformation::setDomain - node " << node_ids(i)
                      << " has " << iU.Size() << " DOFs, while 6 are expected\n";
                  exit(-1);
              }
              size_t index = i * 6;
              for (size_t j = 0; j < 6; j++)
                  m_U0(index + j) = iU(j);
            }
        }
    }

    virtual void revertToLastCommit()
    {
    }

    virtual void commit()
    {
    }

    virtual void update(const VectorType& globalDisplacements)
    {
    }

    virtual ASDShellQ4LocalCoordinateSystem 
    createReferenceCoordinateSystem() const
    {
        // the reference coordinate system in the underformed configuration
        // using the default alignment to the first column of the jacobian at center
        return ASDShellQ4LocalCoordinateSystem(
            m_nodes[0]->getCrds(),
            m_nodes[1]->getCrds(),
            m_nodes[2]->getCrds(),
            m_nodes[3]->getCrds()
        );
    }

    virtual ASDShellQ4LocalCoordinateSystem 
    createLocalCoordinateSystem(const Vector& globalDisplacements) const
    {
        // same as reference
        return createReferenceCoordinateSystem();
    }

    void computeGlobalDisplacements(Vector& globalDisplacements) const
    {
        for (int i = 0; i < 4; i++) {
            int index = i * 6;
            const Vector& iU = m_nodes[i]->getTrialDisp();
            for (int j = 0; j < 6; j++) {
                globalDisplacements(index + j) = iU(j) - m_U0(index + j);
            }
        }
    }

    virtual const Matrix& computeTransformationMatrix(const ASDShellQ4LocalCoordinateSystem& LCS) const
    {
        static Matrix R(24, 24);
        static Matrix T(24, 24);
        static Matrix W(24, 24);
        if (LCS.IsWarped()) {
            R.Zero();
            LCS.ComputeTotalRotationMatrix(R);
            LCS.ComputeTotalWarpageMatrix(W);
            T.addMatrixProduct(0.0, W, R, 1.0);
        }
        else {
            T.Zero();
            LCS.ComputeTotalRotationMatrix(T);
        }
        return T;
    }

    virtual void calculateLocalDisplacements(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        const Matrix& R = computeTransformationMatrix(LCS);
        localDisplacements.addMatrixVector(0.0, R, globalDisplacements, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        Matrix& LHS,
        Vector& RHS,
        bool LHSrequired)
    {
        static MatrixType RT_LHS(24, 24);
        static VectorType RHScopy(24);
        const MatrixType& R = computeTransformationMatrix(LCS);
        RHScopy = RHS;
        RHS.addMatrixTransposeVector(0.0, R, RHScopy, 1.0);
        if (LHSrequired) {
            RT_LHS.addMatrixTransposeProduct(0.0, R, LHS, 1.0);
            LHS.addMatrixProduct(0.0, RT_LHS, R, 1.0);
        }
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType dummy;
        transformToGlobal(LCS, dummy, dummy, LHS, RHS, LHSrequired);
    }

    virtual int internalDataSize() const
    {
        // just the size of the initial displacements
        return 24;
    }

    virtual void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellQ4Transformation - failed to save internal data: vector too small\n";
            exit(-1);
        }
        for (int i = 0; i < 24; i++)
            v(pos++) = m_U0(i);
    }

    virtual void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellQ4Transformation - failed to restore internal data: vector too small\n";
            exit(-1);
        }
        for (int i = 0; i < 24; i++)
            m_U0(i) = v(pos++);
    }

public:

    inline const NodeContainerType& getNodes() const { return m_nodes; }
    inline NodeContainerType& getNodes() { return m_nodes; }

protected:

    NodeContainerType m_nodes = { {nullptr, nullptr, nullptr, nullptr} };
    Vector m_U0 = Vector(24);
};

#endif // !ASDShellQ4Transformation_h

