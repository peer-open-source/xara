//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
#pragma once

#include <array>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>


class Brick02 : public Element {

  public :
    Brick02();
    Brick02(int tag,
            const std::array<int, 8>& node_tags,
            NDMaterial &theMaterial,
            double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);

    // destructor
    virtual ~Brick02();

    const char *getClassType() const final {return "Brick02";}

    FE_Element* createFE_Element(int tag) final;

    void setDomain(Domain *) final;
    int getNumExternalNodes() const final;
    const ID &getExternalNodes() final;
    Node **getNodePtrs() final;
    int getNumDOF() final;

    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // return stiffness matrix
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    void zeroLoad();
    int addLoad(ElementalLoad *, double loadFactor);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();

    // public methods for element output

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &);


    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int parameterID, Information &);

    void Print(OPS_Stream &s, int flag);

private:
    //
    // private methods
    //
    void formInertiaTerms(int tangFlag, Vector &resid);
    void computeBasis();

    const MatrixND<6,3>&
    computeB(int node,
             const double shp[4][8],
             MatrixND<6,3> &B
    ) const noexcept;


    //
    // private attributes
    //
    constexpr static unsigned int
                         NEN = 8,  // number of element nodes
                         NDM = 3,  // Spatial dimensions
                         NDF = 3,  // number of element dof
                         NIP = 8,  // number of integration points
                         NST = 6,  // number of stress components
                         NDOF = NEN*NDF;

    ID conn;       // node tags
    std::array<Node *, NEN> theNodes; // pointers to nodes

    // material information
    NDMaterial *materialPointers[NIP]; // pointers to materials

    double b[3];          // Body forces
    double appliedB[3];   // Body forces applied with load
    int applyLoad;

    Vector *load;
    Matrix *Ki;

    MatrixND<NDOF,NDOF> stiff;
    VectorND<NDOF> resid;
    VectorND<NDOF> response;
    VectorND<NDOF> inertia;

    Matrix K_wrap;
    Vector p_wrap;
    Vector response_wrap;
    Vector inertia_wrap;

    //
    // static attributes
    //

    static Matrix mass;

    // quadrature data
    static constexpr double root3 = 1.73205080757;
    static constexpr double sg[2] = {
      -1.0/root3, 1.0/root3
    };
    static constexpr double wg[NIP] = {
                              1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0 };

    // local nodal coordinates, three coordinates for each node
    double xl[NDM][NEN];

};
