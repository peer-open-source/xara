#pragma once

#include <array>
#include <tuple>
#include <functional>

#include <Element.h>
#include <ID.h>
#include <MatrixND.h>
#include <ShellSection.h>
#include <VectorND.h>

class Domain;
class Matrix;
class Node;
class OPS_Stream;
class Vector;

using namespace OpenSees;


class HeterosisPlate : public Element {
public:
  HeterosisPlate(int tag,
                 const std::array<int, 9>& nodes,
                 ShellSection& section,
                 double drill_scale = 1.0);

  ~HeterosisPlate() override;

  const char* getClassType() const override { return "HeterosisPlate"; }

  int getNumExternalNodes() const override { return nen; }
  const ID& getExternalNodes() override { return node_tags; }
  Node** getNodePtrs() override { return theNodes.data(); }
  int getNumDOF() override { return nxe; }

  void setDomain(Domain*) override;

  int update() override;
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  const Vector& getResistingForce() override;
  const Matrix& getTangentStiff() override;
  const Matrix& getInitialStiff() override;
  const Matrix& getMass() override;

  void Print(OPS_Stream& s, int flag) override;

private:
  // Constant sizes, dimensions, and indices
  constexpr static int nen = 9;         // Number of element nodes
  constexpr static int ndm = 3;         // Number of spatial dimensions
  constexpr static int ndf = 6;         // Number of DOFs per node
  constexpr static int nxe = nen * ndf; // Number of total element DOFs
  constexpr static int nq8 = 8;         //
  constexpr static int nq9 = 9;
  constexpr static int ndof = 26;
  constexpr static int nsr = 8;
  constexpr static int bending_nip = 9;
  constexpr static int shear_nip = 4;
  constexpr static int nip = bending_nip + shear_nip;

  // Ordering of stress resultants in Vectors/Matrices returned by the section
  static constexpr ShellStressLayout scheme = {
      // Membrane
      ShellStress::Fxx,
      ShellStress::Fyy,
      ShellStress::Fxy,
      // Bending
      ShellStress::Mxx,
      ShellStress::Myy,
      ShellStress::Mxy,
      // Shear
      ShellStress::Vxz,
      ShellStress::Vyz,
  };

  using GeometryCoordinates = MatrixND<nq8, 2>;

  static constexpr std::array<std::array<int, 3>, 4> local_edge_nodes = {{
      {0, 4, 1},
      {1, 5, 2},
      {2, 6, 3},
      {3, 7, 0},
  }};

  ID node_tags;                     // array of integer node tags
  std::array<Node*, nen> theNodes;
  std::array<ShellSection*, nip> material;

  VectorND<nxe> p;
  MatrixND<nxe, nxe> K;

  double Ktt;
  double drill_scale;


  std::tuple<int, VectorND<ndof>, MatrixND<ndof, ndof>> compute_stiffness_matrix(
      const GeometryCoordinates& geometry_coordinates,
      const VectorND<ndof>& u_python,
      int tang_flag);

  GeometryCoordinates get_geometry_coordinates() const;
  VectorND<nxe> get_local_displacement_vector() const;
  VectorND<ndof> shell_to_python_displacements(const VectorND<nxe>& u) const;

  static int plate_dof_to_shell_dof(int python_dof);
  static void assemble_python_vector(const VectorND<ndof>& python_vector,
                                     VectorND<nxe>& shell_vector);
  static void assemble_python_matrix(const MatrixND<ndof, ndof>& python_matrix,
                                     MatrixND<nxe, nxe>& shell_matrix);

  static MatrixND<nsr, ndof> assemble_B_matrix(const MatrixND<3, ndof>& B_b);
  static MatrixND<nsr, ndof> assemble_B_matrix(const MatrixND<2, ndof>& B_s);

  void add_drill_stiffness(const VectorND<nxe>& u, int tang_flag);

  VectorND<ndof> compute_edge_force_vector(
      int edge_id,
      const std::function<double(double, double)>& traction) const;

  VectorND<ndof> compute_surface_force_vector(
      const std::function<double(double, double)>& traction) const;

  
  void compute_drill_penalty();

  int form_resid_and_tangent(int tang_flag);
};
