#include "HeterosisPlate.h"

#include <cmath>
#include <functional>

#include <Logging.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Domain.h>
#include <Node.h>
#include <eigen/SymEigDirect3D.h>


using namespace OpenSees;


// TODO: move to separate header
namespace Plate::Quadrature {

template <int n>
struct GaussRule1D {
  VectorND<n> points;
  VectorND<n> weights;
};

template <int order>
GaussRule1D<order>
gauss_legendre_1d();

template <>
GaussRule1D<2>
gauss_legendre_1d<2>()
{
  double a = 1.0 / std::sqrt(3.0);
  return {{-a, a}, {1.0, 1.0}};
}

template <>
GaussRule1D<3>
gauss_legendre_1d<3>()
{
  double a = std::sqrt(3.0 / 5.0);
  return {{-a, 0.0, a}, {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}};
}

template <int order_x, int order_y>
struct TensorProductRule {
  static constexpr int n = order_x * order_y;

  std::array<VectorND<2>, n> points;
  VectorND<n> weights;
};

template <int order_x, int order_y>
TensorProductRule<order_x, order_y>
tensor_product_rule()
{
  TensorProductRule<order_x, order_y> rule;
  const auto rule_x = gauss_legendre_1d<order_x>();
  const auto rule_y = gauss_legendre_1d<order_y>();

  int point_id = 0;
  for (int i = 0; i < order_x; i++) {
    for (int j = 0; j < order_y; j++) {
      rule.points[point_id] = {rule_x.points[i], rule_y.points[j]};
      rule.weights[point_id] = rule_x.weights[i] * rule_y.weights[j];
      point_id++;
    }
  }

  return rule;
}
}


namespace Plate::Interpolation {
VectorND<8>
q8_shape_functions(double xi, double eta)
{
  return {
    0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0),
    0.25 * (1.0 + xi) * (1.0 - eta) * ( xi - eta - 1.0),
    0.25 * (1.0 + xi) * (1.0 + eta) * ( xi + eta - 1.0),
    0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0),
    0.50 * (1.0 - xi * xi) * (1.0 - eta),
    0.50 * (1.0 + xi) * (1.0 - eta * eta),
    0.50 * (1.0 - xi * xi) * (1.0 + eta),
    0.50 * (1.0 - xi) * (1.0 - eta * eta),
  };
}

std::tuple<VectorND<8>, VectorND<8>>
q8_shape_function_gradients_parent(double xi, double eta)
{
  const VectorND<8> dN_dxi = {
      -0.25 * (eta - 1.0) * (eta + 2.0 * xi),
       0.25 * (eta - 1.0) * (eta - 2.0 * xi),
       0.25 * (eta + 1.0) * (eta + 2.0 * xi),
      -0.25 * (eta + 1.0) * (eta - 2.0 * xi),
       xi * (eta - 1.0),
      -0.50 * (eta - 1.0) * (eta + 1.0),
      -xi * (eta + 1.0),
       0.50 * (eta - 1.0) * (eta + 1.0),
  };

  const VectorND<8> dN_deta = {
      -0.25 * (2.0 * eta + xi) * (xi - 1.0),
       0.25 * (2.0 * eta - xi) * (xi + 1.0),
       0.25 * (2.0 * eta + xi) * (xi + 1.0),
      -0.25 * (2.0 * eta - xi) * (xi - 1.0),
       0.50 * (xi - 1.0) * (xi + 1.0),
      -eta * (xi + 1.0),
      -0.50 * (xi - 1.0) * (xi + 1.0),
       eta * (xi - 1.0),
  };

  return std::tuple{dN_dxi, dN_deta};
}

VectorND<9>
q9_shape_functions(double xi, double eta)
{
  const VectorND<3> L_xi = {
    0.50 * xi * (xi - 1.0), 1.0 - xi * xi, 0.50 * xi * (xi + 1.0),
  };

  const VectorND<3> L_eta = {
    0.50 * eta * (eta - 1.0),  1.0 - eta * eta, 0.50 * eta * (eta + 1.0),
  };

  return {
      L_xi[0] * L_eta[0],
      L_xi[2] * L_eta[0],
      L_xi[2] * L_eta[2],
      L_xi[0] * L_eta[2],
      L_xi[1] * L_eta[0],
      L_xi[2] * L_eta[1],
      L_xi[1] * L_eta[2],
      L_xi[0] * L_eta[1],
      L_xi[1] * L_eta[1],
  };
}

std::tuple<VectorND<9>, VectorND<9>>
q9_shape_function_gradients_parent(double xi, double eta)
{
  const VectorND<3> L_xi = {
    0.50 * xi * (xi - 1.0),1.0 - xi * xi, 0.50 * xi * (xi + 1.0)
  };

  const VectorND<3> dL_xi = {
    xi - 0.50, -2.0 * xi, xi + 0.50,
  };

  const VectorND<3> L_eta = {
    0.50 * eta * (eta - 1.0), 1.0 - eta * eta, 0.50 * eta * (eta + 1.0)
  };

  const VectorND<3> dL_eta = {
    eta - 0.50,  -2.0 * eta,  eta + 0.50
  };

  const VectorND<9> dN_dxi = {
      dL_xi[0] * L_eta[0],
      dL_xi[2] * L_eta[0],
      dL_xi[2] * L_eta[2],
      dL_xi[0] * L_eta[2],
      dL_xi[1] * L_eta[0],
      dL_xi[2] * L_eta[1],
      dL_xi[1] * L_eta[2],
      dL_xi[0] * L_eta[1],
      dL_xi[1] * L_eta[1],
  };

  const VectorND<9> dN_deta = {
      L_xi[0] * dL_eta[0],
      L_xi[2] * dL_eta[0],
      L_xi[2] * dL_eta[2],
      L_xi[0] * dL_eta[2],
      L_xi[1] * dL_eta[0],
      L_xi[2] * dL_eta[1],
      L_xi[1] * dL_eta[2],
      L_xi[0] * dL_eta[1],
      L_xi[1] * dL_eta[1],
  };

  return std::tuple{dN_dxi, dN_deta};
}

VectorND<3>
edge_quadratic_shape_functions(double s)
{
  return {
      0.50 * s * (s - 1.0),
      1.0 - s * s,
      0.50 * s * (s + 1.0),
  };
}

VectorND<3>
edge_quadratic_shape_function_derivatives(double s)
{
  return {
    s - 0.50, -2.0 * s, s + 0.50,
  };
}

MatrixND<2, 2>
geometry_jacobian(double xi, double eta,
                  const MatrixND<8, 2>& geometry_coordinates)
{
  const auto [dN_dxi, dN_deta] = q8_shape_function_gradients_parent(xi, eta);

  MatrixND<2, 2> jacobian;
  jacobian.zero();

  for (int i = 0; i < 8; i++) {
    jacobian(0, 0) += dN_dxi[i]  * geometry_coordinates(i, 0);
    jacobian(0, 1) += dN_deta[i] * geometry_coordinates(i, 0);
    jacobian(1, 0) += dN_dxi[i]  * geometry_coordinates(i, 1);
    jacobian(1, 1) += dN_deta[i] * geometry_coordinates(i, 1);
  }

  return jacobian;
}

double
positive_area_jacobian_det(const MatrixND<2, 2>& jacobian,
                           int element_id,
                           const char* context)
{
  double det = jacobian.determinant();

  if (det <= 0.0) {
    opserr << "Element "
           << element_id 
           << " has non-positive Jacobian in " 
           << context << "\n";
  }

  return det;
}


template <int n>
std::tuple<VectorND<n>, VectorND<n>>
parent_to_physical_gradients(const VectorND<n>& dN_dxi,
                             const VectorND<n>& dN_deta,
                             const MatrixND<2, 2>& jacobian)
{
  VectorND<n> dN_dx;
  VectorND<n> dN_dy;

  double det = jacobian.determinant();

  // invert jacobian and store into inv_jacobian.
  MatrixND<2,2> inv_jacobian;
  jacobian.invert(inv_jacobian);

  for (int i = 0; i < n; i++) {
    dN_dx[i] = ( jacobian(1, 1) * dN_dxi[i]
               - jacobian(1, 0) * dN_deta[i]) / det;
    dN_dy[i] = (-jacobian(0, 1) * dN_dxi[i]
               + jacobian(0, 0) * dN_deta[i]) / det;
  }

  return std::tuple{dN_dx, dN_dy};
}

MatrixND<3, 26>
bending_B_matrix(const VectorND<9>& dN_theta_dx,
                 const VectorND<9>& dN_theta_dy)
{
  MatrixND<3, 26> B_b{};

  for (int local_node_id = 0; local_node_id < 9; local_node_id++) {
    const int theta_x_column = 8 + 2 * local_node_id;
    const int theta_y_column = theta_x_column + 1;

    B_b(0, theta_x_column) = dN_theta_dx[local_node_id];
    B_b(1, theta_y_column) = dN_theta_dy[local_node_id];
    B_b(2, theta_x_column) = dN_theta_dy[local_node_id];
    B_b(2, theta_y_column) = dN_theta_dx[local_node_id];
  }

  return B_b;
}

MatrixND<2, 26>
shear_B_matrix(const VectorND<8>& dN_w_dx,
               const VectorND<8>& dN_w_dy,
               const VectorND<9>& N_theta)
{
  MatrixND<2, 26> B_s{};

  for (int local_node_id = 0; local_node_id < 8; local_node_id++) {
    B_s(0, local_node_id) = dN_w_dx[local_node_id];
    B_s(1, local_node_id) = dN_w_dy[local_node_id];
  }

  for (int local_node_id = 0; local_node_id < 9; local_node_id++) {
    const int theta_x_column = 8 + 2 * local_node_id;
    const int theta_y_column = theta_x_column + 1;

    B_s(0, theta_x_column) = -N_theta[local_node_id];
    B_s(1, theta_y_column) = -N_theta[local_node_id];
  }

  return B_s;
}

} // namespace Plate::Interpolation


using namespace Plate::Quadrature;
using namespace Plate::Interpolation;



//
// This is like python's __init__ method.
// it is invoked when an instance is created.
//
HeterosisPlate::HeterosisPlate(int tag,
                               const std::array<int, 9>& nodes,
                               ShellSection& section,
                               double drill_scale)
  : Element(tag, 0),
    node_tags(nen),
    Ktt(0.0),
    drill_scale(drill_scale)
{
  theNodes.fill(nullptr);
  material.fill(nullptr);

  p.zero();
  K.zero();

  // store node ids; the actual node objects are obtained 
  // when setDomain() is called (below).
  for (int i = 0; i < nen; i++)
    node_tags(i) = nodes[i];

  // Each integration point gets its own copy of the section.
  // This is necessary for nonlinear material response
  for (int i = 0; i < nip; i++) {
    material[i] = section.getCopy();
  }
}


//: similar to python's __del__ method; 
HeterosisPlate::~HeterosisPlate()
{
  for (auto& section : material) {
    delete section;
    section = nullptr;
  }

  theNodes.fill(nullptr);
}



// This is invoked immediately after an element is created
void
HeterosisPlate::setDomain(Domain* domain)
{
  if (domain == nullptr)
    return;

  for (int i = 0; i < nen; i++) {
    theNodes[i] = domain->getNode(node_tags(i));
    assert(theNodes[i] != nullptr);
  }

  compute_drill_penalty();
  this->Element::link(*domain);
}

int
HeterosisPlate::commitState()
{
  int success = this->Element::commitState();

  for (auto* section : material) {
    if (section != nullptr)
      success += section->commitState();
  }

  return success;
}

int
HeterosisPlate::revertToLastCommit()
{
  int success = 0;

  for (auto* section : material) {
    if (section != nullptr)
      success += section->revertToLastCommit();
  }

  return success;
}

int
HeterosisPlate::revertToStart()
{
  int success = 0;

  for (auto* section : material) {
    if (section != nullptr)
      success += section->revertToStart();
  }

  p.zero();
  K.zero();

  return success;
}


int
HeterosisPlate::update()
{
  return 0; //form_resid_and_tangent(0);
}

HeterosisPlate::GeometryCoordinates
HeterosisPlate::get_geometry_coordinates() const
{
  GeometryCoordinates geometry_coordinates;
  geometry_coordinates.zero();

  for (int i = 0; i < nq8; i++) {
    const Vector& coordinates = theNodes[i]->getCrds();
    geometry_coordinates(i, 0) = coordinates(0);
    geometry_coordinates(i, 1) = coordinates(1);
  }

  return geometry_coordinates;
}

VectorND<HeterosisPlate::nxe>
HeterosisPlate::get_local_displacement_vector() const
{
  VectorND<nxe> u;
  u.zero();

  for (int node = 0; node < nen; node++) {
    const Vector& trial_disp = theNodes[node]->getTrialDisp();
    for (int dof = 0; dof < ndf; dof++)
      u[node * ndf + dof] = trial_disp(dof);
  }

  return u;
}

VectorND<HeterosisPlate::ndof>
HeterosisPlate::shell_to_python_displacements(const VectorND<nxe>& u) const
{
  VectorND<ndof> u_python;
  u_python.zero();

  for (int i = 0; i < ndof; i++)
    u_python[i] = u[plate_dof_to_shell_dof(i)];

  return u_python;
}

int
HeterosisPlate::plate_dof_to_shell_dof(int python_dof)
{
  if (python_dof < 8)
    return python_dof * ndf + 2;

  int theta_dof = python_dof - 8;
  int local_node_id = theta_dof / 2;
  int component = theta_dof % 2;

  return local_node_id * ndf + 3 + component;
}


void
HeterosisPlate::assemble_python_vector(const VectorND<ndof>& python_vector,
                                       VectorND<nxe>& shell_vector)
{
  for (int i = 0; i < ndof; i++)
    shell_vector[plate_dof_to_shell_dof(i)] += python_vector[i];
}



std::tuple<int, VectorND<HeterosisPlate::ndof>, MatrixND<HeterosisPlate::ndof, HeterosisPlate::ndof>>
HeterosisPlate::compute_stiffness_matrix(const GeometryCoordinates& geometry_coordinates,
                                         const VectorND<ndof>& u_python,
                                         int tang_flag)
{
  // the braces {} here ensure that the vectors and matrices are zero-initialized.
  VectorND<ndof> p_local{};
  VectorND<ndof> p_b{};
  VectorND<ndof> p_s{};
  MatrixND<ndof, ndof> K_local{};
  MatrixND<ndof, ndof> K_b{};
  MatrixND<ndof, ndof> K_s{};

  //
  // Bending
  //
  const auto bending_rule = tensor_product_rule<3, 3>();
  for (int q = 0; q < bending_nip; q++) {
    double xi = bending_rule.points[q][0];
    double eta = bending_rule.points[q][1];
    double weight = bending_rule.weights[q];

    const auto jacobian = geometry_jacobian(xi, eta, geometry_coordinates);
    double det_jacobian = positive_area_jacobian_det(
        jacobian, this->getTag(), "bending quadrature"
    );
    double dvol = weight * det_jacobian;

    const auto [dN_theta_dxi, dN_theta_deta] = q9_shape_function_gradients_parent(xi, eta);
    const auto [dN_theta_dx, dN_theta_dy] = parent_to_physical_gradients(
        dN_theta_dxi,
        dN_theta_deta,
        jacobian
    );

    const auto B_b = bending_B_matrix(dN_theta_dx, dN_theta_dy);

    const auto B_strain = assemble_B_matrix(B_b);
    auto B_equil = B_strain;

    // ShellSection bending equilibrium sign.
    for (int row = 3; row < 6; row++)
      for (int col = 0; col < ndof; col++)
        B_equil(row, col) *= -1.0;

    const VectorND<nsr> strain = B_strain * u_python;

    if (material[q]->setTrialState<nsr, scheme>(strain) != 0)
      return std::tuple{-1, p_local, K_local};

    const VectorND<nsr> stress = material[q]->getResultant<nsr, scheme>();
    const VectorND<ndof> local_force = B_equil ^ stress;

    for (int i = 0; i < ndof; i++)
      p_b[i] += dvol * local_force[i];

    if (tang_flag == 1) {
      const MatrixND<nsr, nsr> dd = material[q]->getTangent<nsr, scheme>();
      K_b.addMatrixTripleProduct(1.0, B_equil, dd, B_strain, dvol);
    }
  }

  //
  // Shear
  //
  const auto shear_rule = tensor_product_rule<2, 2>();
  for (int q = 0; q < shear_nip; q++) {
    double xi     = shear_rule.points[q][0];
    double eta    = shear_rule.points[q][1];
    double weight = shear_rule.weights[q];

    auto* section = material[bending_nip + q];

    const auto jacobian = geometry_jacobian(xi, eta, geometry_coordinates);
    double det_jacobian = positive_area_jacobian_det(
      jacobian, this->getTag(), "shear quadrature"
    );
    double dvol = weight * det_jacobian;

    const auto [dN_w_dxi, dN_w_deta] =
        q8_shape_function_gradients_parent(xi, eta);
    const auto [dN_w_dx, dN_w_dy] = parent_to_physical_gradients(
        dN_w_dxi,
        dN_w_deta,
        jacobian
    );
    const auto N_theta = q9_shape_functions(xi, eta);

    const auto B_s = shear_B_matrix(dN_w_dx, dN_w_dy, N_theta);
    const auto B   = assemble_B_matrix(B_s);
    VectorND<nsr> strain = B * u_python;

    if (section->setTrialState<nsr, scheme>(strain) != 0)
      return std::tuple{-1, p_local, K_local};

    VectorND<nsr> stress = section->getResultant<nsr, scheme>();
    VectorND<ndof> local_force = B^stress;

    for (int i = 0; i < ndof; i++)
      p_s[i] += dvol * local_force[i];

    if (tang_flag == 1) {
      MatrixND<nsr, nsr> dd = section->getTangent<nsr, scheme>();
      K_s.addMatrixTripleProduct(1.0, B, dd, B, dvol);
    }
  }

  for (int i = 0; i < ndof; i++)
    p_local[i] = p_b[i] + p_s[i];

  if (tang_flag == 1) {
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < ndof; j++)
        K_local(i, j) = K_b(i, j) + K_s(i, j);
  }

  return std::tuple{0, p_local, K_local};
}


void
HeterosisPlate::assemble_python_matrix(const MatrixND<ndof, ndof>& python_matrix,
                                       MatrixND<nxe, nxe>& shell_matrix)
{
  for (int i = 0; i < ndof; i++) {
    const int ii = plate_dof_to_shell_dof(i);
    for (int j = 0; j < ndof; j++) {
      const int jj = plate_dof_to_shell_dof(j);
      shell_matrix(ii, jj) += python_matrix(i, j);
    }
  }
}

MatrixND<HeterosisPlate::nsr, HeterosisPlate::ndof>
HeterosisPlate::assemble_B_matrix(const MatrixND<3, ndof>& B_b)
{
  MatrixND<nsr, ndof> B;
  B.zero();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < ndof; j++)
      B(i + 3, j) = B_b(i, j);

  return B;
}

MatrixND<HeterosisPlate::nsr, HeterosisPlate::ndof>
HeterosisPlate::assemble_B_matrix(const MatrixND<2, ndof>& B_s)
{
  MatrixND<nsr, ndof> B;
  B.zero();

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < ndof; j++)
      B(i + 6, j) = B_s(i, j);

  return B;
}

int
HeterosisPlate::form_resid_and_tangent(int tang_flag)
{
  const auto geometry_coordinates = get_geometry_coordinates();
  const auto u = get_local_displacement_vector();
  const auto u_python = shell_to_python_displacements(u);

  p.zero();

  if (tang_flag == 1)
    K.zero();

  const auto [success, p_python, K_python] = compute_stiffness_matrix(
      geometry_coordinates,
      u_python,
      tang_flag
  );

  if (success != 0)
    return success;

  assemble_python_vector(p_python, p);

  if (tang_flag == 1)
    assemble_python_matrix(K_python, K);

  add_drill_stiffness(u, tang_flag);

  return 0;
}


void
HeterosisPlate::compute_drill_penalty()
{
  if (material[0] == nullptr) {
    Ktt = 0.0;
    return;
  }

  VectorND<nsr> zero_strain{};

  auto* section = material[0];
  section->setTrialState<nsr, scheme>(zero_strain);

  static constexpr ShellStressLayout m_layout {
    ShellStress::Fxx,
    ShellStress::Fyy,
    ShellStress::Fxy,
  };

  MatrixND<3, 3> Dm = section->getTangent<3, m_layout>();

  SymEigDirect3D<double,-1> eigensolver;
  std::array<VectorND<3>, 3> eigenvectors;
  VectorND<3> eigenvalues;
  eigensolver(Dm(0,0), Dm(0,1), Dm(0,2), 
              Dm(1,1), Dm(1,2), Dm(2,2), 
              eigenvalues, eigenvectors);


  Ktt = drill_scale * eigenvalues[2];
}


void
HeterosisPlate::add_drill_stiffness(const VectorND<nxe>& u, int tang_flag)
{
  if (Ktt == 0.0)
    return;

  for (int node = 0; node < nen; node++) {
    const int u_dof = node * ndf + 0;
    const int v_dof = node * ndf + 1;
    const int rz_dof = node * ndf + 5;

    p[u_dof] += Ktt * u[u_dof];
    p[v_dof] += Ktt * u[v_dof];
    p[rz_dof] += Ktt * u[rz_dof];

    if (tang_flag == 1) {
      K(u_dof, u_dof) += Ktt;
      K(v_dof, v_dof) += Ktt;
      K(rz_dof, rz_dof) += Ktt;
    }
  }

  const int center_w = 8 * ndf + 2;
  p[center_w] += Ktt * u[center_w];

  if (tang_flag == 1)
    K(center_w, center_w) += Ktt;
}



// Integrate distributed transverse traction on one element edge into local DOFs.
//
// The load contributes only to w DOFs located on the selected quadratic edge
// (3 edge nodes in local Q8 numbering).
VectorND<HeterosisPlate::ndof>
HeterosisPlate::compute_edge_force_vector(
    int edge_id,
    const std::function<double(double, double)>& traction) const
{
  VectorND<ndof> local_force;
  local_force.zero();

  assert(edge_id >= 1 && edge_id <= 4);

  const auto geometry_coordinates = get_geometry_coordinates();
  const auto& edge_node_ids_local = local_edge_nodes[edge_id - 1];
  const auto edge_rule = gauss_legendre_1d<3>();

  for (int q = 0; q < 3; q++) {
    double s = edge_rule.points[q];
    double weight = edge_rule.weights[q];
    const auto N_edge = edge_quadratic_shape_functions(s);
    const auto dN_edge_ds = edge_quadratic_shape_function_derivatives(s);

    double x_q = 0.0;
    double y_q = 0.0;
    double dx_ds = 0.0;
    double dy_ds = 0.0;

    for (int i = 0; i < 3; i++) {
      const int node = edge_node_ids_local[i];
      x_q   += N_edge[i] * geometry_coordinates(node, 0);
      y_q   += N_edge[i] * geometry_coordinates(node, 1);
      dx_ds += dN_edge_ds[i] * geometry_coordinates(node, 0);
      dy_ds += dN_edge_ds[i] * geometry_coordinates(node, 1);
    }

    double jacobian_edge = std::sqrt(dx_ds * dx_ds + dy_ds * dy_ds);
    double q_q = traction(x_q, y_q);

    for (int i = 0; i < 3; i++) {
      const int node = edge_node_ids_local[i];
      local_force[node] += weight * N_edge[i] * q_q * jacobian_edge;
    }
  }

  return local_force;
}

VectorND<HeterosisPlate::ndof>
HeterosisPlate::compute_surface_force_vector(
    const std::function<double(double, double)>& traction) const
{
  VectorND<ndof> local_force;
  local_force.zero();

  const auto geometry_coordinates = get_geometry_coordinates();
  const auto area_rule = tensor_product_rule<3, 3>();

  // Gauss quadrature loop to integrate local_force
  for (int q = 0; q < 9; q++) {
    double xi = area_rule.points[q][0];
    double eta = area_rule.points[q][1];
    double weight = area_rule.weights[q];

    const auto N_w = q8_shape_functions(xi, eta);
    const auto jacobian = geometry_jacobian(xi, eta, geometry_coordinates);
    double det_jacobian = positive_area_jacobian_det(
      jacobian, this->getTag(), "surface pressure quadrature"
    );

    double x_q = 0.0;
    double y_q = 0.0;

    for (int i = 0; i < nq8; i++) {
      x_q += N_w[i] * geometry_coordinates(i, 0);
      y_q += N_w[i] * geometry_coordinates(i, 1);
    }

    double q_q = traction(x_q, y_q);

    for (int i = 0; i < nq8; i++)
      local_force[i] += weight * N_w[i] * q_q * det_jacobian;
  }

  return local_force;
}


//
// Interface methods to OpenSees; 
// these call the internal methods above
//
const Vector&
HeterosisPlate::getResistingForce()
{
  form_resid_and_tangent(0);

  thread_local Vector wrapper(0);
  wrapper.setData(p);
  return wrapper;
}

const Matrix&
HeterosisPlate::getTangentStiff()
{
  form_resid_and_tangent(1);

  thread_local Matrix wrapper(0, 0);
  wrapper.setData(K);
  return wrapper;
}

const Matrix&
HeterosisPlate::getInitialStiff()
{
  return getTangentStiff();
}

const Matrix&
HeterosisPlate::getMass()
{
  thread_local MatrixND<nxe, nxe> M;
  thread_local Matrix wrapper(0, 0);

  M.zero();
  wrapper.setData(M);
  return wrapper;
}

void
HeterosisPlate::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"nodes\": [";
    for (int i = 0; i < nen - 1; i++)
      s << node_tags(i) << ", ";
    s << node_tags(nen - 1) << "]";
    s << "}";
    return;
  }
}
