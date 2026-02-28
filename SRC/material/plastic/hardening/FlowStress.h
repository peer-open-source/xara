#pragma once
#include <cmath>
#include <tuple>

// template <typename Resistance, typename Backstress>
// class MisesFlow {
// public:
//   MisesFlow(double G, Resistance Y, Backstress X)
//   : G(G), Y(Y), X(X) {}

//   std::tuple<double,double>
//   operator()(const VectorND<9>& s, const VectorND<9>& n, double lam) const {
//     const auto [y, dy] = Y(lam);
//     const auto [x, dx] = X(lam);
//     return {f, df};
//   }
// private:
//   double G; // shear modulus
//   Resistance Y; // flow stress
//   Backstress X; // kinematic hardening
// };

namespace Resistance {

class Linear {
public:
  Linear(double H) : H(H) {}
  std::tuple<double,double> operator()(double kappa) const {
    return {H*kappa, H};
  }
private:
  double R(double kappa) const { return H*kappa; }
  double dR(double kappa) const { return H; }
private:
  double H; // hardening modulus
};

class Voce {
public:
  Voce(double b, double Q) : b(b), Q(Q) {}
  std::tuple<double,double> operator()(double kappa) const {
    return {R(kappa), dR(kappa)};
  }
private:
  double R(double kappa) const { return Q*(1.0 - std::exp(-b*kappa)); }
  double dR(double kappa) const { return Q*b*std::exp(-b*kappa); }

private:
  double b; // hardening rate
  double Q; // saturation stress
};



class Mroz {
public:
  Mroz(double b, double Q, double n) : b(b), Q(Q), n(n) {}
  std::tuple<double,double> operator()(double kappa) const {
    return {R(kappa), dR(kappa)};
  }
private:
  double R(double kappa) const { return Q*(1.0 - std::exp(-b*std::pow(kappa, n))); }
  double dR(double kappa) const { return Q*b*n*std::pow(kappa, n-1)*std::exp(-b*std::pow(kappa, n)); }
private:
  double b; // hardening rate
  double Q; // saturation stress
  double n; // initial hardening exponent
};

} // namespace Flow