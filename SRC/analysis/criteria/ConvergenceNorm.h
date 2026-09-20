//
// test {
//   <Quantiy> <iter> <value> <verbosity> [-override] [-sufficient]
// } 
//
// test("NormUnbalance", 10, 1e-6, 2, -override, -sufficient)
// test([
//   xara.Norm("residual|increment|energy|iteration", value, iter, verbosity=9, override, sufficient),
//   xara.Norm(quantity, iter, value, verbosity, override, sufficient, relative="predictor|applied|total"),
// ])
//
#pragma once
#include <array>
#include <vector>
#include <Vector.h>
#include <StandardStream.h>
#include <string>
#include <ConvergenceTest.h>

#include <Logging.h>

class SystemNorm : public ConvergenceTest
{
public:
  enum Status {
    Continue =-1,
    Failure  =-2
  };
  
  enum Quantity: unsigned {
    NormUnbalance              =    1<<0,
    NormDispIncr               =    1<<1,
    EnergyIncr                 =    1<<2,
    // RelativeNormUnbalance      =      40,
    // RelativeNormDispIncr       =      50,
    // RelativeEnergyIncr         =      60, // 
    // RelativeTotalNormDispIncr  =      70, // 9.49
    FixedNumIter               =      80,
  };

  struct NormData 
  {
    double tolerance;
    int    maxIter;
    std::size_t  maxIncrease=0;
    int  norm_order=2;
    bool pass_sufficient = false,
         fail_sufficient = true;
    bool override        = false;
    enum Normalize {
      None,
      NormTotal,
      StepTotal,
      StepStart,
      LoadTotal, // b only
      PathTotal  // x only
    } normalize = None;
    void print(OPS_Stream& s) const {
      s << "tol: " << tolerance
        << ", maxIter: " << maxIter
        << ", maxIncrease: " << maxIncrease
        << ", norm_order: " << norm_order
        << ", pass_sufficient: " << pass_sufficient
        << ", fail_sufficient: " << fail_sufficient
        << ", override: " << override
        << ", normalize: " << normalize
        << "\n";
    }
  };

  static constexpr int types = 10;
  SystemNorm(int  verbosity)
  : ConvergenceTest(-1),
    verbosity(0),
    active(0),
    currentIter(1),
    b0(0.0), x0(0.0), w0(0.0),
    b_increase(0), x_increase(0), w_increase(0),
    b_norms(0), x_norms(0), w_norms(0)
  {

  }

  int getType() { return active; }

  const Vector* getNorms(int type) {
    switch (type) {
      case NormUnbalance: return &b_norms;
      case NormDispIncr:  return &x_norms;
      case EnergyIncr:    return &w_norms;
      default:            return &b_norms;
    }
  }

  ConvergenceTest* getCopy(int iterations) {
    SystemNorm* copy = new SystemNorm(verbosity);
    copy->active = active;
    copy->verbosity = verbosity;
    if (active & NormUnbalance)
      copy->configure(NormUnbalance, b_opt);
    if (active & NormDispIncr)
      copy->configure(NormDispIncr, x_opt);
    if (active & EnergyIncr)
      copy->configure(EnergyIncr, w_opt);
    return copy;
  }

  void setVerbosity(int verbosity) { this->verbosity = verbosity; }

  int configure(int type, NormData& data) {
    if (type < 0 || type >= types)
      return -1;
    active |= type;
    switch (type) {
      case NormUnbalance:
        b_opt = data;
        b_norms.resize(data.maxIter+1);
        if (b_opt.normalize == NormData::PathTotal)
          return -1;
        break;
      case NormDispIncr:
        x_opt = data;
        x_norms.resize(data.maxIter+1);
        break;
      case EnergyIncr:
        w_opt = data;
        w_norms.resize(data.maxIter+1);
        break;
      default:
        return -1;
    }
    return 0;
  }

  //
  //
  int record(const int status) {
    if (verbosity != 0) {
      if (verbosity == 1 && status == Continue)
        pstream << LOG_ITERATE;
      else if (status == Failure)
        pstream << LOG_FAILURE;
      else if (status > 0)
        pstream << LOG_SUCCESS;
      else
        return status;

      pstream << "Iter: "         << pad(currentIter)
              << ", R : "         << pad(b_norms(currentIter-1))
              << ", dX: "         << pad(x_norms(currentIter-1))
              << ", dW: "         << pad(w_norms(currentIter-1)) 
              << "\n";
    }
    return status;
  }

  int getNumTests() { return currentIter; }
  int getMaxNumTests() { return 1; } // TODO
  double getRatioNumToMax() { return 1.0; } // TODO
  const Vector& getNorms() { return b_norms; } // TODO


  int start(LinearSOE& theSOE) {
    const Vector& b = theSOE.getB();

    currentIter = 1;
    b_norms.Zero();
    if (b_opt.normalize == NormData::None)
      b0  = 1.0;
    else if (b_opt.normalize == NormData::StepStart || b_opt.normalize == NormData::StepTotal)
      b0  = b.pNorm(b_opt.norm_order);
    else if (b_opt.normalize == NormData::LoadTotal)
      b0 += b.pNorm(b_opt.norm_order);
    else
      b0 = 1.0;
    return 0;
  }

  int test(const Vector& b, const Vector& x) {
    // The return value upon success is the current iteration number
    const int success = currentIter;

    const int idx = currentIter-1;
    if (currentIter == 1) {
      x_norms.Zero();
      w_norms.Zero();

      // TODO: NormTotal is supposed to normalize by the norm of the sum of x
      // over the full history of analysis. This would require test() to pass
      // dx = (x_n-1)^x 
      double dx = 0;
      if (x_opt.normalize == NormData::None)
        x0 = 1.0;
      else if (x_opt.normalize == NormData::StepStart || x_opt.normalize == NormData::StepTotal)
        x0 = x.pNorm(x_opt.norm_order);
      else if (x_opt.normalize == NormData::PathTotal)
        x0 += x.pNorm(x_opt.norm_order);
      else if (x_opt.normalize == NormData::NormTotal)
        x0 = std::sqrt(x0*x0 + (x^x) + 2.0*dx); // assuming p=2

      if (w_opt.normalize == NormData::None)
        w0 = 1.0;
      else if (w_opt.normalize == NormData::StepStart || w_opt.normalize == NormData::StepTotal)
        w0 = x^b;
      else if (w_opt.normalize == NormData::PathTotal)
        w0 += x^b;
      else if (w_opt.normalize == NormData::NormTotal)
        w0 = 1.0;// std::sqrt(w0*w0 + (x^b) + 2.0*dx);
    }

    int b_status = success;
    if (active & NormUnbalance) {
      b_norms(idx) = b.pNorm(b_opt.norm_order);
      if (b_opt.normalize == NormData::StepTotal || b_opt.normalize == NormData::LoadTotal)
        b0 += b_norms(idx);
      if (idx > 0 && b_norms(idx) > b_norms(idx-1))
        b_increase++;
      //
      //
      if (b_norms(idx) <= b_opt.tolerance*b0) {
        b_status = success;
        if (b_opt.pass_sufficient)
          return record(success);
      }
      else if (b_opt.maxIter <= currentIter) {
        b_status = Failure;
        if (b_opt.override && b_opt.fail_sufficient)
          return success;
        if (b_opt.fail_sufficient)
          return record(Failure);
      }
      else if (b_opt.maxIncrease && b_increase > b_opt.maxIncrease) {
        b_status = Failure;
        if (b_opt.override && b_opt.fail_sufficient)
          return success;
        if (b_opt.fail_sufficient)
          return record(Failure);
      }
      else {
        b_status = Continue;
      }
    }

    int x_status = success;
    if (active & NormDispIncr) {
      x_norms(currentIter) = x.pNorm(x_opt.norm_order);
      if (x_opt.normalize == NormData::StepTotal || x_opt.normalize == NormData::PathTotal)
        x0 += x_norms(currentIter);
      if (x_norms(currentIter) > x_norms(currentIter-1))
        x_increase++;
      //
      //
      if (x_norms(currentIter-1) <= x_opt.tolerance*x0) {
        x_status = success;
        if (x_opt.pass_sufficient)
          return record(success);
      }
      else if (x_opt.maxIter <= currentIter) {
        x_status = Failure;
        if (x_opt.override && x_opt.fail_sufficient)
          return success;
        if (x_opt.fail_sufficient)
          return record(Failure);
      }
      else if (x_opt.maxIncrease && x_increase > x_opt.maxIncrease) {
        x_status = Failure;
        if (x_opt.override && x_opt.fail_sufficient)
          return success;
        if (x_opt.fail_sufficient)
          return record(Failure);
      }
      else {
        x_status = Continue;
      }
    }

    int w_status = success;

    //
    //
    //
    
    if (b_status == success && x_status == success && w_status == success) {
      return record(success);
    }
    else if (b_status == Failure && x_status == Failure && w_status == Failure) {
      return record(Failure);
    }
    else {
      currentIter++;
      return record(Continue);
    }
  }




private:

  int active;
  bool override;

  int verbosity;

  int currentIter;
  double b0, x0, w0;

  NormData b_opt,   x_opt,   w_opt, c_opt;
  int b_increase, x_increase, w_increase;

  Vector b_norms, x_norms, w_norms;

  Status status;


  // std::string pad(double x);
  // std::string pad(int i);
  // StandardStream pstream;
};
