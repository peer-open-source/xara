/* Provided under the Zero-Clause BSD (0BSD) License
 *
 * Copyright (C) 2021 by Claudio Perez
 *
 * Permission to use, copy, modify, and/or distribute this software
 * for any purpose with or without fee is hereby granted.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#define _USE_MATH_DEFINES
#include <tgmath.h>

#ifndef NAN
// #define S_NAN (1.0/0.0)
#define S_NAN sqrt(-1.0)
#else
#define S_NAN NAN
#endif


// Definitions used by Beta func.
#define STOP 1.0e-8 
#define TINY 1.0e-30


/* Weibull PDF
 *
 * f(x|a,b) = (b/a)(x/a)^(b-1) exp(-(x/a)^b)
 */
double wblpdf(double x, double a, double b);
double wblcdf(double x, double a, double b);

/* 
 * Lognormal PDF
 */
double lognpdf(double x, double mean, double stddev);
double logncdf(double x, double mean, double stddev);

/* 
 * PDF for the 2-parameter Beta distribution.
 */
double betapdf(double x, double a, double b);
/* 
 * CDF for the 2-parameter Beta distribution.
 *
 * F(x; a,b) = B(x; a,b) / B(a, b) = I_x(a,b)
 */
double betacdf(double x, double a, double b);
double incbeta(double x, double a, double b);
double beta(double a, double b);

int beta_isvalid(double x, double a, double b) {
  return ((x >= 0) && (x <= 1.0) && (a > 0.0) && (b > 0.0));
}

/* Probability density function for the 2-parameter
 * Beta distribution.
 */
double betapdf(double x, double a, double b) {
  if ((x >= 0) && (x <= 1.0))
    return pow(x, a - 1) * pow(1 - x, b - 1) / beta(a, b);
  else
    // return NAN;
    return 0.0;
}


/* Cumulative distribution function for the 2-parameter
 * Beta distribution.
 *
 * F(x; a,b) = B(x; a,b) / B(a, b) = I_x(a,b)
 */
double betacdf(double x, double a, double b) {
  if (x < 0.0) {
    return 0.0;
  } else if (x <= 1.0) {
    double ans = incbeta(x, a, b); ///beta(a,b);
    return ans > 1.0 ? 1.0 : ans;
  } else {
    return 1.0;
  }
}

double betacdf_0(double x, double a, double b)
// Returns the incomplete beta function Ix(a, b).
{
  double bt;
  if (x < 0.0 || x > 1.0) 
    return NAN;

  if (x == 0.0 || x == 1.0)
    bt=0.0;

  else // Factors in front of the continued fraction.
    bt = exp(lgamma(a+b) - lgamma(a) - lgamma(b) + a*log(x) + b*log(1.0 - x));

  if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly.
    return bt*incbeta(x,a,b)/a;

  else // Use continued fraction after making the symmetry transformation.
    return 1.0-bt*incbeta(1.0-x, b, a)/b; 
}

/* 
 * Beta function
 */
inline double beta(double a, double b) {
 return tgamma(a) * tgamma(b) / tgamma(a + b);
  // return exp(lgamma(a)+lgamma(b)-lgamma(a+b));
}


/* 
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * Obtained under the zlib license at the following link:
 * https://github.com/codeplea/incbeta/blob/master/incbeta.c
 */
double incbeta(double x, double a, double b) {
  if (x < 0.0 || x > 1.0)
    return -1; //NAN;

  /* The continued fraction converges nicely for x < (a+1)/(a+b+2) */
  if (x > (a + 1.0) / (a + b + 2.0)) {
    /* Use the fact that beta is symmetrical. */
    return (1.0 - incbeta(1.0 - x, b, a)); 
  }

  /* Find the first part before the continued fraction. */
  const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
  const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

  /* Use Lentz's algorithm to evaluate the continued fraction. */
  double f = 1.0, c = 1.0, d = 0.0;

  for (int i = 0; i <= 200; ++i) {
    int m = i / 2;

    double numerator;
    if (i == 0) {
      numerator = 1.0; /*First numerator is 1.0.*/
    } else if (i % 2 == 0) {
      numerator = (m * (b - m) * x) /
                  ((a + 2.0 * m - 1.0) * (a + 2.0 * m)); /*Even term.*/
    } else {
      numerator = -((a + m) * (a + b + m) * x) /
                  ((a + 2.0 * m) * (a + 2.0 * m + 1)); /*Odd term.*/
    }

    /* Do an iteration of Lentz's algorithm. */
    d = 1.0 + numerator * d;
    if (fabs(d) < TINY)
      d = TINY;
    d = 1.0 / d;

    c = 1.0 + numerator / c;
    if (fabs(c) < TINY)
      c = TINY;

    const double cd = c * d;
    f *= cd;

    /* Check for stop. */
    if (fabs(1.0 - cd) < STOP) {
      return front * (f - 1.0);
    }
  }
  return NAN; /*Needed more loops, did not converge.*/
}


/* BSD License 
 * http://www.johndcook.com/cpp_phi.html */
double phi(double x) {
  // constants
  double a1 = 0.254829592;
  double a2 = -0.284496736;
  double a3 = 1.421413741;
  double a4 = -1.453152027;
  double a5 = 1.061405429;
  double p = 0.3275911;

  // Save the sign of x
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x) / sqrt(2.0);

  // A&S formula 7.1.26
  double t = 1.0 / (1.0 + p * x);
  double y =
      1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

  return 0.5 * (1.0 + sign * y);
}

/* https://stackoverflow.com/a/3512564 */
double cnd_manual(double x) {
  double L, K, w;
  /* constants */
  double const a1 =  0.31938153,  a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

  L = fabs(x);
  K = 1.0 / (1.0 + 0.2316419 * L);
  w = 1.0 - 1.0 / sqrt(2 * M_PI) * exp(-L * L / 2) *
                (a1 * K + a2 * K * K + a3 * pow(K, 3) + a4 * pow(K, 4) +
                 a5 * pow(K, 5));

  if (x < 0) {
    w = 1.0 - w;
  }
  return w;
}


double std_normpdf(double x) {
  return exp(-pow(x, 2) * 0.5) / sqrt(2.0 * M_PI);
}

double std_normcdf(double x) {
  return 0.5 * erfc(-x * M_SQRT1_2); 
}

double std_lognpdf(double x, double shape) {
  return exp(-pow(log(x), 2) / (2.0 * shape * shape)) /
         (2.0 * x * sqrt(2.0 * M_PI));
}

double normpdf(double x, double mean, double stddev) {
  return std_normpdf(x - mean) / stddev;
}

double normcdf(double x, double mean, double stddev) {
  return std_normcdf((x - mean) / stddev);
}

double lognpdf(double x, double mean, double stddev) {
  double scale = exp(mean);
  return std_lognpdf((x - mean) / scale, stddev) / scale;
}

double logncdf(double x, double mean, double stddev) {
  return normcdf(log(x), mean, stddev);
}

/* ----------------------------------------------------------------------*
 * | WEIBULL
 * ----------------------------------------------------------------------*/

static inline int weibull_isvalid(double x, double scale, double shape) {
  return (int)((x > 0.0) && (shape > 0.0) && (scale > 0.0));
}

/* Weibull PDF
 *
 * f(x|a,b) = (b/a)(x/a)^(b-1) exp(-(x/a)^b)
 *
 * where a is the "scale" parameter, and b the "shape" parameter.
 */
double wblpdf(double x, double scale, double shape) {
  if (weibull_isvalid(x, scale, shape))
    return shape / scale * pow(x / scale, shape - 1) *
           exp(-pow(x / scale, shape));
  else
    return S_NAN;
}

double wblcdf(double x, double a, double b) {
  if (weibull_isvalid(x, a, b))
    return -expm1(-pow(x / a, b));
  else
    return S_NAN;
}
