//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//

#pragma once

template <int nn> void 
constexpr lagrange(const double xi, const double xn[nn], double shp[2][nn])
{
  // constexpr int nn = 2;
  for (int i = 0; i < nn; i++) {
    shp[0][i] = 1.0;
    shp[1][i] = 0.0;

    for (int j = 0; j < nn; j++)
      if (j != i)
        shp[0][i] *= ((xi - xn[j]) / (xn[i] - xn[j]));

    for (int j = 0; j < nn; j++)
      if (j != i)
        shp[1][i] += 1.0 / (xi - xn[j]);

    shp[1][i] *= shp[0][i];

  }
}


template <int ndm, int nn, int deriv, typename T> void 
lagrange(const double xi, T shp[nn]);

// template <int ndm, int nn, int deriv, typename T> void 
// lagrange(const double xi, const double xn[nn], T shp[2][nn]);



template <> void 
lagrange<1, 2, 0>(const double xi, double shp[2])
{
    shp[0] = 0.5*(1. - xi);
    shp[1] = 0.5*(1. + xi);
}

template <> void 
lagrange<1, 2, 1>(const double xi, double shp[2])
{
    shp[0] = -0.5;
    shp[1] =  0.5;
}

template <> void 
lagrange<1, 3, 0>(const double xi, double shp[3])
{
    // one-dimensional quadratic shape functions
    //
    // o------o------o
    // 0      2      1
    shp[0] = 0.5*xi*(xi - 1.0);
    shp[1] = 0.5*xi*(xi + 1.0);
    shp[2] = 1.0 - xi*xi;
}

template <> void 
lagrange<1, 3, 1>(const double xi, double shp[3])
{
    // one-dimensional quadratic shape functions
    //
    // o------o------o
    // 0      2      1
    shp[0] =  0.5 * ( 2.0*xi - 1.0 ) ;
    shp[1] =  0.5 * ( 2.0*xi + 1.0 ) ;
    shp[2] = -2.0*xi ;

}

template <> void 
lagrange<1, 4, 0>(const double xi, double shp[4])
{
     // cubic function
     double xi2 = xi*xi;
     shp[0] = 0.0625*(1. - xi) *(9.*xi2 - 1.);
     shp[1] = 0.5625*(1. - xi2)*(1. - 3.*xi);
     shp[2] = 0.5625*(1. - xi2)*(1.+3.*xi);
     shp[3] = 0.0625*(1.+xi) *(9.*xi2 - 1.);

}

template <> void 
lagrange<1, 4, 1>(const double xi, double shp[4])
{
     // derivative of cubic function
     double xi2 = xi*xi;
     shp[0] = 0.0625*( 1. + 18.*xi - 27.*xi2);
     shp[1] = 0.5625*(-3. -  2.*xi +  9.*xi2);
     shp[2] = 0.5625*( 3. -  2.*xi -  9.*xi2);
     shp[3] = 0.0625*(-1. + 18.*xi + 27.*xi2);

}

