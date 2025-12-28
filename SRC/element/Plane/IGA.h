//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// https://github.com/li-ming-jiang/MyOpenSees
//

namespace IGA {

int    FindSpan(int order, int ncp, int eleId, int nele, Vector& KnotVect);
void   DerBasisFuns(double Idx, Vector Pts, int obf, int n, Vector KnotVect, Matrix** N0n);
void   calcDersBasisFunsAtGPs(int obf, int ncp, Vector KnotVect, int d, int NGPs, int Idx, double* J2, Vector* Xg, Matrix** N0n);
void   Rationalize(Vector WeightsCP, Vector N0, Matrix N1, Vector* R0, Matrix* R1);


void   GaussRule(int NGPs1, Vector* Xg, Vector* WXg);


Vector Dot(const Vector& VecA, const Vector& VecB);
Vector DotDiv(const Vector& VecA, const Vector& VecB);

Vector
DotDiv(const Vector& VecA, const Vector& VecB)
{
    Vector result(VecA.Size());
    for (int i = 1; i <= VecA.Size();i++) {
        result(i - 1) = VecA(i - 1) / VecB(i - 1);
    }
    return result;
}

Vector
Dot(const Vector& VecA, const Vector& VecB)
{
    Vector result(VecA.Size());
    for (int i = 1; i <= VecA.Size();i++) {
        result(i - 1) = VecA(i - 1) * VecB(i - 1);
    }
    return result;
}

#if 0
void 
GaussRule(int NGPs1, Vector* Xg, Vector* WXg)
{
    double NGPs = (double)NGPs1;
    Vector X(NGPs), W(NGPs);
    Matrix result(NGPs, 2);
    //Matrix ptss(NGPs, 2);
    //Vector wtss(NGPs);
    //GaussRule(NGPs);
    int iter = 2;
    double e1 = NGPs * (NGPs + 1);
    int mm = trunc((NGPs + 1) / 2);
    Vector tt(mm);
    int kk = 1;
    double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
    for (int i = 3; i <= 4 * mm - 1;i += 4) {
        tt(kk - 1) = i * PI / (4 * NGPs + 2);
        kk++;
    }
    Vector x0(tt.Size());
    Vector h(tt.Size()), 
           d1(tt.Size()), pk(tt.Size()), dpn(tt.Size()), d2pn(tt.Size()), d3pn(tt.Size()), d4pn(tt.Size());
    for (int i = 1; i <= tt.Size();i++) {
        x0(i - 1) = (1 - (1 - 1 / NGPs) / (8 * NGPs * NGPs)) * cos(tt(i - 1));
    }

    for (int j = 1; j <= iter;j++) {
        Vector pkm1(x0.Size());
        pkm1.Zero();
        pkm1 += 1;
        pk = x0;
        for (int k = 2; k <= NGPs; k++) {
            Vector pkp1 = 2 * Dot(x0, pk) - pkm1 - (Dot(x0, pk) - pkm1) / k;
            pkm1 = pk;
            pk = pkp1;
        }

        Vector den = -1 * Dot(x0, x0) + 1;
        d1   = NGPs * (pkm1 - Dot(x0, pk));
        dpn  = DotDiv(d1, den);
        d2pn = DotDiv((2 * Dot(x0, dpn) - e1 * pk), den);
        d3pn = DotDiv((4 * Dot(x0, d2pn) + (2 - e1) * dpn), den);
        d4pn = DotDiv((6 * Dot(x0, d3pn) + (6 - e1) * d2pn), den);
        Vector uu = DotDiv(pk, dpn);
        Vector vv = DotDiv(d2pn, dpn);

        // Initial approximation H
        h = Dot(-1 * uu, 0.5 * Dot(uu, vv + Dot(uu, DotDiv(Dot(vv, vv) - Dot(uu, d3pn), 3 * dpn))) + 1);

        // Refine H using one step of Newton's method
        Vector p = pk + Dot(h, dpn + 0.5 * Dot(h, d2pn + Dot(h / 3, d3pn + 0.25 * Dot(h, d4pn))));
        x0 += h;
    }
    for (int i = 1; i <= x0.Size();i++) { 
      X(i - 1) = -1 * x0(i - 1) - h(i - 1); 
    }
    Vector fx = d1 - Dot(h, e1 * (pk + 0.5 * Dot(h, dpn + Dot(h / 3, d2pn + 0.25 * Dot(h, d3pn + 0.2 * Dot(h, d4pn))))));
    for (int i = 1;i <= x0.Size();i++) { 
      W(i - 1) = 2 * (1 - X(i - 1) * X(i - 1)) / (fx(i - 1) * fx(i - 1)); 
    }
    if (mm + mm > NGPs) {
      X(mm - 1) = 0;
    }
    if (mm + mm != NGPs) {
      mm = mm - 1;
    }
    for (int i = 1;i <= mm;i++) {
        X(NGPs - i) = -X(i - 1);
        W(NGPs - i) =  W(i - 1);
    }
    for (int i = 1; i <= NGPs; i++) {
        (*Xg)(i - 1)  = X(i - 1);
        (*WXg)(i - 1) = W(i - 1);
    }
}
#endif

// calculate the shape function and derivatives of the quadrature point
int
FindSpan(int order, int ncp, int eleId, int nele, Vector &KnotVect) 
{
    ID Idxs(nele);
    Idxs.Zero();
    int iE = 1;
    for (int i = order; i <= ncp; i++) {
        if (fabs(KnotVect(i - 1) - KnotVect(i)) > 1.4901e-08) {
            Idxs(iE - 1) = i;
            iE = iE + 1;
        }
    }
    int Idx = Idxs(eleId - 1);
    return Idx;
}


// The shape function should be revised for IGA, this is the 1d code
void
DerBasisFuns(double Idx, Vector Pts, int obf, int n, Vector KnotVect, Matrix** N0n)
{
  //
  // Idx: span index of the parametric location; 
  // Pts: parametric points; 
  // obf: order of basis function;
  // n: maximum order of derivatives; 
  // KnotVect: knot vector
  //
    int    PtsSize = obf + 1;
    double iXi;
    double saved, temp;
    double d;

    //N0n.resize(PtsSize,PtsSize);
    //N0n.Zero();
    Matrix Ni(obf + 1, n + 1);
    Matrix ndu(obf + 1, obf + 1);
    Vector left(obf+1); // more 1 space to avoid error

    Vector right(obf+1);
    Matrix a(2, obf + 1);
    for (int i = 1; i <= PtsSize; i++) {
        iXi = Pts(i - 1);
        ndu(0, 0) = 1.0;
        for (int j = 1; j <= obf; j++) {
            left(j) = iXi - KnotVect(Idx - j);
            right(j) = KnotVect(Idx + j - 1) - iXi;
            saved    = 0;
            for (int r = 0; r <= j - 1; r++) {
                // lower triangle
                ndu(j, r) = right(r + 1) + left(j - r);
                temp = ndu(r, j - 1) / ndu(j, r);
                // upper triangle
                ndu(r, j) = saved + right(r + 1) * temp;
                saved = left(j - r) * temp;
            }
            ndu(j, j) = saved;
        }
        for (int r = 1; r <= obf + 1;r++) {
            Ni(r - 1, 0) = ndu(r - 1, obf);
        }

        // Ni.setData(ndu(:, p),[],0); //the syntax is wrong
        //
        for (int r = 0; r <= obf; r++) {
            int s1 = 0,
                s2 = 1;

            a(0, 0) = 1.0;
            for (int k = 1; k <= n; k++) {
                d = 0;
                int rk = r - k;
                int pk = obf - k;
                if (r >= k) {
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                    d = a(s2, 0) * ndu(rk, pk);
                }

                int j1, j2;

                if (rk >= -1)
                    j1 = 1;
                else
                    j1 = -rk;

                if ((r - 1) <= pk)
                    j2 = k - 1;
                else
                    j2 = obf - r;

                for (int j = j1; j <= j2; j++) {
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                    d       += a(s2, j) * ndu(rk + j, pk);
                }
                if (r <= pk) {
                    a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                    d       += a(s2, k) * ndu(r, pk);
                }
                Ni(r, k) = d;
                int j = s1;
                s1 = s2;
                s2 =  j; // switch rows
            }
        }

        int r = obf;
        for (int k = 1;k <= n;k++) {
            for (int rr = 1; rr <= obf + 1; rr++) {
                Ni(rr - 1, k) *= r;// the syntax is also wrong
            }
            r = r * (obf - k);
        }
        (*N0n)[i-1] = Ni;
    }
}

// to calculate the shape functions and derivatives of Gaussian points
// 需要修改成单个单元的版本！
void 
calcDersBasisFunsAtGPs(int obf, int ncp, Vector KnotVect, int d, int NGPs, int Idx, double* J2, Vector* WXg, Matrix** N0n)
{
    // obf: order of basis function, ncp: number of control points;
    // KnotVect: knot vector; d: degree of derivative; NGPs: number of gauss points; Idx: span index of the element;
    Vector Xg(NGPs);
    GaussRule(NGPs,&Xg,WXg);
    //Vector Xg(NGPs);
    
    //for (int i = 1; i <= NGPs; i++) {
    //    Xg(i - 1) = resultGauss(i - 1, 0);
    //    WXg(i - 1) = resultGauss(i - 1, 1);
    //}
    //double J2;
    Vector Xi_e(NGPs);
    //Matrix* N = new Matrix [NGPs];
    //N = new Matrix * [NEl];
    //for (int i = 1; i <= NEl;i++) {
    //    N[i-1] = new Matrix [NGPs];// N = zeros(NEl,NGPs,obf+1,d+1)
    //}
    //int i = Idx[ex - 1];
    *J2 = (KnotVect(Idx) - KnotVect(Idx - 1)) / 2;
    //Vector Xi_e(NGPs);
    for (int ii = 1;ii <= NGPs;ii++) {//有问题！！！
        Xi_e(ii - 1) = (Xg(ii - 1) + 1) * (*J2) + KnotVect(Idx - 1);//Xi(ex,:)=(Xg+1)*J2(ex)+KnotVect(i);
        //Xi_ex(ii - 1) = Xi(ex - 1, ii - 1);
    }
    DerBasisFuns(Idx, Xi_e, obf, d, KnotVect, N0n);//N(ex,:,:,:) = DersBasisFuns(i,Xi(ex,:),obf,d,KnotVect)
    //resultDersBasisFunsAtGPS result;// return a struct
    //result.J2 = J2;
    //result.WXg = WXg;
    //result.Xi = Xi_e;
    //result.N = N0n;
}

void
Rationalize(Vector WeightsCP, Vector N0, Matrix N1, Vector* R0, Matrix* R1)
{
    // convert B-spline to NURBS basis functions for 1D
    Vector N0W(N0.Size());
    double W0 = 0;
    for (int i = 1; i <= N0.Size(); i++) {
        N0W(i - 1) = N0(i - 1) * WeightsCP(i - 1);
        W0 += N0W(i - 1);
    }
    *R0 = N0W / W0;
    //for (int i = 1; i <= R0.Size();i++) {
    //    (RatFunc[0])(i - 1, 0) = R0(i - 1);
    //}
    // First derivatives of NURBS basis functions
    Matrix N1W(N1.noRows(), N1.noCols());
    for (int i = 1; i <= N1.noRows(); i++) {
        for (int j = 1;j <= N1.noCols();j++) {
            N1W(i - 1, j - 1) = N1(i - 1, j - 1) * WeightsCP(j - 1);
        }
    }
    Vector W1(N1.noRows());
    Matrix Temp(N1W.noRows(), N1W.noCols());
    for (int i = 1; i <= N1W.noRows();i++) {
        for (int j = 1; j <= N1W.noCols();j++) {
            W1(i - 1) += N1W(i - 1, j - 1);
        }
    }
    for (int i = 1; i <= N1W.noRows();i++) {
        for (int j = 1; j <= N1W.noCols();j++) {
            Temp(i - 1, j - 1) = (*R0)(j - 1) * W1(i - 1);
        }
    }
    Matrix N1WR0W1 = N1W - Temp;
    *R1 = N1WR0W1 / W0;

    //result.R0 = R0;
    //result.R1 = R1;
}

//double shapeFunction1(int obfs, int ncps, Vector* KnotVects, int nds, int NGPss, int Idxs)
//{
//    //calculate the shape function and derivatives of x and y directions
//    int obf_x = obfs[0]; int obf_y = obfs[1];
//    int ncp_x = ncps[0]; int ncp_y = ncps[1];
//    Vector KnotVect_x = KnotVects[0]; Vector KnotVect_y = KnotVects[1];
//    int nd_x = nds[0]; int nd_y = nds[1];
//    int NGPs_x = NGPss[0]; int NGPs_y = NGPss[1];
//    int Idx_x = Idxs[0]; int Idx_y = Idxs[1];
//
//    resultDersBasisFunsAtGPS ValueBFG_x = calcDersBasisFunsAtGPs(obf_x, ncp_x, KnotVect_x, nd_x, NGPs_x, Idx_x);// return basis function and derivatives of GP
//    resultDersBasisFunsAtGPS ValueBFG_y = calcDersBasisFunsAtGPs(obf_y, ncp_y, KnotVect_y, nd_y, NGPs_y, Idx_y);
//    double Jx = ValueBFG_x.J2; Vector Wx = ValueBFG_x.WXg; Matrix* Nx = ValueBFG_X.N;
//    double Jy = ValueBFG_y.J2; Vector Wy = ValueBFG_y.WXg; Matrix* Ny = ValueBFG_X.N;
//    int NEN = (obf_x + 1) * (obf_y + 1); // number of local basic functions
//    Vector N0(NEN);
//    Matrix N1(2, NEN);
//    double WeightsCP1[] = { 1,1 };
//    Vector WeightsCP(WeightsCP1, 2);
//    const Vector& nd1Crds = theNodes[0]->getCrds();
//    const Vector& nd2Crds = theNodes[1]->getCrds();
//    const Vector& nd3Crds = theNodes[2]->getCrds();
//    const Vector& nd4Crds = theNodes[3]->getCrds();
//    Matrix CtrlPts(4, 3);
//    for (int i = 1; i <= 4; i++) {
//        for (int j = 1; j <= 3;j++) {
//            Vector& ndCrds = theNodes[i - 1]->getCrds();
//            CtrlPts(i - 1, j - 1) = *ndCrds(j - 1);
//        }
//    }
//    Matrix Ke(NEN * 2, NEN * 2);
//    for (int qy = 1; qy <= NGPs_y;qy++) {
//        for (int qx = 1;qx <= NGPs_x;qx++) {
//            int k = 1;
//            for (int j = 1;j <= obf_y + 1;j++) {
//                for (int i = 1; i <= obf_x + 1;i++) {
//                    N0(k - 1) = (Nx[qx - 1])(i - 1, 0) * (Ny[qy - 1])(j - 1, 0);//shape function
//                    N1(0, k - 1) = (Nx[qx - 1])(i - 1, 1) * (Ny[qy - 1])(j - 1, 0);// 1st derivatives
//                    N1(1, k - 1) = (Nx[qx - 1])(i - 1, 0) * (Ny[qy - 1])(j - 1, 1);
//                    k++;
//                }
//            }
//        }
//    }
//    //R1 = Rationalize()
//}

}