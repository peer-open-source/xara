#include <VectorND.h>
#include <MatrixND.h>
#include <cmath>

using namespace OpenSees;

inline double 
Heaviside(double X) { 
  return X > 1e-16 ? 1.0 : (X < 0.0 ? 0.0 : 0.0); // 0.5 
}

template <typename V3, typename M3>
int
EigSY3v0(M3& v, V3& d) 
{
  //.... compute eigenvalues and vectors for a 3 x 3 symmetric matrix
  //
  //.... INPUTS:
  //        M(3,3) - matrix with initial values (only upper half used)
  //
  //.... OUTPUTS
  //        v(3,3) - matrix of eigenvectors (by column)
  //        d(3)   - eigenvalues associated with columns of v
  //        rot    - number of rotations to diagonalize
  //
  //---------------------------------------------------------------eig3==

  //.... Storage done as follows:
  //
  //       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
  //       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
  //       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |
  //
  //        Transformations performed on d(i) and a(i) and v(i,j) become
  //        the eigenvectors.  
  //
  //---------------------------------------------------------------eig3==

    int     i, j, k;
    double  g, h, aij, thresh, t, c, s, tau;

    static const double tol = 1.0e-08;

    //.... move array into one-d arrays
    double a[3];
    a[0] = v(0, 1);
    a[1] = v(1, 2);
    a[2] = v(2, 0);

    double b[3], z[3];
    for (int i = 0; i < 3; i++) {
        d[i] = v(i, i);
        b[i] = v(i, i);
        z[i] = 0.0;

        for (int j = 0; j < 3; j++)
            v(i, j) = 0.0;

        v(i, i) = 1.0;
    }

    int rot = 0;
    int its = 0;

    double sm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);

    while (sm > tol) {
        //.... set convergence test and threshold
        if (its < 3)
            thresh = 0.011 * sm;
        else
            thresh = 0.0;

        //.... perform sweeps for rotations
        for (int i = 0; i < 3; i++) {

            j = (i + 1) % 3;
            k = (j + 1) % 3;

            aij = a[i];

            g = 100.0 * fabs(aij);

            if (fabs(d(i)) + g != fabs(d(i)) ||
                fabs(d(j)) + g != fabs(d(j))) {

                if (fabs(aij) > thresh) {

                    a[i] = 0.0;
                    h = d[j] - d[i];

                    if (fabs(h) + g == fabs(h))
                        t = aij / h;
                    else {
                        //t = 2.0 * sign(h/aij) / ( fabs(h/aij) + sqrt(4.0+(h*h/aij/aij)));
                        double hDIVaij = h / aij;
                        if (hDIVaij > 0.0)
                            t = 2.0 / (hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
                        else
                            t = -2.0 / (-hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
                    }

                    //.... set rotation parameters

                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);

                    //.... rotate diagonal terms

                    h = t * aij;
                    z[i] = z[i] - h;
                    z[j] = z[j] + h;
                    d[i] = d[i] - h;
                    d[j] = d[j] + h;

                    //.... rotate off-diagonal terms

                    h = a[j];
                    g = a[k];
                    a[j] = h + s * (g - h * tau);
                    a[k] = g - s * (h + g * tau);

                    //.... rotate eigenvectors

                    for (int k = 0; k < 3; k++) {
                        g = v(k, i);
                        h = v(k, j);
                        v(k, i) = g - s * (h + g * tau);
                        v(k, j) = h + s * (g - h * tau);
                    }

                    rot = rot + 1;

                } // end if fabs > thresh 
            }
            else
                a[i] = 0.0;

        }

        //.... update the diagonal terms
        for (i = 0; i < 3; i++) {
            b[i] = b[i] + z[i];
            d[i] = b[i];
            z[i] = 0.0;
        }

        its += 1;

        sm = std::fabs(a[0]) + std::fabs(a[1]) + std::fabs(a[2]);

    }

    // sort in descending order (unrolled bubble sort)
    auto sortij = [&d, &v](int i, int j) {
        if (d[i] < d[j]) {
            std::swap(d(i), d(j));
            for (int k = 0; k < 3; ++k)
                std::swap(v(k, i), v(k, j));
        }
    };
    sortij(0, 1);
    sortij(1, 2);
    sortij(0, 1);

    // done
    return 0;
}


static inline void 
addVoightTensorProduct(MatrixND<6,6>& P,
                       const VectorND<3>& n1,
                       const VectorND<3>& n2,
                       const VectorND<3>& n3,
                       const VectorND<3>& n4,
                       const double scale)
{
    const VectorND<6> pi = {
        n1[0]*n2[0], n1[1]*n2[1], n1[2]*n2[2],
        n1[0]*n2[1], n1[1]*n2[2], n1[2]*n2[0]
    };
    const VectorND<6> pj = {
        n3[0]*n4[0], n3[1]*n4[1], n3[2]*n4[2],
        n3[0]*n4[1], n3[1]*n4[2], n3[2]*n4[0]
    };
    P.addTensorProduct(pi, pj, scale);
}

#if 0
template <int nr>
void 
Vector2Tensor(VectorND<nr>& V)
{
//VECTOR2TENSOR transforms vector to 3x3 tensor
//  TENSOR = VECTOR2TENSOR (V)
//  the function transforms vector to a 2x2 or 3x3 tensor with the following ordering:
//  11, 22, 33, 12 (21), 23 (32) , 31 (13) (parentheses for 6x1 vector to symmetric tensor)

//  =========================================================================================
//  FEDEASLab - Release 6.0, July 2025
//  MATLAB Finite Elements for Design, Evaluation and Analysis of Structures
//  Professor Filip C. Filippou (filippou@berkeley.edu)
//  Department of Civil and Environmental Engineering, UC Berkeley
//  Copyright(c) 1998-2025. The Regents of the University of California. All Rights Reserved.
//  =========================================================================================
//  function added                                                                    03-2021
//  -----------------------------------------------------------------------------------------

// arrange in column
V = V(:);

if constexpr (nr == 3) 
    Tensor = [ V(1) V(3) ; V(3) V(2) ];
if constexpr (nr == 6)
    Tensor = [ V(1) V(4) V(6) ; 
               V(4) V(2) V(5) ; 
               V(6) V(5) V(3) ];
if constexpr (nr == 9)
    Tensor = [ V(1) V(4) V(9) ; V(5) V(2) V(6) ; V(8) V(7) V(3) ];

}
#endif

int 
StrsDecA(const VectorND<6> &sig, 
        VectorND<6> &sigpos,
        MatrixND<6,6>* Qpos)
{
    // convert to matrix form
    MatrixND<3,3> sigM; // = Vector2Tensor(sig);

	sigM(0, 0) = sig[0];
	sigM(1, 1) = sig[1];
	sigM(2, 2) = sig[2];
	sigM(0, 1) = sigM(1, 0) = sig[3];
	sigM(1, 2) = sigM(2, 1) = sig[4];
	sigM(0, 2) = sigM(2, 0) = sig[5];

    // find eigenvalues sigI and eigenvectors n
    VectorND<3> sigI;
    if (EigSY3v0(sigM, sigI) < 0)
        return -1;

    // eigenvectors are n = [n1 n2 n3]
    VectorND<3> n1, n2, n3;
    for (int i=0; i<3; i++) {
        n1[i] = sigM(i,0);
        n2[i] = sigM(i,1);
        n3[i] = sigM(i,2);
    }

    // Heaviside function
    double H1 = Heaviside(sigI[0]), 
           H2 = Heaviside(sigI[1]), 
           H3 = Heaviside(sigI[2]);;

    // positive stress

    // McCaulay of principal stresses, M1 = <sigI1>...
    double  M1 = H1*sigI[0],
            M2 = H2*sigI[1],
            M3 = H3*sigI[2];

    MatrixND<3,3> sigposM{};
    sigposM.addTensorProduct(n1, n1, M1);
    sigposM.addTensorProduct(n2, n2, M2);
    sigposM.addTensorProduct(n3, n3, M3);
    sigpos[0] = sigposM(0, 0);
    sigpos[1] = sigposM(1, 1);
    sigpos[2] = sigposM(2, 2);
    sigpos[3] = sigposM(0, 1);
    sigpos[4] = sigposM(1, 2);
    sigpos[5] = sigposM(2, 0);


    // Qpos, Qneg : linearization of stress projections, i.e.
    // d(sig_pos) = Qpos*d(sig)


    if (Qpos != nullptr) {    
        MatrixND<6,6>& Ppos = *Qpos;
        Ppos.zero();

        // Ppos : positive projection operator;
        addVoightTensorProduct(Ppos, n1,n1,  n1,n1,  H1);
        addVoightTensorProduct(Ppos, n2,n2,  n2,n2,  H2);
        addVoightTensorProduct(Ppos, n3,n3,  n3,n3,  H3);

        // p12(x)p12
        double term12 = (fabs(sigI[0] - sigI[1]) <= 1e-16) ? H1 : (M1 - M2)/(sigI[0] - sigI[1]);
        addVoightTensorProduct(Ppos, n1,n2,  n1,n2,  0.25*2*term12); // 1212
        addVoightTensorProduct(Ppos, n1,n2,  n2,n1,  0.25*2*term12); // 1221
        addVoightTensorProduct(Ppos, n2,n1,  n1,n2,  0.25*2*term12); // 2112
        addVoightTensorProduct(Ppos, n2,n1,  n2,n1,  0.25*2*term12); // 2121
        // p13(x)p13
        double term13 = (fabs(sigI[0] - sigI[2]) <= 1e-16) ? H1 : (M1 - M3)/(sigI[0] - sigI[2]);
        addVoightTensorProduct(Ppos, n1,n3,  n1,n3,  0.25*2*term13); // 1313
        addVoightTensorProduct(Ppos, n1,n3,  n3,n1,  0.25*2*term13); // 1321
        addVoightTensorProduct(Ppos, n3,n1,  n1,n3,  0.25*2*term13); // 3113
        addVoightTensorProduct(Ppos, n3,n1,  n3,n1,  0.25*2*term13); // 3121
        // p23(x)p23
        double term23 = (fabs(sigI[1] - sigI[2]) <= 1e-16) ? H2 : (M2 - M3)/(sigI[1] - sigI[2]);
        addVoightTensorProduct(Ppos, n2,n3,  n2,n3,  0.25*2*term23); // 2323
        addVoightTensorProduct(Ppos, n2,n3,  n3,n2,  0.25*2*term23); // 2332
        addVoightTensorProduct(Ppos, n3,n2,  n2,n3,  0.25*2*term23); // 3223
        addVoightTensorProduct(Ppos, n3,n2,  n3,n2,  0.25*2*term23); // 3232

        for (int i=0; i<6; i++)
           for (int j=3; j<6; j++)
               Ppos(i,j) = 2*Ppos(i,j); // 0.5*Ppos(i,j) + 0.5*Ppos(j,i);
    }
    return 0;
    // Qpos = Ppos + 2*(term12*(p12*p12') + term13*(p13*p13') + term23*(p23*p23'));
    // Qpos(:,4:6) = 2*Qpos(:,4:6);
}