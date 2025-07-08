#include <Mate.h>
#include <NosbBase.h>
#include <Logging.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Vector3D.h>
#include <Matrix3D.h>


template <int ndim, int maxfam>
NosbProj<ndim, maxfam>::NosbProj(PeriParticle<ndim> *center, PeriDomain<ndim> &domain, Mate<ndim> *material)
    : center(center), numfam(center->numfam)
{
    delta = domain.delta;
    for (int i = 0; i < numfam; i++)
    {
        neigh[i] = &domain.pts[center->nodefam[i]];
        materials[i] = material->getCopy();
    }
}

template <int ndim, int maxfam>
void
NosbProj<ndim, maxfam>::bond_omega(int i, int w_type)
{
    // in this formulation, omega is constant 
	// for computing all the bond-associated quantities
    double deno = 0.0;
    
    omega.fill(0.0);
    for (int ind=0; ind<numfam; ind++)
    {
        const VectorND<ndim> xi = neigh[ind]->coord - center->coord;
        double dist = xi.norm();
        double x0 = dist / delta;
        if (w_type == 1){
            omega[ind] = 1.0;
        }else if (w_type == 2)
        {
            omega[ind] = std::exp(-std::pow(3.0*x0, 2));
        }
        deno += omega[ind] * center->vol[ind];
    }
    for (int ind=0; ind<numfam; ind++)
    {
        omega[ind] /= deno;
    }
    // return 1.0; // omega[i];
}

template <int ndim, int maxfam>
void 
NosbProj<ndim, maxfam>::init_shape()
{

    bond_omega(0, 1);

    // Integrate K
    // note that Kinv is invariant in this formulation, but may not be in general
    Kinv.zero();
    MatrixND<ndim, ndim> K{};
    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        K.addTensorProduct(xi, xi, omega[i] * center->vol[i]);
        // K += omega[i] * xi.bun(xi) * center->vol[i];
    }

    // Invert K
    K.invert(Kinv);
}

template <int ndim, int maxfam>
MatrixND<ndim, ndim>
NosbProj<ndim, maxfam>::get_A(const VectorND<ndim> &xi)
{
    // Initialize with zeros
    MatrixND<ndim, ndim> A{{{0.0}}};

    A.addDiagonal(1.0);

    A.addTensorProduct(xi, xi, -1.0 / xi.dot(xi));

    return A;
}

template <int ndim, int maxfam>
MatrixND<ndim, ndim>
NosbProj<ndim, maxfam>::get_B(const VectorND<ndim> &xi, const VectorND<ndim> &zeta)
{
    // Initialize with zeros
    MatrixND<ndim, ndim> B{{{0.0}}};

    B.addTensorProduct(zeta, xi, 1.0 / xi.dot(xi));

    return B;
}

template <int ndim, int maxfam>
VectorND<ndim>
NosbProj<ndim, maxfam>::get_T2(int i, const VectorND<ndim> &xi)
{
    const MatrixSD<ndim>& P = materials[i]->getStress();
    VectorND<ndim> B = P * xi;
    B /= xi.dot(xi);
    return B;
}

template <int ndim, int maxfam>
void
NosbProj<ndim, maxfam>::form_trial()
{
    // update materials (materials)

    MatrixND<ndim, ndim> Nmat{{{0.0}}};

    for (int i = 0; i < numfam; i++)
    {
        // PeriParticle<ndim> &other = *neigh[i];
        const VectorND<ndim> eta = neigh[i]->disp - center->disp;
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const VectorND<ndim> zeta = xi + eta;
        Nmat.addTensorProduct(zeta, xi, omega[i] * center->vol[i]);
        // Nmat += omega[i] * zeta.bun(xi) * center->vol[i];
    }
    MatrixND<ndim, ndim> Fmat = Nmat * Kinv;

    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> eta = neigh[i]->disp - center->disp;
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const VectorND<ndim> zeta = xi + eta;

        const MatrixND<ndim, ndim> Amat = this->get_A(xi);
        const MatrixND<ndim, ndim> Bmat = this->get_B(xi, zeta);
        // Fmat <- Fmat * (I - xi\otimes xi / |xi|^2) + zeta\otimes xi / |xi|^2
        materials[i]->setTrialStrain(Fmat * Amat + Bmat);
        // -------- FOR DEBUGGING --------
        // -------- check F --------------
        // MatrixND<ndim, ndim> Fmat2 = Fmat * Amat + Bmat;
        // for (int j = 0; j < ndim; j++)
        //     for (int k = 0; k < ndim; k++) 
        //         printf("%7.2e ", Fmat2(j, k));
        // printf("\n");
        // -------- check E ---------------
        // MatrixND<ndim, ndim> Fmat2 = Fmat * Amat + Bmat;
        // MatrixSD<ndim> E{0.0};
        // for (int j = 0; j < ndim; j++)
        //     for (int k = 0; k < ndim; k++) 
        //         E.ref(j, k) = 0.5 * Fmat2(j, k) + 0.5 * Fmat2(k, j) - (j == k);

        // for (int j = 0; j < ndim; j++)
        //     for (int k = 0; k < ndim; k++)
        //         printf("%7.2e ", E(j, k));
        // printf("\n");
        // -------- check S --------------
        // const MatrixSD<ndim>& Smat = materials[i]->getStress();
        // printf("%7.2e %7.2e %7.2e\n", Smat(0, 0), Smat(1, 1), Smat(0, 1));
        // -------------------------------
        
    }
    
}

template <int ndim, int maxfam>
MatrixND<ndim, ndim>
NosbProj<ndim, maxfam>::sum_PKinv()
{

    MatrixND<ndim, ndim> Qmat{{{0.0}}};
    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const MatrixSD<ndim>& P = materials[i]->getStress();
        const MatrixND<ndim, ndim> A = this->get_A(xi);
        // the correct formula for Qmat is 
        // Qmat += sum_1^numfam { w[i] * P * A * Kinv * vol[i] }
        // Since w[i] and omega[i] are the same in this formulation, we can use omega[i] instead
        Qmat.addMatrix(P * A * Kinv, omega[i] * center->vol[i]);


//      opserr  << Matrix(Qmat) << "\n";
    }
    return Qmat;
}

template <int ndim, int maxfam>
VectorND<ndim>
NosbProj<ndim, maxfam>::bond_force(int i, MatrixND<ndim, ndim> &Qmat)
{
    const VectorND<ndim> xi = neigh[i]->coord - center->coord;
    VectorND<ndim> T2 = this->get_T2(i, xi);

    // printf("T2=(%7.2e, %7.2e)\n", T2[0], T2[1]);

    // the correct formula for T[i] is
    // T[i] = omega[i] * Qmat * xi + w[i]*Pmat*xi/|xi|^2
    // Since w[i] and omega[i] are the same in this formulation and T2 = Pmat*xi/|xi|^2,
    // we can use the following formula instead:
    // T[i] = omega[i] * (Qmat * xi + T2)
    // Because the peridynamics force density at writes as
    // pforce = sum_1^numfam { T[i] * vol[i] }
    // T[i] is returned as 
    // T[i] = omega[i] * vol[i] * correction[i] * (Qmat * xi + T2)
    VectorND<ndim> T_j = Qmat * xi + T2;

    double w = omega[i] * center->vol[i] * center->correction[i];

    return T_j * w;
}
