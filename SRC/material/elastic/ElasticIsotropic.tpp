//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <MatrixND.h>
#include <OPS_Stream.h>

using namespace OpenSees;

template <int ndim, PlaneType type>
ElasticIsotropic<ndim, type>::ElasticIsotropic(int tag, double E, double nu, double rho)
    : Mate<ndim>(tag), E(E), nu(nu), rho(rho)
{
  revertToStart();
}


template <int ndim, PlaneType type>
const char *
ElasticIsotropic<ndim, type>::getClassType() const
{
    return "ElasticIsotropic";
}


template <int ndim, PlaneType type>
int
ElasticIsotropic<ndim, type>::revertToStart()
{
    if constexpr (ndim == 3)
    {
        double tmp = E / (1.0 + nu) / (1.0 - 2.0 * nu);
        ddsdde.ref(0, 0) = (1.0 - nu) * tmp;
        ddsdde.ref(1, 1) = ddsdde.ref(1, 1);
        ddsdde.ref(2, 2) = ddsdde.ref(1, 1);
        ddsdde.ref(0, 1) = nu * tmp;
        ddsdde.ref(0, 2) = ddsdde.ref(1, 2);
        ddsdde.ref(1, 2) = ddsdde.ref(1, 2);
        ddsdde.ref(1, 0) = ddsdde.ref(1, 2);
        ddsdde.ref(2, 0) = ddsdde.ref(1, 2);
        ddsdde.ref(2, 1) = ddsdde.ref(1, 2);
        ddsdde.ref(3, 3) = (1.0 - 2.0 * nu) * tmp;
        ddsdde.ref(4, 4) = ddsdde.ref(4, 4);
        ddsdde.ref(5, 5) = ddsdde.ref(4, 4);
    }
    else if constexpr (ndim == 2 && type == PlaneType::Strain)
    {
        double tmp = E / (1.0+nu) / (1.0-2.0*nu);
        ddsdde.ref(0, 0) = (1.0-nu) * tmp;
        ddsdde.ref(1, 1) = (1.0-nu) * tmp;
        ddsdde.ref(0, 1) = nu * tmp;
//      ddsdde.ref(1, 0) = nu * tmp;
        ddsdde.ref(2, 2) = (1.0 - 2.0 * nu) * 0.5 * tmp;

        /*
        double mu2 = E/(1.0+nu);
        double lam = nu*mu2/(1.0-2.0*nu);
        double mu = 0.50*mu2;
        mu2 += lam;
        
        ddsdde.ref(0,0) = ddsdde.ref(1,1) = mu2;
        ddsdde.ref(0,1) = ddsdde.ref(1,0) = lam;
        ddsdde.ref(2,2) = mu;
	*/
    }
    else if constexpr (ndim == 2 && type == PlaneType::Stress)
    {
        double tmp = E / (1.0 - nu * nu);
        ddsdde.ref(0, 0) = tmp;
        ddsdde.ref(1, 1) = tmp;
        ddsdde.ref(2, 2) = (1.0 - nu)/2.0 * tmp;
        ddsdde.ref(0, 1) = nu * tmp;
//      ddsdde.ref(1, 0) = nu * tmp;

	/*
        double d00 = E/(1.0-nu*nu);
        double d01 = nu*d00;
        double d22 = 0.5*(d00-d01);

        ddsdde.ref(0,0) = ddsdde.ref(1,1) = d00;
        ddsdde.ref(1,0) = ddsdde.ref(0,1) = d01;
        ddsdde.ref(2,2) = d22;
	*/
    }
    return 0;
}


template <int ndim, PlaneType type>
Mate<ndim> *
ElasticIsotropic<ndim, type>::getCopy()
{
    return new ElasticIsotropic<ndim, type>(this->getTag(), E, nu, rho);
}

template <int ndim, PlaneType type>
const MatrixSD<ndim> &
ElasticIsotropic<ndim, type>::getStress()
{
    return Smat;
}

template <int ndim, PlaneType type>
int ElasticIsotropic<ndim, type>::setTrialStrain(const MatrixND<ndim, ndim> &F)
{
    if constexpr (false)
    {
        // Convert deformation gradient F to finite Green-Lagrange tensor e
        MatrixSD<ndim, true> e;
        e.addMatrixTransposeProduct(0.0, F, F, 0.5);
        e.addDiagonal(-0.5);
        return setTrialStrain(e);
    }
    else
    {
        // Convert deformation gradient F to small strain tensor e
        MatrixSD<ndim, true> e;
        for (int i = 0; i < ndim; i++)
            for (int j = 0; j < ndim; j++)
                e.ref(i, j) = 0.5 * (F(i, j) + F(j, i)) - (i == j);
        return setTrialStrain(e);
    }
}

template <int ndim, PlaneType type>
int ElasticIsotropic<ndim, type>::setTrialStrain(const MatrixSD<ndim, true> &e)
{
    Smat.vector = ddsdde * e.vector;
    return 0;
}


template<int ndim, PlaneType type>
auto
ElasticIsotropic<ndim, type>::getTangent()->MatrixSD<ne>
{
    return ddsdde;
}

template<int ndim, PlaneType type>
auto
ElasticIsotropic<ndim, type>::getInitialTangent()->MatrixSD<ne>
{
    return ddsdde;
}

template<int ndim, PlaneType type>
void
ElasticIsotropic<ndim,type>::Print(OPS_Stream& s, int flag) 
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Elastic Isotropic Material Model" << "\n";
        s << "\tE:  "   << E   << "\n";
        s << "\tv:  "   << nu  << "\n";
        s << "\trho:  " << rho << "\n";

    } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_ELEM_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"" << this->getClassType() << "\", ";
        s << "\"E\": "   << E   << ", ";
        s << "\"nu\": "  << nu  << ", ";
        s << "\"rho\": " << rho;
    s << "}";
    }
}
