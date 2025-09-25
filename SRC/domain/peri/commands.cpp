#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <tcl.h>
#include <Parsing.h>
#include <Logging.h>

#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <Logging.h>
#include <PeriDomain.h>
#include <PeriElement.h>
#include <PeriDomainBase.h>

#include <ElasticIsotropic.h>
#include <PeriParticle.h>
#include <NosbProj.h>


Tcl_CmdProc Tcl_Peri;

// Declare the function 'static' to convey that
// it is not expected to be used outside this file.
static int
Tcl_PeriInit(ClientData cd, Tcl_Interp *interp,
             int argc, const char **const argv)
{
    // If init has already been called, clean up (free) the previous
    // domain object.
    if (cd != nullptr)
    {
        PeriDomainBase *domain = static_cast<PeriDomainBase *>(cd);
        delete domain;
    }

    // Parse arguments and create a domain
    // Note that argv[0] holds the name of the command.

    // default values
    int ndim = 0;
    int totnode = 0;
    int maxfam = 0;
    char plane_type = 'x';
    PeriDomainBase *domain;

    int argi = 2;

    // if the size of argv is greater than 2, convert
    // the 3rd element of the array (ie argv[2]) to an integer
    // and store it in totnode
    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ndim) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[2] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Now do the same for the 4th argument (argv[3])
    if (argc > argi)
    {
        if (Tcl_GetInt(interp, argv[argi], &totnode) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[3] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Now do the same for the 5th argument (argv[4])
    if (argc > argi)
    {
        if (Tcl_GetInt(interp, argv[argi], &maxfam) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[4] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }
    // argv[5] is only required for 2D
    if (argc > argi)
    {
        plane_type = char(argv[argi][0]);
    }

    printf("Creating %d D domain with totnode=%d and maxfam=%d...\n", ndim, totnode, maxfam);
    // Allocate a new domain
    if (ndim == 2)
    {
        if (plane_type != 's' && plane_type != 'e')
        {
            printf("ERROR in peri init: if ndim == 2, correct plane_type ('s' or 'e') is required.\n");
            return TCL_ERROR;
        }
        else
        {
            domain = new PeriDomain<2>(totnode, maxfam);
            domain->plane_type = plane_type;
        }
    }
    else if (ndim == 3)
    {
        if (plane_type != 'x')
        {
            printf("EERROR in peri init: if ndim == 3, plane_type is not required.\n");
            return TCL_ERROR;
        }
        else
        {
            domain = new PeriDomain<3>(totnode, maxfam);
        }
    }
    else
    {
        printf("ERROR in peri init: ndim should be 2 or 3\n");
        return TCL_ERROR;
    }

    if (ndim == 3)
    {
        printf("3D dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }
    else if (ndim == 2 && plane_type == 's')
    {
        printf("2D plane stress dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }
    else if (ndim == 2 && plane_type == 'e')
    {
        printf("2D plane strain dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }

    // Store the pointer to the domain in the ClientData
    Tcl_CreateCommand(interp, "peri", Tcl_Peri, static_cast<ClientData>(domain), nullptr);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriSetNode(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                int argc, const char **const argv)
{
    int ind = 0, argi = 2;
    std::array<double, ndim> coord;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "node"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR)
        {
            printf("ERROR in peri node: Couldnt parse argv[2] as the index\n");
            return TCL_ERROR;
        }
        argi++;
    }

    for (int j = 0; j < ndim; j++)
    {
        if (argc > argi)
        {
            if (Tcl_GetDouble(interp, argv[argi], &coord[j]) == TCL_ERROR)
            {
                printf("ERROR in peri node: Couldnt parse argv[%d] as a double\n", argi);
                return TCL_ERROR;
            }
            argi++;
        }
        else
        {
            printf("ERROR in peri node: Not enough arguments\n");
            return -1; // **QUESTION: Is this the correct return value?**
        }
    }
    // Set the coordinates of the particle at index i
    domain->set_coord(ind, coord);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintNode(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                  int argc, const char **const argv)
{
    int ind = 0, argi = 3;

    if (argc == 3)
    {
        // print the coordinates of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: (", i);
            for (int j = 0; j < ndim; j++)
            {
                printf("%f", domain->pts[i].coord[j]);
                if (j < ndim - 1)
                {
                    printf(", ");
                }
            }
            printf(")\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "node"
        // Therefore, the first argument is argv[3]
        for (int i = argi; i < argc; i++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[i], &ind) == TCL_ERROR)
            {
                printf("ERROR in peri prin node: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the coordinates of the particle at index i
            printf("Node %d: (", ind);
            for (int j = 0; j < ndim; j++)
            {
                printf("%f", domain->pts[ind].coord[j]);
                if (j < ndim - 1)
                {
                    printf(", ");
                }
            }
            printf(")\n");
        }
    }
    else
    {
        printf("ERROR in peri prin node: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriCreateFam(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                  int argc, const char **const argv)
{
    int argi = 2;
    double delta = 0.0;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "fam"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &delta) == TCL_ERROR)
        {
            printf("ERROR in peri fam: Couldnt parse argv[2] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Create families for each particle
    domain->create_fam(delta);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriSetVols(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                int argc, const char **const argv)
{
    int ind = 0, argi = 2;
    double vol = 0.0;

    // Parse arguments and set the volume of a particle
    // Note that argv[0] == "peri", argv[1] == "vol"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR)
        {
            printf("ERROR in peri vol: Couldnt parse argv[2] as the index\n");
            return TCL_ERROR;
        }
        argi++;
    }

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[3] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &vol) == TCL_ERROR)
        {
            printf("ERROR in peri vol: Couldnt parse argv[3] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Set the volume of the particle at index i
    domain->set_vols(ind, vol);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintFam(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int ind = 0, argi = 3;

    if (argc == 3)
    {
        // print the number of families and their indices of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: %d neighbors\n", i, domain->pts[i].numfam);
            for (int j = 0; j < domain->pts[i].numfam; j++)
            {
                printf("%d ", domain->pts[i].nodefam[j]);
            }
            printf("\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "fam"
        // Therefore, the first argument is argv[3]
        for (int i = argi; i < argc; i++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[i], &ind) == TCL_ERROR)
            {
                printf("ERROR in peri prin fam: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the families of the particle at index i
            printf("Node %d: %d neighbors\n", ind, domain->pts[ind].numfam);
            for (int j = 0; j < domain->pts[ind].numfam; j++)
            {
                printf("%d ", domain->pts[ind].nodefam[j]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("ERROR in peri prin fam: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriCalcVols(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int argi = 2;
    double space = 0.0;

    // Parse arguments and set the volume of a particle
    // Note that argv[0] == "peri", argv[1] == "cvol"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &space) == TCL_ERROR)
        {
            printf("ERROR in peri cvol: Couldnt parse argv[2] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Calculate the volume of the horizons
    domain->calc_vols(space);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriElem(ClientData clientData, Tcl_Interp *interp,
             int argc, const char **const argv)
{
    PeriDomain<ndim> *pdomain = static_cast<PeriDomain<ndim> *>(clientData);
    PeriElement(1, pdomain);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintVol(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int argi = 3, idx = 0;

    if (argc == 3)
    {
        // print the volume of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: %f\n", i, domain->pts[i].vol_h);
            for (int ind = 0; ind < domain->pts[i].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[i].vol[ind]);
            }
            printf("\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "vol"
        // Therefore, the first argument is argv[3]
        for (int k = argi; k < argc; k++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[k], &idx) == TCL_ERROR)
            {
                printf("ERROR in peri prin vol: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the volume of the particle at index idx
            printf("Node %d: %f\n", idx, domain->pts[idx].vol_h);
            for (int ind = 0; ind < domain->pts[idx].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[idx].vol[ind]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("ERROR in peri prin vol: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintCorr(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                  int argc, const char **const argv)
{
    int argi = 3, idx = 0;

    if (argc == 3)
    {
        // print the volume of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d:\n", i);
            for (int ind = 0; ind < domain->pts[i].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[i].correction[ind]);
            }
            printf("\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "vol"
        // Therefore, the first argument is argv[3]
        for (int k = argi; k < argc; k++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[k], &idx) == TCL_ERROR)
            {
                printf("ERROR in peri prin vol: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the volume of the particle at index idx
            printf("Node %d:\n", idx);
            for (int ind = 0; ind < domain->pts[idx].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[idx].correction[ind]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("ERROR in peri prin vol: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriSetBoun(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                int argc, const char **const argv)
{
    int ndof = 0;
    char btype;
    std::array<double, 2 * ndim + 1> cond;

    // Parse arguments and set the boundary conditions
    // Note that argv[0] == "peri", argv[1] == "boun"
    // Therefore, the first argument is argv[2]
    // the syntax is: `peri boun xlow ylow zlow xhi yhi zhi val ndof btype` for 3D
    // and `peri boun xlow ylow xhi yhi val ndof btype` for 2D
    if (argc != 2 + 2 * ndim + 3)
    {
        printf("ERROR in peri boun: Incorrect number of arguments\n");
        return TCL_ERROR;
    }
    for (int i = 2; i < 2 + 2 * ndim + 1; i++)
    {
        if (Tcl_GetDouble(interp, argv[i], &cond[i - 2]) == TCL_ERROR)
        {
            printf("ERROR in peri boun: Couldnt parse argv[%d] as a double\n", i);
            return TCL_ERROR;
        }
    }
    if (Tcl_GetInt(interp, argv[2 + 2 * ndim + 1], &ndof) == TCL_ERROR)
    {
        printf("ERROR in peri boun: Couldnt parse argv[%d] as an integer\n", 2 + 2 * ndim + 1);
        return TCL_ERROR;
    }
    btype = argv[2 + 2 * ndim + 2][0];
    if (btype != 'f' && btype != 'd')
    {
        printf("ERROR in peri boun: btype should be 'f' or 'd'\n");
        return TCL_ERROR;
    }
    // Create families for specific NOSB type
    //   for (PeriParticle<ndim>& particle : domain.pts) {
    //     nodefam.emplace_back(&particle, domain, new ElasticIsotropic<ndim>(1, 29e3, 0.2, 0.0));
    //   }

    // Set the boundary conditions
    domain->set_bound(cond, ndof, btype);

    return TCL_OK;
}


template <int ndim>
static int
Tcl_PeriIntVerlet(PeriDomain<ndim> &domain, Tcl_Interp *interp, int argc, const char **const argv)
{
    // this is a temporary integrator for the PeriDomain, just for testing purposes
    int tt = 0;
    int part = 0;
    double dens = 0.0;
    double dt_verlet = 0.0;
    // create a 2d vector with the size of (totnode, ndim)
    // and initialize it with zeros
    int totnode = domain.totnode;
    Matrix acc_old(totnode, ndim);
    Matrix acc_new(totnode, ndim);

    acc_old.Zero();
    acc_new.Zero();

    // sparse values from argv
    Tcl_GetInt(interp, argv[2], &tt);
    tt += 1;
    Tcl_GetInt(interp, argv[3], &part);
    Tcl_GetDouble(interp, argv[4], &dens);
    Tcl_GetDouble(interp, argv[5], &dt_verlet);
    // ---------FOR DEBUGGING---------
    // printf("tt=%d, part=%d, dens=%e, dt_verlet=%e\n", tt, part, dens, dt_verlet);
    // -------------------------------

    if (part == 0)
    {
        // x_{n+1} = x_{n} + v_{n} * dt + 0.5 * a_{n} * dt^2
        // compute a_{n} by (pforce + bforce) / dens
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                acc_old(i, j) = (domain.pts[i].pforce[j] + domain.pts[i].bforce[j]) / dens;
        // ---------FOR DEBUGGING---------
        // printf("acc_old[0, 0]=%e\n", acc_old(0, 0));
        // -------------------------------
        // update disp
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                if (domain.pts[i].is_disp_bound[j] == 1)
                {
                    domain.pts[i].disp[j] = domain.pts[i].bdisp[j];
                }
                else
                {
                    domain.pts[i].disp[j] = domain.pts[i].disp[j] + domain.pts[i].vel[j] * dt_verlet + 0.5 * acc_old(i, j) * dt_verlet * dt_verlet;
                }
        // save acc_old[i, j] to domain.pts[i].acc[j]
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                domain.pts[i].acc[j] = acc_old(i, j);
        // ---------FOR DEBUGGING---------
        // for (int i = 0; i < 10; i++){
        // 	printf("Node %d: disp=(%7.2e, %7.2e)\n", i, domain.pts[i].disp[0], domain.pts[i].disp[1]);
        // }
        // -------------------------------
    }
    else if (part == 1)
    {
        // v_{n+1} = v_{n} + 0.5 * (a_{n} + a_{n+1}) * dt
        // compute acc
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                acc_new(i, j) = (domain.pts[i].pforce[j] + domain.pts[i].bforce[j]) / dens;
        // update vel
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                acc_old(i, j) = domain.pts[i].acc[j];
        for (int i = 0; i < totnode; i++)
            for (int j = 0; j < ndim; j++)
                if (domain.pts[i].is_disp_bound[j] == 0)
                    domain.pts[i].vel[j] = domain.pts[i].vel[j] + 0.5 * (acc_old(i, j) + acc_new(i, j)) * dt_verlet;
    }

    return TCL_OK;
}


template <int ndim, int maxfam>
static int
Tcl_PeriForm(PeriDomain<ndim> &domain, Tcl_Interp *interp)
{

    std::vector<NosbProj<ndim, maxfam>> nosb;

    // Create families for specific NOSB type
    if constexpr (ndim == 3)
    {
        for (PeriParticle<ndim> &particle : domain.pts)
        {
            nosb.emplace_back(&particle, domain, new ElasticIsotropic<ndim>(1, 38.4e3, 0.2, 0.0));
        }
    }
    else if (ndim == 2 && domain.plane_type == 's')
    {
        for (PeriParticle<ndim> &particle : domain.pts)
        {
            nosb.emplace_back(&particle, domain, new ElasticIsotropic<ndim, PlaneType::Stress>(1, 38.4e3, 0.2, 0.0));
        }
    }
    else if (ndim == 2 && domain.plane_type == 'e')
    {
        for (PeriParticle<ndim> &particle : domain.pts)
        {
            nosb.emplace_back(&particle, domain, new ElasticIsotropic<ndim, PlaneType::Strain>(1, 38.4e3, 0.2, 0.0));
        }
    }

    // -------- FOR DEBUGGING --------
    // printf("Successfully create families for specific NOSB type\n");
    // -------------------------------

    // Initialize shape tensor
    for (NosbProj<ndim, maxfam> &nosb_i : nosb)
        nosb_i.init_shape();
    // -------- FOR DEBUGGING --------
    // printf("Successfully initialize shape tensor\n");
    // -------------------------------

    // -------- FOR DEBUGGING --------
    // print first 10 shape tensors for debugging
    // for (int i = 0; i < 10; i++)
    // {
    // 	printf("Shape tensor %d:\n", i);
    // 	for (int j = 0; j < ndim; j++)
    // 	{
    // 		for (int k = 0; k < ndim; k++)
    // 		{
    // 			printf("%7.2e ", nosb[i].Kinv(j, k));
    // 		}
    // 		printf("\n");
    // 	}
    // }
    // -------------------------------

    // Form deformation gradients for trial
    for (NosbProj<ndim, maxfam> &nosb_i : nosb)
        nosb_i.form_trial();

    // -------- FOR DEBUGGING --------
    // for (int i = 0; i < 1; i++)
    // {
    // 	printf("Form trial %d\n", i);
    // 	nosb[i].form_trial();
    // }
    // -------------------------------

    // // Form force
    for (NosbProj<ndim, maxfam> &nosb_i : nosb)
    {

    	MatrixND<ndim, ndim> Q = nosb_i.sum_PKinv();

    	for (int j = 0; j < nosb_i.numfam; j++)
    	{
    		const VectorND<ndim> T_j = nosb_i.bond_force(j, Q);

    		nosb_i.center->pforce += T_j;
    		nosb_i.neigh[j]->pforce -= T_j;
    	}
    }

    // -------- FOR DEBUGGING --------

    for (int i = 0; i < 10; i++) {
        // MatrixND<ndim, ndim> Q = nosb[i].sum_PKinv();
        // for (int j = 0; j < ndim; j++)
        //     for (int k = 0; k < ndim; k++)
        //         printf("%7.2e ", Q(j, k));
        // printf("\n");
        // printf("Particle %d: pforce=(%7.2e, %7.2e)\n", i, nosb[i].center->pforce[0], nosb[i].center->pforce[1]);
    }
    

    // //
    // // TODO: Update disp
    // //

    // Create and return force vector
    Tcl_Obj *list = Tcl_NewListObj(nosb.size() * ndim, nullptr);
    for (PeriParticle<ndim> &particle : domain.pts)
    	for (int j = 0; j < ndim; j++)
    		Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(particle.pforce[j]));

    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int Tcl_Peri(ClientData cd, Tcl_Interp *interp,
             int argc, const char **const argv)
{
    // TODO Ensure argv[1] exists
    if (strcmp(argv[1], "init") == 0)
    {
        return Tcl_PeriInit(cd, interp, argc, argv);
    }

    if (cd == nullptr)
    {
        printf("ERROR: ClientData is null.\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }

    PeriDomainBase *domain_base = static_cast<PeriDomainBase *>(cd);

    if (strcmp(argv[1], "node") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriSetNode(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriSetNode(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "fam") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriCreateFam(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriCreateFam(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "svol") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriSetVols(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriSetVols(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "cvol") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriCalcVols(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriCalcVols(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "suco") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            domain->calc_surf_correction();
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            domain->calc_surf_correction();
        }
    }

    if (strcmp(argv[1], "boun") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriSetBoun(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriSetBoun(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "step") == 0)
    {
        // this is a temporary command, will delete it later
        // syntax: peri step tt part dens dt_verlet
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriIntVerlet(*domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriIntVerlet(*domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "form") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            int maxfam = domain->maxfam;
            if (maxfam < 32)
                return Tcl_PeriForm<3,32>(*domain, interp);
            if (maxfam < 64)
                return Tcl_PeriForm<3,64>(*domain, interp);
            if (maxfam < 1024)
                return Tcl_PeriForm<3,1024>(*domain, interp);
            return TCL_ERROR;
        }
        else if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);

            int maxfam = domain->maxfam;
            if (maxfam < 32)
                return Tcl_PeriForm<2,32>(*domain, interp);
            if (maxfam < 64)
                return Tcl_PeriForm<2,64>(*domain, interp);
            if (maxfam < 1024)
                return Tcl_PeriForm<2,1024>(*domain, interp);
            return TCL_ERROR;
        }
    }

    if (strcmp(argv[1], "prin") == 0)
    {

        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            if (strcmp(argv[2], "node") == 0)
            {
                return Tcl_PeriPrintNode(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "fam") == 0)
            {
                return Tcl_PeriPrintFam(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "vol") == 0)
            {
                return Tcl_PeriPrintVol(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "corr") == 0)
            {
                return Tcl_PeriPrintCorr(domain, interp, argc, argv);
            }
            // else if (strcmp(argv[2], "boun") == 0)
            // {
            // 	return Tcl_PeriPrintBoun(domain, interp, argc, argv);
            // }
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            if (strcmp(argv[2], "node") == 0)
            {
                return Tcl_PeriPrintNode(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "fam") == 0)
            {
                return Tcl_PeriPrintFam(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "vol") == 0)
            {
                return Tcl_PeriPrintVol(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "corr") == 0)
            {
                return Tcl_PeriPrintCorr(domain, interp, argc, argv);
            }
            // else if (strcmp(argv[2], "boun") == 0)
            // {
            // 	return Tcl_PeriPrintBoun(domain, interp, argc, argv);
            // }
        }
    }
    return TCL_OK;
}
#pragma clang diagnostic pop
