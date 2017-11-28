/**
 * \file main.cpp
 * \brief Main file (exec) for the computation of the dynamics about the neighborhood of the points \f$ EML_{1,2} \f$ and \f$ SEML_{1,2} \f$.
 *
 * The following models of the Sun-Earth-Moon system are available:
 *   - Quasi-Bicircular Four-Body Problem (QBCP).
 *   - Coupled Circular Restricted Three-Body Problem (CRTBP).
 * The following models of the Sun-Earth-Moon system may be available in the future:
 *   - Bircicular Four-Body Problem (BCP) (work in progress. the real issue being: no clear dynamical equivalent for \f$ EML_{2} \f$.
 *   - Coupled Elliptic Restricted Three-Body Problem (CRTBP) (far from being completed).
 * \author BLB.
 */

//----------------------------------------------------------------------------------------
// Comment & test
//----------------------------------------------------------------------------------------
// STILL a memory leak...
//--------------------
// oftsh.tpp
// ofts.tpp
//
//
// Not currently used:
//--------------------
// ots.tpp
// otsh.tpp
// multimin.c
// multimin_test.cpp
// trajectory.cpp
//--------------------
//
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// TODO:
// CHECKUP COMPLET DES DATA.
// TEST EACH CM/CU/CS
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
// Include
//----------------------------------------------------------------------------------------
//Parallel computing
#include <omp.h>
//Custom
#include "Config.h"
#include "ofs.h"
#include "ofts.h"
#include "poincare.h"
#include "trajectory.h"
#include "pmt.h"
#include "pmode.h"
//Tests
#include "ofs_test.h"
#include "oftsh_test.h"
#include "ofts_test.h"
//C routines and files
extern "C" {
    #include "gnuplot_i.h"
    void fortfunc_(int *ii, float *ff);
    void saxpy_(int *n, double *alpha, double *x, double *y);
    void matvec_(int *n, double **A, double *x, double *y, int *lda);
}

//----------------------------------------------------------------------------------------
// Namespaces
//----------------------------------------------------------------------------------------
using namespace std;


//----------------------------------------------------------------------------------------
// Main
//----------------------------------------------------------------------------------------
/**
 *  \fn int main()
 *  \brief Main routine. It makes of several arguments via argv (see code for details).
 */
int main(int argc, char *argv[])
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "                    OOFTDA                         " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << argc << " arguments have been passed to OOFTDA"       << endl;

    //====================================================================================
    // Initialization of Manip object, which allows for algebraic manipulation
    // of polynomials (exponents in reverse lexicographic order, number of coefficient in
    // homogeneous polynomials...).
    // Note that, for now, the default order of the OFS and OFTS objects are set at
    // compilation (OFS_ORDER_0OFTS_ORDER_0).
    //====================================================================================
    int OFS_ORDER_0  = 30;  //The default OFS_ORDER (Fourier series) is 30
    int OFTS_ORDER_0 = 30;  //The default OFTS_ORDER (Taylor series) is 30
    int nr = max(2*OFS_ORDER_0, OFTS_ORDER_0);
    int nv = Csts::NV;
    Manip::init(nv, nr);


    //====================================================================================
    // Get the user-defined arguments
    //====================================================================================
    // The variable index contains the index of the current argument argv[index]
    // of the main routine. The first index is 1 because index 0 is the application path.
    int index    = 1;

    //------------------------------------------------------------------------------------
    // Get the current orders
    //------------------------------------------------------------------------------------
    OFTS_ORDER   = atoi(argv[index++]); //Order of the Taylor series
    OFS_ORDER    = atoi(argv[index++]); //Order of the Fourier series
    REDUCED_NV   = atoi(argv[index++]); //Number of reduced variables (4, 5 or 6).

    //To check consistency with the desired values
    cout << "Current orders:" << endl;
    cout << "OFTS_ORDER = "   << OFTS_ORDER << endl;
    cout << "OFS_ORDER  = "   << OFS_ORDER  << endl;
    cout << "REDUCED_NV = "   << REDUCED_NV << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // Retrieving the parameters, in this order:
    // 1. compType: Type of computation (QBTBP, NFO2, PM...)
    // 2. model: Type of model (QBCP, CRTBP)
    // 3. coord_sys: Default coordinate system (EM, SEM)
    // 4. dli: Default libration point
    // 5. param_style: Parameterization (PM) style (only used in some computations)
    // 6. mtype: Type of manifold (center, center-stable, center-unstable)
    // 7. storage: Boolean for storage (in txt/bin files) of the results
    //------------------------------------------------------------------------------------
    int compType    = atoi(argv[index++]);
    int model       = atoi(argv[index++]);
    int dli         = atoi(argv[index++]);
    int param_style = atoi(argv[index++]);
    int mtype       = atoi(argv[index++]);
    int storage     = atoi(argv[index++]);

    //Check
    cout << "Current parameters: " << endl;
    cout << "--------------------" << endl;
    cout << "compType  = "   << compType << endl;
    cout << "model     = "   << model << endl;
    cout << "li        = "   << dli << endl;
    cout << "pmstyle   = "   << param_style << endl;
    cout << "mtype     = "   << mtype << endl;
    cout << "storage   = "   << storage << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------------------------------------------------
    // Set the global variable MODEL_TYPE
    //------------------------------------------------------------------------------------
    MODEL_TYPE = model;

    //====================================================================================
    // Initialization of the environnement (the Four-Body Problem at hand).
    // Mandatory to perform any computation, except qbtbp(int)
    //====================================================================================
    init_env(model, dli, mtype, param_style);


    //====================================================================================
    // Master switch on the type of computation required by the user
    // QBTBP   = 0
    // NFO2    = 1
    // PM      = 2
    // PM_TEST = 3
    // PMAP    = 4
    //====================================================================================
    switch(compType)
    {

    //------------------------------------------------------------------------------------
    // Unitary tests of routines for the Fourier-Taylor algebra
    //------------------------------------------------------------------------------------
    case Csts::FT_TEST:
    {
        //Test of the Ofs class (Fourier series)
        ofs_test();
        //Test of the Oftsh class (Fourier-Taylor homogeneous polynomials)
        oftsh_test();
        //Test of the Ofts class (Fourier-Taylor series)
        ofts_test();
        break;
    }


    //------------------------------------------------------------------------------------
    // Sun-Earth-Moon three-body motion resolution up to OFS_ORDER.
    //------------------------------------------------------------------------------------
    case Csts::QBTBP:
    {
        switch(model)
        {
            //----------------------------------------------------------------------------
            // QBTBP & QBCP resolution up to OFS_ORDER.
            // If storage==true, results are stored in the folders
            // data/qbtbp and data/VF/QBCP.
            //----------------------------------------------------------------------------
        case Csts::QBCP:
            qbtbp_and_qbcp(true, storage);
            break;

            //----------------------------------------------------------------------------
            // BCP resolution in OFS format (in order to match QBCP format).
            // Results are stored in the folder data/VF/BCP.
            //----------------------------------------------------------------------------
        case Csts::BCP:
            bcp();
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Compute the complete change of coordinates to:
    // - Get rid of order 1,
    // - Get a normal form for the order 2 of the Hamiltonian of the QBCP.
    // If storage==true, results are stored in the folder data/COC
    //------------------------------------------------------------------------------------
    case Csts::NFO2:
    {
        switch(model)
        {
        case Csts::QBCP:
            nfo2(SEML, storage);
            break;
        default:
            cout << "Error: nfo2 routine only implemented for the QBCP." << endl;
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Parameterization method
    //------------------------------------------------------------------------------------
    case Csts::PM:
    {
        //--------------------------------------------------------------------------------
        //threshold for small divisors = 1e-2 (could be propagated to the user)
        //--------------------------------------------------------------------------------
        double small_div_threshold = 1e-2;

        //--------------------------------------------------------------------------------
        //param method
        //--------------------------------------------------------------------------------
        pmt(SEML.cs.man_type, param_style, small_div_threshold, storage);
        break;
    }

    //------------------------------------------------------------------------------------
    // Test of the parameterization method
    //------------------------------------------------------------------------------------
    case Csts::PM_TEST:
    {
        //--------------------------------------------------------------------------------
        // Additionnal parameters in bash
        //--------------------------------------------------------------------------------
        // 1. Initial conditions
        double si[REDUCED_NV];
        for(int i = 0; i < REDUCED_NV; i++) si[i] = atof(argv[index++]);

        //2. Array of ofts_orders to test
        int n_ofts_order = atoi(argv[index++]);
        int v_ofts_order[n_ofts_order];
        for(int i = 0; i<n_ofts_order; i++)
        {
            v_ofts_order[i] = atoi(argv[index++]);
            if(v_ofts_order[i] > OFTS_ORDER)
            {
                cout << "Warning: a required order is greater than OFTS_ORDER. " << endl;
                cout << "The computation will be cut at OFTS_ORDER." << endl;
            }
        }


        //3. Array of ofs_orders to test
        int n_ofs_order = atoi(argv[index++]);
        int v_ofs_order[n_ofs_order];
        for(int i = 0; i<n_ofs_order; i++)
        {
            v_ofs_order[i] = atoi(argv[index++]);
            if(v_ofs_order[i] > OFS_ORDER)
            {
                cout << "Warning: a required order is greater than OFS_ORDER. " << endl;
                cout << "The computation will be cut at OFS_ORDER." << endl;
            }
        }

        //4. Time interval
        double tvec[2];
        tvec[0] = atof(argv[index++])*SEML.us.T;
        tvec[1] = atof(argv[index++])*SEML.us.T;


        //--------------------------------------------------------------------------------
        // Initialization of the invariant manifold
        //--------------------------------------------------------------------------------
        init_inv_man(SEML);
        init_coc(SEML);

        //--------------------------------------------------------------------------------
        // Test
        //--------------------------------------------------------------------------------
        pm_error_vs_orders_test(n_ofts_order, v_ofts_order, n_ofs_order, v_ofs_order, si, tvec);

        break;
    }

    //------------------------------------------------------------------------------------
    // Poincare, stroboscopic and error maps
    //------------------------------------------------------------------------------------
    case Csts::COMPMAP:
    {
        //--------------------------------------------------------------------------------
        // Additionnal parameters in bash
        //--------------------------------------------------------------------------------
        //1. Array of orders to test
        int nkm = atoi(argv[index++]);
        int km[nkm];
        for(int i = 0; i<nkm; i++)
        {
            km[i] = atoi(argv[index++]);
            if(km[i] > OFTS_ORDER)
            {
                cout << "Warning: a required order is greater than OFTS_ORDER: return." << endl;
                return -1;
            }
        }

        //2. Array of OFS order.
        int nkofs = atoi(argv[index++]);
        int kofs[nkofs];
        for(int i = 0; i<nkofs; i++)
        {
            kofs[i] = atoi(argv[index++]);
            if(kofs[i] > OFS_ORDER)
            {
                cout << "Warning: a required order is greater than OFS_ORDER: return." << endl;
                return -1;
            }
        }

        //3. PMAP parameters
        int  argpmapc  = atoi(argv[index++]);
        cout << "Number of PMAP parameters & settings is " << argpmapc << endl;

        //--------------------------------------------------------------------------------
        // Initialization of the invariant manifold
        //--------------------------------------------------------------------------------
        init_inv_man(SEML);
        init_coc(SEML);


        //--------------------------------------------------------------------------------
        // Pmap init
        // Bash line for parallel computation: export OMP_SET_NUM_THREADS=5
        //--------------------------------------------------------------------------------
        Pmap pmap;
        if(MODEL_TYPE == Csts::CRTBP)
        {
            //Quite stable parameters
            //------------------------
            pmap.pabs       =  1e-14;
            pmap.prel       =  1e-14;
            pmap.proot      =  1e-13;
            pmap.threshold  =  1e-6;
            pmap.tt         =  M_PI;
            pmap.T          =  SEML.us.T;;
            pmap.vdim       =  2;
            pmap.maxRad     =  5.0;
        }
        else
        {
            //Quite stable parameters
            //------------------------
            pmap.pabs       =  1e-14;
            pmap.prel       =  1e-14;
            pmap.proot      =  1e-10;
            pmap.threshold  =  1e-6;
            pmap.tt         =  0.5*SEML.us.T;
            pmap.T          =  SEML.us.T;
            pmap.vdim       =  2;
            pmap.maxRad     =  5.0;
        }


        //--------------------------------------------------------------------------------
        // Parameters that may evolve
        //
        //  - projFreq is the frequency of projection inside the map
        //  - type is the type of map: either PMAP, TMAP, EMAP or IMAP (or IMAPPLANAR)
        //  - tf is the final time of integration (useful only for PMAP/TMAP).
        //  - isGS: force some simplifications in the evaluation
        //    of the manifolds. The parameterization style of the manifold at hand must be
        //    the graph style (Csts::GRAPH) for isGS to be active.
        //  - order: the effective order of the Taylor expansions on the map.
        //  - ofs_order: the effective order of the Fourier expansions on the map.
        //  - max_events: the maximum number of events on the map. Typically,
        //    the number of points on the xy-plane for the PMAP.
        //  - t0 the initial time for all trajectories on the map.
        //    t0 is given as a fraction of the period T of the SEM system.
        //    T is computed in the units (EM or SEM) selected earlier by the user.
        //
        //--------------------------------------------------------------------------------
        pmap.projFreq   =  atoi(argv[index++]);
        pmap.type       =  atoi(argv[index++]);
        pmap.tf         =  atof(argv[index++]); //needs to be a big number of Pmap, not for Tmap

        //isGS and graph style?
        //helps enhance computation
        int isGS     =  atoi(argv[index++]);
        if(isGS == 1 && SEML.param_style == Csts::GRAPH) pmap.isGS = 1;
        else pmap.isGS = 0;

        pmap.order      =  atoi(argv[index++]);
        pmap.ofs_order  =  atoi(argv[index++]);
        pmap.max_events =  atoi(argv[index++]);

        //Specific case of the initial time, that
        //may be initialized as a root of a coefficient
        //of the COC matrix
        pmap.t0         =  atof(argv[index++])*SEML.us.T;
        if(pmap.t0 == -1) pmap.t0 = pij(2, 5); //root of p36 is input is -1

        pmap.dHv        =  atof(argv[index++]);
        pmap.gsize      =  atoi(argv[index++]);
        pmap.gmin       =  atof(argv[index++]);
        pmap.gmax       =  atof(argv[index++]);

        cout << "pmap.type       = "  << pmap.type << endl;
        cout << "pmap.tf         = "  << pmap.tf << endl;
        cout << "pmap.projFreq   = "  << pmap.projFreq << endl;
        cout << "pmap.isGS       = "  << pmap.isGS << endl;
        cout << "pmap.order      = "  << pmap.order << endl;
        cout << "pmap.ofs_order  = "  << pmap.ofs_order  << endl;
        cout << "pmap.max_events = "  << pmap.max_events << endl;
        cout << "pmap.t0         = "  << pmap.t0  << endl;
        cout << "pmap.dHv        = "  << pmap.dHv   << endl;
        cout << "pmap.gsize      = "  << pmap.gsize << endl;
        cout << "pmap.gmin       = "  << pmap.gmin << endl;
        cout << "pmap.gmax       = "  << pmap.gmax << endl;
        cout << "---------------------------------------------------" << endl;

        //--------------------------------------------------------------------------------
        // Additionnal map settings
        //--------------------------------------------------------------------------------
        int append = atoi(argv[index++]);
        int isPlot = atoi(argv[index++]);
        int isPar  = atoi(argv[index++]);
        int method = atoi(argv[index++]);

        cout << "append          = "  << append << endl;
        cout << "isPlot          = "  << isPlot << endl;
        cout << "isPar           = "  << isPar << endl;
        cout << "method          = "  << pmap.type << endl;

        //--------------------------------------------------------------------------------
        // openMP settings
        //--------------------------------------------------------------------------------
        int num_threads = atoi(argv[index++]);
        omp_set_num_threads(num_threads);
        cout << "num_threads     = "  << num_threads << endl;

        //--------------------------------------------------------------------------------
        // Pmap computation
        //--------------------------------------------------------------------------------
        switch(pmap.type)
        {
        case PMAP:
        {
            pmap.tt =  1.0/pmap.projFreq*SEML.us.T; //for pmaps, pmap.tt is the period of projection
            pmap_build(pmap, append, method, isPar);
            break;
        }

        case TMAP:
        {
            tmap_build(pmap, append, method, isPlot, isPar);
            break;
        }

        //--------------------------------------------------------------------------------
        // Precision map
        // for various orders
        //
        //  €€ TODO: add a parameter to take into account
        //   the different possibilities for error maps
        //
        //--------------------------------------------------------------------------------
        case EMAP:
        {
            for(int i = 0; i < nkm; i++)
            {
                pmap.order = km[i];
                pmap_precision(pmap, append, isPar);
            }
            break;
        }


        //--------------------------------------------------------------------------------
        // Invariance error map
        // for various orders
        //
        //  €€ TODO: add a parameter to take into account
        //   the different possibilities for invariance error maps
        //
        //--------------------------------------------------------------------------------
        case IMAP:
        {
            pmap_invariance_error_random(pmap, append, isPar, pmap.dHv, km, nkm);
            break;
        }

        case IMAPPLANAR:
        {
            pmap_invariance_error_random_planar(pmap, append, isPar, pmap.dHv, km, nkm);
            break;
        }


        //--------------------------------------------------------------------------------
        // Energy (hamiltonian) map
        //
        //  €€ TODO: add a parameter to take into account
        //   the different possibilities for invariance error maps
        //
        //--------------------------------------------------------------------------------
        case HMAP:
        {
            //pmap_energy(pmap, append, isPar, pmap.dHv);
            pmap_energy_3d(pmap, append, isPar, pmap.dHv);
            break;
        }

        }
        break; //and of case 4: pmaps
    }

    //------------------------------------------------------------------------------------
    // COCs
    //------------------------------------------------------------------------------------
    case Csts::COC:
    {
        //-----------------------------
        // EM <--> IN <--> SEM
        //-----------------------------
        //GOOD:
        //------------------
        //dynTest_SEMtoEM();
        //dynTest_EMtoSEM();

        //OBSOLETE: €€TODO
        //------------------
        //dynTest_SEtoEM();
        //NCvsEM();
        //QBTBP_IN();
        //SE_DYN_NO_MOON();
        //------------------

        //-----------------------------
        // COC of nfo2:
        //-----------------------------
        //testCOC();
        //testDotCOC();
        //testIntCOC();
        //oftsh_test();

        break;
    }

    //------------------------------------------------------------------------------------
    // Compute the dynamical equivalent of the libration point, and other resonant orbits
    //------------------------------------------------------------------------------------
    case Csts::DYNEQ:
    {
        //-----------------------------
        // Direction computation of the lpdyneq. Works for:
        // - EML1, EML2 in QBCP
        // - SEML1, SEML2 in QBCP
        //-----------------------------
        compute_dyn_eq_lib_point(SEML, storage);

        //-----------------------------
        // Continuation procedures for the computation of the lpdyneq
        // Works for:
        // - EML1 in BCP.
        // - Seems to work for SEML1,2 of the BCP,
        //   but equations of motion need to be checked.
        //-----------------------------
        //continuation_dyn_eq_lib_point(SEML, Csts::CRTBP, Csts::BCP);

        //-----------------------------
        // Continuation procedures for the computation of other resonant orbits
        // Work for:
        // - A 2T/5-resonant orbit at EML1
        // - A T/2-resonant orbit at EML2
        //-----------------------------
        //continuation_res_orbit(SEML, Csts::CRTBP, Csts::QBCP);

        break;
    }


    }


    //------------------------------------------------------------------------------------
    // Note: their a double free or corruption from an unknown source.
    // For now, we do not quit properly!
    //------------------------------------------------------------------------------------
    char ch;
    printf("Press ENTER to finish the program\n");
    scanf("%c",&ch);


    return 0;


    //------------------------------------------------------------------------------------
    // Tests
    //------------------------------------------------------------------------------------
    //testCOC();
    //testDotCOC();
    //testIntCOC();
}



