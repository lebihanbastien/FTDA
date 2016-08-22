/**
 * \file main.cpp
 * \brief Main file for the computation of the neighborhood of the points \f$L_{1,2} \f$ of the Quasi-Bicircular Four-Body Problem (QBCP).
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

//-----------------------------------------------------
// Comment & test
//-----------------------------------------------------
//
//
// STILL a memory leak...
//--------------------
//
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
//-----------------------------------------------------

//-----------------------------------------------------
// TODO:
// CHECKUP COMPLET DES DATA.
// TEST EACH CM/CU/CS
// FIND A WAY TO ADD THE INFO ON THE CORRESPONDING FOLDERS (OTFS_ORDER, etc).
//-----------------------------------------------------



//-----------------------------------------------------
// Include
//-----------------------------------------------------
#include <omp.h>
//Custom
#include "ofs.h"
#include "ofts.h"
#include "poincare.h"
#include "trajectory.h"
#include "ertbp.h"
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

//-----------------------------------------------------
// Namespaces
//-----------------------------------------------------
using namespace std;


// Choice of coefficient type
// for Fourier-Taylor series
//---------------------------
typedef Ofs<cdouble> Ofsc;
typedef Ots<cdouble> Otsc;
typedef Ofts< Ofs<cdouble> > Oftsc;


//-----------------------------------------------------
// Inner routines
//-----------------------------------------------------
/**
    \brief From *char to int using istringstream
**/
double toDigit(char *c)
{
    int res;
    istringstream ss;
    ss.str(c);
    ss >> res;
    return res;
}

//-----------------------------------------------------
// Main
//-----------------------------------------------------
/**
 *  \fn int main()
 *  \brief Main routine.
 */
int main(int argc, char *argv[])
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "                    OOFTDA                         " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;

    //-------------------------------------------------------------------------------------------------------
    // Initialization of static FTDA objects. Mandatory.
    // Note that, for now, the default order of the OFS and OFTS objects are set at compilation
    //-------------------------------------------------------------------------------------------------------
    int OFS_ORDER_0  = 30;
    int OFTS_ORDER_0 = 30;
    int nr = max(2*OFS_ORDER_0, OFTS_ORDER_0);
    int nv = NV;
    FTDA::init(nv, nr);

    //-------------------------------------------------------------------------------------------------------
    // Get the current orders
    //-------------------------------------------------------------------------------------------------------
    int index    = 1; //first index is 1 because index 0 is the application path
    OFTS_ORDER   = toDigit(argv[index++]);
    OFS_ORDER    = toDigit(argv[index++]);
    REDUCED_NV   = toDigit(argv[index++]);
    OTS_ORDER    = 10;

    cout << "OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "OFS_ORDER  = " << OFS_ORDER  << endl;
    cout << "OTS_ORDER  = " << OTS_ORDER  << endl;
    cout << "REDUCED_NV = " << REDUCED_NV << endl;

    //-------------------------------------------------------------------------------------------------------
    // Retrieving the parameters, in this order:
    // 1. Type of computation (QBTBP, NFO2, PM...)
    // 2. Type of model (QBCP, RTBP...)
    // 3. Default coordinate system (EM, SEM...)
    // 4. Boolean for normalization (or not) of the equations of motion
    // 5. Default libration point (EM system)
    // 6. Default libration point (SEM system)
    // 7. PM style (only used in some computations)
    // 8. boolean for storage (in txt files) of the results, for ALL computations !
    //-------------------------------------------------------------------------------------------------------
    int compType  = toDigit(argv[index++]);
    int model     = toDigit(argv[index++]);
    int dcs       = toDigit(argv[index++]);
    int isNorm    = toDigit(argv[index++]);
    int li_EM     = toDigit(argv[index++]);
    int li_SEM    = toDigit(argv[index++]);
    int pms       = toDigit(argv[index++]);
    int mType_EM  = toDigit(argv[index++]);
    int mType_SEM = toDigit(argv[index++]);
    int storage   = toDigit(argv[index++]);


    //Check
    cout << "Current parameters:" << endl;
    cout << "--------------------" << endl;
    cout << "compType  = "   << compType << endl;
    cout << "model     = "   << model << endl;
    cout << "dcs       = "   << dcs << endl;
    cout << "isNorm    = "   << isNorm << endl;
    cout << "li_EM     = "   << li_EM << endl;
    cout << "li_SEM    = "   << li_SEM << endl;
    cout << "pms       = "   << pms << endl;
    cout << "mType_EM  = "   << mType_EM << endl;
    cout << "mType_SEM = "   << mType_SEM << endl;
    cout << "storage   = "   << storage << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------
    // Set the global variable MODEL_TYPE
    //------------------------------------------
    MODEL_TYPE = model;

    //------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, dcs, pms, mType_EM, mType_SEM);

    //------------------------------------------
    // Which default CS?
    //------------------------------------------
    changeDCS(SEML, dcs);

    //-------------------------------------------------------------------------------------------------------
    // Master switch on the type of computation required by the user
    // QBTBP   = 0
    // NFO2    = 1
    // PM      = 2
    // PM_TEST = 3
    // PMAP    = 4
    //-------------------------------------------------------------------------------------------------------
    cout << "T = " << SEML.us_em.T << endl;
    cout << "T_sem = " << SEML.us_sem.T << endl;
    switch(compType)
    {
    //-------------------------------------
    // QBTBP resolution up to OFS_ORDER.
    //-------------------------------------
    case 0:
    {
        switch(model)
        {
            //-------------------------------------
            // QBTBP resolution up to OFS_ORDER.
            // if the testing argument is true,
            // the results are tested on
            // one full period
            //-------------------------------------
        case M_QBCP:
            qbtbp(li_EM, li_SEM, true, dcs);
            break;

            //-------------------------------------
            // BCP resolution in OFS format.
            // Storage in txt file.
            //-------------------------------------
        case M_BCP:
            bcp(li_EM, li_SEM, dcs);
            break;
            //-------------------------------------
            // ERTBP resolution in OFS format.
            // Storage in txt file.
            //-------------------------------------
        case M_ERTBP:
            ertbp(li_EM, li_SEM, dcs);
            break;

        }
        break;
    }

    //------------------------------------------------
    // Compute the complete change of coordinates to:
    // - Get rid of order 1
    // - Get a normal form for the order 2
    // of the Hamiltonian of the QBCP
    // Results are stored in the folder data/COC
    //------------------------------------------------
    case 1:
    {
        switch(model)
        {
        case M_QBCP:
        case M_ERTBP:
            nfo2(SEML, storage);
            break;
        case M_BCP:
//            nfo2_QBP(SEML, storage);
            continuation(SEML, storage);
            break;
        }
        break;
    }


    //------------------------------------------------
    // Parameterization method
    //------------------------------------------------
    case 2:
    {
        pmt(0, storage, pms, SEML.cs.manType);
        break;
    }

    //------------------------------------------------
    // Test of the parameterization method
    //------------------------------------------------
    case 3:
    {
        //------------------------------------------
        // Additionnal parameters in bash
        //------------------------------------------
        // 1. Initial conditions
        double si[REDUCED_NV];
        for(int i = 0; i < REDUCED_NV; i++) si[i] = atof(argv[index++]);

        //2. Array of orders to test
        int nkm = toDigit(argv[index]);
        int km[nkm];
        for(int i = 0; i<nkm; i++)
        {
            km[i] = toDigit(argv[index+1+i]);
            if(km[i] > OFTS_ORDER)
            {
                cout << "Warning: a required order is greater than OFTS_ORDER: the computation will be stopped before this value." << endl;
                nkm = i;
            }
        }

        //------------------------------------------
        // Initialisation of the central manifold
        //------------------------------------------
        initCM(SEML);
        initCOC(SEML);

        //------------------------------------------
        // Test
        //------------------------------------------
        //pmEOvsOrderTest(nkm, km, si);
        pmErrorvsOrderTest(nkm, km, si);


        //------------------------------------------
        // Plot
        //------------------------------------------
        //        gnuplot_ctrl  *h1;
        //        h1 = gnuplot_init();
        //        //coeff_plot(h1, &SEML);
        //        //potential_plot(h1, &SEML);
        //        pij_plot(h1);
        //        //pij_plot(5, 5, h1);
        //        //pij_plot(2, 2, h1);
        //        char ch;
        //        printf("Press ENTER to close the gnuplot window(s)\n");
        //        scanf("%c",&ch);
        //        gnuplot_close(h1);

        //pmProjTest(si);
        //pmSmallDivisors(5e-2);
        //pmTestPrec(0.1);
        //pmErrorvsOrderTest();
        //pmOfsOrderTest(10);
        //pmNorms();
        //pmContributions();
        //pmTestIC();
        break;
    }

    //------------------------------------------------
    // Poincare maps
    //------------------------------------------------
    case 4:
    {
        //------------------------------------------
        // Additionnal parameters in bash
        //------------------------------------------
        //1. Array of orders to test
        int nkm = toDigit(argv[index++]);
        int km[nkm];
        for(int i = 0; i<nkm; i++)
        {
            km[i] = toDigit(argv[index++]);
            if(km[i] > OFTS_ORDER)
            {
                cout << "Warning: a required order is greater than OFTS_ORDER: return." << endl;
                return -1;
            }
        }

        //2. Array of OFS order.
        int nkofs = toDigit(argv[index++]);
        int kofs[nkofs];
        for(int i = 0; i<nkofs; i++)
        {
            kofs[i] = toDigit(argv[index++]);
            if(kofs[i] > OFS_ORDER)
            {
                cout << "Warning: a required order is greater than OFS_ORDER: return." << endl;
                return -1;
            }
        }

        //3. PMAP parameters
        int  argpmapc  = toDigit(argv[index++]);
        cout << "Number of PMAP parameters & settings is " << argpmapc << endl;

        //------------------------------------------
        // Initialisation of the central manifold
        //------------------------------------------
        initCM(SEML);
        initCOC(SEML);


        //------------------------------------------
        // Pmap init
        // Bash line for parallel computation: export OMP_SET_NUM_THREADS=5
        //------------------------------------------
        Pmap pmap;
        if(MODEL_TYPE == M_RTBP)
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


        //------------------------
        //May evolve
        //------------------------
        pmap.projFreq   =  toDigit(argv[index++]);
        pmap.type       =  toDigit(argv[index++]);
        pmap.tf         =     atof(argv[index++]); //needs to be a big number of Pmap, not for Tmap

        //isGS and graph style?
        //helps enhance computation
        int isGS     =  toDigit(argv[index++]);
        if(isGS == 1 && SEML.pms == PMS_GRAPH) pmap.isGS = 1;
        else pmap.isGS = 0;

        pmap.order      =  toDigit(argv[index++]);
        pmap.ofs_order  =  toDigit(argv[index++]);
        pmap.max_events =  toDigit(argv[index++]);

        //Specific case of the initial time, that
        //may be initialized as a root of a coefficient
        //of the COC matrix
        pmap.t0         =  toDigit(argv[index++]);
        if(pmap.t0 == -1) pmap.t0 = pij(2, 5); //root of p36 is input is -1

        pmap.dHv        =  atof(argv[index++]);
        pmap.gsize      =  toDigit(argv[index++]);
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

        //-----------------------
        // Additionnal map settings
        //-----------------------
        int append = toDigit(argv[index++]);
        int isPlot = toDigit(argv[index++]);
        int isPar  = toDigit(argv[index++]);
        int method = toDigit(argv[index++]);

        cout << "append          = "  << append << endl;
        cout << "isPlot          = "  << isPlot << endl;
        cout << "isPar           = "  << isPar << endl;
        cout << "method          = "  << pmap.type << endl;

        //-----------------------
        // openMP settings
        //-----------------------
        int num_threads = toDigit(argv[index++]);
        omp_set_num_threads(num_threads);
        cout << "num_threads     = "  << num_threads << endl;

        //------------------------------------------
        // Pmap computation
        //------------------------------------------
        switch(pmap.type)
        {
        case PMAP:
        {
            pmap.tt =  1.0/pmap.projFreq*SEML.us.T; //for pmaps, pmap.tt is the period of projection
            pmap_build(pmap, append, method, isPlot, isPar);
            break;
        }
        case TMAP:
        {
            tmap_build(pmap, append, method, isPlot, isPar);
            break;
        }


        //-----------------------
        // Precision map
        // for various orders
        //
        //  €€ TODO: add a parameter to take into account
        //   the different possibilities for error maps
        //
        //-----------------------
        case EMAP:
        {
            for(int i = 0; i < nkm; i++)
            {
                pmap.order = km[i];
                pmap_precision(pmap, append, isPar);
            }
            break;
        }
        //-----------------------
        // Invariance error map
        // for various orders
        //
        //  €€ TODO: add a parameter to take into account
        //   the different possibilities for invariance error maps
        //
        //-----------------------
        case IMAP:
        {
            for(int i = 0; i < nkm; i++)
            {
                pmap.order = km[i];
                //pmap_invariance_error(pmap, append, isPar, pmap.dHv);
                pmap_invariance_error_random(pmap, append, isPar, pmap.dHv);
                //pmap_test_error(pmap, append, isPar, pmap.dHv);
            }
            break;
        }

        }
        break; //and of case 4: pmaps
    }

    //------------------------------------------------
    // COCs
    //------------------------------------------------
    case 5:
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
          oftsh_test();

        break;
    }

    //------------------------------------------------
    // Single trajectories
    //------------------------------------------------
    case 6:
    {
        //------------------------------------------
        // Additionnal parameters in bash
        //------------------------------------------
        //3. PMAP parameters
        int  argpmapc  = toDigit(argv[index++]);
        cout << "Number of PMAP parameters & settings is " << argpmapc << endl;

        //------------------------------------------
        // Initialisation of the central manifold
        //------------------------------------------
        initCM(SEML);
        initCOC(SEML);

        //------------------------------------------
        // Pmap init
        // Bash line for parallel computation: export OMP_SET_NUM_THREADS=5
        //------------------------------------------
        Pmap pmap;
        if(MODEL_TYPE == M_RTBP)
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
            pmap.proot      =  1e-13;
            pmap.threshold  =  1e-6;
            pmap.T          =  SEML.us.T;
            pmap.vdim       =  2;
            pmap.maxRad     =  5.0;
        }

        //May evolve
        //------------------------
        pmap.type       =  toDigit(argv[index++]);
        pmap.tf         =  atof(argv[index++]); //needs to be a big number of Pmap, not for Tmap
        pmap.tt         =  atof(argv[index++])*SEML.us.T;


        //isGS and graph style?
        //helps enhance computation
        int isGS     =  toDigit(argv[index++]);
        pmap.isGS = 0;
        if(isGS == 1)
        {
            switch(SEML.cs.manType)
            {
            case MAN_CENTER:
                //if MAN_CENTER, we use
                if(SEML.pms == PMS_GRAPH) pmap.isGS = 1;
                break;
            case MAN_CENTER_S:
            case MAN_CENTER_U:
                pmap.isGS = 1;
                break;
            }
        }
        pmap.order      =  OFTS_ORDER;
        pmap.ofs_order  =  OFS_ORDER;
        pmap.max_events =  toDigit(argv[index++]);

        //Specific case of the initial time, that
        //may be initialized as a root of a coefficient
        //of the COC matrix
        pmap.t0         =  atof(argv[index++]);
        if(pmap.t0 == -1) pmap.t0 = pij(2, 5); //root of p36 is input is -1
        pmap.dHv        =  atof(argv[index++]);

        cout << "pmap.type       = "  << pmap.type << endl;
        cout << "pmap.tf         = "  << pmap.tf << endl;
        cout << "pmap.isGS  = "  << pmap.isGS << endl;
        cout << "pmap.order      = "  << pmap.order << endl;
        cout << "pmap.ofs_order  = "  << pmap.ofs_order  << endl;
        cout << "pmap.max_events = "  << pmap.max_events << endl;
        cout << "pmap.t0         = "  << pmap.t0  << endl;
        cout << "pmap.dHv        = "  << pmap.dHv   << endl;

        //-----------------------
        // openMP settings
        //-----------------------
        int num_threads = toDigit(argv[index++]);
        omp_set_num_threads(num_threads);
        cout << "num_threads       = "  << num_threads << endl;

        //------------------------------------------
        // Initial conditions
        //------------------------------------------
        double si[REDUCED_NV];
        for(int i = 0; i < REDUCED_NV; i++)
        {
            si[i] = atof(argv[index++]);
        }


        //------------------------------------------
        // Trajectory computation
        //------------------------------------------
        trajectory_CU(si, pmap, 1, false);

        break; //and of case 6: trajectories
    }

    }
    return 0;


    //-------------------------------------------------------------------------------------------------------
    // MULTIMIN_TEST
    //-------------------------------------------------------------------------------------------------------
    //test_frprmn();

    //-------------------------------------------------------------------------------------------------------
    // Tests
    //-------------------------------------------------------------------------------------------------------
    //testLegendreRecurrence_OFS();
    //testCOC();
    //testDotCOC();
    //testIntCOC();

    //-------------------------------------------------------------------------------------------------------
    // Unitary tests of routines
    //-------------------------------------------------------------------------------------------------------
    //ofs_test();
    //oftsh_test();
    //ofts_test();

    //-------------------------------------------------------------------------------------------------------
    // Retrieving the Hamiltonian from the VF (not working)
    //-------------------------------------------------------------------------------------------------------
    /*
    vector_fprinf_0(Fh, SEML.cs.F_NF+"fh");

    vector<Oftsc> H = vector<Oftsc>(4);
    vector<Oftsc> Ham = vector<Oftsc>(1);

    Oftsc *F1     = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);

    Oftsc *dF1dp2 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *dF1dq1 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *dF1dq2 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);

    Oftsc *dGdp2  = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *G2     = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *dG2dq1 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *dG2dq2 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);

    Oftsc *dJdq1  = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *J1     = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *dJ1dq2 = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);

    Oftsc *dKdq2  = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);
    Oftsc *K2     = new Oftsc(REDUCED_NV, OFTS_ORDER, OFS_NV, OFS_ORDER);


    //Fh[0] =  dH/ds3 = dH/dp1
    //Fh[1] =  dH/ds4 = dH/dp2
    //Fh[2] = -dH/ds1 = -dH/dq1
    //Fh[3] = -dH/ds2 = -dH/dq2

    //Building F1
    F1->sprim(Fh[0], 3);  //dF1/dp1 = dH/dp1
    dF1dp2->der(*F1, 4);  //dF1/dp2
    dF1dq1->der(*F1, 1);  //dF1/dq1
    dF1dq2->der(*F1, 2);  //dF1/dq2

    //Building dG/dp2 = dH/dp2 - dF1/dp2
    dGdp2->ofts_fsum_u(Fh[1], 1.0+0.0*I, *dF1dp2, -1.0+0.0*I);
    //Building G2
    G2->sprim(*dGdp2, 4);
    dG2dq1->der(*G2, 1);  //dG2/dq1
    dG2dq2->der(*G2, 2);  //dG2/dq2

    //Building dJ/dq1 = dH/dq1 - dF1/dq1 - dG2/dq1
    Fh[2].conjugate();
    Fh[2].conjugate(1);
    dJdq1->ofts_fsum_u(Fh[2], -1.0+0.0*I, *dF1dq1, -1.0+0.0*I);
    dJdq1->ofts_smult_u(*dG2dq1, -1.0+0.0*I);
    //Building J1
    J1->sprim(*dJdq1, 1);
    dJ1dq2->der(*J1, 2);  //dF1/dq2


    //Building dK/dq2 = dJ/dq2 -  dJ1/dq2 = dH/dq2 - dF1/dq2 - dG2/dq2 - dJ1/dq2
    Fh[3].conjugate();
    Fh[3].conjugate(1);
    dKdq2->ofts_fsum_u(Fh[3], -1.0+0.0*I, *dF1dq2, -1.0+0.0*I);
    dKdq2->ofts_smult_u(*dG2dq2, -1.0+0.0*I);
    dKdq2->ofts_smult_u(*dJ1dq2, -1.0+0.0*I);

    //Building K
    K2->sprim(*dKdq2, 2);

    //Building H
    Ham[0].ofts_smult_u(*F1, 1.0+0.0*I);
    Ham[0].ofts_smult_u(*G2, 1.0+0.0*I);
    Ham[0].ofts_smult_u(*J1, 1.0+0.0*I);
    Ham[0].ofts_smult_u(*K2, 1.0+0.0*I);

    H[0].der(Ham[0], 3);
    H[1].der(Ham[0], 4);
    H[2].der(Ham[0], 1);
    H[3].der(Ham[0], 2);



    vector_fprinf_0(H, SEML.cs.F_NF+"H");
    vector_fprinf_0(Fh, SEML.cs.F_NF+"Fh");
    vector_fprinf_0(Ham, SEML.cs.F_NF+"Ham");
    */


    //----------------------------------------------------------------------------------------------------------
    //Continuation methods (backup)
    //----------------------------------------------------------------------------------------------------------
    /*
    char ch;            //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    //Integration tools
    //-------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    if(model.model1.isNormalized)
        sys.function = qbfbp_vfn_cont;
    else
        sys.function = qbfbp_vf_cont;
    sys.jacobian      = NULL;
    sys.dimension     = 42;
    sys.params        = &model;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, 1e-6, 1e-16, 1e-16);

    //Differential correction to get the dynamical equivalent of the libration point (lpdyneq)
    double y0[42];
    lpdyneq_cont(d, y0, h1);

    //Plotting the refined solution
    int Npoints = 5000;
    odePlot(y0,  model.model2.qbcp.T*0.5, d, h1, Npoints, 8);
    odePlot(y0, -model.model2.qbcp.T*0.5, d, h1, Npoints, 8);


    gnuplot_close(h1);
    */

}



