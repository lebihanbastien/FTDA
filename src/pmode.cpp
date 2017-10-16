#include "pmode.h"

/**
 * \file pmode.cpp
 * \brief Integration of the equations of motion for the outputs of the parameterization method.
 * \author BLB.
 * \date 2016
 * \version 1.0
 */


//----------------------------------------------------------------------------------------
//
//          Integration of the pm
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the reduce vector field dot(s) = f(s,t)
 *  \param t the current time
 *  \param y the complex state s (CCM), with real and imag parts separated (CCM8 form)
 *  \param f the complex reduced vector field f(s,t), with real and imag parts separated (CCM8 form)
 *  \param params_void a (void) pointer towards a RVF structure that contains:
 *           (1) the Fourier-Taylor expansions of the rvf
 *           (2) the order of the evaluation
 *           (3) the pulsation n of the system
 *           (4) an OFS temp variable.
 *
 *     Note that the use of void *params_void instead of the direct use of an RVF structure is to comply with GSL ODE integration tools,
 *     that require a vector field of the exact form : int vf(double t, const double y[], double f[], void *params_void)
 *
 **/
int qbfbp_fh(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------
    //Initialization
    //-------------------------------
    RVF *rvf = (RVF*) params_void;

    //-------------------------------
    //Evaluation of the reduced vector field at order rvf->order
    //-------------------------------
    CCM8toRVF8(y, t, rvf->n, rvf->order, rvf->ofs_order, *rvf->fh, *rvf->ofs, f);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------------------------------
//
//          Tests
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Test of the pm of the central manifold of L1,2 on a given orbit, through the computation of various errors (eO, eI, eH), along an orbit initialized by the pm.
 *         Various orders are tested, among the available values.
 *
 *   Note that: the expansions FW, DWf and Wdot are taken from files,
 *   whereas the parameterization itself CM and CMh, and the reduced vector field Fh
 *   are taken from global objects, defined in init.cpp.
 *
 *   Requires initCM().
 *
 **/
void pmErrorvsOrderTest(int nkm, int km[], double si[])
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_PMS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;


    //------------------------------------------------------------------------------------
    // Initialisation of the expansion
    //------------------------------------------------------------------------------------
    vector<Oftsc> FW(6);
    vector<Oftsc> DWf(6);
    vector<Oftsc> Wdot(6);

    readVOFTS_bin(FW,   F_GS+"FW/C_FW",       OFS_ORDER);
    readVOFTS_bin(DWf,  F_GS+"DWf/C_DWf",     OFS_ORDER);
    readVOFTS_bin(Wdot, F_GS+"Wdot/C_Wdot",   OFS_ORDER);

    //------------------------------------------------------------------------------------
    // Initialisation of the COC
    //------------------------------------------------------------------------------------
    matrix<Ofsc> P(6,6);
    matrix<Ofsc> PC(6,6);
    matrix<Ofsc> CQ(6,6);
    matrix<Ofsc> Q(6,6);
    vector<Ofsc> V(6);
    initCOC(P, PC, Q, CQ, V, SEML);

    //------------------------------------------------------------------------------------
    //Integration tools
    //------------------------------------------------------------------------------------
    //For dot(z) = F(z)
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_novar;
    sys.jacobian = NULL;
    sys.dimension = 6;
    sys.params = &SEML;
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());

    //For dot(s) = fh(s)
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh        = &Fh;
    rvf.ofs       = &AUX;
    rvf.order     = OFTS_ORDER;
    rvf.n         = SEML.us.n;
    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*REDUCED_NV;
    sys_fh.params    = &rvf;
    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());


    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the obtained PM             " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Misc variables
    //------------------------------------------------------------------------------------
    //Int variables
    double st0[REDUCED_NV];
    //temp
    Ofsc BUX(OFS_ORDER);


    //------------------------------------------------------------------------------------
    // Orders of the expansions
    //------------------------------------------------------------------------------------
    int order = 5;
    int ofs_order = OFS_ORDER;

    //------------------------------------------------------------------------------------
    // Maximum time
    //------------------------------------------------------------------------------------
    double tmax = (SEML.model == Csts::CRTBP)? 2*M_PI: SEML.us.T;
    tmax *= (SEML.coordsys == Csts::SEM && SEML.model != Csts::CRTBP)? 10.0: 1.0;
    switch(SEML.cs.manType)
    {
        case Csts::MAN_CENTER_S:
            //If we are in the center-stable manifold, we set negative time.
            tmax = -fabs(tmax);
    }

    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    string st0s = "";
    for(int p = 0; p < REDUCED_NV; p++) st0[p] = si[p];

    //------------------------------------------------------------------------------------
    //For plotting
    //------------------------------------------------------------------------------------
    //Gnuplot handlers
    gnuplot_ctrl  **ht = (gnuplot_ctrl**) calloc(6, sizeof(gnuplot_ctrl*));
    for(int i = 0; i <6; i++) ht[i] = gnuplot_init();

    //------------------------------------------------------------------------------------
    //Loop on order
    //------------------------------------------------------------------------------------
    int color = 1;
    for(int ii = 0; ii< nkm; ii++)
    {
        //order of the expansions
        order = km[ii];
        //--------------------------------------------------------------------------------
        // eO, eI, and orbit
        //--------------------------------------------------------------------------------
        rvf.order = order;
        tic();
        errorPlot(st0, CMh, PC, V, FW, DWf, Wdot, 0.5*tmax, d, d_fh, 1000, SEML, order, ofs_order, ht, color++);
        cout << "errorPlot in: " << toc() << endl;
    }

    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);


    //Closing handlers
    for(int i = 0; i <6; i++) gnuplot_close(ht[i]);
    free(ht);


    free(d);
    free(d_fh);

}

/**
 *  \brief Test of the pm of the central manifold of L1,2 on a given orbit, through the computation of the orbital error eO, along an orbit initialized by the pm.
 *         Various orders are tested, among the available values.
 *
 *   Requires initCM().
 *
 **/
void pmEOvsOrderTest(int nkm, int km[], double si[])
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_PMS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;

    //------------------------------------------------------------------------------------
    //Integration tools
    //------------------------------------------------------------------------------------
    //For dot(z) = F(z)
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_novar;
    sys.jacobian = NULL;
    sys.dimension = 6;
    sys.params = &SEML;
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());

    //For dot(s) = fh(s)
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh        = &Fh;
    rvf.ofs       = &AUX;
    rvf.order     = OFTS_ORDER;
    rvf.n         = SEML.us.n;
    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*REDUCED_NV;
    sys_fh.params    = &rvf;
    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());


    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the obtained PM             " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Misc variables
    //------------------------------------------------------------------------------------
    //Int variables
    double st0[REDUCED_NV];
    //temp
    Ofsc BUX(OFS_ORDER);

    //------------------------------------------------------------------------------------
    // Orders of the expansions
    //------------------------------------------------------------------------------------
    int order = 5;
    int ofs_order = OFS_ORDER;

    //------------------------------------------------------------------------------------
    // Maximum time
    //------------------------------------------------------------------------------------
    double tmax = (SEML.model == Csts::CRTBP)? 2*M_PI: SEML.us.T;

    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    string st0s = "";
    for(int p = 0; p < REDUCED_NV; p++) st0[p] = si[p];

    //------------------------------------------------------------------------------------
    //For plotting
    //------------------------------------------------------------------------------------
    //Gnuplot handlers
    gnuplot_ctrl  **ht = (gnuplot_ctrl**) calloc(6, sizeof(gnuplot_ctrl*));
    for(int i = 0; i < 6; i++) ht[i] = gnuplot_init();

    //------------------------------------------------------------------------------------
    //Loop on order
    //------------------------------------------------------------------------------------
    int color = 1;
    for(int ii = 0; ii< nkm; ii++)
    {
        //order of the expansions
        order = km[ii];
        //------------------------------------------------------------------------------------
        // eO, eI, and orbit
        //------------------------------------------------------------------------------------
        rvf.order = order;
        tic();
        eOPlot(st0, CMh, Mcoc, Vcoc, tmax, d, d_fh, 500, SEML, order, ofs_order, ht, color++);
        cout << "eOPlot ended in: " << toc() << endl;
    }


    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);


    //Closing handlers
    for(int i = 0; i <6; i++) gnuplot_close(ht[i]);
    free(ht);


    free(d);
    free(d_fh);
}

/**
 *  \brief Test of the pm of the central manifold of L1,2 on a given orbit, through the computation of the orbital error eO, along an orbit initialized by the pm.
 *         Various OFS orders (of the Fourier expansions) are tested, among the available values.
 *   \param order the order of the Taylor expansions (OFTS objects).
 *   Requires initCM().
 *
 **/
void pmOfsOrderTest(int order)
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;


    //------------------------------------------------------------------------------------
    // Initialisation of the COC
    //------------------------------------------------------------------------------------
    matrix<Ofsc> P(6,6);
    matrix<Ofsc> PC(6,6);
    matrix<Ofsc> CQ(6,6);
    matrix<Ofsc> Q(6,6);
    vector<Ofsc> V(6);
    initCOC(P, PC, Q, CQ, V, SEML);

    //------------------------------------------------------------------------------------
    //Integration tools
    //------------------------------------------------------------------------------------
    //For dot(z) = F(z)
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_novar;
    sys.jacobian = NULL;
    sys.dimension = 6;
    sys.params = &SEML;
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());

    //For dot(s) = fh(s)
    RVF rvf;
    Ofsc AUX(OFS_ORDER);
    rvf.fh        = &Fh;
    rvf.ofs       = &AUX;
    rvf.order     = OFTS_ORDER;
    rvf.ofs_order = OFS_ORDER;
    rvf.n         = SEML.us.n;
    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*REDUCED_NV;
    sys_fh.params    = &rvf;
    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());


    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the obtained PM             " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Misc variables
    //------------------------------------------------------------------------------------
    //Time variables
    double theta1;
    theta1 = 2*M_PI;
    //Int variables
    double st0[4];
    //temp
    Ofsc BUX(OFS_ORDER);
    //Keymap for loop on order
    int km[4];
    km[0] = 10;
    km[1] = 15;
    km[2] = 20;
    km[3] = 30;

    //------------------------------------------------------------------------------------
    // Orders of the Fourier expansions
    //------------------------------------------------------------------------------------
    int ofs_order;

    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    for(int p = 0; p < 4; p++) st0[p] = 0.0;
    st0[0] = 20.0;

    //------------------------------------------------------------------------------------
    //For plotting
    //------------------------------------------------------------------------------------
    //Gnuplot handlers
    gnuplot_ctrl  **ht = (gnuplot_ctrl**) calloc(6, sizeof(gnuplot_ctrl*));
    for(int i = 0; i <6; i++) ht[i] = gnuplot_init();

    //------------------------------------------------------------------------------------
    //Loop on order
    //------------------------------------------------------------------------------------
    int color = 1;
    for(int ii = 0; ii< 4; ii++)
    {
        //order of the expansions
        ofs_order = km[ii];
        //------------------------------------------------------------------------------------
        // eO, eI, and orbit
        //------------------------------------------------------------------------------------
        rvf.order     = order;
        rvf.ofs_order = ofs_order;
        tic();
        eOPlot(st0, CMh, PC, V, theta1/SEML.us.n, d, d_fh, 500, SEML, order, ofs_order, ht, color++);
        cout << "eOPlot in: " << toc() << endl;
    }


    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);


    //Closing handlers
    for(int i = 0; i <6; i++) gnuplot_close(ht[i]);
    free(ht);


    free(d);
    free(d_fh);
}

/**
 *  \brief Evaluates the contributions of each order in W to the computation of W(s,t), with an arbitrary state (s,t)
 **/
void pmContributions()
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);

    //------------------------------------------------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> W(6);
    vector<Oftsc> Wh(6);
    readVOFTS_txt(Wh, F_GS+"W/Wh", OFS_ORDER);
    readVOFTS_txt(W,  F_GS+"W/W", OFS_ORDER);


    //------------------------------------------------------------------------------------
    //Contribution of each order to the evaluation
    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    string st0s;
    cdouble st0[4], s0[4];
    for(int p = 0; p < 4; p++) st0[p] = 5e-0+0.0*I;
    st0s = "5e-0";

    //------------------------------------------------------------------------------------
    //"Realification": s0 = REAL(st0)
    //------------------------------------------------------------------------------------
    s0[0] = 1.0/sqrt(2)*( st0[0]   - st0[2]*I);
    s0[2] = 1.0/sqrt(2)*( st0[2]   - st0[0]*I);
    s0[1] = 1.0/sqrt(2)*( st0[1]   - st0[3]*I);
    s0[3] = 1.0/sqrt(2)*( st0[3]   - st0[1]*I);

    //------------------------------------------------------------------------------------
    // Evaluate W
    //------------------------------------------------------------------------------------
    //z0t = W(s0), z0 = z0t(0.0)
    Ofsc AUX(OFS_ORDER);
    for(int k = 0; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < 6; p++)
        {
            W[p].contribution(s0, AUX, k);
            cout << AUX.l1norm() << " ";
        }
        cout << endl;
    }
}

/**
 *  \brief Evaluates the l1/linf norm of each order in W.
 **/
void pmNorms()
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;


    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    //------------------------------------------------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> W(6);
    vector<Oftsc> Wh(6);
    readVOFTS_bin(Wh, F_GS+"W/Wh", OFS_ORDER);
    readVOFTS_bin(W,  F_GS+"W/W", OFS_ORDER);

    //------------------------------------------------------------------------------------
    //L1-norm
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5, "set logscale y");
    gnuplot_cmd(h5, "set format y \"1e\%%L\"");
    gnuplot_cmd(h5, "set grid");
    gnuplot_cmd(h5, "set title \"l_1(order)\" ");


    double l1n[OFTS_ORDER];
    double l1nMax[OFTS_ORDER+1];
    double kc[OFTS_ORDER+1];
    cout << "---------------------------"<< endl;
    cout << "Linf-norm of W:" << endl;

    for(int k = 0; k <= OFTS_ORDER; k++)
    {
        kc[k] = k;
        l1nMax[k] = 0.0;
        for(int p = 0; p < Csts::NV; p++)
        {
            //l1n[k] = Wh[p].linfnorm(k);
            //l1n[k] = W[p].l1norm(k);
            l1n[k] = W[p].linfnorm(k);
            if(l1n[k] > l1nMax[k]) l1nMax[k] = l1n[k];
        }
        cout << k << " " << l1nMax[k] << endl;
    }

    gnuplot_plot_xy(h5, kc, l1nMax, OFTS_ORDER+1, (char*)"", "points", "1", "2", 1);

    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
}

/**
 *  \brief Small divisors under a certain value
 **/
void pmSmallDivisors(double sdmax)
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_NF    = SEML.cs.F_NF;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;


    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    //------------------------------------------------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------------------------------------------------
    vector<Oftsc> smallDiv(6);
    readVOFTS_txt(smallDiv, F_NF+"W/smallDiv", OFS_ORDER);

    //------------------------------------------------------------------------------------
    //Small divisors
    //------------------------------------------------------------------------------------
    cout << "---------------------------"<< endl;
    cout << "Number of small divisors other than order 0:" << endl;
    for(int k = 2; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < Csts::NV; p++)
        {
            cout <<  smallDiv[p].nsd(k, OFS_ORDER, sdmax) - smallDiv[p].nsd(k, 0, sdmax) << " ";

        }
        cout << endl;

    }



    cout << "---------------------------"<< endl;
    cout << "Number of small divisors at order 0:" << endl;
    for(int k = 2; k <= OFTS_ORDER; k++)
    {
        cout << k << " ";
        for(int p = 0; p < Csts::NV; p++)
        {
            cout <<  smallDiv[p].nsd(k, 0, sdmax) << " ";

        }
        cout << endl;

    }
}

/**
 *  \brief Evaluates the variations of the initial conditions (IC) wrt to the order in W(s,t), with an arbitrary state (s,t)
 **/
void pmTestIC()
{
    //------------------------------------------------------------------------------------
    // Initialisation of the central manifold
    //------------------------------------------------------------------------------------
    initCM(SEML);

    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;


    //Initialization of the configuration
    double st0[4];
    for(int p = 0; p < 4; p++) st0[p] = 0.0;
    st0[0] = 4;

    //"Realification": Xeval = REAL(st0)
    cdouble Xeval[4];
    Xeval[0] = 1.0/sqrt(2)*( st0[0]   - st0[2]*I);
    Xeval[2] = 1.0/sqrt(2)*( st0[2]   - st0[0]*I);
    Xeval[1] = 1.0/sqrt(2)*( st0[1]   - st0[3]*I);
    Xeval[3] = 1.0/sqrt(2)*( st0[3]   - st0[1]*I);

    //------------------------------------------------------------------------------------
    // Evaluate Wc = |W(Xeval, 0)|
    //------------------------------------------------------------------------------------
    vector<Ofsc> z0t(6);
    double z0[4][6];

    //Keymap
    int keyMap[4];
    int m;
    keyMap[0] = 2;
    keyMap[1] = 5;
    keyMap[2] = 10;
    keyMap[3] = OFTS_ORDER;

    for(int i = 0; i < 4; i++)
    {
        m = keyMap[i];
        for(int p = 0; p < 6; p++)
        {
            CM[p].evaluate(Xeval, z0t[p], m, OFS_ORDER);
            z0[i][p] = creal(z0t[p].evaluate(0.0));
        }
    }

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    cout << "IC " << m << ":" << endl;
    cout << "------------------------------" << endl;
    cout << " 2-5        5-10         10-15 " << endl;
    for(int p = 0; p <6; p++)
    {
        for(int i = 0; i < 3; i++)cout << z0[i][p]-z0[i+1][p] << "   ";
        cout << endl;
    }

}


//----------------------------------------------------------------------------------------
//
//          Plot in tests
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Computations of various errors (eO, eI, eH) along an orbit initialized by the pm of the center manifold W(s,t)
 **/
int errorPlot(const double st0[],      //RCM initial conditions
              vector<Oftsc>& Wh,       //TFC manifold
              matrix<Ofsc>& PC,        //COC matrix: z = PC*zh+V
              vector<Ofsc>& V,         //COC vector: z = PC*zh+V
              vector<Oftsc>& FW,       //NC vector field
              vector<Oftsc>& DWf,      //Jacobian of FW
              vector<Oftsc>& Wdot,     //Partial derivative of W wrt time
              double t1,               //Final integration tim
              gsl_odeiv2_driver *dnc,  //driver for NC integration
              gsl_odeiv2_driver *drvf, //drive for RVF integration (reduced vector field)
              int Npoints,             //Number of points on which the errors are estimated
              FBPL& fbpl,          //current QBCP
              int order,               //Order for the eval of the OFTS objects
              int ofs_order,           //Order for the eval of the OFS objects
              gnuplot_ctrl  **ht,      //Gnuplot handlers
              int color)               //Color of the plots
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = fbpl.cs.F_GS;
    string F_PLOT  = fbpl.cs.F_PLOT;
    string F_COC   = fbpl.cs.F_COC;


    //------------------------------------------------------------------------------------
    //Reset integrator
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(dnc);
    gsl_odeiv2_driver_reset(drvf);

    //------------------------------------------------------------------------------------
    //Variables for plotting
    //------------------------------------------------------------------------------------
    //Time
    double tc[Npoints+1];
    //Orbits
    double xc[Npoints+1], yc[Npoints+1], zc[Npoints+1];
    //Orbit in TFC
    double xsc[Npoints+1], ysc[Npoints+1], zsc[Npoints+1];
    //Errors
    double eOc[Npoints+1], eIc[Npoints+1], eHc[Npoints+1];

    //------------------------------------------------------------------------------------
    //Energy
    //------------------------------------------------------------------------------------
    double H, H0;

    //------------------------------------------------------------------------------------
    //Errors
    //------------------------------------------------------------------------------------
    cdouble eId;
    double eIm, eOm;
    double eI[6], eO[6];

    //------------------------------------------------------------------------------------
    // Inner state (NC, TFC...)
    //------------------------------------------------------------------------------------
    double  s1rcm[REDUCED_NV];    //RCM
    cdouble s1ccm[REDUCED_NV];    //CCM
    double  s1ccm8[2*REDUCED_NV]; //CCM8
    double  z1ncr[6], z0ncr[6];   //NC, integrated with NC vector field
    double  z1ncr2[6];            //NC, but computed from z = W(s,t)
    double  z1em[6], z0em[6];     //EM
    double  fnc[6];               //NC vector field
    double  W0;                   //Euclidian norm in NC coordinates
    double  W0km;                 //Euclidian norm in km
    Ofsc AUX;                     //OFS temp variable

    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    //------------------------------------------------------------------------------------
    double n = fbpl.us.n;

    //------------------------------------------------------------------------------------
    //Hamiltonian at the origin
    //------------------------------------------------------------------------------------
    double stLi[REDUCED_NV];
    for(int i = 0; i < REDUCED_NV; i++) stLi[i] = 0.0;
    RCMtoNCbyTFC(stLi, 0.0, n, order, ofs_order, Wh, PC, V, z0ncr, false);
    NCtoSYS(0.0, z0ncr, z0em, (FBPL*) dnc->sys->params);  //Hamiltonian in sys coordinates (for reference)
    double HLi = qbfbp_H(0.0, z0em, dnc->sys->params);

    //------------------------------------------------------------------------------------
    // RCM to NC for NC initial conditions
    //------------------------------------------------------------------------------------
    RCMtoNCbyTFC(st0, 0.0, n, order, ofs_order, Wh, PC, V, z1ncr, false);

    //------------------------------------------------------------------------------------
    // RCM to CCM8 for CCM initial conditions
    //------------------------------------------------------------------------------------
    RCMtoCCM8(st0, s1ccm8);

    //------------------------------------------------------------------------------------
    // Euclidian Norm
    //------------------------------------------------------------------------------------
    W0 = ENorm(z1ncr, 3);
    W0km = W0*fbpl.cs.gamma*fbpl.cs.cr3bp.L;

    //------------------------------------------------------------------------------------
    // Print Initial Conditions
    //------------------------------------------------------------------------------------
    //Initial conditions
    cout << "With current s0 and order, the initial conditions are: " << endl;
    for(int p = 0; p < 6; p++) cout << p << "  " << z1ncr[p] << endl;
    //Euclidian norm in km
    cout << "Approximated distance from Li [km]: " << W0km << endl;

    //------------------------------------------------------------------------------------
    // Initial relative Hamiltonian
    //-----------------------------------------
    //H0 = qbfbp_Hn(0.0, z1ncr, dnc->sys->params)-HLi0;           //Hamiltonian
    NCtoSYS(0.0, z1ncr, z1em, (FBPL*) dnc->sys->params);        //Hamiltonian in EM coordinates (for reference)
    H0 = qbfbp_H(0.0, z1em, dnc->sys->params) -HLi ;
    cout << "H in EM coordinates: " << H0 << endl;

    //------------------------------------------------------------------------------------
    // For plotting (first value)
    //------------------------------------------------------------------------------------
    //Time
    tc[0] = 0.0;
    //NC
    xc[0] = z1ncr[0];
    yc[0] = z1ncr[1];
    zc[0] = z1ncr[2];
    //NC from manifold
    xsc[0] = z1ncr[0];
    ysc[0] = z1ncr[1];
    zsc[0] = z1ncr[2];
    //Errors
    eIc[0]  = 0.0;
    eOc[0]  = 0.0;
    eHc[0]  = sqrt((H0 - HLi)*(H0 - HLi));

    //------------------------------------------------------------------------------------
    // Loop on time
    //------------------------------------------------------------------------------------
    double t = 0.0;
    double t2 = 0.0;
    double ti = 0;
    if(t1 <0)
    {
        dnc->h = -dnc->h;
        drvf->h = -drvf->h;
    }
    for(int i =0; i<= Npoints; i++)
    {
        //------------------------------------------------------------------------------------
        // Apply driver
        //------------------------------------------------------------------------------------
        ti = (double) i * t1 / Npoints;
        gsl_odeiv2_driver_apply (dnc, &t, ti, z1ncr);
        gsl_odeiv2_driver_apply(drvf, &t2, ti, s1ccm8);

        //------------------------------------------------------------------------------------
        // Comparison
        //------------------------------------------------------------------------------------

        //---------------------
        // Using the PM: z = W(s,t)
        // Computed in z1ncr2
        //---------------------
        //CCM8 to RCM
        CCM8toRCM(s1ccm8, s1rcm);
        //CCM8 to CCM
        CCM8toCCM(s1ccm8, s1ccm);
        //z1ncr2 = W(s1rcm, ti)
        RCMtoNCbyTFC(s1rcm, ti, n, order, ofs_order, Wh, PC, V, z1ncr2, false);

        //Evaluating the vector field at z1ncr2 = W(s1rcm, ti)
        qbfbp_vfn_novar(ti, z1ncr2, fnc, dnc->sys->params);

        //Hamiltonian, at the origin
        RCMtoNCbyTFC(stLi, ti, n, order, ofs_order, Wh, PC, V, z0ncr, false);
        NCtoSYS(ti, z0ncr, z0em, (FBPL*) dnc->sys->params);  //Hamiltonian in sys coordinates (for reference)
        HLi = qbfbp_H(ti, z0em, dnc->sys->params);


        //Hamiltonian, taken at z1ncr2 = W(s1rcm, ti)
        NCtoSYS(ti, z1ncr2, z1em, (FBPL*) dnc->sys->params);
        H  = qbfbp_H(ti, z1em, dnc->sys->params);


        //---------------------
        //Error computation
        //---------------------
        for(int p = 0; p < Csts::NV; p++)
        {
            //eO = |z(t) - W(s(t), t)|
            eO[p] = cabs(z1ncr[p] - z1ncr2[p]);
            //eI
            //-------------------------------
            //What is the definition of eI?
            // either :
            // - F(W(s,t)) - FW(s,t)
            // - F(W(s,t)) - DWf(s,t) - Wdot(s,t)
            //------------------------------
            eId = fnc[p]+0.0*I;
            //FW[p].evaluate(s1ccm, AUX, order);
            //eId = AUX.evaluate(n*ti);
            //OR
            DWf[p].evaluate(s1ccm, AUX, order, ofs_order);
            eId -= AUX.evaluate(n*ti);
            Wdot[p].evaluate(s1ccm, AUX, order, ofs_order);
            eId -= AUX.evaluate(n*ti);
            eI[p] = cabs(eId);
        }


        //Taking the maximum (infinity norm)
        eOm = eO[0];
        eIm = eI[0];
        for(int p = 1; p < Csts::NV; p++)
        {
            if(eO[p] > eOm) eOm = eO[p];
            if(eI[p] > eIm) eIm = eI[p];
        }

        //--------------------------------------------------------------------------------
        //Store for plotting
        //--------------------------------------------------------------------------------
        //Timde
        tc[i]  = ti;
        //NC
        xc[i]  = z1ncr[0];
        yc[i]  = z1ncr[1];
        zc[i]  = z1ncr[2];
        //NC from manifold
        xsc[i] = z1ncr2[0];
        ysc[i] = z1ncr2[1];
        zsc[i] = z1ncr2[2];
        //Errors
        eOc[i] = eOm;
        eIc[i] = eIm;
        eHc[i] = sqrt((H - HLi)*(H - HLi));
    }
    if(t1 <0)
    {
        dnc->h = -dnc->h;
        drvf->h = -drvf->h;
    }


    //------------------------------------------------------------------------------------
    //Define the output name
    //------------------------------------------------------------------------------------
    ifstream readStream;
    string ss1, ssW, st0s;
    ss1  = static_cast<ostringstream*>( &(ostringstream() << order) )->str();
    ssW  = static_cast<ostringstream*>( &(ostringstream() << W0) )->str();
    if(SEML.cs.manType == Csts::MAN_CENTER)
    {
        st0s = static_cast<ostringstream*>( &(ostringstream() << creal(st0[0]+st0[1]+st0[2]+st0[3])) )->str();
    }
    else
    {
        st0s = static_cast<ostringstream*>( &(ostringstream() << creal(st0[0]+st0[1]+st0[2]+st0[3])) )->str();
        string s5s = static_cast<ostringstream*>( &(ostringstream() << creal(st0[4])) )->str();
        st0s = st0s+"_s5_"+s5s;
    }
    st0s = st0s+"_test";


    //------------------------------------------------------------------------------------
    //In txt files
    //------------------------------------------------------------------------------------
    string eOstr   = (F_PLOT+"orbits/"+"eO_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eIstr   = (F_PLOT+"orbits/"+"eI_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eHstr   = (F_PLOT+"orbits/"+"eH_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eXYZstr = (F_PLOT+"orbits/"+"XYZ_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string ePMstr  = (F_PLOT+"orbits/"+"XYZ_PM_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    gnuplot_fplot_xy(tc, eOc, Npoints+1, eOstr.c_str());                //eO
    gnuplot_fplot_xy(tc, eIc, Npoints+1, eIstr.c_str());                //eI
    gnuplot_fplot_xy(tc, eHc, Npoints+1, eHstr.c_str());                //eH
    gnuplot_fplot_txyz(tc, xc, yc,   zc,   Npoints, eXYZstr.c_str());   //orbit
    gnuplot_fplot_txyz(tc, xsc, ysc, zsc,  Npoints, ePMstr.c_str());    //PM orbit


    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    //Logscale if necessary (all errors)
    gnuplot_cmd(ht[0], "set logscale y");
    gnuplot_cmd(ht[0], "set format y \"1e\%%L\"");
    gnuplot_cmd(ht[4], "set logscale y");
    gnuplot_cmd(ht[4], "set format y \"1e\%%L\"");
    gnuplot_cmd(ht[5], "set logscale y");
    gnuplot_cmd(ht[5], "set format y \"1e\%%L\"");
    //Grid set by default
    for(int i = 0; i <6; i++) gnuplot_cmd(ht[i],  "set grid");

    //Orbital error
    //----------------
    gnuplot_set_xlabel(ht[0], (char*)"t [-]");
    gnuplot_set_ylabel(ht[0], (char*)"eO [-]");
    gnuplot_plot_xy(ht[0], tc, eOc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "2", color);

    //XY
    //----------------
    gnuplot_set_xlabel(ht[1], (char*)"x [-]");
    gnuplot_set_ylabel(ht[1], (char*)"y [-]");
    gnuplot_plot_xy(ht[1], xc, yc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[1], xsc, ysc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "dashed", "2", color);

    gnuplot_plot_xy(ht[1], &xc[0], &yc[0], 1, "", "points", "1", "2", color);

    //XZ
    //----------------
    gnuplot_set_xlabel(ht[2], (char*)"x [-]");
    gnuplot_set_ylabel(ht[2], (char*)"z [-]");
    gnuplot_plot_xy(ht[2], xc, zc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[2], xsc, zsc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "dashed", "2", color);

    //YZ
    //----------------
    gnuplot_set_xlabel(ht[3], (char*)"y [-]");
    gnuplot_set_ylabel(ht[3], (char*)"z [-]");
    gnuplot_plot_xy(ht[3], yc, zc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[3], ysc, zsc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "dashed", "2", color);


    //Invariance error
    //----------------
    gnuplot_set_xlabel(ht[4], (char*)"t [-]");
    gnuplot_set_ylabel(ht[4], (char*)"eI [-]");
    gnuplot_plot_xy(ht[4], tc, eIc, Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "2", color);


    //Hamiltonian error
    //----------------
    gnuplot_set_xlabel(ht[5], (char*)"t [-]");
    gnuplot_set_ylabel(ht[5], (char*)"eH [-]");
    gnuplot_plot_xy(ht[5], tc, eHc,  Npoints+1, (char*)("Order "+ss1).c_str(), "lines", "1", "3", color);


    //Titles
    //----------------
    gnuplot_cmd(ht[0],  ("set title \"Error in the orbit, |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[1], ("set title \"Orbit in XY plane,  |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[2], ("set title \"Orbit in XZ plane,  |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[3], ("set title \"Orbit in YZ plane,  |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[4],  ("set title \"Error in the invariance equation, |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[5],  ("set title \"Error in the hamiltonian, |W(s,0)| = "+ssW+"\" ").c_str());

    //Additional names
    string eXYstr = (F_PLOT+"orbits/"+"XY_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eXZstr = (F_PLOT+"orbits/"+"XZ_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eYZstr = (F_PLOT+"orbits/"+"YZ_"+"Order_"+ss1+"_Size_"+st0s+".txt");


//    //Save in EPS format
//    gnuplot_cmd(ht[0],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[0], ("set output \""+eOstr+".eps\"").c_str());
//    gnuplot_cmd(ht[0], "replot");
//
//    gnuplot_cmd(ht[1],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[1], ("set output \""+eXYstr+".eps\"").c_str());
//    gnuplot_cmd(ht[1], "replot");
//
//    gnuplot_cmd(ht[2],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[2], ("set output \""+eXZstr+".eps\"").c_str());
//    gnuplot_cmd(ht[2], "replot");
//
//    gnuplot_cmd(ht[3],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[3], ("set output \""+eYZstr+".eps\"").c_str());
//    gnuplot_cmd(ht[3], "replot");
//
//    gnuplot_cmd(ht[4],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[4], ("set output \""+eIstr+".eps\"").c_str());
//    gnuplot_cmd(ht[4], "replot");
//
//    gnuplot_cmd(ht[5],"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
//    gnuplot_cmd(ht[5], ("set output \""+eIstr+".eps\"").c_str());
//    gnuplot_cmd(ht[5], "replot");

    return GSL_SUCCESS;
}

/**
 *  \brief Computations of various the orbial error eO along an orbit initialized by the pm of the center manifold W(s,t)
 **/
int eOPlot(const double st0[],         //RCM initial conditions
              vector<Oftsc>& Wh,       //TFC manifold
              matrix<Ofsc>& PC,        //COC matrix: z = PC*zh+V
              vector<Ofsc>& V,         //COC vector: z = PC*zh+V
              double t1,               //Final integration tim
              gsl_odeiv2_driver *dnc,  //driver for NC integration
              gsl_odeiv2_driver *drvf, //drive for RVF integration (reduced vector field)
              int Npoints,             //Number of points on which the errors are estimated
              FBPL& fbpl,          //current QBCP
              int order,               //Order for the eval of the OFTS objects
              int ofs_order,           //Order for the eval of the OFS objects
              gnuplot_ctrl  **ht,      //Gnuplot handlers
              int color)               //Color of the plots
{
    //------------------------------------------------------------------------------------
    // Strings
    //------------------------------------------------------------------------------------
    string F_GS    = fbpl.cs.F_GS;
    string F_PLOT  = fbpl.cs.F_PLOT;
    string F_COC   = fbpl.cs.F_COC;


    //------------------------------------------------------------------------------------
    //Reset integrator
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(dnc);
    gsl_odeiv2_driver_reset(drvf);

    //------------------------------------------------------------------------------------
    //Variables for plotting
    //------------------------------------------------------------------------------------
    //Time
    double tc[Npoints+1];
    //Orbits
    double xc[Npoints+1], yc[Npoints+1], zc[Npoints+1];
    //Orbit in TFC
    double xsc[Npoints+1], ysc[Npoints+1], zsc[Npoints+1];
    //Errors
    double eOc[Npoints+1];


    //------------------------------------------------------------------------------------
    //Errors
    //------------------------------------------------------------------------------------
    double eOm, eO[6];

    //------------------------------------------------------------------------------------
    // Inner state (NC, TFC...)
    //------------------------------------------------------------------------------------
    double  s1rcm[REDUCED_NV];  //RCM
    cdouble s1ccm[REDUCED_NV];  //CCM
    double  s1ccm8[2*REDUCED_NV]; //CCM8
    double  z1ncr[6];  //NC, integrated with NC vector field
    double  z1ncr2[6]; //NC, but computed from z = W(s,t)
    double  fnc[6];  //NC vector field
    double  W0;         //Euclidian norm in NC coordinates
    double  W0km;       //Euclidian norm in km
    Ofsc AUX;       //OFS temp variable

    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    //------------------------------------------------------------------------------------
    double n = fbpl.us.n;


    //------------------------------------------------------------------------------------
    //Initial time
    //------------------------------------------------------------------------------------
    double t0 = 0.0;//pij(2,5);

    //------------------------------------------------------------------------------------
    // RCM to NC for NC initial conditions
    //------------------------------------------------------------------------------------
    RCMtoNCbyTFC(st0, t0, n, order, ofs_order, Wh, PC, V, z1ncr, false);

    //------------------------------------------------------------------------------------
    // RCM to CCM8 for CCM initial conditions
    //------------------------------------------------------------------------------------
    RCMtoCCM8(st0, s1ccm8);

    //------------------------------------------------------------------------------------
    // Euclidian Norm
    //------------------------------------------------------------------------------------
    W0   = ENorm(z1ncr, 3);
    W0km = W0*fbpl.cs.gamma*fbpl.cs.cr3bp.L;

    //------------------------------------------------------------------------------------
    // Print Initial Conditions
    //------------------------------------------------------------------------------------
    //Initial conditions
    cout << "With current s0 and order, the initial conditions are: " << endl;
    for(int p = 0; p < 6; p++) cout << p << "  " << z1ncr[p] << endl;
    //Euclidian norm in km
    cout << "Approximated distance from Li [km]: " << W0km << endl;

    //------------------------------------------------------------------------------------
    // For plotting (first value)
    //------------------------------------------------------------------------------------
    //Time
    tc[0] = t0;
    //NC
    xc[0] = z1ncr[0];
    yc[0] = z1ncr[1];
    zc[0] = z1ncr[2];
    //NC from manifold
    xsc[0] = z1ncr[0];
    ysc[0] = z1ncr[1];
    zsc[0] = z1ncr[2];


    //------------------------------------------------------------------------------------
    // Loop on time
    //------------------------------------------------------------------------------------
    double t = t0;
    double t2 = t0;
    double ti = t0;
    if(t1 <0)
    {
        dnc->h = -dnc->h;
        drvf->h = -drvf->h;
    }

    for(int i =1; i<= Npoints; i++)
    {
        //------------------------------------------------------------------------------------
        // Apply driver
        //------------------------------------------------------------------------------------
        ti = t0 + (double) i * t1 / Npoints;
        gsl_odeiv2_driver_apply (dnc, &t, ti, z1ncr);
        gsl_odeiv2_driver_apply (drvf, &t2, ti, s1ccm8);

        //------------------------------------------------------------------------------------
        // Comparison
        //------------------------------------------------------------------------------------
        //---------------------
        // Using the PM: z = W(s,t)
        // Computed in z1ncr2
        //---------------------
        //CCM8 to RCM
        CCM8toRCM(s1ccm8, s1rcm);
        //CCM8 to CCM
        CCM8toCCM(s1ccm8, s1ccm);
        //z1ncr2 = W(s1rcm, ti)
        RCMtoNCbyTFC(s1rcm, ti, n, order, ofs_order, Wh, PC, V, z1ncr2, false);
        //Evaluating the vector field at z1ncr2 = W(s1rcm, ti)
        qbfbp_vfn_novar(ti, z1ncr2, fnc, dnc->sys->params);

        //---------------------
        //Error computation
        //---------------------
        for(int p = 0; p < Csts::NV; p++)
        {
            //eO = |z(t) - W(s(t), t)|
            eO[p] = cabs(z1ncr[p] - z1ncr2[p]);
        }


        //Taking the maximum (infinity norm)
        eOm = eO[0];
        for(int p = 1; p < Csts::NV; p++) if(eO[p] > eOm) eOm = eO[p];

        //------------------------------------------------------------------------------------
        //Store for plotting
        //------------------------------------------------------------------------------------
        //Timde
        tc[i]  = ti;
        //NC
        xc[i]  = z1ncr[0];
        yc[i]  = z1ncr[1];
        zc[i]  = z1ncr[2];
        //NC from manifold
        xsc[i] = z1ncr2[0];
        ysc[i] = z1ncr2[1];
        zsc[i] = z1ncr2[2];
        //Errors
        eOc[i] = eOm;
    }
    if(t1 <0)
    {
        dnc->h = -dnc->h;
        drvf->h = -drvf->h;
    }


    cout << "End angle is : " << t*n << endl;


    //------------------------------------------------------------------------------------
    //Define the output name
    //------------------------------------------------------------------------------------
    ifstream readStream;
    string ss1, ssW, st0s, ssofs;
    ss1   = static_cast<ostringstream*>( &(ostringstream() << order) )->str();
    ssofs = static_cast<ostringstream*>( &(ostringstream() << ofs_order) )->str();
    ssW   = static_cast<ostringstream*>( &(ostringstream() << W0) )->str();
    st0s  = static_cast<ostringstream*>( &(ostringstream() << creal(st0[0]+st0[1]+st0[2]+st0[3])) )->str();
    st0s  = st0s+"_test";


    //------------------------------------------------------------------------------------
    //In txt files
    //------------------------------------------------------------------------------------
    string eOstr   = (F_PLOT+"orbits/"+"eO_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string eXYZstr = (F_PLOT+"orbits/"+"XYZ_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    string ePMstr  = (F_PLOT+"orbits/"+"XYZ_PM_"+"Order_"+ss1+"_Size_"+st0s+".txt");
    gnuplot_fplot_xy(tc, eOc, Npoints+1, eOstr.c_str());                //eO
    gnuplot_fplot_txyz(tc, xc, yc,   zc,   Npoints, eXYZstr.c_str());   //orbit
    gnuplot_fplot_txyz(tc, xsc, ysc, zsc,  Npoints, ePMstr.c_str());    //PM orbit


    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    //Logscale if necessary (all errors)
    gnuplot_cmd(ht[0], "set logscale y");
    gnuplot_cmd(ht[0], "set format y \"1e\%%L\"");
    gnuplot_cmd(ht[4], "set logscale y");
    gnuplot_cmd(ht[4], "set format y \"1e\%%L\"");
    gnuplot_cmd(ht[5], "set logscale y");
    gnuplot_cmd(ht[5], "set format y \"1e\%%L\"");
    //Grid set by default
    for(int i = 0; i <6; i++) gnuplot_cmd(ht[i],  "set grid");

    //Orbital error
    //----------------
    gnuplot_set_xlabel(ht[0], (char*)"t [-]");
    gnuplot_set_ylabel(ht[0], (char*)"eO [-]");
    gnuplot_plot_xy(ht[0], tc, eOc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "1", "2", color);

    //XY
    //----------------
    gnuplot_set_xlabel(ht[1], (char*)"x [-]");
    gnuplot_set_ylabel(ht[1], (char*)"y [-]");
    gnuplot_plot_xy(ht[1], xc, yc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[1], xsc, ysc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "dashed", "2", color);

    //XZ
    //----------------
    gnuplot_set_xlabel(ht[2], (char*)"x [-]");
    gnuplot_set_ylabel(ht[2], (char*)"z [-]");
    gnuplot_plot_xy(ht[2], xc, zc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[2], xsc, zsc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "dashed", "2", color);

    //YZ
    //----------------
    gnuplot_set_xlabel(ht[3], (char*)"y [-]");
    gnuplot_set_ylabel(ht[3], (char*)"z [-]");
    gnuplot_plot_xy(ht[3], yc, zc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "1", "2", color);
    gnuplot_plot_xy(ht[3], ysc, zsc, Npoints+1, (char*)("Order "+ss1+ " OFS "+ssofs).c_str(), "lines", "dashed", "2", color);


    //Titles
    //----------------
    gnuplot_cmd(ht[0], ("set title \"Error in the orbit, |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[1], ("set title \"Orbit in XY plane,  |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[2], ("set title \"Orbit in XZ plane,  |W(s,0)| = "+ssW+"\" ").c_str());
    gnuplot_cmd(ht[3], ("set title \"Orbit in YZ plane,  |W(s,0)| = "+ssW+"\" ").c_str());


    return GSL_SUCCESS;
}



