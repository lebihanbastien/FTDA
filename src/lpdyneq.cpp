#include "lpdyneq.h"

/**
 * \file lpdyneq.cpp
 * \brief Computation of the dynamical equivalent to the Libration points for various models.
 * \author BLB.
 * \date August 2016
 * \version 1.1
 */



//========================================================================================
//         Dynamical equivalents to the Libration points
//========================================================================================
/**
 *  \brief Computation of the dynamical equivalents of the libration points. Test function.
 *
 *         - This function makes use of the subroutine lpdyneq that is the true heart of the computation.
 *         In particular, it contains the initialization of the first guess
 *          (the geometrical positions of the CRTBP libration points)
 *          + the differential corrector.
 *         - After lpdyneq, the resulting initial conditions are integrated and plotted on a full orbit.
 *         - Finally, a test of the periodicity of the orbit is performed, via the computation of the error |y(0) - y(T)|.
 *         - If is_stored is true, the results (x(t), y(t)) in synodical coordinates are stored in txt files of the form: "./plot/QBCP/DYNEQ/DYNEQ_QBCP_EM_L1.txt".
 **/
void compute_dyn_eq_lib_point(FBPL &fbpl, int is_stored)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    //------------------------------------------------------------------------------------
    // Integration tools
    // Initial vector field include nonlinear variationnal equations
    // to get the periodic orbit
    //------------------------------------------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function      = fbpl.is_norm? qbfbp_vfn_varnonlin:qbfbp_vf;
    sys.jacobian      = NULL;
    sys.dimension     = 42;
    sys.params        = &fbpl;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
                       Config::configManager().G_PREC_HSTART(),
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL());

    //Building the filename
    string f_model  = init_F_MODEL(fbpl.model);
    string f_csys   = init_F_COORDSYS(fbpl.coord_sys);
    string f_li     = init_F_LI(fbpl.li);
    string f_norm   = fbpl.is_norm?"_NC":"";
    string filename = "plot/"+f_model+"/DYNEQ/DYNEQ_"+f_model+"_"+f_csys+"_"+f_li+f_norm+".txt";

    //====================================================================================
    // 2. Display
    //====================================================================================
    cout << endl;
    cout << "===============================================" << endl;
    cout << "Computation of a dynamical equivalent of a     " << endl;
    cout << "libration point on a full period of the system." << endl;
    cout << "Characteristics:                               " << endl;
    cout << " - Model:  " << f_model << endl;
    cout << " - Number: " << f_li << endl;
    cout << " - Coordinates: " << f_norm+f_csys << endl;
    cout << " - Saved in: " << filename << endl;
    cout << "===============================================" << endl;

    //====================================================================================
    // 3. Single shooting to get the dynamical equivalent of the libration point (lpdyneq)
    //====================================================================================
    cout << "First, a single shooting procedure is used.    " << endl;
    cout << "===============================================" << endl;
    double y0[42];
    lpdyneq(d, y0);

    //====================================================================================
    // 4. Plot & Print in file if necessary
    //====================================================================================
    int Npoints = 500;

    //Building the solution and printing in file
    odePlot2(y0, 42, fbpl.us.T, d, h1, Npoints, 6, fbpl.is_norm, is_stored, "From single shooting", filename.c_str());

    //====================================================================================
    // 5. Periodicity condition for single shooting
    //====================================================================================
    periodicity_condition(y0, 42, 6, fbpl.us.T, d, fbpl.is_norm);

    //====================================================================================
    // 6. Multiple shooting
    //====================================================================================
    int M = 16;
    cout << "Second, a multiple shooting procedure is used,  " << endl;
    cout << "with " << 16 << " patch points.                 " << endl;
    cout << " Note that the multiple shooting procedure      " << endl;
    cout << " is used here only for the sake of completeness." << endl;
    cout << " In particular, it is not used in nfo2.cpp      " << endl;
    cout << " since it did not improve the results.          " << endl;
    cout << "=============================================== " << endl;
    double **ymd  = dmatrix(0, 41, 0, M);
    double *tmd   = dvector(0, M);

    //Building the patch points
    lpdyneq_patch_points(y0, d, fbpl.us.T, ymd, tmd, M);

    //Multiple shooting procedure
    lpdyneq_multiple_shooting(ymd, tmd, ymd, tmd, d, 42, M, 1e-14, false, h1, false);

    //Plot
    odePlotvec(ymd, tmd, 42, M, d, h1, Npoints, 3, true, "From multiple shooting");

    //====================================================================================
    // 7. Close window
    //====================================================================================
    pressEnter(true, "Press ENTER to close the gnuplot window(s)\n");
    gnuplot_close(h1);
}

/**
 *  \brief Main routine for the computation of the dynamical equivalent to the Libration points.
 *         The results are given in the form of the initial conditions (42 dimensions) updated in the array y0[].
 *
 *          €€TODO: BCP case is a work in progress. works only for EML1, in NC coordinates.
 **/
void lpdyneq(gsl_odeiv2_driver *d, double y0[])
{
    //====================================================================================
    // 0. Misc init.
    //====================================================================================
    //Settings for iostd
    Config::configManager().coutlp();

    //Retrieving the parameters
    FBPL* qbp      = (FBPL *) d->sys->params;
    int is_norm =  qbp->is_norm;
    int li           =  qbp->cs.li;
    double tend      =  qbp->us.T;

    //====================================================================================
    // 1. Initial guess
    //====================================================================================
    //Position & Velocity
    double yv[6];

    //------------------------------------------------------------------------------------
    //Approximated position from CR3BP
    //------------------------------------------------------------------------------------
    if(is_norm)
    {
        //Approximated position from CR3BP: null vector
        for(int i =0; i<6; i++) y0[i] = 0.0;

        //€€TODO: BCP case is a work in progress. works only for EML1
        if(qbp->model == Csts::BCP && qbp->coord_sys == Csts::EM && qbp->li_EM == 1)
        {
            //Forced IC
            y0[0] = +4.584016186756542e-03;
            y0[1] = 0.0;
            y0[2] = 0.0;
            y0[3] = 0.0;
            y0[4] = -6.163235628198283e-02;
            y0[5] = 0.0;
            cout << "Force IC in QBP, EML1 case" << endl;
        }
    }
    else
    {
        //Approximated position from CR3BP
        for(int i =0; i<6; i++) yv[i] = 0.0;
        switch(li)
        {
        case 1:
            yv[0] = -SEML.cs.cr3bp.l1.position[0];
            break;
        case 2:
            yv[0] = -SEML.cs.cr3bp.l2.position[0];
            break;
        case 3:
            yv[0] = -SEML.cs.cr3bp.l3.position[0];
            break;
        }

        switch(SEML.coord_sys)
        {
        case Csts::EM:
        {
            //To momenta
            EMvtoEMm(0.0, yv, y0, qbp);
            break;
        }

        case Csts::SEM:
        {
            //To momenta
            SEMvtoSEMm(0.0, yv, y0, qbp);
            break;
        }

        }
    }

    //------------------------------------------------------------------------------------
    // STM is appened to y0 (STM = eye(6) at t = t0)
    //------------------------------------------------------------------------------------
    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_matrixToVector(y0, Id, 6, 6, 6);

    //====================================================================================
    // 2. Differential correction scheme
    //====================================================================================
    //Initial IC
    cout << "lpdyneq. Initial state:" << endl;
    for(int i =0; i < 6; i++) cout << y0[i] << endl;

    double prec = 2e-14;
    cout << "lpdyneq. Starting differential correction..." << endl;
    differential_correction(y0, 0.5*tend, prec, d, 42, 0);
    cout << "lpdyneq. End of differential correction." << endl;

    //====================================================================================
    // 3. Print result
    //====================================================================================
    //Final IC
    cout << "lpdyneq. Final state:" << endl;
    for(int i =0; i<6; i++) cout << y0[i] << endl;
}

//========================================================================================
// Dynamical equivalents via continuation procedures
//========================================================================================
/**
 *  \brief Computation of the dynamical equivalents of the libration points via continuation between two models. WORK IN PROGRESS.
 *         This routine make use of a pseudo-arclength continuation to perform a continuation procedure from the dynamical equivalent
 *         computed in a model_1 to its counterpart in a model_2. The models are nested in the driver *d.
 *
 *         Note that the pseudo-arclength additionnal constraint is not imposed. In place, the minimum norm solution is used.
 **/
void lpdyneq_cont_2(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[])
{
    //============================================================
    // 1. Initialization
    //============================================================
    //Settings for iostd
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------------
    char ch;                //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    //h1 contains the solution in the XY plane.
    //h2 contains the evolution of the X dimension wrt to epsilon
    //gnuplot_cmd(h2, "set xrange [-1.18:-1.14]");
    //gnuplot_cmd(h2, "set yrange [-0.1:0.1]");
    gnuplot_set_xlabel(h2, (char*) "X [synodical coordinates]");
    gnuplot_set_ylabel(h2, (char*) "Epsilon (continuation)");

    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    //------------------------------------------------------------------------------------
    QBCP_I* qbpi     =  (QBCP_I *) d->sys->params;
    FBPL* qbp1     =  &qbpi->model1;              //First model, full model when qpbi->epsilon = 0.0, at the beginning of the process.
    FBPL* qbp2     =  &qbpi->model2;              //First model, full model when qpbi->epsilon = 1.0, at the end of the process.
    int is_norm =  qbp1->is_norm;
    int li           =  qbp1->cs.li;
    //Full period of the Sun
    double tend = 2*M_PI/qbp2->us.n;

    //------------------------------------------------------------------------------------
    // Evaluate the alphas
    //------------------------------------------------------------------------------------
    double alphaf[8];
    eval_array_coef(alphaf, 0.0, qbp1->us.n,  qbp1->n_order_fourier, qbp1->cs.coeffs, 8);

    //============================================================
    // 2. Initialization
    //============================================================
    //CR3BPs
    CR3BP EM, SE;
    init_cr3bp(&EM, Csts::EARTH, Csts::MOON);
    init_cr3bp(&SE, Csts::SUN, Csts::EARTH_AND_MOON);

    //------------------------------------------------------------------------------------
    //Period that we seek: To = fT*Ts
    //------------------------------------------------------------------------------------
    double fT = 0.5;

    //------------------------------------------------------------------------------------
    // Initial conditions
    //------------------------------------------------------------------------------------
    //Position & Velocity
    double yv[6];
    double yf[6];
    for(int i =0; i<6; i++) yv[i] = 0.0;

    //Switch on the coordinate system
    switch(qbp1->coord_sys)
    {
    case Csts::EM:
        switch(li)
        {
        case 1:
            yv[0] = -EM.l1.position[0];
            break;
        case 2:
            yv[0] = -EM.l2.position[0];
            break;
        case 3:
            yv[0] = -EM.l3.position[0];
            break;
        }
        //Position & Momenta
        EMvtoEMm(0.0, yv, yf, qbp1);

        break;
    case Csts::SEM:

        switch(li)
        {
        case 1:
            yv[0] = -SE.l1.position[0];
            break;
        case 2:
            yv[0] = -SE.l2.position[0];
            break;
        case 3:
            yv[0] = -SE.l3.position[0];
            break;
        }
        //Position & Momenta
        SEMvtoSEMm(0.0, yv, yf, qbp1);
        break;
    }


    //------------------------------------------------------------------------------------
    //Transform into y0 if NC coordinates are used
    //------------------------------------------------------------------------------------
    if(is_norm)
    {
        //To NC
        SYStoNC(0.0, yf, y0, qbp1);
    }
    else
    {
        //Simple copy
        for(int i =0; i<6; i++) y0[i] = yf[i];
    }

    //------------------------------------------------------------------------------------
    // 36 variables for variational equations
    //------------------------------------------------------------------------------------
    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_matrixToVector(y0, Id, 6, 6, 6);

    //------------------------------------------------------------------------------------
    //The last 6 components of y0 contain
    // the variational equations for epsilon
    //------------------------------------------------------------------------------------
    for(int i = 42; i< 48; i++) y0[i] = 0.0;

    //============================================================
    // 3. Differential correction scheme
    //============================================================
    //---------------------------------------------------------
    //Arclength stepsize
    //---------------------------------------------------------
    double ds = 8e-3;
    int iterc = 1;

    //---------------------------------------------------------
    //Driver: loose control
    //---------------------------------------------------------
    gsl_odeiv2_control_set_driver (loose_control, d);

    //---------------------------------------------------------
    // First step: solution in the QBCP
    //---------------------------------------------------------
    int Npoints = 500;

    //Diff corr
    cout << "lpdyneq_cont. First step." << endl;
    cout << "Starting diff. corr. with epsilon = " << qbpi->epsilon << endl;
    differential_correction(y0, fT*tend, 1e-14, d, 48, 0);

    //Plotting the refined solution
    odePlot(y0, 48,  0.5*tend, d, h1, Npoints, 8);
    odePlot(y0, 48, -0.5*tend, d, h1, Npoints, 8);

    //---------------------------------------------------------
    // Second step: solution in the QBCP,
    // with epsilon as a parameter
    //---------------------------------------------------------
    cout << "lpdyneq_cont. Second step." << endl;
    double nullvector[7], gv[7];

    //differential_correction_deps(y0, tend, 1e-10, d, 48, 1); //to force the pseudo-arclength constraint;
    differential_correction_deps_mns(y0, nullvector, gv, fT*tend, 1e-10, d, 48, 0, 1);

    //---------------------------------------------------------
    //Plotting the refined solution
    //---------------------------------------------------------
    odePlot(y0, 48,  0.5*fT*tend, d, h1, Npoints, 7);
    odePlot(y0, 48, -0.5*fT*tend, d, h1, Npoints, 7);

    //---------------------------------------------------------
    //Plot in the (x, eps) plane
    //---------------------------------------------------------
    double factor = 1.0;
    double y0n[6], ymn[6], ym[6];

    //Direction of the null vector
    for(int i = 0; i <6; i++) y0n[i] = y0[i];
    y0n[0] = y0[0] + ds*nullvector[0];
    y0n[4] = y0[4] + ds*nullvector[1];
    NCtoSYS(0.0, y0, ym, qbp1);
    NCtoSYS(0.0, y0n, ymn, qbp1);
    string cmd = "set arrow from "+ num_to_string(ym[0]) +","+num_to_string(qbpi->epsilon);
    cmd +=" to " + num_to_string(factor*ymn[0])+","+num_to_string(factor*(qbpi->epsilon +ds*nullvector[2]));
    gnuplot_cmd(h2, cmd.c_str());
    gnuplot_plot_xy(h2, ym, &(qbpi->epsilon), 1, "", "points", "7", "3", 1);

    cout << "New eps = " << qbpi->epsilon << endl;
    printf("Press ENTER to continue\n");
    scanf("%c",&ch);

    //---------------------------------------------------------
    // Third step: solution in the QBCP,
    // with epsilon as a parameter + continuation
    //---------------------------------------------------------
    //Loop
    int Nmax = 3000;
    int iter = 0;
    qbpi->epsilon = 0.0;
    do
    {
        //------------------------
        //  New state, with a correction from null vector
        //------------------------
        y0[0] += ds*nullvector[0];
        y0[4] += ds*nullvector[1];
        qbpi->epsilon += ds*nullvector[2];

        //------------------------
        //  Correct the state with a differential correction procedure
        //------------------------
        // With a square system (pseudo-arclength constraint added)
        //differential_correction_deps_pac(y0, nullvector, gv, ds, fT*tend, 1e-10, d, 48, 1);
        // With the minimum norm solution
        iterc = differential_correction_deps_mns(y0, nullvector, gv, fT*tend, 1e-10, d, 48, (iter%50 == 0), 0);

        //------------------------
        // Plotting the direction of motion
        // in the (x, eps) space
        //------------------------
        if(iter%20 == 0)
        {
            //Direction of the null vector
            for(int i = 0; i <6; i++) y0n[i] = y0[i];
            y0n[0] = y0[0] + ds*nullvector[0];
            y0n[4] = y0[4] + ds*nullvector[1];
            NCtoSYS(0.0, y0, ym, qbp1);
            NCtoSYS(0.0, y0n, ymn, qbp1);
            cmd = "set arrow from "+ num_to_string(ym[0]) +","+num_to_string(qbpi->epsilon)+" to " + num_to_string(factor*ymn[0])+","+num_to_string(factor*(qbpi->epsilon +ds*nullvector[2]));
            gnuplot_cmd(h2, cmd.c_str());
            gnuplot_plot_xy(h2, ym, &(qbpi->epsilon), 1, "", "points", "7", "3", 1);
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);
        }

        //------------------------
        //Plotting the refined solution
        //------------------------
        if(iter%50 == 0)
        {
            NCtoSYS(0.0, y0, ym, &SEML);
            cout << "New state is = " << endl;
            vector_printf(ym, 6);
            cout << "New eps = " << qbpi->epsilon << endl;
            odePlot(y0, 48,  fT*tend, d, h1, Npoints, 8);
            odePlot(y0, 48, -0.5*fT*tend, d, h1, Npoints, 7);
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);
        }
        //Upload iter
        iter++;

        cout << "Old ds = " << ds << endl;
        ds = ds*6/(fabs(iterc)+1);
        cout << "New ds = " << ds << endl;

    }
    while(iter <= Nmax && qbpi->epsilon <= 1.0001);

    //---------------------------------------------------------
    // Fourth step: Last correction much more precise
    //---------------------------------------------------------
    //Set hard control
    gsl_odeiv2_control_set_driver (hard_control, d);

    //Last correction much more precise
    qbpi->epsilon = 1.00;
    if(is_norm) differential_correction(y0, fT*tend, 1e-14, d, 48, 0);
    else differential_correction(y0, fT*tend, 1e-14, d, 48, 0);

    //Final IC
    NCtoSYS(0.0, y0, ymn, qbp1);
    cout << "lpdyneq. IC(EM):" << endl;
    vector_printf(ymn, 6);

    //---------------------------------------------------------
    //Print in file
    //---------------------------------------------------------
    odePlotprint(y0, 48,  fT*tend, d, h1, Npoints, 6, "plot/lpdyneq_from_continuation.txt");

    printf("Press ENTER to continue\n");
    scanf("%c",&ch);
}

/**
 *  \brief Continuation routine for the dynamical equivalents of the libration points. Test function.
 *
 *  This routine makes use of the subroutine lpdyneq_cont_2 to perform a continuation procedure from the dynamical equivalent
 *         computed in a from_model to its counterpart in a to_model.
 */
void continuation_dyn_eq_lib_point(FBPL &fbpl, int from_model, int to_model)
{
    //============================================================
    // 1. Initialization
    //============================================================
    //Model must be normalized.
    if(!fbpl.is_norm)
    {
        cout << "continuation_dyn_eq_lib_point. Warning: selected model is not normalized. Premature ending." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Structures for continuation
    //------------------------------------------------------------------------------------
    QBCP_I model;
    FBPL model1, model2;
    init_fbp_cont(&model, &model1, &model2, Csts::SUN, Csts::EARTH, Csts::MOON, fbpl.is_norm, SEML.li_EM, SEML.li_SEM, 0, from_model, to_model, SEML.coord_sys, SEML.param_style);

    //------------------------------------------------------------------------------------
    // Integration tools
    // Initial vector field include nonlinear variationnal equations
    // to get the periodic orbit
    //------------------------------------------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function      = qbfbp_vfn_cont;
    sys.jacobian      = NULL;
    sys.dimension     = 42+6;
    sys.params        = &model;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
                       Config::configManager().G_PREC_HSTART(),
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL());

    //Differential correction to get the dynamical equivalent of the libration point (lpdyneq)
    gsl_odeiv2_control * loose_control = gsl_odeiv2_control_y_new(1e-10 , 1e-10);
    gsl_odeiv2_control *  hard_control = gsl_odeiv2_control_y_new(1e-16 , 1e-16);

    //============================================================
    // 2. Continuation procedure
    //============================================================
    double y0c[42+6];
    lpdyneq_cont_2(d, loose_control, hard_control, y0c);
}

/**
 *  \brief Computation of lpdyneq with continuation between two models. DEPRECATED
 **/
void lpdyneq_cont(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[], gnuplot_ctrl *h1)
{
    //Settings for iostd
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //Init the plot
    h1 = gnuplot_init();

    //Retrieving the parameters
    QBCP_I* qbpi     =  (QBCP_I *) d->sys->params;
    FBPL* qbp1     =  &qbpi->model1;              //First model, full model when qpbi->epsilon = 0.0, at the beginning of the process.
    FBPL* qbp2     =  &qbpi->model2;              //First model, full model when qpbi->epsilon = 1.0, at the end of the process.
    int is_norm =  qbp1->is_norm;
    int li           =  qbp1->cs.li;

    //Full period of the Sun
    double tend = 2*M_PI/qbp2->us_em.n;

    //----------------------------------------------------------------------------------------------------------
    // Evaluate the alphas
    //----------------------------------------------------------------------------------------------------------
    double alphaf[8];
    eval_array_coef(alphaf, 0.0, qbp1->us.n,  qbp1->n_order_fourier, qbp1->cs.coeffs, 8);

    //----------------------------------------------------------------------------------------------------------
    // Initial conditions
    //----------------------------------------------------------------------------------------------------------
    //CR3BPs
    CR3BP EM;
    init_cr3bp(&EM, Csts::EARTH, Csts::MOON);
    //Position & Velocity
    double yv[6];
    double yf[6];
    for(int i =0; i<6; i++) yv[i] = 0.0;
    switch(li)
    {
    case 1:
        yv[0] = -EM.l1.position[0];
        break;
    case 2:
        yv[0] = -EM.l2.position[0];
        break;
    case 3:
        yv[0] = -EM.l3.position[0];
        break;
    }

    //Position & Momenta
    EMvtoEMm(0.0, yv, yf, qbp1);

    //Transform into y0 if necessary
    if(is_norm)
    {
        //To NC
        SYStoNC(0.0, yf, y0, qbp1);
    }
    else
    {
        //Simple copy
        for(int i =0; i<6; i++) y0[i] = yf[i];
    }

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_matrixToVector(y0, Id, 6, 6, 6);

    //The last 6 components of y0 contain the variational equations for epsilon
    for(int i = 42; i< 48; i++) y0[i] = 0.0;

    //----------------------------------------------------------------------------------------------------------
    // Differential correction scheme
    //----------------------------------------------------------------------------------------------------------
    //Inversion
    double y0s[48];
    double t;

    //Arclength parameters
    gsl_matrix *A = gsl_matrix_calloc(6,6);
    gsl_matrix *ALU = gsl_matrix_calloc(6,6);
    double Aj[7];
    int s;
    gsl_permutation * p6 = gsl_permutation_alloc (6);
    double f[7];
    double ds = 1e-6;

    //Driver: loose control
    gsl_odeiv2_control_set_driver (loose_control, d);


    //-----------------------------------
    //First step:
    //-----------------------------------
    int Npoints = 5000;
    //Diff corr
    cout << "lpdyneq_cont. Starting diff. corr. with epsilon = " << setprecision(15) << qbpi->epsilon << setprecision(15) <<  endl;
    differential_correction(y0, 0.5*tend, 1e-14, d, 48, 0);

    //Plotting the refined solution
    odePlot(y0, 48,  0.5*tend, d, h1, Npoints, 8);
    odePlot(y0, 48, -0.5*tend, d, h1, Npoints, 8);

    char ch;
    printf("Press ENTER to continue\n");
    scanf("%c",&ch);


    //Loop
    int iter    = 0;
    qbpi->epsilon = 0.0;
    do
    {
        //------------------------
        //At the end of the integration, y0[6,...41] contains Df/Dy, y0[42...47] contains Df/Depsilon
        //------------------------
        gsl_odeiv2_driver_reset(d);
        for(int i = 0; i < 48; i++) y0s[i] = y0[i];
        t = 0.0;
        gsl_odeiv2_driver_apply (d, &t, tend, y0s);


        //------------------------
        // Let DP be the 6x7 matrix:
        //  DP = | Df/dy Df/Depsilon|
        //
        // Then build the scalars Aj = determinant of the matrix obtained from DP erasing the j-th column times (-1)^(j+1);
        //------------------------
        for(int j = 0; j < 7; j++)
        {
            for(int k = 0; k < 7; k++)
            {
                if(k < j) for(int p = 0; p < 6; p++) gsl_matrix_set(A, p, k, y0s[6*(k+1)+p]);
                if(k > j) for(int p = 0; p < 6; p++) gsl_matrix_set(A, p, k-1, y0s[6*(k+1)+p]);
            }

            //Copy: ALU = A
            gsl_matrix_memcpy(ALU, A);
            //Lu decomposition
            gsl_linalg_LU_decomp (ALU, p6, &s);
            //Determinant Aj[j] = det(A)
            Aj[j] = gsl_linalg_LU_det(ALU, s);
        }

        //------------------------
        // Compute the arclength equations of motion
        //------------------------
        arclengthvf(f, Aj);

        cout << "Aj = " << endl;
        vector_printf(Aj, 7);

        cout << "f = " << endl;
        vector_printf(f, 7);

        cout << "State before euler step:" << endl;
        vector_printf(y0, 6);
        cout << "epsilon = " << qbpi->epsilon << endl;

        //------------------------
        //  Euler step
        //------------------------
        for(int i = 0; i < 6; i++) y0[i] = y0[i] + ds * f[i];
        qbpi->epsilon = qbpi->epsilon + ds*f[6];

        cout << "New state after euler step:" << endl;
        vector_printf(y0, 6);
        cout << "epsilon = " << qbpi->epsilon << endl;

        char ch;
        printf("Press ENTER to continue\n");
        scanf("%c",&ch);

        //------------------------
        //  Correct the state with a differential correction procedure
        //------------------------
        differential_correction_deps(y0, tend, 1e-10, d, 48, 1);
        //differential_correction(y0, 0.5*tend, 1e-14, d, 48, 1);


        //Upload iter
        iter++;


        printf("Press ENTER to close the gnuplot window(s)\n");
        scanf("%c",&ch);

    }
    while(/*iter <= Nmax &&*/ qbpi->epsilon <= 1.001);

    //Set hard control
    gsl_odeiv2_control_set_driver (hard_control, d);

    cout << "here, after loop." << endl;

    //Last correction much more precise
    if(is_norm) differential_correction(y0, 0.5*tend, 1e-14, d, 48, 0);
    else differential_correction(y0, 0.5*tend, 1e-14, d, 48, 0);


    //Final IC
    cout << "lpdyneq. IC(NC):" << endl;
    for(int i =0; i<6; i++) cout << y0[i] << endl;
}

//========================================================================================
// Continuation procedures for resonant orbits. Ex: 2T/5 resonnance at EML1
//========================================================================================
/**
 *  \brief Computation of resonant orbits via continuation between two models. WORK IN PROGRESS.
 *         This routine make use of a pseudo-arclength continuation to perform a continuation procedure from a resonant orbit
 *         computed in a model_1 to its counterpart in a model_2. The models are nested in the driver *d.
 *
 *         Note that the pseudo-arclength additionnal constraint is not imposed. In place, the minimum norm solution is used.
 *         The available orbits are:
 *                   - A 2T/5-resonant orbit at EML1
 *                   - A T/2-resonant orbit at EML2
 **/
void res_orbit_cont(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[])
{
    //============================================================
    // 1. Initialization
    //============================================================
    //Settings for iostd
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //------------------------------------------------------------------------------------
    //Plotting devices
    //------------------------------------------------------------------------------------
    char ch;                //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    //h1 contains the solution in the XY plane.
    //h2 contains the evolution of the X dimension wrt to epsilon
    //gnuplot_cmd(h2, "set xrange [-1.18:-1.14]");
    //gnuplot_cmd(h2, "set yrange [-0.1:0.1]");
    gnuplot_set_xlabel(h2, (char*) "X [synodical coordinates]");
    gnuplot_set_ylabel(h2, (char*) "Epsilon (continuation)");

    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    //------------------------------------------------------------------------------------
    QBCP_I* qbpi     =  (QBCP_I *) d->sys->params;
    FBPL* qbp1     =  &qbpi->model1;              //First model, full model when qpbi->epsilon = 0.0, at the beginning of the process.
    FBPL* qbp2     =  &qbpi->model2;              //First model, full model when qpbi->epsilon = 1.0, at the end of the process.
    int is_norm =  qbp1->is_norm;
    int li           =  qbp1->cs.li;
    //Full period of the Sun
    double tend = 2*M_PI/qbp2->us.n;

    //------------------------------------------------------------------------------------
    // Evaluate the alphas
    //------------------------------------------------------------------------------------
    double alphaf[8];
    eval_array_coef(alphaf, 0.0, qbp1->us.n,  qbp1->n_order_fourier, qbp1->cs.coeffs, 8);


    //============================================================
    // 2. Initialization
    //============================================================
    //------------------------------------------------------------------------------------
    // Initial conditions + period of the resonant orbit
    //------------------------------------------------------------------------------------
    double yv[6];
    double yf[6];
    for(int i =0; i<6; i++) yv[i] = 0.0;

    //---------------------------------------------------------
    //Period  of the resonant orbit that we seek: To = fT*Ts
    //---------------------------------------------------------
    double fT = 0.5; //by default, but can be changed in the subsequent switch

    //Switch on the coordinate system
    switch(qbp1->coord_sys)
    {
    case Csts::EM:
        //Switch on the libration point
        switch(li)
        {
        case 1:
            //2T/5-periodic lyapunov orbit of the CRTBP around EML1
            yv[0] = -0.827043252705635;
            yv[4] = -0.089404728400948;
            fT = 2.0/5;
            break;
        case 2:
            //T/2-periodic lyapunov orbit of the CRTBP around EML2 (022 in Andreu)
            yv[0] = -1.130364852389014;
            yv[4] = -0.128328493319532;
            fT = 0.5;
            break;
        }

        //Position & Momenta
        EMvtoEMm(0.0, yv, yf, qbp1);
        break;

    case Csts::SEM:
        cout << "res_orbit_cont. Csts::SEM coordinate system is currently not supported. end." << endl;
        return;
        break;
    }

    //------------------------------------------------------------------------------------
    //Transform into y0 if NC coordinates are used
    //------------------------------------------------------------------------------------
    if(is_norm)
    {
        //To NC
        SYStoNC(0.0, yf, y0, qbp1);
    }
    else
    {
        //Simple copy
        for(int i =0; i<6; i++) y0[i] = yf[i];
    }

    //------------------------------------------------------------------------------------
    // 36 variables for variational equations
    //------------------------------------------------------------------------------------
    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_matrixToVector(y0, Id, 6, 6, 6);

    //------------------------------------------------------------------------------------
    //The last 6 components of y0 contain
    // the variational equations for epsilon
    //------------------------------------------------------------------------------------
    for(int i = 42; i< 48; i++) y0[i] = 0.0;

    //============================================================
    // 3. Differential correction scheme
    //============================================================
    //---------------------------------------------------------
    //Arclength stepsize
    //---------------------------------------------------------
    double ds = 5e-3;

    //---------------------------------------------------------
    //Driver: loose control to fasten the computation
    //---------------------------------------------------------
    gsl_odeiv2_control_set_driver (loose_control, d);

    //---------------------------------------------------------
    // First step: solution in the QBCP
    //---------------------------------------------------------
    int Npoints = 500;

    //Diff corr
    cout << "res_orbit_cont. First step." << endl;
    cout << "Starting diff. corr. with epsilon = " << qbpi->epsilon << endl;
    differential_correction(y0, fT*tend, 1e-14, d, 48, 0);

    //Plotting the refined solution
    odePlot(y0, 48,  0.5*tend, d, h1, Npoints, 8);
    odePlot(y0, 48, -0.5*tend, d, h1, Npoints, 8);

    //---------------------------------------------------------
    // Second step: solution in the QBCP,
    // with epsilon as a parameter
    //---------------------------------------------------------
    cout << "res_orbit_cont. Second step." << endl;
    double nullvector[7], gv[7];

    //differential_correction_deps(y0, tend, 1e-10, d, 48, 1);
    differential_correction_deps_mns(y0, nullvector, gv, fT*tend, 1e-10, d, 48, 0, 1);

    //---------------------------------------------------------
    //Plotting the refined solution
    //---------------------------------------------------------
    odePlot(y0, 48,  0.5*fT*tend, d, h1, Npoints, 7);
    odePlot(y0, 48, -0.5*fT*tend, d, h1, Npoints, 7);

    //---------------------------------------------------------
    //Plot in the (x, eps) plane
    //---------------------------------------------------------
    double factor = 1.0;
    double y0n[6], ymn[6], ym[6];

    //Direction of the null vector
    for(int i = 0; i <6; i++) y0n[i] = y0[i];
    y0n[0] = y0[0] + ds*nullvector[0];
    y0n[4] = y0[4] + ds*nullvector[1];
    NCtoSYS(0.0, y0, ym, qbp1);
    NCtoSYS(0.0, y0n, ymn, qbp1);
    string cmd = "set arrow from "+ num_to_string(ym[0]) +","+num_to_string(qbpi->epsilon);
    cmd +=" to " + num_to_string(factor*ymn[0])+","+num_to_string(factor*(qbpi->epsilon +ds*nullvector[2]));
    gnuplot_cmd(h2, cmd.c_str());
    gnuplot_plot_xy(h2, ym, &(qbpi->epsilon), 1, "", "points", "7", "3", 1);


    cout << "New eps = " << qbpi->epsilon << endl;
    printf("Press ENTER to continue\n");
    scanf("%c",&ch);

    //---------------------------------------------------------
    // Third step: solution in the QBCP,
    // with epsilon as a parameter + continuation
    //---------------------------------------------------------
    //Loop
    int Nmax = 3000;
    int iter = 0;
    qbpi->epsilon = 0.0;
    do
    {
        //------------------------
        //  New state, with a correction from null vector
        //------------------------
        y0[0] += ds*nullvector[0];
        y0[4] += ds*nullvector[1];
        qbpi->epsilon += ds*nullvector[2];

        //------------------------
        //  Correct the state with a differential correction procedure
        //------------------------
        // With a square system (pseudo-arclength constraint added)
        //differential_correction_deps_pac(y0, nullvector, gv, ds, fT*tend, 1e-10, d, 48, 1);
        // With the minimum norm solution
        differential_correction_deps_mns(y0, nullvector, gv, fT*tend, 1e-10, d, 48, 0, 0);

        //------------------------
        // Plotting the direction of motion
        // in the (x, eps) space
        //------------------------
        if(iter%20 == 0)
        {
            //Direction of the null vector
            for(int i = 0; i <6; i++) y0n[i] = y0[i];
            y0n[0] = y0[0] + ds*nullvector[0];
            y0n[4] = y0[4] + ds*nullvector[1];
            NCtoSYS(0.0, y0, ym, qbp1);
            NCtoSYS(0.0, y0n, ymn, qbp1);
            cmd = "set arrow from "+ num_to_string(ym[0]) +","+num_to_string(qbpi->epsilon)+" to " + num_to_string(factor*ymn[0])+","+num_to_string(factor*(qbpi->epsilon +ds*nullvector[2]));
            gnuplot_cmd(h2, cmd.c_str());
            gnuplot_plot_xy(h2, ym, &(qbpi->epsilon), 1, "", "points", "7", "3", 1);
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);
        }

        //------------------------
        //Plotting the refined solution
        //------------------------
        if(iter%100 == 0)
        {
            NCtoSYS(0.0, y0, ym, &SEML);
            cout << "New state is = " << endl;
            vector_printf(ym, 6);
            cout << "New eps = " << qbpi->epsilon << endl;
            odePlot(y0, 48,  0.5*fT*tend, d, h1, Npoints, 8);
            odePlot(y0, 48, -0.5*fT*tend, d, h1, Npoints, 8);
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);
        }
        //Upload iter
        iter++;

    }
    while(iter <= Nmax && qbpi->epsilon <= 1.0001);

    //---------------------------------------------------------
    // Fourth step: Last correction much more precise
    //---------------------------------------------------------
    //Set hard control
    gsl_odeiv2_control_set_driver (hard_control, d);

    //Last correction much more precise
    qbpi->epsilon = 1.00;
    if(is_norm) differential_correction(y0, fT*tend, 1e-14, d, 48, 0);
    else differential_correction(y0, fT*tend, 1e-14, d, 48, 0);

    //Final IC
    NCtoSYS(0.0, y0, ymn, qbp1);
    cout << "lpdyneq. IC(EM):" << endl;
    vector_printf(ymn, 6);

    //---------------------------------------------------------
    //Print in file
    //---------------------------------------------------------
    switch(li)
    {
        case 1:
        odePlotprint(y0, 48,  fT*tend, d, h1, Npoints, 6, "plot/eml1_2T5_resonance.txt");
        break;

        case 2:
        odePlotprint(y0, 48,  fT*tend, d, h1, Npoints, 6, "plot/eml2_halfT_resonance.txt");
        break;
    }


    printf("Press ENTER to continue\n");
    scanf("%c",&ch);
}

/**
 *  \brief Continuation routine for the computation of resonant orbits. Test function.
 *
 *  This routine makes use of the subroutine res_orbit_cont to perform a continuation procedure from a resonant orbit
 *         computed in a model from_model to its counterpart in a model to_model.
 *         The available orbits are:
 *                   - A 2T/5-resonant orbit at EML1
 *                   - A T/2-resonant orbit at EML2
 */
void continuation_res_orbit(FBPL &fbpl, int from_model, int to_model)
{
    //============================================================
    // 1. Initialization
    //============================================================
    //------------------------------------------------------------------------------------
    //Model must be normalized.
    //------------------------------------------------------------------------------------
    if(!fbpl.is_norm)
    {
        cout << "continuation_res_orbit. Warning: selected model is not normalized. Premature ending." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Structures for continuation
    //------------------------------------------------------------------------------------
    QBCP_I model;
    FBPL model1, model2;
    init_fbp_cont(&model, &model1, &model2, Csts::SUN, Csts::EARTH, Csts::MOON, fbpl.is_norm, SEML.li_EM, SEML.li_SEM, 0, from_model, to_model, SEML.coord_sys, SEML.param_style);

    //------------------------------------------------------------------------------------
    // Integration tools
    // Initial vector field include nonlinear variationnal equations
    // to get the periodic orbit
    //------------------------------------------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function      = qbfbp_vfn_cont;
    sys.jacobian      = NULL;
    sys.dimension     = 42+6;
    sys.params        = &model;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, Config::configManager().G_PREC_HSTART(),
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL());

    //Differential correction to get the dynamical equivalent of the libration point (lpdyneq)
    gsl_odeiv2_control * loose_control = gsl_odeiv2_control_y_new(1e-10 , 1e-10);
    gsl_odeiv2_control *  hard_control = gsl_odeiv2_control_y_new(1e-16 , 1e-16);

    //============================================================
    // 2. Continuation procedure
    //============================================================
    double y0c[42+6];
    res_orbit_cont(d, loose_control, hard_control, y0c);
}

//====================================================================================
// Periodicity Condition
//====================================================================================
/**
 *  \brief Test the periodicity condition y(0) = y(t1), with initial condition y, for the vector field contained in the driver d.
 *
 *  Note that the integer N is the number of variables associated to the driver d, and should be also the number of variables in y. However, the periodicity condition is tested only on
 *  NvarTest variables (a usual example is Nvar = 42 but NvarTest = 6).
 **/
int periodicity_condition(const double y[], int Nvar, int NvarTest, double t1, gsl_odeiv2_driver *d, int is_norm)
{
    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(d);

    //Retrieving the parameters
    FBPL* qbp = (FBPL *) d->sys->params;

    //Initial conditions
    double ys[Nvar], y0[Nvar], y1[Nvar], y0_nat[Nvar], y1_nat[Nvar];
    for(int i=0; i< Nvar; i++) ys[i] = y[i];

    //First point in system coordinates
    if(is_norm) NCtoSYS(0.0, ys, y0, qbp);
    else for(int i=0; i< Nvar; i++) y0[i] = ys[i];

    //First point in native coordinates
    for(int i=0; i< Nvar; i++) y0_nat[i] = ys[i];

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    double t = 0.0;
    gsl_odeiv2_driver_apply (d, &t, t1, ys);

    //------------------------------------------------------------------------------------
    //Final state
    //------------------------------------------------------------------------------------
    // In system coordinates
    if(is_norm) NCtoSYS(t, ys, y1, qbp);
    else for(int i=0; i< Nvar; i++) y1[i] = ys[i];

    //In native coordinates
    for(int i=0; i< Nvar; i++) y1_nat[i] = ys[i];

    if(t1 <0) d->h = -d->h;

    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    // Norm |y(0) - y(T)|
    double yNorm = 0.0, yNorm_nat = 0.0;
    for(int i = 0; i < NvarTest; i++) yNorm += (y0[i] - y1[i])*(y0[i] - y1[i]);
    for(int i = 0; i < NvarTest; i++) yNorm_nat += (y0_nat[i] - y1_nat[i])*(y0_nat[i] - y1_nat[i]);
    yNorm = sqrt(yNorm);
    yNorm_nat = sqrt(yNorm_nat);

    cout << "===============================================" << endl;
    cout << "Periodicity condition after single shooting:   " << endl;
    cout << "Norm |y(0) - y(T)| in system coordinates = " << yNorm << endl;
    cout << "Norm |y(0) - y(T)| in native coordinates = " << yNorm_nat << endl;
    cout << "===============================================" << endl;

    return GSL_SUCCESS;
}



//====================================================================================
// Subroutines
//====================================================================================
/**
 *   \brief Computes M+1 patch points, indexed between 0 and M,
 *          along the T-periodic orbit starting at y0 = [x0, y0, z0, vx0, vy0, vz0, eye(6)]
 *          Each STM is initialized with the identity matrix.
 **/
void lpdyneq_patch_points(const double *y0, gsl_odeiv2_driver *d, double T, double **ym, double *tm, int M)
{
    //====================================================================================
    //Initial conditions + eye(6)
    //====================================================================================
    double ys[42], t = 0.0;
    for(int i = 0; i < 42; i++) ys[i] = y0[i];

    //First patch point is y0
    for(int i = 0; i < 42; i++) ym[i][0] = y0[i];
    tm[0] = 0.0;

    //====================================================================================
    //Loop on the grid [1...M]
    //====================================================================================
    for(int k = 1; k <= M; k++)
    {
        //Time
        tm[k] = 1.0*k*T/M;
        gsl_odeiv2_driver_reset(d);

        //Integration
        gsl_odeiv2_driver_apply (d, &t, tm[k], ys);

        //Storage in ym - each STM is initialized with the identity matrix
        for(int i= 0; i < 6;  i++)  ym[i][k] = ys[i];
        for(int i= 6; i < 42;  i++) ym[i][k] = y0[i];

        //The STM is restarded: STM(tk) = eye(6);
        for(int j=6; j<42; j++) ys[j] = y0[j];
    }
}


/**
 *  \brief Inverse a 2x2 matrix. Used in lpdyneq_cont.
 **/
void invMat22(double A[2][2], double invA[2][2])
{
    double deter      = +A[0][0]*A[1][1] - A[0][1]*A[1][0];
    invA[0][0] = +A[1][1]/deter;
    invA[0][1] = -A[0][1]/deter;
    invA[1][0] = -A[1][0]/deter;
    invA[1][1] = +A[0][0]/deter;
}

/**
 *  \brief Arclength vector field. For all i in [[1,8]], fi = Ai/sum(Aj^2).
 **/
void arclengthvf(double *f, double *Aj)
{
    //div = sqrt(sum(Aj^2)), j = 0 to 7
    double div = 0;
    for(int i = 0; i < 7; i++) div += Aj[i]*Aj[i];
    div = sqrt(div);

    //Equations of motions
    for(int i = 0; i < 7; i++)
    {
        f[i] = +Aj[i]/div;
    }
}
