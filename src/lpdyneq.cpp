#include "qbcp.h"

/**
 * \file lpdyneq.cpp
 * \brief Computation of the dynamical equivalent to the Libration points for various models. Declaration are in qbcp.h for now.
 * \author BLB.
 * \date November 2015
 * \version 1.0
 */


//-----------------------------------------------------------------------------
//Dynamical equivalents to the Libration points
//-----------------------------------------------------------------------------
/**
 *  \brief Main routine for the computation of the dynamical equivalent to the Libration points.
 **/
void lpdyneq(gsl_odeiv2_driver *d, double y0[], gnuplot_ctrl *h1)
{
    //Settings for iostd
    coutlp();

    //Retrieving the parameters
    QBCP_L* qbp      = (QBCP_L *) d->sys->params;
    int isNormalized =  qbp->isNormalized;
    int li           =  qbp->cs.li;
    double tend      =  qbp->us.T;

    //----------------------------------------------------------------------------------------------------------
    // Initial correction
    //----------------------------------------------------------------------------------------------------------
    //Position & Velocity
    double yv[6], yf[6];

    //--------------------
    //Approximated position
    //--------------------
    if(isNormalized)
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

        switch(SEML.fwrk)
        {
        case F_EM:
        {
            //To momenta
            EMvtoEMm(0.0, yv, yf, &SEML);
            //To NC
            EMtoNC(0.0, yf, y0, qbp);
            break;
        }

        case F_SEM:
        {
            //To momenta
            SEMvtoSEMm(0.0, yv, yf, &SEML);
            //To NC
            SEMtoNC(0.0, yf, y0, qbp);
            break;
        }

        }


        if(qbp->model == M_BCP && qbp->li_EM == 1)
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

        switch(SEML.fwrk)
        {
        case F_EM:
        {
            //To momenta
            EMvtoEMm(0.0, yv, y0, &SEML);
            break;
        }

        case F_SEM:
        {
            //To momenta
            SEMvtoSEMm(0.0, yv, y0, &SEML);
            break;
        }

        }
    }

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_alloc(6,6);
    gsl_matrix_set_identity (Id);

    //Storing eye(6) into the initial vector
    gslc_matrixToVector(y0, Id, 6, 6, 6);

    //----------------------------------------------------------------------------------------------------------
    // Differential correction scheme
    //----------------------------------------------------------------------------------------------------------
    if(qbp->model != M_ERTBP)
    {
        cout << "lpdyneq. Starting differential correction..." << endl;
        if(isNormalized) differential_correction(y0, 0.5*tend, 1e-14, d, 42, 0);
        else differential_correction(y0, 0.5*tend, 1e-15, d, 42, 0);
        cout << "lpdyneq. End of differential correction." << endl;
    }

    //----------------------------------------------------------------------------------------------------------
    // Print result
    //----------------------------------------------------------------------------------------------------------
    //Final IC
    cout << "lpdyneq. final IC(NC):" << endl;
    for(int i =0; i<6; i++) cout << y0[i] << endl;
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


void arclengthvf(double *f, double *Aj)
{
    //div = sqrt(sum(Aj^2)), j = 0 to 6
    double div = 0;
    for(int i = 0; i < 7; i++) div += Aj[i]*Aj[i];
    div = sqrt(div);

    //Equations of motions
    for(int i = 0; i < 7; i++)
    {
        f[i] = +Aj[i]/div;
    }
}

/**
 *  \brief Computation of lpdyneq with continuation between two models. WORK IN PROGRESS.
 **/
void lpdyneq_cont(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[], gnuplot_ctrl *h1)
{
    //Settings for iostd
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //Init the plot
    h1 = gnuplot_init();

    //Retrieving the parameters
    QBCP_I* qbpi     =  (QBCP_I *) d->sys->params;
    QBCP_L* qbp1     =  &qbpi->model1;              //First model, full model when qpbi->epsilon = 0.0, at the beginning of the process.
    QBCP_L* qbp2     =  &qbpi->model2;              //First model, full model when qpbi->epsilon = 1.0, at the end of the process.

    int isNormalized =  qbp1->isNormalized;
    int li           =  qbp1->cs.li;

    //Full period of the Sun
    double tend = 2*M_PI/qbp2->us_em.n;

    //----------------------------------------------------------------------------------------------------------
    // Evaluate the alphas
    //----------------------------------------------------------------------------------------------------------
    double alphaf[8];
    evaluateCoef(alphaf, 0.0, qbp1->us.n,  qbp1->nf, qbp1->cs.coeffs, 8);

    //----------------------------------------------------------------------------------------------------------
    // Initial conditions
    //----------------------------------------------------------------------------------------------------------
    //CR3BPs
    CR3BP EM;
    init_CR3BP(&EM, EARTH, MOON);
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
    if(isNormalized)
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
    gsl_matrix *Id = gsl_matrix_alloc(6,6);
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
    gsl_matrix *A = gsl_matrix_alloc(6,6);
    gsl_matrix *ALU = gsl_matrix_alloc(6,6);
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


    //Loop
    int iter    = 0;
    qbpi->epsilon = 0.0;
    do
    {
        //------------------------
        //At the end of the integration, y0[6,...41] contains Df/Dy, y0[42...47] contains Df/Depsilon
        //------------------------
        gsl_odeiv2_driver_reset(d);
        for(int i = 0; i <48; i++) y0s[i] = y0[i];
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

    }while(/*iter <= Nmax &&*/ qbpi->epsilon <= 1.001);

    //Set hard control
    gsl_odeiv2_control_set_driver (hard_control, d);

    cout << "here, after loop." << endl;

    //Last correction much more precise
    if(isNormalized) differential_correction(y0, 0.5*tend, 1e-14, d, 48, 0);
    else differential_correction(y0, 0.5*tend, 1e-14, d, 48, 0);


    //Final IC
    cout << "lpdyneq. IC(NC):" << endl;
    for(int i =0; i<6; i++) cout << y0[i] << endl;
}

