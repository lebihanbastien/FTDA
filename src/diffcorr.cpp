/**
 * \file  diffcorr.cpp
 * \brief Contains all the routines that perform differential correction procedures.
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "diffcorr.h"

//------------------------------------------------------------------------------------------------------------
//Differential correction
//------------------------------------------------------------------------------------------------------------

/**
 *  \brief Performs a differential correction procedure on ystart in order to get a periodic orbit of period t1.
 *         The algorithm assumes that the orbit is in the xy plane, is symmetric wrt to the x-axis, and has a period T = t1.
 *
 *  \param ystart the initial conditions that needs to be corrected.
 *  \param t1 the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param N:  the number of points in the plot if the steps of the diffcorr are plotted (see isPlotted).
 *  \param isPlotted:  if true, the steps of the diffcorr procedure are plotted in a temporary gnuplot window.
 *
 *   In brief: the algorithm corrects  [ ystart[0] ystart[4] ] in order to get [y[1] y[3]](t1) = 0
 *             i.e. when the trajectory crosses the line t = t1.
 **/
int differential_correction(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //Init
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 500;

    double dyf[2];
    double dy0[2];
    double t;

    double A[2][2];
    double invA[2][2];
    double deter;

    //Optionnal plotting
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //Correction loop
    do
    {
        //Update the starting point (STM (0) = Id included)
        for(i=0; i< N; i++) y[i] = ystart[i];
        //Integration until t=t1 is reached
        t = 0.0;


        //Optionnal plotting
        if(isPlotted)
        {
            odePlotGen(y, N, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        }

        //Integration
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //Update the matrix A
        A[0][0] = y[5+7];  //Phi21
        A[0][1] = y[5+11]; //Phi25

        A[1][0] = y[5+19]; //Phi41
        A[1][1] = y[5+23]; //Phi45

        //Inverse A
        deter      = +A[0][0]*A[1][1] - A[0][1]*A[1][0];
        invA[0][0] = +A[1][1]/deter;
        invA[0][1] = -A[0][1]/deter;
        invA[1][0] = -A[1][0]/deter;
        invA[1][1] = +A[0][0]/deter;

        //Update the final error
        dyf[0] = -y[1];
        dyf[1] = -y[3];

        //Inversion of the error
        dy0[0] = invA[0][0]*dyf[0] + invA[0][1]*dyf[1];
        dy0[1] = invA[1][0]*dyf[0] + invA[1][1]*dyf[1];


        //Update the state
        ystart[0] += dy0[0];
        ystart[4] += dy0[1];
    }
    while(fabs(dyf[0])> eps_diff && fabs(dyf[1])> eps_diff && (++iter) < itermax);


    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dyf[0]));
        return GSL_FAILURE;
    }

    gnuplot_close(hc);


    return GSL_SUCCESS;
}

/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //Retrieve the model
    QBCP_I* qbpi =  (QBCP_I *) d->sys->params;

    //Init
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double dx;

    //Gsl matrices and vectors
    gsl_matrix *DP   = gsl_matrix_alloc(6,7);
    gsl_matrix *Prod = gsl_matrix_alloc(6,6);
    gsl_vector *P    = gsl_vector_alloc(6);
    gsl_vector *P1   = gsl_vector_alloc(6);
    gsl_vector *Pn   = gsl_vector_alloc(7);

    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (6);

    //Optionnal plotting
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //Correction loop
    do
    {
        //Update the starting point (STM (0) = Id included)
        for(i=0; i< N; i++) y[i] = ystart[i];

        //Integration until t=t1 is reached
        t = 0.0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, N, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //Update the matrix DP = DP(x^k), a 6x7 matrix
        gslc_vectorToMatrix(DP, y, 6, 7, 6);
        //maybe need to subtract the identity here??
        //for(int i = 0; i < 6; i++) gsl_matrix_set(DP, i, i, gsl_matrix_get(DP, i, i)-1.0);

        //Update P(x^k) = x^k(T) - x^k(0) a 6x1 vector
        for(int i = 0; i < 6; i++) gsl_vector_set(P, i, y[i] - ystart[i]);

        //Compute Prod = DP(x^k)*DP(x^k)^T, a 6x6 matrix
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DP , DP, 0.0, Prod);

        //Inverse and product
        int s;
        gsl_linalg_LU_decomp (Prod, p , &s);
        //P1(x^k) = Prod^{-1}*P(x^k), a 6x1 vector
        gsl_linalg_LU_solve(Prod, p, P, P1);

        //Pn(x^k) = DP(x^k)^T*P1(x^k), a 7x1 vector
        gsl_blas_dgemv (CblasTrans, 1.0, DP, P1, 0.0, Pn);

        //Update the state
        for(int i = 0; i < 6; i++) ystart[i] -= gsl_vector_get(Pn, i);
        qbpi->epsilon -= gsl_vector_get(Pn, 6);

        //Norm of the correction
        dx = gsl_blas_dnrm2(Pn);

        cout << "Diffcorr: " << iter << "  dx = " << dx << endl;
    }
    while((fabs(dx)> eps_diff) && (++iter) < itermax);


    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dx));
        return GSL_FAILURE;
    }

    gsl_matrix_free(DP);
    gsl_matrix_free(Prod);
    gsl_vector_free(P);
    gsl_vector_free(P1);
    gsl_vector_free(Pn);
    gnuplot_close(hc);


    return GSL_SUCCESS;
}

/*
 Differential correction - Minimum norm in all 6 dimensions
 */
int differential_correction_MN(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int isPlotted)
{
    //Init
    double y[42];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double dyf[2];
    double t;

    //Gsl matrices
    gsl_matrix *DP   = gsl_matrix_alloc(6,6);
    gsl_matrix *id   = gsl_matrix_alloc(6,6);

    //Gsl vector
    gsl_vector *v0n  = gsl_vector_alloc(6);
    gsl_vector *v1n  = gsl_vector_alloc(6);

    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (6);

    //Identity
    gsl_matrix_set_identity(id);

    //Optionnal plotting
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //Correction loop
    do
    {
        //Update the starting point (STM (0) = Id included)
        for(i=0; i<42; i++) y[i] = ystart[i];

        //Integration until t=t1 is reached
        t = 0.0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, 42, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //Update the matrix A = Df
        gslc_vectorToMatrix(DP, y, 6, 6, 6);
        //Substract the identity: A = Df -id
        gsl_matrix_sub(DP, id);

        //Update v0n
        for(int i = 0; i < 6; i++) gsl_vector_set(v0n, i, y[i]-ystart[i]);

        //Inverse A-Id
        int s;
        gsl_linalg_LU_decomp (DP, p , &s);

        //v1n = DP^{-1}*v0n
        gsl_linalg_LU_solve(DP, p, v0n, v1n);


        //Update the state
        for(int i = 0; i < 6; i++) ystart[i] -= gsl_vector_get(v1n, i);
    }
    while(fabs(gsl_vector_max(v1n)) > eps_diff && (++iter) < itermax);


    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dyf[0]));
        return GSL_FAILURE;
    }

    gnuplot_close(hc);

    gsl_matrix_free(DP);
    gsl_matrix_free(id);
    gsl_permutation_free(p);
    gsl_vector_free(v0n);
    gsl_vector_free(v1n);


    return GSL_SUCCESS;
}

/*
    Differential correction in the case of minimum norm (MN)
    FOR LYAPUNOV ORBITS ONLY
*/
int differential_correction_MN_planar(double ystart[], double yhalf[], double *tcross, double t1, OdeStruct *ode_s)
{

    //Optionnal plotting
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //Current error on -pxf
    double dyf[1];
    double y[42];
    int iter = 0;
    int status;
    do
    {
        //Update the starting point
        for(int i=0; i<42; i++) y[i] = ystart[i];

        //Integration until y=0 is reached
        status = custom_odezero(y, tcross, t1, ode_s);
        odePlotGen(ystart, 42, *tcross, ode_s->d, hc, 5000, "", "lines", "1", "2", 5);
        gsl_odeiv2_driver_reset(ode_s->d);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: odezero failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //Af matrix
        double Af[1][2];
        Af[0][0]= y[5+19];  //Phi41
        Af[0][1]= y[5+23];  //Phi45
        //Bf vector
        double Bf[2];
        Bf[0] = y[5+7];  //Phi21
        Bf[1] = y[5+11]; //Phi25

        //Derivative of y @ y= t12 in ydot
        double ydot[42];
        qbfbp_vf(*tcross, y, ydot, ode_s->d->sys->params);

        //ppf vector
        double ppf = ydot[3]/y[4]; //dot(px)/dot(y)

        //B
        double B[1][2];
        B[0][0] = Af[0][0] - ppf*Bf[0];
        B[0][1] = Af[0][1] - ppf*Bf[1];


        //Inversion of the error (minimum norm solution)
        //-----------------------------------------------
        //Final error
        dyf[0] = -y[3];

        //B*Az_new^T
        double BBT;
        BBT = B[0][0]*B[0][0] + B[0][1]*B[0][1];

        //Af*inv(B*Az_new^T)*dyf = dy0
        double dy0[3];
        dy0[0] = B[0][0]/BBT*dyf[0];
        dy0[1] = B[0][1]/BBT*dyf[0];

        //Updating the state
        ystart[0] = ystart[0] + dy0[0];  //x0
        ystart[4] = ystart[4] + dy0[1];  //py0

        cout << "iter : " << iter << endl;
        cout << "fabs(dyf[0]) = " << fabs(dyf[0]) << endl;

        cout << "dy0[0] = " << dy0[0] << endl;
        cout << "dy0[1] = " << dy0[1] << endl;

    }
    while(fabs(dyf[0])> 1e-11 && (++iter) < 50);

    if(iter>=50)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Premature ending\n");
        return GSL_FAILURE;
    }

    //Update the half-period state
    for(int i=0; i<42; i++) yhalf[i] = y[i];

    return GSL_SUCCESS;
}

/*
    Differential correction in the case of x0 fixed
    FOR LYAPUNOV ORBITS ONLY
*/
int differential_correction_x0_fixed(double ystart[], double yhalf[], double *tcross, double t1, OdeStruct *ode_s)
{
    //Optionnal plotting
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");


    //Current error on -pxf
    double dyf[1];
    double y[42];
    int iter = 0;
    int status;
    do
    {
        //Update the starting point
        for(int i=0; i<42; i++) y[i] = ystart[i];
        //Integration until y=0 is reached
        status = custom_odezero(y, tcross, t1, ode_s);
        if(iter%30 == 0)
        {
            //odePlotGen(ystart, *tcross, ode_s->d, hc, 5000, "lines", "");
            //gsl_odeiv2_driver_reset(ode_s->d);
        }

        if(status != GSL_SUCCESS)
        {
            char ch;
            printf("WARNING: odezero failed to converge in differential_correction. Premature ending.\n");
            odePlotGen(ystart, 42, *tcross, ode_s->d, hc, 5000, "LastDiffCorr", "lines", "1", "2", 5);
            gsl_odeiv2_driver_reset(ode_s->d);
            printf("Press ENTER to close the gnuplot window(s)\n");
            scanf("%c",&ch);
            return GSL_FAILURE;
        }

        //Af matrix
        double Af = y[5+23];  //Phi45
        //Bf vector
        double Bf = y[5+11];  //Phi25
        //Derivative of y @ y= t12 in ydot
        double ydot[42];
        qbfbp_vf(*tcross, y, ydot, ode_s->d->sys->params);
        //ppf vector
        double ppf = ydot[3]/y[4]; //dot(px)/dot(y)
        //B
        double B = Af - ppf*Bf;

        //Inversion of the error
        //-----------------------------------------------
        //Final error
        dyf[0] = -y[3];
        double dy0[1];
        dy0[0] = 1.0/B*dyf[0];

        //Updating the state
        //-----------------------------------------------
        if(dy0[0] > 1e-5) dy0[0] *= 1e-1;
        ystart[4] = ystart[4] + dy0[0];  //py0

        cout << "iter : " << iter << endl;
        cout << "fabs(dyf[0]) = " << fabs(dyf[0]) << endl;
        cout << "dy0[0] = " << dy0[0] << endl;

    }
    while(fabs(dyf[0])> 1e-11 && (++iter) < 100);

    if(iter>= 100)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Premature ending\n");
        return GSL_FAILURE;
    }

    //Update the half-period state
    for(int i=0; i<42; i++) yhalf[i] = y[i];

    return GSL_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------
// ODE PLOT
//------------------------------------------------------------------------------------------------------------

/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1.
 **/
int odePlot(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color)
{
    gsl_odeiv2_driver_reset(d);

    // EM state on a Npoints grid
    double xEM[Npoints];
    double yEM[Npoints];

    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) d->sys->params;

    //Initial conditions
    double ys[N], ye[N];
    for(int i=0; i<N;i++) ys[i] = y[i];
    double ti = 0;

    //First point
    NCtoSYS(0.0, ys, ye, qbp);
    xEM[0] = ye[0];
    yEM[0] = ye[1];

    //Loop and integration
    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i =1; i<= Npoints; i++)
    {
        ti = i * t1 / Npoints;
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Update the SYS state
        NCtoSYS(ti, ys, ye, qbp);
        xEM[i] = ye[0];
        yEM[i] = ye[1];
    }
    if(t1 <0) d->h = -d->h;

    //Plotting
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)"SYS coordinates", "lines", "1", "4", 4);

    return GSL_SUCCESS;

}

/**
 *  \brief Same as odePlot, but with more flexibility on the parameters: can choose the title, the line style, the line type and the line color.
 **/
int odePlotGen(const double y[],
               int N,
               double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc)
{
    gsl_odeiv2_driver_reset(d);

    double xc[Npoints];
    double yc[Npoints];

    //Initial conditions
    double ys[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
    double ti = 0;

    xc[0] = ys[0];
    yc[0] = ys[1];
    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i =1; i<= Npoints; i++)
    {
        ti = i * t1 / Npoints;
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        xc[i] = ys[0];
        yc[i] = ys[1];
    }
    if(t1 <0) d->h = -d->h;

    //------------------------------------------
    //Plot
    //------------------------------------------
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xc, yc, Npoints, title, ls, lt, lw, lc);

    return GSL_SUCCESS;
}

/**
 *  \brief Same as odePlotGen, but the 3D result are also printed in a txt file, with the name (folder+"orbits/"+title+".txt").
 **/
int odePrintGen(const double y[],
               int N,
               double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc,
               string folder)
{
    gsl_odeiv2_driver_reset(d);

    double xc[Npoints];
    double yc[Npoints];
    double zc[Npoints];
    double tc[Npoints];

    //Initial conditions
    double ys[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
    double ti = 0;

    //Loop and integration
    xc[0] = ys[0];
    yc[0] = ys[1];
    zc[0] = ys[2];
    tc[0] = 0.0;
    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i =1; i<= Npoints; i++)
    {
        ti = i * t1 / Npoints;
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        xc[i] = ys[0];
        yc[i] = ys[1];
        zc[i] = ys[2];
        tc[i] = ti;
    }
    if(t1 <0) d->h = -d->h;

    //------------------------------------------
    //Plot
    //------------------------------------------
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xc, yc, Npoints, title, ls, lt, lw, lc);

    //------------------------------------------
    //In txt files
    //------------------------------------------
    string eXYZstr = (folder+"orbits/"+title+".txt");
    gnuplot_fplot_txyz(tc, xc, yc, zc, Npoints, eXYZstr.c_str());   //orbit

    return GSL_SUCCESS;
}


