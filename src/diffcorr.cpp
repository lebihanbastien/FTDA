/**
 * \file  diffcorr.cpp
 * \brief Contains all the routines that perform differential correction procedures.
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "diffcorr.h"


//----------------------------------------------------------------------------------------
//          Multiple shooting
//----------------------------------------------------------------------------------------
/**
 *  \brief Solve a definite-positive linear system:
 *         - Decompose a ncs x ncs matrix M that is definite-positive, using GSL routines.
 *         - Then inverse the system Fv = M*K3.
 **/
int ftc_inv_dfls(gsl_matrix* M, gsl_vector*Fv, gsl_vector* K3, int ncs)
{
    //Name of the routine
    string fname = "ftc_inv_dfls";

    //====================================================================================
    //We start by turning off the handler in the next chunk of code, so that we can test
    //if M is truly definite-positive.
    //We will restore it at the end of the decomposition + inversion process.
    //====================================================================================
    gsl_error_handler_t* error_handler = gsl_set_error_handler_off();

    //====================================================================================
    //Since M is in theory definite-positive,
    //a Cholesky decomposition can be used, instead of a LU decomposition.
    //====================================================================================
    int status = gsl_linalg_cholesky_decomp (M);

    //====================================================================================
    //If the decomposition went bad, status is different from GSL_SUCCESS. For example,
    //if M is not strictly definite-positive, the constant GSL_EDOM is returned by
    //gsl_linalg_cholesky_decomp. In any case, we try a LU decomposition.
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Cholesky decomposition failed. A LU decomposition is tested. ref_errno = " << gsl_strerror(status) << endl;
        int s;
        gsl_permutation* p = gsl_permutation_alloc (ncs);
        status = gsl_linalg_LU_decomp (M, p, &s);
        if(!status) status = gsl_linalg_LU_solve(M, p, Fv, K3);
        gsl_permutation_free (p);

        //--------------------------------------------------------------------------------
        //If the decomposition or solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". The LU decomposition failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return GSL_FAILURE;
        }

    }
    else  //if not, the decomposition went fine, and we can go on with the inversion.
    {
        status = gsl_linalg_cholesky_solve(M, Fv, K3);
        //--------------------------------------------------------------------------------
        //If the solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". Cholesky solving failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return GSL_FAILURE;
        }
    }

    //We check that the result is not undefined
    if(gsl_isnan(gsl_vector_max(K3)))
    {
        cerr << fname << ". The decomposition + solving failed, result contains NaN" << endl;
        //Finally, we return an error to the user
        gsl_set_error_handler (error_handler);
        return GSL_FAILURE;
    }


    //====================================================================================
    //We restore the error handler at the end
    //====================================================================================
    gsl_set_error_handler (error_handler);
    return GSL_SUCCESS;
}

/**
 *  \brief Computes the correction vector associated to the minimum norm solution.
 *         Given:
 *              - an ncs x 1   error vector Fv
 *              - an nfv x ncs Jacobian DF,
 *         This routine computes the correction vector associated to
 *         the minimum norm solution:
 *
 *              DQv = DF^T x (DF x DF^T)^{-1} Fv.
 **/
int ftc_corrvec_mn(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int nfv, int ncs)
{
    //Name of the routine
    string fname = "ftc_corrvec_mn";

    //====================================================================================
    // Initialize the local GSL objects
    //====================================================================================
    gsl_matrix *M   = gsl_matrix_calloc(ncs, ncs);
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //====================================================================================
    // Update M = DF x DF^T
    //====================================================================================
    int status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);
    if(status)
    {
        cerr << fname << ". The computation of M failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Update correction vector DQv = DF^T x (DF x DF^T)^{-1} Fv
    //====================================================================================
    //Inverse Fv = M*K3 via Cholesky decomposition of M
    status = ftc_inv_dfls(M, Fv, K3, ncs);
    //We check that the inversing went well
    if(status)
    {
        cerr << fname << ". The decomposition + solving failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return GSL_FAILURE;
    }
    //Then, DQv = DF^T*K3
    status = gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);
    if(status)
    {
        cerr << fname << ". The computation of DQv failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Kill the local GSL objects and return GSL_SUCCESS
    //====================================================================================
    gsl_matrix_free(M);
    gsl_vector_free(K3);
    return GSL_SUCCESS;
}

/**
 *  \brief Computes the correction vector for a square system.
 *         Given:
 *              - an ncs x 1   error vector Fv
 *              - an ncs x ncs Jacobian DF,
 *         This routine computes the correction vector given by:
 *
 *              DQv = DF^{-1} x Fv.
 **/
int ftc_corrvec_square(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int ncs)
{
    //Name of the routine
    string fname = "ftc_corrvec_square";

    //====================================================================================
    // Initialize the local GSL objects
    //====================================================================================
    gsl_matrix *M   = gsl_matrix_calloc(ncs, ncs);
    gsl_permutation* p = gsl_permutation_alloc (ncs);
    int s;


    //M = DF
    gsl_matrix_memcpy(M, DF);

    //====================================================================================
    // Update correction vector DQv = DF^{-1} x Fv
    //====================================================================================
    //Inverse Fv = M*DQv via LU decomposition of M
    int status = gsl_linalg_LU_decomp (M, p, &s);
    if(!status) status = gsl_linalg_LU_solve(M, p, Fv, DQv);

    //--------------------------------------------------------------------------------
    //If the decomposition or solving went bad, status is different from GSL_SUCCESS.
    //--------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". The LU decomposition failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_permutation_free (p);
        gsl_matrix_free(M);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Kill the local GSL objects and return FTC_SUCCESS
    //====================================================================================
    gsl_permutation_free (p);
    gsl_matrix_free(M);
    return GSL_SUCCESS;
}


/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int lpdyneq_multiple_shooting(double** ymd, double* tmd, double** ymdn, double* tmdn,
                              gsl_odeiv2_driver* d, int N, int mgs, double prec,
                              int isPlotted, gnuplot_ctrl* h1, int strConv)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Name of the routine
    string fname = "lpdyneq_multiple_shooting";

    //Cumulated norm of the error
    double normC = 0.0, normC_old = 1e5;
    int status = GSL_SUCCESS;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*(mgs+1);  //free variables
    int ncs = 6*(mgs);  //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Jacobian at patch points
    gsl_matrix **Ji = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF  = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M   = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    gsl_odeiv2_driver_reset(d);
    while(iter < itermax)
    {
        //--------------------------------------------------------------------------------
        //Trajectory plot
        //--------------------------------------------------------------------------------
        //gnuplot_plot_xyz(h1, ymdn[0], ymdn[1], ymdn[2], mgs+1, (char*)"", "points", "1", "4", 4);

        //--------------------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //--------------------------------------------------------------------------------
        //#pragma omp parallel for if(isPar)
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Initialization for // computing
            //----------------------------------------------------------------------------
            double ye[N], tv = 0;

            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ymdn[i][k];
            tv = tmdn[k];

            //Storing eye(6) into the initial vector
            gslc_matrixToVector(ye, Id, 6, 6, 6);

            //Integration in [tmdn[k], tmdn[k+1]]
            gsl_odeiv2_driver_apply (d, &tv, tmdn[k+1], ye);

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                if(k != mgs-1) gsl_vector_set(Fv, i + 6*k, ye[i] - ymdn[i][k+1]);
                else gsl_vector_set(Fv, i + 6*k, ye[i] - ymdn[i][0]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                if(k != mgs-1)
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 6*(k+1), -gsl_matrix_get(Id, i, j));
                    }
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j      ,     -gsl_matrix_get(Id, i, j));
                    }


                }
            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        normC  = gsl_blas_dnrm2(Fv);


        //--------------------------------------------------------------------------------
        // Check that all points are under a given threshold
        //--------------------------------------------------------------------------------
        cout << fname << ". n° " << iter+1 << "/" << itermax << ". nerror = " << normC << endl;
        if(normC < prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }


        //--------------------------------------------------------------------------------
        // If a strong convergence is desired, check that we are always decreasing the error
        //--------------------------------------------------------------------------------
        if(strConv && normC_old < normC)
        {
            cout << fname << ". error[n+1] < error[n]. " << endl;
            return GSL_FAILURE;
        }
        normC_old = normC;


        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        //status = ftc_corrvec_square(DQv, Fv, DF, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return GSL_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        for(int k = 0; k <= mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 6*k+i);
        }

        //--------------------------------------------------------------------------------
        // Update number of iterations
        //--------------------------------------------------------------------------------
        iter++;
    }

    //====================================================================================
    // Check that the desired condition was met
    //====================================================================================
    if(normC > prec)
    {
        cerr << fname << ". The desired precision was not reached."  << endl;
        return GSL_FAILURE;
    }


    //====================================================================================
    // Free
    //====================================================================================
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}


//----------------------------------------------------------------------------------------
//Differential correction
//----------------------------------------------------------------------------------------
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
    //====================================================================================
    //Init
    //====================================================================================
    double y[N], yc[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 1000;

    double dyf[2];
    double dy0[2];
    dy0[0] = 0.0;
    dy0[1] = 0.0;
    double t;

    double A[2][2];
    double invA[2][2];
    double deter;
    double eps_c = 10;
    double eps_m = 10;

    //====================================================================================
    //Optionnal plotting
    //====================================================================================
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //====================================================================================
    //Update the starting point (STM (0) = Id included)
    //====================================================================================
    for(i=0; i< N; i++) yc[i] = ystart[i];

    //====================================================================================
    //Correction loop
    //====================================================================================
    do
    {
        //================================================================================
        //Update the starting point (STM (0) = Id included)
        //================================================================================
        for(i=0; i< N; i++) y[i] = yc[i];
        //Integration until t=t1 is reached
        t = 0.0;

        //Optionnal plotting
        if(isPlotted)
        {
            odePlotGen(y, N, t1, d, hc, 500, "DiffCorr", "lines", "1", "2", 8);
        }

        //================================================================================
        //Integration
        //================================================================================
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //================================================================================
        //Update the matrix A
        //================================================================================
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

        //================================================================================
        //Update the final error
        //================================================================================
        dyf[0] = -y[1];
        dyf[1] = -y[3];

        //================================================================================
        //Check that yc leads to a better precision, and update the output if so.
        //================================================================================
        //Current precision
        eps_c = sqrt(dyf[0]*dyf[0] + dyf[1]*dyf[1]);

        //Compare to the best precision yet
        if(eps_c < eps_m)
        {
            eps_m = eps_c;
            for(i=0; i< N; i++) ystart[i] = yc[i];
        }

        //================================================================================
        //Break the loop if precision is good
        //================================================================================
        if(eps_m < eps_diff) break;

        //================================================================================
        //Inversion of the error
        //================================================================================
        dy0[0] = invA[0][0]*dyf[0] + invA[0][1]*dyf[1];
        dy0[1] = invA[1][0]*dyf[0] + invA[1][1]*dyf[1];

        //================================================================================
        //Update the state
        //================================================================================
        yc[0] += 0.1*dy0[0];
        yc[4] += 0.1*dy0[1];
    }
    while((++iter) < itermax);

    //====================================================================================
    // Warning message if no convergence
    //====================================================================================
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", eps_m);
        return GSL_FAILURE;
    }

    //====================================================================================
    // Display info
    //====================================================================================
    cout << "differential_correction. Required eps = " << eps_diff << endl;
    cout << "differential_correction. Obtained eps = " << eps_m    << endl;

    gnuplot_close(hc);
    return GSL_SUCCESS;
}

/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps_mns(double *ystart, double *nullvector, double *fvv, double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted, int isFirst)
{
    //-------------------------------------------------------------------------------
    //Retrieve the model
    //-------------------------------------------------------------------------------
    QBCP_I* qbpi =  (QBCP_I *) d->sys->params;

    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double dx;


    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *DP   = gsl_matrix_calloc(2,3);
    gsl_matrix *Prod = gsl_matrix_calloc(2,2);
    gsl_vector *P    = gsl_vector_calloc(2);
    gsl_vector *P1   = gsl_vector_calloc(2);
    gsl_vector *Pn   = gsl_vector_calloc(3);
    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (2);
    int s;

    //-------------------------------------------------------------------------------
    //Optionnal plotting
    //-------------------------------------------------------------------------------
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //-------------------------------------------------------------------------------
    //Correction loop
    //-------------------------------------------------------------------------------
    do
    {
        //-------------------------------------------------------------------------------
        //Update the starting point (STM (0) = Id included)
        //-------------------------------------------------------------------------------
        for(i=0; i< N; i++) y[i] = ystart[i];

        //-------------------------------------------------------------------------------
        //Integration until t=t1 is reached
        //-------------------------------------------------------------------------------
        t = 0.0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, N, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        //Integration
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //-------------------------------------------------------------------------------
        // Update the Jacobian
        //-------------------------------------------------------------------------------
        //Update the matrix DP = DP(x^k), a 2x3 matrix
        gsl_matrix_set(DP, 0, 0, y[5+7]);  //Phi21
        gsl_matrix_set(DP, 0, 1, y[5+11]); //Phi25
        gsl_matrix_set(DP, 1, 0, y[5+19]); //Phi41
        gsl_matrix_set(DP, 1, 1, y[5+23]); //Phi45
        gsl_matrix_set(DP, 0, 2, y[43]);   //dy/deps
        gsl_matrix_set(DP, 1, 2, y[45]);   //dpx/deps

        //-------------------------------------------------------------------------------
        // Update the error vector = [y px]
        //-------------------------------------------------------------------------------
        gsl_vector_set(P, 0, y[1]);
        gsl_vector_set(P, 1, y[3]);

        //-------------------------------------------------------------------------------
        //Minimum norm solution
        //-------------------------------------------------------------------------------
        //Update the matrix
        //Compute Prod = DP(x^k)*DP(x^k)^T, a 2x2 matrix
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DP , DP, 0.0, Prod);
        //Inverse and product
        gsl_linalg_LU_decomp (Prod, p , &s);
        //P1(x^k) = Prod^{-1}*P(x^k), a 2x1 vector
        gsl_linalg_LU_solve(Prod, p, P, P1);
        //Pn(x^k) = DP(x^k)^T*P1(x^k), a 3x1 vector
        gsl_blas_dgemv (CblasTrans, 1.0, DP, P1, 0.0, Pn);

        //-------------------------------------------------------------------------------
        //Update the state
        //-------------------------------------------------------------------------------
        ystart[0]     -= gsl_vector_get(Pn, 0);
        ystart[4]     -= gsl_vector_get(Pn, 1);
        qbpi->epsilon -= gsl_vector_get(Pn, 2);

        //-------------------------------------------------------------------------------
        //Norm of the error
        //-------------------------------------------------------------------------------
        dx = gsl_blas_dnrm2(P);
    }
    while((fabs(dx)> eps_diff) && (++iter) < itermax);

    cout << "dx = " << dx << endl;
    //-------------------------------------------------------------------------------
    //End the computation if there was no convergence
    //-------------------------------------------------------------------------------
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dx));
        return GSL_FAILURE;
    }

    //-------------------------------------------------------------------------------
    //Compute the "good" vector of free variables
    //-------------------------------------------------------------------------------
    fvv[0] = ystart[0];
    fvv[1] = ystart[4];
    fvv[2] = qbpi->epsilon;

    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(2);
    gsl_matrix *Q     = gsl_matrix_calloc(3,3);
    gsl_matrix *R     = gsl_matrix_calloc(3,2);
    gsl_matrix *DPT   = gsl_matrix_calloc(3,2);

    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DPT, DP);
    //QR decomposition
    gsl_linalg_QR_decomp (DPT, work);
    gsl_linalg_QR_unpack (DPT, work, Q, R);

    //----------------------------------
    // The distinction is made between
    // the first use of the routine and the other uses
    //----------------------------------
    double nv[3];
    double dotNV = 0;
    int sign = 1;

    if(isFirst)
    {
        //If it is the first use, we go in the direction so that depsilon > 0
        sign = gsl_matrix_get(Q, 2, 2) > 0? -1:+1;
        for(int i = 0; i < 3; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, 2);
    }
    else
    {
        //If not, we follow the last direction of the null vector.
        for(int i = 0; i < 3; i++)
        {
            nv[i] = gsl_matrix_get(Q, i, 2);
            dotNV += nv[i]*nullvector[i];
        }
        sign = dotNV > 0? 1:-1;
        for(int i = 0; i < 3; i++) nullvector[i] = sign*nv[i];
    }


    //-------------------------------------------------------------------------------
    //Free memory
    //-------------------------------------------------------------------------------
    gsl_matrix_free(DP);
    gsl_matrix_free(DPT);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_matrix_free(Prod);
    gsl_vector_free(work);
    gsl_vector_free(P);
    gsl_vector_free(P1);
    gsl_vector_free(Pn);
    gnuplot_close(hc);


    return iter;
}

/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //-------------------------------------------------------------------------------
    //Retrieve the model
    //-------------------------------------------------------------------------------
    QBCP_I* qbpi =  (QBCP_I *) d->sys->params;

    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double dx;

    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *DP   = gsl_matrix_calloc(6,7);
    gsl_matrix *Prod = gsl_matrix_calloc(6,6);
    gsl_vector *P    = gsl_vector_calloc(6);
    gsl_vector *P1   = gsl_vector_calloc(6);
    gsl_vector *Pn   = gsl_vector_calloc(7);
    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (6);

    //-------------------------------------------------------------------------------
    //Optionnal plotting
    //-------------------------------------------------------------------------------
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //-------------------------------------------------------------------------------
    //Correction loop
    //-------------------------------------------------------------------------------
    do
    {
        //-------------------------------------------------------------------------------
        //Update the starting point (STM (0) = Id included)
        //-------------------------------------------------------------------------------
        for(i=0; i< N; i++) y[i] = ystart[i];

        //-------------------------------------------------------------------------------
        //Integration until t=t1 is reached
        //-------------------------------------------------------------------------------
        t = 0.0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, N, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        //Integration
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //-------------------------------------------------------------------------------
        // Update the Jacobian
        //-------------------------------------------------------------------------------
        //Update the matrix DP = DP(x^k), a 6x7 matrix
        gslc_vectorToMatrix(DP, y, 6, 7, 6);
        //maybe need to subtract the identity here??
        //for(int i = 0; i < 6; i++) gsl_matrix_set(DP, i, i, gsl_matrix_get(DP, i, i)-1.0);

        //-------------------------------------------------------------------------------
        // Update the error vector
        //-------------------------------------------------------------------------------
        //Update P(x^k) = x^k(T) - x^k(0) a 6x1 vector
        for(int i = 0; i < 6; i++) gsl_vector_set(P, i, y[i] - ystart[i]);


        //-------------------------------------------------------------------------------
        //Minimum norm solution
        //-------------------------------------------------------------------------------
        //Compute Prod = DP(x^k)*DP(x^k)^T, a 6x6 matrix
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DP , DP, 0.0, Prod);
        //Inverse and product
        int s;
        gsl_linalg_LU_decomp (Prod, p , &s);
        //P1(x^k) = Prod^{-1}*P(x^k), a 6x1 vector
        gsl_linalg_LU_solve(Prod, p, P, P1);
        //Pn(x^k) = DP(x^k)^T*P1(x^k), a 7x1 vector
        gsl_blas_dgemv (CblasTrans, 1.0, DP, P1, 0.0, Pn);

        //-------------------------------------------------------------------------------
        //Update the state
        //-------------------------------------------------------------------------------
        for(int i = 0; i < 6; i++) ystart[i] -= gsl_vector_get(Pn, i);
        qbpi->epsilon -= gsl_vector_get(Pn, 6);

        //-------------------------------------------------------------------------------
        //Norm of the correction
        //-------------------------------------------------------------------------------
        dx = gsl_blas_dnrm2(Pn);
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

/**
 *  \brief Performs a differential correction procedure on ystart in order to get a periodic orbit of period t1.
 *
 *  \param ystart the initial conditions that needs to be corrected.
 *  \param t1 the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param N:  the number of points in the plot if the steps of the diffcorr are plotted (see isPlotted).
 *  \param isPlotted:  if true, the steps of the diffcorr procedure are plotted in a temporary gnuplot window.
 *
 **/
int differential_correction_T(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //=====================================
    //Init
    //=====================================
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 100;
    double t;

    //For GSL inversions
    int s;
    gsl_matrix* A = gsl_matrix_calloc(6,6);
    gsl_vector *dyf = gsl_vector_calloc(6);
    gsl_vector *dy0 = gsl_vector_calloc(6);
    gsl_permutation * p = gsl_permutation_alloc (6);

    //Identity
    gsl_matrix *id   = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity(id);

    //=====================================
    //Optionnal plotting
    //=====================================
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //=====================================
    //Correction loop
    //=====================================
    double prec = 0.0;
    do
    {
        //=====================================
        //Update the starting point (STM (0) = Id included)
        //=====================================
        for(i=0; i< N; i++) y[i] = ystart[i];
        t = 0.0;

        //Optionnal plotting
        if(isPlotted)
        {
            odePlotGen(y, N, t1, d, hc, 5000, "DiffCorr", "lines", "1", "2", 8);
        }

        //=====================================
        //Integration until t = t1
        //=====================================
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , t1 , y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //=====================================
        //Update the matrix A
        //=====================================
        gslc_vectorToMatrix(A, y, 6, 6 ,6);
        //Substract the identity: A = Df -id
        gsl_matrix_sub(A, id);

        //=====================================
        //Update the final error
        //=====================================
        for(int i = 0; i < 6; i++) gsl_vector_set(dyf, i, y[i]-ystart[i]);

        //=====================================
        //Break the loop if precision is good
        //=====================================
        prec = gsl_blas_dnrm2(dyf);
        if(prec < eps_diff) break;
        cout << "differential_correction. Obtained eps = " << prec << endl;

        //=====================================
        //Inversion of the error
        //=====================================
        gsl_linalg_LU_decomp (A, p, &s);
        gsl_linalg_LU_solve(A, p, dyf, dy0);

        //=====================================
        //Update the state
        //=====================================
        for(int i = 0; i < 6; i++) ystart[i] -= gsl_vector_get(dy0, i);

        //char ch;
        //printf("Press ENTER\n");
        //scanf("%c",&ch);
    }
    while((++iter) < itermax);

    //=====================================
    // Warning message if no convergence
    //=====================================
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", prec);
        return GSL_FAILURE;
    }

    //=====================================
    // Display info
    //=====================================
    cout << "differential_correction. Required eps = " << eps_diff << endl;
    cout << "differential_correction. Obtained eps = " << prec     << endl;

    gnuplot_close(hc);

    //=====================================
    // Free
    //=====================================
    gsl_matrix_free(A);
    gsl_vector_free(dyf);
    gsl_vector_free(dy0);
    gsl_permutation_free(p);


    return GSL_SUCCESS;
}

//----------------------------------------------------------------------------------------
// ODE PLOT - NEW VERSION
//----------------------------------------------------------------------------------------
/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in either NC or SYS coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1. Print in txt files is included vis \c isStored integer.
 **/
int odePlot2(const double y[], int N, double t1, gsl_odeiv2_driver *d,
             gnuplot_ctrl  *h1, int Npoints, int color, int isNormalized, int isStored,
             string legend, string filename)
{
    //------------------------------------------
    // Init
    //------------------------------------------
    gsl_odeiv2_driver_reset(d);

    //EM state on a Npoints grid
    double xEM[Npoints];
    double yEM[Npoints];

    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) d->sys->params;

    //Initial conditions
    double ys[N], ye[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
    double ti = 0;

    //First point
    if(isNormalized)
    {
        NCtoSYS(0.0, ys, ye, qbp);
        xEM[0] = ye[0];
        yEM[0] = ye[1];
    }
    else
    {
        xEM[0] = ys[0];
        yEM[0] = ys[1];
    }

    //------------------------------------------
    //Loop and integration
    //------------------------------------------
    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i =1; i<= Npoints; i++)
    {
        ti = i * t1 / Npoints;
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        //Update the SYS state
        if(isNormalized)
        {
            NCtoSYS(ti, ys, ye, qbp);
            xEM[i] = ye[0];
            yEM[i] = ye[1];
        }
        else
        {
            xEM[i] = ys[0];
            yEM[i] = ys[1];
        }
    }
    if(t1 <0) d->h = -d->h;

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)legend.c_str(), "lines", "1", "4", color);


    //------------------------------------------
    //In txt files
    //------------------------------------------
    if(isStored)
    gnuplot_fplot_xy(xEM, yEM, Npoints, filename.c_str());   //orbit

    return GSL_SUCCESS;
}


/**
 *  \brief Integrate the states (ymdn**, tmdn*), in either NC or SYS coordinates,
 *         and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1.
 **/
int odePlotvec(double **ymdn, double *tmdn, int N, int mgs, gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1, int Npoints, int color, int isNormalized, string legend)
{
    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(d);

    //EM state on a Npoints grid
    double xEM[Npoints];
    double yEM[Npoints];

    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) d->sys->params;

    //Initial conditions
    double ys[N], ye[N];

    //------------------------------------------------------------------------------------
    //Loop and integration
    //------------------------------------------------------------------------------------
    double ts, ti;
    for(int k = 0; k < mgs; k++)
    {
        for(int i= 0; i < N; i++) ys[i] = ymdn[i][k];
        ts  = tmdn[k];

        //Update the SYS state
        if(isNormalized)
        {
            NCtoSYS(ts, ys, ye, qbp);
            xEM[0] = ye[0];
            yEM[0] = ye[1];
        }
        else
        {
            xEM[0] = ys[0];
            yEM[0] = ys[1];
        }

        if(tmdn[mgs] < 0) d->h = -d->h;
        for(int i = 1; i <= Npoints; i++)
        {
            ti = tmdn[k] + 1.0*i/Npoints*(tmdn[k+1] - tmdn[k]);

            //Integration
            gsl_odeiv2_driver_apply (d, &ts, ti, ys);

            //Update the SYS state
            if(isNormalized)
            {
                NCtoSYS(ts, ys, ye, qbp);
                xEM[i] = ye[0];
                yEM[i] = ye[1];
            }
            else
            {
                xEM[i] = ys[0];
                yEM[i] = ys[1];
            }
        }
        if(tmdn[mgs] < 0) d->h = -d->h;

        //--------------------------------------------------------------------------------
        //Plotting
        //--------------------------------------------------------------------------------
        // Of each leg
        if(k == 0) gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)legend.c_str(), "lines", "2", "2", color);
        else gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)"", "lines", "2", "2", color);

        // Of each point
        gnuplot_plot_xy(h1, &xEM[0], &yEM[0], 1, (char*)"", "points", "2", "2", color);
    }

    //------------------------------------------------------------------------------------
    // Aesthetics
    //------------------------------------------------------------------------------------
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");

    return GSL_SUCCESS;
}



//----------------------------------------------------------------------------------------
// ODE PLOT: Kept for consistency
//----------------------------------------------------------------------------------------
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
    for(int i=0; i<N; i++) ys[i] = y[i];
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
    gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)"", "lines", "1", "4", color);

    return GSL_SUCCESS;

}

/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1.
 **/
int odePlot3D(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color)
{
    gsl_odeiv2_driver_reset(d);

    // EM state on a Npoints grid
    double xEM[Npoints];
    double yEM[Npoints];
    double zEM[Npoints];

    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) d->sys->params;

    //Initial conditions
    double ys[N], ye[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
    double ti = 0;

    //First point
    NCtoSYS(0.0, ys, ye, qbp);
    xEM[0] = ye[0];
    yEM[0] = ye[1];
    zEM[0] = ye[2];

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
        zEM[i] = ye[2];
    }
    if(t1 <0) d->h = -d->h;

    //Plotting
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, xEM, yEM, zEM, Npoints, (char*)"", "lines", "1", "4", color);

    return GSL_SUCCESS;

}

/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1. Print in txt files is included
 **/
int odePlotprint(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color, string filename)
{
    gsl_odeiv2_driver_reset(d);

    // EM state on a Npoints grid
    double xEM[Npoints];
    double yEM[Npoints];

    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) d->sys->params;

    //Initial conditions
    double ys[N], ye[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
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

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xEM, yEM, Npoints, (char*)"", "lines", "1", "4", color);

    //------------------------------------------
    //In txt files
    //------------------------------------------
    gnuplot_fplot_xy(xEM, yEM, Npoints, filename.c_str());   //orbit

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

//----------------------------------------------------------------------------------------
//Differential correction: not used
//----------------------------------------------------------------------------------------
/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps_pac(double y0[], double nullvector[], double fvv[], double ds, double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //-------------------------------------------------------------------------------
    //Retrieve the model
    //-------------------------------------------------------------------------------
    QBCP_I* qbpi =  (QBCP_I *) d->sys->params;

    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double dx;
    double X0[3];

    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *DG   = gsl_matrix_calloc(3,3);
    gsl_vector *P    = gsl_vector_calloc(3);
    gsl_vector *Pn   = gsl_vector_calloc(3);

    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (3);

    //-------------------------------------------------------------------------------
    //Optionnal plotting
    //-------------------------------------------------------------------------------
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //-------------------------------------------------------------------------------
    //Correction loop
    //-------------------------------------------------------------------------------
    do
    {
        //-------------------------------------------------------------------------------
        //Update the starting point (STM (0) = Id included)
        //-------------------------------------------------------------------------------
        for(i=0; i< N; i++) y[i] = y0[i];

        //Integration until t = t1 is reached
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

        //-------------------------------------------------------------------------------
        //Build the free variables vector = [x0 py0 eps]^T
        //-------------------------------------------------------------------------------
        X0[0] = y0[0];
        X0[1] = y0[4];
        X0[2] = qbpi->epsilon;

        //-------------------------------------------------------------------------------
        // Update the Jacobian
        //-------------------------------------------------------------------------------
        //Update the matrix DP = DP(x^k), a 2x3 matrix
        gsl_matrix_set(DG, 0, 0, y[5+7]);  //Phi21
        gsl_matrix_set(DG, 0, 1, y[5+11]); //Phi25
        gsl_matrix_set(DG, 1, 0, y[5+19]); //Phi41
        gsl_matrix_set(DG, 1, 1, y[5+23]); //Phi45
        gsl_matrix_set(DG, 0, 2, y[43]);   //dy/deps
        gsl_matrix_set(DG, 1, 2, y[45]);   //dpx/deps
        //The last row is equal to the null vector of the previous step
        for(int i = 0; i < 3; i++) gsl_matrix_set(DG, 2, i, nullvector[i]);

        //-------------------------------------------------------------------------------
        //Build the error vector = [y px pseudo-arclength constraint]^T
        //-------------------------------------------------------------------------------
        //Update the first 2 rows of the error: P(x^k) = x^k(T)
        gsl_vector_set(P, 0, y[1]);
        gsl_vector_set(P, 1, y[3]);
        //The last row is equal to: (X0 - fvv)^T * nullvector - ds
        double res = 0;
        for(int i = 0; i < 3; i++)
        {
            X0[i] -= fvv[i];
            res += X0[i]*nullvector[i];
        }
        gsl_vector_set(P, 2, res - ds);

        //-------------------------------------------------------------------------------
        //Get the first-order correction
        //-------------------------------------------------------------------------------
        int s;
        gsl_linalg_LU_decomp (DG, p , &s);
        //P1 = DG^{-1}*P, 3x1 vector
        gsl_linalg_LU_solve(DG, p, P, Pn);

        //-------------------------------------------------------------------------------
        //Update the state
        //-------------------------------------------------------------------------------
        y0[0]         -= gsl_vector_get(Pn, 0);
        y0[4]         -= gsl_vector_get(Pn, 1);
        qbpi->epsilon -= gsl_vector_get(Pn, 2);

        //-------------------------------------------------------------------------------
        //Norm of the correction
        //-------------------------------------------------------------------------------
        dx = gsl_blas_dnrm2(Pn);

    }
    while((fabs(dx)> eps_diff) && (++iter) < itermax);


    //-------------------------------------------------------------------------------
    //End the computation if there was no convergence
    //-------------------------------------------------------------------------------
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dx));
        return GSL_FAILURE;
    }

    //-------------------------------------------------------------------------------
    //Compute the "good" vector of free variables
    //-------------------------------------------------------------------------------
    fvv[0] = y0[0];
    fvv[1] = y0[4];
    fvv[2] = qbpi->epsilon;

    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(2);
    gsl_matrix *Q     = gsl_matrix_calloc(3,3);
    gsl_matrix *R     = gsl_matrix_calloc(3,2);
    gsl_matrix *DPT   = gsl_matrix_calloc(3,2);

    //Get the matrix DP as a submatrix of DG
    gsl_matrix_view DP = gsl_matrix_submatrix(DG, 0, 0, 2, 3);

    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DPT, &DP.matrix);
    //QR decomposition
    gsl_linalg_QR_decomp (DPT, work);
    gsl_linalg_QR_unpack (DPT, work, Q, R);

    //Null vector is the last column of Q
    double nv[3];
    double dotNV = 0;
    for(int i = 0; i < 3; i++)
    {
        nv[i] = gsl_matrix_get(Q, i, 2);
        dotNV += nv[i]*nullvector[i];
    }
    int sign = gsl_matrix_get(Q, 2, 2) > 0? 1:-1;
    //int sign = dotNV > 0? 1:-1;
    for(int i = 0; i < 3; i++) nullvector[i] = sign*nv[i];

    //-------------------------------------------------------------------------------
    //Free memory
    //-------------------------------------------------------------------------------
    gsl_matrix_free(DG);
    gsl_matrix_free(DPT);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_vector_free(Pn);
    gsl_vector_free(P);
    gsl_vector_free(work);
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
    double t;

    //Gsl matrices
    gsl_matrix *DP   = gsl_matrix_calloc(6,6);
    gsl_matrix *id   = gsl_matrix_calloc(6,6);

    //Gsl vector
    gsl_vector *v0n  = gsl_vector_calloc(6);
    gsl_vector *v1n  = gsl_vector_calloc(6);

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
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(gsl_vector_max(v1n)));
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

