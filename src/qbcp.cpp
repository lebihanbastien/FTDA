#include "qbcp.h"

extern "C"{
    #include "nrutil.h"
}


/**
 * \file qbcp.cpp
 * \brief Implementation of various vector fields in the Sun-Earth-Moon QBCP,
 *        as well as additional routines for inner changes of coordinates
 *       (Earth-Moon <-> Inertial <-> Sun-(Earth+Moon).
 * \author BLB.
 */


//----------------------------------------------------------------------------------------
//
// Derivatives for ODE integration
//
//----------------------------------------------------------------------------------------
//-------------------------------------------------------
// Vector fields
//-------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates (no variationnal equations).
 **/
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Memory allocation
    //------------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *Q    = gsl_matrix_calloc(6,6);

    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp  = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //------------------------------------------------------------------------------------
    vfn_stm(y, Q, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //------------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);

    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbfbp_vf(double t, const double y[], double f[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Memory allocation
    //------------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *B    = gsl_matrix_calloc(6,6);

    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, 8);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double Ps[3];
    eval_array_coef(Ps, t, n, qbp->n_order_fourier, qbp->cs.Ps, 3);
    double Pe[3];
    eval_array_coef(Pe, t, n, qbp->n_order_fourier, qbp->cs.Pe, 3);
    double Pm[3];
    eval_array_coef(Pm, t, n, qbp->n_order_fourier, qbp->cs.Pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    vf_state(y, f, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //------------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //------------------------------------------------------------------------------------
    vf_stm(y, B, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //------------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //------------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = B * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, STM, 0.0, STMd);

    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(B);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varlin_trans(double t, const double y[], double f[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Memory allocation
    //------------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *Q    = gsl_matrix_calloc(6,6);

    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    // STM derivatives, linearized case
    //------------------------------------------------------------------------------------
    vfn_stm_lin_trans(y, Q, alpha, ps, pe, pm, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //------------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);

    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}

/**
 *  \brief Variationnal matrix of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_Dfn_varnonlin(double t, const double y[], double **Df, void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp  = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //------------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //------------------------------------------------------------------------------------
    vfn_stm(y, Df, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return GSL_SUCCESS;
}

//-------------------------------------------------------
// Inner routines for the computation of the vector fields
//-------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
void set_vareq_matrix(gsl_matrix *Q, double b[], double alpha[])
{
    //------------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //------------------------------------------------------------------------------------
    gsl_matrix_set(Q, 0, 0,  alpha[1]);
    gsl_matrix_set(Q, 0, 1,  alpha[2]);
    gsl_matrix_set(Q, 0, 3,  alpha[0]);
    gsl_matrix_set(Q, 1, 0, -alpha[2]);
    gsl_matrix_set(Q, 1, 1,  alpha[1]);
    gsl_matrix_set(Q, 1, 4,  alpha[0]);
    gsl_matrix_set(Q, 2, 2,  alpha[1]);
    gsl_matrix_set(Q, 2, 5,  alpha[0]);
    gsl_matrix_set(Q, 3, 0,  b[0]+alpha[14]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]+alpha[14]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]+alpha[14]);
    gsl_matrix_set(Q, 3, 3, -alpha[1]);
    gsl_matrix_set(Q, 3, 4,  alpha[2]);
    gsl_matrix_set(Q, 4, 3, -alpha[2]);
    gsl_matrix_set(Q, 4, 4, -alpha[1]);
    gsl_matrix_set(Q, 5, 5, -alpha[1]);
}

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized). Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm )
{
    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0]*y[3] + alpha[1]*y[0] + alpha[2]*y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0]*y[4] + alpha[1]*y[1] - alpha[2]*y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0]*y[5] + alpha[1]*y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1]*y[3] + alpha[2]*y[4] + alpha[14]*y[0]
           - alpha[3];
    if(me != 0) f[3] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);


    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1]*y[4] - alpha[2]*y[3] + alpha[14]*y[1]
           - alpha[4];
    if(me != 0) f[4] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);

    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1]*y[5] + alpha[14]*y[2];
    if(me != 0) f[5] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);


    return GSL_SUCCESS;
}

/**
 *  \brief Update the vector field of the state in Normalized-Centered (NC) coordinates. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma)
{
    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0]*y[3] + alpha[1]*y[0] + alpha[2]*y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0]*y[4] + alpha[1]*y[1]-alpha[2]*y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0]*y[5] + alpha[1]*y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1]*y[3] + alpha[2]*y[4] + alpha[14]*y[0]
           + alpha[12];
    if(me != 0) f[3] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);

    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1]*y[4] - alpha[2]*y[3] + alpha[14]*y[1]
           + alpha[13];
    if(me != 0) f[4] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);

    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1]*y[5] + alpha[14]*y[2];
    if(me != 0) f[5] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);

    return GSL_SUCCESS;
}


/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm(const double y[], gsl_matrix *Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma)
{
    double b[6];

    //------------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //------------------------------------------------------------------------------------
    double factor = alpha[5]/pow(gamma,3.0);

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
           + Uij(y, pm, qpm2, mm, factor, 0, 0)
           + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
           + Uij(y, pm, qpm2, mm, factor, 0, 1)
           + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
           + Uij(y, pm, qpm2, mm, factor, 1, 1)
           + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
           + Uij(y, pm, qpm2, mm, factor, 0, 2)
           + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
           + Uij(y, pm, qpm2, mm, factor, 1, 2)
           + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
           + Uij(y, pm, qpm2, mm, factor, 2, 2)
           + Uij(y, ps, qps2, ms, factor, 2, 2);

    //------------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //------------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}


/**
 *  \brief Update the Normalized-Centered variational equation matrix - case of a double **Df, instead of a gsl_matrix. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vfn_stm(const double y[], double **Df, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma)
{
    double b[6];

    //------------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //------------------------------------------------------------------------------------
    double factor = alpha[5]/pow(gamma,3.0);

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
           + Uij(y, pm, qpm2, mm, factor, 0, 0)
           + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
           + Uij(y, pm, qpm2, mm, factor, 0, 1)
           + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
           + Uij(y, pm, qpm2, mm, factor, 1, 1)
           + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
           + Uij(y, pm, qpm2, mm, factor, 0, 2)
           + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
           + Uij(y, pm, qpm2, mm, factor, 1, 2)
           + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
           + Uij(y, pm, qpm2, mm, factor, 2, 2)
           + Uij(y, ps, qps2, ms, factor, 2, 2);

    //------------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //------------------------------------------------------------------------------------
    for(int i = 0; i <6; i++) for(int j = 0;j <6; j++) Df[i][j] = 0.0;
    Df[0][0] = alpha[1];
    Df[0][1] = alpha[2];
    Df[0][3] = alpha[0];
    Df[1][0] =-alpha[2];
    Df[1][1] = alpha[1];
    Df[1][4] = alpha[0];
    Df[2][2] = alpha[1];
    Df[2][5] = alpha[0];
    Df[3][0] = b[0]+alpha[14];
    Df[3][1] = b[1];
    Df[3][2] = b[3];
    Df[4][0] = b[1];
    Df[4][1] = b[2]+alpha[14];
    Df[4][2] = b[4];
    Df[5][0] = b[3];
    Df[5][1] = b[4];
    Df[5][2] = b[5]+alpha[14];
    Df[3][3] =-alpha[1];
    Df[3][4] = alpha[2];
    Df[4][3] =-alpha[2];
    Df[4][4] =-alpha[1];
    Df[5][5] =-alpha[1];

    return GSL_SUCCESS;
}


/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double y[], gsl_matrix *Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm)
{
    double b[6];

    //------------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //------------------------------------------------------------------------------------
    double factor = alpha[5];

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
           + Uij(y, pm, qpm2, mm, factor, 0, 0)
           + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
           + Uij(y, pm, qpm2, mm, factor, 0, 1)
           + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
           + Uij(y, pm, qpm2, mm, factor, 1, 1)
           + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
           + Uij(y, pm, qpm2, mm, factor, 0, 2)
           + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
           + Uij(y, pm, qpm2, mm, factor, 1, 2)
           + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
           + Uij(y, pm, qpm2, mm, factor, 2, 2)
           + Uij(y, ps, qps2, ms, factor, 2, 2);

    //------------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //------------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}

/**
 *  \brief Update the Normalized-Centered linearized variational equation matrix.
 **/
int vfn_stm_lin_trans(const double y[], gsl_matrix *Q, double alpha[],
                     double ps[], double pe[], double pm[],
                     double ms, double me, double mm,
                     double gamma)
{
    double b[6];

    //------------------------------------------------------------------------------------
    // First, we need to take into account the fact that the positions
    // of the primaries need to be translated by the dyn. eq. of the libration point
    //------------------------------------------------------------------------------------
    double xe[3], xm[3], xs[3];
    double re2, rm2, rs2;
    //Factor for the normalized potential
    double factor = alpha[5]/pow(gamma, 3.0);

    //Position of the Earth/Moon/Sun, translated by the dyneq of Li
    //------------------------------------
    xe[0] = pe[0] - y[0];
    xe[1] = pe[1] - y[1];
    xe[2] = pe[2] - y[2];

    xs[0] = ps[0] - y[0];
    xs[1] = ps[1] - y[1];
    xs[2] = ps[2] - y[2];

    xm[0] = pm[0] - y[0];
    xm[1] = pm[1] - y[1];
    xm[2] = pm[2] - y[2];


    //Euclidian distance from the primaries
    //------------------------------------
    re2 = xe[0]*xe[0] + xe[1]*xe[1] + xe[2]*xe[2];
    rm2 = xm[0]*xm[0] + xm[1]*xm[1] + xm[2]*xm[2];
    rs2 = xs[0]*xs[0] + xs[1]*xs[1] + xs[2]*xs[2];


    //------------------------------------------------------------------------------------
    // Linearized variational equations, with translated positions of the primaries
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //------------------------------------------------------------------------------------
    //Q1 matrix
    b[0] =  Uijlin(xe, re2, me, factor, 0, 0)
          + Uijlin(xm, rm2, mm, factor, 0, 0)
          + Uijlin(xs, rs2, ms, factor, 0, 0);

    b[1] =  Uijlin(xe, re2, me, factor, 0, 1)
          + Uijlin(xm, rm2, mm, factor, 0, 1)
          + Uijlin(xs, rs2, ms, factor, 0, 1);

    b[2] =  Uijlin(xe, re2, me, factor, 1, 1)
          + Uijlin(xm, rm2, mm, factor, 1, 1)
          + Uijlin(xs, rs2, ms, factor, 1, 1);

    b[3] =  Uijlin(xe, re2, me, factor, 0, 2)
          + Uijlin(xm, rm2, mm, factor, 0, 2)
          + Uijlin(xs, rs2, ms, factor, 0, 2);

    b[4] =  Uijlin(xe, re2, me, factor, 1, 2)
          + Uijlin(xm, rm2, mm, factor, 1, 2)
          + Uijlin(xs, rs2, ms, factor, 1, 2);

    b[5] =  Uijlin(xe, re2, me, factor, 2, 2)
          + Uijlin(xm, rm2, mm, factor, 2, 2)
          + Uijlin(xs, rs2, ms, factor, 2, 2);


    //------------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //------------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;

}


/**
 *  \brief Computes the second-order derivatives of the linearized potential of one given primary
 **/
double Uijlin(double pc[], double qpc2, double mc, double factor, int i, int j)
{
    if(mc != 0.0)
    {
        if(i == j) return factor*mc/pow(qpc2, 5.0/2)*(3*pc[i]*pc[i] - qpc2);
        else return factor*3*mc/pow(qpc2, 5.0/2)*pc[i]*pc[j];
    }else return 0.0;
}

/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j)
{
    if(mc != 0.0)
    {
        if(i == j) return factor*(3*mc/pow(qpc2,5.0/2)*pow(y[i]-pc[i],2.0) - mc/pow(qpc2,3.0/2));
        else return factor*3*mc/pow(qpc2,5.0/2)*(y[i]-pc[i])*(y[j]-pc[j]);
    }else return 0.0;
}


//----------------------------------------------------------------------------------------
//
// Derivatives for ODE integration in continuation procedures
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. Normalized-Centered coordinates. Variational equations included.
 **/
int qbfbp_vfn_cont(double t, const double y[], double f[], void *params_void)
{
    //Inner variables
    double ys[42], f1[42], f2[42];

    //ys = y[1:42]
    for(int i = 0; i < 42; i++) ys[i] = y[i];

    //Retrieving the parameters
    QBCP_I* qbpi = (QBCP_I *) params_void;

    //VF for model 1 & 2
    qbfbp_vfn_varnonlin(t, ys, f1, &qbpi->model1);
    qbfbp_vfn_varnonlin(t, ys, f2, &qbpi->model2);

    //f = (1-eps)*f1 + eps*f2 for the first 42 components
    for(int i = 0; i < 42; i++) f[i] = (1.0-qbpi->epsilon)*f1[i] + qbpi->epsilon*f2[i];

    //Jacobian matrices
    double **Df1 = dmatrix(0, 5, 0, 5);
    double **Df2 = dmatrix(0, 5, 0, 5);
    qbfbp_Dfn_varnonlin(t, ys, Df1, &qbpi->model1);
    qbfbp_Dfn_varnonlin(t, ys, Df2, &qbpi->model2);

    //Last 6 components: variational equations of epsilon
    for(int i = 0; i < 6; i++)
    {
        f[42+i] = 0.0;
        //omega[i]  = sum(j) Dfij*j
        for(int j = 0; j < 6; j++) f[42+i] += ((1.0-qbpi->epsilon)*Df1[i][j] + qbpi->epsilon*Df2[i][j])*y[42+j];
        //omega[i] += df/depsilon = f2 - f1
        f[42+i] += f2[i] - f1[i];
    }

    //Memory release
    free_dmatrix(Df1, 0, 5, 0, 5);
    free_dmatrix(Df2, 0, 5, 0, 5);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. System coordinates (non-normalized).
 **/
int qbfbp_vf_cont(double t, const double y[], double f[], void *params_void)
{
    double f1[42], f2[42];
    //Retrieving the parameters
    QBCP_I* qbpi = (QBCP_I *) params_void;

    //VF for model 1 & 2
    qbfbp_vf(t, y, f1, &qbpi->model1);
    qbfbp_vf(t, y, f2, &qbpi->model2);

    //f = (1-eps)*f1 + eps*f2
    for(int i = 0; i <42; i++) f[i] = (1.0-qbpi->epsilon)*f1[i] + qbpi->epsilon*f2[i];

    return GSL_SUCCESS;
}


//----------------------------------------------------------------------------------------
//
// Hamiltonians
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian of the QBCP with SYS units and SYS coordinates
 **/
double qbfbp_H(double t, const double y[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double Ps[3];
    eval_array_coef(Ps, t, n, qbp->n_order_fourier, qbp->cs.Ps, 3);
    double Pe[3];
    eval_array_coef(Pe, t, n, qbp->n_order_fourier, qbp->cs.Pe, 3);
    double Pm[3];
    eval_array_coef(Pm, t, n, qbp->n_order_fourier, qbp->cs.Pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //------------------------------------------------------------------------------------
    // Hamiltonian
    //------------------------------------------------------------------------------------
    double H = 0.5*alpha[0]*(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]) + alpha[1]*(y[3]*y[0] + y[4]*y[1] + y[5]*y[2])
             + alpha[2]*(y[3]*y[1] - y[4]*y[0])
             + alpha[3]*y[0] + alpha[4]*y[1]
             - 0.5*alpha[14]*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2])
             - alpha[5]*( me/pow(qPe2, 1.0/2)
                        + mm/pow(qPm2, 1.0/2)
                        + ms/pow(qPs2, 1.0/2) );

    return H;
}


/**
 *  \brief Hamiltonian of the QBCP with EM units and Normalized-Centered coordinates. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
double qbfbp_Hn(double t, const double y[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);


    //------------------------------------------------------------------------------------
    // Hamiltonian
    //------------------------------------------------------------------------------------
    double H = 0.5*alpha[0]*(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])
             + alpha[1]*(y[3]*y[0] + y[4]*y[1] + y[5]*y[2])
             + alpha[2]*(y[3]*y[1] - y[4]*y[0])
             - y[0]*alpha[12]
             - y[1]*alpha[13]
             - 0.5*alpha[14]*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2])
             - alpha[5]/pow(gamma,3.0)*( me/pow(qpe2, 1.0/2) + mm/pow(qpm2, 1.0/2) + ms/pow(qps2, 1.0/2) );
    return H;
}

//----------------------------------------------------------------------------------------
//
// Floquet analysis
//
//----------------------------------------------------------------------------------------
//
/**
 *  \brief Same as the routine qbfbp_vfn + Integration of the matrix Pbar that appears in the Floquet analysis of the hamiltonian of order 2 in the neighborhood of L1,2.
 *         (see nfo2.h)
 **/
int qbfbp_vfn_varlin_trans_P(double t, const double y[], double f[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Memory allocation
    //------------------------------------------------------------------------------------
    gsl_matrix *P    = gsl_matrix_calloc(6,6);
    gsl_matrix *Pd   = gsl_matrix_calloc(6,6);
    gsl_matrix *Q    = gsl_matrix_calloc(6,6);
    gsl_matrix *B    = gsl_matrix_calloc(6,6);

    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //------------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    // STM derivatives, linearized case
    //------------------------------------------------------------------------------------
    vfn_stm_lin_trans(y, Q, alpha, ps, pe, pm, ms, me, mm, gamma);

    //------------------------------------------------------------------------------------
    // Derivatives of P
    //------------------------------------------------------------------------------------
    //P(t) update
    gslc_vectorToMatrix(P, y, 6, 6, 6);

    //Matrix B such that P'(t) = Q(t)P(t) - P(t)*B
    gslc_vectorToMatrix(B, qbp->B, 6, 6, 0);

    //Pd = Q*P
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,  1.0, Q, P, 0.0, Pd);
    //Pd -= P*B = Q*P-P*B = dot(P)
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, P, B, 1.0, Pd);

    //Matrix to Vector
    gslc_matrixToVector(f, Pd, 6, 6 ,6);

    //------------------------------------------------------------------------------------
    //Memory release
    //------------------------------------------------------------------------------------
    gsl_matrix_free(P);
    gsl_matrix_free(Q);
    gsl_matrix_free(B);
    gsl_matrix_free(Pd);

    return GSL_SUCCESS;
}


/**
 *  \brief QBCP vector field in the form of the matrix Q1(t), Q2(t) and Q3(t) defined by:
 *              x' =  Q2^T*x + 2Q3*Y
 *              y' = -2Q1 *x -  Q2*y
 *          Used in the Floquet analysis (see nfo2.h)
 **/
int qbfbp_Q(double t, const double y[], gsl_matrix *Q1, gsl_matrix *Q2, gsl_matrix *Q3, void *params_void)
{
    //------------------------------------------------------------------------------------
    // Memory allocation
    //------------------------------------------------------------------------------------
    gsl_matrix *Q = gsl_matrix_calloc(6,6);

    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    FBPL* qbp = (FBPL *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //------------------------------------------------------------------------------------
    double alpha[noc];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, noc);

    //------------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //------------------------------------------------------------------------------------
    double ps[3];
    eval_array_coef(ps, t, n, qbp->n_order_fourier, qbp->cs.ps, 3);
    double pe[3];
    eval_array_coef(pe, t, n, qbp->n_order_fourier, qbp->cs.pe, 3);
    double pm[3];
    eval_array_coef(pm, t, n, qbp->n_order_fourier, qbp->cs.pm, 3);

    //------------------------------------------------------------------------------------
    // Distances to 2nd power
    //------------------------------------------------------------------------------------
    //    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    //    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    //    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);


    //------------------------------------------------------------------------------------
    // STM derivatives, linearized case
    //------------------------------------------------------------------------------------
    vfn_stm_lin_trans(y, Q, alpha, ps, pe, pm, ms, me, mm, gamma);


    //------------------------------------------------------------------------------------
    // Q1, Q2, Q3 from Q:
    //
    // Q  =  | Q2^T  +2Q3 |
    //       | -2Q1   -Q2 |
    //------------------------------------------------------------------------------------
    gsl_matrix_view Q12 = gsl_matrix_submatrix (Q , 0 , 3 , 3 , 3 );   //Q12
    gsl_matrix_view Q21 = gsl_matrix_submatrix (Q , 3 , 0 , 3 , 3 );   //Q21
    gsl_matrix_view Q22 = gsl_matrix_submatrix (Q , 3 , 3 , 3 , 3 );   //Q22

    //Q1
    gsl_matrix_scale(&Q21.matrix, -0.5);
    gsl_matrix_memcpy(Q1, &Q21.matrix);

    //Q2
    gsl_matrix_scale(&Q22.matrix, -1.0);
    gsl_matrix_memcpy(Q2, &Q22.matrix);

    //Q3
    gsl_matrix_scale(&Q12.matrix, 0.5);
    gsl_matrix_memcpy(Q3, &Q12.matrix);


    //------------------------------------------------------------------------------------
    // Memory release
    //------------------------------------------------------------------------------------
    gsl_matrix_free(Q);
    return GSL_SUCCESS;
}

//----------------------------------------------------------------------------------------
//
// Change of coordinates: SEM <-> IN <-> EM
//
//----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Normalized-Centered coordinates to Earth-Moon coordinates
 **/
void NCtoEM(double t, const double yNC[], double yEM[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];

}

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}


//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------
/**
 *  \brief COC: Normalized-Centered coordinates to system coordinates. Use in priority instead of NCtoEM or NCtoSEM.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;
    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from system coordinates to Normalized-Centered coordinates
 **/
void SYStoNC(double t, const double yEM[], double yNC[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -ySEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -ySEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +ySEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -ySEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -ySEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +ySEM[5]/gamma;
}

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], FBPL *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    eval_array_coef(alpha, t, n, qbp->n_order_fourier, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    ySEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    ySEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    ySEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    ySEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    ySEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    ySEM[5] = +gamma*yNC[5];

}


//-----------------------------------------------------------------------------
// COC: tests & plots
//-----------------------------------------------------------------------------
/**
 *  \brief Plot the Primaries' Three-Body motion of the QBCP.
 **/
void QBTBP_IN()
{
    cout << "----------------------------------------------" << endl;
    cout << "Plot of the QBTBP_IN:                         " << endl;
    cout << "----------------------------------------------" << endl;

    //Param in EM units
    double n  = SEML.us_em.n;
    double ni = SEML.us_em.ni;
    double mu = SEML.us_em.mu_EM;
    double ns = SEML.us_em.ns;
    double as = SEML.us_em.as;
    double ai = SEML.us_em.ai;
    double ms = SEML.us_em.ms;
    double tf = 0.5*SEML.us_em.T;

    //z and Z
    int fftN = OFS_ORDER;
    string filename = "data/qbtbp/";
    Ofsc zt(fftN);
    Ofsc Zt(fftN);
    //zc
    read_ofs_txt(zt, filename+"bjc");
    //Zc
    read_ofs_txt(Zt, filename+"cjc");
    //r & R
    double r1, r2, R1, R2, f1, f2;

    //----------------------------------------------------------------------------------------------------------
    //Plotting devices
    //----------------------------------------------------------------------------------------------------------
    char ch; //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2, *h3, *h4, *h5;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();
    h5 = gnuplot_init();


    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");

    gnuplot_setstyle(h2, (char*)"lines");
    gnuplot_set_xlabel(h2, (char*)"x [-]");
    gnuplot_set_ylabel(h2, (char*)"y [-]");

    gnuplot_setstyle(h3, (char*)"lines");
    gnuplot_set_xlabel(h3, (char*)"x [-]");
    gnuplot_set_ylabel(h3, (char*)"y [-]");

    gnuplot_setstyle(h4, (char*)"lines");
    gnuplot_set_xlabel(h4, (char*)"x [-]");
    gnuplot_set_ylabel(h4, (char*)"y [-]");

    gnuplot_setstyle(h5, (char*)"lines");
    gnuplot_set_xlabel(h5, (char*)"x [-]");
    gnuplot_set_ylabel(h5, (char*)"y [-]");

    //gnuplot_cmd(h4, "set xrange [-1.004:-0.996]");
    //gnuplot_cmd(h4, "set yrange [-0.004:+0.004]");

    //----------------------------------------------------------------------------------------------------------
    //Plot in the IN/SE/EM systems
    //----------------------------------------------------------------------------------------------------------
    int N = 5000;
    double t;
    //Earth, Moon and Sun
    double yINe[6];
    double yINm[6];
    double yINs[6];
    //Velocity to zero (not necessary in this context) + z coord to zero
    for(int i = 2; i<6; i++)
    {
        yINe[i] = 0.0;
        yINm[i] = 0.0;
        yINs[i] = 0.0;
    }

    double yEMe[6];
    double yEMm[6];
    double yEMs[6];
    //Velocity to zero (not necessary in this context) + z coord to zero
    for(int i = 2; i<6; i++)
    {
        yEMe[i] = 0.0;
        yEMm[i] = 0.0;
        yEMs[i] = 0.0;
    }

    double ySEe[6];
    double ySEm[6];
    double ySEs[6];
    //Velocity to zero (not necessary in this context) + z coord to zero
    for(int i = 2; i<6; i++)
    {
        ySEe[i] = 0.0;
        ySEm[i] = 0.0;
        ySEs[i] = 0.0;
    }

    //For plotting
    //------------------------------------
    //EM
    double x_EMe[N];
    double y_EMe[N];
    double x_EMm[N];
    double y_EMm[N];
    double x_EMs[N];
    double y_EMs[N];
    //IN
    double x_INe[N];
    double y_INe[N];
    double x_INm[N];
    double y_INm[N];
    double x_INs[N];
    double y_INs[N];
    //SE
    double x_SEe[N];
    double y_SEe[N];
    double x_SEm[N];
    double y_SEm[N];
    double x_SEs[N];
    double y_SEs[N];

    double x_SEb[N];
    double y_SEb[N];

    //Time vector
    double tvec[N];
    //r/h length scale ratio
    double hr_ratio[N];
    //Loop on time
    for(int i = 0 ; i < N; i++)
    {
        //update time in EM units
        t = tf*i/N;
        //Time vector
        tvec[i] = t;
        //r
        r1 = creal(eval_z(zt, t, n, ni, ai));
        r2 = cimag(eval_z(zt, t, n, ni, ai));
        //R
        R1 = creal(eval_z(Zt, t, n, ns, as));
        R2 = cimag(eval_z(Zt, t, n, ns, as));
        //h
        f1 = R1 - mu*r1;
        f2 = R2 - mu*r2;
        //r/h length scale ratio
        //f = sqrt(f1*f1 + f2*f2);
        hr_ratio[i] = sqrt(r1*r1 + r2*r2)/sqrt(f1*f1 + f2*f2);


        //Inertial
        //----------------------
        //Earth
        yINe[0] = -ms/(1.0+ms)*R1 + mu*r1;
        yINe[1] = -ms/(1.0+ms)*R2 + mu*r2;
        //Moon
        yINm[0] = -ms/(1.0+ms)*R1 - (1-mu)*r1;
        yINm[1] = -ms/(1.0+ms)*R2 - (1-mu)*r2;
        //Sun
        yINs[0] = 1.0/(1.0+ms)*R1;
        yINs[1] = 1.0/(1.0+ms)*R2;

        //Earth
        x_INe[i] = yINe[0];
        y_INe[i] = yINe[1];
        //Moon
        x_INm[i] = yINm[0];
        y_INm[i] = yINm[1];
        //Sun
        x_INs[i] = yINs[0];
        y_INs[i] = yINs[1];

        //EM
        //----------------------
        INtoEM(t, yINe, yEMe, &SEML);
        INtoEM(t, yINm, yEMm, &SEML);
        INtoEM(t, yINs, yEMs, &SEML);

        //Earth
        x_EMe[i] = yEMe[0];
        y_EMe[i] = yEMe[1];
        //Moon
        x_EMm[i] = yEMm[0];
        y_EMm[i] = yEMm[1];
        //Sun
        x_EMs[i] = yEMs[0];
        y_EMs[i] = yEMs[1];


        //SEM
        //----------------------
        //From EM[EM] to SE[SE]
        EMmtoSEMm(t, yEMe, ySEe, &SEML);
        EMmtoSEMm(t, yEMm, ySEm, &SEML);
        EMmtoSEMm(t, yEMs, ySEs, &SEML);

        //Earth
        x_SEe[i] = ySEe[0];
        y_SEe[i] = ySEe[1];
        //Moon
        x_SEm[i] = ySEm[0];
        y_SEm[i] = ySEm[1];
        //Sun
        x_SEs[i] = ySEs[0];
        y_SEs[i] = ySEs[1];

        //Barycenter of Earth and Moon
        //----------------------
        x_SEb[i] = mu*x_SEm[i] + (1-mu)*x_SEe[i];
        y_SEb[i] = mu*y_SEm[i] + (1-mu)*y_SEe[i];

    }

    int color = 1;

    //Sun, Earth & Moon in SE
    //gnuplot_plot_xy(h4, x_SEe, y_SEe, N, (char*)"Earth in SE", "lines",  "1", "2", color++);
    //gnuplot_plot_xy(h4, x_SEm, y_SEm, N, (char*)"Moon  in SE", "lines",  "1", "2",  color++);
    gnuplot_plot_xy(h4, x_SEb, y_SEb, N, (char*)"Bar   in SE", "points", "1", "2",  color++);
    gnuplot_plot_xy(h4, x_SEs, y_SEs, N, (char*)"Sun   in SE", "points", "1", "2", color++);
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //Sun, Earth & Moon in IN
    gnuplot_plot_xy(h1, x_INe, y_INe, N, (char*)"Earth in IN", "lines", "1", "2", color++);
    gnuplot_plot_xy(h1, x_INm, y_INm, N, (char*)"Moon  in IN", "lines", "1", "2", color++);
    gnuplot_plot_xy(h1, x_INs, y_INs, N, (char*)"Sun   in IN", "lines", "1", "2", color++);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //Sun, Earth & Moon in EM
    gnuplot_plot_xy(h2, x_EMe, y_EMe, N, (char*)"Earth in EM", "points", "1", "2", color++);
    gnuplot_plot_xy(h2, x_EMm, y_EMm, N, (char*)"Moon  in EM", "points", "1", "2", color++);
    gnuplot_plot_xy(h2, x_EMs, y_EMs, N, (char*)"Sun   in EM", "lines", "1", "2", color++);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //Earth & Moon in EM
    gnuplot_plot_xy(h3, x_EMe, y_EMe, N, (char*)"Earth in EM", "lines", "1", "2", color++);
    gnuplot_plot_xy(h3, x_EMm, y_EMm, N, (char*)"Moon  in EM", "lines", "1", "2", color++);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //r/h ratio
    gnuplot_plot_xy(h5, tvec, hr_ratio, N, (char*)"r/h ratio", "lines", "1", "2", 1);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);

}

/**
 *  \brief Computes a solution in NC EM coordinates and plot it in NC SEM coordinates in a temporary gnuplot window.
 **/
int odePlot_EMtoSEM(const double y[],
                      double t1,
                      FBPL &qbcp_l,
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
    double yENC[42];
    double ySNC[42];

    for(int i=0; i<42; i++) yENC[i] = y[i];
    double ti = 0;

    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i = 0; i < Npoints; i++)
    {
        ti = i * t1 / Npoints;
        if(i > 0) gsl_odeiv2_driver_apply (d, &t, ti, yENC);
        //ENC to SNC
        NCEMmtoNCSEMm(t, yENC, ySNC, &qbcp_l);
        //Storage
        xc[i] = ySNC[0];
        yc[i] = ySNC[1];
        zc[i] = ySNC[2];
        tc[i] = ti;
    }
    if(t1 <0) d->h = -d->h;
    //gnuplot_setstyle(h1, (char*)plotType.c_str());
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xc, yc, Npoints, title, ls, lt, lw, lc);

    //------------------------------------------
    //In txt files
    //------------------------------------------
    string str_z_nc = (folder+"orbits/"+title+".txt");
    gnuplot_fplot_txyz(tc, xc, yc, zc, Npoints, str_z_nc.c_str());   //orbit

    return GSL_SUCCESS;

}

/**
 *  \brief Computes a solution in NC SEM coordinates and plot it in NC EM coordinates in a temporary gnuplot window.
 **/
int odePlot_SEMtoEM(const double y[],
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
    double ySNC[42];
    double yENC[42];

    for(int i=0; i<42; i++) ySNC[i] = y[i];
    double ti = 0;

    double t = 0.0;
    if(t1 <0) d->h = -d->h;
    for(int i = 0; i < Npoints; i++)
    {
        ti = i * t1 / Npoints;
        if(i > 0) gsl_odeiv2_driver_apply (d, &t, ti, ySNC);
        //SNC to ENC
        NCSEMmtoNCEMm(ti, ySNC, yENC, &SEML);
        //Storage
        xc[i] = yENC[0];
        yc[i] = yENC[1];
        zc[i] = yENC[2];
        tc[i] = ti;
    }
    if(t1 <0) d->h = -d->h;
    //gnuplot_setstyle(h1, (char*)plotType.c_str());
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xc, yc, Npoints, title, ls, lt, lw, lc);


    //------------------------------------------
    //In txt files
    //------------------------------------------
    string str_z_nc = (folder+"orbits/"+title+".txt");
    gnuplot_fplot_txyz(tc, xc, yc, zc, Npoints, str_z_nc.c_str());   //orbit

    return GSL_SUCCESS;

}


/**
 *  \brief Test of the dynamics: SEM to EM system.
 **/
void dynTest_SEMtoEM()
{
    cout << "----------------------------------------------" << endl;
    cout << "Test of the SE dynamics                       " << endl;
    cout << "----------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific) << setprecision(15);

    //-------------------------------------------------
    //Initialization
    //-------------------------------------------------
    //Physical params in EM units
    //--------------------------
    double T_EM  = SEML.us_em.T;     //Period
    //Physical params in SE units
    //--------------------------
    double T_SE = SEML.us_sem.T;     //Period

    //CR3BPs
    CR3BP SE, EM;
    init_cr3bp(&SE, Csts::SUN, Csts::EARTH_AND_MOON);
    init_cr3bp(&EM, Csts::EARTH, Csts::MOON);

    //-------------------------------------------------
    //Plotting devices
    //-------------------------------------------------
    char ch; //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    //-------------------------------------------------
    // SEML focused on the SEM system
    //-------------------------------------------------
    change_coord(SEML, Csts::SEM);

    //-------------------------------------------------
    //Integration tools
    //-------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_varnonlin;
    sys.jacobian = NULL;
    sys.dimension = 42;
    sys.params = &SEML;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());

    //-------------------------------------------------
    //Initial conditions
    //-------------------------------------------------
    double ySNC[42];    //in NC units, momenta
    double ySEM0[42];   //in SEM units, velocities
    double ySEMm[42];   //in SEM units, momenta

    //Li Sun-Earth in SE units
    for(int i =0; i< 42; i++) ySEM0[i] = 0.0;
    switch(SEML.cs_sem.li)
    {
        case 1:
        {
            ySEM0[0] = -SE.l1.position[0];
            break;
        }

        case 2:
        {
            ySEM0[0] = -SE.l2.position[0];
            break;
        }

        case 3:
        {
            ySEM0[0] = -SE.l3.position[0];
            break;
        }
    }


    //To momenta
    SEMvtoSEMm(0.0, ySEM0, ySEMm, &SEML);
    //To NC
    SEMtoNC(0.0, ySEMm, ySNC, &SEML);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);
    //Storing eye(6) into the initial vector
    gslc_matrixToVector(ySNC, Id, 6, 6, 6);
    double tfac = 1.0;


    //-------------------------------------------------
    //Plot Initial conditions (L point)
    //-------------------------------------------------
    gnuplot_plot_xy(h1, &ySNC[0], &ySNC[1], 1, "Initial Li point", "points", "1", "3", 4);

    //-------------------------------------------------
    //Diff Corr
    //-------------------------------------------------
    cout << "----------------------------------------------" << endl;
    cout << "Diff Corr                                     " << endl;
    cout << "----------------------------------------------" << endl;
    differential_correction(ySNC, 0.5*T_SE, 1e-15, d, 42, 0);

    //-------------------------------------------------
    //Plot
    //-------------------------------------------------
    string Li, title;
    if(SEML.cs_sem.li == 2) Li = "L2";
    else Li = "L1";
    title = "DYNEQ of "+Li+" in SEM";
    odePrintGen(ySNC, 42, tfac*T_SE, d, h1, 5000, title.c_str(), "lines", "1", "3", 4, SEML.cs_sem.F_PLOT);

    //-------------------------------------------------
    //Using EM frame
    //-------------------------------------------------
    change_coord(SEML, Csts::EM);

    //-------------------------------
    //SNC to ENC
    //-------------------------------
    double yENC[42];        //in NC units, momenta
    NCSEMmtoNCEMm(0.0, ySNC, yENC, &SEML);
    //Storing eye(6) into the initial vector
    gslc_matrixToVector(yENC, Id, 6, 6, 6);
    //differential_correction(yENC, 0.5*T_EM, 1e-14, d, 0); //NOT NEEDED, ALREADY VERY PRECISE!
    title = "DYNEQ of "+Li+" prop. in EM back in SEM";
    odePlot_EMtoSEM(yENC, tfac*T_EM, SEML, d, h1, 5000, title.c_str(), "lines", "dashed", "3", 2, SEML.cs_sem.F_PLOT);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    gnuplot_close(h1);
    gnuplot_close(h2);
}


/**
 *  \brief Test of the dynamics: EM to SEM system.
 **/
void dynTest_EMtoSEM()
{
    cout << "----------------------------------------------" << endl;
    cout << "Test of the SE dynamics                       " << endl;
    cout << "----------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific) << setprecision(15);

    //-------------------------------------------------
    //Initialization
    //-------------------------------------------------
    //Physical params in EM units
    //--------------------------
    double T_EM  = SEML.us_em.T;     //Period
    //Physical params in SE units
    //--------------------------
    double T_SE = SEML.us_sem.T;     //Period

    //CR3BPs
    CR3BP SE, EM;
    init_cr3bp(&SE, Csts::SUN, Csts::EARTH_AND_MOON);
    init_cr3bp(&EM, Csts::EARTH, Csts::MOON);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //-------------------------------------------------
    //Plotting devices
    //-------------------------------------------------
    char ch; //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    //-------------------------------------------------
    //Using EM frame
    //-------------------------------------------------
    change_coord(SEML, Csts::EM);

    //-------------------------------------------------
    //Integration tools
    //-------------------------------------------------
    //System
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_varnonlin;
    sys.jacobian = NULL;
    sys.dimension = 42;
    sys.params = &SEML;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,
    Config::configManager().G_PREC_HSTART(),
    Config::configManager().G_PREC_ABS(),
    Config::configManager().G_PREC_REL());

    //Initial conditions
    //-------------------------------------------------
    double yEMm[42];    //in EM units, momenta
    double yEM0[42];    //in EM units, velocities
    double yENC[42];    //in NC units, momenta

    //Li Earth-Moon in EM units
    for(int i =0; i< 42; i++) yEM0[i] = 0.0;

    switch(SEML.cs_em.li)
    {
        case 1:
        {
            yEM0[0] = -SEML.cs.cr3bp.l1.position[0];
            break;
        }

        case 2:
        {
            yEM0[0] = -SEML.cs.cr3bp.l2.position[0];
            break;
        }

        case 3:
        {
            yEM0[0] = -SEML.cs.cr3bp.l3.position[0];
            break;
        }
    }


    //Computing the momenta
    EMvtoEMm(0.0, yEM0, yEMm, &SEML);
    //EM to NC
    EMtoNC(0.0, yEMm, yENC, &SEML);
    //Storing eye(6) into the initial vector
    gslc_matrixToVector(yENC, Id, 6, 6, 6);

    //-------------------------------------------------
    //Plot Initial conditions (L point)
    //-------------------------------------------------
    gnuplot_plot_xy(h1, &yENC[0], &yENC[1], 1, "Initial Li point", "points", "1", "3", 4);

    double tfac = 1.0;
    //-------------------------------------------------
    //Diff Corr
    //-------------------------------------------------
    cout << "----------------------------------------------" << endl;
    cout << "Diff Corr                                     " << endl;
    cout << "----------------------------------------------" << endl;
    differential_correction(yENC, 0.5*T_EM, 1e-14, d, 42, 0);
    //-------------------------------------------------
    //Plot
    //-------------------------------------------------
    string Li, title;
    if(SEML.cs_em.li == 2) Li = "L2";
    else Li = "L1";
    title = "DYNEQ of "+Li+" in EM";
    odePrintGen(yENC, 42, tfac*T_EM, d, h1, 5000, title.c_str(), "lines", "1", "3", 1, SEML.cs_em.F_PLOT);


    //-------------------------------------------------
    //Using SEM frame
    //-------------------------------------------------
    change_coord(SEML, Csts::SEM);

    //-----------------
    //ENC to SNC
    //-----------------
    double ySNC[42];    //in NC units, momenta
    NCEMmtoNCSEMm(0.0, yENC, ySNC, &SEML);
    //Storing eye(6) into the initial vector
    gslc_matrixToVector(ySNC, Id, 6, 6, 6);

    //-----------------
    // Plotting
    //-----------------
    title = "DYNEQ of "+Li+" prop. in SEM back in EM";
    odePlot_SEMtoEM(ySNC, tfac*T_SE, d, h1, 5000,  title.c_str(), "lines", "1", "3", 4, SEML.cs_em.F_PLOT);

    //-----------------
    // With Diff Corr in SNC
    //-----------------
    //Diff Corr in SNC
    differential_correction(ySNC, T_SE, 1e-14, d, 42, 0);
    //Change framework for integration
    change_coord(SEML, Csts::EM);
    //SNC to ENC coordinates for IC
    NCSEMmtoNCEMm(0.0, ySNC, yENC, &SEML);
    //Plot
    title = "DYNEQ of "+Li+" diffcor in SEM, prop. in SEM, back in EM";
    odePrintGen(yENC, 42, tfac*T_EM, d, h1, 5000, title.c_str(), "lines", "1", "3", 1, SEML.cs_em.F_PLOT);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    gnuplot_close(h1);
    gnuplot_close(h2);
}
