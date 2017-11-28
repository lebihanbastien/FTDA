#ifndef COC_H_INCLUDED
#define COC_H_INCLUDED


/**
 * \file coc.h
 * \brief Implements the complete change of coordinates for the Normalized-Centered Hamiltonian of the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  Change of coordinates for the hamiltonian equations of the QBCP
 *  This coc allows to get:
 *  - The dynamical equivalent of Li as a fixed point (get rid of order one of the hamiltonian)
 *  - "normal form" of the order 2 of the hamiltonian
 *
 */

//Custom
#include "qbcp.h"
#include "define_env.h"
#include "timec.h"
#include "matrix.h"

#include "init.h"
#include "nfo2.h"

#include <stdio.h>
#include <iostream>

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFTS version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief To TFS format for the whole COC
 **/
void tfs_from_ofs(matrix<Ofsc> &P,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &PCdot,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &Xe,
             vector<Ofsc> &Xm,
             vector<Ofsc> &Xs,
             vector<Ofsc> &V,
             vector<Ofsc> &Vdot,
             Ofsc &ILe,
             Ofsc &ILm,
             Ofsc &ILs);

/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c fbpl.F_COC.
 *      * Vdot = dot(V) is computed
 *      * The products P*C and C*Q are computed in PC and QC
 *      * dot(PC) = PCdot is computed
 *
 *  2. Second step
 *  The vectors Xe[3:5], Xm[3:5], and X[3:5] contains the time-dependent positions of the primaries in the xy-plane.
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[3]^2 + Xc[4]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 **/
void tfts_init_coc(matrix<Ofsc> &P,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &PCdot,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &Xe,
             vector<Ofsc> &Xm,
             vector<Ofsc> &Xs,
             vector<Ofsc> &V,
             vector<Ofsc> &Vdot,
             Ofsc &ILe,
             Ofsc &ILm,
             Ofsc &ILs,
             FBPL& fbpl);
/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c fbpl.F_COC.
 *      * Vdot = dot(V) is computed
 *      * The products P*C and C*Q are computed in PC and QC
 *      * dot(PC) = PCdot is computed
 *
 *  2. Second step
 *  The vectors Xe[2], Xm[2], and X[s] contains the time-dependent positions of the primaries in the xy-plane.
 *  Note that this is the TRANSLATED versions of these positions i.e. they are taking into account the motion
 *  of the periodic orbit, dynamical equivalent of the libration point L1,2.
 *
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[0]^2 + Xc[1]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 *
 **/
void init_coc(matrix<Ofsc> &P,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &PCdot,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &Xe,
             vector<Ofsc> &Xm,
             vector<Ofsc> &Xs,
             vector<Ofsc> &V,
             vector<Ofsc> &Vdot,
             Ofsc &ILe,
             Ofsc &ILm,
             Ofsc &ILs,
             FBPL& fbpl);

//----------------------------------
// COC
//----------------------------------
/**
 *  \brief Apply the change of variables at every orders in zIN/zOut.
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Oftsc> &zIn,
              vector<Oftsc> &zOut);
/**
 *  \brief Apply the change of variables at order m in zIN/zOut, with or without the zero order V
 *         The change of variables is of the form: zOut = PC * zIN (+ V)
 *         Warning: If flag == 1, the zero order is added
 **/
void tfts_applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Oftsc> &zIn,
              vector<Oftsc> &zOut,
              int m,
              int flag);


/**
 *  \brief Apply the change of variables at order m in zIN/zOut, with or without the zero order V
 *         The change of variables is of the form: zOut = PC * zIN (+ V)
 *         Warning: If flag == 1, the zero order is added
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Oftsc> &zIn,
              vector<Oftsc> &zOut,
              int m,
              int flag);

/**
 *  \brief Apply the inverse change of variables at every orders in zIN/zOut.
 *         The change of variables is of the form: zOut = CQ*(zIn - V).
 *
 *   Note: stp is a temporary variable required to store zIn-V
 **/
void applyInvCOC(matrix<Ofsc> &CQ,
                 vector<Ofsc> &V,
                 vector<Oftsc> &zIn,
                 vector<Oftsc> &zOut,
                 vector<Oftsc> &ztp);
//----------------------------------
// dot(COC)
//----------------------------------
/**
 *  \brief Apply the derivative of the change of variables at every orders in zIN/zOut.
 *         The change of variables is of the form: zdot = PC zh + PC zhdot + Vdot
 **/
void applyDotCOC(matrix<Ofsc> &PC,
                 matrix<Ofsc> &PCdot,
                 vector<Ofsc> &Vdot,
                 vector<Oftsc> &zh,
                 vector<Oftsc> &zhdot,
                 vector<Oftsc> &zdot);
/**
 *  \brief Apply the derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zdot = PC zh + PC zhdot + Vdot
 **/
void applyDotCOC(matrix<Ofsc> &PC,
                 matrix<Ofsc> &PCdot,
                 vector<Ofsc> &Vdot,
                 vector<Oftsc> &zh,
                 vector<Oftsc> &zhdot,
                 vector<Oftsc> &zdot,
                 int m);

/**
 *  \brief Apply the inverse derivative of the change of variables at every orders m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 *
 *         Warning: zhdot is used as temporary variable within the routine
 *         However, the final result is good.
 **/
void applyInvDotCOC(matrix<Ofsc> &CQ,
                    matrix<Ofsc> &PCdot,
                    vector<Ofsc> &Vdot,
                    vector<Oftsc> &zh,
                    vector<Oftsc> &zdot,
                    vector<Oftsc> &zhdot,
                    vector<Oftsc> &ztp_1,
                    vector<Oftsc> &ztp_2);

/**
 *  \brief Apply the inverse derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 **/
void applyInvDotCOC(matrix<Ofsc> &CQ,
                    matrix<Ofsc> &PCdot,
                    vector<Ofsc> &Vdot,
                    vector<Oftsc> &zh,
                    vector<Oftsc> &zdot,
                    vector<Oftsc> &zhdot,
                    vector<Oftsc> &ztp_1,
                    vector<Oftsc> &ztp_2,
                    int m);

/**
 *  \brief Apply the inverse derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 **/
void tfts_apply_inv_coc_der(matrix<Ofsc> &CQ,
                    matrix<Ofsc> &PCdot,
                    vector<Ofsc> &Vdot,
                    vector<Oftsc> &zh,
                    vector<Oftsc> &zdot,
                    vector<Oftsc> &zhdot,
                    vector<Oftsc> &ztp_1,
                    vector<Oftsc> &ztp_2,
                    int m);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFS version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c fbpl.F_COC.
 *
 *  2. Second step
 *  The vectors Xe[2], Xm[2], and X[s] contains the time-dependent positions of the primaries in the xy-plane.
 *  Note that this is the TRANSLATED versions of these positions i.e. they are taking into account the motion
 *  of the periodic orbit, dynamical equivalent of the libration point L1,2.
 *
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[0]^2 + Xc[1]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 *  NOTE: This routine is a simplified version of the routine of the same name, used for OFTS computations. It is only used for OFS version of the coc
 *  which itself is only used as a test version of the OFTS coc.
 *
 **/
void initCOC_OFS(matrix<Ofsc> &P,
             matrix<Ofsc> &Q,
             vector<Ofsc> &V,
             vector<Ofsc> &Xe,
             vector<Ofsc> &Xm,
             vector<Ofsc> &Xs,
             Ofsc &ILe,
             Ofsc &ILm,
             Ofsc &ILs,
             FBPL& fbpl);

/**
 *  \brief Apply the change of variables in zIN/zOut (OFS version).
 *       The change of variables is of the form: zOut = (P*C) zIN + V
 *
 *   This routine shows that using P as it is, without computing PC in the preprocess,
 *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX in the src code).
 **/
void applyCOC_OFS(matrix<Ofsc> &P, vector<Ofsc> &V, vector<Ofsc> &zIn, vector<Ofsc> &zOut);

/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC_OFS(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Ofsc> &zIn,
              vector<Ofsc> &zOut);

/**
 *  \brief Apply the change of variables (without the order zero) in zIN/zOut (OFS version).
 *       The change of variables is of the form: zOut = (P*C) zIN
 *
 *   As applyCOC_OFS, this routine shows that using P as it is, without computing PC in the preprocess,
 *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX in the src code).
 **/
void applyModifiedCOC_OFS(matrix<Ofsc> &P, vector<Ofsc> &V, vector<Ofsc> &zIn, vector<Ofsc> &zOut);

/**
 *  \brief Apply the inverse change of variables in zIN/zOut (OFS version).
 *         The change of variables is of the form: zOut = CQ*(zIn - V).
 *
 *   This routine shows that using Q as it is, without computing CQ in the preprocess,
 *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX, EUX in the src code).
 **/
void applyInvCOC_OFS(matrix<Ofsc>& Q, vector<Ofsc>& V, vector<Ofsc>& zIn, vector<Ofsc>& zOut);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Scalar version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Apply the change of variables in zIN/zOut (scalar version).
 *       The change of variables is of the form: zOut = (P*C) zIN + V
 */
void applyCOC(matrix<Ofsc> PC, vector<Ofsc> V, vector<cdouble> &zIn, vector<cdouble> &zOut,  double t);


/**
 *  \brief Apply the derivative change of variables at every order in zIN/zOut (scalar version).
 */
void applyDotCOC(matrix<Ofsc> PC,
                 matrix<Ofsc> PCdot,
                 vector<Ofsc> Vdot,
                 vector<cdouble> &zh,
                 vector<cdouble> &zhdot,
                 vector<cdouble> &zdot,
                 double t);

/**
 *  \brief Apply the inverse of the derivative change of variables in zIN/zOut (scalar version).
 */
void applyInvDotCOC(matrix<Ofsc> CQ,
                    matrix<Ofsc> PCdot,
                    vector<Ofsc> Vdot,
                    vector<cdouble> &zh,
                    vector<cdouble> &zdot,
                    vector<cdouble> &zhdot,
                    vector<cdouble> &ztp,
                    double t);


/**
 *  \brief Apply the inverse of the change of variables at every order in zIN/zOut (scalar version).
 */
void applyInvCOC(matrix<Ofsc> CQ, vector<Ofsc> V, vector<cdouble> &zIn, vector<cdouble> &zOut, vector<cdouble> &ztp, double t);




//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Tests
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test routine for the change of coordinates: comparison between an object x and COC^{-1}(COC(x)). Tested on OFS and OFTS objects.
 **/
void testCOC();

/**
 *  \brief Test routine for the derivatives of change of coordinates. Tested OFTS objects.
 **/
void testDotCOC();

/**
 *  \brief Test routine for the change of coordinates on a trajectory integrated on a full period.
 *
 *      The same initial conditions are integrated both in NC and TFC coordinates on a full period T, and the results are compared.
 *      Note that the vector field in TFC makes use of the scalar version of the COC, so this test only proves that the scalar COC is self-consistent
 *      but no more than that.
 **/
void testIntCOC();

/**
 *  \brief Gives the vector field f[12] in TFC coordinates from the state y[12] in TFC coordinates.
 *
 *  The process is the following:
 *  - The state zh is given in the form y[i] = re(zh[i%2]), y[i+1] = Im(zh[i%2]) (real and imag part of the state are stored in y)
 *  - zh is build from y and the scalar version of the COC is applied : z = COC(zh)
 *  - The classical vector field of the QBCP is applied on z: zdot = qbfbp_vfn_novar(z)
 *  - The inverse of the derivatives of the coc is applied on the corresponding vector field: f = invCOCdot(zdot)
 *
 *  This routine is essentially used in testIntCOC to test the COC on a trajectory integrated on a full period.
 **/
int qbfbp_Fh(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Finds the first positive root t0 of the equation Pij(nt) = 0. Note that ii and jj must be given between 0 and 5
 **/
double pij(int ii, int jj);

/**
 *  \brief Plot Pij(nt)
 **/
void pij_plot(int ii, int jj, gnuplot_ctrl  *h1);

/**
 *  \brief Plot Pij(nt)
 **/
void pij_plot(gnuplot_ctrl  *h1);

/**
 *  \brief Plot the coefficient of qbp
 **/
void coeff_plot(gnuplot_ctrl  *h1, FBPL* qbp);

/**
 *  \brief Plot the potential of the primaries when the s/c is fixed at the Lagrange point.
 **/
void potential_plot(gnuplot_ctrl  *h1, FBPL* qbp);


 /**
 * \struct COC
 * \brief Structure to store the change of coordinate elements used in the vector field routine qbfbp_Fh.
 **/
typedef struct COC COC;
struct COC
{
    //PC=P(theta)*C
    matrix<Ofsc>* PC;
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc>* PCdot;
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc>* CQ;
    //V
    vector<Ofsc>* V;
    //Vdot=dot(V)
    vector<Ofsc>* Vdot;
    //FBPL
    FBPL* fbpl;
    vector<cdouble>* zh;  //current state in TFC variables
    vector<cdouble>* z;   //current state in Nc  variables
    vector<cdouble>* zd;  //temporary variable used in COC
    vector<cdouble>* zhd; //temporary variable used in COC
    vector<cdouble>* ztp; //temporary variable used in COC
};

#endif // COC_H_INCLUDED
