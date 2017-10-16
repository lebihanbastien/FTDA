#ifndef PM_H_INCLUDED
#define PM_H_INCLUDED

#include "pmcoc.h"
#include "matrix.h"

#include "odezero.h"
#include "init.h"

//Gnuplot
extern "C" {
#include "gnuplot_i.h"
}

/**
 * \file pm.h
 * \brief Implements the parameterization method for both autonomous and non-autonomous Hamiltonian vector fields of the Sun-Earth-Moon problem. Deprecated.
           Use routines of pmt.h instead.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  More precisely, the parameterization method is applied to compute a parameterization of the central manifold of the Earth-Moon libration point L1,2:
 *
 *   - The parameterization is given as Fourier-Taylor expansions in the case of non-autonomous Hamiltonian vector fields (QBCP, BCP).
 *   - The parameterization is given as pure    Taylor expansions in the case of     autonomous Hamiltonian vector fields (CRTBP).
 */


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Main routine
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Compute the parameterization of the central manifold of the dynamical equivalent of the Earth-Moon libration points L1,2, up to a given order.
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_PM...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.F_PM):
 *
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 *
 **/
void pm(int OutputEachOrder, int Output);

/**
 *  \brief Same routine as pm, with additional test.
 *
 *  Invariance Equation test:
 *
 *  - IE_ots: checks that the invariance equation is satisfied by W(s,t) at each order.
 *  - IEh_ots: checks that the invariance equation is satisfied by Wh(s,t) at each order.
 *  - xi5_ots; consistency test of the cohomological equations in the coordinates s (CCM).
 *  - Wh6_ots; consistency test of the cohomological equations in the coordinates Wh (TFC).
 *  - D1_ots & D2_ots: splitted consistency test  of the cohomological equations in the coordinates Wh.
 *    More precisly:
 *
 *      D1 = T_FWh - DF(0)Wh - FWh, with T_FWh the complete vector field, FWh the vector field without the contribution of the last order in Wh
 *      D2 = T_DWfh - LaLs*dWh[p]/dx[j] - L*fh - DWfh, with DWfh the complete derivatives of the vector field and DWfh the derivatives of the
 *                                                     vector field without the contribution of the last order in Wh
 *
 */
void pmTested();

/**
 *  \brief Compute the parameterization of the central manifold of the dynamical equivalent of the Earth-Moon libration points L1,2, up to a given order.
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_PM...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_PM):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 **/
void pm_normalform(int OutputEachOrder, int Output);

/**
 *  \brief Compute the 6D parameterization of the transit/non-transit trajectories about the Earth-Moon libration points L1,2, up to a given order.
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_PM...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_PM):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 *
 *  Notes:
 *     The routine apply a mixed style parameterization in order to separate as many submanifolds as possible. Here is a list, with the corresponding sets of indices
 *     used for the mixed style (see Haro 2014, section 2.2.3).
 *
 *            set I          |  submanifold described by si = 0, for all i in I
 *          ----------------------------------------------------------------
 *           {1,2,3,4,6}     |  unstable manifold
 *           {1,2,3,4,5}     |  stable manifold
 *           {1,2,3,4}       |  hyperbolic normal part
 *           {2,4,5,6}       |  planar family
 *           {2,4,6}         |  unstable manifold of the planar family
 *           {2,4,5}         |  stable manifold of the planar family
 *           {2,4}           |  hyperbolic normal part of the planar family
 *           {1,3,5,6}       |  vertical family
 *           {1,3,6}         |  unstable manifold of the vertical family
 *           {1,3,5}         |  stable manifold of the vertical family
 *           {1,3}           |  hyperbolic normal part of the vertical family
 *           {5,6}           |  center manifold
 *           {6}             |  center-unstable manifold
 *           {5}             |  center-stable manifold
 *
 **/
void pm_mixedstyle(int OutputEachOrder, int Output);


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFTS version of the recurrence
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Applying the Normalized-Centered Vector Field on a Ofts object at order m
 *        WARNING: This routine assumes that the potential has been udpated up to order m
 **/
void applyVF(vector<Ofsc>& alpha,
             vector<Oftsc>& W,
             vector<Oftsc>& FW,
             vector<Oftsc>& Unx,  //vector of size OFTS_ORDER
             vector<Oftsc>& Uny,  //vector of size OFTS_ORDER
             vector<Oftsc>& Unz,  //vector of size OFTS_ORDER
             int m);

/**
 * \brief Applying the Normalized-Centered Vector Field on a Ofts object at order m
 *        WARNING: This routine assumes that the potential has been udpated up to order m
 *        Case with Earth, Moon and Sun potential separated.
 **/
void applyVF(vector<Ofsc>  &alpha,
             vector<Oftsc> &W,
             vector<Oftsc> &FW,
             vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
             vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
             int m);

/**
 *   \brief Initializes the Legendre-derived expansions of one given primary at a given order m.
 *          Note that it is better to reset the order m of En[1] and Enx[2] in this routine.
 *          Since they are linear in the variable W[i], we always have the good value in it.
 **/
void initLegPrim(vector<Oftsc> &Wt,   //vector of size 6
                 vector<Oftsc> &En,   //vector of size OFTS_ORDER
                 vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                 vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                 vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                 vector<Ofsc> &Xe,    //modified position of the Earth
                 Ofsc const &ILe,     //modified inverse orbit radius of the Earth
                 double Ke,
                 Ofsc &AUX,
                 Ofsc &BUX,
                 Ofsc &CUX,
                 Ofsc &DUX,
                 int m);


/**
 *   \brief Initializes the Legendre-derived expansions of the Primaries at a given order m.
 *          Note that it is better to reset the order m of En[1] and Enx[2] in this routine.
 *          Since they are linear in the variable W[i], we always have the good value in it.
 **/
void initLegendrePoly( vector<Oftsc> &Wt,   //vector of size 6
                       vector<Oftsc> &En,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                       vector<Ofsc>  &Xe,   //modified position of the Earth
                       vector<Ofsc>  &Xm,   //modified position of the Moon
                       vector<Ofsc>  &Xs,   //modified position of the Sun
                       Ofsc &ILe,           //modified inverse orbit radius of the Earth
                       Ofsc &ILm,           //modified inverse orbit radius of the Moon
                       Ofsc &ILs,           //modified inverse orbit radius of the Sun
                       QBCP_L& qbcp_l,
                       int m);

/**
 *   \brief Apply the Legendre-derived recurrence for the potential of one given primary.
 *          No reset here, since we may add other terms afterwards.
 **/
void applyLegRec(Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                 vector<Oftsc> &En,   //vector of size OFTS_ORDER
                 Ofsc &ILe,           //modified inverse orbit radius of the Primary
                 Oftsc &xse,          //=xe*xt+ye*yt @order m
                 Ofsc &AUX,           //spare OFS object
                 Ofsc &temp,          //spare OFS object
                 int m,                  //current order of the W expansion
                 int k);                 //current En[k] that needs to be updated

/**
 *   \brief Apply the Legendre-derived recurrence for the derivatives of the potential of one given primary.
 *          No reset here, since we may add other terms afterwards.
 **/
void applyLegRecDer(vector<Oftsc> &Wt,   //tilde state vector (xt, yt,..)
                    Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                    vector<Oftsc> &En,   //vector of size OFTS_ORDER
                    vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                    vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                    vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                    vector<Ofsc>  &Xe,   //modified position of the primary
                    Ofsc &ILe,           //modified inverse orbit radius of the Primary
                    Oftsc &xse,          //=xe*xt+ye*yt @order m
                    Ofsc &AUX,           //spare OFS object
                    Ofsc &BUX,           //spare OFS object
                    Ofsc &temp,          //spare OFS object
                    int m,                  //current order of the W expansion
                    int k);                 //current En[k] that needs to be updated

/**
 *   \brief Update the order m of the kth coefficient in the Legendre-derived objects En, Mn, etc (k < m),
 *          Provided that the (k-2t)h and (k-1)th coefficient have been updated.
 **/
void updateLegPoly(vector<Oftsc> &Wt,   //vector of size 6: tilde state vector (xt, yt,..)
                   Oftsc &Rho2,         //Ofts containing Wt[0]^2+Wt[1]^2+Wt[2]^2
                   vector<Oftsc> &En,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                   vector<Ofsc> &Xe,    //modified position of the Earth
                   vector<Ofsc> &Xm,    //modified position of the Moon
                   vector<Ofsc> &Xs,    //modified position of the Sun
                   Ofsc &ILe,           //modified inverse orbit radius of the Earth
                   Ofsc &ILm,           //modified inverse orbit radius of the Moon
                   Ofsc &ILs,           //modified inverse orbit radius of the Sun
                   Oftsc &xse,          //=xe*xt+ye*yt
                   Oftsc &xsm,          //=xm*xt+ym*yt
                   Oftsc &xss,          //=xs*xt+ys*yt
                   Ofsc &AUX,           //Temporary variables
                   Ofsc &BUX,           //Temporary variables
                   Ofsc &temp,          //Temporary variables
                   int m,                  //current order of the W expansion
                   int k);                 //current En[k] (and Enx[k+1]) that needs to be updated


/**
 *   \brief Update the whole Legendre-derived objects at order m > 1.
 **/
void updatePotentialExpansion(vector<Oftsc> &Wt,    //vector of size 6
                              Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                              vector<Oftsc> &En,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                              vector<Ofsc> &Xe,    //modified position of the Earth
                              vector<Ofsc> &Xm,    //modified position of the Moon
                              vector<Ofsc> &Xs,    //modified position of the Sun
                              Ofsc &ILe,           //modified inverse orbit radius of the Earth
                              Ofsc &ILm,           //modified inverse orbit radius of the Moon
                              Ofsc &ILs,           //modified inverse orbit radius of the Sun
                              Oftsc &xse,          //=xe*xt+ye*yt
                              Oftsc &xsm,          //=xm*xt+ym*yt
                              Oftsc &xss,          //=xs*xt+ys*yt
                              QBCP_L& qbcp_l,        //current QBCP
                              int m);                 //targeted order of the W expansion


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluate the potential (for testing)
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test the validity of the expansions of the potential (En, Sn, Mn, Un) and the vector field (FW).
 **/
void pmTestPrec(double smax);

/**
 *   \brief Evaluate the different potentials at order m at time n*t: expansion vs real implementation
 **/
void evaluatePotential(cdouble X[],
                       cdouble& Ez,
                       cdouble& Mz,
                       cdouble& Sz,
                       cdouble& Uz,
                       cdouble& Ez2,
                       cdouble& Mz2,
                       cdouble& Sz2,
                       cdouble& Uz2,
                       vector<Oftsc> &En,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                       vector<Oftsc> &W,   //vector of size OFTS_ORDER
                       vector<Ofsc> &Xe,    //modified position of the Earth
                       vector<Ofsc> &Xm,    //modified position of the Moon
                       vector<Ofsc> &Xs,    //modified position of the Sun
                       double t,
                       QBCP_L& qbcp_l,
                       int m);

/**
 *   \brief Evaluate the derivatives of the potentials at order m at time n*t: expansion vs real implementation
 **/
void evaluatePotentialDer(cdouble X[],
                          cdouble& Ex, cdouble& Mx, cdouble& Sx, cdouble& Ux,
                          cdouble& Ey, cdouble& My, cdouble& Sy, cdouble& Uy,
                          cdouble& Ez, cdouble& Mz, cdouble& Sz, cdouble& Uz,
                          cdouble& Ex2, cdouble& Mx2, cdouble& Sx2, cdouble& Ux2,
                          cdouble& Ey2, cdouble& My2, cdouble& Sy2, cdouble& Uy2,
                          cdouble& Ez2, cdouble& Mz2, cdouble& Sz2, cdouble& Uz2,
                          vector<Oftsc> &Enx, vector<Oftsc> &Mnx, vector<Oftsc> &Snx, vector<Oftsc> &Unx,
                          vector<Oftsc> &Eny, vector<Oftsc> &Mny, vector<Oftsc> &Sny, vector<Oftsc> &Uny,
                          vector<Oftsc> &Enz, vector<Oftsc> &Mnz, vector<Oftsc> &Snz, vector<Oftsc> &Unz,
                          vector<Oftsc> &W,   vector<Ofsc> &Xe,   vector<Ofsc> &Xm,   vector<Ofsc> &Xs,
                          double t,
                          QBCP_L& qbcp_l,
                          int m);


/**
 *   \brief Evaluate the the vector field at order m at time n*t: expansion vs real implementation
 **/
void evaluateVectorField(cdouble s0[],
                         cdouble VF1[],
                         cdouble VF2[],
                         vector<Oftsc> &W,    //vector of size OFTS_ORDER
                         vector<Oftsc> &FW,   //vector of size OFTS_ORDER
                         double t,
                         QBCP_L& qbcp_l,
                         int m);

/**
 *  \brief Earth potential = (1-mu)/(gamma^3*||qpe||)
 **/
cdouble EarthPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xe, QBCP_L& qbcp_l);

/**
 *  \brief Moon potential = mu/(gamma^3*||qpm||)
 **/
cdouble MoonPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xm, QBCP_L& qbcp_l);

/**
 *  \brief Sun potential = ms/(gamma^3*||qps||)
 **/
cdouble SunPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xs, QBCP_L& qbcp_l);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//         Inner routines
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Prints a given vector W of type \c T in a txt files of the form "filename+i.txt", with i = 0, length(W)-1. The type \c T must have an implementation of the operator \c <<.
 **/
template <typename T> void  vector_fprinf(T& W, string filename);


//---------------------------------------------------------------------------------------------------------------------------------------
//          Reading
//---------------------------------------------------------------------------------------------------------------------------------------
void readVOFTS_txt(vector<Oftsc> &W, string filename, int fftN);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFS version of the recurrence
//
//---------------------------------------------------------------------------------------------------------------------------------------
//Test the recurrence routines
void testLegendreRecurrence_OFS();

//Init the recurrence on Legendre-derived objects
void initLegendreRecurrence_OFS( vector<Ofsc> &W,    //vector of size 6
                                 Ofsc &Rho2,          //Ofts containing W[0]^2+W[1]^2+W[2]^2
                                 vector<Ofsc> &En,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mn,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Sn,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Enx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Eny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Enz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mnx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mnz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Snx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Sny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Snz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Xe,    //modified position of the Earth
                                 vector<Ofsc> &Xm,    //modified position of the Moon
                                 vector<Ofsc> &Xs,    //modified position of the Sun
                                 Ofsc &ILe,           //modified inverse orbit radius of the Earth
                                 Ofsc &ILm,           //modified inverse orbit radius of the Moon
                                 Ofsc &ILs,           //modified inverse orbit radius of the Sun
                                 QBCP_L& qbcp_l);


//Apply recurrence on Legendre-derived objects for m > 1
//Update the order m of the kth coefficient in En, Mn, etc (k < m)
//Provided that the (k-2t)h and (k-1)th coefficient have been updated
void applyLegendreRecurrence_OFS(vector<Ofsc> &W,    //vector of size 6
                                 Ofsc &Rho2,          //Ofts containing W[0]^2+W[1]^2+W[2]^2
                                 vector<Ofsc> &En,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mn,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Sn,   //vector of size OFTS_ORDER
                                 vector<Ofsc> &Enx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Eny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Enz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mnx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Mnz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Snx,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Sny,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Snz,  //vector of size OFTS_ORDER
                                 vector<Ofsc> &Xe,    //modified position of the Earth
                                 vector<Ofsc> &Xm,    //modified position of the Moon
                                 vector<Ofsc> &Xs,    //modified position of the Sun
                                 Ofsc &ILe,           //modified inverse orbit radius of the Earth
                                 Ofsc &ILm,           //modified inverse orbit radius of the Moon
                                 Ofsc &ILs,           //modified inverse orbit radius of the Sun
                                 QBCP_L& qbcp_l,       //QBFPB on a given Li
                                 int k);


//Update the potential up to order nPot in the expansions
void updatePotentialExpansion_OFS(vector<Ofsc> &W,    //vector of size 6
                             Ofsc &Rho2,          //Ofts containing W[0]^2+W[1]^2+W[2]^2
                             vector<Ofsc> &En,   //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Mn,   //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Sn,   //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Enx,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Eny,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Mnx,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Mny,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Mnz,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Snx,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Sny,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Snz,  //vector of size OFTS_ORDER+1
                             vector<Ofsc> &Xe,    //modified position of the Earth
                             vector<Ofsc> &Xm,    //modified position of the Moon
                             vector<Ofsc> &Xs,    //modified position of the Sun
                             Ofsc &ILe,           //modified inverse orbit radius of the Earth
                             Ofsc &ILm,           //modified inverse orbit radius of the Moon
                             Ofsc &ILs,           //modified inverse orbit radius of the Sun
                             QBCP_L& qbcp_l);     //QBFPB on a given Li

//---------------------------------------------------------------------------------------------------------------------------------------
//          FTDA version
//---------------------------------------------------------------------------------------------------------------------------------------


#endif // PM_H_INCLUDED
