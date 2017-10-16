#ifndef PMT_H_INCLUDED
#define PMT_H_INCLUDED

/**
 * \file pmt.h
 * \brief Implements the parameterization method for both autonomous and non-autonomous Hamiltonian vector fields of the Sun-Earth-Moon problem.
 *        uses time form of the Fourier-Taylor series to perform algebraic operations.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  More precisely, the parameterization method is applied to compute a parameterization of the central manifold of the Earth-Moon libration point L1,2:
 *
 *   - The parameterization is given as Fourier-Taylor expansions in the case of non-autonomous Hamiltonian vector fields (QBCP, BCP).
 *   - The parameterization is given as pure    Taylor expansions in the case of     autonomous Hamiltonian vector fields (CRTBP).
 */

#include "pmcoc.h"
#include "matrix.h"
#include "init.h"

//---------------------------------------------------------------------------------------------------------------------------------------
//         Main routine
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Compute the parameterization of the central manifold of the dynamical equivalent of the Earth-Moon libration points L1,2, up to a given order.
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *  \param pms: parameterization style (graph, normal form or mixed)
 *  \param manType: type of manifold (center, center-stable, center-unstable or center-hyperbolic).
 *
 *   Inputs:
 *      - The FBPL structure SEML, which contains the addresses of the data folders (F_COC, F_PMS...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_PMS, in binary format):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *      - FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *      - DWf(s,t), the product DW * fh
 *      - Wdotc(s,t), the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 **/
void pmt(int OutputEachOrder, int Output, int pms, int manType);

//---------------------------------------------------------------------------------------------------------------------------------------
//         Init routines
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Initialization of the vector fields (alpha coefficients) used in pm.
 **/
void tfts_initVF(vector<Ofsc> &alpha);

/**
 *  \brief Initialization of various matrices that helps build the order one of the parameterization. Namely:
 *
 *      | iw1 0    0    0    0   0   |
 *      | 0   w2   0    0    0   0   |
 *      | 0   0    iw3  0    0   0   |
 * DB = | 0   0    0   -iw1  0   0   |
 *      | 0   0    0    0   -w2  0   |
 *      | 0   0    0    0    0  -iw3 |
 *
 * DF0 = DB in matrix<Ofsc> format
 *
 *     |  iw1 0    0    0    0   0  |
 *     |  0   0    0    0    w2  0  |
 *     |  0   iw3  0    0    0   0  |
 * H = |  0   0   -iw1  0    0   0  | = (L N)
 *     |  0   0    0    0    0  -w2 |
 *     |  0   0    0   -iw3  0   0  |
 *
 * H2 = H in matrix<Ofsc> format
 *
 *     |  iw1 0    0    0    |
 *     |  0   0    0    0    |
 *     |  0   iw3  0    0    |
 * L = |  0   0   -iw1  0    |
 *     |  0   0    0    0    |
 *     |  0   0    0   -iw3  |
 *
 * L2 = L in matrix<Ofsc> format
 *
 *     |  0   0  |
 *     |  w2  0  |
 *     |  0   0  |
 * N = |  0   0  |
 *     |  0  -w2 |
 *     |  0   0  |
 *
 *        |1/iw1  0    0      0     0     0    |
 *        | 0     0   1/iw3   0     0     0    |
 *        | 0     0    0    -1/iw1  0     0    |
 * Hinv = | 0     0    0      0     0  -1/iw3  |
 *        | 0    1/w2  0      0     0     0    |
 *        | 0     0    0      0   -1/w2   0    |
 *
 * with (w1, w2, w3) the frequency associated to the xy elliptic motion, hyperbolic motion and z elliptic motion, respectively.
 *
 *  WARNING: all these objects are initialized in OFS format, and NOT in TFS format, since they are only used in pure OFS computations.
 **/
void tfts_initOrderOne(gsl_matrix_complex *DB,
                       gsl_matrix_complex *H,
                       gsl_matrix_complex *Hinv,
                       gsl_matrix_complex *La,
                       matrix<Ofsc>& DF0,
                       matrix<Ofsc>& H2,
                       matrix<Ofsc>& Hinv2,
                       int manType);
/**
 *  \brief Initialization of the PM objects (W, Wh, and fh). TFS formatting is included.
 **/
void tfts_initPM(vector<Oftsc>& W,
                 vector<Oftsc>& Wh,
                 vector<Oftsc>& fh,
                 matrix<Ofsc>& PC,
                 vector<Ofsc>& V,
                 gsl_matrix_complex *L,
                 gsl_matrix_complex *La,
                 int manType);

//---------------------------------------------------------------------------------------------------------------------------------------
//         Cohomological equations
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Resolution of the cohomological equations at order m.
 *
 *         The right-hand side of the cohomological equation are in the expansion xi. The expansions eta and fh are updated.
 **/
void cohomEq(vector<Oftsc>& eta,
             vector<Oftsc>& xi,
             vector<Oftsc>& fh,
             gsl_matrix_complex *La,
             int m,
             int pms,
             double threshold,
             int **VIs,
             int  *VIn,
             int ims);

//---------------------------------------------------------------------------------------------------------------------------------------
//         Recurrence
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Adds the contribution of the terms of order m to the order m of Un. These terms are part of the potential of the primary with position Xe and factor Ke.
 **/
void updateDerPotprim_zero(vector<Oftsc>& W,        //parameterization in NC coordinates
                           vector<Ofsc>& Xe,        //position of the primary
                           Oftsc& E1, Oftsc& E2,    //intermediate steps
                           Oftsc& Et1, Oftsc& Et2,  //intermediate steps
                           double Ke,               //primary factor
                           vector<Oftsc>& Un,       //potential to update
                           int m);                  //order

/**
 *  \brief Adds the contribution of the terms of order < m to the order m of Un. These terms are part of the the potential of the primary with position Xe and factor Ke.
 **/
void updateDerPotprim(vector<Oftsc>& W,     //parameterization in NC coordinates
                      vector<Ofsc>& Xe,     //position of the primary
                      Oftsc& E1, Oftsc& E2, //intermediate steps
                      double Ke,            //primary factor
                      vector<Oftsc>& Un,    //potential to update
                      int m);               //order

/**
 *  \brief Adds the contribution of the terms of order m to the order m of the potential Un.
 **/
void updateDerPot_zero(vector<Oftsc>& W,        //parameterization in NC coordinates
                       vector<Ofsc>& Xe,        //Earth position
                       vector<Ofsc>& Xm,        //Moon position
                       vector<Ofsc>& Xs,        //Sun position
                       vector<Oftsc>& PrimPt,   //Intermediate steps: contr. of order < m to order m
                       vector<Oftsc>& PrimPt2,  //Intermediate steps: contr. of order = m to order m
                       double Ke,               //Earth factor
                       double Km,               //Moon factor
                       double Ks,               //Sun factor
                       vector<Oftsc>& Un,       //potential to update
                       int m);                  //order

/**
 *  \brief Adds the contribution of the terms of order < m to the order m of the potential Un.
 **/
void updateDerPot(vector<Oftsc>& W,        //parameterization in NC coordinates
                  vector<Ofsc>& Xe,        //Earth position
                  vector<Ofsc>& Xm,        //Moon position
                  vector<Ofsc>& Xs,        //Sun position
                  vector<Oftsc>& PrimPt,   //Intermediate steps: contr. of order < m to order m
                  double Ke,               //Earth factor
                  double Km,               //Moon factor
                  double Ks,               //Sun factor
                  vector<Oftsc>& Un,       //potential to update
                  int m);                  //order

/**
 *  \brief Update the vector field FWc at order m, knowing the parameterization W and the potential Un.
 **/
void applyVF(vector<Ofsc>& alpha,
             vector<Oftsc>& W,
             vector<Oftsc>& FWc,
             vector<Oftsc>& Un,
             int m);

//---------------------------------------------------------------------------------------------------------------------------------------
//         Mixed style
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Initialize the arrays that contains the segregated families of solutions in the mixed style.
 **/
void initMS(int **VIs,
            int  *VIn,
            int manType);
/**
 *  \brief Free the memory allocated by initMS.
 **/
void freeMS(int **VIs,
            int  *VIn,
            int ims);


/**
 * \brief Inner routine. Checks that the indix \c ind is in the array \c VI, of size \c sVI.
 *
 *  Note that there is a shif of indices: VI contains indices in the form 1, 2, 3, ..., whereas ind is an indix of the form 0, 1, 2, ...
 *
 **/
bool inVI(int ind, int *VI, int sVI);
/**
 * \brief Inner routine. Checks that all indices of the form \c kv[VI[i]-1] are null, where \c VI is an array of integers of size \c sVI.
 *
 *  Note that there is a shif of indices: VI contains indices in the form 1, 2, 3, ..., whereas ind is an indix of the form 0, 1, 2, ...
 *
 **/
bool nullIndices(int *kv, int *VI, int sVI);
/**
 * \brief Inner routine. For all arrays VIs[vi] of size VIn[vi], checks the boolean:
 *                                      inVI(p, VIs[vi], VIn[vi]) && nullIndices(kv, VIs[vi], VIn[vi])
 *        This the standard test for mixed style parameterization (see Haro 2014, section 2.2.3).
 *
 **/
bool isNF(int p, int *kv, int **VIs, int *VIn, int numberOfSets);

#endif // PMT_H_INCLUDED
