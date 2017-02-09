#ifndef PMODE_H_INCLUDED
#define PMODE_H_INCLUDED

/**
 * \file pmode.h
 * \brief Integration of the equations of motion for the outputs of the parameterization method.
 * \author BLB.
 * \date 2016
 * \version 1.0
 */


#include "pmt.h"
#include "pmcoc.h"


// Structure for the reduced vector field
typedef struct RVF RVF;
struct RVF
{
    vector<Oftsc>* fh;   //Reduced vector field
    Ofsc* ofs;           //Auxiliary Ofs object
    int order;              //Order of the evaluation     (<= OFTS_ORDER)
    int ofs_order;          //Order of the ofs evaluation (<= OFS_ORDER)
    double n;               //Pulsation of the QBCP
};


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Integration of the pm
//
//---------------------------------------------------------------------------------------------------------------------------------------
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
int qbfbp_fh(double t, const double y[], double f[], void *params_void);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Tests
//
//---------------------------------------------------------------------------------------------------------------------------------------
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
void pmErrorvsOrderTest(int nkm, int km[], double si[]);

/**
 *  \brief Test of the pm of the central manifold of L1,2 on a given orbit, through the computation of the orbital error eO, along an orbit initialized by the pm.
 *         Various orders are tested, among the available values.
 *
 *   Requires initCM().
 *
 **/
void pmEOvsOrderTest(int nkm, int km[], double si[]);

/**
 *  \brief Test of the pm of the central manifold of L1,2 on a given orbit, through the computation of the orbital error eO, along an orbit initialized by the pm.
 *         Various OFS orders (of the Fourier expansions) are tested, among the available values.
 *   \param order the order of the Taylor expansions (OFTS objects).
 *   Requires initCM().
 *
 **/
void pmOfsOrderTest(int order);

/**
 *  \brief Evaluates the contributions of each order in W to the computation of W(s,t), with an arbitrary state (s,t)
 **/
void pmContributions();

/**
 *  \brief Evaluates the l1/linf norm of each order in W.
 **/
void pmNorms();

/**
 *  \brief Small divisors under a certain value
 **/
void pmSmallDivisors(double sdmax);

/**
 *  \brief Evaluates the variations of the initial conditions (IC) wrt to the order in W(s,t), with an arbitrary state (s,t)
 **/
void pmTestIC();



//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Plot in tests
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computations of various errors (eO, eI, eH) along an orbit initialized by the pm of the center manifold W(s,t)
 **/
int errorPlot(const double st0[],      //RCM initial conditions
              vector<Oftsc>& Wh,    //TFC manifold
              matrix<Ofsc>& PC,     //COC matrix: z = PC*zh+V
              vector<Ofsc>& V,      //COC vector: z = PC*zh+V
              vector<Oftsc>& FW,    //NC vector field
              vector<Oftsc>& DWf,   //Jacobian of FW
              vector<Oftsc>& Wdot,  //Partial derivative of W wrt time
              double t1,               //Final integration tim
              gsl_odeiv2_driver *dnc,  //driver for NC integration
              gsl_odeiv2_driver *drvf, //drive for RVF integration (reduced vector field)
              int Npoints,             //Number of points on which the errors are estimated
              QBCP_L& qbcp_l,         //current QBCP
              int order,               //Order for the eval of the OFTS objects
              int ofs_order,           //Order for the eval of the OFS objects
              gnuplot_ctrl  **ht,      //Gnuplot handlers
              int color);              //Color of the plots

/**
 *  \brief Computations of various the orbial error eO along an orbit initialized by the pm of the center manifold W(s,t)
 **/
int eOPlot(const double st0[],         //RCM initial conditions
              vector<Oftsc>& Wh,    //TFC manifold
              matrix<Ofsc>& PC,     //COC matrix: z = PC*zh+V
              vector<Ofsc>& V,      //COC vector: z = PC*zh+V
              double t1,               //Final integration tim
              gsl_odeiv2_driver *dnc,  //driver for NC integration
              gsl_odeiv2_driver *drvf, //drive for RVF integration (reduced vector field)
              int Npoints,             //Number of points on which the errors are estimated
              QBCP_L& qbcp_l,         //current QBCP
              int order,               //Order for the eval of the OFTS objects
              int ofs_order,           //Order for the eval of the OFS objects
              gnuplot_ctrl  **ht,      //Gnuplot handlers
              int color);              //Color of the plots



#endif // PMODE_H_INCLUDED
