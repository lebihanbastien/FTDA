#ifndef POINCARE_H_INCLUDED
#define POINCARE_H_INCLUDED

extern "C"
{
    #include "nrutil.h"
    #include "multimin.h"
}


#include "pmode.h"
#include "init.h"
#include "parameters.h"

//Integration methods
#define DUAL_INT 1
#define DUAL_INT_NO_RESET 2
#define DUAL_INT_STEPPED  3
#define SINGLE_INT 4

//Types of map
#define PMAP  1
#define TMAP  2
#define EMAP  3
#define IMAP  4
#define HMAP  5

/**
 * \file poincare.h
 * \brief Computation of Poincare, Error, and Period map for the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 */

/**
 *  \struct Pmap
 *  \brief  Poincare map parameters.
 *
 *  Contains all necessary parameters and constants to define a suitable Poincare map, Stroboscopic map, or Precision map
 *  (energy level, order of expansions, maximum number of events, initial time...).
 *
 **/
typedef struct Pmap Pmap;
struct Pmap
{
    //Expansions
    int order;           //the order of the expansions
    int ofs_order;       //the order of the expansions

    //Integration & event
    double tf;           //final time
    double t0;           //init time
    double tt;           //test time
    int projFreq;        //test projection frequency
    double T;            //period
    int max_events;      //maximum z=0 events allowed
    double threshold;    //the threshold above which a reset of z(t) is needed
    double pabs;         //absolute tolerance during integration
    double prel;         //relative tolerance during integration
    double proot;        //the precision on the roots z = 0

    //Energy
    double Hv;           //value of the Hamiltonian
    double H0;           //value of the Hamiltonian @Li
    double dHv;          //Hv-H0

    //Grid
    int  gsize;          //size of the grid
    double gmin;         //min of the grid
    double gmax;         //max of the grid

    //Checking divergence
    double maxRad;

    //Dimension
    int vdim;

    //Model bool
    bool isGS;

    //Type
    int type;
};

/**
 *  \struct Orbit
 *  \brief  Defines a given orbit with proper arrays to store results.
 **/
typedef struct Orbit Orbit;
struct Orbit
{
    //------------------------------------------------------------------------------------
    //Parent
    //------------------------------------------------------------------------------------
    Pmap     *pmap;         //Poincare map (parent)
    QBCP_L   *qbcp_l;       //QBCP around a given Li point (parent)

    //------------------------------------------------------------------------------------
    //Parameterization (common to all orbits)
    //------------------------------------------------------------------------------------
    vector<Oftsc>*  W;      //z(t) = W(s(t), t)
    vector<Oftsc>*  Wh;     //zh(t) = Wh(s(t), t)
    matrix<Oftsc>*  DW;     //Jacobian of W
    Ofsc* ofs;              //Auxiliary Ofs object
    double  n;              //Pulsation of the QBCP
    //------------------------------------------------------------------------------------
    //COC (common to all orbits)
    //------------------------------------------------------------------------------------
    matrix<Ofsc>* PC;       //COC matrix
    vector<Ofsc>* V;        //COC vector

    //------------------------------------------------------------------------------------
    //For event detection
    //------------------------------------------------------------------------------------
    value_params *val_par;  //Event parameters
    value_function *fvalue; //fvalue for event detection
    double** z0_mat;        //Pointer towards the stored position of events (NC)
    double** zh0_mat;       //Pointer towards the stored position of events (TFR)
    double** s0_mat;        //Pointer towards the stored position of events (RCM)
    double*  te_mat;        //Pointer towards the stored time of events
    double*  hz;            //Energy z(t) at each event
    double*  hw;            //Energy W(s(t),t) at each event
    double*  ePm;           //Projection error at each event
    int last_indix;         //indix of the last computed event in z0_mat
    int reset_number;       //number of reset during computation (for dual integration)
    int      *nevent;       //the label of the events

    //------------------------------------------------------------------------------------
    //Characteristics
    //------------------------------------------------------------------------------------
    double   *z0;           //Initial position in NC coordinates dim = 6
    double   *zh0;          //Initial position in TFR coordinates dim = 6
    double   *si;           //Initial RCM configuration dim = REDUCED_NV
    double   *s0d;          //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
    cdouble  *s0;           //Initial position in CCM8 coordinates (real+imag part) dim = 4
    double   *xf;           //Final position dim = 6
    double    tf;           //final time after computation
    double    eOm;          //Orbit error
    int       int_method;   //integration method used to compute the orbit; -1 if not computed
    int       label;        //label of the orbit
    int       vdim;         //the dimension along wich we guarantee a certain initial energy value


    //------------------------------------------------------------------------------------
    //ODE integration
    //------------------------------------------------------------------------------------
    OdeStruct *ode_s_6;     //NC ode struct
    OdeStruct *ode_s_8;     //TFC ode struct

    OdeStruct *ode_s_6_root; //NC ode struct
    OdeStruct *ode_s_8_root; //TFC ode struct
};


//----------------------------------------------------------------------------------------
//
//  Poincare map
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Computes a Poincare map
 *   \param pmap a reference to the Poincare map parameters
 *   \param isPlot if true, the Poincare map is plotted during the computation
 *
 *    Requires initCM and initCOC
 **/
void pmap_build(Pmap &pmap, int append, int return_method,  bool isPlot, bool isPar);

/**
 *   \brief Computes a Stroboscopic map
 *   \param pmap a reference to the Stroboscopic map parameters
 *   \param isPlot if true, the Stroboscopic map is plotted during the computation
 *
 *    Requires initCM and initCOC
 **/
void tmap_build(Pmap &pmap, int append, int method, bool isPlot, bool isPar);

/**
 *   \brief Precision on a Poincare map. Parallelized version
 *   \param pmap a reference to the Poincare maps parameters
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_precision(Pmap &pmap, int append, bool isPar);

/**
 *   \brief Invariance error computation
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_invariance_error(Pmap &pmap, int append, bool isPar, double hzmax);

/**
 *   \brief Invariance error computation
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_invariance_error_random(Pmap &pmap, int append, bool isPar, double hzmax);

/**
 *   \brief Test error computation
 *   \param pmap a reference to the Poincare maps parameters
 *   \param hzmax the maximum energy value allowed.
 *          Note that only positive dhz are selected
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_test_error(Pmap &pmap, int append, bool isPar, double hzmax);

/**
 *   \brief Energy of the initial conditions on a Poincare map. Parallelized version
 *   \param pmap a reference to the Poincare maps parameters
 *
 *    Requires initCM and initCOC
 *
 *   REMARK: may be good to "force" p36(t0= 0.0) (e.g. pmap.t0 =  +1.044814582930593 for L2) so that each IC begins on z = 0 plane
 *   If so, the way H(0) = cst is guaranteed must be changed because we need also to ensure that s4 = 0.0 (which is not the case, since it is the variable
 *   that ensures H(0) = cst
 **/
void pmap_energy(Pmap &pmap, int append, bool isPar, double hzmax);


//----------------------------------------------------------------------------------------
//
//  Init of one orbit
//
//----------------------------------------------------------------------------------------

/**
    \brief Initialize an orbit wrt a Poincare map so that H(orbit.s0) = H(Pmap)
 **/
int orbit_init_pmap(Orbit &orbit, double st0[]);

/**
    \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(Orbit &orbit, const double si[], double t0);

/**
 *   \brief Computes the hamiltonian at the position st0
 *   \param pmap a reference to the pmap that carries a set of useful parameters
 *   \param st0 the input state
 *
 *   WARNING: Direct use of CM inside this
 **/
double orbit_ham(Pmap &pmap, double st0[]);

/**
    \brief Computes the difference between an given Ham value and the current configuration: (orbit.s1, orbit.s2, sr, 0.0) (RTBB) or (orbit.s1, orbit.s2, orbit.s3, sr)
    (QBCP)
    \param  sr a double to complete the current tested configuration
    \param  params a pointer to the orbit with a given Ham value (why void? the idea is to a have a generic function, but might be useless at this point)
    \return the difference between the two hamiltonians
 **/
double deltham(double sr, void *params);


//----------------------------------------------------------------------------------------
//
//  Computation of one orbit
//
//----------------------------------------------------------------------------------------
/**
    \brief Computes the poincare map of one given orbit
    \param orbit a pointer to the orbit
**/
int orbit_compute_pmap(Orbit *orbit, int method);

/**
 *   \brief Computes the poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *          The precise roots of z = 0 are computed at the end of the integration.
 *   \param orbit a pointer to the orbit
 **/
int dual_pmap(Orbit *orbit);

/**
 *   \brief Computes the poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *          The precise roots of z = 0 are computed during the integration.
 *   \param orbit a pointer to the orbit
 **/
int dual_pmap_stepped(Orbit *orbit);

/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *   \param orbit a pointer to the orbit
 **/
int dual_tmap(Orbit *orbit);

/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a single method: only the NC vector field is integrated.
 *          The current state is projected on the center manifold after each half-period.
 *   \param orbit a pointer to the orbit
 **/
int single_tmap(Orbit *orbit);

/**
 *   \brief Computes the period/stroboscopic/poincare map of one given orbit, with a single method: only the NC vector field is integrated.
 *          The current state is projected on the center manifold after each half-period.
 *   \param orbit a pointer to the orbit
 **/
int single_tmap_tested(Orbit *orbit);

/**
 *   \brief Computes the error/precision map of one given orbit, with a dual method: both vector fields (NC and reduced) are integrated in parallel.
 *   \param orbit a pointer to the orbit
 **/
int dual_emap(Orbit *orbit);

/**
 *   \brief Refine the root z = 0 for a dual integration NC+rvf. Used in dual_pmap_stepped.
 **/
int refine_root_dual(Orbit *orbit,
                     double *yv,
                     double *si,
                     cdouble *s1,
                     double *t,
                     double *yv_mat,
                     double *sv_mat,
                     double *t_mat,
                     int events);

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap(Orbit *orbit);
int single_pmap_proj(Orbit *orbit);

/**
 *   \brief Computes the poincare map of one given orbit, with a single method: only the Nc vector field is computed. The distance with respect to the
 *          central manifold is evaluated each time that z = 0.
 *   \param orbit a pointer to the orbit
 **/
int single_pmap_plot(Orbit *orbit,gnuplot_ctrl  *h1, int color);

/**
 *   \brief Refine the root z = 0 for a single integration NC+rvf. Used in single_pmap.
 **/
int refine_root(Orbit *orbit,
                double *yv,
                double *t,
                double *s1,
                double *yv_mat,
                double *t_mat,
                double omega1,
                double omega3,
                int events);



//----------------------------------------------------------------------------------------
//
//  Print & Read
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Print the poincare map of and orbit in a txt file
 **/
void orbit_pmap_fprint(Orbit *orbit, string filename, int append);

/**
 *   \brief Print the poincare map of and orbit in a txt file. Lighter version
 **/
void orbit_pmap_fprint_small(Orbit *orbit, string filename, int append);

/**
    \brief Print the precision map of and orbit in a txt file
**/
void orbit_precision_fprint(Orbit *orbit, string filename, int append);

/**
 *   \brief Writing the precision map header in the txt file
 **/
void header_precision_fprint(string filename);

/**
 *   \brief Writing the poincare map header in the txt file
 **/
void header_pmap_fprint(string filename);

/**
 *   \brief Writing the poincare map header in the txt file. Lighter version
 **/
void header_pmap_fprint_small(string filename);

/**
    \brief Print the energy map of and orbit in a txt file
**/
void orbit_energy_fprint(Orbit *orbit, string filename, double hzmax, int append);

/**
 *   \brief Writing the energy map header in the txt file
 **/
void header_energy_fprint(string filename);

/**
    \brief Print the energy map of and orbit in a bin file
**/
void orbit_imap_fprint_bin(Orbit *orbit, string filename, int append);


/**
    \brief Print the energy map of and orbit in a bin file
**/
void orbit_energy_fprint_bin(Orbit *orbit, string filename, int append);

//----------------------------------------------------------------------------------------
//
// Steppers
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Integrates the current state yv/sv one step, using both NC and rvf vector fields.
 **/
double gslc_dual_step( Orbit *orbit,
                       double yv[],
                       double sv[],
                       double eO[],
                       double zr1[],
                       double *t,
                       double *tr,
                       double t1,
                       int resetOn);


/**
 *   \brief Integrates the current state yv/sv up to t = t1, using both NC and rvf vector fields.
 **/
int gslc_dual_evolve(Orbit *orbit,
                     double yv[],
                     double sv[],
                     double eO[],
                     double z1[],
                     double *t,
                     double t1,
                     double threshold);


//----------------------------------------------------------------------------------------
//
// Plotting (deprecated)
//
//----------------------------------------------------------------------------------------
//Plot one orbit
void orbit_plot(Orbit *orbit, gnuplot_ctrl *h1, int type, int points, OdeStruct *ode_s_6, OdeStruct *ode_s_8);
//Plot one orbit (3D)
void orbit_plot_3d(Orbit *orbit, gnuplot_ctrl *h1, int points, OdeStruct *ode_s_6, OdeStruct *ode_s_8);
//Plot poincare map for one orbit
void orbit_poincare_plot(Orbit *orbit, gnuplot_ctrl *h1, gnuplot_ctrl *h2, int color);
//Plot T map for one orbit
void orbit_Tmap_plot(Orbit *orbit, gnuplot_ctrl *h1, gnuplot_ctrl *h2, int color);


//----------------------------------------------------------------------------------------
//
//  Square distance & minimization of it (deprecated)
//
//----------------------------------------------------------------------------------------
/**
    \brief Given an initial guess st0, computes the min argument of square distance between a given configuraton z1 and z = W(g(st), t) (see sqdist).
    i.e. st1 = argmin sqdist(st, z1, t, &orbit)

    \param Orbit a reference to the current orbit
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param st0 the initial guess in real TFC configuration (dim REDUCED_NV)
    \param st1 the min argument to update real TFC configuration (dim REDUCED_NV)
    \param multimin_method an integer to select a multidimensionnal minimization method
**/
void argmin_sqdist(Orbit &orbit, double z1[], double t, double st0[], double st1[], int multimin_method);


/**
    \brief Computes the square distance between a given configuraton z1 and z0 = W(g(st0), t)

    \param st0 an array of REDUCED_NV double which gives the configuration to study in real TFC coordinates
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param params a set of parameters needed to perform evaluations of W and g
    \return a double, the square distance between a given configuraton z1 and z0 = W(g(st0), t)
**/
double sqdist(double st0[], double z1[], double t, void *params);


/**
    \brief Computes the gradient square distance between a given configuraton z1 and z0 = W(g(st0), t) (see sqdist)

    \param st0 an array of REDUCED_NV double which gives the configuration to study in real TFC coordinates
    \param df the gradient to update
    \param z1 the target NC configuration (dim 6)
    \param t the current time
    \param params a set of parameters needed to perform evaluations of W and g
**/
void dsqdist(double st0[], double df[], double z1[], double t, void *params);


//----------------------------------------------------------------------------------------
//
//  Orbit C structure handling
//
//----------------------------------------------------------------------------------------
/**
    \brief Initialize one orbit structure
 **/
void init_orbit(Orbit *orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Oftsc>* DW,
                matrix<Ofsc>*  PC,
                vector<Ofsc>*  V,
                Pmap *pmap,
                QBCP_L *qbcp_l,
                Ofsc* orbit_ofs,
                int vdim,
                int label);
/**
    \brief Initialize one orbit structure
 **/
void init_orbit(Orbit *orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Oftsc>* DW,
                matrix<Ofsc>*  PC,
                vector<Ofsc>*  V,
                value_params   *val_par,
                value_function *fvalue,
                OdeStruct *ode_s_6,
                OdeStruct *ode_s_8,
                OdeStruct *ode_s_6_root,
                OdeStruct *ode_s_8_root,
                Pmap *pmap,
                QBCP_L *qbcp_l,
                Ofsc* orbit_ofs,
                int vdim,
                int label);

/**
    \brief Free one orbit
 **/
void free_orbit(Orbit *orbit);


//----------------------------------------------------------------------------------------
//
//  Inner routines
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Initialize the grid on which the poincare map will be evaluated
 **/
void init_grid(double *grid, double gmin, double gmax, int gsize);

/**
 *   \brief Initialize the ode structure for the NC vector field
 **/
void init_ode_NC(OdeStruct &ode_s_6, Pmap &pmap);

/**
 *   \brief Initialize the ode structure for the reduced vector field
 **/
void init_ode_CCM(OdeStruct &ode_s_8, RVF &rvf, Ofsc &rvf_ofs, Pmap &pmap);


#endif // POINCARE_H_INCLUDED
