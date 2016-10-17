#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

/**
 * \file config.h
 * \brief Configuration file. It allows to initialize the environment of the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

#include <vector>
#include "parameters.h"
#include "ofts.h"
#include "matrix.h"
#include "define_env.h"

//---------------------------------------
//QBCP
//---------------------------------------
extern QBCP SEM;              //global structure that describes the Sun-Earth-Moon system
extern QBCP_L SEML;           //global structure that describes the Sun-Earth-Moon system around Li (initialized on the EM  framework)
extern QBCP_L SEML_SEM;       //global structure that describes the Sun-Earth-Moon system around Li (initialized on the SEM framework)

/**
 *   \brief Initialization of the environnement (Sun, Earth, Moon, Li...).
 *
 *    The global structures SEM, SEML and SEML_SEML are initialized in order to describe the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem about a given libration point. More precisely:
 *      - SEM is a QBCP structure that contains the numerical constants of the three primaries + the numerical constants of the associated CRTBPs.
 *
 *      - SEML is a QBCP_L structure that contains the numerical constants and parameters of the model (defined by the integer \c model) in the four possible coordinates systems: Earth-Moon (EM), Earth-Moon Normalized-Centered (EMNC), Sun-(Earth+Moon) (SEM), and Sun-(Earth+Moon) Normalized-Centered (SEMN). The default coordinate system is the EMNC about the libration point \c li_EM.
 *
 *      - SEML_SEM is a QBCP_L structure that contains the numerical constants and parameters of the model (defined by the integer \c model) in the four possible coordinates systems. The default coordinate system is the SEMNC about the libration point \c li_SEM.
 *
 *    The structures are also initialized in terms of manifolds. The integers \c manType_EM and \c manType_SEM define the type of manifold that is attached to the model (stable, center-stable...), while the integer \c pmStyle defines the style of parameterization that will be used (graph, normal form...).
 *
 *    Note that the routine qbtbp(true) must have been run at least *once*.
 *
 *    Note that, although everything is pretty much called "QBCP", the structures can harbour other models, such as the coupled CRTBP and the BCP (probably better to use more generic names in the future).
 *
 *    Note that the parameter \c isNormalized is probably deprecated, but is kept for compatibility with old code. Setting isNormalized = true = 1 is always preferable.
 **/
void init_env(int li_EM, int li_SEM, int isNormalized, int model, int pmStyle, int manType_EM, int manType_SEM);

//---------------------------------------
//Center manifold
//---------------------------------------
extern vector<Oftsc>  CM;    ///center manifold in NC coordinates
extern vector<Oftsc>  CMdot; ///center manifold in NC coordinates
extern vector<Oftsc> CMh;    ///center manifold in TFC coordinates
extern matrix<Oftsc> JCM;    ///jacobian of CM
extern vector<Oftsc>  Fh;    ///reduced vector field
extern vector<Oftsc>  DWFh;   ///JCM * Fh

/**
 *   \brief Initialization of the parameterization center manifold around a given libration point, encoded in the \c qbcp structure.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    These data files must have been previously computed through the use of the routine pm(int, int) or pmt(int, int, int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void initCM(QBCP_L &qbcp);
/**
 *   \brief Update of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    These data files must have been previously computed through the use of the routine pm(int, int) or pmt(int, int, int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables updated by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void updateCM(QBCP_L &qbcp);

//---------------------------------------
//COC
//---------------------------------------
extern matrix<Ofsc>  Mcoc;    ///COC matrix
extern matrix<Ofsc>  Pcoc;    ///COC matrix (Mcoc = Pcoc*Complex matrix)
extern vector<Ofsc>  Vcoc;    ///COC vector
extern matrix<Ofsc>  MIcoc;   ///COC matrix = inv(Mcoc)
extern matrix<Ofsc>  PIcoc;   ///COC matrix = inv(Pcoc)

/**
 *  Main routine for the update of the Change of Coordinates (COC), from TFC to NC coordinates.
 **/
void initCOC(QBCP_L &qbcp);

/**
 *  Main routine for the update of the COC
 **/
void updateCOC(QBCP_L &qbcp);

/**
 *  \brief Update a complex Fourier series given as the component (k,p) of a matrix matrix, from a given txt file.
 *  \param xFFT: the \c ofs<cdouble> to update
 *  \param matrix: the beginning of the name of the source txt file. e.g. "alpha"
 *  \param k the lign indix of the desired component in the matrix \c matrix.
 *  \param p the column indix of the desired component in the matrix \c matrix.
 *
 *  As an example, the call readCOC(xFFT, "alpha", 2, 1), will update xFFT with the file "alpha21.txt".
 **/
void readCOC(Ofsc& xFFT, string matrix, int k, int p);

/**
 *  \brief Lighter version of the init routine, with only P, PC and V retrieved.
 **/
void initCOC(matrix<Ofsc> &P,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &V,
             QBCP_L& qbcp_l);


/**
 *   \brief Number to string inner routine
 **/
string numTostring(double num);

#endif // CONFIG_H_INCLUDED
