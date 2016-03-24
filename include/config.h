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
extern QBCP SEM;               ///Sun-Earth-Moon system
extern QBCP_L SEML;            ///Sun-Earth-Moon system around Li
extern QBCP_L SEML_SEM;        ///Sun-Earth-Moon system around Li

/**
 *   \brief Initialization of the environnement (Sun, Earth, Moon, Li...).
 *
 *    The global variables SEM and SEML are initialized in order to describe the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem around the Lagrangian point LI,
 *    LI being given as an integer in the file parameters.h. The default coordinates system is the normalized-centered (NC) coordinates of the Earth-Moon system.
 *    Note that, in order for the initialization to be complete - Sun-(Earth+Moon) motion given as Fourier series within SEML -
 *    the routine qbtbp(true) must have been run *once.
 **/
void init_env(int li_EM, int li_SEM, int isNormalized, int model, int fwrk, int pmStyle, int manType_EM, int manType_SEM);

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
 *   \brief Initialization of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 **/
void initCM(QBCP_L &qbcp);
/**
 *   \brief Update of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
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
 *  Main routine for the initialization of the COC
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
