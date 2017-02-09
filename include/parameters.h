#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/**
 * \file parameters.h
 * \brief Configuration file. User-defined constants of the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

#include <string>
#include <complex.h>

using namespace std;

//------------------------------------------------------------------------------------
// Environment
//------------------------------------------------------------------------------------
/// Numerotation of the Solar System planets and objects, consistent with JPL's HORIZON numerotation. Here, the Sun.
#define SUN 10
/// Mercury
#define MERCURY 199
/// Venus
#define VENUS 299
/// Earth
#define EARTH 399
/// Mars
#define MARS 499
/// Jupiter
#define JUPITER 599
/// Saturn
#define SATURN 699
/// Uranus
#define URANUS 799
/// Neptune
#define NEPTUNE 899
/// Pluto
#define PLUTO 999
/// The Moon
#define MOON 301
/// Custom indix for the Sun+Earth system
#define EARTH_AND_MOON 700

//------------------------------------------------------------------------------------
//   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
//------------------------------------------------------------------------------------
//Model
#define M_RTBP  0 ///< RTBP model indix
#define M_QBCP  1 ///< QBCP model indix
#define M_BCP   2 ///< BCP model indix
#define M_ERTBP 3 ///< ERTBP model indix
//Frameworks
#define F_EM 0  ///< Earth-Moon framework indix
#define F_SEM 1 ///< Sun-(Earth+Moon) framework indix
#define F_SE 2  ///< Sun-Earth framework indix

//Misc
#define POTENTIAL_ORDER 60 ///< Maximum order of the potential of the primaries
#define NV 6               ///< Number of state variables
#define OFS_NV 1           ///< Number of variables in the OFS object (a priori always 1)

//------------------------------------------------------------------------------------
//   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
//------------------------------------------------------------------------------------
#define XY 1  ///< XY plane, for plotting purposes
#define YZ 2  ///< YZ plane, for plotting purposes
#define XZ 3  ///< XZ plane, for plotting purposes

#define MAX_EVENTS 100 ///< Maximum events in Poincare & Period maps

//------------------------------------------------------------------------------------
//   Global constants
//------------------------------------------------------------------------------------
extern int OFTS_ORDER;  ///< Order of the Taylor series in Fourier-Taylor series
extern int OFS_ORDER;   ///< Order of the Fourier series in Fourier-Taylor series
extern int OTS_ORDER;   ///< Order of the Taylor series in pure Taylor series
extern int MODEL_TYPE;  ///< Type of model (chosen in the constants beginning by "M_")
extern int REDUCED_NV;  ///< Number of reduced variables (e.g. 4 for a center manifold, 5 for a center-stable...)

//------------------------------------------------------------------------------------
//   Misc
//------------------------------------------------------------------------------------
extern const std::string SSEPR; ///< Small separation mark
extern const std::string SEPR;  ///< Big separation mark

//------------------------------------------------------------------------------------
//   typedef
//------------------------------------------------------------------------------------
typedef complex double cdouble;

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp();

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp();

/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp();

//------------------------------------------------------------------------------------
//   Print
//------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double *y, int n);

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble *y, int n);

//------------------------------------------------------------------------------------
//   Norm
//------------------------------------------------------------------------------------
/**
 *  Euclidian norm computed on the first k components of a complex vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(cdouble z0[], int k);

/**
 *  Euclidian norm computed on the first k components of a double vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(double z0[], int k);


#endif // PARAMETERS_H_INCLUDED
