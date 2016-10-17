#ifndef ERTBP_H_INCLUDED
#define ERTBP_H_INCLUDED

#include <math.h>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
//Custom
#include "ofts.h"
#include "define_env.h"
//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_sf_bessel.h>

/**
 * \file ertbp.h
 * \brief Routines of the Elliptic Restricted Three-Body Problem (ERTBP) (header)
 * \author BLB.
 * \date November 2015
 * \version 1.0
 */

 /**
 *  \brief Main routine to compute the Elliptic Three-Body Problem in Ofs format.
 */
void ertbp(int li_EM, int li_SEM, int coordsys);


#endif // ERTBP_H_INCLUDED
