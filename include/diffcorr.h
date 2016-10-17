#ifndef DIFFCORR_H_INCLUDED
#define DIFFCORR_H_INCLUDED

/**
 * \file  diffcorr.h
 * \brief Contains all the routines that perform differential correction procedures.
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>

#include "qbtbp.h"
#include "qbcp.h"
#include "odezero.h"

extern "C"{
#include "gnuplot_i.h"
#include "nrutil.h"
}

//------------------------------------------------------------------------------------------------------------
//Differential correction
//------------------------------------------------------------------------------------------------------------
/**
 *  \brief Performs a differential correction procedure on ystart in order to get a periodic orbit of period t1.
 *         The algorithm assumes that the orbit is in the xy plane, is symmetric wrt to the x-axis, and has a period T = t1.
 *
 *  \param ystart the initial conditions that needs to be corrected.
 *  \param t1 the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param N:  the number of points in the plot if the steps of the diffcorr are plotted (see isPlotted).
 *  \param isPlotted:  if true, the steps of the diffcorr procedure are plotted in a temporary gnuplot window.
 *
 *   In brief: the algorithm corrects  [ ystart[0] ystart[4] ] in order to get [y[1] y[3]](t1)
 *             i.e. when the trajectory crosses the line t = t1.
 **/
int differential_correction(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);
/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);
/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps_mns(double ystart[], double nullvector[], double fvv[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted, int isFirst);

/**
 *  \brief Performs a differential correction procedure on ystart in order to get a periodic orbit of period t1.
 *
 *  \param ystart the initial conditions that needs to be corrected.
 *  \param t1 the final integration time (period of the orbit).
 *  \param eps_diff: the desired precision (magnitude of the final error).
 *  \param d:  a gsl_odeiv2_driver object used to integrate the equations of motion.
 *  \param N:  the number of points in the plot if the steps of the diffcorr are plotted (see isPlotted).
 *  \param isPlotted:  if true, the steps of the diffcorr procedure are plotted in a temporary gnuplot window.
 *
 **/
int differential_correction_T(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);

//------------------------------------------------------------------------------------------------------------
// ODE PLOT - NEW VERSION
//------------------------------------------------------------------------------------------------------------
/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in either NC or SYS coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1. Print in txt files is included vis \c isStored integer.
 **/
int odePlot2(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color, int isNormalized, int isStored, string filename);

//------------------------------------------------------------------------------------------------------------
// ODE PLOT
//------------------------------------------------------------------------------------------------------------
/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1.
 **/
int odePlot(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color);

/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1.
 **/
int odePlot3D(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color);

/**
 *  \brief Integrate the state y[] up to t = t1 on a Npoints grid, in NC coordinates, and plot the corresponding result in SYS coordinates (e.g. EM or SEM coord.).
 *         Then the results are plotted on a temporary gnuplot window via the handle *h1. Print in txt files is included
 **/
int odePlotprint(const double y[], int N, double t1, gsl_odeiv2_driver *d, gnuplot_ctrl  *h1, int Npoints, int color, string filename);

/**
 *  \brief Same as odePlot, but with more flexibility on the parameters: can choose the title, the line style, the line type and the line color.
 **/
int odePlotGen(const double y[],
               int N, double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc);

/**
 *  \brief Same as odePlotGen, but the 3D result are also printed in a txt file, with the name (folder+"orbits/"+title+".txt").
 **/
int odePrintGen(const double y[],
               int N,
               double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc,
               string folder);

//------------------------------------------------------------------------------------------------------------
//Differential correction: not used
//------------------------------------------------------------------------------------------------------------
/**
 * \brief Differential correction with a fixed final time used in the continutation procedures from a model to another (e.g CRTBP to BCP).
 *        The state is 7-dimensional. The last component is deps, where eps is the continuation parameter. DEPRECATED.
 **/
int differential_correction_deps_pac(double y0[], double nullvector[], double fvv[], double ds, double t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);
/*
  Differential correction - Minimum norm in all 6 dimensions
 */
int differential_correction_MN(double ystart[], double t1, double eps_diff, gsl_odeiv2_driver *d, int isPlotted);
/*
 Differential correction in the case of minimum norm (MN)
 FOR LYAPUNOV ORBITS ONLY
*/
int differential_correction_MN_planar(double ystart[], double yhalf[], double *tcross, double t1, OdeStruct *ode_s);
/*
 Differential correction in the case of x0 fixed
 FOR LYAPUNOV ORBITS ONLY
*/
int differential_correction_x0_fixed(double ystart[], double yhalf[], double *tcross, double t1, OdeStruct *ode_s);

#endif // DIFFCORR_H_INCLUDED
