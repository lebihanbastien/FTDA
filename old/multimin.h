#ifndef MULTIMIN_H_INCLUDED
#define MULTIMIN_H_INCLUDED

/**
 * \file multimin.h
 * \brief Function minimization in \f$ R^n \f$.
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

#include "nrutil.h"
#include "math.h"

//Constants to switch between various multmin methods
#define MULTIMIN_SIMPLEX  1  //Simplex method
#define MULTIMIN_FRPRMN   2  //Fletcher-Reeves-Polak-Ribiere minimization
#define MULTIMIN_DFRPRMN  3  //Fletcher-Reeves-Polak-Ribiere minimization with gradient

//----------------------------------------------------------------------------------------
//
//  Simplex method
//
//----------------------------------------------------------------------------------------
/*
Multidimensional minimization of the function funk(x) where x[0..ndim-1] is a vector in ndim
dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[0..ndim]
[0..ndim-1] is input. Its ndim+1 rows are ndim -dimensional vectors which are the vertices of
the starting simplex. Also input is the vector y[0..ndim] , whose components must be pre-
initialized to the values of funk evaluated at the ndim+1 vertices (rows) of p ; and ftol the
fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
nfunk gives the number of function evaluations taken.
*/
void amoeba(double z1[], double t, void *params,
            double **p, double y[], int ndim, double ftol,
            double (*funk)(double [], double[], double, void *),
            int *nfunk);

/*
Extrapolates by a factor fac through the face of the simplex across from the high point, tries
it, and replaces the high point if the new point is better.
*/
double amotry(double z1[], double t, void *params, double **p, double y[], double psum[], int ndim,
             double (*funk)(double [], double[], double, void *), int ihi, double fac);

//----------------------------------------------------------------------------------------
//
//  Fletcher-Reeves-Polak-Ribiere minimization
//
//----------------------------------------------------------------------------------------
/*
Given a starting point p[0..n-1] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations.
*/
void frprmn(double z1[], double t, void *params, double p[], int n, double ftol, int *iter, double *fret,
            double (*funk)(double [], double[], double, void *), void (*dfunk)(double [], double[], double[], double, void *));

/*
Given a starting point p[0..n-1] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine dlinmin is called to perform line minimizations.
*/
void dfrprmn(double z1[], double t, void *params, double p[], int n, double ftol, int *iter, double *fret,
            double (*funk)(double [], double[], double, void *), void (*dfunk)(double [], double[], double[], double, void *));

/*
Given an n -dimensional point p[0..n-1] and an n -dimensional direction xi[0..n-1] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and brent .
*/
void linmin(double z1[], double t, void *params, double p[], double xi[], int n, double *fret, double (*funk)(double [], double[], double, void *));


/*
Given an n -dimensional point p[0..n-1] and an n -dimensional direction xi[0..n-1] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and dbrent .
*/
void dlinmin(double z1[], double t, void *params, double p[], double xi[], int n, double *fret, double (*funk)(double [], double[], double, void *),
             void (*dfunk)(double [], double[], double[], double, void *));

/*
  Artificial function of one variable along one direction of minimization in the routine linmin
*/
double f1dim(double z1[], double t, void *params, double x);

/*
  Gradient of an artificial function of one variable along one direction of minimization in the routine linmin
*/
double df1dim(double z1[], double t, void *params, double x);


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  1D min
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/*
Brent minimization

Given a function f , and given a bracketing triplet of abscissas ax , bx , cx (such that bx is
between ax and cx , and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value.
*/
double brent(double z1[], double t, void *params, double ax, double bx, double cx, double (*f)(double[], double, void *, double), double tol, double *xmin);

/*
Brent minimization with gradient

Given a function f and its derivative function df , and given a bracketing triplet of abscissas ax ,
bx , cx [such that bx is between ax and cx , and f(bx) is less than both f(ax) and f(cx) ],
this routine isolates the minimum to a fractional precision of about tol using a modification of
Brent’s method that uses derivatives. The abscissa of the minimum is returned as xmin , and
the minimum function value is returned as dbrent , the returned function value.
*/
double dbrent(double z1[], double t, void *params, double ax, double bx, double cx, double (*f)(double[], double, void *, double),
              double (*df)(double[], double, void *, double), double tol, double *xmin);
/*
Given a function func , and given distinct initial points ax and bx , this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax , bx , cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa , fb , and fc .
*/
void mnbrak(double z1[], double t, void *params, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double[], double, void *, double));


#endif // MULTIMIN_H_INCLUDED
