#ifndef FTDA_H_INCLUDED
#define FTDA_H_INCLUDED

/**
 * \file Manip.h
 * \brief Basic operations for manipulation of polynomials
 *       (exponents in reverse lexicographic order, number of coefficient in
 *       homogeneous polynomials...).
 *       Note that the algebraic operations on series (sum, product) are not done here
 *       but in ofts.h/ofts.tpp.
 * \author BLB using code by Angel Jorba.
 *
 *      Based on Jorba 1999 (http://www.maia.ub.es/~angel/soft.html).
 */

#include <iostream>
#include <iomanip>
#include <limits.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>


/**
 *  \brief Basic operations for manipulation of polynomials.
 *         Manip::psi(i,j) for i=2..nov, j=1...nor, is the number of monomials of
 *         degree j with i variables.
 *         Note that the algebraic operations on series (sum, product) are not done here
 *         but in ofts.h/ofts.tpp.
 **/
class Manip;
class Manip
{
private:
    ///static maximum degree
    static int nor;
    ///static number of variables
    static int nov;
    ///contains psi(i,j) for i=2..nov, which is the number of monomials of degree j with i variables.
    static long  int **psi;

public:

    /**
     *   \brief Initializes the table psi(i,j), which contains the number of
     *          homogeneous polynomials of order j with i variables.
     *
     *   parameters:
     *   nr: maximum degree we are going to work with. it can not be greater than 63.
     *   nv: number of variables.
     *
     *   returned value: number of kbytes allocated by the internal tables.
     *   Based on a routine by Angel Jorba, 1999.
     **/
    static int init(int nv, int nr);
    /**
     *  \brief Frees the space allocated by Manip::init.
     **/
    static void free();
    /**
     *   \brief Returns the number of monomials of degree nr with nv variables,
     *          making use of the table Manip::psi.
     *
     *  parameters:
     *  nv: number of variables
     *  nr: order we are interested in (input).
     *  returned value: number of monomials of order no.
     **/
    static long int nmon(int nv, int nr);
    /**
     *  \brief given a multiindex k, this routine computes the next one
     *  according to the (reverse) lexicographic order.
     *
     *  parameters:
     *  k: array of nv components containing the multiindex. It is overwritten on exit
     *   (input and output).
     **/
    static void prxkt(int k[], int nv);
};



//----------------------------------------------------------------------------------------
//Binomial coefficients. These routines do not make use of Manip::psi.
//----------------------------------------------------------------------------------------
/**
 * \brief Computes the binomial coefficient (x y).
 **/
unsigned long binomial(unsigned long n, unsigned long k);

/**
 * \brief Computes the binomial coefficient (n k).
 **/
unsigned long gcd_ui(unsigned long x, unsigned long y);

#endif // FTDA_H_INCLUDED
