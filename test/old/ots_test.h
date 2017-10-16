#ifndef OTS_TEST_H_INCLUDED
#define OTS_TEST_H_INCLUDED

//std
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <sstream>

//custom
#include "oftsh.h"
#include "ofts.h"

using namespace std;

/**
 * \fn void ofts_test()
 * \brief Main routine for Ofts class test.
 */
void ots_test();

/**
 * \fn void ots_test_conjugate()
 * \brief Test of the routine: Ots<T>::conjugate().
 */
void ots_test_conjugate();

/**
 * \fn void ots_test_sprod()
 * \brief Test of the routine: Ots<T>::sprod(Ots<T> const& a, Ots<T> const& b).
 */
void ots_test_sprod();

/**
 * \fn void ots_test_smprod()
 * \brief Test of the routine: Ots<T>& Ots<T>::smprod(Ots<T> const& a, Ots<T> const& b, T const& m).
 */
void ots_test_smprod();

/**
 * \fn void ots_test_smult()
 * \brief Test of the routine: Ots<T>& Ots<T>::smult(Ots<T> const& a, T const& m).
 */
void ots_test_smult();

/**
 * \fn void ots_test_divs()
 * \brief Test of the routine: Ots<T>& Ots<T>::divs(Ots<T> const& a, Ots<T> const& b).
 */
void ots_test_divs();

/**
 * \fn void ots_test_pows()
 * \brief Test of the routine: Ots<T>::pows(Ots<T> const& a,  T const& alpha).
 */
void ots_test_pows();

#endif // OTS_TEST_H_INCLUDED
