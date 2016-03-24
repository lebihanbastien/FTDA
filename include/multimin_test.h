#ifndef MULTIMIN_TEST_H_INCLUDED
#define MULTIMIN_TEST_H_INCLUDED

/**
 * \file multimin_test.h
 * \brief Test of the function minimization routines defined in multimin.h
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "timec.h"

extern "C"
{
    #include "nrutil.h"
    #include "multimin.h"
}

using namespace std;

//Test function of the routine frprmn of multimin.c
void test_frprmn();

#endif // MULTIMIN_TEST_H_INCLUDED
