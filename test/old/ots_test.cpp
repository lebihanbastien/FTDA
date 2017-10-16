/**
 * \file ots_test.cpp
 * \brief Test file for the Ots (Taylor series) class.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

#include "timec.h"
#include "ots_test.h"
#define VAR 2.0

/**
 * \fn void ots_test()
 * \brief Main routine for Ots class test.
 */
void ots_test()
{
    char ch;
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "         Test routine for the Oftsh class          " << endl;
    cout << "    specialized in the form Oftsh<Ofs <U> >        " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "The following public functions are tested:         " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "void Ots<T>::evaluate(U X[], T& z)" << endl;
    cout << "void Ots<T>::evaluate_conjugate(U X[], T& z)" << endl;

    cout << "See source code for details. " << endl;
    coutlp();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);

    //------------------------
    // Tests.
    //------------------------
    ots_test_pows();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;


    ots_test_divs();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ots_test_smult();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ots_test_sprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ots_test_smprod();
    cout << "Press Enter to proceed with the tests." << endl;
    scanf("%c",&ch);
    cout << endl;

    ots_test_conjugate();


    cout << "---------------------------------------------------" << endl;
    cout << "End of test session.                               " << endl;
    cout << "---------------------------------------------------" << endl;
}

/**
 * \fn void ots_test_conjugate()
 * \brief Test of the routine: Ots<T>::conjugate().
 */
void ots_test_conjugate()
{
    cout << "--------------------------------------" << endl;
    cout << " Test of the routines:                " << endl;
    cout << " Ots<T>::conjugate()                  " << endl;
    cout << "--------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2;
    double complex resd1, resd2;

    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX;

    //Set random coefs in xd and xdc
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();

    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //------------------------------------------------------------------------
    //Set results
    //------------------------------------------------------------------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate_conjugate(Xvard, resd1);
    res1  = res1;
    resd1 = conj(resd1);

    //------------------------------------------------------------------------
    //Take operation
    //------------------------------------------------------------------------
    xd.conjugate();
    xdc.conjugate();

    //------------------------------------------------------------------------
    //Set results
    //------------------------------------------------------------------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);

    //------------------------------------------------------------------------
    //Comparison. double case
    //------------------------------------------------------------------------
    cout << "1. double case.  " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;

    //------------------------------------------------------------------------
    //Comparison. double complex case
    //------------------------------------------------------------------------
    cout << "2. double complex case. " << endl;
    cout << "Expected result:        " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result:        " << creal(resd2)  << cimag(resd2) << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ots_test_smult()
 * \brief Test of the routine: Ots<T>& Ots<T>::smult(Ots<T> const& a, T const& m).
 */
void ots_test_smult()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>& Ots<T>::smult(Ots<T> const& a, T const& m)                      " << endl;
    cout << "-------------------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2;
    double complex resd1, resd2;

    //coeff
    double c = 2.0;
    double complex cd = 2.0+I;


    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    res1  = res1  + res2*c;
    resd1 = resd1 + resd2*cd;

    //Take operation
    //--------------
    xd.smult(xd2, c);
    xdc.smult(xdc2, cd);

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>& Ots<T>::mult(Ots<T> const& a, T const& m)                       " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;
    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    res1  = res2*c;
    resd1 = resd2*cd;

    //Take operation
    //--------------
    xd.mult(xd2, c);
    xdc.mult(xdc2, cd);

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ots_test_sprod()
 * \brief Test of the routine: Ots<T>::sprod(Ots<T> const& a, Ots<T> const& b).
 */
void ots_test_sprod()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>& Ots<T>::sprod(Ots<T> const& a, Ots<T> const& b)                 " << endl;
    cout << "-------------------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);
    Ots< double > xd3(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc3(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2, res3;
    double complex resd1, resd2, resd3;


    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    xd3.evaluate(Xvar, res3);
    xdc3.evaluate(Xvard, resd3);

    res1  = res1  + res2*res3;
    resd1 = resd1 + resd2*resd3;

    //Take product
    //--------------
    tic();
    xd.sprod(xd2, xd3);
    double trots = toc();
    cout << "Product in : " << trots << "s." << endl;

    tic();
    xdc.sprod(xdc2, xdc3);
    double tcots = toc();
    cout << "Product in : " << tcots << "s." << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;
    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    xd3.evaluate(Xvar, res3);
    xdc3.evaluate(Xvard, resd3);

    res1  = res1  + res2*res3;
    resd1 = resd1 + resd2*resd3;

    //Take product
    //--------------
    tic();
    for(int k = 0; k<= xd.getOrder(); k++) xd.sprod(xd2, xd3, k);
    cout << "Product in : " << toc() << "s." << endl;

    tic();
    for(int k = 0; k<= xd.getOrder(); k++) xdc.sprod(xdc2, xdc3, k);
    cout << "Product in : " << toc() << "s." << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;



    cout << "-----------------------" << endl;
    cout << " Comparison with OFTS  " << endl;
    cout << "-----------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Ofts< Ofs<double> > xd_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double> > xd2_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc2_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double> > xd3_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc3_fts(REDUCED_NV, OTS_ORDER, 1, 0);

    //Set random coefs
    xd_fts.setRandomCoefs();
    xdc_fts.setRandomCoefs();
    xd2_fts.setRandomCoefs();
    xdc2_fts.setRandomCoefs();
    xd3_fts.setRandomCoefs();
    xdc3_fts.setRandomCoefs();


    //Take real product
    //--------------
    tic();
    xd_fts.ofts_sprod(xd2_fts, xd3_fts);
    double trofts = toc();
    cout << "Real OFTS product in : " << trofts << "s." << endl;
    coutsp();
    cout << "which is " << trofts/trots << " times longer than with OTS objects." << endl;
    coutlp();

    //Take complex product
    //--------------
    tic();
    xdc_fts.ofts_sprod(xdc2_fts, xdc3_fts);
    double tcofts = toc();
    cout << "Complex OFTS product in : " << tcofts << "s." << endl;
    coutsp();
    cout << "which is " << tcofts/tcots << " times longer than with OTS objects." << endl;
    coutlp();
}

/**
 * \fn void ots_test_smprod()
 * \brief Test of the routine: Ots<T>& Ots<T>::smprod(Ots<T> const& a, Ots<T> const& b, T const& m).
 */
void ots_test_smprod()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>& Ots<T>::smprod(Ots<T> const& a, Ots<T> const& b, T const& m)    " << endl;
    cout << "-------------------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);
    Ots< double > xd3(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc3(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2, res3;
    double complex resd1, resd2, resd3;

    double c = 2.0;
    double complex c2 = 2.0+1.0*I;

    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    xd3.evaluate(Xvar, res3);
    xdc3.evaluate(Xvard, resd3);

    res1  = res1  + c*res2*res3;
    resd1 = resd1 + c2*resd2*resd3;

    //Take product
    //--------------
    tic();
    xd.smprod(xd2, xd3, c);
    double trots = toc();
    cout << "Product in : " << trots << "s." << endl;

    tic();
    xdc.smprod(xdc2, xdc3, c2);
    double tcots = toc();
    cout << "Product in : " << tcots << "s." << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;


    cout << "-----------------------" << endl;
    cout << " 2. At order n         " << endl;
    cout << "-----------------------" << endl;
    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    //Set results
    //--------------
    xd.evaluate(Xvar, res1);
    xdc.evaluate(Xvard, resd1);

    xd2.evaluate(Xvar, res2);
    xdc2.evaluate(Xvard, resd2);

    xd3.evaluate(Xvar, res3);
    xdc3.evaluate(Xvard, resd3);

    res1  = res1  + c*res2*res3;
    resd1 = resd1 + c2*resd2*resd3;

    //Take product
    //--------------
    tic();
    for(int k = 0; k<= xd.getOrder(); k++) xd.smprod(xd2, xd3, c, k);
    cout << "Product in : " << toc() << "s." << endl;

    tic();
    for(int k = 0; k<= xd.getOrder(); k++) xdc.smprod(xdc2, xdc3, c2, k);
    cout << "Product in : " << toc() << "s." << endl;

    //Set results
    //--------------
    xd.evaluate(Xvar, res2);
    xdc.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;



    cout << "-----------------------" << endl;
    cout << " Comparison with OFTS  " << endl;
    cout << "-----------------------" << endl;
    //Initialization of OFTS objects, used to initialize the OFTSH objects
    //--------------
    Ofts< Ofs<double> > xd_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double> > xd2_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc2_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double> > xd3_fts(REDUCED_NV, OTS_ORDER, 1, 0);
    Ofts< Ofs<double complex> > xdc3_fts(REDUCED_NV, OTS_ORDER, 1, 0);

    //Set random coefs
    xd_fts.setRandomCoefs();
    xdc_fts.setRandomCoefs();
    xd2_fts.setRandomCoefs();
    xdc2_fts.setRandomCoefs();
    xd3_fts.setRandomCoefs();
    xdc3_fts.setRandomCoefs();

    //Coeffs to mult
    Ofs<double> ofs(0), temp(0);
    Ofs<double complex> ofsd(0), tempd(0);
    ofs.setCoef(c, 0);
    ofsd.setCoef(c2, 0);

    //Take real product
    //--------------
    tic();
    xd_fts.ofts_smprod_t(xd2_fts, xd3_fts, ofs, temp);
    double trofts = toc();
    cout << "Real OFTS product in : " << trofts << "s." << endl;
    coutsp();
    cout << "which is " << trofts/trots << " times longer than with OTS objects." << endl;
    coutlp();

    //Take complex product
    //--------------
    tic();
    xdc_fts.ofts_smprod_t(xdc2_fts, xdc3_fts, ofsd, tempd);
    double tcofts = toc();
    cout << "Complex OFTS product in : " << tcofts << "s." << endl;
    coutsp();
    cout << "which is " << tcofts/tcots << " times longer than with OTS objects." << endl;
    coutlp();
}

/**
 * \fn void ots_test_divs()
 * \brief Test of the routine: Ots<T>& Ots<T>::divs(Ots<T> const& a, Ots<T> const& b).
 */
void ots_test_divs()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>& Ots<T>::divs(Ots<T> const& a, Ots<T> const& b)                  " << endl;
    cout << "-------------------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);
    Ots< double > xd3(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc3(REDUCED_NV, OTS_ORDER);
    Ots< double > xd4(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc4(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2;
    double complex resd1, resd2;

    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    //Take product
    //--------------
    tic();
    xd.divs(xd2, xd3);  //xd = xd2/xd3
    xd4.sprod(xd, xd3); //xd4 = xd*xd3 = xd2
    double trots = toc();
    cout << "Division in : " << trots << "s." << endl;

    tic();
    xdc.divs(xdc2, xdc3);  //xdc = xdc2/xdc3
    xdc4.sprod(xdc, xdc3); //xdc4 = xdc*xdc3 = xdc2
    double tcots = toc();
    cout << "Division in : " << tcots << "s." << endl;

    //Set results
    //--------------
    xd2.evaluate(Xvar, res1);
    xdc2.evaluate(Xvard, resd1);

    xd4.evaluate(Xvar, res2);
    xdc4.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}

/**
 * \fn void ots_test_pows()
 * \brief Test of the routine: Ots<T>::pows(Ots<T> const& a,  T const& alpha).
 */
void ots_test_pows()
{
    cout << "-------------------------------------------------------------------------" << endl;
    cout << " Test of the routines:                                                   " << endl;
    cout << " Ots<T>::pows(Ots<T> const& a,  T const& alpha)                          " << endl;
    cout << "-------------------------------------------------------------------------" << endl;

    //------------------------------------------------------------------------
    //Initialization of OTS objects, used to initialize the OTSH objects
    //------------------------------------------------------------------------
    Ots< double > xd(REDUCED_NV, OTS_ORDER);
    Ots< double complex> xdc(REDUCED_NV, OTS_ORDER);
    Ots< double > xd2(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc2(REDUCED_NV, OTS_ORDER);
    Ots< double > xd3(REDUCED_NV, OTS_ORDER);
    Ots< double complex > xdc3(REDUCED_NV, OTS_ORDER);

    //------------------------------------------------------------------------
    //result storage
    //------------------------------------------------------------------------
    double res1, res2;
    double complex resd1, resd2;

    //------------------------------------------------------------------------
    //Variables
    //------------------------------------------------------------------------
    double Xvar[REDUCED_NV];
    double complex Xvard[REDUCED_NV];

    //array of coordinates (randomly init)
    for(int i = 0; i < REDUCED_NV; i++) Xvar[i]  = 1.0*rand()/RAND_MAX;
    for(int i = 0; i < REDUCED_NV; i++) Xvard[i] = 1.0*rand()/RAND_MAX + I*rand()/RAND_MAX ;


    cout << "-----------------------" << endl;
    cout << " 1. At all terms       " << endl;
    cout << "-----------------------" << endl;

    //Set random coefs
    xd.setRandomCoefs();
    xdc.setRandomCoefs();
    xd2.setRandomCoefs();
    xdc2.setRandomCoefs();
    xd3.setRandomCoefs();
    xdc3.setRandomCoefs();

    double alpha = -3.5;

    //Take product
    //--------------
    tic();
    xd.pows(xd2, alpha);        //xd = xd2^alpha
    xd3.pows(xd, 1.0/alpha);    //xd3 = xd2
    double trots = toc();
    cout << "Division in : " << trots << "s." << endl;

    tic();
    xdc.pows(xdc2, alpha);      //xdc = xdc2^alpha
    xdc3.pows(xdc, 1.0/alpha);  //xdc3 = xdc2
    double tcots = toc();
    cout << "Division in : " << tcots << "s." << endl;

    //Set results
    //--------------
    xd2.evaluate(Xvar, res1);
    xdc2.evaluate(Xvard, resd1);

    xd3.evaluate(Xvar, res2);
    xdc3.evaluate(Xvard, resd2);


    //Comparison. double case
    //--------------
    cout << "1. double case.         " << endl;
    cout << "Expected result: " << creal(res1) << cimag(res1) << "i" << endl;
    cout << "Obtained result: " << creal(res2) << cimag(res2) << "i" << endl;
    res2 -= res1;
    cout << "Delta: " << cabs(res2) << endl;


    //Comparison. double complex case
    //--------------
    cout << "2. double complex case.         " << endl;
    cout << "Expected result: " << creal(resd1)  << cimag(resd1) << "i" << endl;
    cout << "Obtained result: " << creal(resd2)  << cimag(resd2)  << "i" << endl;
    resd2 -= resd1;
    cout << "Delta: " << cabs(resd2) << endl;
}
