#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>

#include "tests.h"
#include "ofts.h"
#include "ofs.h"

using namespace std;

//Tests of OTS routines
void oneVariableTest(int order, int corder)
{
    //Fourier-Taylor serie
    Ofts< Ots<double> > ofts1(1, order,   2, corder);
    Ofts< Ots<double> > ofts2(1, order,   2, corder);
    Ofts< Ots<double> > ofts3(1, order,   2, corder);

    //Initialisation
    //Note that we only initialize non-null value for k<=order/2 and |l|<= nf/2
    //So that the final results contains ALL non-null coefficients.
    int nf = corder/2;
    Ofsd ofs1(nf);

    for(int ord= 0; ord<= order/2; ord++)
    {
        for(int l=-nf/2; l<=nf/2; l++) ofs1.setCoef((double)rand()/(10.0*(l*l+1)*RAND_MAX),l);

        //Fourier to Taylor
        fs2ts(ofts1.getCoef(ord,0), ofs1);
    }

    double evalOp1, evalOp2;
    double t = 0.2;
    double epsilon = 1.0/389;


    //Product
    //--------------------------------------------------------------------------
    //for(int k =0; k<=order; k++) ofts3.sprod(ofts1, ofts1, k);
    //evalOp1 = creal(fts2scalar(ofts1, epsilon, t)*fts2scalar(ofts1, epsilon, t));
    //evalOp2 = creal(fts2scalar(ofts3, epsilon, t));

    //Smult
    //--------------------------------------------------------------------------
    Ots<double> ots1(2, corder);
    fs2ts(&ots1, ofs1);
    double alpha = 2.534;
    for(int k =0; k<=order; k++) ofts3.ofts_smult_tu(ofts1, ots1, alpha, k);
    evalOp1 = creal(fts2scalar(ofts1, epsilon, t)*ofs1.evaluate(t)*alpha);
    evalOp2 = creal(fts2scalar(ofts3, epsilon, t));


    //Power: for power function, one needs a lot of orders to get a good approximation.
    //--------------------------------------------------------------------------
    //Pure taylor series seem to work fine. However, Fourier-Taylor series lead to divergence.
    //Does this test have a real meaning? Maybe better to test in situ i.e. with real series
    //For which the Fourier coefficients have a meaning and with decreasing coefficients when |l| increases.
    //--------------------------------------------------------------------------
    //double alpha = 1.5;
    //for(int k =0; k<=order; k++) ofts3.pows(ofts1, alpha, k);
    //evalOp1 = creal(cpow(fts2scalar(ofts1, epsilon, t), alpha)); //creal(cpow(fts2scalar(ofts1, epsilon, t), alpha));
    //evalOp2 = creal(fts2scalar(ofts3, epsilon, t));


    cout << "F(eval(ofs1)): " << evalOp1 << endl;
    cout << "eval(F(ofs1)): " << evalOp2 << endl;
}

//Tests of OFS routines
void oneVariableTest_Ofs(int order, int corder)
{
    //Fourier-Taylor serie
    Ofts< Ofsd > ofts1(1, order,   2, corder);
    Ofts< Ofsd > ofts2(1, order,   2, corder);
    Ofts< Ofsd > ofts3(1, order,   2, corder);

    //Initialisation
    //Note that we only initialize non-null value for k<=order/2 and |l|<= nf/2
    //So that the final results contains ALL non-null coefficients.
    int nf = corder;
    Ofsd ofs1(nf);

    for(int ord= 0; ord<= order/2; ord++)
    {
        for(int l=-nf/2; l<=nf/2; l++) ofs1.setCoef((double)rand()/(10.0*(l*l+1)*RAND_MAX),l);

        //Fourier to Taylor
        ofts1.getCoef(ord,0)->ccopy(ofs1);
    }

    double evalOp1, evalOp2;
    double t = 0.2;
    double epsilon = 1.0/389;

    //Product
    //--------------------------------------------------------------------------
    //for(int k =0; k<=order; k++) ofts3.sprod(ofts1, ofts1, k);
    //evalOp1 = creal(fts2scalar(ofts1, epsilon, t)*fts2scalar(ofts1, epsilon, t));
    //evalOp2 = creal(fts2scalar(ofts3, epsilon, t));

    //Smult
    //--------------------------------------------------------------------------
    //double alpha = 2.534;
    //for(int k =0; k<=order; k++) ofts3.smult(ofts1, ofs1, alpha, k);
    //evalOp1 = creal(fts2scalar(ofts1, epsilon, t)*ofs1.evaluate(t)*alpha);
    //evalOp2 = creal(fts2scalar(ofts3, epsilon, t));


    //Power: for power function, one needs a lot of orders to get a good approximation.
    //--------------------------------------------------------------------------
    //Pure taylor series seem to work fine. However, Fourier-Taylor series lead to divergence.
    //Does this test have a real meaning? Maybe better to test in situ i.e. with real series
    //For which the Fourier coefficients have a meaning and with decreasing coefficients when |l| increases.
    //--------------------------------------------------------------------------
    //Guarantee that ofts1[0] = 1.0!!
    ofts1.getCoef(0,0)->zero();
    ofts1.setCoef0(1.0,0,0);
    cout << ofts1 << endl;

    double alpha = -3.0/2;
    for(int k =0; k<=order; k++) ofts3.ofts_pows(ofts1, alpha, k);
    evalOp1 = creal(cpow(fts2scalar(ofts1, epsilon, t), alpha+0.0*I)); //creal(cpow(fts2scalar(ofts1, epsilon, t), alpha));
    evalOp2 = creal(fts2scalar(ofts3, epsilon, t));


    cout << "F(eval(ofs1)): " << evalOp1 << endl;
    cout << "eval(F(ofs1)): " << evalOp2 << endl;
}


//Tests of OTS derivatives routines
void derivTest()
{
    cout << "-------------------------------" << endl;
    cout << "Testing the partial derivatives" << endl;
    cout << "       of OTS objects.         " << endl;
    cout << "-------------------------------" << endl;
    cout << "Test on series of 3 variables and of order 3. " << endl;
    Ots<cdouble> ots(3, 3);
    Ots<cdouble> ots2(3, 3);
    ots.setRandomCoefs();

    cout << "Initial serie:" << endl;
    cout << ots << endl;
    cout << "------------------" << endl;

    ots2.der(ots, 1);
    cout << "Derivative wrt to x1:" << endl;
    cout << ots2 << endl;
    cout << "------------------" << endl;

    ots2.der(ots, 2);
    cout << "Derivative wrt to x2:" << endl;
    cout << ots2 << endl;
    cout << "------------------" << endl;

    ots2.der(ots, 3);
    cout << "Derivative wrt to x3:" << endl;
    cout << ots2 << endl;
    cout << "------------------" << endl;

}
