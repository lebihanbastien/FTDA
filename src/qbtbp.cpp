#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>

#include "qbtbp.h"

using namespace std;

//n angular velocity
//ms sun mass
//as orbit radius of the sun
void qbtbp (Ofts< Ots<double> > &ofts_z, Ofts< Ots<double> > &ofts_Z, double n, double ms, double as, double ns)
{
    int order = ofts_z.getOrder();
    int nv = ofts_z.getVariables();
    int corder = ofts_z.getCOrder();
    int cnv = ofts_z.getCVariables();
    //Here nf = corder/2 because nf is the order of the Ofs coefficients (Fourier form), not the Ots (Taylor form)
    int nf = ofts_z.getCOrder()/2;

    double mu = 0.0121505816;
    double eps = 1/as;


    //---------------------------------------------------------
    //Order 0
    //---------------------------------------------------------
    ofts_z.setCoef(1.0,0);
    ofts_Z.setCoef(1.0,0);


    //---------------------------------------------------------
    //Order 1
    //One may choose to set it at the beginning or to compute it in the loop.
    //---------------------------------------------------------
    //Fourier series
    Ofs<double> u1(nf);
    Ofs<double> v1(nf);

    u1.setCoef(-ms/(as*as)/6.0, 0);
    u1.setCoef(-ms/(as*as)*(24.0*n*n+24*n+9.0)/(64*pow(n,4.0)-16*n*n), -2);
    u1.setCoef(+ms/(as*as)*9.0/(64*pow(n,4.0)-16*n*n), 2);

    //Fourier to Taylor: this is the solution @order 1.
    //-----------------------------------
    //fs2ts(ofts_z.getCoef(1,0), u1);
    //fs2ts(ofts_Z.getCoef(1,0), v1);

    //---------------------------------------------------------
    //Order n
    //---------------------------------------------------------
    Ofts< Ots<double> > z1(nv, order, cnv, corder);
    Ofts< Ots<double> > z2(nv, order, cnv, corder);
    Ofts< Ots<double> > z3(nv, order, cnv, corder);
    Ofts< Ots<double> > z4(nv, order, cnv, corder);


    Ofts< Ots<double> > a1(nv, order, cnv, corder);
    Ofts< Ots<double> > a2(nv, order, cnv, corder);
    Ofts< Ots<double> > a3(nv, order, cnv, corder);
    Ofts< Ots<double> > a4(nv, order, cnv, corder);
    Ofts< Ots<double> > a5(nv, order, cnv, corder);
    Ofts< Ots<double> > a6(nv, order, cnv, corder);
    Ofts< Ots<double> > a7(nv, order, cnv, corder);


    Ofts< Ots<double> > b1(nv, order, cnv, corder);
    Ofts< Ots<double> > b2(nv, order, cnv, corder);
    Ofts< Ots<double> > b3(nv, order, cnv, corder);
    Ofts< Ots<double> > b4(nv, order, cnv, corder);
    Ofts< Ots<double> > b5(nv, order, cnv, corder);
    Ofts< Ots<double> > b6(nv, order, cnv, corder);
    Ofts< Ots<double> > b7(nv, order, cnv, corder);

    Ofts< Ots<double> > Pm(nv, order, cnv, corder);
    Ofts< Ots<double> > Qm(nv, order, cnv, corder);

    //sigma1 = exp(-itheta)
    Ofs<double> sigma1(nf);
    sigma1.setCoef(1.0, -1);
    Ots<double> s1(cnv, corder);
    fs2ts(&s1, sigma1);


    //sigma2 = exp(+itheta)
    Ofs<double> sigma2(nf);
    sigma2.setCoef(1.0, 1);
    Ots<double> s2(cnv, corder);
    fs2ts(&s2, sigma2);

    //epsilon = 0+epsilon+0 (order 1 of Taylor serie)
    Ofts< Ots<double> > epsilon(nv, order, cnv, corder);
    epsilon.setCoef(1.0,1);

    //Order 0
    int m;
    for(m = 0; m<=0 ; m++)
    {
        z1.ccopy(ofts_z, m);
        //z1 = \bar{z1}
        z1.conjugate(m);
        //z2 = \bar{z}^^(-3/2)
        z2.pows(z1, -3.0/2, m);
        //z3 = z^(-1/2)
        z3.pows(ofts_z, -1.0/2, m);
        //z4 = z2*z3;
        z4.sprod(z2, z3, m);


        //a1 = mu*exp(itheta)*z
        a1.smult(ofts_z, s2, mu, m);
        //a2 = epsilon*a1
        a2.sprod(a1, epsilon, m);
        //a3+= Z
        a3.smult(ofts_Z, 1.0, m);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.smult(a2, -1.0, m);
        //a4 = a3^(-1/2)
        a4.pows(a3, -1.0/2, m);
        //a5 = \bar{a3}
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2)
        a6.pows(a5, -3.0/2, m);
        //a7 = a4*a6;
        a7.sprod(a4,a6, m);


        //b1 = (1-mu)*exp(ithetb)*z
        b1.smult(ofts_z, s2, (1.0-mu), m);
        //b2 = epsilon*b1
        b2.sprod(b1, epsilon, m);
        //b3+= Z
        b3.smult(ofts_Z, 1.0, m);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.smult(b2, 1.0, m);
        //b4 = b3^(-1/2)
        b4.pows(b3, -1.0/2, m);
        //b5 = \bbr{b3}
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2)
        b6.pows(b5, -3.0/2, m);
        //b7 = b4*b6;
        b7.sprod(b4,b6, m);


        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.smult(b7, s1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.smult(a7, s1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.smult(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.smult(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.smult(a7, -ns*ns*(1-mu), m);

    }

    Ofs<double>  Pfm(nf);
    Ofs<double>  Qfm(nf);

    Ofs<double> ufm(nf);
    Ofs<double> vfm(nf);

    //-------------------------
    //Recurrence
    //-------------------------
    for(m = 1; m <= ofts_z.getOrder(); m++)
    {

        //z1 =ofts_z at order m-1
        z1.ccopy(ofts_z, m-1);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m-1);
        //z2 = \bar{z}^^(-3/2)
        if(m>1) z2.pows(z1, -3.0/2, m-1);
        z2.pows(z1, -3.0/2, m);
        //z3 = z^(-1/2)
        if(m>1) z3.pows(ofts_z, -1.0/2, m-1);
        z3.pows(ofts_z, -1.0/2, m);
        //z4 = z2*z3;
        z4.sprod(z2, z3, m);

        //cout << "Order: " << m << endl;
        //cout << "z1: \n" << z1 << endl;
        //cout << "z2: \n" << z2 << endl;



        //a1 = mu*exp(itheta)*z at order m-1
        if(m>1) a1.smult(ofts_z, s2, mu, m-1);
        //a2 = epsilon*a1 at order m
        a2.sprod(a1, epsilon, m);
        //a3+= Z
        if(m>1) a3.smult(ofts_Z, 1.0, m-1);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.smult(a2, -1.0, m);
        //a4 = a3^(-1/2)
        if(m>1) a4.pows(a3, -1.0/2, m-1);
        a4.pows(a3, -1.0/2, m);
        //a5 = \bar{a3}
        if(m>1)
        {
            a5.ccopy(a3, m-1);
            a5.conjugate(m-1);
        }
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2)
        if(m>1) a6.pows(a5, -3.0/2, m-1);
        a6.pows(a5, -3.0/2, m);
        //a7 = a4*a6;
        a7.sprod(a4,a6, m);

//        cout << "\n---------------------------\n" << endl;
//        cout << "Order: " << m << endl;
//        cout << "a3: \n" << a3 << endl;
//        cout << "a4: \n" << a4 << endl;

        //b1 = (1-mu)*exp(ithetb)*z
        if(m>1) b1.smult(ofts_z, s2, (1.0-mu), m-1);
        //b2 = epsilon*b1
        b2.sprod(b1, epsilon, m);
        //b3+= Z
        if(m>1) b3.smult(ofts_Z, 1.0, m-1);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.smult(b2, 1.0, m);
        //b4 = b3^(-1/2)
        if(m>1) b4.pows(b3, -1.0/2, m-1);
        b4.pows(b3, -1.0/2, m);
        //b5 = \bbr{b3}
        if(m>1)
        {
            b5.ccopy(b3, m-1);
            b5.conjugate(m-1);
        }
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2)
        if(m>1) b6.pows(b5, -3.0/2, m-1);
         b6.pows(b5, -3.0/2, m);
        //b7 = b4*b6;
        b7.sprod(b4,b6, m);

        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.smult(b7, s1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.smult(a7, s1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.smult(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.smult(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.smult(a7, -ns*ns*(1-mu), m);

        //-------------------------
        // Solving equations
        //-------------------------
        Pfm.ts2fs(Pm.getTerm(m)->getCoef(0));
        Qfm.ts2fs(Qm.getTerm(m)->getCoef(0));

        //Order 0
        ufm.setCoef(-1.0/3*Pfm.getCoef(0), 0);  //u0 = -1/3*p0
        vfm.setCoef(-1.0/(3*ns*ns)*Qfm.getCoef(0), 0); //v0 = -1/(3*ns^2)*p0

        //Order 1 to nf
        //solving a 2*2 system
        double k1, k2, k3, l1, l2, l3, uj, vj;
        for(int j = 1; j<= nf; j++)
        {
            k1 = -pow(j*n,2.0) + 2*j*n - 3.0/2;
            k2 = -pow(j*n,2.0) - 2*j*n - 3.0/2;
            k3 = -3.0/2;

            l1 = -pow(j*n,2.0) + 2*j*n*ns - 3.0/2*ns*ns;
            l2 = -pow(j*n,2.0) - 2*j*n*ns - 3.0/2*ns*ns;
            l3 = -3.0/2*ns*ns;


            //u-j
            uj = ( k2*Pfm.getCoef(-j) - k3*Pfm.getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, -j);
            //uj
            uj = (-k3*Pfm.getCoef(-j) + k1*Pfm.getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, j);

            //v-j
            vj = ( l2*Qfm.getCoef(-j) - l3*Qfm.getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, -j);
            //vj
            vj = (-l3*Qfm.getCoef(-j) + l1*Qfm.getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, j);

        }


        //Update z and Z
        fs2ts(ofts_z.getCoef(m,0), ufm);
        fs2ts(ofts_Z.getCoef(m,0), vfm);

        cout << "Order " << m << " completed." << endl;
    }

    //bj and cj
    Ofs<double> bj(nf);
    fts2fs(&bj, ofts_z, eps);
    Ofs<double> cj(nf);
    fts2fs(&cj, ofts_Z, eps);


    //cout << ofts_z << endl;

    ofstream bjstream("data/bj.txt");
    bjstream << "bj: \n" << bj << endl;
    bjstream.close();
    ofstream cjstream("data/cj.txt");
    cjstream << "cj: \n" << cj << endl;
    cjstream.close();

}


