#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <vector>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

//Custom
#include "custom_ode.h"
#include "qbtbp.h"

using namespace std;

//n angular velocity
//ms sun mass
//as orbit radius of the sun
void qbtbp_ofs (Ofs<double complex> &bjc, Ofs<double complex> &cjc, Ofts< Ofs<double> > &ofts_z, Ofts< Ofs<double> > &ofts_Z, double n, double ms, double as, double ns)
{
    int order = ofts_z.getOrder();
    int nv = ofts_z.getVariables();
    int corder = ofts_z.getCOrder();
    int cnv = ofts_z.getCVariables();
    int nf = ofts_z.getCOrder();

    double mu = 0.0121505816;
    double eps = 1/as;

    //---------------------------------------------------------
    //Order 0
    //---------------------------------------------------------
    ofts_z.setCoef(1.0, 0);
    ofts_Z.setCoef(1.0, 0);


    //---------------------------------------------------------
    //Order 1
    //---------------------------------------------------------
    //Fourier series
    Ofs<double> u1(nf);
    Ofs<double> v1(nf);

    u1.setCoef(-ms/(as*as)/6.0, 0);
    u1.setCoef(-ms/(as*as)*(24.0*n*n+24*n+9.0)/(64*pow(n,4.0)-16*n*n), -2);
    u1.setCoef(+ms/(as*as)*9.0/(64*pow(n,4.0)-16*n*n), 2);

    //Fourier to Taylor
    //fs2ts(ofts_z.getCoef(1,0), u1);
    //fs2ts(ofts_Z.getCoef(1,0), v1);

    //---------------------------------------------------------
    //Order n
    //---------------------------------------------------------
    Ofts< Ofs<double> > z1(nv, order, cnv, corder);
    Ofts< Ofs<double> > z2(nv, order, cnv, corder);
    Ofts< Ofs<double> > z3(nv, order, cnv, corder);
    Ofts< Ofs<double> > z4(nv, order, cnv, corder);
    Ofts< Ofs<double> > zspare(nv, order, cnv, corder);

    Ofts< Ofs<double> > a1(nv, order, cnv, corder);
    Ofts< Ofs<double> > a2(nv, order, cnv, corder);
    Ofts< Ofs<double> > a3(nv, order, cnv, corder);
    Ofts< Ofs<double> > a4(nv, order, cnv, corder);
    Ofts< Ofs<double> > a5(nv, order, cnv, corder);
    Ofts< Ofs<double> > a6(nv, order, cnv, corder);
    Ofts< Ofs<double> > a7(nv, order, cnv, corder);


    Ofts< Ofs<double> > b1(nv, order, cnv, corder);
    Ofts< Ofs<double> > b2(nv, order, cnv, corder);
    Ofts< Ofs<double> > b3(nv, order, cnv, corder);
    Ofts< Ofs<double> > b4(nv, order, cnv, corder);
    Ofts< Ofs<double> > b5(nv, order, cnv, corder);
    Ofts< Ofs<double> > b6(nv, order, cnv, corder);
    Ofts< Ofs<double> > b7(nv, order, cnv, corder);

    Ofts< Ofs<double> > Pm(nv, order, cnv, corder);
    Ofts< Ofs<double> > Qm(nv, order, cnv, corder);

    //sigma1 = exp(-itheta)
    Ofs<double> sigma1(nf);
    sigma1.setCoef(1.0, -1);
    //sigma2 = exp(+itheta)
    Ofs<double> sigma2(nf);
    sigma2.setCoef(1.0, 1);

    Ofts< Ofs<double> > epsilon(nv, order, cnv, corder);
    epsilon.setCoef(1.0, 1);

    int m;
    for(m = 0; m<=0 ; m++)
    {
        z1.ccopy(ofts_z, m);
        //z1 = \bar{z1}
        z1.conjugate(m);
        //z2 = \bar{z}^^(-3/2). At order 0, z1 = z2
        z2.ccopy(z1, m);
        //z3 = z^(-1/2). At order 0, z3 = z
        z3.ccopy(ofts_z, m);
        //z4 = z2*z3;
        z4.sprod(z2, z3, m);
        //z2 = \bar{z}^^(-1/2). At order 0, z1 = z2
        zspare.ccopy(z1, m);

        //a1 = mu*exp(itheta)*z
        a1.smult(ofts_z, sigma2, mu, m);
        //a2 = epsilon*a1
        a2.sprod(a1, epsilon, m);
        //a3+= Z
        a3.smult(ofts_Z, 1.0, m);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.smult(a2, -1.0, m);
        //a4 = a3^(-1/2). At order 0, a4 = a3
        a4.ccopy(a3, m);

        //a5 = \bar{a3}
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2). At order 0, a6 = a5
        a6.ccopy(a5, m);
        //a7 = a4*a6;
        a7.sprod(a4,a6, m);


        //b1 = (1-mu)*exp(ithetb)*z
        b1.smult(ofts_z, sigma2, (1.0-mu), m);
        //b2 = epsilon*b1
        b2.sprod(b1, epsilon, m);
        //b3+= Z
        b3.smult(ofts_Z, 1.0, m);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.smult(b2, 1.0, m);
        //b4 = b3^(-1/2). At order 0, b4 = b3
        b4.ccopy(b3, m);
        //b5 = \bbr{b3}
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2). At order 0, b6 = b5
        b6.ccopy(b5, m);
        //b7 = b4*b6;
        b7.sprod(b4,b6, m);


        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.smult(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.smult(a7, sigma1,  ms/(as*as),m);
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
        //z1 = ofts_z at order m-1
        z1.ccopy(ofts_z, m-1);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m-1);

        //z2 = \bar{z}^^(-3/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) z2.pows(z1, -3.0/2, m-1);
        z2.pows(z1, -3.0/2, m);

        //zspare = \bar{z}^^(-1/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) zspare.pows(z1, -1.0/2, m-1);
        zspare.pows(z1, -1.0/2, m);

        //z3 = z^(-1/2)
        //CAN BE USED BECAUSE z[0] = 1.0 !!
        if(m>1) z3.pows(ofts_z, -1.0/2, m-1);
        z3.pows(ofts_z, -1.0/2, m);

        //z4 = z2*z3;
        z4.sprod(z2, z3, m);

        //a1 = mu*exp(itheta)*z at order m-1
        if(m>1) a1.smult(ofts_z, sigma2, mu, m-1);
        //a2 = epsilon*a1 at order m
        a2.sprod(a1, epsilon, m);
        //a3+= Z at order m-1
        if(m>1) a3.smult(ofts_Z, 1.0, m-1);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.smult(a2, -1.0, m);

        //a4 = a3^(-1/2)
        //CAN BE USED BECAUSE a3[0] = 1.0 !!
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
        //CAN BE USED BECAUSE a5[0] = 1.0 !!
        if(m>1) a6.pows(a5, -3.0/2, m-1);
        a6.pows(a5, -3.0/2, m);

        //a7 = a4*a6;
        a7.sprod(a4,a6, m);

        //b1 = (1-mu)*exp(itheta)*z
        if(m>1) b1.smult(ofts_z, sigma2, (1.0-mu), m-1);
        //b2 = epsilon*b1
        b2.sprod(b1, epsilon, m);
        //b3+= Z
        if(m>1) b3.smult(ofts_Z, 1.0, m-1);
        //b3+= (1-mu)*exp(itheta)*epsilon*z
        b3.smult(b2, 1.0, m);

        //b4 = b3^(-1/2)
        //CAN BE USED BECAUSE b3[0] = 1.0 !!
        if(m>1) b4.pows(b3, -1.0/2, m-1);
        b4.pows(b3, -1.0/2, m);

        //b5 = \bar{b3}
        if(m>1)
        {
            b5.ccopy(b3, m-1);
            b5.conjugate(m-1);
        }
        b5.ccopy(b3, m);
        b5.conjugate(m);

        //b6 = b5^(-3/2)
        //CAN BE USED BECAUSE b5[0] = 1.0 !!
        if(m>1) b6.pows(b5, -3.0/2, m-1);
        b6.pows(b5, -3.0/2, m);

        //b7 = b4*b6;
        b7.sprod(b4,b6, m);

        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.smult(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.smult(a7, sigma1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.smult(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.smult(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.smult(a7, -ns*ns*(1-mu), m);

        //-------------------------
        // Solving equations
        //-------------------------
        Pfm.ccopy(Pm.getTerm(m)->getCoef(0));
        Qfm.ccopy(Qm.getTerm(m)->getCoef(0));


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
        ofts_z.getCoef(m,0)->ccopy(ufm);
        ofts_Z.getCoef(m,0)->ccopy(vfm);



    }

    cout << "QBTBP: End of recursive computation." << endl;

    //----------------------------------------------------------------------
    // bj and cj
    //----------------------------------------------------------------------
    //bj and cj
    Ofs<double> bj(nf);
    Ofs<double> cj(nf);

    fts2fs(&bj, ofts_z, eps);
    fts2fs(&cj, ofts_Z, eps);

    ofstream curentStream;
    curentStream.open("data/bj_ofs.txt");
    curentStream << "bj: \n" << bj << endl;
    curentStream.close();
    curentStream.open("data/cj_ofs.txt");
    curentStream << "cj: \n" << cj << endl;
    curentStream.close();


    //----------------------------------------------------------------------
    // z and Z in final complex form
    //----------------------------------------------------------------------
    for(int l = -nf; l<=nf;l++) bjc.setCoef(bj.getCoef(l), l);
    for(int l = -nf; l<=nf;l++) cjc.setCoef(cj.getCoef(l), l);


    cout << "QBTBP: z and Z have been obtained in final complex form." << endl;

    //----------------------------------------------------------------------
    // The alphas
    //----------------------------------------------------------------------
    //alpha6
    //-----------------------------
    Ofs<double> alpha6(nf);
    //z4 = 1/r2 = z^(-1/2)*zbar^(-1/2)
    z4.zero();
    z4.sprod(z3, zspare);
    fts2fs(&alpha6, z4, eps);

    curentStream.open("data/alpha6.txt");
    curentStream << "alpha6: \n" << alpha6 << endl;
    curentStream.close();
    curentStream.open("data/alpha6_even.txt");
    curentStream << 0 << " " << alpha6.getCoef(0) << endl;
    for(int l = 1; l<=nf;l++)
    curentStream << l << " " << alpha6.getCoef(-l) + alpha6.getCoef(l) << endl;
    curentStream.close();


    //alpha1
    //-----------------------------
    Ofs<double> alpha1(nf);

    //alpha1 = alpha6*alpha6
    alpha1.sprod(alpha6,alpha6);

    curentStream.open("data/alpha1.txt");
    curentStream << "alpha1: \n" << alpha1 << endl;
    curentStream.close();
    curentStream.open("data/alpha1_even.txt");
    curentStream << 0 << " " << alpha1.getCoef(0) << endl;
    for(int l = 1; l<=nf;l++)
    curentStream << l << " " << alpha1.getCoef(-l) + alpha1.getCoef(l) << endl;
    curentStream.close();



    //alpha2 & alpha3
    //-----------------------------
    Ofs<double complex > alpha1c(nf);
    Ofs<double complex > alpha2(nf);
    Ofs<double complex > alpha3(nf);

    Ofs<double complex > zdot(nf);
    Ofs<double complex > zbar(nf);
    Ofs<double complex > zdzb(nf);
    Ofs<double complex > zdzbbar(nf);
    Ofs<double complex > minus_rezdzb(nf);
    Ofs<double complex > imzdzb(nf);

    //zdot over i
    for(int l = -nf; l<=nf;l++) zdot.setCoef(I*(1+l*n)*bj.getCoef(l), l);
    //zbar
    for(int l = -nf; l<=nf;l++) zbar.setCoef(bj.getCoef(l), l);
    zbar.conjugate();
    zdzb.sprod(zdot, zbar);
    //alpha1c
    for(int l = -nf; l<=nf;l++) alpha1c.setCoef(alpha1.getCoef(l), l);
    //zdzbbar = conj(zdzb)
    zdzbbar.ccopy(zdzb);
    zdzbbar.conjugate();
    //minus_rezdzb = -0.5*(zdzb + conj(zdzb))
    for(int l = -nf; l<=nf;l++) minus_rezdzb.setCoef(-0.5*(zdzb.getCoef(l) + zdzbbar.getCoef(l)), l);
    //imzdzb = -0.5i*(zdzb - conj(zdzb))
    for(int l = -nf; l<=nf;l++) imzdzb.setCoef(-0.5*I*(zdzb.getCoef(l) - zdzbbar.getCoef(l)), l);

    //alpha2 = minus_rezdzb*1/r2
    alpha2.sprod(alpha1c, minus_rezdzb);

    curentStream.open("data/alpha2.txt");
    curentStream << "alpha2: \n" << alpha2 << endl;
    curentStream.close();

    curentStream.open("data/alpha2_odd.txt");
    for(int l = 0; l<=nf;l++)
    curentStream << l << " " << creal(-alpha2.getCoef(-l) + alpha2.getCoef(l)) << " " <<  cimag(-alpha2.getCoef(-l) + alpha2.getCoef(l)) << endl;
    curentStream.close();

    //alpha3 = imzdzb*1/r2
-   alpha3.sprod(alpha1c, imzdzb);

    curentStream.open("data/alpha3.txt");
    curentStream << "alpha3: \n" << alpha3 << endl;
    curentStream.close();

    curentStream.open("data/alpha3_even.txt");
    curentStream << 0 << " " << creal(alpha3.getCoef(0)) << " " <<  cimag(alpha3.getCoef(0)) << endl;
    for(int l = 1; l<=nf;l++)
    curentStream << l << " " << creal(alpha3.getCoef(-l) + alpha3.getCoef(l)) << " " <<  cimag(alpha3.getCoef(-l) + alpha3.getCoef(l)) << endl;
    curentStream.close();


    //alpha4 & alpha5
    //-----------------------------
    Ofs<double complex > ddZt(nf);
    Ofs<double complex > ddZ(nf);
    Ofs<double complex > sigma1c(nf);

    Ofs<double complex > ddZzb(nf);
    Ofs<double complex > ddZzbbar(nf);
    Ofs<double complex > re_ddZzb(nf);
    Ofs<double complex > im_ddZzb(nf);


    Ofs<double complex > alpha4(nf);
    Ofs<double complex > alpha5(nf);


    //double dot{Z} temp
    for(int l = -nf; l<=nf;l++) ddZt.setCoef(-as*(ns+l*n)*(ns+l*n)*cj.getCoef(l), l);
     //sigma1c
    for(int l = -nf; l<=nf;l++) sigma1c.setCoef(sigma1.getCoef(l), l);
    //Final double dot{Z}
    ddZ.sprod(ddZt, sigma1c);
    //ddZzb = ddZ*zbar
    ddZzb.sprod(ddZ, zbar);
    //Conjugate
    ddZzbbar.ccopy(ddZzb);
    ddZzbbar.conjugate();
    //Imaginary and real part
    for(int l = -nf; l<=nf;l++) re_ddZzb.setCoef(0.5*(ddZzb.getCoef(l) + ddZzbbar.getCoef(l)), l);
    for(int l = -nf; l<=nf;l++) im_ddZzb.setCoef(-0.5*I*(ddZzb.getCoef(l) - ddZzbbar.getCoef(l)), l);

    //alpha4 & alpha5
    alpha4+= (complex double)-ms/(1+ms)*re_ddZzb;
    alpha5+= (complex double)-ms/(1+ms)*im_ddZzb;

    curentStream.open("data/alpha4.txt");
    curentStream << "alpha4: \n" << alpha4 << endl;
    curentStream.close();

    curentStream.open("data/alpha5.txt");
    curentStream << "alpha5: \n" << alpha5 << endl;
    curentStream.close();

    curentStream.open("data/alpha4_even.txt");
    curentStream << 0 << " " << creal(alpha4.getCoef(0)) << " " <<  cimag(alpha4.getCoef(0)) << endl;
    for(int l = 1; l<=nf;l++)
    curentStream << l << " " << creal(alpha4.getCoef(-l) + alpha4.getCoef(l)) << " " <<  cimag(alpha4.getCoef(-l) + alpha4.getCoef(l)) << endl;
    curentStream.close();

    curentStream.open("data/alpha5_odd.txt");
    for(int l = 0; l<=nf;l++)
    curentStream << l << " " << creal(-alpha5.getCoef(-l) + alpha5.getCoef(l)) << " " <<  cimag(-alpha5.getCoef(-l) + alpha5.getCoef(l)) << endl;
    curentStream.close();


    //alpha7 & alpha8
    //-----------------------------
    Ofs<double complex > Zct(nf);
    Ofs<double complex > Zc(nf);


    Ofs<double complex > Zczb(nf);
    Ofs<double complex > Zczbbar(nf);
    Ofs<double complex > re_Zczb(nf);
    Ofs<double complex > im_Zczb(nf);


    Ofs<double complex > alpha7(nf);
    Ofs<double complex > alpha8(nf);


    //Z temp
    for(int l = -nf; l<=nf;l++) Zct.setCoef(as*cj.getCoef(l), l);
    //Final Z
    Zc.sprod(Zct, sigma1c);
    //Zczb = Zc*zbar
    Zczb.sprod(Zc, zbar);
    //Conjugate
    Zczbbar.ccopy(Zczb);
    Zczbbar.conjugate();
    //Imaginary and real part
    for(int l = -nf; l<=nf;l++) re_Zczb.setCoef(0.5*(Zczb.getCoef(l) + Zczbbar.getCoef(l)), l);
    for(int l = -nf; l<=nf;l++) im_Zczb.setCoef(-0.5*I*(Zczb.getCoef(l) - Zczbbar.getCoef(l)), l);

    //alpha7 & alpha8
    alpha7.sprod(alpha1c, re_Zczb);
    alpha8.sprod(alpha1c, im_Zczb);

    curentStream.open("data/alpha7.txt");
    curentStream << "alpha7: \n" << alpha7 << endl;
    curentStream.close();

    curentStream.open("data/alpha8.txt");
    curentStream << "alpha8: \n" << alpha8 << endl;
    curentStream.close();

    curentStream.open("data/alpha7_even.txt");
    curentStream << 0 << " " << creal(alpha7.getCoef(0)) << " " <<  cimag(alpha7.getCoef(0)) << endl;
    for(int l = 1; l<=nf;l++)
    curentStream << l << " " << creal(alpha7.getCoef(-l) + alpha7.getCoef(l)) << " " <<  cimag(alpha7.getCoef(-l) + alpha7.getCoef(l)) << endl;
    curentStream.close();

    curentStream.open("data/alpha8_odd.txt");
    for(int l = 0; l<=nf;l++)
    curentStream << l << " " << creal(-alpha8.getCoef(-l) + alpha8.getCoef(l)) << " " <<  cimag(-alpha8.getCoef(-l) + alpha8.getCoef(l)) << endl;
    curentStream.close();


    cout << "QBTBP: end of the computation of the alphas" << endl;
}

// Note: the use of gsl livrary forces us to use double variables
// As a consequence, in the case of z and Z, we need to use real and imaginary parts
int qbtbp_derivatives(double t, const double y[], double f[], void *params)
{
    //Parameters for the qbtbp
    double mu = * (double * ) params;
    double ms = * ((double * ) (params)+1);

    //Reconstruction of z and Z
    double complex z = y[0] + I*y[1];
    double complex Z = y[2] + I*y[3];

    double complex temp1 = Z-mu*z;
    temp1 = temp1/pow(cabs(temp1), 3.0);

    double complex temp2 = Z+(1-mu)*z;
    temp2 = temp2/pow(cabs(temp2), 3.0);

    double complex zdd = -z/pow(cabs(z), 3.0) + ms*(temp1-temp2);
    double complex Zdd = -(1+ms)*(mu*temp2 + (1-mu)*temp1);

    //-------------------------------------------------------------------------------
    //Phase space derivatives
    //-------------------------------------------------------------------------------
    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];
    f[4] = creal(zdd);
    f[5] = cimag(zdd);
    f[6] = creal(Zdd);
    f[7] = cimag(Zdd);

    return GSL_SUCCESS;
}

void analyticVsnumeric(double t1, Ofs<double complex> &bjc, Ofs<double complex> &cjc, custom_ode_structure ode_s, double as, double n, double ns)
{
    int nf = bjc.getOrder();

    //z(0) and Z(0)
    double complex z0 = bjc.evaluate(0.0);
    double complex Z0 = as*cjc.evaluate(0.0);

    //zdot(0) and Zdot(0)
    Ofs<double complex> zdot(nf);
    for(int l = -nf; l<=nf; l++) zdot.setCoef(I*(1+l*n)*bjc.getCoef(l), l);
    double complex zdot0 = zdot.evaluate(0.0);

    Ofs<double complex> Zdot(nf);
    for(int l = -nf; l<=nf; l++) Zdot.setCoef(I*(ns+l*n)*cjc.getCoef(l), l);
    double complex Zdot0 = as*Zdot.evaluate(0.0);


    //Initial conditions
    double yv[8];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    cout << "-------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << "Initial positions: " << endl;
    cout << "z(t=0.0): " << creal(z0) << " + " << cimag(z0) << "i" <<  endl;
    cout << "Z(t=0.0): " << creal(Z0) << " + " << cimag(Z0) << "i" <<  endl;


    //Loop
    double h = 1e-6; //first guess for stepper
    int status;
    double t = 0.0;
    while (t < t1)
    {
        //Evolve one step
        status = gsl_odeiv2_evolve_apply (ode_s.e, ode_s.c, ode_s.s, &ode_s.sys, &t, t1, &h, yv);
        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;
    }


    //Final state
    cout << "-------------------------------------------" << endl;
    cout << "End of integration." << endl;
    cout << "Final t: " << t << endl << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Numerical results:" << endl;
    cout << "z(t=t1): " << yv[0] << " + " << yv[1] << "i" << endl;
    cout << "Z(t=t1): " << yv[2] << " + " << yv[3] << "i" << endl;


    //Analytical final state
    double complex zfinal = (cos(t1)+I*sin(t1))*bjc.evaluate(n*t1);
    double complex Zfinal = as*(cos(ns*t1)+I*sin(ns*t1))*cjc.evaluate(n*t1);

    double complex zdotfinal = (cos(t1)+I*sin(t1))*zdot.evaluate(n*t1);
    double complex Zdotfinal = as*(cos(ns*t1)+I*sin(ns*t1))*Zdot.evaluate(n*t1);

    cout << "Analytical results" << endl;
    cout << "z(t=t1): " << creal(zfinal) << " + " << cimag(zfinal) << "i" <<  endl;
    cout << "Z(t=t1): " << creal(Zfinal) << " + " << cimag(Zfinal) << "i" <<  endl;

    cout << "Discrepancy between analytical and numerical results: " << endl;
    cout << "dz    = " << cabs(zfinal - yv[0]-I*yv[1])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dzdot = " << cabs(zdotfinal - yv[4]-I*yv[5])/*/cabs(zdotfinal)*100 << " %" */<< endl;
    cout << "dZ    = " << cabs(Zfinal - yv[2]-I*yv[3])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dZdot = " << cabs(Zdotfinal - yv[6]-I*yv[7])/*/cabs(Zdotfinal)*100 << " %" */<< endl;




}
