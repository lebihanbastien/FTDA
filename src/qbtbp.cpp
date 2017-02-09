/**
 * \file qbtbp.cpp
 * \brief Routines of the Quasi-Bicircular Three-Body Problem (QBTBP) (src)
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */
#include "qbtbp.h"

//-----------------------------------------------------------------------------
// Main routine: computation of the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Main routine to compute the QBTBP in Ofs format. The iteratives equation of the QBTBP are solved using qbtbp_ofs.
 *   If necessary, the result is tested using qbtbp_test.
 *  \param isTestOn: if true, the following tests are performed:
 *                   - a test on a full period is performed vs the numerical integration of the QBTBP.
 *                   - a comparison between the coefficients obtained with FFT and the ones obtained via algebraic manipulations on Ofs objects.
 *                   - a test of the results in Earth-Moon/Inertial/Sun-Earth frameworks.
 */
void qbtbp(int li_EM, int li_SEM, int isTestOn, int coordsys)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Resolution of the                   " << endl;
    cout << "   Quasi-Bicircular Three-Body Problem (QBTBP)     " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << " qbtbp. The qbtbp will be computed up to order " << OFS_ORDER << endl;

    //-------------------------------------
    //Initialization the integration tools
    //-------------------------------------
    //Init the qbcp
    QBCP qbcp;
    init_QBCP(&qbcp, SUN, EARTH, MOON, M_QBCP);

    //Init the qbcp focused on one libration point
    QBCP_L qbcp_l;
    init_QBCP_L(&qbcp_l, &qbcp, true, li_EM, li_SEM, true, M_QBCP, coordsys, PMS_GRAPH, MAN_CENTER, MAN_CENTER); //note PM style and type are not used

    //Init of the internal/external motions
    Ofts<Ofsd> ofts_z(1,OFS_ORDER,2,OFS_ORDER);
    Ofts<Ofsd> ofts_Z(1,OFS_ORDER,2,OFS_ORDER);

    cout << " qbtbp. end of initialization." << endl;

    //-------------------------------------
    //Computation of the QBTBP
    //-------------------------------------
    tic();
    qbtbp_ofs(ofts_z, ofts_Z, qbcp_l, coordsys);
    cout << " qbtbp. end of computation in: " << toc() << "s." << endl;

    //If the user wants to test the results
    if(isTestOn)
    {
        //----------------------------
        //Init the integration tools
        //----------------------------
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        OdeStruct ode_s;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method
        //Ode solver parameters
        double param[2];
        param[0] = qbcp_l.us_em.mu_EM; //note that the qbtbp is computed in EM framework
        param[1] = qbcp_l.us_em.ms;    //note that the qbtbp is computed in EM framework
        //General structures
        init_ode_structure(&ode_s, T, T_root, Config::configManager().G_PREC_ABS(),
        Config::configManager().G_PREC_REL(), Config::configManager().G_PREC_ROOT(),
        8, Config::configManager().G_PREC_HSTART(),
        qbtbp_derivatives, NULL, param);

        //----------------------------
        //Init the external/internal
        //motion in Ofs format
        //----------------------------
        Ofsd bj(OFS_ORDER);
        Ofsd cj(OFS_ORDER);
        Ofsc bjc(OFS_ORDER);
        Ofsc cjc(OFS_ORDER);

        //The epsilon paramater
        double eps = 1.0/qbcp_l.us_em.as;  //note that the qbtbp is computed in EM framework

        //From ots to ofs
        fts2fs(&bj, ofts_z, eps);
        fts2fs(&cj, ofts_Z, eps);

        //From double to cdouble
        doubleToComplex(bj, bjc);
        doubleToComplex(cj, cjc);

        //----------------------------
        //Analytical Vs Numerical results
        //----------------------------
        char ch;
        cout << " qbtbp. The test is performed on one full period by default." << endl;
        qbtbp_test(2*M_PI, bjc, cjc, ode_s, qbcp_l);
        cout << "Press Enter to proceed with the tests." << endl;
        scanf("%c",&ch);
        qbtbp_test_FFT_vs_OFS(bjc, cjc, OFS_ORDER, 500, 1, ode_s, qbcp_l);
        cout << "Press Enter to proceed with the tests." << endl;
        scanf("%c",&ch);
        qbtbp_test_IN_EM_SEM(0.0*M_PI, bjc, cjc);
        cout << "Press Enter to end the test session." << endl;
        scanf("%c",&ch);

    }
}

//-----------------------------------------------------------------------------
// Main routine: computation of the BCP
//-----------------------------------------------------------------------------
void bcp_alpha(QBCP_L &qbcp_l)
{
    //--------------------------------------------------------------------
    //Parameters for one specific libration point
    //--------------------------------------------------------------------
    double mu_EM = qbcp_l.us_em.mu_EM;
    double gamma = qbcp_l.cs_em.gamma;
    double c1    = qbcp_l.cs_em.c1;
    double ms    = qbcp_l.us_em.ms;
    double as    = qbcp_l.us_em.as;
    int nf       = qbcp_l.nf;

    //--------------------------------------------------------------------
    // 3. Set all alpha functions
    //--------------------------------------------------------------------
    Ofs<cdouble > alpha1c(nf);
    Ofs<cdouble > alpha2c(nf);
    Ofs<cdouble > alpha3c(nf);
    Ofs<cdouble > alpha4c(nf);
    Ofs<cdouble > alpha5c(nf);
    Ofs<cdouble > alpha6c(nf);
    Ofs<cdouble > alpha7c(nf);
    Ofs<cdouble > alpha8c(nf);
    Ofs<cdouble > alpha9c(nf);
    Ofs<cdouble > alpha10c(nf);
    Ofs<cdouble > alpha11c(nf);
    Ofs<cdouble > alpha12c(nf);
    Ofs<cdouble > alpha13c(nf);
    Ofs<cdouble > alpha14c(nf);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(nf);
    Ofs<cdouble > Ye(nf);
    Ofs<cdouble > Ze(nf);
    Ofs<cdouble > Xm(nf);
    Ofs<cdouble > Ym(nf);
    Ofs<cdouble > Zm(nf);
    Ofs<cdouble > Xs(nf);
    Ofs<cdouble > Ys(nf);
    Ofs<cdouble > Zs(nf);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(nf);
    Ofs<cdouble > ye(nf);
    Ofs<cdouble > ze(nf);
    Ofs<cdouble > xm(nf);
    Ofs<cdouble > ym(nf);
    Ofs<cdouble > zm(nf);
    Ofs<cdouble > xs(nf);
    Ofs<cdouble > ys(nf);
    Ofs<cdouble > zs(nf);

    //-------------------------
    // The alphas
    //-------------------------
    //alpha1 = 1
    alpha1c.setCoef(1.0+0.0*I, 0);
    //alpha2 = 0
    //alpha3 = 1
    alpha3c.setCoef(1.0+0.0*I, 0);
    //alpha4 = ms/as^2*cos(theta) = ms/(2*as^2)*(exp(itheta) + exp(-itheta))
    alpha4c.setCoef(0.5*ms/(as*as)+0.0*I, -1);
    alpha4c.setCoef(0.5*ms/(as*as)+0.0*I, +1);
    //alpha5 = -ms/as^2*sin(theta) = +I*ms/(2*as^2)*(exp(itheta) - exp(-itheta))
    alpha5c.setCoef(-0.5*ms/(as*as)*I, -1);
    alpha5c.setCoef(+0.5*ms/(as*as)*I, +1);
    //alpha6 = 1
    alpha6c.setCoef(1.0+0.0*I, 0);
    //alpha13 = alpha4/gamma-c1
    alpha13c.ofs_mult(alpha4c, 1.0/gamma+0.0*I);
    alpha13c.addCoef(-c1+0.0*I, 0);
    //alpha14 = alpha5/gamma
    alpha14c.ofs_mult(alpha5c, 1.0/gamma+0.0*I);

    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    alpha9c.setCoef(+mu_EM+0.0*I, 0);
    //Moon
    alpha11c.setCoef(+mu_EM-1.0+0.0*I, 0);
    //Sun
    //alpha7 = as*cos(theta) = as/2*(exp(itheta) + exp(-itheta))
    alpha7c.setCoef(0.5*as+0.0*I, -1);
    alpha7c.setCoef(0.5*as+0.0*I, +1);
    //alpha8 = -as*sin(theta) = +I*as/2*(exp(itheta) - exp(-itheta))
    alpha8c.setCoef(-0.5*as*I, -1);
    alpha8c.setCoef(+0.5*as*I, +1);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    Xe.setCoef(+mu_EM+0.0*I, 0);
    //Earth, NC coordinates
    xe.setCoef(c1 - mu_EM/gamma + 0.0*I, 0);
    //Moon, EM coordinates
    Xm.setCoef(+mu_EM-1.0+0.0*I, 0);
    //Moon, NC coordinates
    xm.setCoef(c1 - (mu_EM-1.0)/gamma + 0.0*I, 0);

    //Sun, EM coordinates
    //---------------
    //Xs = as*cos(theta) = as/2*(exp(itheta) + exp(-itheta))
    Xs.setCoef(0.5*as+0.0*I, -1);
    Xs.setCoef(0.5*as+0.0*I, +1);
    //Ys = -as*sin(theta) = +I*as/2*(exp(itheta) - exp(-itheta))
    Ys.setCoef(-0.5*as*I, -1);
    Ys.setCoef(+0.5*as*I, +1);
    //Sun, NC coordinates;
    //---------------
    xs.ofs_smult(Xs, -1.0/gamma);
    xs.addCoef(c1, 0);
    ys.ofs_smult(Ys, -1.0/gamma);

    //--------------------------
    //Put in data file
    //--------------------------
    cout << "bcp. storage in txt files. " << endl;
    ofs_sst(alpha1c, qbcp_l.cs_em.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(alpha2c, qbcp_l.cs_em.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(alpha3c, qbcp_l.cs_em.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(alpha4c, qbcp_l.cs_em.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(alpha5c, qbcp_l.cs_em.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(alpha6c, qbcp_l.cs_em.F_COEF+"alpha6", 1, "_fft");
     //Sun
    ofs_sst(alpha7c, qbcp_l.cs_em.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(alpha8c, qbcp_l.cs_em.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(alpha9c,  qbcp_l.cs_em.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(alpha10c, qbcp_l.cs_em.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(alpha11c, qbcp_l.cs_em.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(alpha12c, qbcp_l.cs_em.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(alpha13c, qbcp_l.cs_em.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(alpha14c, qbcp_l.cs_em.F_COEF+"alpha14", 0, "_fft");

    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //---------------
    ofs_sst(Xs, qbcp_l.cs_em.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, qbcp_l.cs_em.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, qbcp_l.cs_em.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, qbcp_l.cs_em.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, qbcp_l.cs_em.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, qbcp_l.cs_em.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, qbcp_l.cs_em.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, qbcp_l.cs_em.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, qbcp_l.cs_em.F_COEF+"Pe3", 1, "_fft");


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //---------------
    ofs_sst(xs, qbcp_l.cs_em.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, qbcp_l.cs_em.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, qbcp_l.cs_em.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, qbcp_l.cs_em.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, qbcp_l.cs_em.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, qbcp_l.cs_em.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, qbcp_l.cs_em.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, qbcp_l.cs_em.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, qbcp_l.cs_em.F_COEF+"pe3", 1, "_fft");
}

void bcp_delta(QBCP_L &qbcp_l)
{
    //--------------------------------------------------------------------
    //Parameters for one specific libration point
    //--------------------------------------------------------------------
    double mu_SE  = qbcp_l.us_sem.mu_SEM;
    double gamma  = qbcp_l.cs_sem.gamma;
    double c1     = qbcp_l.cs_sem.c1;
    //double mm     = qbcp_l.us_sem.mm;
    double am     = qbcp_l.us_sem.ai;
    int nf        = qbcp_l.nf;

    //--------------------------------------------------------------------
    // 3. Set all delta functions
    //--------------------------------------------------------------------
    Ofs<cdouble > delta1c(nf);
    Ofs<cdouble > delta2c(nf);
    Ofs<cdouble > delta3c(nf);
    Ofs<cdouble > delta4c(nf);
    Ofs<cdouble > delta5c(nf);
    Ofs<cdouble > delta6c(nf);
    Ofs<cdouble > delta7c(nf);
    Ofs<cdouble > delta8c(nf);
    Ofs<cdouble > delta9c(nf);
    Ofs<cdouble > delta10c(nf);
    Ofs<cdouble > delta11c(nf);
    Ofs<cdouble > delta12c(nf);
    Ofs<cdouble > delta13c(nf);
    Ofs<cdouble > delta14c(nf);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(nf);
    Ofs<cdouble > Ye(nf);
    Ofs<cdouble > Ze(nf);
    Ofs<cdouble > Xm(nf);
    Ofs<cdouble > Ym(nf);
    Ofs<cdouble > Zm(nf);
    Ofs<cdouble > Xs(nf);
    Ofs<cdouble > Ys(nf);
    Ofs<cdouble > Zs(nf);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(nf);
    Ofs<cdouble > ye(nf);
    Ofs<cdouble > ze(nf);
    Ofs<cdouble > xm(nf);
    Ofs<cdouble > ym(nf);
    Ofs<cdouble > zm(nf);
    Ofs<cdouble > xs(nf);
    Ofs<cdouble > ys(nf);
    Ofs<cdouble > zs(nf);

    //-------------------------
    // The deltas
    //-------------------------
    //delta1 = 1
    delta1c.setCoef(1.0+0.0*I, 0);
    //delta2 = 0
    //delta3 = 1
    delta3c.setCoef(1.0+0.0*I, 0);
    //delta4 = -mm/am^2*cos(theta) = -mm/(2*am^2)*(exp(itheta) + exp(-itheta))
    //delta4c.setCoef(-0.5*mm/(am*am)+0.0*I, -1);
    //delta4c.setCoef(-0.5*mm/(am*am)+0.0*I, +1);
    //delta5 = -mm/am^2*sin(theta) = +I*mm/(2*am^2)*(exp(itheta) - exp(-itheta))
    //delta5c.setCoef(-0.5*mm/(am*am)*I, -1);
    //delta5c.setCoef(+0.5*mm/(am*am)*I, +1);
    //delta6 = 1
    delta6c.setCoef(1.0+0.0*I, 0);
    //delta13 = delta4/gamma-c1
    delta13c.ofs_mult(delta4c, 1.0/gamma+0.0*I);
    delta13c.addCoef(-c1+0.0*I, 0);
    //delta14 = delta5/gamma
    delta14c.ofs_mult(delta5c, 1.0/gamma+0.0*I);

    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    delta9c.setCoef(+mu_SE-1.0+0.0*I, 0);
    //Sun
    delta7c.setCoef(+mu_SE+0.0*I, 0);
    //Moon
    //delta11 = -am*cos(theta) + mu_SE-1.0 = -am/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    delta11c.setCoef(mu_SE-1.0 +0.0*I, 0);
    delta11c.setCoef(-0.5*am+0.0*I, -1);
    delta11c.setCoef(-0.5*am+0.0*I, +1);
    //delta12 = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    delta12c.setCoef(-0.5*am*I, -1);
    delta12c.setCoef(+0.5*am*I, +1);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    Xe.setCoef(mu_SE-1.0+0.0*I, 0);
    //Earth, NC coordinates
    xe.setCoef(c1 - (mu_SE-1.0)/gamma + 0.0*I, 0);
    //Sun, EM coordinates
    Xs.setCoef(+mu_SE+0.0*I, 0);
    //Sun, NC coordinates
    xs.setCoef(c1 - mu_SE/gamma + 0.0*I, 0);

    //Moon, EM coordinates
    //---------------
    //Xm = -am*cos(theta) + mu_SE - 1 = am/2*(exp(itheta) + exp(-itheta)) + mu_SE - 1
    Xm.setCoef(mu_SE-1.0 +0.0*I, 0);
    Xm.setCoef(-0.5*am+0.0*I, -1);
    Xm.setCoef(-0.5*am+0.0*I, +1);
    //Ym = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    Ym.setCoef(-0.5*am*I, -1);
    Ym.setCoef(+0.5*am*I, +1);
    //Moon, NC coordinates;
    //---------------
    xm.ofs_smult(Xm, -1.0/gamma);
    xm.addCoef(c1, 0);
    ym.ofs_smult(Ym, -1.0/gamma);

    //--------------------------
    //Put in data file
    //--------------------------
    cout << "bcp. storage in txt files. " << endl;
    ofs_sst(delta1c, qbcp_l.cs_sem.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(delta2c, qbcp_l.cs_sem.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(delta3c, qbcp_l.cs_sem.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(delta4c, qbcp_l.cs_sem.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(delta5c, qbcp_l.cs_sem.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(delta6c, qbcp_l.cs_sem.F_COEF+"alpha6", 1, "_fft");
     //Sun
    ofs_sst(delta7c, qbcp_l.cs_sem.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(delta8c, qbcp_l.cs_sem.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(delta9c,  qbcp_l.cs_sem.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(delta10c, qbcp_l.cs_sem.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(delta11c, qbcp_l.cs_sem.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(delta12c, qbcp_l.cs_sem.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(delta13c, qbcp_l.cs_sem.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(delta14c, qbcp_l.cs_sem.F_COEF+"alpha14", 0, "_fft");

    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //---------------
    ofs_sst(Xs, qbcp_l.cs_sem.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, qbcp_l.cs_sem.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, qbcp_l.cs_sem.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, qbcp_l.cs_sem.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, qbcp_l.cs_sem.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, qbcp_l.cs_sem.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, qbcp_l.cs_sem.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, qbcp_l.cs_sem.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, qbcp_l.cs_sem.F_COEF+"Pe3", 1, "_fft");


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //---------------
    ofs_sst(xs, qbcp_l.cs_sem.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, qbcp_l.cs_sem.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, qbcp_l.cs_sem.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, qbcp_l.cs_sem.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, qbcp_l.cs_sem.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, qbcp_l.cs_sem.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, qbcp_l.cs_sem.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, qbcp_l.cs_sem.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, qbcp_l.cs_sem.F_COEF+"pe3", 1, "_fft");
}

void bcp_delta_2(QBCP_L &qbcp_l)
{
    //--------------------------------------------------------------------
    //Parameters for one specific libration point
    //--------------------------------------------------------------------
    double mu_EM  = qbcp_l.us_sem.mu_EM;
    double mu_SE  = qbcp_l.us_sem.mu_SEM;
    double gamma  = qbcp_l.cs_sem.gamma;
    double c1     = qbcp_l.cs_sem.c1;
    //double mm     = qbcp_l.us_sem.mm;
    //double me     = qbcp_l.us_sem.me;
    double ai     = qbcp_l.us_sem.ai;
    int nf        = qbcp_l.nf;
    //Distance of Earth & Moon from their barycenter
    double am = (1- mu_EM)*ai;
    double ae = mu_EM*ai;
    //Factor that appears in the definition of the deltas.
    //double fem = mm/(am*am) - me/(ae*ae);

    //--------------------------------------------------------------------
    // 3. Set all delta functions
    //--------------------------------------------------------------------
    Ofs<cdouble > delta1c(nf);
    Ofs<cdouble > delta2c(nf);
    Ofs<cdouble > delta3c(nf);
    Ofs<cdouble > delta4c(nf);
    Ofs<cdouble > delta5c(nf);
    Ofs<cdouble > delta6c(nf);
    Ofs<cdouble > delta7c(nf);
    Ofs<cdouble > delta8c(nf);
    Ofs<cdouble > delta9c(nf);
    Ofs<cdouble > delta10c(nf);
    Ofs<cdouble > delta11c(nf);
    Ofs<cdouble > delta12c(nf);
    Ofs<cdouble > delta13c(nf);
    Ofs<cdouble > delta14c(nf);

    //Redundancy for the positions of the primaries
    Ofs<cdouble > Xe(nf);
    Ofs<cdouble > Ye(nf);
    Ofs<cdouble > Ze(nf);
    Ofs<cdouble > Xm(nf);
    Ofs<cdouble > Ym(nf);
    Ofs<cdouble > Zm(nf);
    Ofs<cdouble > Xs(nf);
    Ofs<cdouble > Ys(nf);
    Ofs<cdouble > Zs(nf);

    //Positions of the primaries in NC coordinates
    Ofs<cdouble > xe(nf);
    Ofs<cdouble > ye(nf);
    Ofs<cdouble > ze(nf);
    Ofs<cdouble > xm(nf);
    Ofs<cdouble > ym(nf);
    Ofs<cdouble > zm(nf);
    Ofs<cdouble > xs(nf);
    Ofs<cdouble > ys(nf);
    Ofs<cdouble > zs(nf);

    //-------------------------
    // The deltas
    //-------------------------
    //delta1 = 1
    delta1c.setCoef(1.0+0.0*I, 0);
    //delta2 = 0
    //delta3 = 1
    delta3c.setCoef(1.0+0.0*I, 0);
    //delta4 = 0; // OR -fem*cos(theta) = -0.5*fem*(exp(itheta) + exp(-itheta))
    //delta4c.setCoef(-0.5*fem+0.0*I, -1);
    //delta4c.setCoef(-0.5*fem+0.0*I, +1);
    //delta5 = 0; // OR -fem*sin(theta) = +I*0.5*fem*(exp(itheta) - exp(-itheta))
    //delta5c.setCoef(-0.5*fem*I, -1);
    //delta5c.setCoef(+0.5*fem*I, +1);
    //delta6 = 1
    delta6c.setCoef(1.0+0.0*I, 0);
    //delta13 = delta4/gamma-c1
    delta13c.ofs_mult(delta4c, 1.0/gamma+0.0*I);
    delta13c.addCoef(-c1+0.0*I, 0);
    //delta14 = delta5/gamma
    delta14c.ofs_mult(delta5c, 1.0/gamma+0.0*I);

    //-------------------------
    // Primaries
    //-------------------------
    //Earth
    //delta9 = +ae*cos(theta) + mu_SE-1.0  = ae/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    delta9c.setCoef(mu_SE-1.0 +0.0*I, 0);
    delta9c.setCoef(0.5*ae+0.0*I, -1);
    delta9c.setCoef(0.5*ae+0.0*I, +1);
    //delta10 = +ae*sin(theta) = -I*ae/2*(exp(itheta) - exp(-itheta))
    delta10c.setCoef(+0.5*ae*I, -1);
    delta10c.setCoef(-0.5*ae*I, +1);

    //Sun
    delta7c.setCoef(+mu_SE+0.0*I, 0);

    //Moon
    //delta11 = -am*cos(theta) + mu_SE-1.0  = -am/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    delta11c.setCoef(mu_SE-1.0 +0.0*I, 0);
    delta11c.setCoef(-0.5*am+0.0*I, -1);
    delta11c.setCoef(-0.5*am+0.0*I, +1);
    //delta12 = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    delta12c.setCoef(-0.5*am*I, -1);
    delta12c.setCoef(+0.5*am*I, +1);

    //-------------------------
    // The primaries, again
    //-------------------------
    //Earth, EM coordinates
    //---------------
    //Xe = +ae*cos(theta) + mu_SE-1.0  = ae/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    Xe.setCoef(mu_SE-1.0 +0.0*I, 0);
    Xe.setCoef(0.5*ae+0.0*I, -1);
    Xe.setCoef(0.5*ae+0.0*I, +1);
    //Ye = +ae*sin(theta) = -I*ae/2*(exp(itheta) - exp(-itheta))
    Ye.setCoef(+0.5*ae*I, -1);
    Ye.setCoef(-0.5*ae*I, +1);
    //Earth, NC coordinates;
    //---------------
    xe.ofs_smult(Xe, -1.0/gamma);
    xe.addCoef(c1, 0);
    ye.ofs_smult(Ye, -1.0/gamma);

    //Sun, EM coordinates
    Xs.setCoef(+mu_SE+0.0*I, 0);
    //Sun, NC coordinates
    xs.setCoef(c1 - mu_SE/gamma + 0.0*I, 0);

    //Moon, EM coordinates
    //---------------
    //Xm = -am*cos(theta) + mu_SE-1.0  = am/2*(exp(itheta) + exp(-itheta)) + mu_SE-1.0
    Xm.setCoef(mu_SE-1.0 +0.0*I, 0);
    Xm.setCoef(-0.5*am+0.0*I, -1);
    Xm.setCoef(-0.5*am+0.0*I, +1);
    //Ym = -am*sin(theta) = +I*am/2*(exp(itheta) - exp(-itheta))
    Ym.setCoef(-0.5*am*I, -1);
    Ym.setCoef(+0.5*am*I, +1);
    //Moon, NC coordinates;
    //---------------
    xm.ofs_smult(Xm, -1.0/gamma);
    xm.addCoef(c1, 0);
    ym.ofs_smult(Ym, -1.0/gamma);

    //--------------------------
    //Put in data file
    //--------------------------
    cout << "bcp. storage in txt files. " << endl;
    ofs_sst(delta1c, qbcp_l.cs_sem.F_COEF+"alpha1", 1, "_fft");
    ofs_sst(delta2c, qbcp_l.cs_sem.F_COEF+"alpha2", 0, "_fft");
    ofs_sst(delta3c, qbcp_l.cs_sem.F_COEF+"alpha3", 1, "_fft");
    ofs_sst(delta4c, qbcp_l.cs_sem.F_COEF+"alpha4", 1, "_fft");
    ofs_sst(delta5c, qbcp_l.cs_sem.F_COEF+"alpha5", 0, "_fft");
    ofs_sst(delta6c, qbcp_l.cs_sem.F_COEF+"alpha6", 1, "_fft");
     //Sun
    ofs_sst(delta7c, qbcp_l.cs_sem.F_COEF+"alpha7", 1, "_fft");
    ofs_sst(delta8c, qbcp_l.cs_sem.F_COEF+"alpha8", 0, "_fft");
    //Earth
    ofs_sst(delta9c,  qbcp_l.cs_sem.F_COEF+"alpha9",  1, "_fft");
    ofs_sst(delta10c, qbcp_l.cs_sem.F_COEF+"alpha10", 0, "_fft");
    //Moon
    ofs_sst(delta11c, qbcp_l.cs_sem.F_COEF+"alpha11", 1, "_fft");
    ofs_sst(delta12c, qbcp_l.cs_sem.F_COEF+"alpha12", 0, "_fft");
    //NC additional coeffs
    ofs_sst(delta13c, qbcp_l.cs_sem.F_COEF+"alpha13", 1, "_fft");
    ofs_sst(delta14c, qbcp_l.cs_sem.F_COEF+"alpha14", 0, "_fft");

    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_Zc without much difference
    //---------------
    ofs_sst(Xs, qbcp_l.cs_sem.F_COEF+"Ps1", 1, "_fft");
    ofs_sst(Ys, qbcp_l.cs_sem.F_COEF+"Ps2", 0, "_fft");
    ofs_sst(Zs, qbcp_l.cs_sem.F_COEF+"Ps3", 1, "_fft");

    ofs_sst(Xm, qbcp_l.cs_sem.F_COEF+"Pm1", 1, "_fft");
    ofs_sst(Ym, qbcp_l.cs_sem.F_COEF+"Pm2", 0, "_fft");
    ofs_sst(Zm, qbcp_l.cs_sem.F_COEF+"Pm3", 1, "_fft");

    ofs_sst(Xe, qbcp_l.cs_sem.F_COEF+"Pe1", 1, "_fft");
    ofs_sst(Ye, qbcp_l.cs_sem.F_COEF+"Pe2", 0, "_fft");
    ofs_sst(Ze, qbcp_l.cs_sem.F_COEF+"Pe3", 1, "_fft");


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the ofs_sst option of gsl_zc without much difference
    //---------------
    ofs_sst(xs, qbcp_l.cs_sem.F_COEF+"ps1", 1, "_fft");
    ofs_sst(ys, qbcp_l.cs_sem.F_COEF+"ps2", 0, "_fft");
    ofs_sst(zs, qbcp_l.cs_sem.F_COEF+"ps3", 1, "_fft");

    ofs_sst(xm, qbcp_l.cs_sem.F_COEF+"pm1", 1, "_fft");
    ofs_sst(ym, qbcp_l.cs_sem.F_COEF+"pm2", 0, "_fft");
    ofs_sst(zm, qbcp_l.cs_sem.F_COEF+"pm3", 1, "_fft");

    ofs_sst(xe, qbcp_l.cs_sem.F_COEF+"pe1", 1, "_fft");
    ofs_sst(ye, qbcp_l.cs_sem.F_COEF+"pe2", 0, "_fft");
    ofs_sst(ze, qbcp_l.cs_sem.F_COEF+"pe3", 1, "_fft");
}


/**
 *  \brief Main routine to compute the Bicircular Three-Body Problem in Ofs format.
 */
void bcp(int li_EM, int li_SEM, int coordsys)
{

    //--------------------------------------------------------------------
    // 1. Init
    //--------------------------------------------------------------------
    //Init the fbp
    QBCP fbp;
    init_QBCP(&fbp, SUN, EARTH, MOON, M_BCP);

    //Init the fbp focused on one libration point
    QBCP_L qbcp_l;
    init_QBCP_L(&qbcp_l, &fbp, 1, li_EM, li_SEM, true, M_BCP, coordsys, PMS_GRAPH, MAN_CENTER, MAN_CENTER);  //Note: PM style is NOT used

    //--------------------------------------------------------------------
    // 2. Splash
    //--------------------------------------------------------------------
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "              Storage of the                       " << endl;
    cout << "    Bicircular Three-Body Problem (BCP)            " << endl;
    cout << "                                                   " << endl;
    cout << "      save in " << qbcp_l.cs_em.F_COEF               << endl;
    cout << "      save in " << qbcp_l.cs_sem.F_COEF              << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //--------------------------------------------------------------------
    // 3. Alphas
    //--------------------------------------------------------------------
    bcp_alpha(qbcp_l);

    //--------------------------------------------------------------------
    // 4. Deltas
    //--------------------------------------------------------------------
    //bcp_delta(qbcp_l); //libration point on the Earth-Sun line
    bcp_delta_2(qbcp_l); //libration points on the Sun-Bem line
}

//-----------------------------------------------------------------------------
// Computing the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Computes the QBTBP in Ofs format.
 *
 *   The iteratives equation of the QBTBP are solved. The periodic functions \f$ \alpha_i \f$ (from the Earth-Moon p.o.v) and
 *   \f$ \beta_i \f$ (from the Sun-(Earth+Moon) p.o.v) are computed in stored in txt files.
 *   These functions are computed both through operations on Fourier series and FFT of the integrated motion.
 */
void qbtbp_ofs (Ofts< Ofsd > &zr_ofts, Ofts< Ofsd > &Zr_ofts, QBCP_L& qbcp_l, int coordsys)
{
    //---------------------------------------------------------
    //Parameters
    //---------------------------------------------------------
    int order  = zr_ofts.getOrder();         //order of the Taylor expansion
    int nv     = zr_ofts.getNV();            //number of variables of the Taylor expansion
    int nf     = zr_ofts.getCOrder();        //order of the Fourier coefficient
    int cnv    = zr_ofts.getCVariables();    //number of variables of the Fourier coefficient

    //Physical param
    double n  = qbcp_l.us_em.n;  //mean angular motion in EM units
    double as = qbcp_l.us_em.as; //Sun-(Earth+Moon barycenter) average distance in EM units
    double ms = qbcp_l.us_em.ms; //Sun mass in EM units

    //---------------------------------------------------------
    //Declaration of the variables
    //---------------------------------------------------------
    //Order 1
    Ofsd u1(nf);
    Ofsd v1(nf);
    //Order n of zr and Zr
    Ofts< Ofsd > z1(nv, order, cnv, nf);     //z1 = \bar{zr_ofts}
    Ofts< Ofsd > z2(nv, order, cnv, nf);     //z2 = \bar{zr_ofts}^(-3/2)
    Ofts< Ofsd > z3(nv, order, cnv, nf);     //z3 = zr_ofts^(-1/2)
    Ofts< Ofsd > z4(nv, order, cnv, nf);     //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
    Ofts< Ofsd > z5(nv, order, cnv, nf);     //z5 = \bar{zr_ofts}^(-1/2)
    Ofts< Ofsd > z6(nv, order, cnv, nf);     //z6 = z3*z5 = (zr_ofts*\bar{zr_ofts})^(-1/2)

    Ofts< Ofsd > a1(nv, order, cnv, nf);     //a1 = mu*exp(itheta)*zr_ofts
    Ofts< Ofsd > a2(nv, order, cnv, nf);     //a2 = epsilon*a1
    Ofts< Ofsd > a3(nv, order, cnv, nf);     //a3 = Zr_ofts-mu*exp(itheta)*epsilon*zr_ofts
    Ofts< Ofsd > a4(nv, order, cnv, nf);     //a4 = a3^(-1/2).
    Ofts< Ofsd > a5(nv, order, cnv, nf);     //a5 = \bar{a3}
    Ofts< Ofsd > a6(nv, order, cnv, nf);     //a6 = a5^(-3/2).
    Ofts< Ofsd > a7(nv, order, cnv, nf);     //a7 = a4*a6;

    Ofts< Ofsd > b1(nv, order, cnv, nf);     //b1 = (1-mu)*exp(ithetb)*zr_ofts
    Ofts< Ofsd > b2(nv, order, cnv, nf);     //b2 = epsilon*b1
    Ofts< Ofsd > b3(nv, order, cnv, nf);     //b3 = Zr_ofts+(1-mu)*exp(ithetb)*epsilon*zr_ofts
    Ofts< Ofsd > b4(nv, order, cnv, nf);     //b4 = b3^(-1/2)
    Ofts< Ofsd > b5(nv, order, cnv, nf);     //b5 = \bar{b3}
    Ofts< Ofsd > b6(nv, order, cnv, nf);     //b6 = b5^(-3/2)
    Ofts< Ofsd > b7(nv, order, cnv, nf);     //b7 = b4*b6;

    Ofts< Ofsd > Pm(nv, order, cnv, nf);     //Pm = -ms/as^2*exp(-itheta)*b7 + ms/as^2*exp(-itheta)*a7 - z4
    Ofts< Ofsd > Qm(nv, order, cnv, nf);     //Qm = -ns^2*mu*b7 - ns^2*(1-mu)*a7

    //Ofs used to solve the system at each order n
    Ofsd  Pfm(nf);
    Ofsd  Qfm(nf);
    Ofsd  ufm(nf);
    Ofsd  vfm(nf);

    //---------------------------------------------------------
    //Specific Ofs and Ofts objects
    //---------------------------------------------------------
    //sigma1 = exp(-itheta)
    Ofsd sigma1(nf);
    sigma1.setCoef(1.0, -1);
    //sigma2 = exp(+itheta)
    Ofsd sigma2(nf);
    sigma2.setCoef(1.0, 1);
    //epsilon = Ofts with just 1.0 at order 1
    Ofts< Ofsd > epsilon(nv, order, cnv, nf);
    epsilon.setCoef0(1.0, 1, 0);

    //---------------------------------------------------------
    //Order 0 of zr and Zr
    //---------------------------------------------------------
    zr_ofts.setCoef0(1.0, 0, 0);
    Zr_ofts.setCoef0(1.0, 0, 0);

    //---------------------------------------------------------
    //Order 1 (optionnal)
    //---------------------------------------------------------
    u1.setCoef(-ms/(as*as)/6.0, 0);
    u1.setCoef(-ms/(as*as)*(24.0*n*n+24*n+9.0)/(64*pow(n,4.0)-16*n*n), -2);
    u1.setCoef(+ms/(as*as)*9.0/(64*pow(n,4.0)-16*n*n), 2);

    //Fourier to Taylor
    //fs2ts(zr_ofts.getCoef(1,0), u1);
    //fs2ts(Zr_ofts.getCoef(1,0), v1);

    //---------------------------------------------------------
    //Order 0 of the recurrence
    //---------------------------------------------------------
    int m = 0;
    qbtbp_ofs_recurrence(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, z6,
                         a1, a2, a3, a4, a5, a6, a7,
                         b1, b2, b3, b4, b5, b6, b7,
                         Pm, Qm, epsilon, Pfm, Qfm, ufm, vfm,
                         sigma1, sigma2, m, nf, qbcp_l);

    //---------------------------------------------------------
    //Recurrence: m = 1 to zr_ofts.getOrder()
    //---------------------------------------------------------
    for(m = 1; m <= zr_ofts.getOrder(); m++)
    {
        qbtbp_ofs_recurrence(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, z6,
                             a1, a2, a3, a4, a5, a6, a7,
                             b1, b2, b3, b4, b5, b6, b7,
                             Pm, Qm, epsilon, Pfm, Qfm, ufm, vfm,
                             sigma1, sigma2, m, nf, qbcp_l);
    }
    cout << "   qbtbp_ofs. end of recursive computation." << endl;

    //---------------------------------------------------------
    //Manipulation and Storage in txt files of the OFS
    //---------------------------------------------------------

    //---------------------
    //Computation of the vector field coefficient for the EMLi libration point of qbcp_l
    //---------------------
    qbtbp_ofs_fft_alpha(zr_ofts, Zr_ofts, nf, qbcp_l);

    //---------------------
    //Computation of the vector field coefficient for the SEMLi libration point of qbcp_l
    //---------------------
    qbtbp_ofs_fft_delta(zr_ofts, Zr_ofts, nf, qbcp_l);

    //---------------------
    // DEPRECATED.
    // - Computation of the vector field coefficient
    // for the EMLi libration point of qbcp_l through OFS manipulation.
    // - Storage in data files of the general QBTBP.
    // - Note that ONLY the element stored in qbtbp folder can be used
    // the rest is DEPRECATED. Use FFT coefficients instead.
    // TO BE ENHANCED: isolate the computation in qbtbp folder
    //---------------------
    qbtbp_ofs_storage(zr_ofts, Zr_ofts, z1, z2, z3, z4, z5, z6,
                      a1, a2, a3, a4, a5, a6, a7,
                      b1, b2, b3, b4, b5, b6, b7,
                      Pm, Qm, Pfm, Qfm, ufm, vfm,
                      sigma1, sigma2, nf, qbcp_l);


    cout << "   qbtbp_ofs. end of the computation of the alphas and betas" << endl;
}

//-----------------------------------------------------------------------------
// Storing the results.
//-----------------------------------------------------------------------------
/**
 *  \brief Computes one specific FFT and stores it txt file.
 *  \param xGSL0: a GSL vector of size N which contains the discrete evaluations of the function f on which to perform the FFT.
 *  \param filename: a string. The final result is stored in the file filename+"_fft.txt"
 *  \param nf: an integer. The order of the Fourier series obtained after FFT.
 *  \param flag: a boolean. If true, the function f is supposed even (sum of cosinus). If false, f is supposed odd (sum of sinus).
 */
void qbtbp_ofs_fft_unpack(gsl_vector *xGSL0, string filename, int nf, int N, int flag)
{
    //Copy of xGSL0 into temp vector
    gsl_vector *xGSL = gsl_vector_calloc(N);
    gsl_vector_memcpy(xGSL, xGSL0);

    gsl_vector_complex *data_complex = gsl_vector_complex_calloc(N);
    gsl_fft_real_wavetable * wavetable = gsl_fft_real_wavetable_alloc (N);
    gsl_fft_real_workspace * workspace = gsl_fft_real_workspace_alloc (N);
    Ofsc xFFT(nf);
    ofstream curentStream;

    gsl_fft_real_transform (xGSL->data, 1, N, wavetable, workspace);
    gsl_fft_halfcomplex_unpack(xGSL->data , data_complex->data ,  xGSL->stride ,xGSL->size);
    //Order 0
    if(flag) //even case
        xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, 0))/N+I*GSL_IMAG(gsl_vector_complex_get(data_complex, 0))/N,  0);
    else //odd case
        xFFT.setCoef(0.0,  0);

    //Order n
    for(int i = 1; i<= nf; i++)
    {
        if(flag) //even case
        {
            //Negative frequencies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, N-i))/N, -i);
            //Positive frequencies
            xFFT.setCoef(+GSL_REAL(gsl_vector_complex_get(data_complex, i))/N,  i);
        }
        else //odd case
        {
            //Negative frequencies
            xFFT.setCoef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, N-i))/N, -i);
            //Positive frequencies
            xFFT.setCoef(I*GSL_IMAG(gsl_vector_complex_get(data_complex, i))/N,  i);
        }

    }

    //Storage in txt file
    curentStream.open((filename+"_fft.txt").c_str());
    curentStream <<  xFFT << endl;
    curentStream.close();


    if(flag) //even case
    {
        //Cosinus expansion version
        curentStream.open((filename+"c_fft.txt").c_str());
        curentStream << 0 << " " << creal(xFFT.ofs_getCoef(0)) << endl;
        for(int l = 1; l<=nf; l++) curentStream << l << " " << creal(xFFT.ofs_getCoef(-l) + xFFT.ofs_getCoef(l))  <<  endl;
        curentStream.close();
    }
    else //odd case
    {
        //Sinus expansion version
        curentStream.open((filename+"c_fft.txt").c_str());
        for(int l = 0; l<=nf; l++)
            curentStream << l << " " <<  cimag(xFFT.ofs_getCoef(-l) - xFFT.ofs_getCoef(l)) << endl;
        curentStream.close();
    }


    gsl_vector_complex_free(data_complex);
    gsl_vector_free(xGSL);
    gsl_fft_real_wavetable_free(wavetable);
    gsl_fft_real_workspace_free(workspace);
}



/**
 *  \brief FFT and Storage in txt files of OFS objects alpha_i in qbtbp_ofs. Makes use of FFT routines via qbtbp_ofs_fft_unpack.
 *         The following coefficients are computed and stored:
 *              - The alpha (1 to 14).
 *              - The position of the primaries in EM coordinates system (redundant with alpha7-12): Xe, Ye, Ze
 *              - The position of the primaries in NC coordinates system: xe, ye, ze
 */
void qbtbp_ofs_fft_alpha( Ofts< Ofsd > &zt,     //zt = normalized Earth-Moon motion
                          Ofts< Ofsd > &Zt,     //Zt = normalized Sun-(Earth+Moon) motion
                          int nf,               //order of the Fourier expansions
                          QBCP_L& qbcp_l)       //QBCP
{
    //--------------------------
    //Physical params specific to the QBCP_L
    //--------------------------
    double c1    = qbcp_l.cs_em.c1;
    double gamma = qbcp_l.cs_em.gamma;

    //--------------------------
    //Mass ratio
    //--------------------------
    double mu_EM = qbcp_l.us_em.mu_EM;

    //--------------------------
    //Physical params in EM units
    //--------------------------
    double ns = qbcp_l.us_em.ns;  //Sun-(Earth+Moon) mean angular motion
    double ni = qbcp_l.us_em.ni;  //Earth-Moon mean angular motion
    double n  = qbcp_l.us_em.n;   //n = ni - ns
    double as = qbcp_l.us_em.as;  //Sun-(Earth+Moon) mean distance
    double ai = qbcp_l.us_em.ai;  //Earth-Moon mean distance
    double ms = qbcp_l.us_em.ms;  //Sun mass

    //--------------------------
    //FFT params
    //--------------------------
    int N      = pow(2,12);
    double tf  = 2*M_PI/n;
    double t   = 0;
    double eps = 1.0/as;

    //--------------------------
    //QBTBP init
    //--------------------------
    Ofsd bj(nf);
    Ofsd cj(nf);
    Ofsc ztc(nf);
    Ofsc Ztc(nf);

    //--------------------------
    //Building z(t) and Z(t)
    //--------------------------
    //From ots to ofs
    fts2fs(&bj, zt, eps);
    fts2fs(&cj, Zt, eps);

    //From double to cdouble in ztc and Ztc
    doubleToComplex(bj, ztc);
    doubleToComplex(cj, Ztc);

    //Derivatives
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t)
    cdouble zi;
    cdouble Zi;
    cdouble zidot;
    cdouble ziddot;
    cdouble Ziddot;

    //--------------------------
    // GSL objects
    //--------------------------
    gsl_vector *gsl_alpha1  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha2  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha3  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha4  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha5  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha6  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha7  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha8  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha9  = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha10 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha11 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha12 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha13 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha14 = gsl_vector_calloc(N);
    //rq: alpha15 is let equal to zero because of ERTBP vector field
    gsl_vector *gsl_alpha16 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha17 = gsl_vector_calloc(N);
    gsl_vector *gsl_alpha18 = gsl_vector_calloc(N);

    //Redundancy for the positions of the primaries
    gsl_vector *gsl_Xe = gsl_vector_calloc(N);
    gsl_vector *gsl_Ye = gsl_vector_calloc(N);
    gsl_vector *gsl_Ze = gsl_vector_calloc(N);
    gsl_vector *gsl_Xm = gsl_vector_calloc(N);
    gsl_vector *gsl_Ym = gsl_vector_calloc(N);
    gsl_vector *gsl_Zm = gsl_vector_calloc(N);
    gsl_vector *gsl_Xs = gsl_vector_calloc(N);
    gsl_vector *gsl_Ys = gsl_vector_calloc(N);
    gsl_vector *gsl_Zs = gsl_vector_calloc(N);


    //Positions of the primaries in NC coordinates
    gsl_vector *gsl_xe = gsl_vector_calloc(N);
    gsl_vector *gsl_ye = gsl_vector_calloc(N);
    gsl_vector *gsl_ze = gsl_vector_calloc(N);
    gsl_vector *gsl_xm = gsl_vector_calloc(N);
    gsl_vector *gsl_ym = gsl_vector_calloc(N);
    gsl_vector *gsl_zm = gsl_vector_calloc(N);
    gsl_vector *gsl_xs = gsl_vector_calloc(N);
    gsl_vector *gsl_ys = gsl_vector_calloc(N);
    gsl_vector *gsl_zs = gsl_vector_calloc(N);

    //--------------------------
    //Inner double variables
    //--------------------------
    double rv2;
    double alpha1i;
    double alpha2i;
    double alpha3i;
    double alpha4i;
    double alpha5i;
    double alpha6i;
    double alpha7i;
    double alpha8i;
    double alpha9i;
    double alpha10i;
    double alpha11i;
    double alpha12i;
    double alpha13i;
    double alpha14i;
    //rq: alpha15 is let equal to zero because of ERTBP vector field
    double alpha16i;
    double alpha17i;
    double alpha18i;

    double Ze[3], Zm[3], Zs[3];
    double ze[3], zm[3], zs[3];

    //Temp
    double alphati;

    //Derivatives
    double alpha1doti;
    double alpha2doti;
    double alpha3doti;

    //--------------------------
    //Time loop to prepare the FFT
    //--------------------------
    for(int i = 0 ; i < N; i++)
    {
        //---------------
        //update the time
        //---------------
        t = tf*i/N;

        //---------------
        //z(t), zdot(t), Z(t) and Zdot(t)
        //---------------
        zi     =  evz(ztc, t, n, ni, ai);
        Zi     =  evz(Ztc, t, n, ns, as);
        zidot  =  evzdot(ztc, ztcdot, t, n, ni, ai);
        ziddot =  evzddot(ztc, ztcdot, ztcddot, t, n, ni, ai);
        Ziddot =  evzddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
        //rv2= 1/(z*zb)
        rv2    =  creal(1.0/(zi*conj(zi)));

        //---------------
        // The alphas
        //---------------
        alpha1i  = +rv2;
        alpha2i  = -rv2*creal(zidot*conj(zi));
        alpha3i  = +rv2*cimag(zidot*conj(zi));
        alpha4i  = -ms/(1.0+ms)*creal(Ziddot*conj(zi));
        alpha5i  = -ms/(1.0+ms)*cimag(Ziddot*conj(zi));
        alpha6i  = +creal(cpow(zi*conj(zi), -1.0/2+0.0*I));
        //Sun
        alpha7i  = +rv2*creal(Zi*conj(zi));
        alpha8i  = +rv2*cimag(Zi*conj(zi));
        //Earth
        alpha9i  = +mu_EM;
        alpha10i = +0.0;
        //Moon
        alpha11i = +mu_EM-1.0;
        alpha12i = +0.0;
        //Temp
        alphati = +creal(zi*conj(zi));

        //Derivatives of alpha1,2,3
        alpha1doti = creal(+0.0*I-( zidot*conj(zi)+zi*conj(zidot) ) / cpow(zi*conj(zi), 2.0+0.0*I));
        alpha2doti = creal(+0.0*I-( creal(ziddot*conj(zi) + zidot*conj(zidot))*zi*conj(zi)
                                    - (zidot*conj(zi)+zi*conj(zidot))*creal(zidot*conj(zi)) ) / cpow(zi*conj(zi), 2.0+0.0*I));
        alpha3doti = creal(+0.0*I+( cimag(ziddot*conj(zi) + zidot*conj(zidot))*zi*conj(zi)
                                    - (zidot*conj(zi)+zi*conj(zidot))*cimag(zidot*conj(zi)) ) / cpow(zi*conj(zi), 2.0+0.0*I));

        //Alpha13 and 14, for NC computation
        alpha13i = -(alphati*alphati*(alpha2doti*alpha1i - alpha1doti*alpha2i)+ alphati*(alpha2i*alpha2i+alpha3i*alpha3i))*c1 + alpha4i/gamma;
        alpha14i =   alphati*alphati*(alpha3doti*alpha1i - alpha1doti*alpha3i)*c1 + alpha5i/gamma;

        //Alpha16i
        alpha16i = alpha2doti - alpha2i*alpha1doti/alpha1i + alpha2i*alpha2i + alpha3i*alpha3i;
        //Alpha17i
        alpha17i = alpha3doti - alpha3i*alpha1doti/alpha1i;
        //Alpha18i
        alpha18i = alpha1doti/alpha1i;

        //---------------
        // The primaries, again
        //---------------
        //Sun, EM coordinates
        Zs[0] = alpha7i;
        Zs[1] = alpha8i;
        Zs[2] = 0.0;
        //Sun, NC coordinates
        SYStoNC_prim(Zs, zs, c1, gamma);
        //Earth, EM coordinates
        Ze[0] = alpha9i;
        Ze[1] = alpha10i;
        Ze[2] = 0.0;
        //Earth, NC coordinates
        SYStoNC_prim(Ze, ze, c1, gamma);
        //Moon, EM coordinates
        Zm[0] = alpha11i;
        Zm[1] = alpha12i;
        Zm[2] = 0.0;
        //Moon, NC coordinates
        SYStoNC_prim(Zm, zm, c1, gamma);

        //---------------
        //Storage in GSL objects at each time
        //---------------
        gsl_vector_set(gsl_alpha1,  i, alpha1i);
        gsl_vector_set(gsl_alpha2,  i, alpha2i);
        gsl_vector_set(gsl_alpha3,  i, alpha3i);
        gsl_vector_set(gsl_alpha4,  i, alpha4i);
        gsl_vector_set(gsl_alpha5,  i, alpha5i);
        gsl_vector_set(gsl_alpha6,  i, alpha6i);
        gsl_vector_set(gsl_alpha7,  i, alpha7i);
        gsl_vector_set(gsl_alpha8,  i, alpha8i);
        gsl_vector_set(gsl_alpha9,  i, alpha9i);
        gsl_vector_set(gsl_alpha10, i, alpha10i);
        gsl_vector_set(gsl_alpha11, i, alpha11i);
        gsl_vector_set(gsl_alpha12, i, alpha12i);
        gsl_vector_set(gsl_alpha13, i, alpha13i);
        gsl_vector_set(gsl_alpha14, i, alpha14i);

        gsl_vector_set(gsl_alpha16, i, alpha16i);
        gsl_vector_set(gsl_alpha17, i, alpha17i);
        gsl_vector_set(gsl_alpha18, i, alpha18i);

        //---------------
        //Primary, EM coordinates
        //---------------
        gsl_vector_set(gsl_Xs,  i, Zs[0]);
        gsl_vector_set(gsl_Ys,  i, Zs[1]);
        gsl_vector_set(gsl_Zs,  i, Zs[2]);

        gsl_vector_set(gsl_Xe,  i, Ze[0]);
        gsl_vector_set(gsl_Ye,  i, Ze[1]);
        gsl_vector_set(gsl_Ze,  i, Ze[2]);

        gsl_vector_set(gsl_Xm,  i, Zm[0]);
        gsl_vector_set(gsl_Ym,  i, Zm[1]);
        gsl_vector_set(gsl_Zm,  i, Zm[2]);

        //---------------
        //Primary, NC coordinates
        //---------------
        gsl_vector_set(gsl_xs,  i, zs[0]);
        gsl_vector_set(gsl_ys,  i, zs[1]);
        gsl_vector_set(gsl_zs,  i, zs[2]);

        gsl_vector_set(gsl_xe,  i, ze[0]);
        gsl_vector_set(gsl_ye,  i, ze[1]);
        gsl_vector_set(gsl_ze,  i, ze[2]);

        gsl_vector_set(gsl_xm,  i, zm[0]);
        gsl_vector_set(gsl_ym,  i, zm[1]);
        gsl_vector_set(gsl_zm,  i, zm[2]);

    }

    //--------------------------
    //Unpack the FFT and put in data file
    //--------------------------
    qbtbp_ofs_fft_unpack(gsl_alpha1, qbcp_l.cs_em.F_COEF+"alpha1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha2, qbcp_l.cs_em.F_COEF+"alpha2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha3, qbcp_l.cs_em.F_COEF+"alpha3", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha4, qbcp_l.cs_em.F_COEF+"alpha4", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha5, qbcp_l.cs_em.F_COEF+"alpha5", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha6, qbcp_l.cs_em.F_COEF+"alpha6", nf, N, 1);
    //Sun
    qbtbp_ofs_fft_unpack(gsl_alpha7, qbcp_l.cs_em.F_COEF+"alpha7", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha8, qbcp_l.cs_em.F_COEF+"alpha8", nf, N, 0);
    //Earth
    qbtbp_ofs_fft_unpack(gsl_alpha9,  qbcp_l.cs_em.F_COEF+"alpha9",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha10, qbcp_l.cs_em.F_COEF+"alpha10", nf, N, 0);
    //Moon
    qbtbp_ofs_fft_unpack(gsl_alpha11, qbcp_l.cs_em.F_COEF+"alpha11", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha12, qbcp_l.cs_em.F_COEF+"alpha12", nf, N, 0);
    //NC additional coef
    qbtbp_ofs_fft_unpack(gsl_alpha13, qbcp_l.cs_em.F_COEF+"alpha13",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha14, qbcp_l.cs_em.F_COEF+"alpha14",  nf, N, 0);
    //VF with velocity additional coef
    qbtbp_ofs_fft_unpack(gsl_alpha16, qbcp_l.cs_em.F_COEF+"alpha16",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_alpha17, qbcp_l.cs_em.F_COEF+"alpha17",  nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_alpha18, qbcp_l.cs_em.F_COEF+"alpha18",  nf, N, 0);


    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_Zc without much difference
    //---------------
    qbtbp_ofs_fft_unpack(gsl_Xs, qbcp_l.cs_em.F_COEF+"Ps1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ys, qbcp_l.cs_em.F_COEF+"Ps2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zs, qbcp_l.cs_em.F_COEF+"Ps3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xe, qbcp_l.cs_em.F_COEF+"Pe1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ye, qbcp_l.cs_em.F_COEF+"Pe2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Ze, qbcp_l.cs_em.F_COEF+"Pe3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xm, qbcp_l.cs_em.F_COEF+"Pm1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ym, qbcp_l.cs_em.F_COEF+"Pm2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zm, qbcp_l.cs_em.F_COEF+"Pm3", nf, N, 1);


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_zc without much difference
    //---------------
    qbtbp_ofs_fft_unpack(gsl_xs, qbcp_l.cs_em.F_COEF+"ps1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ys, qbcp_l.cs_em.F_COEF+"ps2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zs, qbcp_l.cs_em.F_COEF+"ps3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xe, qbcp_l.cs_em.F_COEF+"pe1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ye, qbcp_l.cs_em.F_COEF+"pe2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_ze, qbcp_l.cs_em.F_COEF+"pe3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xm, qbcp_l.cs_em.F_COEF+"pm1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ym, qbcp_l.cs_em.F_COEF+"pm2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zm, qbcp_l.cs_em.F_COEF+"pm3", nf, N, 1);


    //--------------------------
    //Umemory release
    //--------------------------
    gsl_vector_free(gsl_alpha1);
    gsl_vector_free(gsl_alpha2);
    gsl_vector_free(gsl_alpha3);
    gsl_vector_free(gsl_alpha4);
    gsl_vector_free(gsl_alpha5);
    gsl_vector_free(gsl_alpha6);
    gsl_vector_free(gsl_alpha7);
    gsl_vector_free(gsl_alpha8);
    gsl_vector_free(gsl_alpha9);
    gsl_vector_free(gsl_alpha10);
    gsl_vector_free(gsl_alpha11);
    gsl_vector_free(gsl_alpha12);
    gsl_vector_free(gsl_alpha13);
    gsl_vector_free(gsl_alpha14);
    gsl_vector_free(gsl_alpha16);
    gsl_vector_free(gsl_alpha17);
    gsl_vector_free(gsl_alpha18);
    gsl_vector_free(gsl_Xe);
    gsl_vector_free(gsl_Ye);
    gsl_vector_free(gsl_Ze);
    gsl_vector_free(gsl_Xm);
    gsl_vector_free(gsl_Ym);
    gsl_vector_free(gsl_Zm);
    gsl_vector_free(gsl_Xs);
    gsl_vector_free(gsl_Ys);
    gsl_vector_free(gsl_Zs);
    gsl_vector_free(gsl_xe);
    gsl_vector_free(gsl_ye);
    gsl_vector_free(gsl_ze);
    gsl_vector_free(gsl_xm);
    gsl_vector_free(gsl_ym);
    gsl_vector_free(gsl_zm);
    gsl_vector_free(gsl_xs);
    gsl_vector_free(gsl_ys);
    gsl_vector_free(gsl_zs);
}


/**
 *  \brief FFT and Storage in txt files of OFS objects delta_i in qbtbp_ofs. Makes use of FFT routines via qbtbp_ofs_fft_unpack.
 *         The following coefficients are computed and stored:
 *              - The deltas (1 to 14).
 *              - The position of the primaries in SE coordinates system (redundant with alpha7-12): Xe, Ye, Ze
 *              - The position of the primaries in SENC coordinates system: xe, ye, ze
 */
void qbtbp_ofs_fft_delta(Ofts< Ofsd > &zt,     //zt = normalized Earth-Moon motion
                         Ofts< Ofsd > &Zt,     //Zt = normalized Sun-(Earth+Moon) motion
                         int nf,               //order of the Fourier expansions
                         QBCP_L &qbcp_l)       //QBCP
{
    //--------------------------
    //Physical params in EM units
    //--------------------------
    double as_EM  = qbcp_l.us_em.as;     //Sun-(Earth+Moon) mean distance
    double ms_EM  = qbcp_l.us_em.ms;     //Sun mass
    double mu_EM  = qbcp_l.us_em.mu_EM;  //Mass ratio
    double mu_SEM = qbcp_l.us_em.mu_SEM; //Mass ratio

    //--------------------------
    //Physical params in SE units
    //--------------------------
    double ns = qbcp_l.us_sem.ns;    //Sun-(Earth+Moon) mean angular motion
    double ni = qbcp_l.us_sem.ni;    //Earth-Moon mean angular motion
    double n  = qbcp_l.us_sem.n;     //n = ni - ns
    double as = qbcp_l.us_sem.as;    //Sun-(Earth+Moon) mean distance
    double ai = qbcp_l.us_sem.ai;    //Earth-Moon mean distance

    //--------------------------
    //Physical params specific to the QBCP_L
    //--------------------------
    double c1    = qbcp_l.cs_sem.c1;
    double gamma = qbcp_l.cs_sem.gamma;

    //--------------------------
    //FFT params
    //--------------------------
    int N      = pow(2,12);
    double tf  = 2*M_PI/n;
    double t   = 0;
    double eps = 1.0/as_EM;

    //--------------------------
    // Building z(t) and Z(t)
    //--------------------------
    Ofsd bj(nf);
    Ofsd cj(nf);
    Ofsc ztc(nf);
    Ofsc Ztc(nf);

    //From ots to ofs
    fts2fs(&bj, zt, eps);
    fts2fs(&cj, Zt, eps);

    //From double to cdouble in ztc and Ztc
    doubleToComplex(bj, ztc);
    doubleToComplex(cj, Ztc);

    //Derivatives.
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives.
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t) and derivatives
    cdouble zi;
    cdouble Zi;
    cdouble Zidot;
    cdouble Ziddot;

    //--------------------------
    // GSL objects
    //--------------------------
    gsl_vector *gsl_delta1  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta2  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta3  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta4  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta5  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta6  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta7  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta8  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta9  = gsl_vector_calloc(N);
    gsl_vector *gsl_delta10 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta11 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta12 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta13 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta14 = gsl_vector_calloc(N);
    //rq: delta15 is let equal to zero because of BCP vector field
    gsl_vector *gsl_delta16 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta17 = gsl_vector_calloc(N);
    gsl_vector *gsl_delta18 = gsl_vector_calloc(N);


    //Redundancy for the positions of the primaries
    gsl_vector *gsl_Xe = gsl_vector_calloc(N);
    gsl_vector *gsl_Ye = gsl_vector_calloc(N);
    gsl_vector *gsl_Ze = gsl_vector_calloc(N);
    gsl_vector *gsl_Xm = gsl_vector_calloc(N);
    gsl_vector *gsl_Ym = gsl_vector_calloc(N);
    gsl_vector *gsl_Zm = gsl_vector_calloc(N);
    gsl_vector *gsl_Xs = gsl_vector_calloc(N);
    gsl_vector *gsl_Ys = gsl_vector_calloc(N);
    gsl_vector *gsl_Zs = gsl_vector_calloc(N);


    //Positions of the primaries in NC coordinates
    gsl_vector *gsl_xe = gsl_vector_calloc(N);
    gsl_vector *gsl_ye = gsl_vector_calloc(N);
    gsl_vector *gsl_ze = gsl_vector_calloc(N);
    gsl_vector *gsl_xm = gsl_vector_calloc(N);
    gsl_vector *gsl_ym = gsl_vector_calloc(N);
    gsl_vector *gsl_zm = gsl_vector_calloc(N);
    gsl_vector *gsl_xs = gsl_vector_calloc(N);
    gsl_vector *gsl_ys = gsl_vector_calloc(N);
    gsl_vector *gsl_zs = gsl_vector_calloc(N);


    //--------------------------
    //Inner double variables
    //--------------------------
    double delta1i;
    double delta2i;
    double delta3i;
    double delta4i;
    double delta5i;
    double delta6i;
    double delta7i;
    double delta8i;
    double delta9i;
    double delta10i;
    double delta11i;
    double delta12i;
    double delta13i;
    double delta14i;
    //rq: delta15 is let equal to zero because of BCP vector field
    double delta16i;
    double delta17i;
    double delta18i;

    double Ze[3], Zm[3], Zs[3];
    double ze[3], zm[3], zs[3];

    //Temp
    double deltati;

    //Derivatives
    double delta1doti;
    double delta2doti;
    double delta3doti;


    //Temp variables
    //double g1;
    //double g2;
    double Rv2;

    //Time loop to prepare the FFT
    for(int i = 0 ; i < N; i++)
    {
        //---------------
        //update time
        //---------------
        t = tf*i/N;

        //---------------
        //z(t),  Z(t) and derivatives
        //---------------
        zi     =  evz(ztc, t, n, ni, ai);
        Zi     =  evz(Ztc, t, n, ns, as);
        Zidot  =  evzdot(Ztc, Ztcdot, t, n, ns, as);
        Ziddot =  evzddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
        //rv2= 1/(Z*Zb)
        Rv2 = creal(1.0/(Zi*conj(Zi)));

        //---------------
        // The deltas
        //---------------
        delta1i  = +Rv2;
        delta2i  = -Rv2*creal(Zidot * conj(Zi));
        delta3i  = +Rv2*cimag(Zidot * conj(Zi));
        delta4i  = +0.0;
        delta5i  = +0.0;
        delta6i  = sqrt(Rv2);
        //Sun
        delta7i  = +mu_SEM;
        delta8i  = +0.0;
        //Earth
        delta9i  = (mu_EM)*Rv2*creal(zi*conj(Zi)) - ms_EM/(1.0+ms_EM);
        delta10i = (mu_EM)*Rv2*cimag(zi*conj(Zi));
        //Moon
        delta11i = (mu_EM-1)*Rv2*creal(zi*conj(Zi)) - ms_EM/(1.0+ms_EM);
        delta12i = (mu_EM-1)*Rv2*cimag(zi*conj(Zi));


        //Temp
        deltati = +creal(Zi*conj(Zi));

        //Derivatives of delta1,2,3
        delta1doti = creal(+0.0*I-( Zidot*conj(Zi)+Zi*conj(Zidot) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));
        delta2doti = creal(+0.0*I-( creal(Ziddot*conj(Zi) + Zidot*conj(Zidot))*Zi*conj(Zi)
                                    - (Zidot*conj(Zi)+Zi*conj(Zidot))*creal(Zidot*conj(Zi)) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));
        delta3doti = creal(+0.0*I+( cimag(Ziddot*conj(Zi) + Zidot*conj(Zidot))*Zi*conj(Zi)
                                    - (Zidot*conj(Zi)+Zi*conj(Zidot))*cimag(Zidot*conj(Zi)) ) / cpow(Zi*conj(Zi), 2.0+0.0*I));

        //Delta13 and 14, for NC computation
        delta13i = -(deltati*deltati*(delta2doti*delta1i - delta1doti*delta2i)+ deltati*(delta2i*delta2i+delta3i*delta3i))*c1 + delta4i/gamma;
        delta14i =   deltati*deltati*(delta3doti*delta1i - delta1doti*delta3i)*c1 + delta5i/gamma;

        //Delta16i
        delta16i = delta2doti - delta2i*delta1doti/delta1i + delta2i*delta2i + delta3i*delta3i;
        //Delta17i
        delta17i = delta3doti - delta3i*delta1doti/delta1i;
        //Delta18i
        delta18i = delta1doti/delta1i;

        //---------------
        // The primaries, again
        //---------------
        //Sun, EM coordinates
        Zs[0] = delta7i;
        Zs[1] = delta8i;
        Zs[2] = 0.0;
        //Sun, NC coordinates
        SYStoNC_prim(Zs, zs, c1, gamma);
        //Earth, EM coordinates
        Ze[0] = delta9i;
        Ze[1] = delta10i;
        Ze[2] = 0.0;
        //Earth, NC coordinates
        SYStoNC_prim(Ze, ze, c1, gamma);
        //Moon, EM coordinates
        Zm[0] = delta11i;
        Zm[1] = delta12i;
        Zm[2] = 0.0;
        //Moon, NC coordinates
        SYStoNC_prim(Zm, zm, c1, gamma);

        //---------------
        //Storage in GSL objects at each time
        //---------------
        gsl_vector_set(gsl_delta1,  i, delta1i);
        gsl_vector_set(gsl_delta2,  i, delta2i);
        gsl_vector_set(gsl_delta3,  i, delta3i);
        gsl_vector_set(gsl_delta4,  i, delta4i);
        gsl_vector_set(gsl_delta5,  i, delta5i);
        gsl_vector_set(gsl_delta6,  i, delta6i);
        gsl_vector_set(gsl_delta7,  i, delta7i);
        gsl_vector_set(gsl_delta8,  i, delta8i);
        gsl_vector_set(gsl_delta9,  i, delta9i);
        gsl_vector_set(gsl_delta10, i, delta10i);
        gsl_vector_set(gsl_delta11, i, delta11i);
        gsl_vector_set(gsl_delta12, i, delta12i);
        gsl_vector_set(gsl_delta13, i, delta13i);
        gsl_vector_set(gsl_delta14, i, delta14i);

        gsl_vector_set(gsl_delta16, i, delta16i);
        gsl_vector_set(gsl_delta17, i, delta17i);
        gsl_vector_set(gsl_delta18, i, delta18i);

        //---------------
        //Primary, EM coordinates
        //---------------
        gsl_vector_set(gsl_Xs,  i, Zs[0]);
        gsl_vector_set(gsl_Ys,  i, Zs[1]);
        gsl_vector_set(gsl_Zs,  i, Zs[2]);

        gsl_vector_set(gsl_Xe,  i, Ze[0]);
        gsl_vector_set(gsl_Ye,  i, Ze[1]);
        gsl_vector_set(gsl_Ze,  i, Ze[2]);

        gsl_vector_set(gsl_Xm,  i, Zm[0]);
        gsl_vector_set(gsl_Ym,  i, Zm[1]);
        gsl_vector_set(gsl_Zm,  i, Zm[2]);

        //---------------
        //Primary, NC coordinates
        //---------------
        gsl_vector_set(gsl_xs,  i, zs[0]);
        gsl_vector_set(gsl_ys,  i, zs[1]);
        gsl_vector_set(gsl_zs,  i, zs[2]);

        gsl_vector_set(gsl_xe,  i, ze[0]);
        gsl_vector_set(gsl_ye,  i, ze[1]);
        gsl_vector_set(gsl_ze,  i, ze[2]);

        gsl_vector_set(gsl_xm,  i, zm[0]);
        gsl_vector_set(gsl_ym,  i, zm[1]);
        gsl_vector_set(gsl_zm,  i, zm[2]);

    }

    //--------------------------
    //Unpack the FFT and put in data file
    //--------------------------
    qbtbp_ofs_fft_unpack(gsl_delta1, qbcp_l.cs_sem.F_COEF+"alpha1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta2, qbcp_l.cs_sem.F_COEF+"alpha2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta3, qbcp_l.cs_sem.F_COEF+"alpha3", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta4, qbcp_l.cs_sem.F_COEF+"alpha4", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta5, qbcp_l.cs_sem.F_COEF+"alpha5", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta6, qbcp_l.cs_sem.F_COEF+"alpha6", nf, N, 1);
    //Sun
    qbtbp_ofs_fft_unpack(gsl_delta7, qbcp_l.cs_sem.F_COEF+"alpha7", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta8, qbcp_l.cs_sem.F_COEF+"alpha8", nf, N, 0);
    //Earth
    qbtbp_ofs_fft_unpack(gsl_delta9,  qbcp_l.cs_sem.F_COEF+"alpha9",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta10, qbcp_l.cs_sem.F_COEF+"alpha10", nf, N, 0);
    //Moon
    qbtbp_ofs_fft_unpack(gsl_delta11, qbcp_l.cs_sem.F_COEF+"alpha11", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta12, qbcp_l.cs_sem.F_COEF+"alpha12", nf, N, 0);
    //NC additional coef
    qbtbp_ofs_fft_unpack(gsl_delta13, qbcp_l.cs_sem.F_COEF+"alpha13",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta14, qbcp_l.cs_sem.F_COEF+"alpha14",  nf, N, 0);
    //VF with velocity additionnal coef
    qbtbp_ofs_fft_unpack(gsl_delta16, qbcp_l.cs_sem.F_COEF+"alpha16",  nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_delta17, qbcp_l.cs_sem.F_COEF+"alpha17",  nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_delta18, qbcp_l.cs_sem.F_COEF+"alpha18",  nf, N, 0);

    //---------------
    //Primary, EM coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_Zc without much difference
    //---------------
    qbtbp_ofs_fft_unpack(gsl_Xs, qbcp_l.cs_sem.F_COEF+"Ps1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ys, qbcp_l.cs_sem.F_COEF+"Ps2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zs, qbcp_l.cs_sem.F_COEF+"Ps3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xe, qbcp_l.cs_sem.F_COEF+"Pe1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ye, qbcp_l.cs_sem.F_COEF+"Pe2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Ze, qbcp_l.cs_sem.F_COEF+"Pe3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_Xm, qbcp_l.cs_sem.F_COEF+"Pm1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_Ym, qbcp_l.cs_sem.F_COEF+"Pm2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_Zm, qbcp_l.cs_sem.F_COEF+"Pm3", nf, N, 1);


    //---------------
    //Primary, NC coordinates
    //Note that, at this step, the vertical motion of the primaries is undefined,
    //so we can put either Even or Odd in the qbtbp_ofs_fft_unpack option of gsl_zc without much difference
    //---------------
    qbtbp_ofs_fft_unpack(gsl_xs, qbcp_l.cs_sem.F_COEF+"ps1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ys, qbcp_l.cs_sem.F_COEF+"ps2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zs, qbcp_l.cs_sem.F_COEF+"ps3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xe, qbcp_l.cs_sem.F_COEF+"pe1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ye, qbcp_l.cs_sem.F_COEF+"pe2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_ze, qbcp_l.cs_sem.F_COEF+"pe3", nf, N, 1);

    qbtbp_ofs_fft_unpack(gsl_xm, qbcp_l.cs_sem.F_COEF+"pm1", nf, N, 1);
    qbtbp_ofs_fft_unpack(gsl_ym, qbcp_l.cs_sem.F_COEF+"pm2", nf, N, 0);
    qbtbp_ofs_fft_unpack(gsl_zm, qbcp_l.cs_sem.F_COEF+"pm3", nf, N, 1);


    //--------------------------
    //Memory release
    //--------------------------
    gsl_vector_free(gsl_delta1);
    gsl_vector_free(gsl_delta2);
    gsl_vector_free(gsl_delta3);
    gsl_vector_free(gsl_delta4);
    gsl_vector_free(gsl_delta5);
    gsl_vector_free(gsl_delta6);
    gsl_vector_free(gsl_delta7);
    gsl_vector_free(gsl_delta8);
    gsl_vector_free(gsl_delta9);
    gsl_vector_free(gsl_delta10);
    gsl_vector_free(gsl_delta11);
    gsl_vector_free(gsl_delta12);
    gsl_vector_free(gsl_delta13);
    gsl_vector_free(gsl_delta14);
    gsl_vector_free(gsl_delta16);
    gsl_vector_free(gsl_delta17);
    gsl_vector_free(gsl_delta18);
    gsl_vector_free(gsl_Xe);
    gsl_vector_free(gsl_Ye);
    gsl_vector_free(gsl_Ze);
    gsl_vector_free(gsl_Xm);
    gsl_vector_free(gsl_Ym);
    gsl_vector_free(gsl_Zm);
    gsl_vector_free(gsl_Xs);
    gsl_vector_free(gsl_Ys);
    gsl_vector_free(gsl_Zs);
    gsl_vector_free(gsl_xe);
    gsl_vector_free(gsl_ye);
    gsl_vector_free(gsl_ze);
    gsl_vector_free(gsl_xm);
    gsl_vector_free(gsl_ym);
    gsl_vector_free(gsl_zm);
    gsl_vector_free(gsl_xs);
    gsl_vector_free(gsl_ys);
    gsl_vector_free(gsl_zs);
}

//-----------------------------------------------------------------------------
// Reccurence for computing the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief One step of the recurrence scheme in qbtbp_ofs. See comments in src/h for the signification of the parameters.
 */
void qbtbp_ofs_recurrence(Ofts< Ofsd > &zr_ofts, //zr_ofts = normalized Earth-Moon motion
                          Ofts< Ofsd > &Zr_ofts, //Zr_ofts = normalized Sun-(Earth+Moon) motion
                          Ofts< Ofsd > &z1,      //z1 = \bar{zr_ofts}
                          Ofts< Ofsd > &z2,      //z2 = \bar{zr_ofts}^(-3/2)
                          Ofts< Ofsd > &z3,      //z3 = zr_ofts^(-1/2)
                          Ofts< Ofsd > &z4,      //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
                          Ofts< Ofsd > &z5,      //z5 = \bar{zr_ofts}^(-1/2)
                          Ofts< Ofsd > &z6,      //z6 = z3*z5 = (zr_ofts*\bar{zr_ofts})^(-1/2)
                          Ofts< Ofsd > &a1,      //a1 = mu*exp(itheta)*zr_ofts
                          Ofts< Ofsd > &a2,      //a2 = epsilon*a1
                          Ofts< Ofsd > &a3,      //a3 = Zr_ofts-mu*exp(itheta)*epsilon*zr_ofts
                          Ofts< Ofsd > &a4,      //a4 = a3^(-1/2).
                          Ofts< Ofsd > &a5,      //a5 = \bar{a3}
                          Ofts< Ofsd > &a6,      //a6 = a5^(-3/2).
                          Ofts< Ofsd > &a7,      //a7 = a4*a6,
                          Ofts< Ofsd > &b1,      //b1 = (1-mu)*exp(ithetb)*zr_ofts
                          Ofts< Ofsd > &b2,      //b2 = epsilon*b1
                          Ofts< Ofsd > &b3,      //b3 = Zr_ofts+(1-mu)*exp(ithetb)*epsilon*zr_ofts
                          Ofts< Ofsd > &b4,      //b4 = b3^(-1/2)
                          Ofts< Ofsd > &b5,      //b5 = \bar{b3}
                          Ofts< Ofsd > &b6,      //b6 = b5^(-3/2)
                          Ofts< Ofsd > &b7,      //b7 = b4*b6,
                          Ofts< Ofsd > &Pm,      //Pm = -ms/as^2*exp(-itheta)*b7 + ms/as^2*exp(-itheta)*a7 - z4
                          Ofts< Ofsd > &Qm,      //Qm = -ns^2*mu*b7 - ns^2*(1-mu)*a7
                          Ofts< Ofsd > &epsilon, //epsilon = Ofts with just 1.0 at order 1
                          Ofsd  &Pfm,            //Pfm = Pm(j) at order n
                          Ofsd  &Qfm,            //Qfm = Qm(j) at order n
                          Ofsd  &ufm,            //ufm = zr_ofts(j) at order n
                          Ofsd  &vfm,            //vfm = Zr_ofts(j) at order n
                          Ofsd &sigma1,          //sigma1 = exp(-itheta)
                          Ofsd &sigma2,          //sigma2 = exp(+itheta)
                          int m,                 //order
                          int nf,                //order of the Fourier expansions
                          QBCP_L &qbcp_l)        //QBCP

{
    //Physical param
    double n  = qbcp_l.us_em.n;
    double mu = qbcp_l.us_em.mu_EM;
    double ns = qbcp_l.us_em.ns;
    double as = qbcp_l.us_em.as;
    double ms = qbcp_l.us_em.ms;

    if(m == 0)
    {
        //---------------------------------------------------------
        //Order 0 of the recurrence
        //---------------------------------------------------------
        z1.ccopy(zr_ofts, m);
        //z1 = \bar{zr_ofts}
        z1.conjugate(m);
        //z2 = \bar{zr_ofts}^(-3/2). At order 0, z1 = z2
        z2.ccopy(z1, m);
        //z3 = zr_ofts^(-1/2). At order 0, z3 = z
        z3.ccopy(zr_ofts, m);
        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);
        //z5 = \bar{zr_ofts}^(-1/2). At order 0, z5 = z1
        z5.ccopy(z1, m);
        //z6 = z3*z5
        z6.ofts_sprod(z3, z5, m);

        //a1 = mu*exp(itheta)*z
        a1.ofts_smult_tu(zr_ofts, sigma2, mu, m);
        //a2 = epsilon*a1
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z
        a3.ofts_smult_u(Zr_ofts, 1.0, m);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);
        //a4 = a3^(-1/2). At order 0, a4 = a3
        a4.ccopy(a3, m);

        //a5 = \bar{a3}
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2). At order 0, a6 = a5
        a6.ccopy(a5, m);
        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);


        //b1 = (1-mu)*exp(itheta)*z
        b1.ofts_smult_tu(zr_ofts, sigma2, (1.0-mu), m);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        b3.ofts_smult_u(Zr_ofts, 1.0, m);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);
        //b4 = b3^(-1/2). At order 0, b4 = b3
        b4.ccopy(b3, m);
        //b5 = \bar{b3}
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2). At order 0, b6 = b5
        b6.ccopy(b5, m);
        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);


        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, sigma1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);
    }
    else
    {
        //---------------------------------------------------------
        //Order m of the recurrence
        //---------------------------------------------------------
        //z1 = zr_ofts at order m-1
        z1.ccopy(zr_ofts, m-1);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m-1);

        //z2 = \bar{z}^^(-3/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) z2.ofts_pows(z1, -3.0/2, m-1);
        z2.ofts_pows(z1, -3.0/2, m);

        //z5 = \bar{z}^^(-1/2)
        //CAN BE USED BECAUSE z1[0] = 1.0 !!
        if(m>1) z5.ofts_pows(z1, -1.0/2, m-1);
        z5.ofts_pows(z1, -1.0/2, m);

        //z3 = z^(-1/2)
        //CAN BE USED BECAUSE z[0] = 1.0 !!
        if(m>1) z3.ofts_pows(zr_ofts, -1.0/2, m-1);
        z3.ofts_pows(zr_ofts, -1.0/2, m);

        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);

        //z6 = z3*z5
        z6.ofts_sprod(z3, z5, m);

        //a1 = mu*exp(itheta)*z at order m-1
        if(m>1) a1.ofts_smult_tu(zr_ofts, sigma2, mu, m-1);
        //a2 = epsilon*a1 at order m
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z at order m-1
        if(m>1) a3.ofts_smult_u(Zr_ofts, 1.0, m-1);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);

        //a4 = a3^(-1/2)
        //CAN BE USED BECAUSE a3[0] = 1.0 !!
        if(m>1) a4.ofts_pows(a3, -1.0/2, m-1);
        a4.ofts_pows(a3, -1.0/2, m);

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
        if(m>1) a6.ofts_pows(a5, -3.0/2, m-1);
        a6.ofts_pows(a5, -3.0/2, m);

        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);

        //b1 = (1-mu)*exp(itheta)*z
        if(m>1) b1.ofts_smult_tu(zr_ofts, sigma2, (1.0-mu), m-1);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        if(m>1) b3.ofts_smult_u(Zr_ofts, 1.0, m-1);
        //b3+= (1-mu)*exp(itheta)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);

        //b4 = b3^(-1/2)
        //CAN BE USED BECAUSE b3[0] = 1.0 !!
        if(m>1) b4.ofts_pows(b3, -1.0/2, m-1);
        b4.ofts_pows(b3, -1.0/2, m);

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
        if(m>1) b6.ofts_pows(b5, -3.0/2, m-1);
        b6.ofts_pows(b5, -3.0/2, m);

        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);

        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, sigma1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, sigma1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);

        //-------------------------
        // Solving equations
        //-------------------------
        Pfm.ccopy(Pm.getTerm(m)->getCoef(0));
        Qfm.ccopy(Qm.getTerm(m)->getCoef(0));


        //Order 0
        ufm.setCoef(-1.0/3*Pfm.ofs_getCoef(0), 0);          //u0 = -1/3*p0
        vfm.setCoef(-1.0/(3*ns*ns)*Qfm.ofs_getCoef(0), 0);  //v0 = -1/(3*ns^2)*p0

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
            uj = ( k2*Pfm.ofs_getCoef(-j) - k3*Pfm.ofs_getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, -j);
            //uj
            uj = (-k3*Pfm.ofs_getCoef(-j) + k1*Pfm.ofs_getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, j);
            //v-j
            vj = ( l2*Qfm.ofs_getCoef(-j) - l3*Qfm.ofs_getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, -j);
            //vj
            vj = (-l3*Qfm.ofs_getCoef(-j) + l1*Qfm.ofs_getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, j);

        }

        //Update z and Z
        zr_ofts.getCoef(m,0)->ccopy(ufm);
        Zr_ofts.getCoef(m,0)->ccopy(vfm);

    }

    if(m == zr_ofts.getOrder())
    {
        //Last update to account for final order in all series
        //---------------------------------
        //z1 = zr_ofts at order m-1
        z1.ccopy(zr_ofts, m);
        //z1 = \bar{z1} at order m-1
        z1.conjugate(m);
        //z2 = \bar{z}^^(-3/2)
        if(m>1) z2.ofts_pows(z1, -3.0/2, m);
        //z5 = \bar{z}^^(-1/2)
        if(m>1) z5.ofts_pows(z1, -1.0/2, m);
        //z3 = z^(-1/2)
        if(m>1) z3.ofts_pows(zr_ofts, -1.0/2, m);
    }

}

//-----------------------------------------------------------------------------
// Integrating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Derivatives of the QBTBP. To plug into GSL integrator.
 *
 * Note: the use of gsl livrary forces us to use double variables
 * As a consequence, in the case of z and Z, we need to use real and imaginary parts
 * as separate variables
 */
int qbtbp_derivatives(double t, const double y[], double f[], void *params)
{
    //Parameters for the qbtbp
    double mu = * (double * ) params;
    double ms = * ((double * ) (params)+1);

    //Reconstruction of z and Z
    cdouble z = y[0] + I*y[1];
    cdouble Z = y[2] + I*y[3];

    cdouble temp1 = Z-mu*z;
    temp1 = temp1/pow(cabs(temp1), 3.0);

    cdouble temp2 = Z+(1-mu)*z;
    temp2 = temp2/pow(cabs(temp2), 3.0);

    cdouble zdd = +0.0*I-z/pow(cabs(z), 3.0) + ms*(temp1-temp2);
    cdouble Zdd = -(1+ms)*(mu*temp2 + (1-mu)*temp1);

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

//-----------------------------------------------------------------------------
// Testing the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, Ofsc &bjc, Ofsc &cjc, OdeStruct ode_s, QBCP_L &qbcp_l)
{
    cout << "---------------------------------------------------" << endl;
    cout << "               Test of the                         " << endl;
    cout << "   Quasi-Bicircular Three-Body Problem (QBTBP)     " << endl;
    cout << "---------------------------------------------------" << endl;
    //Initialization
    double as = qbcp_l.us_em.as;
    double n  = qbcp_l.us_em.n;
    double ns = qbcp_l.us_em.ns;
    int nf    = bjc.getOrder();

    //z(0) and Z(0)
    cdouble z0 = bjc.evaluate(0.0);
    cdouble Z0 = as*cjc.evaluate(0.0);

    //zdot(0) and Zdot(0)
    Ofsc zdot(nf);
    for(int l = -nf; l<=nf; l++) zdot.setCoef(I*(1.0+l*n)*bjc.ofs_getCoef(l), l);
    cdouble zdot0 = zdot.evaluate(0.0);

    Ofsc Zdot(nf);
    for(int l = -nf; l<=nf; l++) Zdot.setCoef(I*(ns+l*n)*cjc.ofs_getCoef(l), l);
    cdouble Zdot0 = as*Zdot.evaluate(0.0);

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

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << "Initial positions: " << endl;
    cout << "Internal motion: z(t=0.0) = " << creal(z0) << " + " << cimag(z0) << "i" <<  endl;
    cout << "External motion: Z(t=0.0) = " << creal(Z0) << " + " << cimag(Z0) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;

    //Loop
    double h = Config::configManager().G_PREC_HSTART(); //first guess for stepper
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
    cout << "Final t: " << t << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Results from integration" << endl;
    cout << "Internal motion: z(t=t1) = " << yv[0] << " + " << yv[1] << "i" << endl;
    cout << "External motion: Z(t=t1) = " << yv[2] << " + " << yv[3] << "i" << endl;
    cout << "-------------------------------------------" << endl;

    //Analytical final state
    cdouble zfinal = (cos(t1)+I*sin(t1))*bjc.evaluate(n*t1);
    cdouble Zfinal = as*(cos(ns*t1)+I*sin(ns*t1))*cjc.evaluate(n*t1);

    cdouble zdotfinal = (cos(t1)+I*sin(t1))*zdot.evaluate(n*t1);
    cdouble Zdotfinal = as*(cos(ns*t1)+I*sin(ns*t1))*Zdot.evaluate(n*t1);


    cout << "Analytical results" << endl;
    cout << "Internal motion: z(t=t1) = " << creal(zfinal) << " + " << cimag(zfinal) << "i" <<  endl;
    cout << "External motion: Z(t=t1) = " << creal(Zfinal) << " + " << cimag(Zfinal) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;


    cout << "Absolute delta between analytical and numerical results: " << endl;
    cout << "dz    = " << cabs(zfinal - yv[0]-I*yv[1])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dzdot = " << cabs(zdotfinal - yv[4]-I*yv[5])/*/cabs(zdotfinal)*100 << " %" */<< endl;
    cout << "dZ    = " << cabs(Zfinal - yv[2]-I*yv[3])/*/cabs(zfinal)*100 << " %" */<< endl;
    cout << "dZdot = " << cabs(Zdotfinal - yv[6]-I*yv[7])/*/cabs(Zdotfinal)*100 << " %" */<< endl;
    cout << "-------------------------------------------" << endl;
}

/**
 *  \brief Test function used inside qbtbp_test_IN_EM_SEM.
 */
void qbpcomp(Ofsc &bjc,
             Ofsc &cjc,
             const double xe0_SEM[],
             const double xe0_EM[],
             const double ye0_IN[],
             double t0,
             double t0c)
{
    //State vectors
    double xem_SEM[6];
    double xem_EM[6];
    double xe0_IN_from_SEM[6];
    double xe0_IN_from_EM[6];
    double xe0_EM2[6];
    double xe0_SEM2[6];

    //-----------------------------
    // SEM -> EM -> IN
    //-----------------------------
    //Planet position & momenta in SEM ref and SEM units
    SEMvtoSEMm(t0c, xe0_SEM, xem_SEM, &SEML);
    //Planet position & momenta in EM ref and EM units
    SEMmtoEMm(t0c, xem_SEM, xem_EM, &SEML);
    //Planet position & velocity in EM ref and EM units
    EMmtoEMv(t0, xem_EM, xe0_EM2, &SEML);
    //Planet position & velocity in IN ref and EM units
    EMtoIN(t0, xe0_EM2, xe0_IN_from_SEM,  &SEML);

    //-----------------------------
    // EM -> IN
    //-----------------------------
    //Planet position & velocities in IN ref and EM units
    EMtoIN(t0, xe0_EM, xe0_IN_from_EM,  &SEML);
    cout << "-------------------------------------------" << endl;
    cout << "Error in initial planet position & velocity wrt IN" << endl;
    cout << "SEM -> EM -> IN         EM->IN  " << endl;
    for(int i = 0; i < 6; i++)
    {
        cout <<  xe0_IN_from_SEM[i]-ye0_IN[i] << "    " << xe0_IN_from_EM[i]-ye0_IN[i] << endl;
    }

    //-----------------------------
    // EM -> IN -> SE
    //-----------------------------
    //Planet position & momenta in EM ref and EM units
    EMvtoEMm(t0, xe0_EM, xem_EM, &SEML);
    //Planet position & momenta in SEM ref and SEM units
    EMmtoSEMm(t0, xem_EM, xem_SEM, &SEML);
    //Planet position & velocity in SEM ref and SEM units
    SEMmtoSEMv(t0c, xem_SEM, xe0_SEM2, &SEML);

    cout << "-------------------------------------------" << endl;
    cout << "SEM          EM -> IN -> SEM     Delta " << endl;
    for(int i = 0; i < 6; i++)
    {
        cout << xe0_SEM[i] << "    " <<  xe0_SEM2[i] << "    " <<  xe0_SEM[i]-xe0_SEM2[i] << endl;
    }
    cout << "-------------------------------------------" << endl;
}

/**
 *  \brief Test function to compare the analytical solutions found in IN, EM, and SEM coordinates
 */
void qbtbp_test_IN_EM_SEM(double t1, Ofsc &bjc, Ofsc &cjc)
{
    //Initialization
    USYS us_em  = SEML.us_em;
    USYS us_sem = SEML.us_sem;
    int nf      = SEML.nf;
    //Derivatives of z and Z
    Ofsc zdot(bjc);
    Ofsc Zdot(cjc);
    zdot.dot(us_em.n);
    Zdot.dot(us_em.n);

    cout << "---------------------------------------------------" << endl;
    cout << "               Test in SEM units                   " << endl;
    cout << "---------------------------------------------------" << endl;
    double t0 = t1;             //new t0 in EM units
    double t0c =t0*us_em.ns;    //new t0 in SEM units

    //z(t0) and Z(t0)
    cdouble z0 = evz(bjc, t0, us_em.n, us_em.ni, us_em.ai);
    cdouble Z0 = evz(cjc, t0, us_em.n, us_em.ns, us_em.as);
    cdouble zdot0 = evzdot(bjc, zdot, t0, us_em.n, us_em.ni, us_em.ai);
    cdouble Zdot0 = evzdot(cjc, Zdot, t0, us_em.n, us_em.ns, us_em.as);

    //-------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //-------------------------------------------------------------------------------
    double delta[15];
    evaluateCoef(delta, t0c, us_sem.n, nf, SEML.cs_sem.coeffs, 15);
    double deltad[15];
    evaluateCoefDerivatives(deltad, t0c, us_sem.n, nf, SEML.cs_sem.coeffs, 15);

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //-------------------------------------------------------------------------------
    double alpha[15];
    evaluateCoef(alpha, t0, us_em.n, nf, SEML.cs_em.coeffs, 15);
    double alphad[15];
    evaluateCoefDerivatives(alphad, t0, us_em.n, nf, SEML.cs_em.coeffs, 15);

    cout << "---------------------------------------------------" << endl;
    cout << "Earth  " << endl;
    //Earth position & momenta in SEM ref and SEM units
    double xe0_SEM[6];
    xe0_SEM[0] = delta[8];
    xe0_SEM[1] = delta[9];
    xe0_SEM[2] = 0.0;
    xe0_SEM[3] = deltad[8];
    xe0_SEM[4] = deltad[9];
    xe0_SEM[5] = 0.0;
    //Earth position & momenta in EM ref and EM units
    double xe0_EM[6];
    xe0_EM[0] = us_em.mu_EM;
    xe0_EM[1] = 0.0;
    xe0_EM[2] = 0.0;
    xe0_EM[3] = 0.0;
    xe0_EM[4] = 0.0;
    xe0_EM[5] = 0.0;
    //Earth position & momenta in IN ref and EM units
    double ye0_IN[6];
    ye0_IN[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) + us_em.mu_EM*creal(z0);
    ye0_IN[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) + us_em.mu_EM*cimag(z0);
    ye0_IN[2] = 0.0;
    ye0_IN[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) + us_em.mu_EM*creal(zdot0);
    ye0_IN[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) + us_em.mu_EM*cimag(zdot0);
    ye0_IN[5] = 0.0;
    //Comparison
    qbpcomp(bjc, cjc, xe0_SEM, xe0_EM, ye0_IN, t0, t0c);

    cout << "---------------------------------------------------" << endl;
    cout << "Moon  " << endl;
    //Moon position & velocity in SEM ref and SEM units
    double xm0_SEM[6];
    xm0_SEM[0] = delta[10];
    xm0_SEM[1] = delta[11];
    xm0_SEM[2] = 0.0;
    xm0_SEM[3] = deltad[10];
    xm0_SEM[4] = deltad[11];
    xm0_SEM[5] = 0.0;
    //Moon position & velocity in EM ref and EM units
    double xm0_EM[6];
    xm0_EM[0] = us_em.mu_EM-1;
    xm0_EM[1] = 0.0;
    xm0_EM[2] = 0.0;
    xm0_EM[3] = 0.0;
    xm0_EM[4] = 0.0;
    xm0_EM[5] = 0.0;
    //Moon position & momenta in IN ref and EM units
    double ym0_IN[6];
    ym0_IN[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) - (1-us_em.mu_EM)*creal(z0);
    ym0_IN[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) - (1-us_em.mu_EM)*cimag(z0);
    ym0_IN[2] = 0.0;
    ym0_IN[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) - (1-us_em.mu_EM)*creal(zdot0);
    ym0_IN[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) - (1-us_em.mu_EM)*cimag(zdot0);
    ym0_IN[5] = 0.0;
    //Comparison
    qbpcomp(bjc, cjc, xm0_SEM, xm0_EM, ym0_IN, t0, t0c);


    cout << "---------------------------------------------------" << endl;
    cout << "Sun  " << endl;
    //Sun position & velocity in SEM ref and SEM units
    double xs0_SEM[6];
    xs0_SEM[0] = us_sem.mu_SEM;
    xs0_SEM[1] = 0.0;
    xs0_SEM[2] = 0.0;
    xs0_SEM[3] = 0.0;
    xs0_SEM[4] = 0.0;
    xs0_SEM[5] = 0.0;
    //Sun position & velocity in EM ref and EM units
    double xs0_EM[6];
    xs0_EM[0] = alpha[6];
    xs0_EM[1] = alpha[7];
    xs0_EM[2] = 0.0;
    xs0_EM[3] = alphad[6];
    xs0_EM[4] = alphad[7];
    xs0_EM[5] = 0.0;
    //Sun position & momenta in IN ref and EM units
    double ys0_IN[6];
    ys0_IN[0] = 1.0/(1.0+us_em.ms)*creal(Z0);
    ys0_IN[1] = 1.0/(1.0+us_em.ms)*cimag(Z0);
    ys0_IN[2] = 0.0;
    ys0_IN[3] = 1.0/(1.0+us_em.ms)*creal(Zdot0);
    ys0_IN[4] = 1.0/(1.0+us_em.ms)*cimag(Zdot0);
    ys0_IN[5] = 0.0;
    //Comparison
    qbpcomp(bjc, cjc, xs0_SEM, xs0_EM, ys0_IN, t0, t0c);

}

/**
 *  \brief Comparison of the QBTBP computed via FFT or via OFS computations
 */
void qbtbp_test_FFT_vs_OFS(Ofsc &bjc,       //zt(t)
                           Ofsc &cjc,       //Zt(t)
                           int nf,          //order of the Fourier expansions
                           int N,           //Number of points
                           int type,        //Type of reference
                           OdeStruct ode_s, //ode structure
                           QBCP_L& qbcp_l)  //QBCP
{
    reset_ode_structure(&ode_s);
    //For print
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //Creation of the matrices of results
    double **alpha_INT_VS_OFS = dmatrix(0, 7, 0, N-1);
    double **alpha_INT_VS_FFT = dmatrix(0, 7, 0, N-1);
    double **alpha_INT        = dmatrix(0, 7, 0, N-1);

    //Retrieving from txt files
    double *alpha_OFS = dvector(0, 8*(nf+1)-1);
    double *alpha_FFT = dvector(0, 8*(nf+1)-1);

    cout << "Retrieving OFS coefficients " << endl;
    coefRetrieving(qbcp_l.cs.F_COEF+"alpha", alpha_OFS, nf, 0, 0, 8);
    cout << "Retrieving FFT coefficients " << endl;
    coefRetrieving(qbcp_l.cs.F_COEF+"alpha", alpha_FFT, nf, 0, 1, 8);

    //Physical params in EM units
    //--------------------------
    double ns = qbcp_l.us_em.ns;  //Sun-(Earth+Moon) mean angular motion
    double ni = qbcp_l.us_em.ni;  //Earth-Moon mean angular motion
    double n  = qbcp_l.us_em.n;   //n = ni - ns
    double as = qbcp_l.us_em.as;  //Sun-(Earth+Moon) mean distance
    double ai = qbcp_l.us_em.ai;  //Earth-Moon mean distance
    double ms = qbcp_l.us_em.ms;  //Sun mass

    //QBTBP init
    //--------------------------
    Ofsc ztc(bjc);
    Ofsc Ztc(cjc);
    //Derivatives
    Ofsc ztcdot(ztc);
    Ofsc Ztcdot(Ztc);
    //Derivation
    ztcdot.dot(n);
    Ztcdot.dot(n);

    //Double derivatives
    Ofsc ztcddot(ztcdot);
    Ofsc Ztcddot(Ztcdot);
    //Derivation
    ztcddot.dot(n);
    Ztcddot.dot(n);

    //z(t) and Z(t)
    cdouble zi;
    cdouble Zi;
    cdouble zidot;
    cdouble Ziddot;

    double rv2;
    double alpha1i;
    double alpha2i;
    double alpha3i;
    double alpha4i;
    double alpha5i;
    double alpha6i;
    double eval_OFS[8];
    double eval_FFT[8];

    //z(0) and Z(0)
    cdouble z0 = bjc.evaluate(0.0);
    cdouble Z0 = as*cjc.evaluate(0.0);

    //zdot(0) and Zdot(0)
    Ofsc zdot(nf);
    for(int l = -nf; l<=nf; l++) zdot.setCoef(I*(1.0+l*n)*bjc.ofs_getCoef(l), l);
    cdouble zdot0 = zdot.evaluate(0.0);

    Ofsc Zdot(nf);
    for(int l = -nf; l<=nf; l++) Zdot.setCoef(I*(ns+l*n)*cjc.ofs_getCoef(l), l);
    cdouble Zdot0 = as*Zdot.evaluate(0.0);

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

    //Loop
    double ti;
    double tf  = 2*M_PI/n;
    double t   = 0.0;
    int status;
    double f[8];
    for(int i = 0 ; i < N; i++)
    {
        //update
        ti = tf*i/N;

        //Computing the reference values
        if(type == 0)
        {
            //z(t), zdot(t), Z(t) and Zdot(t) if the type of reference is the numerical integration  (type = 0)
            //Evolve one step
            status = gsl_odeiv2_driver_apply (ode_s.d, &t, ti, yv);
            //Compute z(t) and Z(t)
            zi     =  yv[0] + yv[1]*I;
            Zi     =  yv[2] + yv[3]*I;
            zidot  =  yv[4] + yv[5]*I;
            qbtbp_derivatives(t, yv, f, ode_s.sys.params);
            Ziddot =  f[6] + f[7]*I;
            rv2    =  creal(1.0/(zi*conj(zi)));
        }
        else
        {
            //z(t), zdot(t), Z(t) and Zdot(t) if the type of reference is the Fourier expansion (type = 1)
            zi     =  evz(ztc, t, n, ni, ai);
            Zi     =  evz(Ztc, t, n, ns, as);
            zidot  =  evzdot(ztc, ztcdot, t, n, ni, ai);
            Ziddot =  evzddot(Ztc, Ztcdot, Ztcddot, t, n, ns, as);
            rv2    =  creal(1.0/(zi*conj(zi)));
        }


        //Evalutation of FFT and OFS results
        evaluateCoef(eval_OFS, t, n, nf, alpha_OFS, 8);
        evaluateCoef(eval_FFT, t, n, nf, alpha_FFT, 8);

        alpha1i  = rv2;
        alpha2i  = -rv2*creal(zidot*conj(zi));
        alpha3i  = +rv2*cimag(zidot*conj(zi));
        alpha4i  = -ms/(1.0+ms)*creal(Ziddot*conj(zi));
        alpha5i  = -ms/(1.0+ms)*cimag(Ziddot*conj(zi));
        alpha6i  = creal(cpow(zi*conj(zi), -1.0/2+0.0*I));

        //Results in matrices
        alpha_INT_VS_OFS[0][i] = fabs(alpha1i - eval_OFS[0]);
        alpha_INT_VS_OFS[1][i] = fabs(alpha2i - eval_OFS[1]);
        alpha_INT_VS_OFS[2][i] = fabs(alpha3i - eval_OFS[2]);
        alpha_INT_VS_OFS[3][i] = fabs(alpha4i - eval_OFS[3]);
        alpha_INT_VS_OFS[4][i] = fabs(alpha5i - eval_OFS[4]);
        alpha_INT_VS_OFS[5][i] = fabs(alpha6i - eval_OFS[5]);
        alpha_INT_VS_OFS[6][i] = fabs(rv2*creal(Zi*conj(zi)) - eval_OFS[6]);
        alpha_INT_VS_OFS[7][i] = fabs(rv2*cimag(Zi*conj(zi)) - eval_OFS[7]);

        alpha_INT_VS_FFT[0][i] = fabs(alpha1i - eval_FFT[0]);
        alpha_INT_VS_FFT[1][i] = fabs(alpha2i - eval_FFT[1]);
        alpha_INT_VS_FFT[2][i] = fabs(alpha3i - eval_FFT[2]);
        alpha_INT_VS_FFT[3][i] = fabs(alpha4i - eval_FFT[3]);
        alpha_INT_VS_FFT[4][i] = fabs(alpha5i - eval_FFT[4]);
        alpha_INT_VS_FFT[5][i] = fabs(alpha6i - eval_FFT[5]);
        alpha_INT_VS_FFT[6][i] = fabs(rv2*creal(Zi*conj(zi)) - eval_FFT[6]);
        alpha_INT_VS_FFT[7][i] = fabs(rv2*cimag(Zi*conj(zi)) - eval_FFT[7]);

        alpha_INT[0][i] = alpha1i;
        alpha_INT[1][i] = alpha2i;
        alpha_INT[2][i] = alpha3i;
        alpha_INT[3][i] = alpha4i;
        alpha_INT[4][i] = alpha5i;
        alpha_INT[5][i] = alpha6i;
        alpha_INT[6][i] = rv2*creal(Zi*conj(zi));
        alpha_INT[7][i] = rv2*cimag(Zi*conj(zi));

        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;
    }

    //Maximum
    double max_INT_VS_OFS[8];
    double max_INT_VS_FFT[8];
    double max_INT[8];
    for(int j = 0; j < 8; j ++)
    {
        max_INT_VS_OFS[j] = alpha_INT_VS_OFS[j][0];
        max_INT_VS_FFT[j] = alpha_INT_VS_FFT[j][0];
        max_INT[j]        = alpha_INT[j][0];
    }

    for(int i = 1 ; i < N; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            if(max_INT_VS_OFS[j] < alpha_INT_VS_OFS[j][i]) max_INT_VS_OFS[j] = alpha_INT_VS_OFS[j][i];
            if(max_INT_VS_FFT[j] < alpha_INT_VS_FFT[j][i]) max_INT_VS_FFT[j] = alpha_INT_VS_FFT[j][i];
            if(max_INT[j] < alpha_INT[j][i]) max_INT[j] = alpha_INT[j][i];
        }
    }

    //Print
    cout << "---------------------------------------" << endl;
    cout << "       FFT vs OFS                      " << endl;
    cout << "---------------------------------------" << endl;
    if(type == 0) cout << "The numerical integration is used to compute the reference values " << endl;
    else          cout << "The Fourier expansions computed in qbtbp(int) are used to compute the reference values " << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(3);
    cout << "i         max REF           max REF vs OFS       max REF vs FFT      " << endl;
    for(int j = 0; j < 8; j++)
    {
        cout << j << "      " << max_INT[j] << "           " << max_INT_VS_OFS[j]  << "            " <<  max_INT_VS_FFT[j]  << endl;
    }
    cout << "---------------------------------------" << endl;

    //Conclusion
    cout << "Conclusion: " << endl;
    cout << "- OFS and FFT coefficients both represent quite well the real integrated system (see comparison with type == 0). " << endl;
    cout << "- However, the FFT coefficients are better suited to represent the approximated system given by the Fourier expansions computed in qbtbp(int) (see comparison with type == 1). " << endl;
    cout << "- This is ought to the approximations made to compute z(OFS)^{-a}, a > 0, when computing the OFS coefficients." << endl;

    //Memory release
    free_dvector(alpha_OFS, 0, 8*(nf+1)-1);
    free_dvector(alpha_FFT, 0, 8*(nf+1)-1);
    free_dmatrix(alpha_INT_VS_OFS, 0, 7, 0, N-1);
    free_dmatrix(alpha_INT       , 0, 7, 0, N-1);
    free_dmatrix(alpha_INT_VS_FFT, 0, 7, 0, N-1);
}

//-----------------------------------------------------------------------------
// Evaluating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(Ofsc& zt, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*zt.evaluate(n*t);
}

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(Ofsc& zt, Ofsc& ztdot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*(ztdot.evaluate(n*t) + I*ni*zt.evaluate(n*t));
}

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(Ofsc& zt, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*( 2*I*ni*ztdot.evaluate(n*t) - ni*ni*zt.evaluate(n*t) + ztddot.evaluate(n*t));
}


//-----------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< order; i++)
    {
        //From trigo formulae (fastest)
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
        //From GSL
        //cR[i] =  gsl_sf_cos((i+1)*omega*t);
        //sR[i] =  gsl_sf_sin((i+1)*omega*t);
        //Native C
        //cR[i] =  cos((i+1)*omega*t);
        //sR[i] =  sin((i+1)*omega*t);
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13) alpha[l] = evaluateOdd(t, omega , order, header, sR);    //Odd funtions (alpha_2,5,8,10,12,14)
        else  alpha[l] = evaluateEven(t, omega, order, header, cR);                                                  //Even functions
    }
}

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13)
        {
            alpha[l] = evaluateOddDerivative(t, omega , order, header, cR); //Odd funtions (alpha_2,5,8,10,12,14)
        }
        else  alpha[l] = evaluateEvenDerivative(t, omega, order, header, sR);   //Even functions
    }
}


//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(double t, double omega, int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*cR[i-1];//even type
    result += coef[0];
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double t, double omega,  int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += -omega*i*coef[i]*sR[i-1];//even type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(double t, double omega,  int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*sR[i-1]; //odd type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double t, double omega,  int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += omega*i*coef[i]*cR[i-1];//odd type
    return result;
}



//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Deprecated routine for the computation of the alphas and betas
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computation of the coefficients alphas and betas that appear in the QBCP vector field. NOT USED ANYMORE
 *         These computations are performed via algebraic manipulation on the Fourier series, contrary to other similar routines
 *         Such as qbtbp_ofs_fft_alpha (for the alpha coefficients) that makes use of FFT procedures.
 **/
void qbtbp_ofs_storage(   Ofts< Ofsd > &zr_ofts,     //zr_ofts = normalized Earth-Moon motion
                          Ofts< Ofsd > &Zr_ofts,     //Zr_ofts = normalized Sun-(Earth+Moon) motion
                          Ofts< Ofsd > &z1,     //z1 = \bar{zr_ofts}
                          Ofts< Ofsd > &z2,     //z2 = \bar{zr_ofts}^(-3/2)
                          Ofts< Ofsd > &z3,     //z3 = zr_ofts^(-1/2)
                          Ofts< Ofsd > &z4,     //z4 = z2*z3 = zr_ofts^(-1/2)*\bar{zr_ofts}^(-3/2)
                          Ofts< Ofsd > &z5,     //z5 = \bar{zr_ofts}^(-1/2)
                          Ofts< Ofsd > &z6,     //z6 = z3*z5 = (zr_ofts*\bar{zr_ofts})^(-1/2)
                          Ofts< Ofsd > &a1,     //a1 = mu*exp(itheta)*zr_ofts
                          Ofts< Ofsd > &a2,     //a2 = epsilon*a1
                          Ofts< Ofsd > &a3,     //a3 = Zr_ofts-mu*exp(itheta)*epsilon*zr_ofts
                          Ofts< Ofsd > &a4,     //a4 = a3^(-1/2).
                          Ofts< Ofsd > &a5,     //a5 = \bar{a3}
                          Ofts< Ofsd > &a6,     //a6 = a5^(-3/2).
                          Ofts< Ofsd > &a7,     //a7 = a4*a6,
                          Ofts< Ofsd > &b1,     //b1 = (1-mu)*exp(ithetb)*zr_ofts
                          Ofts< Ofsd > &b2,     //b2 = epsilon*b1
                          Ofts< Ofsd > &b3,     //b3 = Zr_ofts+(1-mu)*exp(ithetb)*epsilon*zr_ofts
                          Ofts< Ofsd > &b4,     //b4 = b3^(-1/2)
                          Ofts< Ofsd > &b5,     //b5 = \bar{b3}
                          Ofts< Ofsd > &b6,     //b6 = b5^(-3/2)
                          Ofts< Ofsd > &b7,     //b7 = b4*b6,
                          Ofts< Ofsd > &Pm,     //Pm = -ms/as^2*exp(-itheta)*b7 + ms/as^2*exp(-itheta)*a7 - z4
                          Ofts< Ofsd > &Qm,     //Qm = -ns^2*mu*b7 - ns^2*(1-mu)*a7
                          Ofsd  &Pfm,            //Pfm = Pm(j) at order n
                          Ofsd  &Qfm,            //Qfm = Qm(j) at order n
                          Ofsd  &ufm,            //ufm = zr_ofts(j) at order n
                          Ofsd  &vfm,            //vfm = Zr_ofts(j) at order n
                          Ofsd &sigma1,          //sigma1 = exp(-itheta)
                          Ofsd &sigma2,          //sigma2 = exp(+itheta)
                          int nf,                       //order of the Fourier expansions
                          QBCP_L& qbcp_l)             //QBCP
{
    //Physical param
    double n  = qbcp_l.us_em.n;
    double mu = qbcp_l.us_em.mu_EM;
    double ns = qbcp_l.us_em.ns;
    double as = qbcp_l.us_em.as;
    double ms = qbcp_l.us_em.ms;
    double muSE = qbcp_l.us_em.mu_SE;

    double eps = 1.0/as;                     //length epsilon
    int order  = zr_ofts.getOrder();         //order of the Taylor expansion
    int nv     = zr_ofts.getNV();            //number of variables of the Taylor expansion
    int cnv    = zr_ofts.getCVariables();    //number of variables of the Fourier coefficient

    //----------------------------------------------------------------------
    // Init
    //----------------------------------------------------------------------
    Ofsd bj(nf);
    Ofsd cj(nf);
    Ofsc bjc(nf);
    Ofsc cjc(nf);

    //r^2  = 1/alpha1
    Ofsd r2(nf);                             //Ofs double version
    Ofts< Ofsd >  zr2(nv, order, cnv, nf);   //Ofts version


    //The pre computations
    Ofsc d1(nf);
    Ofsc d2(nf);
    Ofsc d3(nf);
    Ofsc d31(nf);
    Ofsc d4(nf);
    Ofsc d41(nf);
    Ofsc d5(nf);
    Ofsc d51(nf);
    Ofsc d6(nf);
    Ofsc d61(nf);
    Ofsc d7(nf);
    Ofsc d71(nf);
    Ofsc d8(nf);
    Ofsc d81(nf);
    Ofsc d9(nf);
    Ofsc d91(nf);
    Ofsc d10(nf);
    Ofsc d101(nf);
    Ofsc d11(nf);
    Ofsc d111(nf);
    Ofsc d12(nf);
    Ofsc d121(nf);


    Ofsc d1r(nf);
    Ofsc d1i(nf);
    Ofsc d3r(nf);
    Ofsc d3i(nf);
    Ofsc d4r(nf);
    Ofsc d4i(nf);
    Ofsc d5r(nf);
    Ofsc d5i(nf);
    Ofsc d6r(nf);
    Ofsc d6i(nf);
    Ofsc d7r(nf);
    Ofsc d7i(nf);
    Ofsc d8r(nf);
    Ofsc d8i(nf);

    Ofsc d9r(nf);
    Ofsc d9i(nf);
    Ofsc d10r(nf);
    Ofsc d10i(nf);
    Ofsc d11r(nf);
    Ofsc d11i(nf);
    Ofsc d12r(nf);
    Ofsc d12i(nf);

    //The alphas
    Ofsd alpha1(nf);
    Ofsd alpha6(nf);
    Ofs<cdouble > alpha1c(nf);
    Ofs<cdouble > alpha2c(nf);
    Ofs<cdouble > alpha3c(nf);
    Ofs<cdouble > alpha4c(nf);
    Ofs<cdouble > alpha5c(nf);
    Ofs<cdouble > alpha6c(nf);
    Ofs<cdouble > alpha7c(nf);
    Ofs<cdouble > alpha8c(nf);
    //Note that alpha9 is not defined in Andreu (the numerotation seems to have jumped one term)
    //Note alpha10...19 are defined in Andreu but not used here
    //For consistency, their name have not been used for other variables
    Ofs<cdouble > alpha20c(nf);
    Ofs<cdouble > alpha21c(nf);
    Ofs<cdouble > alpha22c(nf);
    Ofs<cdouble > alpha23c(nf);
    Ofs<cdouble > alpha24c(nf);
    Ofs<cdouble > alpha25c(nf);
    Ofs<cdouble > alpha26c(nf);
    Ofs<cdouble > alpha27c(nf);
    Ofs<cdouble > alpha28c(nf);

    //sigma1c & sigma2c
    Ofs<cdouble > sigma1c(nf);
    Ofs<cdouble > sigma2c(nf);
    doubleToComplex(sigma1, sigma1c);
    doubleToComplex(sigma2, sigma2c);

    //----------------------------------------------------------------------
    // bj and cj
    //----------------------------------------------------------------------
    //From ots to ofs
    fts2fs(&bj, zr_ofts, eps);
    fts2fs(&cj, Zr_ofts, eps);
    //Store in txt files
    ofstream curentStream;
    curentStream.open("data/qbtbp/bj.txt");
    curentStream << "bj: \n" << bj << endl;
    curentStream.close();
    curentStream.open("data/qbtbp/cj.txt");
    curentStream << "cj: \n" << cj << endl;
    curentStream.close();

    //----------------------------------------------------------------------
    // bj and cj in complex form
    //----------------------------------------------------------------------
    //From double to cdouble
    doubleToComplex(bj, bjc);
    doubleToComplex(cj, cjc);
    //Store in txt files
    curentStream.open("data/qbtbp/bjc.txt");
    curentStream <<  bjc << endl;
    curentStream.close();
    curentStream.open("data/qbtbp/cjc.txt");
    curentStream <<  cjc << endl;
    curentStream.close();


    //----------------------------------------------------------------------
    // zr and Zr in complex form and derivatives. We recall that:
    //sigma1 = exp(-itheta)
    //sigma2 = exp(+itheta
    //----------------------------------------------------------------------
    Ofsc zr(bjc);
    Ofsc Zr(cjc);
    //dot
    Ofsc zrdot(zr);
    zrdot.dot(n);
    Ofsc Zrdot(Zr);
    Zrdot.dot(n);
    //ddot
    Ofsc zrddot(zrdot);
    zrddot.dot(n);
    Ofsc Zrddot(Zrdot);
    Zrddot.dot(n);
    //conj
    Ofsc zrc(zr);
    zrc.conjugate();
    Ofsc Zrc(Zr);
    Zrc.conjugate();

    //Order 0:
    //----------------------------------
    //d1 = zr*conj(zr)
    d1.ofs_prod(zr, zrc);
    //d2 = Z*conj(Z) = as^2*Zr*conj(Zr)
    d2.ofs_mprod(Zr, Zrc, as*as+0.0*I);
    //d31 = zr*conj(Zr)
    d31.ofs_prod(zr, Zrc);
    //d3 = z*conj(Z) = as*sigma2*zr*conj(Zr)
    d3.ofs_mprod(sigma2c, d31, as+0.0*I);
    //d41 = Zr*conj(zr)
    d41.ofs_prod(Zr, zrc);
    //d4 = Z*conj(z) = as*sigma1*Zr*conj(zr)
    d4.ofs_mprod(sigma1c, d41, as+0.0*I);
    //d1r = real(d1), d1i = imag(d1)
    realPart(d1, d1r);
    imagPart(d1, d1i);
    //d3r = real(d3), d3i = imag(d3)
    realPart(d3, d3r);
    imagPart(d3, d3i);
    //d4r = real(d4), d4i = imag(d4)
    realPart(d4, d4r);
    imagPart(d4, d4i);


    //First order:
    //----------------------------------
    //d51 = dot(zr) + i*zr
    d51.ofs_fsum(zrdot, 1.0+0.0*I, zr, I);
    //d5 = dot(z)*conj(z)
    d5.ofs_prod(d51, zrc);

    //d61 = sigma2*d51 = exp(int)*(dot(zr) + i*zr)
    d61.ofs_prod(sigma2c, d51);
    //d6 = dot(z)*conj(Z) = as*exp(int)*(dot(zr) + i*zr)*conj(Zr)
    d6.ofs_mprod(d61, Zrc, as+0.0*I);

    //d71 = dot(Zr) + i*ns*Zr
    d71.ofs_fsum(Zrdot, 1.0+0.0*I, Zr, I*ns);
    //d7 = dot(Z)*conj(Z) = as*as*(dot(Zr) + i*ns*Zr)*conj(Zr)
    d7.ofs_mprod(d71, Zrc, as*as+0.0*I);

    //d81 = sigma1*d71 = exp(-int)*(dot(Zr) + i*ns*Zr)
    d81.ofs_sprod(sigma1c, d71);
    //d8 = dot(Z)*conj(z) = as*exp(-int)*(dot(Zr) + i*ns*Zr)*conj(z)
    d8.ofs_mprod(d81, zrc, as+0.0*I);

    //d5r = real(d5), d5i = imag(d5)
    realPart(d5, d5r);
    imagPart(d5, d5i);
    //d6r = real(d6), d6i = imag(d6)
    realPart(d6, d6r);
    imagPart(d6, d6i);
    //d7r = real(d7), d7i = imag(d7)
    realPart(d7, d7r);
    imagPart(d7, d7i);
    //d8r = real(d8), d8i = imag(d8)
    realPart(d8, d8r);
    imagPart(d8, d8i);

    //Second order:
    //----------------------------------
    //d91 = ddot(zr) +2*i*dot(zr) - zr
    d91.ofs_fsum(zrddot, 1.0+0.0*I, zrdot, 2.0*I);
    d91.ofs_fsum(d91, 1.0+0.0*I, zr, -1.0+0.0*I);
    //d9 = ddot(z)*conj(z)
    d9.ofs_prod(d91, zrc);

    //d101 = sigma2*d91 = exp(int)*(ddot(zr) +2*i*dot(zr) - zr)
    d101.ofs_prod(sigma2c, d91);
    //d10 = ddot(z)*conj(Z) = as*exp(int)*(ddot(zr) +2*i*dot(zr) - zr)*conj(Zr)
    d10.ofs_mprod(d101, Zrc, as+0.0*I);

    //d111 = ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr
    d111.ofs_fsum(Zrddot, 1.0+0.0*I, Zrdot, 2*ns*I);
    d111.ofs_fsum(d111, 1.0+0.0*I, Zr, -ns*ns+0.0*I);
    //d11 = ddot(Z)*conj(Z) = as*as*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)*conj(Zr)
    d11.ofs_mprod(d111, Zrc, as*as+0.0*I);

    //d121 = sigma1*d111 = exp(-int)*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)
    d121.ofs_prod(sigma1c, d111);
    //d12 = ddot(Z)*conj(z) = as*exp(-int)*(ddot(Zr) +2*i*ns*dot(Zr) - ns*ns*Zr)*conj(zr)
    d12.ofs_mprod(d121, zrc, as+0.0*I);


    //d9r = real(d9), d9i = imag(d9)
    realPart(d9, d9r);
    imagPart(d9, d9i);

    //d10r = real(d10), d10i = imag(d10)
    realPart(d10, d10r);
    imagPart(d10, d10i);

    //d11r = real(d11), d11i = imag(d11)
    realPart(d11, d11r);
    imagPart(d11, d11i);

    //d12r = real(d12), d12i = imag(d12)
    realPart(d12, d12r);
    imagPart(d12, d12i);


    //----------------------------------------------------------------------
    // The alphas
    //----------------------------------------------------------------------
    //r^2 (square radius) = 1/alpha1 = alpha27
    //----------------------------------------
    alpha27c.ccopy(d1);
    //Storage in txt file
    curentStream.open((qbcp_l.cs.F_COEF+"r2.txt").c_str());
    curentStream << alpha27c << endl;
    curentStream.close();
    //Storage in txt file (alpha27 !!)
    ofs_sst(alpha27c, qbcp_l.cs.F_COEF+"alpha27", 1, "");
    //-----------------------------

    //alpha6
    //-----------------------------
    //z4 = 1/r2 = z^(-1/2)*zbar^(-1/2)
    z4.zero();
    z4.ofts_sprod(z3, z5);
    //alpha6
    fts2fs(&alpha6, z4, eps);
    //alpha6c
    doubleToComplex(alpha6, alpha6c);
    //Storage in txt file
    ofs_sst(alpha6c, qbcp_l.cs.F_COEF+"alpha6", 1, "");
    //-----------------------------

    //alpha25 = c2 * alpha6
    //-----------------------------
    alpha25c.ofs_mult(alpha6c, cn(qbcp_l, 2)+0.0*I);
    //Storage in txt file
    ofs_sst(alpha25c, qbcp_l.cs.F_COEF+"alpha25", 1, "");
    //-----------------------------

    //alpha1
    //-----------------------------
    //alpha1 = alpha6*alpha6
    alpha1.ofs_sprod(alpha6,alpha6);
    //alpha1c
    doubleToComplex(alpha1, alpha1c);
    //Storage in txt file
    ofs_sst(alpha1c, qbcp_l.cs.F_COEF+"alpha1", 1, "");
    //-----------------------------

    //alpha2
    //-----------------------------
    //alpha2c = -real(dot(z)*conj(z))*1/r2
    alpha2c.ofs_mprod(alpha1c, d5r, -1.0+0.0*I);
    //force the zero order to zero (odd function)
    alpha2c.setCoef(0.0, 0);
    //force the odd nature
    for(int l = -nf; l< 0; l++) alpha2c.setCoef(+0.0*I-alpha2c.ofs_getCoef(-l), l);
    //Storage in txt file
    ofs_sst(alpha2c, qbcp_l.cs.F_COEF+"alpha2", 0, "");
    //-----------------------------

    //alpha3
    //-----------------------------
    //alpha3c = +imag(dot(z)*conj(z))*1/r2
    alpha3c.ofs_sprod(alpha1c, d5i);
    //Storage in txt file
    ofs_sst(alpha3c, qbcp_l.cs.F_COEF+"alpha3", 1, "");
    //-----------------------------

    //alpha4
    //-----------------------------
    alpha4c+= (+0.0*I-ms/(1+ms))*d12r;
    //Storage in txt file
    ofs_sst(alpha4c, qbcp_l.cs.F_COEF+"alpha4", 1, "");
    //-----------------------------

    //alpha5
    //-----------------------------
    alpha5c+= (+0.0*I-ms/(1+ms))*d12i;
    //Storage in txt file
    ofs_sst(alpha5c, qbcp_l.cs.F_COEF+"alpha5", 0, "");
    //-----------------------------

    //alpha7
    //-----------------------------
    //alpha7
    alpha7c.ofs_sprod(alpha1c, d4r);
    //Storage in txt file
    ofs_sst(alpha7c, qbcp_l.cs.F_COEF+"alpha7", 1, "");
    //-----------------------------

    //alpha8
    //-----------------------------
    alpha8c.ofs_sprod(alpha1c, d4i);
    alpha8c.setCoef(0.0, 0);
    //Storage in txt file
    ofs_sst(alpha8c, qbcp_l.cs.F_COEF+"alpha8", 0, "");
    //-----------------------------


    //----------------------------------------------------------------------
    // alpha28 = 1/sqrt(tilde(alpha7)^2 + tilde(alpha8)^2)
    //----------------------------------------------------------------------
    Ofs<cdouble > alpha78c(nf);
    alpha78c.ofs_smprod(alpha7c, alpha7c, 1.0/(as*as)+0.0*I);
    alpha78c.ofs_smprod(alpha8c, alpha8c, 1.0/(as*as)+0.0*I);

    //Normalized version
    Ofs<cdouble > e1(nf);
    e1.ofs_smult(alpha78c, 1.0/alpha78c.ofs_getCoef(0));

    //e1^(-3/2)
    Ofs<cdouble > e2(nf);
    e2.ofs_epspow(e1, -1.0/2+0.0*I);
    //e2*cpow(alpha78c.getCoef(0), -1.0/2)
    alpha28c.ofs_smult(e2, cpow(alpha78c.ofs_getCoef(0), -1.0/2+0.0*I));

    //Test: uncomment to check
    //---------------------------------------------------------------------
//    double alpha = -1.0/2;
//    double t = 1.0;
//    cdouble a78i   = alpha78c.evaluate(n*t);
//    cdouble a78i32 = cpow(a78i, alpha);
//    cdouble alpha28ci    = alpha28c.evaluate(n*t);
//    cout <<  setiosflags(ios::scientific) << setprecision(15);
//    cout << "------------------------------------------------" << endl;
//    cout << " Test of the inverse alpha28 = pow(a78, -1/2)" << endl;
//    cout << "------------------------------------------------" << endl;
//    cout << "alpha28(t) = " << creal(alpha28ci) << "  " << cimag(alpha28ci) << endl;
//    cout << "pow(a78, -1/2) = " << creal(a78i32) << "  " << cimag(a78i32) << endl;
//    cout << "alpha28(t)-pow(a78, -1/2) = " << creal(alpha28ci)-creal(a78i32)
//         << "  " << cimag(alpha28ci)-cimag(a78i32) << endl;
//    cout << "------------------------------------------------" << endl;
    //---------------------------------------------------------------------

    //Storage in txt file
    ofs_sst(alpha28c, qbcp_l.cs.F_COEF+"alpha28", 1, "");
    //-----------------------------

    //alpha20 = alpha6*alpha28*ms/as^2
    //-----------------------------
    //alpha20
    alpha20c.ofs_smprod(alpha6c, alpha28c, ms/(as*as)+0.0*I);
    //Storage in txt file
    ofs_sst(alpha20c, qbcp_l.cs.F_COEF+"alpha20", 1, "");
    //-----------------------------


    //alpha21 = alpha28*alpha28
    //-----------------------------
    //alpha21
    alpha21c.ofs_prod(alpha28c, alpha28c);
    //Storage in txt file
    ofs_sst(alpha21c, qbcp_l.cs.F_COEF+"alpha21", 1, "");
    //-----------------------------

    //alpha22 = alpha28*alpha28*tilde(alpha7) = alpha21*alpha7/as
    //-----------------------------
    //alpha22
    alpha22c.ofs_mprod(alpha7c, alpha21c, 1.0/as+0.0*I);
    //Storage in txt file
    ofs_sst(alpha22c, qbcp_l.cs.F_COEF+"alpha22", 1, "");
    //-----------------------------

    //alpha23 = alpha28*alpha28*tilde(alpha8) = alpha21*alpha8/as
    //-----------------------------
    //alpha23
    alpha23c.ofs_mprod(alpha8c, alpha21c, 1.0/as+0.0*I);
    //Storage in txt file
    ofs_sst(alpha23c, qbcp_l.cs.F_COEF+"alpha23", 0, "");
    //-----------------------------


    //Derivatives useful in the forthcoming computations
    Ofs<cdouble > alpha1dotc(alpha1c);
    Ofs<cdouble > alpha2dotc(alpha2c);
    Ofs<cdouble > alpha3dotc(alpha3c);
    alpha1dotc.dot(n);
    alpha2dotc.dot(n);
    alpha3dotc.dot(n);

    //Spares
    Ofsc sp1(nf);
    Ofsc sp2(nf);
    Ofsc sp3(nf);
    Ofsc sp4(nf);
    Ofsc sp5(nf);
    Ofsc sp6(nf);
    Ofsc sp7(nf);
    Ofsc sp8(nf);
    //sp1 = alpha2dotc*alpha1c - alpha1dotc*alpha2c
    sp1.ofs_prod(alpha2dotc, alpha1c);
    sp2.ofs_prod(alpha1dotc, alpha2c);
    sp1.ofs_fsum(sp1,  1.0+0.0*I, sp2,  -1.0+0.0*I); //fsum used in place
    //sp4 = alpha27^2*(alpha2dotc*alpha1c - alpha1dotc*alpha2c)
    sp3.ofs_prod(alpha27c, alpha27c);
    sp4.ofs_prod(sp3, sp1);
    //sp5 = alpha2^2+alpha3^2
    sp5.ofs_prod(alpha2c, alpha2c);
    sp6.ofs_prod(alpha3c, alpha3c);
    sp5.ofs_fsum(sp5, 1.0+0.0*I, sp6, 1.0+0.0*I); //fsum used in place
    //sp7 = alpha27*(alpha2^2+alpha3^2)
    sp7.ofs_prod(alpha27c, sp5);
    //sp8 = - sp4 - sp7
    //BEWARE: the definition of alpha24 has changed, so alpha6c*c1 is not added anymore!
    sp8.ofs_smult(sp7,  -1.0+0.0*I);
    sp8.ofs_smult(sp4,  -1.0+0.0*I);

    //alpha24
    //-----------------------------
    alpha24c.ofs_fsum(sp8, qbcp_l.cs.c1+0.0*I, alpha4c, 1.0/qbcp_l.cs.gamma+0.0*I);
    //Storage in txt file
    ofs_sst(alpha24c, qbcp_l.cs.F_COEF+"alpha24", 1, "");
    //WARNING: alpha9 is set equal to alpha24 for now (needs improvement, avoid double definition!)
    ofs_sst(alpha24c, qbcp_l.cs.F_COEF+"alpha9", 1, "");
    //-----------------------------

    //Spares are used again
    //sp1 = alpha3dotc*alpha1c - alpha1dotc*alpha3c
    sp1.ofs_prod(alpha3dotc, alpha1c);
    sp2.ofs_prod(alpha1dotc, alpha3c);
    sp1.ofs_fsum(sp1, 1.0+0.0*I, sp2, -1.0+0.0*I); //fsum used in place
    //sp4 = alpha27^2*(alpha3dotc*alpha1c - alpha1dotc*alpha3c)
    sp4.ofs_prod(sp3, sp1);

    //alpha26
    //-----------------------------
    alpha26c.ofs_fsum(sp4, qbcp_l.cs.c1+0.0*I, alpha5c, 1.0/qbcp_l.cs.gamma+0.0*I);
    //Storage in txt file
    ofs_sst(alpha26c, qbcp_l.cs.F_COEF+"alpha26", 0, "");
    //WARNING: alpha10 is set equal to alpha26 for now (needs improvement, avoid double definition!)
    ofs_sst(alpha26c, qbcp_l.cs.F_COEF+"alpha10", 0, "");
    //-----------------------------


    //---------------------------------------------------------------------------------------------------------------------------------
    // The gammas
    //---------------------------------------------------------------------------------------------------------------------------------
    Ofsc g1(nf);
    Ofsc g2(nf);
    Ofsc g3(nf);
    Ofsc g31(nf);
    Ofsc g32(nf);
    Ofsc g4(nf);
    Ofsc g41(nf);
    Ofsc g42(nf);
    Ofsc g51(nf);
    Ofsc g61(nf);

    //g1 = d7r + mu*mu*d5r - mu*(d6r+d8r);
    g1.ofs_fsum(d7r, 1.0+0.0*I, d5r,  mu*mu+0.0*I);
    g1.ofs_fsum(g1,  1.0+0.0*I, d6r, -mu+0.0*I);
    g1.ofs_fsum(g1,  1.0+0.0*I, d8r, -mu+0.0*I);

    //g2 = d7i + mu*mu*d5i - mu*(d6i+d8i);
    g2.ofs_fsum(d7i, 1.0+0.0*I, d5i,  mu*mu+0.0*I);
    g2.ofs_fsum(g2 , 1.0+0.0*I, d6i, -mu+0.0*I);
    g2.ofs_fsum(g2 , 1.0+0.0*I, d8i, -mu+0.0*I);

    //g31 = re(ddot(Z)*conj(zh)) = re(ddot(Z)*conj(Z)) - mu*re(ddot(Z)*conj(z)) = d11r - mu*d12r
    g31.ofs_fsum(d11r, 1.0+0.0*I, d12r, -mu+0.0*I);
    //g32 = re(ddot(z)*conj(zh)) = re(ddot(z)*conj(Z)) - mu*re(ddot(z)*conj(z)) = d10r - mu*d9r
    g32.ofs_fsum(d10r, 1.0+0.0*I, d9r, -mu+0.0*I);
    //g3 = mu/(1-mu+ms)*(ms/(1+ms)*g31 + (1-mu)*g32)
    g3.ofs_fsum(g31, mu/(1-mu+ms)*ms/(1+ms)+0.0*I, g32, mu/(1-mu+ms)*(1-mu)+0.0*I);

    //g41 = im(ddot(Z)*conj(zh)) = im(ddot(Z)*conj(Z)) - mu*im(ddot(Z)*conj(z)) = d11i - mu*d12i
    g41.ofs_fsum(d11i, 1.0+0.0*I, d12i, -mu+0.0*I);
    //g42 = im(ddot(z)*conj(zh)) = im(ddot(z)*conj(Z)) - mu*im(ddot(z)*conj(z)) = d10i - mu*d9i
    g42.ofs_fsum(d10i, 1.0+0.0*I, d9i, -mu+0.0*I);
    //g4 = mu/(1-mu+ms)*(ms/(1+ms)*g41 + (1-mu)*g42)
    g4.ofs_fsum(g41, mu/(1-mu+ms)*ms/(1+ms)+0.0*I, g42, mu/(1-mu+ms)*(1-mu)+0.0*I);

    //g51 = re(z*zh) = re(z*conj(Z)) - mu*Re(z*conj(z))
    g51.ofs_fsum(d3r, 1.0+0.0*I, d1r, -mu+0.0*I);
    //g61 = im(z*zh) = im(z*conj(Z)) - mu*im(z*conj(z))
    g61.ofs_fsum(d3i, 1.0+0.0*I, d1i, -mu+0.0*I);
    //Note: g5 and g6 are not computed since:
    // - the expansion 1/h2 = beta1 is needed
    // - beta7 = gamma5 and beta8 = gamma6, so they can be computed as betas instead of gammas

    //----------------------------------------------------------------------
    // The betas
    //----------------------------------------------------------------------
    Ofts< Ofsd >  zh(nv, order, cnv, nf);       //zh = Z-mu*z
    Ofsd zh0(nf);                               //order 0

    //zh+= -mu*zr_ofts
    zh.ofts_smult_u(zr_ofts, -mu);
    //zh+= Z*as*sigma1
    zh.ofts_smult_tu(Zr_ofts, sigma1, as);

    //Take the order 0
    zh0.ccopy(zh.getTerm(0)->getCoef(0));

    //zp = 1 - mu/as*exp(int)
    Ofsd zp(nf);
    zp.setCoef(1.0, 0);
    zp.setCoef(-mu/as, 1);

    //zp1= 1/zp
    Ofsd zp1(nf);
    zp1.ofs_epspow(zp, -1.0);

    //zh0inv = 1/as*exp(int)*1/zp = inverse of order 0 of zh
    Ofsd zh0inv(nf);
    zh0inv.ofs_sprod(zp1, sigma2);
    zh0inv *= 1.0/as;

    //zhinv= zh^-1
    Ofts< Ofsd >  zhinv(nv, order, cnv, nf);
    zhinv.pows(zh,  zh0inv, zh0inv, -1.0);

    //To Ofs
    Ofsd zh_ofs(nf);
    Ofsd zhinv_ofs(nf);
    fts2fs(&zh_ofs, zh, eps);
    fts2fs(&zhinv_ofs, zhinv, eps);


    //Test: uncomment to check
    //---------------------------------------------------------------------
    //    alpha = -1.0;
    //    t = 1.0;
    //    cout <<  setiosflags(ios::scientific) << setprecision(15);
    //    cout << "------------------------------------------------" << endl;
    //    cout << " Test of the inverse of zh = Z - mu * z         " << endl;
    //    cout << "------------------------------------------------" << endl;
    //    cout << "zhinv_ofs(t) = " << creal(zhinv_ofs.evaluate(t)) << "  " << cimag(zhinv_ofs.evaluate(t)) << endl;
    //    cout << "1/zh_ofs(t) = " << creal(cpow(zh_ofs.evaluate(t), alpha)) << "  " << cimag(cpow(zh_ofs.evaluate(t), alpha)) << endl;
    //    cout << "zhinv_ofs(t)-1/zh_ofs(t) = " << creal(zhinv_ofs.evaluate(t))-creal(cpow(zh_ofs.evaluate(t), alpha))
    //         << "  " << cimag(zhinv_ofs.evaluate(t))-cimag(cpow(zh_ofs.evaluate(t), alpha)) << endl;
    //    cout << "------------------------------------------------" << endl;
    //---------------------------------------------------------------------


    //zh is divided by as*exp(-int) i.e. is multiplied by sigma2/as
    Ofts< Ofsd >  zhn(nv, order, cnv, nf);
    zhn.ofts_smult_tu(zh, sigma2, 1.0/as);

    //zp1= pow(zp, -1/2)
    Ofsd zp12(nf);
    zp12.ofs_epspow(zp, -1.0/2);

    //zh2inv= zhn^-1/2
    Ofts< Ofsd >  zh2inv(nv, order, cnv, nf);
    zh2inv.pows(zhn, zp1, zp12, -1.0/2);

    //To Ofs
    Ofsd zh2inv_ofs(nf);
    fts2fs(&zh2inv_ofs, zh2inv, eps);

    //beta6
    //-----------------------------
    Ofsc zh2invc(nf);
    Ofsc zh2invconj(nf);
    doubleToComplex(zh2inv_ofs, zh2invc);
    doubleToComplex(zh2inv_ofs, zh2invconj);
    zh2invconj.conjugate();

    //beta6 = zh2invc*zh2invc
    Ofsc beta6c(nf);
    beta6c.ofs_prod(zh2invc, zh2invconj);
    beta6c.ofs_mult(beta6c, 1.0/as+0.0*I);

    //Test: uncomment to check
    //---------------------------------------------------------------------
    //alpha = -1.0/2;
    //t = 1.0;
    //cdouble zhi = zh_ofs.evaluate(t);
    //cdouble hi2  = zhi*conj(zhi);
    //cdouble hi = cpow(hi2, alpha);
    //cout <<  setiosflags(ios::scientific) << setprecision(15);
    //cout << "------------------------------------------------" << endl;
    //cout << " Test of the inverse beta6 = 1/h" << endl;
    //cout << "------------------------------------------------" << endl;
    //cout << "beta6(t) = " << creal(beta6c.evaluate(t)) << "  " << cimag(beta6c.evaluate(t)) << endl;
    //cout << "pow(zh_ofs(t), -1/2) = " << creal(hi) << "  " << cimag(hi) << endl;
    //cout << "beta6(t)-pow(zh_ofs(t), -1/2) = " << creal(beta6c.evaluate(t))-creal(hi)
    //     << "  " << cimag(beta6c.evaluate(t))-cimag(hi) << endl;
    //cout << "------------------------------------------------" << endl;
    //---------------------------------------------------------------------


    //Storage in txt file
    ofs_sst(beta6c, qbcp_l.cs.F_COEF+"beta6", 1, "");
    //-----------------------------


    //beta1
    //-----------------------------
    Ofsc zhinvc(nf);
    Ofsc zhinvconj(nf);
    doubleToComplex(zhinv_ofs, zhinvc);
    doubleToComplex(zhinv_ofs, zhinvconj);
    zhinvconj.conjugate();

    //beta1 = zhinvc*zhinvc
    Ofsc beta1c(nf);
    beta1c.ofs_prod(zhinvc, zhinvconj);

    //Storage in txt file
    ofs_sst(beta1c, qbcp_l.cs.F_COEF+"beta1", 1, "");
    //-----------------------------


    //beta2
    //-----------------------------
    Ofsc beta2c(nf);
    //beta2 = - gamma1/h2 = -beta1*gamma1
    beta2c.ofs_mprod(g1, beta1c, -1.0+0.0*I);
    beta2c.setCoef(0.0, 0);
    //Storage in txt file
    ofs_sst(beta2c, qbcp_l.cs.F_COEF+"beta2", 0, "");
    //-----------------------------


    //beta3
    //-----------------------------
    Ofsc beta3c(nf);
    //beta3 = + gamma2/h2 = +beta1*gamma2
    beta3c.ofs_prod(g2, beta1c);
    //Storage in txt file
    ofs_sst(beta3c, qbcp_l.cs.F_COEF+"beta3", 1, "");
    //-----------------------------


    //beta4
    //-----------------------------
    Ofsc beta4c(g3);
    //Storage in txt file
    ofs_sst(beta4c, qbcp_l.cs.F_COEF+"beta4", 1, "");
    //-----------------------------

    //beta5
    //-----------------------------
    Ofsc beta5c(g4);
    //Storage in txt file
    ofs_sst(beta5c, qbcp_l.cs.F_COEF+"beta5", 0, "");
    //-----------------------------

    //beta7
    //-----------------------------
    //beta7 = muSE - 1 + beta1*g51
    Ofsc beta7c(nf);
    beta7c.setCoef(muSE - 1, 0);
    beta7c.ofs_sprod(beta1c, g51);
    //Storage in txt file
    ofs_sst(beta7c, qbcp_l.cs.F_COEF+"beta7", 1, "");
    //-----------------------------

    //beta8
    //-----------------------------
    //beta8 = beta1*g61
    Ofsc beta8c(nf);
    beta8c.ofs_sprod(beta1c, g61);
    beta8c.setCoef(0.0, 0);
    //Storage in txt file
    ofs_sst(beta8c, qbcp_l.cs.F_COEF+"beta8", 0, "");
    //-----------------------------

}



//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Backup void qbtbp_ofs and qbtbp_ots
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
///
//   \brief Solution of the qbtbp in Ots format (obsolete)
//
void qbtbp_ots (Ofts< Ots<double> > &ofts_z, Ofts< Ots<double> > &ofts_Z, double n, double ms, double as, double ns)
{
    int order = ofts_z.getOrder();
    int nv = ofts_z.getNV();
    int corder = ofts_z.getCOrder();
    int cnv = ofts_z.getCVariables();
    //Here nf = corder/2 because nf is the order of the Ofs coefficients (Fourier form), not the Ots (Taylor form)
    int nf = ofts_z.getCOrder()/2;

    double mu = 0.0121505816;
    double eps = 1.0/as;


    //---------------------------------------------------------
    //Order 0
    //---------------------------------------------------------
    ofts_z.setCoef0(1.0,0, 0);
    ofts_Z.setCoef0(1.0,0, 0);


    //---------------------------------------------------------
    //Order 1
    //One may choose to set it at the beginning or to compute it in the loop.
    //---------------------------------------------------------
    //Fourier series
    Ofsd u1(nf);
    Ofsd v1(nf);

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
    Ofsd sigma1(nf);
    sigma1.setCoef(1.0, -1);
    Ots<double> s1(cnv, corder);
    fs2ts(&s1, sigma1);


    //sigma2 = exp(+itheta)
    Ofsd sigma2(nf);
    sigma2.setCoef(1.0, 1);
    Ots<double> s2(cnv, corder);
    fs2ts(&s2, sigma2);

    //epsilon = 0+epsilon+0 (order 1 of Taylor serie)
    Ofts< Ots<double> > epsilon(nv, order, cnv, corder);
    epsilon.setCoef0(1.0,1,0);

    //Order 0
    int m;
    for(m = 0; m<=0 ; m++)
    {
        z1.ccopy(ofts_z, m);
        //z1 = \bar{z1}
        z1.conjugate(m);
        //z2 = \bar{z}^^(-3/2)
        z2.ofts_pows(z1, -3.0/2, m);
        //z3 = z^(-1/2)
        z3.ofts_pows(ofts_z, -1.0/2, m);
        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);


        //a1 = mu*exp(itheta)*z
        a1.ofts_smult_tu(ofts_z, s2, mu, m);
        //a2 = epsilon*a1
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z
        a3.ofts_smult_u(ofts_Z, 1.0, m);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);
        //a4 = a3^(-1/2)
        a4.ofts_pows(a3, -1.0/2, m);
        //a5 = \bar{a3}
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2)
        a6.ofts_pows(a5, -3.0/2, m);
        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);


        //b1 = (1-mu)*exp(ithetb)*z
        b1.ofts_smult_tu(ofts_z, s2, (1.0-mu), m);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        b3.ofts_smult_u(ofts_Z, 1.0, m);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);
        //b4 = b3^(-1/2)
        b4.ofts_pows(b3, -1.0/2, m);
        //b5 = \bbr{b3}
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2)
        b6.ofts_pows(b5, -3.0/2, m);
        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);


        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, s1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, s1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);

    }

    Ofsd  Pfm(nf);
    Ofsd  Qfm(nf);

    Ofsd ufm(nf);
    Ofsd vfm(nf);

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
        if(m>1) z2.ofts_pows(z1, -3.0/2, m-1);
        z2.ofts_pows(z1, -3.0/2, m);
        //z3 = z^(-1/2)
        if(m>1) z3.ofts_pows(ofts_z, -1.0/2, m-1);
        z3.ofts_pows(ofts_z, -1.0/2, m);
        //z4 = z2*z3;
        z4.ofts_sprod(z2, z3, m);

        //cout << "Order: " << m << endl;
        //cout << "z1: \n" << z1 << endl;
        //cout << "z2: \n" << z2 << endl;



        //a1 = mu*exp(itheta)*z at order m-1
        if(m>1) a1.ofts_smult_tu(ofts_z, s2, mu, m-1);
        //a2 = epsilon*a1 at order m
        a2.ofts_sprod(a1, epsilon, m);
        //a3+= Z
        if(m>1) a3.ofts_smult_u(ofts_Z, 1.0, m-1);
        //a3+= -mu*exp(itheta)*epsilon*z
        a3.ofts_smult_u(a2, -1.0, m);
        //a4 = a3^(-1/2)
        if(m>1) a4.ofts_pows(a3, -1.0/2, m-1);
        a4.ofts_pows(a3, -1.0/2, m);
        //a5 = \bar{a3}
        if(m>1)
        {
            a5.ccopy(a3, m-1);
            a5.conjugate(m-1);
        }
        a5.ccopy(a3, m);
        a5.conjugate(m);
        //a6 = a5^(-3/2)
        if(m>1) a6.ofts_pows(a5, -3.0/2, m-1);
        a6.ofts_pows(a5, -3.0/2, m);
        //a7 = a4*a6;
        a7.ofts_sprod(a4,a6, m);

//        cout << "\n---------------------------\n" << endl;
//        cout << "Order: " << m << endl;
//        cout << "a3: \n" << a3 << endl;
//        cout << "a4: \n" << a4 << endl;

        //b1 = (1-mu)*exp(ithetb)*z
        if(m>1) b1.ofts_smult_tu(ofts_z, s2, (1.0-mu), m-1);
        //b2 = epsilon*b1
        b2.ofts_sprod(b1, epsilon, m);
        //b3+= Z
        if(m>1) b3.ofts_smult_u(ofts_Z, 1.0, m-1);
        //b3+= (1-mu)*exp(ithetb)*epsilon*z
        b3.ofts_smult_u(b2, 1.0, m);
        //b4 = b3^(-1/2)
        if(m>1) b4.ofts_pows(b3, -1.0/2, m-1);
        b4.ofts_pows(b3, -1.0/2, m);
        //b5 = \bbr{b3}
        if(m>1)
        {
            b5.ccopy(b3, m-1);
            b5.conjugate(m-1);
        }
        b5.ccopy(b3, m);
        b5.conjugate(m);
        //b6 = b5^(-3/2)
        if(m>1) b6.ofts_pows(b5, -3.0/2, m-1);
        b6.ofts_pows(b5, -3.0/2, m);
        //b7 = b4*b6;
        b7.ofts_sprod(b4,b6, m);

        //Pm += -ms/as^2*exp(-itheta)*b7
        Pm.ofts_smult_tu(b7, s1, -ms/(as*as),m);
        //Pm +=  ms/as^2*exp(-itheta)*a7
        Pm.ofts_smult_tu(a7, s1,  ms/(as*as),m);
        //Pm +=  -z4
        Pm.ofts_smult_u(z4, -1.0, m);


        //Qm += -ns^2*mu*b7
        Qm.ofts_smult_u(b7, -ns*ns*mu, m);
        //Qm += -ns^2*(1-mu)*a7
        Qm.ofts_smult_u(a7, -ns*ns*(1-mu), m);

        //-------------------------
        // Solving equations
        //-------------------------
        Pfm.ts2fs(Pm.getTerm(m)->getCoef(0));
        Qfm.ts2fs(Qm.getTerm(m)->getCoef(0));

        //Order 0
        ufm.setCoef(-1.0/3*Pfm.ofs_getCoef(0), 0);  //u0 = -1/3*p0
        vfm.setCoef(-1.0/(3*ns*ns)*Qfm.ofs_getCoef(0), 0); //v0 = -1/(3*ns^2)*p0

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
            uj = ( k2*Pfm.ofs_getCoef(-j) - k3*Pfm.ofs_getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, -j);
            //uj
            uj = (-k3*Pfm.ofs_getCoef(-j) + k1*Pfm.ofs_getCoef(j))/(k1*k2-k3*k3);
            ufm.setCoef(uj, j);

            //v-j
            vj = ( l2*Qfm.ofs_getCoef(-j) - l3*Qfm.ofs_getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, -j);
            //vj
            vj = (-l3*Qfm.ofs_getCoef(-j) + l1*Qfm.ofs_getCoef(j))/(l1*l2-l3*l3);
            vfm.setCoef(vj, j);

        }


        //Update z and Z
        fs2ts(ofts_z.getCoef(m,0), ufm);
        fs2ts(ofts_Z.getCoef(m,0), vfm);

        cout << "Order " << m << " completed." << endl;
    }

    //bj and cj
    Ofsd bj(nf);
    fts2fs(&bj, ofts_z, eps);
    Ofsd cj(nf);
    fts2fs(&cj, ofts_Z, eps);


    //cout << ofts_z << endl;

    ofstream bjstream("data/bj.txt");
    bjstream << "bj: \n" << bj << endl;
    bjstream.close();
    ofstream cjstream("data/cj.txt");
    cjstream << "cj: \n" << cj << endl;
    cjstream.close();

}

*/
//-----------------------------------------------------------------------------
// FFT from NR in C
//-----------------------------------------------------------------------------
/**
 * \brief Calculates the Fourier transform of a set of n real-valued data points. Not currently used.
 *
 * Replaces this data (which
 * is stored in array data[0..n-1] ) by the positive frequency half of its complex Fourier transform.
 * The real-valued first and last components of the complex transform are returned as elements
 * data[1] and data[2] , respectively. n must be a power of 2. This routine also calculates the
 * inverse transform of a complex data array if it is the transform of real data. (Result in this case
 * must be multiplied by 2/n .)
 **/
void realft(double data[], unsigned long n, int isign)
{
    void four1(double data[], unsigned long nn, int isign);
    unsigned long i,i1,i2,i3,i4,np3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta; //Double precision for the trigonometric recurrences.
    theta=3.141592653589793/(double) (n>>1);//Initialize the recurrence.
    if (isign == 1)
    {
        c2 = -0.5;
        four1(data,n>>1,1);  //The forward transform is here.
    }
    else
    {
        c2=0.5; //Otherwise set up for an inverse transform.
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2; i<=(n>>2); i++)
    {
        //Case i=1 done separately below.
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]); //The two separate transforms are separated out of data.
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i; //Here they are recombined to form
        data[i2]=h1i+wr*h2i+wi*h2r; //the true transform of the original real data.
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr; //The recurrence.
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1)
    {
        //Squeeze the first and last data to-
        //gether to get them all within the
        //original array.
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    }
    else
    {
        //This is the inverse transform for the
        //case isign=-1.
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n>>1,-1);
    }

}
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/**
 * \brief Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
 * data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
 * data is a complex array of length nn or, equivalently, a real array of length 2*nn . nn MUST
 * be an integer power of 2 (this is not checked for!).
 **/
void four1(double data[], unsigned long nn, int isign)
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;
    n=nn << 1;
    j=1;
    for (i=1; i<n; i+=2)
    {
        //This is the bit-reversal section of the routine.
        if (j > i)
        {

            SWAP(data[j],data[i]);   //Exchange the two complex numbers.
            SWAP(data[j+1],data[i+1]);
        }
        m=nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    //Here begins the Danielson-Lanczos section of the routine.
    mmax=2;
    while (n > mmax)  //Outer loop executed log 2 nn times.
    {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax); // Initialize the trigonometric recurrence.
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1; m<mmax; m+=2)  //Here are the two nested inner loops.
        {
            for (i=m; i<=n; i+=istep)
            {
                //This is the Danielson-Lanczos for-mula:
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr; //Trigonometric recurrence.
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

