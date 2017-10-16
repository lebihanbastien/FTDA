#include "pm.h"
#include "pmt.h"
#include <iostream>
#include <fstream>
extern "C" {
    #include "nrutil.h"
}

/**
 * \file pm.cpp
 * \brief Implements the parameterization method for both autonomous and non-autonomous Hamiltonian vector fields of the Sun-Earth-Moon problem.
 *         Deprecated. Use routines in pmt.cpp instead
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  More precisely, the parameterization method is applied to compute a parameterization of the central manifold of the Earth-Moon libration point L1,2:
 *
 *   - The parameterization is given as Fourier-Taylor expansions in the case of non-autonomous Hamiltonian vector fields (QBCP, BCP).
 *   - The parameterization is given as pure    Taylor expansions in the case of     autonomous Hamiltonian vector fields (CRTBP).
 */

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Main routine (deprecated, use pmt instead, from pmt.h)
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Compute the parameterization of the central manifold of the dynamical equivalent of the Earth-Moon libration points L1,2, up to a given order.
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_GS...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_GS):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 **/
void pm(int OutputEachOrder, int Output)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Application of the parameterization method     " << endl;
    cout << "  on the non-autonomous vector field of the QBCP   " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(4);
    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_COC   = SEML.cs.F_COC;
    string F_COEF  = SEML.cs.F_COEF;

    cout << "F_COEF = " << F_COEF << endl;
    cout << "F_COC = "  << F_COC << endl;
    cout << "F_GS = "   << F_GS << endl;

    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    matrix<Ofsc> P(6,6);      //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> Q(6,6);      //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    vector<Ofsc> V(6);        //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> Xe(3);       //The vectors Xe, Xm, Xs, modified position
    vector<Ofsc> Xm(3);       //of the primaries used to compute the potential:
    vector<Ofsc> Xs(3);       //Xf = True Xf - V, for f = e, m, s
    Ofsc ILe(OFS_ORDER);      //The Ofs ILe, ILm and ILs, modified inverse orbit radii
    Ofsc ILm(OFS_ORDER);      //of the primaries used to compute the potential
    Ofsc ILs(OFS_ORDER);      //ILf = 1.0/||Xf||, for f = e, m, s
    matrix<Ofsc> PC(6,6);     //PC=P(theta)*C
    matrix<Ofsc> PCdot(6,6);  //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> CQ(6,6);     //CQ=Cinv*Pinv=inv(PC)
    vector<Ofsc> Vdot(6);     //Vdot=dot(V)

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Initialization of the alphas
    //------------------------------------------
    int noc = SEML.numberOfCoefs;
    vector<Ofsc> alpha(noc);
    //Switch on the model to initialize the alphas
    switch(SEML.model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    case Csts::ERTBP:
    {
        //Read from txt files
        ifstream readStream;
        string ss1;
        for(int i = 0; i < noc; i++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
            readOFS_txt(alpha[i], F_COEF+"alpha"+ss1+"_fft");
        }

        break;
    }

    case Csts::CRTBP:
    {
        cout << "pm. The use of the RTBP has been detected when initializing the alphas." << endl;
        //All coeffs to zero
        for(int j = 0; j < noc;  j++) alpha[j].setCoef(0.0, 0);
        //Non-null coeffs
        alpha[0].setCoef(1.0, 0);               //alpha1  = 1.0
        alpha[2].setCoef(1.0, 0);               //alpha3  = 1.0
        alpha[5].setCoef(1.0, 0);               //alpha6  = 1.0
        alpha[8].setCoef(SEML.cs.mu, 0);  //alpha9  = Xe
        alpha[10].setCoef(SEML.cs.mu, 0); //alpha11 = Xm
        alpha[12].setCoef(-SEML.cs.c1, 0);   //alpha13 = -c1
        break;
    }

    default:
    {
        cout << "error in pm. Unknown model." << endl;
        return;
    }

    }



    //------------------------------------------
    // Getting DB from files
    //------------------------------------------
    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |
    gsl_matrix_complex *DB = gsl_matrix_complex_calloc(6, 6);
    if(SEML.model == Csts::QBCP || SEML.model == Csts::BCP)
    {
        //Reading an OFS from a text file
        string filename = F_COC+"DB";
        glsc_matrix_complex_read(DB, filename);
    }
    else //RTBP
    {
        cout << "pm. The use of the RTBP has been detected when initializing DB." << endl;
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, c2;
        c2 = SEML.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);

        cout << "omega1 = " << om1 << endl;
        cout << "omega2 = " << om2 << endl;
        //--------------------
        //Init DB
        //--------------------
        gsl_matrix_complex_set(DB, 0, 0, gslc_complex(+om1*I));
        gsl_matrix_complex_set(DB, 1, 1, gslc_complex(+la1+0.0*I));
        gsl_matrix_complex_set(DB, 2, 2, gslc_complex(+om2*I));
        gsl_matrix_complex_set(DB, 3, 3, gslc_complex(-om1*I));
        gsl_matrix_complex_set(DB, 4, 4, gslc_complex(-la1+0.0*I));
        gsl_matrix_complex_set(DB, 5, 5, gslc_complex(-om2*I));
    }
    cout << setprecision(15) << endl;
    //------------------------------------------
    //DF0 = DB in matrix<Ofsc> format
    //------------------------------------------
    Ofsc BUX(OFS_ORDER);
    matrix<Ofsc> DF0(6, 6);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    DF0.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    DF0.setCoef(BUX, 1, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    DF0.setCoef(BUX, 2, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    DF0.setCoef(BUX, 3, 3);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    DF0.setCoef(BUX, 4, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    DF0.setCoef(BUX, 5, 5);


    //------------------------------------------
    // Setting H, Hinv, Lambda, LamdaL, LambdaN from DB
    // A this step, the central manifold is selected.
    //
    // TO-DO: something more generic to be able
    // to compute the central-stable and central-unstable manifold
    //
    //------------------------------------------
    gsl_matrix_complex *H     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Hinv  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *La    = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *LaL   = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex *LaN   = gsl_matrix_complex_calloc(2, 2);
    //H = (L N)
    gsl_matrix_complex *L    = gsl_matrix_complex_calloc(6, 4);
    gsl_matrix_complex *N    = gsl_matrix_complex_calloc(6, 2);

    //------------------------------------------
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L N)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //------------------------------------------
    gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
    gsl_matrix_complex_set(H, 1, 4, gsl_matrix_complex_get(DB, 1, 1));
    gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
    gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
    gsl_matrix_complex_set(H, 4, 5, gsl_matrix_complex_get(DB, 4, 4));
    gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

    //------------------------------------------
    //H2 = H in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> H2(6, 6);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    H2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    H2.setCoef(BUX, 1, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    H2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    H2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    H2.setCoef(BUX, 4, 5);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    H2.setCoef(BUX, 5, 3);

    //------------------------------------------
    //    |  iw1 0    0    0    |
    //    |  0   0    0    0    |
    //    |  0   iw3  0    0    |
    //L = |  0   0   -iw1  0    |
    //    |  0   0    0    0    |
    //    |  0   0    0   -iw3  |
    //------------------------------------------
    gsl_matrix_complex_view H11 = gsl_matrix_complex_submatrix (H , 0 , 0 , 6 , 4 );
    gsl_matrix_complex_memcpy(L, &H11.matrix);

    //matrix version of L
    matrix<Ofsc> L2(6, 4);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    L2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    L2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    L2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    L2.setCoef(BUX, 5, 3);


    //------------------------------------------
    //    |  0   0  |
    //    |  w2  0  |
    //    |  0   0  |
    //N = |  0   0  |
    //    |  0  -w2 |
    //    |  0   0  |
    //------------------------------------------
    H11 = gsl_matrix_complex_submatrix (H , 0 , 4 , 6 , 2 );
    gsl_matrix_complex_memcpy(N, &H11.matrix);


    //------------------------------------------
    //       |1/iw1  0    0      0     0     0    |
    //       | 0     0   1/iw3   0     0     0    |
    //       | 0     0    0    -1/iw1  0     0    |
    //Hinv = | 0     0    0      0     0  -1/iw3  |
    //       | 0    1/w2  0      0     0     0    |
    //       | 0     0    0      0   -1/w2   0    |
    //------------------------------------------
    gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
    gsl_matrix_complex_set(Hinv, 4, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
    gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
    gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
    gsl_matrix_complex_set(Hinv, 5, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
    gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

    //------------------------------------------
    //Hinv2 = Hinv in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> Hinv2(6, 6);
    BUX.zero();
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    Hinv2.setCoef(BUX, 0, 0);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    Hinv2.setCoef(BUX, 4, 1);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    Hinv2.setCoef(BUX, 1, 2);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    Hinv2.setCoef(BUX, 2, 3);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    Hinv2.setCoef(BUX, 5, 4);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    Hinv2.setCoef(BUX, 3, 5);

    //------------------------------------------
    //                     | LaL  T  |
    //Lambda = Hinv*DB*H = | 0  LaN  | with T = 0 in this context
    //
    //------------------------------------------
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(6, 6);
    //AUX = DB*H
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , H , zero_c , AUX );
    //B = Hinv*AUX
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Hinv , AUX , zero_c , La);
    //------------------------------------------
    gsl_matrix_complex_view La11 = gsl_matrix_complex_submatrix (La , 0 , 0 , 4 , 4 );   //S11
    gsl_matrix_complex_memcpy(LaL, &La11.matrix);
    gsl_matrix_complex_view La12 = gsl_matrix_complex_submatrix (La , 4 , 4 , 2 , 2 );   //S12
    gsl_matrix_complex_memcpy(LaN, &La12.matrix);

    //------------------------------------------
    // Initialisation of the TFC manifold ("W hat")
    //------------------------------------------
    cdouble db;
    gsl_complex dbc;
    //vector of size 6
    vector<Oftsc> Wh(6);
    //The order 0 is equal to zero, the order 1 is equal to L, which is diagonal
    //----------------------
    //     |   iw1*s1 |
    //     |     0    |
    //     |   iw3*s2 |
    //Wh = |  -iw1*s3 |
    //     |     0    |
    //     |  -iw3*s4 |
    //----------------------
    // Wh[0] = L11*s1 = L[0][0]*s[0]
    dbc = gsl_matrix_complex_get(L, 0, 0);
    db = gslc_complex(dbc);
    Wh[0].setCoef0(db, 1, 0);
    // Wh[1] = 0
    // Wh[2] = L32*s2 = L[2][1]*s[1]
    dbc = gsl_matrix_complex_get(L, 2, 1);
    db = gslc_complex(dbc);
    Wh[2].setCoef0(db, 1, 1);
    // Wh[3] = L43*s3 = L[3][2]*s[2]
    dbc = gsl_matrix_complex_get(L, 3, 2);
    db = gslc_complex(dbc);
    Wh[3].setCoef0(db, 1, 2);
    // Wh[4] = 0
    // Wh[5] = L64*s4 = L[5][3]*s[3]
    dbc = gsl_matrix_complex_get(L, 5, 3);
    db = gslc_complex(dbc);
    Wh[5].setCoef0(db, 1, 3);

    //------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------
    vector<Oftsc> W(6);   //W = TFC(Wh) (in NC coordinates)
    vector<Oftsc> Wt(6);  //Wt = W - V  (in NC coordinates)
    //Applying the COC (with and without zero order) @order 0 and 1
    for(int i = 0; i < 2; i++)
    {
        applyCOC(PC,  V,  Wh, W,  i, 1);   //order zero is added
        applyCOC(PC,  V,  Wh, Wt, i, 0);   //order zero is NOT added
    }

    //------------------------------------------
    // Differential of Wh (6*4)
    //------------------------------------------
    //Contribution of order < k to k
    matrix<Oftsc> DWh(6,4);  //TFC
    //Contribution of order <=k to k
    matrix<Oftsc> C_DW(6,4);   //NC

    //------------------------------------------
    // Initialisation of the complete VF
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> FWh(6);  //TFC
    vector<Oftsc> FW(6);   //NC
    vector<Oftsc> FW2(6);  //Temp
    vector<Oftsc> FW3(6);  //Temp
    //Contribution of order <=k to k
    vector<Oftsc> C_FW(6);   //NC

    //------------------------------------------
    // Initialisation of the reduced VF
    //------------------------------------------
    vector<Oftsc> fh(4);     //TFC
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        fh[i].setCoef0(db, 1, i);
    }
    //LaLs = fh @ order 0 and 1
    vector<Oftsc> LaLs(4);
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        //db = 1.0;
        LaLs[i].setCoef0(db, 1, i);
    }

    //------------------------------------------
    // Initialisation of DW times fh
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> DWfh(6);  //TFC
    //Contribution of order <=k to k
    vector<Oftsc> C_DWf(6);   //NC

    //------------------------------------------
    // Initialisation of dot(W)
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> Whdot(6);      //TFC
    vector<Oftsc> Wdot(6);       //NC
    //Contribution of order <=k to k
    vector<Oftsc> C_Wdot(6);     //NC


    //------------------------------------------
    // Initialisation of the Legendre objects
    //------------------------------------------
    //Potential
    vector<Oftsc> En(OFTS_ORDER+2);     //Earth
    vector<Oftsc> Sn(OFTS_ORDER+2);     //Moon
    vector<Oftsc> Mn(OFTS_ORDER+2);     //Sun
    vector<Oftsc> Un(OFTS_ORDER+2);     //Earth+Moon+Sun
    //Potential derivatives
    vector<Oftsc> Enx(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Eny(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Enz(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Mnx(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mny(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mnz(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Snx(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Sny(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Snz(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Unx(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Uny(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Unz(OFTS_ORDER+2);    //Earth+Moon+Sun

    //------------------------------------------
    // Initialisation of additionnal OFTS objects (see comments)
    // Function of Wh ("W hat")
    // Updated in updatePotentialExpansion routine
    //------------------------------------------
    Oftsc Rho2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);   //=xt^2+yt^2+zt^2
    Oftsc xse(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xe*xt+ye*yt
    Oftsc xsm(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xm*xt+ym*yt
    Oftsc xss(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xs*xt+ys*yt

    //------------------------------------------
    // Initialisation of various objects
    // Used to solve the pm equations
    //------------------------------------------
    vector<Oftsc> E(6);    //[DWh<m fh<m]m-[Fh(Wh<m)]m
    vector<Oftsc> eta(6);  //eta = Hinv*E
    vector<Oftsc> xi(6);   //xi = Hinv*Wh at order >=2
    cdouble lap;       //contains LaN
    cdouble lar;       //contains LaL*r for each monomial
    cdouble aux;       //temp obj
    cdouble bux;       //temp obj
    int kv[REDUCED_NV];       //contains the indices for each monomial

    cout << "  pm. Initialization of the recurrence.  " << endl;

    //-------------------------------------------------------------
    // Initialization of the recurrence
    //-------------------------------------------------------------
    //The potential is updated: Un, Unx, Uny, Unz
    // @order 0, 1 and 2
    //Wt is used to update the Potential
    updatePotentialExpansion(Wt, Rho2, En,
                             Mn, Sn, Un,
                             Enx, Eny, Enz,
                             Mnx, Mny, Mnz,
                             Snx, Sny, Snz,
                             Unx, Uny, Unz,
                             Xe, Xm, Xs,
                             ILe, ILm, ILs,
                             xse, xsm, xss,
                             SEML, 2);

    cout << "  pm. End of initialisation.  " << endl;

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ order 0 and 1
    //------------------------------------------
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    for(int k = 0; k < 2 ; k++)
    {
        //vector field
        //-----------------
        //applyVF(alpha, W, FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);
        //Copy for consistency testing
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);

        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, k);

        //dot(W)
        //-----------------
        for(int i = 0; i < 6; i++) Wdot[i].dot(W[i],   SEML.us.n, k);     //in NC
        for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], SEML.us.n, k);     //in TFC
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, k);     //in NC, complete contribution

        //----------------------
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, k);   //in NC
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, k);   //in NC, complete contribution

        // DWh times f @order k
        //----------------------
        //smvprod_t(DW,  fh, DWf,  k);     //in NC
        smvprod_t(DWh, fh, DWfh, k);       //in TFC
        smvprod_t(C_DW,  fh, C_DWf,  k);   //in NC, complete contribution
    }

    cout << "  pm. End of order 0/1.  " << endl;

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ orders >= 2
    //------------------------------------------
    for(int m = 2; m <= OFTS_ORDER; m++)
    {

        cout << "pm. start of order " << m << endl;
        tic();
        //------------------------------------------
        //The potential is updated: Un, Unx, Uny, Unz
        //Wt is used to update the Potential
        //What is updated: contributions of order < m to order m
        //------------------------------------------
        if(m > 2) updatePotentialExpansion(Wt, Rho2, En,
                                               Mn, Sn, Un,
                                               Enx, Eny, Enz,
                                               Mnx, Mny, Mnz,
                                               Snx, Sny, Snz,
                                               Unx, Uny, Unz,
                                               Xe, Xm, Xs,
                                               ILe, ILm, ILs,
                                               xse, xsm, xss,
                                               SEML, m);

        //------------------------------------------
        // Applying the VF: [F(W<m)]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);

        //------------------------------------------
        //FWh = COCinv(FW)
        //------------------------------------------
        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, m);

        //cout << "pm. end of applyInvDotCOC." << endl;

        //------------------------------------------
        // dot(W<m)
        // It is NOT updated here: the terms of order < m
        // dot not influence dot(W) at order m
        //------------------------------------------
        //for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], qbcp.n, m);
        //for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], qbcp.n, m);
        //cout << "pm. end of dot." << " in " << toc() << " s. " << endl;

        //------------------------------------------
        // Differential of Wh
        // It is NOT updated here, but at the end of the loop.
        // Indeed, Wh & W contain no order m at this point!
        //------------------------------------------

        //------------------------------------------
        // [DW<m times f<m]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        smvprod_t(DWh, fh, DWfh, m);
        //smvprod_t(DW, fh, DWf, m);


        //------------------------------------------
        // We have [Fh(Wh<m)]m and [DWh<m fh<m]m
        // Update of eta at order m
        //------------------------------------------
        // [E]m = [DWh<m fh<m]m - [Fh(Wh<m)]m
        for(int i = 0 ; i< 6; i++) E[i].ofts_fsum_u(FWh[i], -1.0+0.0*I, DWfh[i], 1.0+0.0*I, m);
        //[eta]m = [Hinv*E]m;
        smvprod_u(Hinv2, E, eta, m);



        //------------------------------------------
        // Equation in the normal space:
        // Last two components of eta
        //------------------------------------------
        //Sum on the normal components of eta (last n-d components)
        for(int p = REDUCED_NV; p<= Csts::NV-1; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
            //Reset kv
            kv[0] = m;
            for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Reset aux
                aux = 0.0+0.0*I;
                //Update aux = laL.kv
                for(int ii = 0; ii < REDUCED_NV; ii++)
                {
                    lar = gslc_complex(gsl_matrix_complex_get(LaL, ii, ii));
                    aux += lar*(kv[ii]+0.0*I);
                }

                //Update xi
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    bux = 0.0*I-eta[p].getCoef(m, k)->ofs_getCoef(j)/(SEML.us.n*I*j - (lap - aux));
                    xi[p].getCoef(m, k)->setCoef(bux, j);
                }
                //Update kv
                if(k < FTDA::nmon(REDUCED_NV, m)-1)  FTDA::prxkt(kv, REDUCED_NV);
            }
        }

        //------------------------------------------
        // Equation in the tangential space:
        // Graph style is used by default
        // As a consequence:
        // - xi^p_m = 0
        // - fh^p_m = -eta^p_m
        //------------------------------------------
        //Sum on the tangential components of eta (first d components)
        for(int p = 0; p< REDUCED_NV; p++)
        {
            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Update fh
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    bux = eta[p].getCoef(m, k)->ofs_getCoef(j);
                    fh[p].getCoef(m, k)->setCoef(0.0*I-bux, j);
                }
            }
        }

        //------------------------------------------
        // Wh = H2*xi @ order m
        // Recall that H2 = H in matrix format
        //------------------------------------------
        smvprod_u(H2, xi, Wh, m);

        //------------------------------------------
        // Update, W, Wt and W2 @ order m
        //------------------------------------------
        applyCOC(PC,  V,  Wh, W,  m, 1);     //order zero is added
        applyCOC(PC,  V,  Wh, Wt, m, 0);     //order zero is NOT added


        //------------------------------------------
        // Update the differential to take into account the new terms @ order m
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, m);


        //------------------------------------------
        // Update the potential expansions at order m
        //------------------------------------------
        initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, SEML, m);

        //------------------------------------------------------------------------------------
        //
        // Note: after this point, only necessary for testing, can be discarded for speed
        //
        //------------------------------------------------------------------------------------
        //------------------------------------------
        // update the VF
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, m);

        //------------------------------------------
        // dot(W)
        //------------------------------------------
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, m);

        //------------------------------------------
        // Differential of Wh
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, m);

        //------------------------------------------
        // DWh times f @order m
        //------------------------------------------
        smvprod_t(C_DW,  fh, C_DWf,  m);

        cout << "pm. end of update testing." << " in " << toc() << " s. "<< endl;
        cout << "  pm. End of order " << m << " in " << toc() << " s. " << endl;

        //If txt files have to be updated at each order
        if(OutputEachOrder)
        {
            tic();
            //print W & FW
            vector_fprinf(W,     F_GS+"W/W");
            vector_fprinf(Wt,    F_GS+"W/Wt");
            vector_fprinf(Wh,    F_GS+"W/Wh");
            vector_fprinf(FW,    F_GS+"FW/FW");
            vector_fprinf(C_FW,  F_GS+"FW/C_FW");
            vector_fprinf(FWh,   F_GS+"FW/FWh");
            //fprintf Wdot
            vector_fprinf(C_Wdot,  F_GS+"Wdot/C_Wdot");
            //Reduced vector field
            vector_fprinf(fh, F_GS+"rvf/fh");
            //Potential
            vector_fprinf(Unx, F_GS+"Potential/Unx");
            vector_fprinf(Uny, F_GS+"Potential/Uny");
            vector_fprinf(Unz, F_GS+"Potential/Unz");
            vector_fprinf(En, F_GS+"Potential/En");
            vector_fprinf(Mn, F_GS+"Potential/Mn");
            vector_fprinf(Sn, F_GS+"Potential/Sn");
            vector_fprinf(Un, F_GS+"Potential/Un");
            //DW*vf
            vector_fprinf(C_DWf,  F_GS+"DWf/C_DWf");
            cout << "  pm. End of printing of order " << m  <<  " in " << toc() << " s. " << endl;
        }
    }


    //------------------------------------------
    // Printing
    //------------------------------------------
    //If txt files have to be updated at the very last order
    if(Output && !OutputEachOrder)
    {
        tic();

        //------------------------------------------
        // Txt
        //------------------------------------------
        //print W & FW
        vector_fprinf(W,     F_GS+"W/W");
        vector_fprinf(Wt,    F_GS+"W/Wt");
        vector_fprinf(Wh,    F_GS+"W/Wh");
        vector_fprinf(FW,    F_GS+"FW/FW");
        vector_fprinf(C_FW,  F_GS+"FW/C_FW");
        vector_fprinf(FWh,   F_GS+"FW/FWh");
        //fprintf Wdot
        vector_fprinf(C_Wdot,  F_GS+"Wdot/C_Wdot");
        //Reduced vector field
        vector_fprinf(fh, F_GS+"rvf/fh");
        //Potential
        vector_fprinf(Unx, F_GS+"Potential/Unx");
        vector_fprinf(Uny, F_GS+"Potential/Uny");
        vector_fprinf(Unz, F_GS+"Potential/Unz");
        vector_fprinf(En,  F_GS+"Potential/En");
        vector_fprinf(Mn, F_GS+"Potential/Mn");
        vector_fprinf(Sn, F_GS+"Potential/Sn");
        vector_fprinf(Un, F_GS+"Potential/Un");
        //DW*vf
        vector_fprinf(C_DWf,  F_GS+"DWf/C_DWf");



        //------------------------------------------
        // Binary
        //------------------------------------------
        writeVOFTS_bin(W,     F_GS+"W/W");
        writeVOFTS_bin(Wt,    F_GS+"W/Wt");
        writeVOFTS_bin(Wh,    F_GS+"W/Wh");
        writeVOFTS_bin(FW,    F_GS+"FW/FW");
        writeVOFTS_bin(C_FW,  F_GS+"FW/C_FW");
        writeVOFTS_bin(FWh,   F_GS+"FW/FWh");
        //fprintf Wdot
        writeVOFTS_bin(C_Wdot,  F_GS+"Wdot/C_Wdot");
        //Reduced vector field
        writeVOFTS_bin(fh, F_GS+"rvf/fh");
        //Potential
        writeVOFTS_bin(Unx, F_GS+"Potential/Unx");
        writeVOFTS_bin(Uny, F_GS+"Potential/Uny");
        writeVOFTS_bin(Unz, F_GS+"Potential/Unz");
        writeVOFTS_bin(En,  F_GS+"Potential/En");
        writeVOFTS_bin(Mn, F_GS+"Potential/Mn");
        writeVOFTS_bin(Sn, F_GS+"Potential/Sn");
        writeVOFTS_bin(Un, F_GS+"Potential/Un");
        //DW*vf
        writeVOFTS_bin(C_DWf,  F_GS+"DWf/C_DWf");


        cout << "  pm. End of printing" <<  " in " << toc() << " s. " << endl;
    }

}

/**
 *  \brief Same routine as pm, with additional test.
 *
 *  Invariance Equation test:
 *
 *  - IE_ots: checks that the invariance equation is satisfied by W(s,t) at each order.
 *  - IEh_ots: checks that the invariance equation is satisfied by Wh(s,t) at each order.
 *  - xi5_ots; consistency test of the cohomological equations in the coordinates s (CCM).
 *  - Wh6_ots; consistency test of the cohomological equations in the coordinates Wh (TFC).
 *  - D1_ots & D2_ots: splitted consistency test  of the cohomological equations in the coordinates Wh.
 *    More precisly:
 *
 *      D1 = T_FWh - DF(0)Wh - FWh, with T_FWh the complete vector field, FWh the vector field without the contribution of the last order in Wh
 *      D2 = T_DWfh - LaLs*dWh[p]/dx[j] - L*fh - DWfh, with DWfh the complete derivatives of the vector field and DWfh the derivatives of the
 *                                                     vector field without the contribution of the last order in Wh
 *
 */
void pmTested()
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Application of the parameterization method     " << endl;
    cout << "  on the non-autonomous vector field of the QBCP  " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(4);
    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_GS;
    string F_COEF  = SEML.cs.F_COEF;
    string F_COC   = SEML.cs.F_COC;


    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(6,6);
    //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    matrix<Ofsc> Q(6,6);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(6);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential:
    //Xf = True Xf - V, for f = e, m, s
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);
    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    //ILf = 1.0/||Xf||, for f = e, m, s
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);
    //PC=P(theta)*C
    matrix<Ofsc> PC(6,6);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(6,6);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(6,6);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(6);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Initialization of the alphas
    //------------------------------------------
    int noc = SEML.numberOfCoefs;
    vector<Ofsc> alpha(noc);
    //Switch on the model to initialize the alphas
    switch(SEML.model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    case Csts::ERTBP:
    {
        //Read from txt files
        ifstream readStream;
        string ss1;
        for(int i = 0; i < noc; i++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
            readOFS_txt(alpha[i], F_COEF+"alpha"+ss1+"_fft");
        }

        break;
    }

    case Csts::CRTBP:
    {
        cout << "pm. The use of the RTBP has been detected when initializing the alphas." << endl;
        //All coeffs to zero
        for(int j = 0; j < noc;  j++) alpha[j].setCoef(0.0, 0);
        //Non-null coeffs
        alpha[0].setCoef(1.0, 0);               //alpha1  = 1.0
        alpha[2].setCoef(1.0, 0);               //alpha3  = 1.0
        alpha[5].setCoef(1.0, 0);               //alpha6  = 1.0
        alpha[8].setCoef(SEML.us_em.mu_EM, 0);  //alpha9  = Xe
        alpha[10].setCoef(SEML.us_em.mu_EM, 0); //alpha11 = Xm
        alpha[12].setCoef(-SEML.cs_em.c1, 0);   //alpha13 = -c1
        break;
    }

    default:
    {
        cout << "error in pm. Unknown model." << endl;
        return;
    }

    }


    //------------------------------------------
    // Getting DB from files
    //------------------------------------------
    gsl_matrix_complex *DB = gsl_matrix_complex_calloc(6, 6);
    if(SEML.model == Csts::QBCP || SEML.model == Csts::BCP)
    {
        //Reading an OFS from a text file
        string filename = F_COC+"DB";
        glsc_matrix_complex_read(DB, filename);
    }
    else //RTBP
    {
        cout << "pm. The use of the RTBP has been detected when initializing DB." << endl;
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, c2;
        c2 = SEML.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        //cout << "omega1 = " << om1 << endl;
        //cout << "omega2 = " << om2 << endl;
        //--------------------
        //Init DB
        //--------------------
        gsl_matrix_complex_set(DB, 0, 0, gslc_complex(+om1*I));
        gsl_matrix_complex_set(DB, 1, 1, gslc_complex(+la1+0.0*I));
        gsl_matrix_complex_set(DB, 2, 2, gslc_complex(+om2*I));
        gsl_matrix_complex_set(DB, 3, 3, gslc_complex(-om1*I));
        gsl_matrix_complex_set(DB, 4, 4, gslc_complex(-la1+0.0*I));
        gsl_matrix_complex_set(DB, 5, 5, gslc_complex(-om2*I));
    }

    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |

    //cout << "real(DB)" << endl;
    //gslc_matrix_complex_printf_real(DB);
    //cout << "imag(DB)" << endl;
    //gslc_matrix_complex_printf_imag(DB);

    //------------------------------------------
    //DF0 = DB in matrix<Ofsc> format
    //------------------------------------------
    Ofsc BUX(OFS_ORDER);
    matrix<Ofsc> DF0(6, 6);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    DF0.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    DF0.setCoef(BUX, 1, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    DF0.setCoef(BUX, 2, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    DF0.setCoef(BUX, 3, 3);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    DF0.setCoef(BUX, 4, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    DF0.setCoef(BUX, 5, 5);


    //------------------------------------------
    // Setting H, Hinv, Lambda, LamdaL, LambdaN from DB
    // A this step, the central manifold is selected.
    //
    // TO-DO: something more generic to be able
    // to compute the central-stable and central-unstable manifold
    //
    //------------------------------------------
    gsl_matrix_complex *H     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Hinv  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *La    = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *LaL   = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex *LaN   = gsl_matrix_complex_calloc(2, 2);
    //H = (L N)
    gsl_matrix_complex *L    = gsl_matrix_complex_calloc(6, 4);
    gsl_matrix_complex *N    = gsl_matrix_complex_calloc(6, 2);

    //------------------------------------------
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L N)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //------------------------------------------
    gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
    gsl_matrix_complex_set(H, 1, 4, gsl_matrix_complex_get(DB, 1, 1));
    gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
    gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
    gsl_matrix_complex_set(H, 4, 5, gsl_matrix_complex_get(DB, 4, 4));
    gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

    //cout << "real(H)" << endl;
    //gslc_matrix_complex_printf_real(H);
    //cout << "imag(H)" << endl;
    //gslc_matrix_complex_printf_imag(H);


    //------------------------------------------
    //H2 = H in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> H2(6, 6);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    H2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    H2.setCoef(BUX, 1, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    H2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    H2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    H2.setCoef(BUX, 4, 5);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    H2.setCoef(BUX, 5, 3);

    //------------------------------------------
    //    |  iw1 0    0    0    |
    //    |  0   0    0    0    |
    //    |  0   iw3  0    0    |
    //L = |  0   0   -iw1  0    |
    //    |  0   0    0    0    |
    //    |  0   0    0   -iw3  |
    //------------------------------------------
    gsl_matrix_complex_view H11 = gsl_matrix_complex_submatrix (H , 0 , 0 , 6 , 4 );
    gsl_matrix_complex_memcpy(L, &H11.matrix);

    //cout << "imag(L)" << endl;
    //gslc_matrix_complex_printf_imag(L);
    matrix<Ofsc> L2(6, 4);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    L2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    L2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    L2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    L2.setCoef(BUX, 5, 3);


    //------------------------------------------
    //    |  0   0  |
    //    |  w2  0  |
    //    |  0   0  |
    //N = |  0   0  |
    //    |  0  -w2 |
    //    |  0   0  |
    //------------------------------------------
    H11 = gsl_matrix_complex_submatrix (H , 0 , 4 , 6 , 2 );
    gsl_matrix_complex_memcpy(N, &H11.matrix);


    //------------------------------------------
    //       |1/iw1  0    0      0     0     0    |
    //       | 0     0   1/iw3   0     0     0    |
    //       | 0     0    0    -1/iw1  0     0    |
    //Hinv = | 0     0    0      0     0  -1/iw3  |
    //       | 0    1/w2  0      0     0     0    |
    //       | 0     0    0      0   -1/w2   0    |
    //------------------------------------------
    gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
    gsl_matrix_complex_set(Hinv, 4, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
    gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
    gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
    gsl_matrix_complex_set(Hinv, 5, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
    gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

    //------------------------------------------
    //Hinv2 = Hinv in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> Hinv2(6, 6);
    BUX.zero();
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    Hinv2.setCoef(BUX, 0, 0);

    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    Hinv2.setCoef(BUX, 4, 1);

    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    Hinv2.setCoef(BUX, 1, 2);

    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    Hinv2.setCoef(BUX, 2, 3);

    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    Hinv2.setCoef(BUX, 5, 4);

    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    Hinv2.setCoef(BUX, 3, 5);

    //------------------------------------------
    //                     | LaL  T  |
    //Lambda = Hinv*DB*H = | 0  LaN  | with T = 0 in this context
    //
    //------------------------------------------
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(6, 6);
    //AUX = DB*H
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , H , zero_c , AUX );
    //B = Hinv*AUX
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Hinv , AUX , zero_c , La);

    //------------------------------------------
    gsl_matrix_complex_view La11 = gsl_matrix_complex_submatrix (La , 0 , 0 , 4 , 4 );   //S11
    gsl_matrix_complex_memcpy(LaL, &La11.matrix);

    gsl_matrix_complex_view La12 = gsl_matrix_complex_submatrix (La , 4 , 4 , 2 , 2 );   //S12
    gsl_matrix_complex_memcpy(LaN, &La12.matrix);

    //cout << "real(La)" << endl;
    //gslc_matrix_complex_printf_real(La);
    //cout << "imag(La)" << endl;
    //gslc_matrix_complex_printf_imag(La);



    //------------------------------------------
    // Initialisation of the TFC manifold ("W hat")
    //------------------------------------------
    cdouble db;
    gsl_complex dbc;
    //vector of size 6
    vector<Oftsc> Wh(6);
    //The order 0 is equal to zero, the order 1 is equal to L, which is diagonal
    //----------------------
    //     |   iw1*s1 |
    //     |     0    |
    //     |   iw3*s2 |
    //Wh = |  -iw1*s3 |
    //     |     0    |
    //     |  -iw3*s4 |
    //----------------------
    // Wh[0] = L11*s1 = L[0][0]*s[0]
    dbc = gsl_matrix_complex_get(L, 0, 0);
    db = gslc_complex(dbc);
    Wh[0].setCoef0(db, 1, 0);
    // Wh[1] = 0
    // Wh[2] = L32*s2 = L[2][1]*s[1]
    dbc = gsl_matrix_complex_get(L, 2, 1);
    db = gslc_complex(dbc);
    Wh[2].setCoef0(db, 1, 1);
    // Wh[3] = L43*s3 = L[3][2]*s[2]
    dbc = gsl_matrix_complex_get(L, 3, 2);
    db = gslc_complex(dbc);
    Wh[3].setCoef0(db, 1, 2);
    // Wh[4] = 0
    // Wh[5] = L64*s4 = L[5][3]*s[3]
    dbc = gsl_matrix_complex_get(L, 5, 3);
    db = gslc_complex(dbc);
    Wh[5].setCoef0(db, 1, 3);

    //------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------
    vector<Oftsc> W(6);   //W = TFC(Wh) (in NC coordinates)
    vector<Oftsc> Wt(6);  //Wt = W - V  (in NC coordinates)

    //Applying the COC (with and without zero order) @order 0 and 1
    for(int i = 0; i < 2; i++)
    {
        applyCOC(PC,  V,  Wh, W,  i, 1);   //order zero is added
        applyCOC(PC,  V,  Wh, Wt, i, 0);   //order zero is NOT added
    }

    //------------------------------------------
    // Differential of Wh (6*4)
    //------------------------------------------
    //Contribution of order < k to k
    matrix<Oftsc> DWh(6,4);  //TFC
    //matrix<Oftsc> DW(6,4);   //NC
    //Contribution of order <=k to k
    matrix<Oftsc> C_DW(6,4);   //NC

    //------------------------------------------
    // Initialisation of the complete VF
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> FWh(6);  //TFC
    vector<Oftsc> FW(6);   //NC
    vector<Oftsc> FW2(6);  //Temp
    vector<Oftsc> FW3(6);  //Temp
    //Contribution of order <=k to k
    vector<Oftsc> C_FW(6);   //NC

    //------------------------------------------
    // Initialisation of the reduced VF
    //------------------------------------------
    vector<Oftsc> fh(4);     //TFC
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        fh[i].setCoef0(db, 1, i);
    }
    //LaLs = fh @ order 0 and 1
    vector<Oftsc> LaLs(4);
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        //db = 1.0;
        LaLs[i].setCoef0(db, 1, i);
    }

    //------------------------------------------
    // Initialisation of DW times fh
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> DWfh(6);  //TFC
    //vector<Oftsc> DWf(6);   //NC
    //Contribution of order <=k to k
    vector<Oftsc> C_DWf(6);   //NC

    //------------------------------------------
    // Initialisation of dot(W)
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> Whdot(6);      //TFC
    vector<Oftsc> Wdot(6);       //NC
    //Contribution of order <=k to k
    vector<Oftsc> C_Wdot(6);     //NC


    //------------------------------------------
    // Initialisation of the Legendre objects
    //------------------------------------------
    //Potential
    vector<Oftsc> En(OFTS_ORDER+2);     //Earth
    vector<Oftsc> Sn(OFTS_ORDER+2);     //Moon
    vector<Oftsc> Mn(OFTS_ORDER+2);     //Sun
    vector<Oftsc> Un(OFTS_ORDER+2);     //Earth+Moon+Sun
    //Potential derivatives
    vector<Oftsc> Enx(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Eny(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Enz(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Mnx(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mny(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mnz(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Snx(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Sny(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Snz(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Unx(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Uny(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Unz(OFTS_ORDER+2);    //Earth+Moon+Sun

    //------------------------------------------
    // Initialisation of additionnal OFTS objects (see comments)
    // Function of Wh ("W hat")
    // Updated in updatePotentialExpansion routine
    //------------------------------------------
    Oftsc Rho2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);   //=xt^2+yt^2+zt^2
    Oftsc xse(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xe*xt+ye*yt
    Oftsc xsm(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xm*xt+ym*yt
    Oftsc xss(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xs*xt+ys*yt

    //------------------------------------------
    // Initialisation of various objects
    // Used to solve the pm equations
    //------------------------------------------
    vector<Oftsc> E(6);    //[DWh<m fh<m]m-[Fh(Wh<m)]m
    vector<Oftsc> eta(6);  //eta = Hinv*E
    vector<Oftsc> xi(6);   //xi = Hinv*Wh at order >=2
    cdouble lap;       //contains LaN
    cdouble lar;       //contains LaL*r for each monomial
    cdouble aux;       //temp obj
    cdouble bux;       //temp obj
    int kv[REDUCED_NV];       //contains the indices for each monomial


    //------------------------------------------
    // For testing
    //------------------------------------------
    //Copy for consistency testing
    matrix<Oftsc> T_DWh(6,4);  //TFC
    vector<Oftsc> T_FWh(6);    //TFC
    vector<Oftsc> T_DWfh(6);   //TFC
    vector<Oftsc> T_Whdot(6);  //TFC
    vector<Oftsc> T_FW2(6);    //Temp
    vector<Oftsc> T_FW3(6);    //Temp
    //Error in the invariance equation (global object)
    vector<Oftsc> IE(6);  //NC
    vector<Oftsc> IEh(6); //TFC
    //Error in the invariance equation (ots object)
    vector<Otsc>  IE_ots(6);  //NC
    vector<Otsc>  IEh_ots(6); //TFC
    //Error in the invariance equation (eta variable)
    vector<Oftsc> xi1(6);
    matrix<Oftsc> xi2(6,4);
    vector<Oftsc> xi3(6);
    vector<Oftsc> xi4(6);
    vector<Oftsc> xi5(6);
    vector<Otsc>  xi5_ots(6);
    //Error in the invariance equation (Wh variable)
    vector<Oftsc> Wh1(6);
    matrix<Oftsc> Wh2(6,4);
    vector<Oftsc> Wh3(6);
    vector<Oftsc> Wh4(6);
    vector<Oftsc> Wh5(6);
    vector<Oftsc> Wh6(6);
    vector<Otsc>  Wh6_ots(6);
    //Error on the vector field
    vector<Oftsc> D1(6);
    vector<Otsc>  D1_ots(6);
    //Error on the DW*f product
    vector<Oftsc> D2(6);
    vector<Otsc>  D2_ots(6);

    cout << "  pm. Initialization of the recurrence.  " << endl;

    //-------------------------------------------------------------
    // Initialization of the recurrence
    //-------------------------------------------------------------
    //The potential is updated: Un, Unx, Uny, Unz
    // @order 0, 1 and 2
    //Wt is used to update the Potential
    updatePotentialExpansion(Wt, Rho2, En,
                             Mn, Sn, Un,
                             Enx, Eny, Enz,
                             Mnx, Mny, Mnz,
                             Snx, Sny, Snz,
                             Unx, Uny, Unz,
                             Xe, Xm, Xs,
                             ILe, ILm, ILs,
                             xse, xsm, xss,
                             SEML, 2);

    cout << "  pm. End of initialisation.  " << endl;

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ order 0 and 1
    //------------------------------------------
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    for(int k = 0; k < 2 ; k++)
    {
        //vector field
        //-----------------
        //applyVF(alpha, W, FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);
        //Copy for consistency testing
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);

        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, k);

        //dot(W)
        //-----------------
        for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], SEML.us.n, k);       //in NC
        for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], SEML.us.n, k);     //in TFC
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, k);     //in NC, complete contribution


        //----------------------
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, k);   //in NC
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, k);   //in NC, complete contribution

        // DWh times f @order k
        //----------------------
        //smvprod_t(DW,  fh, DWf,  k);     //in NC
        smvprod_t(DWh, fh, DWfh, k);     //in TFC
        smvprod_t(C_DW,  fh, C_DWf,  k); //in NC, complete contribution

        // Invariance equation test
        //-------------------------
        //FWh = COCinv(FW)
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, T_FWh, T_FW2, T_FW3, k);
        //dot(W)
        for(int i = 0; i < 6; i++) T_Whdot[i].dot(Wh[i], SEML.us.n, k);
        // Differential of Wh
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) T_DWh.der(Wh[i], j+1, i, j, k);
        // DWh times f @order k
        smvprod_t(T_DWh, fh, T_DWfh, k);

        //IE
        for(int i = 0; i<6; i++)
        {
            IE[i].ofts_smult_u(C_FW[i],     1.0+0.0*I, k);
            IE[i].ofts_smult_u(C_DWf[i],   -1.0+0.0*I, k);
            IE[i].ofts_smult_u(C_Wdot[i],  -1.0+0.0*I, k);
        }

        //IEh
        for(int i = 0; i<6; i++)
        {
            IEh[i].ofts_smult_u(T_FWh[i],     1.0+0.0*I, k);
            IEh[i].ofts_smult_u(T_DWfh[i],   -1.0+0.0*I, k);
            IEh[i].ofts_smult_u(T_Whdot[i],  -1.0+0.0*I, k);
        }

    }


    cout << "  pm. End of order 0/1.  " << endl;

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ orders >= 2
    //------------------------------------------
    for(int m = 2; m <= OFTS_ORDER; m++)
    {
        cout << "pm. start of order " << m << endl;
        tic();
        //------------------------------------------
        //The potential is updated: Un, Unx, Uny, Unz
        //Wt is used to update the Potential
        //What is updated: contributions of order < m to order m
        //------------------------------------------
        if(m > 2) updatePotentialExpansion(Wt, Rho2, En,
                                               Mn, Sn, Un,
                                               Enx, Eny, Enz,
                                               Mnx, Mny, Mnz,
                                               Snx, Sny, Snz,
                                               Unx, Uny, Unz,
                                               Xe, Xm, Xs,
                                               ILe, ILm, ILs,
                                               xse, xsm, xss,
                                               SEML, m);

        //------------------------------------------
        // Applying the VF: [F(W<m)]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);

        //------------------------------------------
        //FWh = COCinv(FW)
        //------------------------------------------
        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, m);

        //cout << "pm. end of applyInvDotCOC." << endl;

        //------------------------------------------
        // dot(W<m)
        // It is NOT updated here: the terms of order < m
        // dot not influence dot(W) at order m
        //------------------------------------------
        //for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], qbcp.n, m);
        //for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], qbcp.n, m);
        //cout << "pm. end of dot." << " in " << toc() << " s. " << endl;

        //------------------------------------------
        // Differential of Wh
        // It is NOT updated here, but at the end of the loop.
        // Indeed, Wh & W contain no order m at this point!
        //------------------------------------------

        //------------------------------------------
        // [DW<m times f<m]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        smvprod_t(DWh, fh, DWfh, m);
        //smvprod_t(DW, fh, DWf, m);

        //------------------------------------------
        // We have [Fh(Wh<m)]m and [DWh<m fh<m]m
        // Update of eta at order m
        //------------------------------------------
        // [E]m = [DWh<m fh<m]m - [Fh(Wh<m)]m
        for(int i = 0 ; i< 6; i++) E[i].ofts_fsum_u(FWh[i], -1.0+0.0*I, DWfh[i], 1.0+0.0*I, m);
        //[eta]m = [Hinv*E]m;
        smvprod_u(Hinv2, E, eta, m);

        //------------------------------------------
        // Equation in the normal space:
        // Last two components of eta
        //------------------------------------------
        //Sum on the normal components of eta (last n-d components)
        for(int p = REDUCED_NV; p<= Csts::NV-1; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
            //Reset kv
            kv[0] = m;
            for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Reset aux
                aux = 0.0+0.0*I;
                //Update aux = laL.kv
                for(int ii = 0; ii < REDUCED_NV; ii++)
                {
                    lar = gslc_complex(gsl_matrix_complex_get(LaL, ii, ii));
                    aux += lar*(kv[ii]+0.0*I);
                }

                //Update xi
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    bux = 0.0*I-eta[p].getCoef(m, k)->ofs_getCoef(j)/(SEML.us.n*I*j - (lap - aux));
                    xi[p].getCoef(m, k)->setCoef(bux, j);
                }
                //Update kv
                if(k < FTDA::nmon(REDUCED_NV, m)-1)  FTDA::prxkt(kv, REDUCED_NV);
            }



            //---------------------------------
            //Normal cohomological equations
            //Consistency Test
            //---------------------------------
            //xi1 = lap*xi[p] at order m
            xi1[p].ofts_smult_u(xi[p], (cdouble) lap, m);
            //xi2(p,j) = dxi[p]/dx[j] at order m
            for(int j = 0 ; j < REDUCED_NV; j++) xi2.der(xi[p], j+1, p, j, m);
            //xi3 = xi2*LaLs
            smvprod_t(xi2, LaLs, xi3, m);
            //xi4 = dot(xi)
            xi4[p].dot(xi[p], SEML.us.n, m);
            //xi5 = xi1 - xi3 - xi4 - eta (needs to be zero!)
            xi5[p].ofts_smult_u(xi1[p],   1.0+0.0*I, m);
            xi5[p].ofts_smult_u(xi3[p],  -1.0+0.0*I, m);
            xi5[p].ofts_smult_u(xi4[p],  -1.0+0.0*I, m);
            xi5[p].ofts_smult_u(eta[p],  -1.0+0.0*I, m);
        }

        //------------------------------------------
        // Equation in the tangential space:
        // Graph style is used by default
        // As a consequence:
        // - xi^p_m = 0
        // - fh^p_m = -eta^p_m
        //------------------------------------------
        //Sum on the tangential components of eta (first d components)
        for(int p = 0; p< REDUCED_NV; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));

            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Update fh
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    bux = eta[p].getCoef(m, k)->ofs_getCoef(j);
                    fh[p].getCoef(m, k)->setCoef(+0.0*I-bux, j);
                }
            }

            //---------------------------------
            //Tangential cohomological equations
            //Consistency Test
            //---------------------------------
            //xi1 = lap*xi[p] at order m
            xi1[p].ofts_smult_u(xi[p], (cdouble) lap, m);
            //xi2(p,j) = dxi[p]/dx[j] at order m
            for(int j = 0 ; j < REDUCED_NV; j++) xi2.der(xi[p], j+1, p, j, m);
            //xi3 = xi2*LaLs
            smvprod_t(xi2, LaLs, xi3, m);
            //xi4 = dot(xi)
            xi4[p].dot(xi[p], SEML.us.n, m);
            //xi5 = xi1 - xi3 - xi4 - eta - fh (needs to be zero!)
            xi5[p].ofts_smult_u(xi1[p],   1.0+0.0*I, m);
            xi5[p].ofts_smult_u(xi3[p],  -1.0+0.0*I, m);
            xi5[p].ofts_smult_u(xi4[p],  -1.0+0.0*I, m);
            xi5[p].ofts_smult_u(eta[p],  -1.0+0.0*I, m);
            xi5[p].ofts_smult_u(fh[p],   -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) xi5_ots[p].setCoef(xi5[p].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);
        }

        //------------------------------------------
        // Wh = H2*xi @ order m
        // Recall that H2 = H in matrix format
        //------------------------------------------
        smvprod_u(H2, xi, Wh, m);

        //---------------------------------
        //Tangential & Tangential cohomological equations
        //Consistency Test in Wh variables
        //---------------------------------
        //Wh1 = DF(0)Wh @ order m
        smvprod_u(DF0, Wh, Wh1, m);
        //Wh2(p,j) = dWh[p]/dx[j] at order m
        for(int p = 0; p < Csts::NV; p++) for(int j = 0 ; j < REDUCED_NV; j++) Wh2.der(Wh[p], j+1, p, j, m);
        //Wh3 = Wh2*LaLs
        smvprod_t(Wh2, LaLs, Wh3, m);
        //Wh4 = L*fh
        smvprod_u(L2, fh, Wh4, m);
        //Wh5 = dot(Wh)
        for(int p = 0; p < Csts::NV; p++) Wh5[p].dot(Wh[p], SEML.us.n, m);
        //Wh6 = Wh1 - Wh3 - - Wh4 - Wh5 - E
        for(int p = 0; p < Csts::NV; p++)
        {
            Wh6[p].ofts_smult_u(Wh1[p],   1.0+0.0*I, m);
            Wh6[p].ofts_smult_u(Wh3[p],  -1.0+0.0*I, m);
            Wh6[p].ofts_smult_u(Wh4[p],  -1.0+0.0*I, m);
            Wh6[p].ofts_smult_u(Wh5[p],  -1.0+0.0*I, m);
            Wh6[p].ofts_smult_u(E[p],    -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) Wh6_ots[p].setCoef(Wh6[p].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);
        }


        //------------------------------------------
        // Update, W, Wt and W2 @ order m
        //------------------------------------------
        applyCOC(PC,  V,  Wh, W,  m, 1);     //order zero is added
        applyCOC(PC,  V,  Wh, Wt, m, 0);     //order zero is NOT added

        //------------------------------------------
        // Update the differential to take into account the new terms @ order m
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, m);
//        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DW.der(W[i], j+1, i, j, m);

        //------------------------------------------------------------------------------------
        //
        // Note: after this point, only necessary for testing, can be discarded for speed
        //
        //------------------------------------------------------------------------------------
        initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, SEML, m);

        //------------------------------------------
        // update the VF
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, m);

        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ, PCdot, Vdot, Wh, C_FW, T_FWh, T_FW2, T_FW3, m);

        //------------------------------------------
        // dot(W)
        //------------------------------------------
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, m);
        for(int i = 0; i < 6; i++) T_Whdot[i].dot(Wh[i], SEML.us.n, m);

        //------------------------------------------
        // Differential of Wh
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) T_DWh.der(Wh[i], j+1, i, j, m);
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, m);

        //------------------------------------------
        // DWh times f @order m
        //------------------------------------------
        smvprod_t(T_DWh, fh, T_DWfh, m);
        smvprod_t(C_DW,  fh, C_DWf,  m);

        cout << "pm. end of update testing." << " in " << toc() << " s. "<< endl;

        //------------------------------------------
        // Invariance equation test
        //------------------------------------------

        for(int i = 0; i<6; i++)
        {
            IE[i].ofts_smult_u(C_FW[i],     1.0+0.0*I, m);
            IE[i].ofts_smult_u(C_DWf[i],   -1.0+0.0*I, m);
            IE[i].ofts_smult_u(C_Wdot[i],  -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) IE_ots[i].setCoef(IE[i].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);
        }

        //IEh
        for(int i = 0; i<6; i++)
        {
            IEh[i].ofts_smult_u(T_FWh[i],     1.0+0.0*I, m);
            IEh[i].ofts_smult_u(T_DWfh[i],   -1.0+0.0*I, m);
            IEh[i].ofts_smult_u(T_Whdot[i],  -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) IEh_ots[i].setCoef(IEh[i].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);
        }

        //D1 & D2
        for(int i = 0; i<6; i++)
        {
            D1[i].ofts_smult_u(T_FWh[i],   1.0+0.0*I, m);
            D1[i].ofts_smult_u(Wh1[i],    -1.0+0.0*I, m);
            D1[i].ofts_smult_u(FWh[i],    -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) D1_ots[i].setCoef(D1[i].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);

            D2[i].ofts_smult_u(T_DWfh[i],   1.0+0.0*I, m);
            D2[i].ofts_smult_u(Wh3[i],     -1.0+0.0*I, m);
            D2[i].ofts_smult_u(Wh4[i],     -1.0+0.0*I, m);
            D2[i].ofts_smult_u(DWfh[i],    -1.0+0.0*I, m);
            for(int jj = 0; jj < FTDA::nmon(REDUCED_NV, m); jj++) D2_ots[i].setCoef(D2[i].getCoef(m, jj)->getCoefMaxNorm()+0.0*I, m, jj);
        }


        cout << "pm. end of invariance equation test." << " in " << toc() << " s. "<< endl;

        cout << "  pm. End of order " << m << " in " << toc() << " s. " << endl;
    }

    //------------------------------------------
    // Printing
    //------------------------------------------
    tic();
    //Invariance Equation test
    vector_fprinf(IE_ots,  F_GS+"Test/IE_ots");
    vector_fprinf(IEh_ots, F_GS+"Test/IEh_ots");
    //print W & FW
    vector_fprinf(W,     F_GS+"W/W");
    vector_fprinf(Wt,     F_GS+"W/Wt");
    vector_fprinf(Wh,    F_GS+"W/Wh");
    vector_fprinf(FW,    F_GS+"FW/FW");
    vector_fprinf(C_FW,  F_GS+"FW/C_FW");
    vector_fprinf(FWh,   F_GS+"FW/FWh");
    //fprintf Whdot & Wdot
    vector_fprinf(C_Wdot,  F_GS+"Wdot/C_Wdot");
    vector_fprinf(T_Whdot, F_GS+"Wdot/T_Whdot");
    //Reduced vector field
    vector_fprinf(fh, F_GS+"rvf/fh");
    //Secondary testing
    vector_fprinf(xi5_ots, F_GS+"Test/IE_xi");
    vector_fprinf(Wh6_ots, F_GS+"Test/IE_Wh");
    vector_fprinf(D1_ots,  F_GS+"Test/D1");
    vector_fprinf(D2_ots,  F_GS+"Test/D2");
    vector_fprinf(D1,      F_GS+"Test/D1v");
    //Potential
    vector_fprinf(Unx, F_GS+"Potential/Unx");
    vector_fprinf(Uny, F_GS+"Potential/Uny");
    vector_fprinf(Unz, F_GS+"Potential/Unz");
    vector_fprinf(En, F_GS+"Potential/En");
    vector_fprinf(Mn, F_GS+"Potential/Mn");
    vector_fprinf(Sn, F_GS+"Potential/Sn");
    vector_fprinf(Un, F_GS+"Potential/Un");
    //DW*vf
    vector_fprinf(T_DWfh, F_GS+"DWf/T_DWfh");
    vector_fprinf(C_DWf,  F_GS+"DWf/C_DWf");
    cout << "  pm. End of printing" <<  " in " << toc() << " s. " << endl;

}

/**
 *  \brief Compute the parameterization of the central manifold of the dynamical equivalent of the Earth-Moon libration points L1,2, up to a given order (normal form style).
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_GS...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_GS):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 **/
void pm_normalform(int OutputEachOrder, int Output)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Application of the parameterization method     " << endl;
    cout << "  on the non-autonomous vector field of the QBCP  " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(4);
    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_NF    = SEML.cs.F_NF;
    string F_COC   = SEML.cs.F_COC;
    string F_COEF  = SEML.cs.F_COEF;

    cout << "F_COEF = " << F_COEF << endl;
    cout << "F_COC = "  << F_COC << endl;
    cout << "F_NF = "   << F_NF << endl;


    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    matrix<Ofsc> P(6,6);      //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> Q(6,6);      //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    vector<Ofsc> V(6);        //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> Xe(3);       //The vectors Xe, Xm, Xs, modified position
    vector<Ofsc> Xm(3);       //of the primaries used to compute the potential:
    vector<Ofsc> Xs(3);       //Xf = True Xf - V, for f = e, m, s
    Ofsc ILe(OFS_ORDER);      //The Ofs ILe, ILm and ILs, modified inverse orbit radii
    Ofsc ILm(OFS_ORDER);      //of the primaries used to compute the potential
    Ofsc ILs(OFS_ORDER);      //ILf = 1.0/||Xf||, for f = e, m, s
    matrix<Ofsc> PC(6,6);     //PC=P(theta)*C
    matrix<Ofsc> PCdot(6,6);  //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> CQ(6,6);     //CQ=Cinv*Pinv=inv(PC)
    vector<Ofsc> Vdot(6);     //Vdot=dot(V)

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Initialization of the alphas
    //------------------------------------------
    int noc = SEML.numberOfCoefs;
    vector<Ofsc> alpha(noc);
    //Switch on the model to initialize the alphas
    switch(SEML.model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    case Csts::ERTBP:
    {
        //Read from txt files
        ifstream readStream;
        string ss1;
        for(int i = 0; i < noc; i++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
            readOFS_txt(alpha[i], F_COEF+"alpha"+ss1+"_fft");
        }

        break;
    }

    case Csts::CRTBP:
    {
        cout << "pm. The use of the RTBP has been detected when initializing the alphas." << endl;
        //All coeffs to zero
        for(int j = 0; j < noc;  j++) alpha[j].setCoef(0.0, 0);
        //Non-null coeffs
        alpha[0].setCoef(1.0, 0);               //alpha1  = 1.0
        alpha[2].setCoef(1.0, 0);               //alpha3  = 1.0
        alpha[5].setCoef(1.0, 0);               //alpha6  = 1.0
        alpha[8].setCoef(SEML.us_em.mu_EM, 0);  //alpha9  = Xe
        alpha[10].setCoef(SEML.us_em.mu_EM, 0); //alpha11 = Xm
        alpha[12].setCoef(-SEML.cs_em.c1, 0);   //alpha13 = -c1
        break;
    }

    default:
    {
        cout << "error in pm. Unknown model." << endl;
        return;
    }

    }



    //------------------------------------------
    // Getting DB from files
    //------------------------------------------
    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |
    gsl_matrix_complex *DB = gsl_matrix_complex_calloc(6, 6);
    if(SEML.model == Csts::QBCP || SEML.model == Csts::BCP)
    {
        //Reading an OFS from a text file
        string filename = F_COC+"DB";
        glsc_matrix_complex_read(DB, filename);
    }
    else //RTBP
    {
        cout << "pm. The use of the RTBP has been detected when initializing DB." << endl;
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, c2;
        c2 = SEML.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        //cout << "omega1 = " << om1 << endl;
        //cout << "omega2 = " << om2 << endl;
        //--------------------
        //Init DB
        //--------------------
        gsl_matrix_complex_set(DB, 0, 0, gslc_complex(+om1*I));
        gsl_matrix_complex_set(DB, 1, 1, gslc_complex(+la1+0.0*I));
        gsl_matrix_complex_set(DB, 2, 2, gslc_complex(+om2*I));
        gsl_matrix_complex_set(DB, 3, 3, gslc_complex(-om1*I));
        gsl_matrix_complex_set(DB, 4, 4, gslc_complex(-la1+0.0*I));
        gsl_matrix_complex_set(DB, 5, 5, gslc_complex(-om2*I));
    }
    cout << setprecision(15) << endl;
    //------------------------------------------
    //DF0 = DB in matrix<Ofsc> format
    //------------------------------------------
    Ofsc BUX(OFS_ORDER);
    matrix<Ofsc> DF0(6, 6);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    DF0.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    DF0.setCoef(BUX, 1, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    DF0.setCoef(BUX, 2, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    DF0.setCoef(BUX, 3, 3);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    DF0.setCoef(BUX, 4, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    DF0.setCoef(BUX, 5, 5);


    //------------------------------------------
    // Setting H, Hinv, Lambda, LamdaL, LambdaN from DB
    // A this step, the central manifold is selected.
    //
    // TO-DO: something more generic to be able
    // to compute the central-stable and central-unstable manifold
    //
    //------------------------------------------
    gsl_matrix_complex *H     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Hinv  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *La    = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *LaL   = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex *LaN   = gsl_matrix_complex_calloc(2, 2);
    //H = (L N)
    gsl_matrix_complex *L    = gsl_matrix_complex_calloc(6, 4);
    gsl_matrix_complex *N    = gsl_matrix_complex_calloc(6, 2);

    //------------------------------------------
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L N)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //------------------------------------------
    gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
    gsl_matrix_complex_set(H, 1, 4, gsl_matrix_complex_get(DB, 1, 1));
    gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
    gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
    gsl_matrix_complex_set(H, 4, 5, gsl_matrix_complex_get(DB, 4, 4));
    gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

    //------------------------------------------
    //H2 = H in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> H2(6, 6);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    H2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    H2.setCoef(BUX, 1, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    H2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    H2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    H2.setCoef(BUX, 4, 5);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    H2.setCoef(BUX, 5, 3);

    //------------------------------------------
    //    |  iw1 0    0    0    |
    //    |  0   0    0    0    |
    //    |  0   iw3  0    0    |
    //L = |  0   0   -iw1  0    |
    //    |  0   0    0    0    |
    //    |  0   0    0   -iw3  |
    //------------------------------------------
    gsl_matrix_complex_view H11 = gsl_matrix_complex_submatrix (H , 0 , 0 , 6 , 4 );
    gsl_matrix_complex_memcpy(L, &H11.matrix);

    //matrix version of L
    matrix<Ofsc> L2(6, 4);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    L2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    L2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    L2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    L2.setCoef(BUX, 5, 3);


    //------------------------------------------
    //    |  0   0  |
    //    |  w2  0  |
    //    |  0   0  |
    //N = |  0   0  |
    //    |  0  -w2 |
    //    |  0   0  |
    //------------------------------------------
    H11 = gsl_matrix_complex_submatrix (H , 0 , 4 , 6 , 2 );
    gsl_matrix_complex_memcpy(N, &H11.matrix);


    //------------------------------------------
    //       |1/iw1  0    0      0     0     0    |
    //       | 0     0   1/iw3   0     0     0    |
    //       | 0     0    0    -1/iw1  0     0    |
    //Hinv = | 0     0    0      0     0  -1/iw3  |
    //       | 0    1/w2  0      0     0     0    |
    //       | 0     0    0      0   -1/w2   0    |
    //------------------------------------------
    gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
    gsl_matrix_complex_set(Hinv, 4, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
    gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
    gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
    gsl_matrix_complex_set(Hinv, 5, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
    gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

    //------------------------------------------
    //Hinv2 = Hinv in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> Hinv2(6, 6);
    BUX.zero();
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    Hinv2.setCoef(BUX, 0, 0);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    Hinv2.setCoef(BUX, 4, 1);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    Hinv2.setCoef(BUX, 1, 2);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    Hinv2.setCoef(BUX, 2, 3);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    Hinv2.setCoef(BUX, 5, 4);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    Hinv2.setCoef(BUX, 3, 5);

    //------------------------------------------
    //                     | LaL  T  |
    //Lambda = Hinv*DB*H = | 0  LaN  | with T = 0 in this context
    //
    //------------------------------------------
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(6, 6);
    //AUX = DB*H
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , H , zero_c , AUX );
    //B = Hinv*AUX
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Hinv , AUX , zero_c , La);
    //------------------------------------------
    gsl_matrix_complex_view La11 = gsl_matrix_complex_submatrix (La , 0 , 0 , 4 , 4 );   //S11
    gsl_matrix_complex_memcpy(LaL, &La11.matrix);
    gsl_matrix_complex_view La12 = gsl_matrix_complex_submatrix (La , 4 , 4 , 2 , 2 );   //S12
    gsl_matrix_complex_memcpy(LaN, &La12.matrix);

    //------------------------------------------
    // Initialisation of the TFC manifold ("W hat")
    //------------------------------------------
    cdouble db;
    gsl_complex dbc;
    //vector of size 6
    vector<Oftsc> Wh(6);
    //The order 0 is equal to zero, the order 1 is equal to L, which is diagonal
    //----------------------
    //     |   iw1*s1 |
    //     |     0    |
    //     |   iw3*s2 |
    //Wh = |  -iw1*s3 |
    //     |     0    |
    //     |  -iw3*s4 |
    //----------------------
    // Wh[0] = L11*s1 = L[0][0]*s[0]
    dbc = gsl_matrix_complex_get(L, 0, 0);
    db = gslc_complex(dbc);
    Wh[0].setCoef0(db, 1, 0);
    // Wh[1] = 0
    // Wh[2] = L32*s2 = L[2][1]*s[1]
    dbc = gsl_matrix_complex_get(L, 2, 1);
    db = gslc_complex(dbc);
    Wh[2].setCoef0(db, 1, 1);
    // Wh[3] = L43*s3 = L[3][2]*s[2]
    dbc = gsl_matrix_complex_get(L, 3, 2);
    db = gslc_complex(dbc);
    Wh[3].setCoef0(db, 1, 2);
    // Wh[4] = 0
    // Wh[5] = L64*s4 = L[5][3]*s[3]
    dbc = gsl_matrix_complex_get(L, 5, 3);
    db = gslc_complex(dbc);
    Wh[5].setCoef0(db, 1, 3);

    //------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------
    vector<Oftsc> W(6);   //W = TFC(Wh) (in NC coordinates)
    vector<Oftsc> Wt(6);  //Wt = W - V  (in NC coordinates)
    //Applying the COC (with and without zero order) @order 0 and 1
    for(int i = 0; i < 2; i++)
    {
        applyCOC(PC,  V,  Wh, W,  i, 1);   //order zero is added
        applyCOC(PC,  V,  Wh, Wt, i, 0);   //order zero is NOT added
    }

    //------------------------------------------
    // Differential of Wh (6*4)
    //------------------------------------------
    //Contribution of order < k to k
    matrix<Oftsc> DWh(6,4);  //TFC
    //Contribution of order <=k to k
    matrix<Oftsc> C_DW(6,4);   //NC

    //------------------------------------------
    // Initialisation of the complete VF
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> FWh(6);  //TFC
    vector<Oftsc> FW(6);   //NC
    vector<Oftsc> FW2(6);  //Temp
    vector<Oftsc> FW3(6);  //Temp
    //Contribution of order <=k to k
    vector<Oftsc> C_FW(6);   //NC

    //------------------------------------------
    // Initialisation of the reduced VF
    //------------------------------------------
    vector<Oftsc> fh(4);     //TFC
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        fh[i].setCoef0(db, 1, i);
    }
    //LaLs = fh @ order 0 and 1
    vector<Oftsc> LaLs(4);
    for(int i = 0; i < 4; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        //db = 1.0;
        LaLs[i].setCoef0(db, 1, i);
    }

    //------------------------------------------
    // Initialisation of DW times fh
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> DWfh(6);  //TFC
    //Contribution of order <=k to k
    vector<Oftsc> C_DWf(6);   //NC

    //------------------------------------------
    // Initialisation of dot(W)
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> Whdot(6);      //TFC
    vector<Oftsc> Wdot(6);       //NC
    //Contribution of order <=k to k
    vector<Oftsc> C_Wdot(6);     //NC


    //------------------------------------------
    // Initialisation of the Legendre objects
    //------------------------------------------
    //Potential
    vector<Oftsc> En(OFTS_ORDER+2);     //Earth
    vector<Oftsc> Sn(OFTS_ORDER+2);     //Moon
    vector<Oftsc> Mn(OFTS_ORDER+2);     //Sun
    vector<Oftsc> Un(OFTS_ORDER+2);     //Earth+Moon+Sun
    //Potential derivatives
    vector<Oftsc> Enx(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Eny(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Enz(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Mnx(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mny(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mnz(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Snx(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Sny(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Snz(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Unx(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Uny(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Unz(OFTS_ORDER+2);    //Earth+Moon+Sun

    //------------------------------------------
    // Initialisation of additionnal OFTS objects (see comments)
    // Function of Wh ("W hat")
    // Updated in updatePotentialExpansion routine
    //------------------------------------------
    Oftsc Rho2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);   //=xt^2+yt^2+zt^2
    Oftsc xse(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xe*xt+ye*yt
    Oftsc xsm(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xm*xt+ym*yt
    Oftsc xss(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xs*xt+ys*yt

    //------------------------------------------
    // Initialisation of various objects
    // Used to solve the pm equations
    //------------------------------------------
    vector<Oftsc> E(6);    //[DWh<m fh<m]m-[Fh(Wh<m)]m
    vector<Oftsc> eta(6);  //eta = Hinv*E
    vector<Oftsc> xi(6);   //xi = Hinv*Wh at order >=2
    cdouble lap;       //contains LaN
    cdouble lar;       //contains LaL*r for each monomial
    cdouble aux;       //temp obj
    cdouble bux;       //temp obj
    int kv[REDUCED_NV];       //contains the indices for each monomial
    cdouble divisor;          //divisor in the cohomological equations
    double threshold = 1e-2;  //threshold for small divisors = 1e-2


    //------------------------------------------
    // Initialisation of the small divisor container
    //------------------------------------------
    vector<Oftsc> smallDiv(6);


    cout << "  pm. Initialization of the recurrence.  " << endl;



    //-------------------------------------------------------------
    // Initialization of the recurrence
    //-------------------------------------------------------------
    //The potential is updated: Un, Unx, Uny, Unz
    // @order 0, 1 and 2
    //Wt is used to update the Potential
    updatePotentialExpansion(Wt, Rho2, En,
                             Mn, Sn, Un,
                             Enx, Eny, Enz,
                             Mnx, Mny, Mnz,
                             Snx, Sny, Snz,
                             Unx, Uny, Unz,
                             Xe, Xm, Xs,
                             ILe, ILm, ILs,
                             xse, xsm, xss,
                             SEML, 2);

    cout << "  pm. End of initialisation.  " << endl;

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ order 0 and 1
    //------------------------------------------
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    for(int k = 0; k < 2 ; k++)
    {
        //vector field
        //-----------------
        //applyVF(alpha, W, FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);
        //Copy for consistency testing
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);

        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, k);

        //dot(W)
        //-----------------
        for(int i = 0; i < 6; i++) Wdot[i].dot(W[i],   SEML.us.n, k);     //in NC
        for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], SEML.us.n, k);     //in TFC
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, k);     //in NC, complete contribution

        //----------------------
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, k);   //in NC
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, k);   //in NC, complete contribution

        // DWh times f @order k
        //----------------------
        //smvprod_t(DW,  fh, DWf,  k);     //in NC
        smvprod_t(DWh, fh, DWfh, k);       //in TFC
        smvprod_t(C_DW,  fh, C_DWf,  k);   //in NC, complete contribution


    }

    cout << "  pm. End of order 0/1.  " << endl;

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ orders >= 2
    //------------------------------------------
    for(int m = 2; m <= OFTS_ORDER; m++)
    {

        cout << "pm. start of order " << m << endl;
        tic();
        //------------------------------------------
        //The potential is updated: Un, Unx, Uny, Unz
        //Wt is used to update the Potential
        //What is updated: contributions of order < m to order m
        //------------------------------------------
        if(m > 2) updatePotentialExpansion(Wt, Rho2, En,
                                               Mn, Sn, Un,
                                               Enx, Eny, Enz,
                                               Mnx, Mny, Mnz,
                                               Snx, Sny, Snz,
                                               Unx, Uny, Unz,
                                               Xe, Xm, Xs,
                                               ILe, ILm, ILs,
                                               xse, xsm, xss,
                                               SEML, m);

        //------------------------------------------
        // Applying the VF: [F(W<m)]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);

        //------------------------------------------
        //FWh = COCinv(FW)
        //------------------------------------------
        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, m);

        //cout << "pm. end of applyInvDotCOC." << endl;

        //------------------------------------------
        // dot(W<m)
        // It is NOT updated here: the terms of order < m
        // dot not influence dot(W) at order m
        //------------------------------------------
        //for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], qbcp.n, m);
        //for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], qbcp.n, m);
        //cout << "pm. end of dot." << " in " << toc() << " s. " << endl;

        //------------------------------------------
        // Differential of Wh
        // It is NOT updated here, but at the end of the loop.
        // Indeed, Wh & W contain no order m at this point!
        //------------------------------------------

        //------------------------------------------
        // [DW<m times f<m]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        smvprod_t(DWh, fh, DWfh, m);
        //smvprod_t(DW, fh, DWf, m);


        //------------------------------------------
        // We have [Fh(Wh<m)]m and [DWh<m fh<m]m
        // Update of eta at order m
        //------------------------------------------
        // [E]m = [DWh<m fh<m]m - [Fh(Wh<m)]m
        for(int i = 0 ; i< 6; i++) E[i].ofts_fsum_u(FWh[i], -1.0+0.0*I, DWfh[i], 1.0+0.0*I, m);
        //[eta]m = [Hinv*E]m;
        smvprod_u(Hinv2, E, eta, m);



        //------------------------------------------
        // Equation in the normal space:
        // Last two components of eta
        //------------------------------------------
        //Sum on the normal components of eta (last n-d components)
        for(int p = REDUCED_NV; p<= Csts::NV-1; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
            //Reset kv
            kv[0] = m;
            for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Reset aux
                aux = 0.0+0.0*I;
                //Update aux = laL.kv
                for(int ii = 0; ii < REDUCED_NV; ii++)
                {
                    lar = gslc_complex(gsl_matrix_complex_get(LaL, ii, ii));
                    aux += lar*(kv[ii]+0.0*I);
                }

                //Update xi
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    divisor = SEML.us.n*I*j - (lap - aux);
                    bux = 0.0*I-eta[p].getCoef(m, k)->ofs_getCoef(j)/divisor;
                    xi[p].getCoef(m, k)->setCoef(bux, j);
                    //Small divisor
                    smallDiv[p].getCoef(m, k)->setCoef(divisor, j);
                }
                //Update kv
                if(k < FTDA::nmon(REDUCED_NV, m)-1)  FTDA::prxkt(kv, REDUCED_NV);
            }
        }

        //------------------------------------------
        // Equation in the tangent space:
        // First four components of eta
        //------------------------------------------
        //Sum on the tangential components of eta (first d components)
        for(int p = 0; p< REDUCED_NV; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
            //Reset kv
            kv[0] = m;
            for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;

            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Reset aux
                aux = 0.0+0.0*I;
                //Update aux = laL.kv
                for(int ii = 0; ii < REDUCED_NV; ii++)
                {
                    lar = gslc_complex(gsl_matrix_complex_get(LaL, ii, ii));
                    aux += lar*(kv[ii]+0.0*I);
                }

                //Update xi
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    //divisor = jnI - (lap - laL.kv)
                    divisor = SEML.us.n*I*j - (lap - aux);
                    //Small divisor
                    smallDiv[p].getCoef(m, k)->setCoef(divisor, j);

                    if(cabs(divisor) < threshold) //if small divisor, we set the rvf
                    {
                        bux = eta[p].getCoef(m, k)->ofs_getCoef(j);
                        fh[p].getCoef(m, k)->setCoef(0.0*I-bux, j);

                        if(true)
                        {
                            cout <<  setprecision(1);
                            cout << "Small divisor: p = " << p << ", kv = ("<< kv[0] << ", " << kv[1] << ", " << kv[2]  << ", " << kv[3]  << "), j = " << j << "." << endl;
                            cout <<  setprecision(4);
                            cout << "Divisor = " << creal(divisor) << "  " << cimag(divisor)  << ", -eta    = " << creal(0.0*I-bux) << "  " << cimag(0.0*I-bux) << endl;
                            cout << "-------------" << endl;
                        }

                    }
                    else //normal form
                    {
                        bux = 0.0*I-eta[p].getCoef(m, k)->ofs_getCoef(j)/divisor;
                        xi[p].getCoef(m, k)->setCoef(bux, j);
                    }
                }
                //Update kv
                if(k < FTDA::nmon(REDUCED_NV, m)-1)  FTDA::prxkt(kv, REDUCED_NV);
            }
        }

        //------------------------------------------
        // Wh = H2*xi @ order m
        // Recall that H2 = H in matrix format
        //------------------------------------------
        smvprod_u(H2, xi, Wh, m);

        //------------------------------------------
        // Update, W, Wt and W2 @ order m
        //------------------------------------------
        applyCOC(PC,  V,  Wh, W,  m, 1);     //order zero is added
        applyCOC(PC,  V,  Wh, Wt, m, 0);     //order zero is NOT added


        //------------------------------------------
        // Update the differential to take into account the new terms @ order m
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, m);


        //------------------------------------------------------------------------------------
        //
        // Note: after this point, only necessary for testing, can be discarded for speed
        //
        //------------------------------------------------------------------------------------
        initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, SEML, m);

        //------------------------------------------
        // update the VF
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, m);

        //------------------------------------------
        // dot(W)
        //------------------------------------------
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, m);

        //------------------------------------------
        // Differential of Wh
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, m);

        //------------------------------------------
        // DWh times f @order m
        //------------------------------------------
        smvprod_t(C_DW,  fh, C_DWf,  m);

        cout << "pm. end of update testing." << " in " << toc() << " s. "<< endl;
        cout << "  pm. End of order " << m << " in " << toc() << " s. " << endl;

        //If txt files have to be updated at each order
        if(OutputEachOrder)
        {
            tic();
            //print W & FW
            vector_fprinf(W,     F_NF+"W/W");
            vector_fprinf(smallDiv, F_NF+"W/smallDiv");
            vector_fprinf(Wt,    F_NF+"W/Wt");
            vector_fprinf(Wh,    F_NF+"W/Wh");
            vector_fprinf(FW,    F_NF+"FW/FW");
            vector_fprinf(C_FW,  F_NF+"FW/C_FW");
            vector_fprinf(FWh,   F_NF+"FW/FWh");
            //fprintf Wdot
            vector_fprinf(C_Wdot,  F_NF+"Wdot/C_Wdot");
            //Reduced vector field
            vector_fprinf(fh, F_NF+"rvf/fh");
            //Potential
            vector_fprinf(Unx, F_NF+"Potential/Unx");
            vector_fprinf(Uny, F_NF+"Potential/Uny");
            vector_fprinf(Unz, F_NF+"Potential/Unz");
            vector_fprinf(En, F_NF+"Potential/En");
            vector_fprinf(Mn, F_NF+"Potential/Mn");
            vector_fprinf(Sn, F_NF+"Potential/Sn");
            vector_fprinf(Un, F_NF+"Potential/Un");
            //DW*vf
            vector_fprinf(C_DWf,  F_NF+"DWf/C_DWf");
            cout << "  pm. End of printing of order " << m  <<  " in " << toc() << " s. " << endl;
        }
    }


    //------------------------------------------
    // Printing
    //------------------------------------------
    //If txt files have to be updated at the very last order
    if(Output && !OutputEachOrder)
    {
        tic();
        //print W & FW
        vector_fprinf(W,     F_NF+"W/W");
        vector_fprinf(smallDiv, F_NF+"W/smallDiv");
        vector_fprinf(Wt,    F_NF+"W/Wt");
        vector_fprinf(Wh,    F_NF+"W/Wh");
        vector_fprinf(FW,    F_NF+"FW/FW");
        vector_fprinf(C_FW,  F_NF+"FW/C_FW");
        vector_fprinf(FWh,   F_NF+"FW/FWh");
        //fprintf Wdot
        vector_fprinf(C_Wdot,  F_NF+"Wdot/C_Wdot");
        //Reduced vector field
        vector_fprinf(fh, F_NF+"rvf/fh");
        //Potential
        vector_fprinf(Unx, F_NF+"Potential/Unx");
        vector_fprinf(Uny, F_NF+"Potential/Uny");
        vector_fprinf(Unz, F_NF+"Potential/Unz");
        vector_fprinf(En,  F_NF+"Potential/En");
        vector_fprinf(Mn, F_NF+"Potential/Mn");
        vector_fprinf(Sn, F_NF+"Potential/Sn");
        vector_fprinf(Un, F_NF+"Potential/Un");
        //DW*vf
        vector_fprinf(C_DWf,  F_NF+"DWf/C_DWf");
        cout << "  pm. End of printing" <<  " in " << toc() << " s. " << endl;
    }

}

/**
 *  \brief Compute the 6D parameterization of the transit/non-transit trajectories about the Earth-Moon libration points L1,2, up to a given order (mixed style).
 *  \param OutputEachOrder: if true, the results are printed in txt files after each computed order.
 *  \param Output: if true and OutputEachOrder = false, the results are printed at the very end of the computation.
 *
 *   Inputs:
 *      - The QBCP_L structure SEML, which contains the addresses of the data folders (F_COC, F_GS...),
 *        the order of the expansions (both Fourier and Fourier-Taylor),
 *        and some useful parameters (gamma, c1)
 *
 *      - In the case of the QBCP and BCP, the COC in data txt files, computed with the routines contained in nfo2.cpp.
 *        In the case of autonomous vector fields, such as the CRTBP, the COC is computed within the routine pm.
 *
 *  Outputs (stored in the folder SEML.cs.F_GS):
 *      - W(s,t),  the parameterization of the central manifold in Normalized-Centered coordinates.
 *      - Wt(s,t), the parameterization of the central manifold in Normalized-Centered coordinates, translated to the periodic orbit of L1,2.
 *      - Wh(s,t), the parameterization of the central manifold in TFC coordinates.
 *
 *      - FW(s,t),   partial the vector field dot(W), without the term of current order n in W(s,t).
 *      - FWh(s,t),   partial the vector field dot(Wh), without the term of current order n in Wh(s,t).
 *      - C_FW(s,t), the complete vector field dot(W).
 *      - fh(s,t), the reduced vector field: dot(s) = fh(s,t).
 *
 *      - En(s,t), the Earth potential up to order the order of the parameterization.
 *      - Mn(s,t), the Moon  potential up to order the order of the parameterization.
 *      - Sn(s,t), the Sun   potential up to order the order of the parameterization.
 *      - Un(s,t), the Sun+Earth+Moon potential up to order the order of the parameterization.
 *      - Unx(s,t) = dr(Un)/drx, with dr = partial derivative
 *      - Uny(s,t) = dr(Un)/dry, with dr = partial derivative
 *      - Unz(s,t) = dr(Un)/drz, with dr = partial derivative
 *
 *      - Wdot(s,t) : the partial derivative of W(s,t) with respect to time (not equal to dot(W) = dW/dt)
 *
 *  Notes:
 *     The routine apply a mixed style parameterization in order to separate as many submanifolds as possible. Here is a list, with the corresponding sets of indices
 *     used for the mixed style (see Haro 2014, section 2.2.3).
 *
 *            set I          |  submanifold described by si = 0, for all i in I
 *          ----------------------------------------------------------------
 *           {1,2,3,4,6}     |  unstable manifold
 *           {1,2,3,4,5}     |  stable manifold
 *           {1,2,3,4}       |  hyperbolic normal part
 *           {2,4,5,6}       |  planar family
 *           {2,4,6}         |  unstable manifold of the planar family
 *           {2,4,5}         |  stable manifold of the planar family
 *           {2,4}           |  hyperbolic normal part of the planar family
 *           {1,3,5,6}       |  vertical family
 *           {1,3,6}         |  unstable manifold of the vertical family
 *           {1,3,5}         |  stable manifold of the vertical family
 *           {1,3}           |  hyperbolic normal part of the vertical family
 *           {5,6}           |  center manifold
 *           {6}             |  center-unstable manifold
 *           {5}             |  center-stable manifold
 *
 **/
void pm_mixedstyle(int OutputEachOrder, int Output)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "    Application of the parameterization method     " << endl;
    cout << "  on the non-autonomous vector field of the QBCP  " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout <<  setw(5) << setprecision(4);
    //Misc complex numbers
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_MS;
    string F_COC   = SEML.cs.F_COC;
    string F_COEF  = SEML.cs.F_COEF;

    cout << "F_COEF = " << F_COEF << endl;
    cout << "F_COC = "  << F_COC << endl;
    cout << "F_GS = "   << F_GS << endl;


    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    matrix<Ofsc> P(6,6);      //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> Q(6,6);      //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    vector<Ofsc> V(6);        //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> Xe(3);       //The vectors Xe, Xm, Xs, modified position
    vector<Ofsc> Xm(3);       //of the primaries used to compute the potential:
    vector<Ofsc> Xs(3);       //Xf = True Xf - V, for f = e, m, s
    Ofsc ILe(OFS_ORDER);      //The Ofs ILe, ILm and ILs, modified inverse orbit radii
    Ofsc ILm(OFS_ORDER);      //of the primaries used to compute the potential
    Ofsc ILs(OFS_ORDER);      //ILf = 1.0/||Xf||, for f = e, m, s
    matrix<Ofsc> PC(6,6);     //PC=P(theta)*C
    matrix<Ofsc> PCdot(6,6);  //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> CQ(6,6);     //CQ=Cinv*Pinv=inv(PC)
    vector<Ofsc> Vdot(6);     //Vdot=dot(V)

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Initialization of the alphas
    //------------------------------------------
    int noc = SEML.numberOfCoefs;
    vector<Ofsc> alpha(noc);
    //Switch on the model to initialize the alphas
    switch(SEML.model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    case Csts::ERTBP:
    {
        //Read from txt files
        ifstream readStream;
        string ss1;
        for(int i = 0; i < noc; i++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
            readOFS_txt(alpha[i], F_COEF+"alpha"+ss1+"_fft");
        }

        break;
    }

    case Csts::CRTBP:
    {
        cout << "pm. The use of the RTBP has been detected when initializing the alphas." << endl;
        //All coeffs to zero
        for(int j = 0; j < noc;  j++) alpha[j].setCoef(0.0, 0);
        //Non-null coeffs
        alpha[0].setCoef(1.0, 0);               //alpha1  = 1.0
        alpha[2].setCoef(1.0, 0);               //alpha3  = 1.0
        alpha[5].setCoef(1.0, 0);               //alpha6  = 1.0
        alpha[8].setCoef(SEML.us_em.mu_EM, 0);  //alpha9  = Xe
        alpha[10].setCoef(SEML.us_em.mu_EM, 0); //alpha11 = Xm
        alpha[12].setCoef(-SEML.cs_em.c1, 0);   //alpha13 = -c1
        break;
    }

    default:
    {
        cout << "error in pm. Unknown model." << endl;
        return;
    }

    }



    //------------------------------------------
    // Getting DB from files
    //------------------------------------------
    //     | iw1 0    0    0    0   0   |
    //     | 0   w2   0    0    0   0   |
    //     | 0   0    iw3  0    0   0   |
    //DB = | 0   0    0   -iw1  0   0   |
    //     | 0   0    0    0   -w2  0   |
    //     | 0   0    0    0    0  -iw3 |
    gsl_matrix_complex *DB = gsl_matrix_complex_calloc(6, 6);
    if(SEML.model == Csts::QBCP || SEML.model == Csts::BCP)
    {
        //Reading an OFS from a text file
        string filename = F_COC+"DB";
        glsc_matrix_complex_read(DB, filename);
    }
    else //RTBP
    {
        cout << "pm. The use of the RTBP has been detected when initializing DB." << endl;
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, c2;
        c2 = SEML.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);

        cout << "omega1 = " << om1 << endl;
        cout << "la1    = " << la1 << endl;
        cout << "omega2 = " << om2 << endl;

        //--------------------
        //Init DB
        //--------------------
        gsl_matrix_complex_set(DB, 0, 0, gslc_complex(+om1*I));
        gsl_matrix_complex_set(DB, 1, 1, gslc_complex(+la1+0.0*I));
        gsl_matrix_complex_set(DB, 2, 2, gslc_complex(+om2*I));
        gsl_matrix_complex_set(DB, 3, 3, gslc_complex(-om1*I));
        gsl_matrix_complex_set(DB, 4, 4, gslc_complex(-la1+0.0*I));
        gsl_matrix_complex_set(DB, 5, 5, gslc_complex(-om2*I));
    }
    cout << setprecision(15) << endl;


    //------------------------------------------
    //DF0 = DB in matrix<Ofsc> format
    //------------------------------------------
    Ofsc BUX(OFS_ORDER);
    matrix<Ofsc> DF0(6, 6);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    DF0.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    DF0.setCoef(BUX, 1, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    DF0.setCoef(BUX, 2, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    DF0.setCoef(BUX, 3, 3);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    DF0.setCoef(BUX, 4, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    DF0.setCoef(BUX, 5, 5);


    //------------------------------------------
    // Setting H, Hinv, Lambda, LamdaL, LambdaN from DB
    // A this step, the central manifold is selected.
    //
    // TO-DO: something more generic to be able
    // to compute the central-stable and central-unstable manifold
    //
    //------------------------------------------
    gsl_matrix_complex *H     = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *Hinv  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *La    = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex *LaL   = gsl_matrix_complex_calloc(6, 6);
    //H = (L N), with N null, H = L
    gsl_matrix_complex *L    = gsl_matrix_complex_calloc(6, 6);

    //------------------------------------------
    //    |  iw1 0    0    0    0   0  |
    //    |  0   0    0    0    w2  0  |
    //    |  0   iw3  0    0    0   0  |
    //H = |  0   0   -iw1  0    0   0  | = (L)
    //    |  0   0    0    0    0  -w2 |
    //    |  0   0    0   -iw3  0   0  |
    //------------------------------------------
    gsl_matrix_complex_set(H, 0, 0, gsl_matrix_complex_get(DB, 0, 0));
    gsl_matrix_complex_set(H, 1, 4, gsl_matrix_complex_get(DB, 1, 1));
    gsl_matrix_complex_set(H, 2, 1, gsl_matrix_complex_get(DB, 2, 2));
    gsl_matrix_complex_set(H, 3, 2, gsl_matrix_complex_get(DB, 3, 3));
    gsl_matrix_complex_set(H, 4, 5, gsl_matrix_complex_get(DB, 4, 4));
    gsl_matrix_complex_set(H, 5, 3, gsl_matrix_complex_get(DB, 5, 5));

    //------------------------------------------
    //H2 = H in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> H2(6, 6);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    H2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    H2.setCoef(BUX, 1, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    H2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    H2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    H2.setCoef(BUX, 4, 5);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    H2.setCoef(BUX, 5, 3);

    //------------------------------------------
    //L = H, here
    //------------------------------------------
    gsl_matrix_complex_view H11 = gsl_matrix_complex_submatrix (H , 0 , 0 , 6 , 6 );
    gsl_matrix_complex_memcpy(L, &H11.matrix);

    //matrix version of L
    matrix<Ofsc> L2(6, 6);
    BUX.zero();
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    L2.setCoef(BUX, 0, 0);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    L2.setCoef(BUX, 1, 4);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    L2.setCoef(BUX, 2, 1);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    L2.setCoef(BUX, 3, 2);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    L2.setCoef(BUX, 4, 5);
    BUX.setCoef(gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    L2.setCoef(BUX, 5, 3);


    //------------------------------------------
    //       |1/iw1  0    0      0     0     0    |
    //       | 0     0   1/iw3   0     0     0    |
    //       | 0     0    0    -1/iw1  0     0    |
    //Hinv = | 0     0    0      0     0  -1/iw3  |
    //       | 0    1/w2  0      0     0     0    |
    //       | 0     0    0      0   -1/w2   0    |
    //------------------------------------------
    gsl_matrix_complex_set(Hinv, 0, 0, gsl_complex_inverse(gsl_matrix_complex_get(DB, 0, 0)));
    gsl_matrix_complex_set(Hinv, 4, 1, gsl_complex_inverse(gsl_matrix_complex_get(DB, 1, 1)));
    gsl_matrix_complex_set(Hinv, 1, 2, gsl_complex_inverse(gsl_matrix_complex_get(DB, 2, 2)));
    gsl_matrix_complex_set(Hinv, 2, 3, gsl_complex_inverse(gsl_matrix_complex_get(DB, 3, 3)));
    gsl_matrix_complex_set(Hinv, 5, 4, gsl_complex_inverse(gsl_matrix_complex_get(DB, 4, 4)));
    gsl_matrix_complex_set(Hinv, 3, 5, gsl_complex_inverse(gsl_matrix_complex_get(DB, 5, 5)));

    //------------------------------------------
    //Hinv2 = Hinv in matrix<Ofsc> format
    //------------------------------------------
    matrix<Ofsc> Hinv2(6, 6);
    BUX.zero();
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 0, 0)), 0);
    Hinv2.setCoef(BUX, 0, 0);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 1, 1)), 0);
    Hinv2.setCoef(BUX, 4, 1);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 2, 2)), 0);
    Hinv2.setCoef(BUX, 1, 2);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 3, 3)), 0);
    Hinv2.setCoef(BUX, 2, 3);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 4, 4)), 0);
    Hinv2.setCoef(BUX, 5, 4);
    BUX.setCoef(1.0/gslc_complex(gsl_matrix_complex_get(DB, 5, 5)), 0);
    Hinv2.setCoef(BUX, 3, 5);

    //------------------------------------------
    //                     | LaL  T  |
    //Lambda = Hinv*DB*H = | 0  LaN  | with T = 0 in this context
    //
    //------------------------------------------
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(6, 6);
    //AUX = DB*H
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , DB , H , zero_c , AUX );
    //B = Hinv*AUX
    gsl_blas_zgemm (CblasNoTrans , CblasNoTrans , one_c , Hinv , AUX , zero_c , La);
    //------------------------------------------
    gsl_matrix_complex_view La11 = gsl_matrix_complex_submatrix (La , 0 , 0 , 6 , 6 );   //S11
    gsl_matrix_complex_memcpy(LaL, &La11.matrix);

    cout << "real(Lambda) = " << endl;
    gslc_matrix_complex_printf_real(La);
    cout << "imag(Lambda) = " << endl;
    gslc_matrix_complex_printf_imag(La);


    cout << "  pm. Initialisation of the TFC manifold.  " << endl;
    //------------------------------------------
    // Initialisation of the TFC manifold ("W hat")
    //------------------------------------------
    cdouble db;
    gsl_complex dbc;
    //vector of size 6
    vector<Oftsc> Wh(6);
    //The order 0 is equal to zero, the order 1 is equal to H, which is diagonal
    //----------------------
    //     |   iw1*s1 |
    //     |    w2*s5 |
    //     |   iw3*s2 |
    //Wh = |  -iw1*s3 |
    //     |  - w2*s6 |
    //     |  -iw3*s4 |
    //----------------------
    // Wh[0] = L11*s1 = L[0][0]*s[0]
    dbc = gsl_matrix_complex_get(L, 0, 0);
    db  = gslc_complex(dbc);
    Wh[0].setCoef0(db, 1, 0);
    // Wh[1] = L25*s5 = L[1][4]*s[4]
    dbc = gsl_matrix_complex_get(L, 1, 4);
    db  = gslc_complex(dbc);
    Wh[1].setCoef0(db, 1, 4);
    // Wh[2] = L32*s2 = L[2][1]*s[1]
    dbc = gsl_matrix_complex_get(L, 2, 1);
    db  = gslc_complex(dbc);
    Wh[2].setCoef0(db, 1, 1);
    // Wh[3] = L43*s3 = L[3][2]*s[2]
    dbc = gsl_matrix_complex_get(L, 3, 2);
    db  = gslc_complex(dbc);
    Wh[3].setCoef0(db, 1, 2);
    // Wh[4] = L56*s6 = L[4][5]*s[5]
    dbc = gsl_matrix_complex_get(L, 4, 5);
    db  = gslc_complex(dbc);
    Wh[4].setCoef0(db, 1, 5);
    // Wh[5] = L64*s4 = L[5][3]*s[3]
    dbc = gsl_matrix_complex_get(L, 5, 3);
    db  = gslc_complex(dbc);
    Wh[5].setCoef0(db, 1, 3);

    cout << "  pm. Initialisation of the NC manifolds.  " << endl;
    //------------------------------------------
    // Initialisation of the NC manifolds
    //------------------------------------------
    vector<Oftsc> W(6);   //W = TFC(Wh) (in NC coordinates)
    vector<Oftsc> Wt(6);  //Wt = W - V  (in NC coordinates)
    //Applying the COC (with and without zero order) @order 0 and 1
    for(int i = 0; i < 2; i++)
    {
        applyCOC(PC,  V,  Wh, W,  i, 1);   //order zero is added
        applyCOC(PC,  V,  Wh, Wt, i, 0);   //order zero is NOT added
    }

    //------------------------------------------
    // Differential of Wh (6*6)
    //------------------------------------------
    //Contribution of order < k to k
    matrix<Oftsc> DWh(6,6);  //TFC
    //Contribution of order <=k to k
    matrix<Oftsc> C_DW(6,6);   //NC

    cout << "  pm. Initialisation of the complete VF.  " << endl;
    //------------------------------------------
    // Initialisation of the complete VF
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> FWh(6);  //TFC
    vector<Oftsc> FW(6);   //NC
    vector<Oftsc> FW2(6);  //Temp
    vector<Oftsc> FW3(6);  //Temp
    //Contribution of order <=k to k
    vector<Oftsc> C_FW(6);   //NC

    cout << "  pm. Initialisation of the reduced VF.  " << endl;
    //------------------------------------------
    // Initialisation of the reduced VF
    //------------------------------------------
    vector<Oftsc> fh(6);     //TFC
    for(int i = 0; i < 6; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        fh[i].setCoef0(db, 1, i);
    }

    cout << "  pm. Initialisation of LaLs.  " << endl;
    //LaLs = fh @ order 0 and 1
    vector<Oftsc> LaLs(6);
    for(int i = 0; i < 6; i++)
    {
        // fh[i] =  LaL[i][i]*s[i]
        dbc = gsl_matrix_complex_get(LaL, i, i);
        db = gslc_complex(dbc);
        LaLs[i].setCoef0(db, 1, i);
    }

    //------------------------------------------
    // Initialisation of DW times fh
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> DWfh(6);  //TFC
    //Contribution of order <=k to k
    vector<Oftsc> C_DWf(6);   //NC

    //------------------------------------------
    // Initialisation of dot(W)
    //------------------------------------------
    //Contribution of order < k to k
    vector<Oftsc> Whdot(6);      //TFC
    vector<Oftsc> Wdot(6);       //NC
    //Contribution of order <=k to k
    vector<Oftsc> C_Wdot(6);     //NC


    //------------------------------------------
    // Initialisation of the Legendre objects
    //------------------------------------------
    //Potential
    vector<Oftsc> En(OFTS_ORDER+2);     //Earth
    vector<Oftsc> Sn(OFTS_ORDER+2);     //Moon
    vector<Oftsc> Mn(OFTS_ORDER+2);     //Sun
    vector<Oftsc> Un(OFTS_ORDER+2);     //Earth+Moon+Sun
    //Potential derivatives
    vector<Oftsc> Enx(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Eny(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Enz(OFTS_ORDER+2);    //Earth
    vector<Oftsc> Mnx(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mny(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Mnz(OFTS_ORDER+2);    //Moon
    vector<Oftsc> Snx(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Sny(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Snz(OFTS_ORDER+2);    //Sun
    vector<Oftsc> Unx(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Uny(OFTS_ORDER+2);    //Earth+Moon+Sun
    vector<Oftsc> Unz(OFTS_ORDER+2);    //Earth+Moon+Sun

    //------------------------------------------
    // Initialisation of additionnal OFTS objects (see comments)
    // Function of Wh ("W hat")
    // Updated in updatePotentialExpansion routine
    //------------------------------------------
    Oftsc Rho2(REDUCED_NV, OFTS_ORDER, 1, OFS_ORDER);   //=xt^2+yt^2+zt^2
    Oftsc xse(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xe*xt+ye*yt
    Oftsc xsm(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xm*xt+ym*yt
    Oftsc xss(REDUCED_NV,  OFTS_ORDER, 1, OFS_ORDER);   //=xs*xt+ys*yt

    //------------------------------------------
    // Initialisation of various objects
    // Used to solve the pm equations
    //------------------------------------------
    vector<Oftsc> E(6);    //[DWh<m fh<m]m-[Fh(Wh<m)]m
    vector<Oftsc> eta(6);  //eta = Hinv*E
    vector<Oftsc> xi(6);   //xi = Hinv*Wh at order >=2
    cdouble lap;       //contains LaN
    cdouble lar;       //contains LaL*r for each monomial
    cdouble aux;       //temp obj
    cdouble bux;       //temp obj
    int kv[REDUCED_NV];       //contains the indices for each monomial
    cdouble divisor;          //divisor in the cohomological equations
    double threshold = 1e-2;  //threshold for small divisors = 1e-2

    //------------------------------------------
    // Indices for mixed style
    //------------------------------------------
    int **VIs;
    int  *VIn;
    VIs = (int **) calloc(14, sizeof(int*));
    VIn = (int *)  calloc(14, sizeof(int));

    //--------------------------
    // Unstable manifold
    //--------------------------
    VIn[0] = 5;
    VIs[0] = (int *) calloc(5, sizeof(int));
    VIs[0][0] = 1;
    VIs[0][1] = 2;
    VIs[0][2] = 3;
    VIs[0][3] = 4;
    VIs[0][4] = 6;

    //--------------------------
    // Stable manifold
    //--------------------------
    VIn[1] = 5;
    VIs[1] = (int *) calloc(5, sizeof(int));
    VIs[1][0] = 1;
    VIs[1][1] = 2;
    VIs[1][2] = 3;
    VIs[1][3] = 4;
    VIs[1][4] = 5;

    //--------------------------
    // Hyperbolic normal part (stable + unstable)
    //--------------------------
    VIn[2] = 4;
    VIs[2] = (int *) calloc(4, sizeof(int));
    VIs[2][0] = 1;
    VIs[2][1] = 2;
    VIs[2][2] = 3;
    VIs[2][3] = 4;

    //--------------------------
    // Planar family
    //--------------------------
    VIn[3] = 4;
    VIs[3] = (int *) calloc(4, sizeof(int));
    VIs[3][0] = 2;
    VIs[3][1] = 4;
    VIs[3][2] = 5;
    VIs[3][3] = 6;

    //--------------------------
    // Unstable manifold of the planar family
    //--------------------------
    VIn[4] = 3;
    VIs[4] = (int *) calloc(3, sizeof(int));
    VIs[4][0] = 2;
    VIs[4][1] = 4;
    VIs[4][2] = 6;

    //--------------------------
    // Stable manifold of the planar family
    //--------------------------
    VIn[5] = 3;
    VIs[5] = (int *) calloc(3, sizeof(int));
    VIs[5][0] = 2;
    VIs[5][1] = 4;
    VIs[5][2] = 5;

    //--------------------------
    // Hyperbolic normal part of the planar family
    //--------------------------
    VIn[6] = 2;
    VIs[6] = (int *) calloc(2, sizeof(int));
    VIs[6][0] = 2;
    VIs[6][1] = 4;

    //--------------------------
    // Vertical family
    //--------------------------
    VIn[7] = 4;
    VIs[7] = (int *) calloc(4, sizeof(int));
    VIs[7][0] = 1;
    VIs[7][1] = 3;
    VIs[7][2] = 5;
    VIs[7][3] = 6;

    //--------------------------
    // Unstable manifold of the vertical family
    //--------------------------
    VIn[8] = 3;
    VIs[8] = (int *) calloc(3, sizeof(int));
    VIs[8][0] = 1;
    VIs[8][1] = 3;
    VIs[8][2] = 6;

    //--------------------------
    // Stable manifold of the vertical family
    //--------------------------
    VIn[9] = 3;
    VIs[9] = (int *) calloc(3, sizeof(int));
    VIs[9][0] = 1;
    VIs[9][1] = 3;
    VIs[9][2] = 5;

    //--------------------------
    // Hyperbolic normal part of the vertical family
    //--------------------------
    VIn[10] = 2;
    VIs[10] = (int *) calloc(2, sizeof(int));
    VIs[10][0] = 1;
    VIs[10][1] = 3;

    //--------------------------
    // Center manifold
    //--------------------------
    VIn[11] = 2;
    VIs[11] = (int *) calloc(2, sizeof(int));
    VIs[11][0] = 5;
    VIs[11][1] = 6;

    //--------------------------
    // Center-unstable manifold
    //--------------------------
    VIn[12] = 1;
    VIs[12] = (int *) calloc(1, sizeof(int));
    VIs[12][0] = 6;

    //--------------------------
    // Center-stable manifold
    //--------------------------
    VIn[13] = 1;
    VIs[13] = (int *) calloc(1, sizeof(int));
    VIs[13][0] = 5;


    cout << "  pm. Initialization of the recurrence.  " << endl;

    //-------------------------------------------------------------
    // Initialization of the recurrence
    //-------------------------------------------------------------
    //The potential is updated: Un, Unx, Uny, Unz
    // @order 0, 1 and 2
    //Wt is used to update the Potential
    updatePotentialExpansion(Wt, Rho2, En,
                             Mn, Sn, Un,
                             Enx, Eny, Enz,
                             Mnx, Mny, Mnz,
                             Snx, Sny, Snz,
                             Unx, Uny, Unz,
                             Xe, Xm, Xs,
                             ILe, ILm, ILs,
                             xse, xsm, xss,
                             SEML, 2);

    cout << "  pm. End of initialisation.  " << endl;

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ order 0 and 1
    //------------------------------------------
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    for(int k = 0; k < 2 ; k++)
    {
        //vector field
        //-----------------
        //applyVF(alpha, W, FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);
        //Copy for consistency testing
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, k);
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, k);

        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, k);

        //dot(W)
        //-----------------
        for(int i = 0; i < 6; i++) Wdot[i].dot(W[i],   SEML.us.n, k);     //in NC
        for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], SEML.us.n, k);     //in TFC
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, k);     //in NC, complete contribution

        //----------------------
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, k);   //in NC
        if(k > 0) for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, k);   //in NC, complete contribution

        // DWh times f @order k
        //----------------------
        //smvprod_t(DW,  fh, DWf,  k);     //in NC
        smvprod_t(DWh, fh, DWfh, k);       //in TFC
        smvprod_t(C_DW,  fh, C_DWf,  k);   //in NC, complete contribution


    }

    cout << "  pm. End of order 0/1.  " << endl;

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Computing @ orders >= 2
    //------------------------------------------
    for(int m = 2; m <= OFTS_ORDER; m++)
    {

        cout << "pm. start of order " << m << endl;
        tic();
        //------------------------------------------
        //The potential is updated: Un, Unx, Uny, Unz
        //Wt is used to update the Potential
        //What is updated: contributions of order < m to order m
        //------------------------------------------
        if(m > 2) updatePotentialExpansion(Wt, Rho2, En,
                                               Mn, Sn, Un,
                                               Enx, Eny, Enz,
                                               Mnx, Mny, Mnz,
                                               Snx, Sny, Snz,
                                               Unx, Uny, Unz,
                                               Xe, Xm, Xs,
                                               ILe, ILm, ILs,
                                               xse, xsm, xss,
                                               SEML, m);

        //------------------------------------------
        // Applying the VF: [F(W<m)]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);

        //------------------------------------------
        //FWh = COCinv(FW)
        //------------------------------------------
        //FWh = COCinv(FW)
        //-----------------
        applyInvDotCOC(CQ,PCdot,Vdot, Wh, FW, FWh, FW2, FW3, m);

        //cout << "pm. end of applyInvDotCOC." << endl;

        //------------------------------------------
        // dot(W<m)
        // It is NOT updated here: the terms of order < m
        // dot not influence dot(W) at order m
        //------------------------------------------
        //for(int i = 0; i < 6; i++) Wdot[i].dot(W[i], qbcp.n, m);
        //for(int i = 0; i < 6; i++) Whdot[i].dot(Wh[i], qbcp.n, m);
        //cout << "pm. end of dot." << " in " << toc() << " s. " << endl;

        //------------------------------------------
        // Differential of Wh
        // It is NOT updated here, but at the end of the loop.
        // Indeed, Wh & W contain no order m at this point!
        //------------------------------------------

        //------------------------------------------
        // [DW<m times f<m]m
        // What is updated: contributions of order < m to order m
        //------------------------------------------
        smvprod_t(DWh, fh, DWfh, m);
        //smvprod_t(DW, fh, DWf, m);


        //------------------------------------------
        // We have [Fh(Wh<m)]m and [DWh<m fh<m]m
        // Update of eta at order m
        //------------------------------------------
        // [E]m = [DWh<m fh<m]m - [Fh(Wh<m)]m
        for(int i = 0 ; i< 6; i++) E[i].ofts_fsum_u(FWh[i], -1.0+0.0*I, DWfh[i], 1.0+0.0*I, m);
        //[eta]m = [Hinv*E]m;
        smvprod_u(Hinv2, E, eta, m);


        //------------------------------------------
        // Equation in the tangent space:
        // First four components of eta
        //------------------------------------------
        //Sum on the tangential components of eta (first d components)
        for(int p = 0; p< Csts::NV; p++)
        {
            //Update lap
            lap = gslc_complex(gsl_matrix_complex_get(La, p, p));
            //Reset kv
            kv[0] = m;
            for(int i = 1; i< REDUCED_NV; i++) kv[i] = 0;
            //Sum on all the coefficients of eta^p_m
            for(int k = 0; k< FTDA::nmon(REDUCED_NV, m); k++)
            {
                //Reset aux
                aux = 0.0+0.0*I;
                //Update aux = laL.kv
                for(int ii = 0; ii < REDUCED_NV; ii++)
                {
                    lar = gslc_complex(gsl_matrix_complex_get(LaL, ii, ii));
                    aux += lar*(kv[ii]+0.0*I);
                }

                //Update xi
                for(int j = -OFS_ORDER; j<= OFS_ORDER; j++)
                {
                    //divisor = jnI - (lap - laL.kv)
                    divisor = SEML.us.n*I*j - (lap - aux);

                    //if i is in VIs[vi] and kv[j] == 0 for all j in VIs[vi]
                    if(isNF(p, kv, VIs, VIn, 14))
                    {
                        //Normal form
                        bux = 0.0*I-eta[p].getCoef(m, k)->ofs_getCoef(j)/divisor;
                        xi[p].getCoef(m, k)->setCoef(bux, j);

                        //if small divisor, send a warning. Should not happend if the set of indices is good
                        if(cabs(divisor) < threshold)
                        {

                            cout <<  setprecision(1);
                            cout << "Small divisor: p = " << p << ", kv = ("<< kv[0] << ", " << kv[1] << ", " << kv[2]
                                                                            << ", " << kv[3] << ", " << kv[4] << ", " << kv[5] << "), j = "
                                                                            << j << "." << endl;
                            cout <<  setprecision(4);
                            cout << "Divisor = " << creal(divisor) << "  " << cimag(divisor)  << ", -eta    = " << creal(0.0*I-bux) << "  " << cimag(0.0*I-bux) << endl;
                            cout << "-------------" << endl;
                        }
                    }
                    else
                    {
                        //Graph style
                        bux = eta[p].getCoef(m, k)->ofs_getCoef(j);
                        fh[p].getCoef(m, k)->setCoef(0.0*I-bux, j);
                    }
                }
                //Update kv
                if(k < FTDA::nmon(REDUCED_NV, m)-1)  FTDA::prxkt(kv, REDUCED_NV);
            }
        }

        //------------------------------------------
        // Wh = H2*xi @ order m
        // Recall that H2 = H in matrix format
        //------------------------------------------
        smvprod_u(H2, xi, Wh, m);

        //------------------------------------------
        // Update, W, Wt and W2 @ order m
        //------------------------------------------
        applyCOC(PC,  V,  Wh, W,  m, 1);     //order zero is added
        applyCOC(PC,  V,  Wh, Wt, m, 0);     //order zero is NOT added


        //------------------------------------------
        // Update the differential to take into account the new terms @ order m
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) DWh.der(Wh[i], j+1, i, j, m);


        //------------------------------------------------------------------------------------
        //
        // Note: after this point, only necessary for testing, can be discarded for speed
        //
        //------------------------------------------------------------------------------------
        initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, SEML, m);

        //------------------------------------------
        // update the VF
        //------------------------------------------
        //W is used in the vector field
        applyVF(alpha, W, C_FW, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, m);
        //applyVF(alpha, W, C_FW, Unx, Uny, Unz, m);

        //------------------------------------------
        // dot(W)
        //------------------------------------------
        for(int i = 0; i < 6; i++) C_Wdot[i].dot(W[i], SEML.us.n, m);

        //------------------------------------------
        // Differential of Wh
        //------------------------------------------
        for(int i = 0 ; i < Csts::NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) C_DW.der(W[i], j+1, i, j, m);

        //------------------------------------------
        // DWh times f @order m
        //------------------------------------------
        smvprod_t(C_DW,  fh, C_DWf,  m);

        cout << "pm. end of update testing." << " in " << toc() << " s. "<< endl;
        cout << "  pm. End of order " << m << " in " << toc() << " s. " << endl;

        //If txt files have to be updated at each order
        if(OutputEachOrder)
        {
            tic();
            //print W & FW
            vector_fprinf(W,     F_GS+"W/W");
            vector_fprinf(Wt,    F_GS+"W/Wt");
            vector_fprinf(Wh,    F_GS+"W/Wh");
            vector_fprinf(FW,    F_GS+"FW/FW");
            vector_fprinf(C_FW,  F_GS+"FW/C_FW");
            vector_fprinf(FWh,   F_GS+"FW/FWh");
            //fprintf Wdot
            vector_fprinf(C_Wdot,  F_GS+"Wdot/C_Wdot");
            //Reduced vector field
            vector_fprinf(fh, F_GS+"rvf/fh");
            //Potential
            vector_fprinf(Unx, F_GS+"Potential/Unx");
            vector_fprinf(Uny, F_GS+"Potential/Uny");
            vector_fprinf(Unz, F_GS+"Potential/Unz");
            vector_fprinf(En, F_GS+"Potential/En");
            vector_fprinf(Mn, F_GS+"Potential/Mn");
            vector_fprinf(Sn, F_GS+"Potential/Sn");
            vector_fprinf(Un, F_GS+"Potential/Un");
            //DW*vf
            vector_fprinf(C_DWf,  F_GS+"DWf/C_DWf");
            cout << "  pm. End of printing of order " << m  <<  " in " << toc() << " s. " << endl;
        }
    }


    //------------------------------------------
    // Printing
    //------------------------------------------
    //If txt files have to be updated at the very last order
    if(Output && !OutputEachOrder)
    {
        tic();

        //------------------------------------------
        // Txt
        //------------------------------------------
        //print W & FW
        vector_fprinf(W,     F_GS+"W/W");
        vector_fprinf(Wt,    F_GS+"W/Wt");
        vector_fprinf(Wh,    F_GS+"W/Wh");
        vector_fprinf(FW,    F_GS+"FW/FW");
        vector_fprinf(C_FW,  F_GS+"FW/C_FW");
        vector_fprinf(FWh,   F_GS+"FW/FWh");
        //fprintf Wdot
        vector_fprinf(C_Wdot,  F_GS+"Wdot/C_Wdot");
        //Reduced vector field
        vector_fprinf(fh, F_GS+"rvf/fh");
        //Potential
        vector_fprinf(Unx, F_GS+"Potential/Unx");
        vector_fprinf(Uny, F_GS+"Potential/Uny");
        vector_fprinf(Unz, F_GS+"Potential/Unz");
        vector_fprinf(En,  F_GS+"Potential/En");
        vector_fprinf(Mn, F_GS+"Potential/Mn");
        vector_fprinf(Sn, F_GS+"Potential/Sn");
        vector_fprinf(Un, F_GS+"Potential/Un");
        //DW*vf
        vector_fprinf(C_DWf,  F_GS+"DWf/C_DWf");



        //------------------------------------------
        // Binary
        //------------------------------------------
        writeVOFTS_bin(W,     F_GS+"W/W");
        writeVOFTS_bin(Wt,    F_GS+"W/Wt");
        writeVOFTS_bin(Wh,    F_GS+"W/Wh");
        writeVOFTS_bin(FW,    F_GS+"FW/FW");
        writeVOFTS_bin(C_FW,  F_GS+"FW/C_FW");
        writeVOFTS_bin(FWh,   F_GS+"FW/FWh");
        //fprintf Wdot
        writeVOFTS_bin(C_Wdot,  F_GS+"Wdot/C_Wdot");
        //Reduced vector field
        writeVOFTS_bin(fh, F_GS+"rvf/fh");
        //Potential
        writeVOFTS_bin(Unx, F_GS+"Potential/Unx");
        writeVOFTS_bin(Uny, F_GS+"Potential/Uny");
        writeVOFTS_bin(Unz, F_GS+"Potential/Unz");
        writeVOFTS_bin(En,  F_GS+"Potential/En");
        writeVOFTS_bin(Mn, F_GS+"Potential/Mn");
        writeVOFTS_bin(Sn, F_GS+"Potential/Sn");
        writeVOFTS_bin(Un, F_GS+"Potential/Un");
        //DW*vf
        writeVOFTS_bin(C_DWf,  F_GS+"DWf/C_DWf");


        cout << "  pm. End of printing" <<  " in " << toc() << " s. " << endl;
    }

}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFTS version of the recurrence
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Applying the Normalized-Centered Vector Field on a Ofts object at order m
 *        WARNING: This routine assumes that the potential has been udpated up to order m
 **/
void applyVF(vector<Ofsc>& alpha,
             vector<Oftsc>& W,
             vector<Oftsc>& FW,
             vector<Oftsc>& Unx,  //vector of size OFTS_ORDER
             vector<Oftsc>& Uny,  //vector of size OFTS_ORDER
             vector<Oftsc>& Unz,  //vector of size OFTS_ORDER
             int m)
{
    //-----------------------------------------
    //Configuration variables
    //-----------------------------------------
    //FW[0] = alpha1*px + alpha2*x + alpha3*y
    FW[0].ofts_sfsum_tt(W[3], alpha[0], W[0], alpha[1], W[1],  alpha[2], m);
    //FW[1] = alpha1*py + alpha2*y - alpha3*x
    FW[1].ofts_sfsum_tt(W[4], alpha[0], W[1], alpha[1], W[0], -alpha[2], m);
    //FW[2] = alpha1*pz + alpha2*z
    FW[2].ofts_sfsum_t(W[5], alpha[0], W[2], alpha[1], m);

    //-----------------------------------------
    //Momentum variables
    //-----------------------------------------
    //FW[3] = -alpha2*px + alpha3*py + alpha24 + order m of the derivative of the potential
    FW[3].ofts_sfsum_t(W[3], -alpha[1], W[4], alpha[2], m);
    //At order 0, need to add alpha13
    if(m == 0) FW[3].addCoef(alpha[12], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++) FW[3].ofts_smult_t(Unx[i+1], alpha[5], m);

    //FW[4] = -alpha2*py - alpha3*px + alpha26 + order m of the derivative of the potential
    FW[4].ofts_sfsum_t(W[4], -alpha[1], W[3], -alpha[2], m);
    //At order 0, need to add alpha14
    if(m == 0) FW[4].addCoef(alpha[13], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++) FW[4].ofts_smult_t(Uny[i+1], alpha[5], m);

    //FW[5] = -alpha2*pz  + order m of the derivative of the potential
    FW[5].ofts_smult_t(W[5], -alpha[1], m);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++) FW[5].ofts_smult_t(Unz[i+1], alpha[5], m);
}

/**
 * \brief Applying the Normalized-Centered Vector Field on a Ofts object at order m
 *        WARNING: This routine assumes that the potential has been udpated up to order m
 *        Case with Earth, Moon and Sun potential separated.
 **/
void applyVF(vector<Ofsc>  &alpha,
             vector<Oftsc> &W,
             vector<Oftsc> &FW,
             vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
             vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
             vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
             vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
             int m)
{
    //-----------------------------------------
    //Configuration variables
    //-----------------------------------------
    //FW[0] = alpha1*px + alpha2*x + alpha3*y
    FW[0].ofts_sfsum_tt(W[3], alpha[0], W[0], alpha[1], W[1],  alpha[2], m);
    //FW[1] = alpha1*py + alpha2*y - alpha3*x
    FW[1].ofts_sfsum_tt(W[4], alpha[0], W[1], alpha[1], W[0], -alpha[2], m);
    //FW[2] = alpha1*pz + alpha2*z
    FW[2].ofts_sfsum_t(W[5], alpha[0], W[2], alpha[1], m);

    //-----------------------------------------
    //Momentum variables
    //-----------------------------------------
    //FW[3] = -alpha2*px + alpha3*py + alpha24 + order m of the derivative of the potential
    FW[3].ofts_sfsum_t(W[3], -alpha[1], W[4], alpha[2], m);
    //At order 0, need to add alpha13
    if(m == 0) FW[3].addCoef(alpha[12], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++)
    {
        FW[3].ofts_smult_t(Enx[i+1], alpha[5], m);
        FW[3].ofts_smult_t(Mnx[i+1], alpha[5], m);
        FW[3].ofts_smult_t(Snx[i+1], alpha[5], m);
    }

    //FW[4] = -alpha2*py - alpha3*px + alpha26 + order m of the derivative of the potential
    FW[4].ofts_sfsum_t(W[4], -alpha[1], W[3], -alpha[2], m);
    //At order 0, need to add alpha14
    if(m == 0) FW[4].addCoef(alpha[13], 0, 0);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++)
    {
        FW[4].ofts_smult_t(Eny[i+1], alpha[5], m);
        FW[4].ofts_smult_t(Mny[i+1], alpha[5], m);
        FW[4].ofts_smult_t(Sny[i+1], alpha[5], m);

    }

    //FW[5] = -alpha2*pz  + order m of the derivative of the potential
    FW[5].ofts_smult_t(W[5], -alpha[1], m);
    //Sum of the Earth+Moon+Sun potential up to degree m+1 (order m)
    for(int i = 0; i <= m ; i++)
    {
        FW[5].ofts_smult_t(Enz[i+1], alpha[5], m);
        FW[5].ofts_smult_t(Mnz[i+1], alpha[5], m);
        FW[5].ofts_smult_t(Snz[i+1], alpha[5], m);
    }
}

/**
 *   \brief Initializes the Legendre-derived expansions of one given primary at a given order m.
 *          Note that it is better to reset the order m of En[1] and Enx[2] in this routine.
 *          Since they are linear in the variable W[i], we always have the good value in it.
 **/
void initLegPrim(vector<Oftsc> &Wt,   //vector of size 6
                 vector<Oftsc> &En,   //vector of size OFTS_ORDER
                 vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                 vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                 vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                 vector<Ofsc> &Xe,    //modified position of the Earth
                 Ofsc const &ILe,     //modified inverse orbit radius of the Earth
                 double Ke,
                 Ofsc &AUX,
                 Ofsc &BUX,
                 Ofsc &CUX,
                 Ofsc &DUX,
                 int m)
{
    //En[0] = Ke*ILe
    //------------------------------------------
    AUX.ofs_mult(ILe, Ke+0.0*I);
    if(m == 0) En[0].setCoef(AUX, 0, 0);    //reset here

    //En[1] = Ke*ILe^3*(Xe*W[0] + Ye*W[1]);
    //------------------------------------------
    //BUX = ILe*AUX = Ke*ILe^2
    BUX.ofs_prod(AUX, ILe);
    //CUX = ILe*BUX = Ke*ILe^3
    CUX.ofs_prod(BUX, ILe);
    //AUX = CUX*Xe = Ke*ILe^3*Xe
    AUX.ofs_prod(CUX, Xe[0]);
    //BUX = CUX*Ye = Ke*ILe^3*Ye
    BUX.ofs_prod(CUX, Xe[1]);
    //At every order!
    En[1].ofts_fsum_t(Wt[0], AUX, Wt[1], BUX, m);  //reset here

    //------------------------------------------
    //Derivatives of the  potential
    //------------------------------------------
    //Enx[0] = Eny[0] = Enz[0] = 0
    if(m == 0)
    {
        //Enx[1] = Ke*ILe^3*Xe = AUX
        Enx[1].setCoef(AUX, 0, 0);      //reset here
        //Eny[1] = Ke*ILe^3*Ye = BUX
        Eny[1].setCoef(BUX, 0, 0);      //reset here
    }

    //------------------------------------------
    //Enx
    //------------------------------------------
    //AUX = 1/pe^2
    AUX.ofs_prod(ILe, ILe);
    //BUX = 1/pe^3
    BUX.ofs_prod(AUX, ILe);
    //CUX = 1/pe^4
    CUX.ofs_prod(BUX, ILe);
    //DUX = 1/pe^5
    DUX.ofs_prod(CUX, ILe);

    //BUX = xf/pe^5
    BUX.ofs_prod(DUX, Xe[0]);
    //CUX = xf^2/pe^5
    CUX.ofs_prod(BUX, Xe[0]);

    //AUX = 1/pe^2
    AUX.ofs_prod(ILe, ILe);
    //BUX = 1/pe^3
    BUX.ofs_prod(AUX, ILe);


    //AUX = 3*Ke*CUX - Ke*BUX = 3*Ke*xf^2/pe^5 - Ke/pe^3
    AUX.ofs_fsum(CUX, 3.0*Ke+0.0*I, BUX, -Ke+0.0*I);

    //BUX = xf/pe^5
    BUX.ofs_prod(DUX, Xe[0]);
    //CUX = xf*yf/pe^5
    CUX.ofs_prod(BUX, Xe[1]);
    //CUX= 3*Kf*xf*yf/pe^5
    CUX.ofs_mult(CUX, 3.0*Ke+0.0*I);

    Enx[2].ofts_fsum_t(Wt[0], AUX, Wt[1], CUX, m);  //reset here


    //------------------------------------------
    //Eny
    //------------------------------------------
    //AUX = 1/pe^2
    AUX.ofs_prod(ILe, ILe);
    //BUX = 1/pe^3
    BUX.ofs_prod(AUX, ILe);
    //CUX = 1/pe^4
    CUX.ofs_prod(BUX, ILe);
    //DUX = 1/pe^5
    DUX.ofs_prod(CUX, ILe);

    //BUX = yf/pe^5
    BUX.ofs_prod(DUX, Xe[1]);
    //CUX = yf^2/pe^5
    CUX.ofs_prod(BUX, Xe[1]);

    //AUX = 1/pe^2
    AUX.ofs_prod(ILe, ILe);
    //BUX = 1/pe^3
    BUX.ofs_prod(AUX, ILe);

    //AUX = 3*Ke*CUX - Ke*BUX = 3*Ke*yf^2/pe^5 - Ke/pe^3
    AUX.ofs_fsum(CUX, 3.0*Ke+0.0*I, BUX, -Ke+0.0*I);

    //BUX = xf/pe^5
    BUX.ofs_prod(DUX, Xe[0]);
    //CUX = xf*yf/pe^5
    CUX.ofs_prod(BUX, Xe[1]);
    //CUX= 3*Kf*xf*yf/pe^5
    CUX.ofs_mult(CUX, 3.0*Ke+0.0*I);

    Eny[2].ofts_fsum_t(Wt[0], CUX, Wt[1], AUX, m);  //reset here


    //---------
    //Enz
    //---------
    //AUX = 1/pe^2
    AUX.ofs_prod(ILe, ILe);
    //BUX = 1/pe^3
    BUX.ofs_prod(AUX, ILe);
    //BUX = -Ke/pe^3
    BUX.ofs_mult(BUX, -Ke+0.0*I);

    Enz[2].ofts_mult_t(Wt[2], BUX, m);    //reset here
}


/**
 *   \brief Initializes the Legendre-derived expansions of the Primaries at a given order m.
 *          Note that it is better to reset the order m of En[1] and Enx[2] in this routine.
 *          Since they are linear in the variable W[i], we always have the good value in it.
 **/
void initLegendrePoly( vector<Oftsc> &Wt,   //vector of size 6
                       vector<Oftsc> &En,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                       vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                       vector<Ofsc>  &Xe,   //modified position of the Earth
                       vector<Ofsc>  &Xm,   //modified position of the Moon
                       vector<Ofsc>  &Xs,   //modified position of the Sun
                       Ofsc &ILe,           //modified inverse orbit radius of the Earth
                       Ofsc &ILm,           //modified inverse orbit radius of the Moon
                       Ofsc &ILs,           //modified inverse orbit radius of the Sun
                       QBCP_L& qbcp_l,
                       int m)
{
    //------------------------------------------
    // Initialization
    //------------------------------------------
    double gamma =  qbcp_l.cs.gamma;
    double Ke = qbcp_l.us.me/pow(gamma, 3.0);
    double Km = qbcp_l.us.mm/pow(gamma, 3.0);
    double Ks = qbcp_l.us.ms/pow(gamma, 3.0);

    //------------------------------------------
    //Temporary variables
    //------------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);
    Ofsc DUX(OFS_ORDER);

    //------------------------------------------
    //Earth potential
    //------------------------------------------
    initLegPrim(Wt, En,  Enx,  Eny, Enz, Xe, ILe, Ke, AUX, BUX, CUX, DUX,  m);

    //------------------------------------------
    //Moon potential
    //------------------------------------------
    initLegPrim(Wt, Mn,  Mnx,  Mny, Mnz, Xm, ILm, Km, AUX, BUX, CUX, DUX, m);

    //------------------------------------------
    //Sun potential
    //------------------------------------------
    initLegPrim(Wt, Sn,  Snx,  Sny, Snz, Xs, ILs, Ks, AUX, BUX, CUX, DUX, m);

    //--------------------------------
    //Earth+Moon+Sum
    //--------------------------------
    for(int i =  0; i <= 1; i++) //sum up to i = 1
    {
        Un[i].ofts_mult_u(En[i],   1.0+0.0*I, m);        //reset here
        Un[i].ofts_smult_u(Mn[i],  1.0+0.0*I, m);
        Un[i].ofts_smult_u(Sn[i],  1.0+0.0*I, m);
    }



    for(int i =  0; i <= 2; i++) //sum up to i = 2
    {
        Unx[i].ofts_mult_u(Enx[i],   1.0+0.0*I, m);     //reset here
        Unx[i].ofts_smult_u(Mnx[i],  1.0+0.0*I, m);
        Unx[i].ofts_smult_u(Snx[i],  1.0+0.0*I, m);

        Uny[i].ofts_mult_u(Eny[i],   1.0+0.0*I, m);     //reset here
        Uny[i].ofts_smult_u(Mny[i],  1.0+0.0*I, m);
        Uny[i].ofts_smult_u(Sny[i],  1.0+0.0*I, m);

        Unz[i].ofts_mult_u(Enz[i],   1.0+0.0*I, m);     //reset here
        Unz[i].ofts_smult_u(Mnz[i],  1.0+0.0*I, m);
        Unz[i].ofts_smult_u(Snz[i],  1.0+0.0*I, m);
    }
}

/**
 *   \brief Apply the Legendre-derived recurrence for the potential of one given primary.
 *          No reset here, since we may add other terms afterwards.
 **/
void applyLegRec(Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                 vector<Oftsc> &En,   //vector of size OFTS_ORDER
                 Ofsc &ILe,           //modified inverse orbit radius of the Primary
                 Oftsc &xse,          //=xe*xt+ye*yt @order m
                 Ofsc &AUX,           //spare OFS object
                 Ofsc &temp,          //spare OFS object
                 int m,               //current order of the W expansion
                 int k)               //current En[k] that needs to be updated
{
    //AUX = ILe^2*(2*k-1)/k
    AUX.ofs_mprod(ILe, ILe, (2.0*k - 1.0)/(double)k+0.0*I);
    //smprod at order m
    En[k].ofts_smprod_t(xse, En[k-1], AUX, m, temp);
    //AUX = ILe^2*(2*k-1)/k
    AUX.ofs_mprod(ILe, ILe, -(k-1.0)/(double)k+0.0*I);
    //smprod at order m
    En[k].ofts_smprod_t(Rho2, En[k-2], AUX, m, temp);
}

/**
 *   \brief Apply the Legendre-derived recurrence for the derivatives of the potential of one given primary.
 *          No reset here, since we may add other terms afterwards.
 **/
void applyLegRecDer(vector<Oftsc> &Wt,   //tilde state vector (xt, yt,..)
                    Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                    vector<Oftsc> &En,   //vector of size OFTS_ORDER
                    vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                    vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                    vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                    vector<Ofsc>  &Xe,   //modified position of the primary
                    Ofsc &ILe,           //modified inverse orbit radius of the Primary
                    Oftsc &xse,          //=xe*xt+ye*yt @order m
                    Ofsc &AUX,           //spare OFS object
                    Ofsc &BUX,           //spare OFS object
                    Ofsc &temp,          //spare OFS object
                    int m,               //current order of the W expansion
                    int k)               //current En[k] that needs to be updated
{
    //AUX = ILe^2*(2*k-1)/k
    AUX.ofs_mprod(ILe, ILe, (2.0*k - 1.0)/k+0.0*I);

    //-----------------------------------------------
    // First part of the recursive equations
    //-----------------------------------------------
    //Enx
    BUX.ofs_prod(AUX, Xe[0]);  //BUX = xe*ILe^2*(2*k-1)/k
    Enx[k].ofts_smult_t(En[k-1], BUX, m);
    Enx[k].ofts_smprod_t(xse, Enx[k-1], AUX, m, temp);
    //Eny
    BUX.ofs_prod(AUX, Xe[1]);  //BUX = ye*ILe^2*(2*k-1)/k
    Eny[k].ofts_smult_t(En[k-1], BUX, m);
    Eny[k].ofts_smprod_t(xse, Eny[k-1], AUX, m, temp);
    //Enz
    Enz[k].ofts_smprod_t(xse, Enz[k-1], AUX, m, temp);

    //-----------------------------------------------
    // Second part of the recursive equations
    //-----------------------------------------------
    //AUX = -ILe^2*(k-1)/k
    AUX.ofs_mprod(ILe, ILe, -(k - 1.0)/k+0.0*I);
    //BUX = 2*AUX
    BUX.ofs_mult(AUX, 2.0+0.0*I);

    //Enx
    Enx[k].ofts_smprod_t(En[k-2], Wt[0]    , BUX, m, temp);
    Enx[k].ofts_smprod_t(Rho2   , Enx[k-2],  AUX, m, temp);

    //Eny
    Eny[k].ofts_smprod_t(En[k-2], Wt[1]    , BUX, m, temp);
    Eny[k].ofts_smprod_t(Rho2   , Eny[k-2],  AUX, m, temp);
    //Enz
    Enz[k].ofts_smprod_t(En[k-2], Wt[2]    , BUX, m, temp);
    Enz[k].ofts_smprod_t(Rho2   , Enz[k-2],  AUX, m, temp);
}


/**
 *   \brief Update the order m of the kth coefficient in the Legendre-derived objects En, Mn, etc (k < m),
 *          Provided that the (k-2t)h and (k-1)th coefficient have been updated.
 **/
void updateLegPoly(vector<Oftsc> &Wt,   //vector of size 6: tilde state vector (xt, yt,..)
                   Oftsc &Rho2,         //Ofts containing Wt[0]^2+Wt[1]^2+Wt[2]^2
                   vector<Oftsc> &En,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                   vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                   vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                   vector<Ofsc> &Xe,    //modified position of the Earth
                   vector<Ofsc> &Xm,    //modified position of the Moon
                   vector<Ofsc> &Xs,    //modified position of the Sun
                   Ofsc &ILe,           //modified inverse orbit radius of the Earth
                   Ofsc &ILm,           //modified inverse orbit radius of the Moon
                   Ofsc &ILs,           //modified inverse orbit radius of the Sun
                   Oftsc &xse,          //=xe*xt+ye*yt
                   Oftsc &xsm,          //=xm*xt+ym*yt
                   Oftsc &xss,          //=xs*xt+ys*yt
                   Ofsc &AUX,           //Temporary variables
                   Ofsc &BUX,           //Temporary variables
                   Ofsc &temp,          //Temporary variables
                   int m,               //current order of the W expansion
                   int k)               //current En[k] (and Enx[k+1]) that needs to be updated
{
    //At this step:
    // -  En[0, .., k-1] are at order m
    // -  Enx[0, .., k]  are at order m
    //We want to update:
    // -  En[k] at order m
    // -  Enx[k+1] at order m

    //------------------------------------------
    //For the Earth potential
    //------------------------------------------
    //Update En[k] @order m
    applyLegRec(Rho2, En, ILe, xse, AUX, temp, m, k);
    //Update Enx, Eny, Enz[k+1] @order m
    applyLegRecDer(Wt, Rho2, En, Enx, Eny, Enz, Xe, ILe, xse, AUX, BUX, temp, m, k+1);


    //------------------------------------------
    //For the Moon potential
    //------------------------------------------
    //Update Mn[k] @order m
    applyLegRec(Rho2, Mn, ILm, xsm, AUX, temp, m, k);
    //Update Mnx, Mny, Mnz[k+1] @order m
    applyLegRecDer(Wt, Rho2, Mn, Mnx, Mny, Mnz, Xm, ILm, xsm, AUX, BUX, temp, m, k+1);

    //------------------------------------------
    //For the Sun potential
    //------------------------------------------
    //Update Sn[k] @order m
    applyLegRec(Rho2, Sn, ILs, xss, AUX, temp, m, k);
    //Update Snx, Sny, Snz[k+1] @order m
    applyLegRecDer(Wt, Rho2, Sn, Snx, Sny, Snz, Xs, ILs, xss, AUX, BUX, temp, m, k+1);

    //------------------------------------------
    //For the Earth+Moon+Sun
    //------------------------------------------
    Un[k].ofts_smult_u(En[k],  1.0+0.0*I, m);
    Un[k].ofts_smult_u(Mn[k],  1.0+0.0*I, m);
    Un[k].ofts_smult_u(Sn[k],  1.0+0.0*I, m);

    Unx[k+1].ofts_smult_u(Enx[k+1],  1.0+0.0*I, m);
    Unx[k+1].ofts_smult_u(Mnx[k+1],  1.0+0.0*I, m);
    Unx[k+1].ofts_smult_u(Snx[k+1],  1.0+0.0*I, m);

    Uny[k+1].ofts_smult_u(Eny[k+1],  1.0+0.0*I, m);
    Uny[k+1].ofts_smult_u(Mny[k+1],  1.0+0.0*I, m);
    Uny[k+1].ofts_smult_u(Sny[k+1],  1.0+0.0*I, m);

    Unz[k+1].ofts_smult_u(Enz[k+1],  1.0+0.0*I, m);
    Unz[k+1].ofts_smult_u(Mnz[k+1],  1.0+0.0*I, m);
    Unz[k+1].ofts_smult_u(Snz[k+1],  1.0+0.0*I, m);

    //Now, En[k]...Un[k] are at order m
    //And  Enx[k+1]...Unz[k+1] are at order m
}



/**
 *   \brief Update the whole Legendre-derived objects at order m > 1.
 **/
void updatePotentialExpansion(vector<Oftsc> &Wt,    //vector of size 6
                              Oftsc &Rho2,         //Ofts containing W[0]^2+W[1]^2+W[2]^2
                              vector<Oftsc> &En,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                              vector<Oftsc> &Enx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Eny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Enz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mnx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Mnz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Snx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Sny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Snz,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Unx,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Uny,  //vector of size OFTS_ORDER
                              vector<Oftsc> &Unz,  //vector of size OFTS_ORDER
                              vector<Ofsc> &Xe,    //modified position of the Earth
                              vector<Ofsc> &Xm,    //modified position of the Moon
                              vector<Ofsc> &Xs,    //modified position of the Sun
                              Ofsc &ILe,           //modified inverse orbit radius of the Earth
                              Ofsc &ILm,           //modified inverse orbit radius of the Moon
                              Ofsc &ILs,           //modified inverse orbit radius of the Sun
                              Oftsc &xse,          //=xe*xt+ye*yt
                              Oftsc &xsm,          //=xm*xt+ym*yt
                              Oftsc &xss,          //=xs*xt+ys*yt
                              QBCP_L& qbcp_l,       //current QBCP
                              int m)                  //targeted order of the W expansion
{
    //------------------------------------------
    //Update Rho2 at order m
    //------------------------------------------
    Rho2.ofts_sprod(Wt[0], Wt[0], m);
    Rho2.ofts_sprod(Wt[1], Wt[1], m);
    Rho2.ofts_sprod(Wt[2], Wt[2], m);

    //------------------------------------------
    //Update En[0, 1], Mn[0, 1] and Sn[0, 1] at order m-1
    //Update Enx[0, 1, 2],... and Snz[0, 1, 2] at order m-1
    //Note: the update is at order m-1 because En[1] and Enx[2] are linear in the W[i]
    //------------------------------------------
    if(m == 2) initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, qbcp_l, m-2);
    initLegendrePoly(Wt, En, Mn, Sn, Un, Enx, Eny, Enz, Mnx, Mny, Mnz, Snx, Sny, Snz, Unx, Uny, Unz, Xe, Xm, Xs, ILe, ILm, ILs, qbcp_l, m-1);

    if(m == 2)
    {
        //Update xse @ order m-2
        xse.ofts_sfsum_t(Wt[0], Xe[0], Wt[1], Xe[1], m-2);
        //Update xsm @ order m-2
        xsm.ofts_sfsum_t(Wt[0], Xm[0], Wt[1], Xm[1], m-2);
        //Update xss @ order m-2
        xss.ofts_sfsum_t(Wt[0], Xs[0], Wt[1], Xs[1], m-2);
    }

    //Update xse @ order m-1
    xse.ofts_sfsum_t(Wt[0], Xe[0], Wt[1], Xe[1], m-1);
    //Update xsm @ order m-1
    xsm.ofts_sfsum_t(Wt[0], Xm[0], Wt[1], Xm[1], m-1);
    //Update xss @ order m-1
    xss.ofts_sfsum_t(Wt[0], Xs[0], Wt[1], Xs[1], m-1);

    //At this step:
    // - W[i] are at order m-1
    // - Rho2 is at order m
    // - En[0, 1] are at order m-1
    // - En[2,...,m-1] are at order m-1
    // - Enx[0, 1, 2] are at order m-1
    // - Enx[3, ..., m] are at order m-1

    //------------------------------------------
    //Temporary variables
    //------------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);
    Ofsc temp(OFS_ORDER);

    //-----------------------------------------------------------------
    //Update En[2,...,m-1], Mn[2,...,m-1]  and Sn[2,...,m-1] at order m
    //Update Enx[3,...,m],  Mnx[3,...,m-1] and Snx[3,...,m]  at order m
    //WARNING: only the contribution of order m-1 to order m is computed!
    //-----------------------------------------------------------------
    for(int k = 2; k < m; k++)
    {
        //cout << "updatePotentialExpansion. Order is " << m  << endl;
        updateLegPoly(Wt, Rho2, En,
                      Mn, Sn, Un,
                      Enx, Eny, Enz,
                      Mnx, Mny, Mnz,
                      Snx, Sny, Snz,
                      Unx, Uny, Unz,
                      Xe, Xm, Xs,
                      ILe, ILm, ILs,
                      xse, xsm, xss,
                      AUX, BUX, temp,
                      m, k);
    }

    //At this step:
    // - W[i] are at order m-1
    // - Rho2 is at order m
    // - En[0, 1] are at order m-1
    // - En[2,...,m-1] are at order m
    // - Enx[0, 1, 2] are at order m-1
    // - Enx[3, ..., m] are at order m

    //-----------------------------------------------------------------
    //Update En[m], Mn[m] and Sn[m] at order m (first non null coefficient)
    //Update Enx[m+1], Mnx[m+1] and Snx[m+1] at order m (first non null coefficient)
    //WARNING: only the contribution of order m-1 to order m is computed!
    //-----------------------------------------------------------------
    updateLegPoly(Wt, Rho2, En,
                  Mn, Sn, Un,
                  Enx, Eny, Enz,
                  Mnx, Mny, Mnz,
                  Snx, Sny, Snz,
                  Unx, Uny, Unz,
                  Xe, Xm, Xs,
                  ILe, ILm, ILs,
                  xse, xsm, xss,
                  AUX, BUX, temp,
                  m, m);

    cout << "updatePotentialExpansion. end of updateLegPoly @order m." << " in " << toc() << " s. " << endl;

    //At this step:
    // - W[i] are at order m-1
    // - Rho2 is at order m
    // - En[0, 1] are at order m-1
    // - En[2,...,m] are at order m
    // - Enx[0, 1, 2] are at order m-1
    // - Enx[3, ..., m+1] are at order m

    //i.e: the contribution of the potential to the order m of the vector field F(W(s, t)) is the sum of the order m of En[2,...,m]
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluate the potential (for testing)
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test the validity of the expansions of the potential (En, Sn, Mn, Un) and the vector field (FW).
 **/
void pmTestPrec(double smax)
{
    //------------------------------------------
    // Strings
    //------------------------------------------
    string F_GS    = SEML.cs.F_PMS;
    string F_PLOT  = SEML.cs.F_PLOT;
    string F_COC   = SEML.cs.F_COC;

    //----------------------------------------------------------------------------------------
    // Initialization of the COC
    //----------------------------------------------------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(6,6);
    //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    matrix<Ofsc> Q(6,6);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(6);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential:
    //Xf = True Xf - V, for f = e, m, s
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);
    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    //ILf = 1.0/||Xf||, for f = e, m, s
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);
    //Complement
    matrix<Ofsc> PC(6,6);    //PC=P(theta)*C
    matrix<Ofsc> PCdot(6,6); //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> CQ(6,6);    //CQ=Cinv*Pinv=inv(PC)
    vector<Ofsc> Vdot(6);   //Vdot=dot(V)
    //---------------
    //Init routine
    //---------------
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //----------------------------------------------------------------------------------------
    // Initialisation of OFTS objects
    //----------------------------------------------------------------------------------------
    //Manifolds and vector fields
    vector<Oftsc> W(6);
    vector<Oftsc> Wt(6);
    vector<Oftsc> Wh(6);
    vector<Oftsc> FW(6);

    //Potential
    vector<Oftsc> En(OFTS_ORDER+2);     //Earth
    vector<Oftsc> Sn(OFTS_ORDER+2);     //Moon
    vector<Oftsc> Mn(OFTS_ORDER+2);     //Sun
    vector<Oftsc> Un(OFTS_ORDER+2);     //Earth+Moon+Sun

    readVOFTS_bin(W,    F_GS+"W/W",             OFS_ORDER);
    readVOFTS_bin(Wt,   F_GS+"W/Wt",            OFS_ORDER);
    readVOFTS_bin(Wh,   F_GS+"W/Wh",            OFS_ORDER);
    readVOFTS_bin(FW,   F_GS+"FW/C_FW",         OFS_ORDER);
    readVOFTS_bin(En,   F_GS+"Potential/En",    OFS_ORDER);
    readVOFTS_bin(Mn,   F_GS+"Potential/Mn",    OFS_ORDER);
    readVOFTS_bin(Sn,   F_GS+"Potential/Sn",    OFS_ORDER);
    readVOFTS_bin(Un,   F_GS+"Potential/Un",    OFS_ORDER);

    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the obtained PM             " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(6);

    //----------------------------------------------------------------------------------------
    //For plotting purposes
    //----------------------------------------------------------------------------------------
    int plotSize = 500;
    double Wc[plotSize];       //PM
    double VF1c[6][plotSize];  //Vector field
    double VF1cmax[plotSize];  //Vector field
    double Uc[plotSize];       //Potential

    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    gnuplot_cmd(h1, "set logscale y");
    gnuplot_cmd(h1, "set format y \"1e\%%L\"");
    gnuplot_cmd(h2, "set logscale y");
    gnuplot_cmd(h2, "set format y \"1e\%%L\"");

    gnuplot_cmd(h1, "set grid");
    gnuplot_cmd(h2, "set grid");

    gnuplot_cmd(h1, "set key right bottom");
    gnuplot_cmd(h2, "set key right bottom");

    //----------------------------------------------------------------------------------------
    //For evaluation purposes
    //----------------------------------------------------------------------------------------
    cdouble Xeval[4], st0[4];
    double st0v[plotSize];
    int keyMap[5];
    //Vector field
    cdouble VF1[6], VF2[6];
    //Potential
    cdouble Ea, M, S, U;
    cdouble Ea2, M2, S2, U2;

    //Fill the grid
    for(int i = 1; i <= plotSize; i++) st0v[i-1] = (double) i*smax/((double)plotSize);

    //Keymap
    keyMap[0] = 2;
    keyMap[1] = 5;
    keyMap[2] = 8;
    keyMap[3] = 10;
    keyMap[4] = OFTS_ORDER;

    //----------------------------------------------------------------------------------------
    //Evaluations: expansion vs real implementation
    //----------------------------------------------------------------------------------------
    int color = 1;
    ifstream readStream;
    string ss1;
    int m;

    cdouble z0[6];
    vector<Ofsc> z0t(6);

    for(int ml = 0 ; ml < 4; ml++)
    {
        m = keyMap[ml];

        for(int k = 0 ; k < plotSize; k++)
        {
            //Initialization of the configuration
            for(int p = 0; p < 4; p++) st0[p] = st0v[k]+0.0*I;
            //st0s = "1e-4";
            //"Realification": Xeval = REAL(st0)
            Xeval[0] = 1.0/sqrt(2)*( st0[0]   - st0[2]*I);
            Xeval[2] = 1.0/sqrt(2)*( st0[2]   - st0[0]*I);
            Xeval[1] = 1.0/sqrt(2)*( st0[1]   - st0[3]*I);
            Xeval[3] = 1.0/sqrt(2)*( st0[3]   - st0[1]*I);

            //------------------------------------------
            // Evaluate Wc = |W(Xeval, 0)|
            //------------------------------------------

            for(int p = 0; p < 6; p++)
            {
                W[p].evaluate(Xeval, z0t[p], m, OFS_ORDER);
                z0[p] = z0t[p].evaluate(0.0, OFS_ORDER);
            }
            Wc[k] = ENorm(z0,3)*SEML.cs.cr3bp.L*SEML.cs.gamma;



            //Evaluate the potential: expansion vs real implementation
            evaluatePotential(Xeval, Ea, M, S, U,
                              Ea2, M2, S2, U2,
                              En, Mn, Sn, Un,
                              Wt, Xe,  Xm,  Xs,
                              0.0, SEML, m);

            //Evaluate the vector field: expansion vs real implementation
            evaluateVectorField(Xeval, VF1, VF2, W, FW, 1.0, SEML, m);

            //------------------------------------------
            //For plotting purposes
            //------------------------------------------
            VF1cmax[k] = cabs(VF1[0] - VF2[0]);
            for(int p = 0; p< Csts::NV; p++)
            {
                VF1c[p][k] = cabs(VF1[p] - VF2[p]);
                if(VF1cmax[k] < VF1c[p][k]) VF1cmax[k] = VF1c[p][k];
            }

            Uc[k] = cabs(U2 - U);
        }

        //----------------------------------------------------------------------------------------
        // Plotting
        //----------------------------------------------------------------------------------------
        ss1 = static_cast<ostringstream*>( &(ostringstream() << m) )->str();
        gnuplot_set_xlabel(h1, (char*)"|W(s,0)| [-]");
        gnuplot_set_ylabel(h1, (char*)"Error on the potential [-]");
        gnuplot_cmd(h1,"set yrange [1e-16:1e-0] ");
        gnuplot_plot_xy(h1, Wc, Uc, plotSize, (char*)("Order. "+ss1).c_str(), "lines", "4", "3", color);
        //In txt files
        gnuplot_fplot_xy(Wc, Uc, plotSize, (F_PLOT+"PotentialvsW_"+ss1+".txt").c_str());

        gnuplot_set_xlabel(h2, (char*)"|W(s,0)| [-]");
        gnuplot_set_ylabel(h2, (char*)"Error on the vector field [-]");
        gnuplot_cmd(h1,"set yrange [1e-16:1e-0] ");
        gnuplot_plot_xy(h2, Wc, VF1cmax, plotSize, (char*)("Order. "+ss1).c_str(), "lines", "4", "3", color++);
        //In txt files
        gnuplot_fplot_xy(Wc, VF1cmax, plotSize, (F_PLOT+"VFvsW_"+ss1+".txt").c_str());
    }


    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //----------------------------------------------------------------------------------------
    //Save in EPS format
    //----------------------------------------------------------------------------------------
    //gnuplot_cmd(h1,"set terminal postscript eps size 3.5,2.62 enhanced solid color font 'Helvetica,20' linewidth 2");
    //gnuplot_cmd(h1, "set terminal postscript eps solid color enhanced");

    if(SEML.model == Csts::CRTBP)
    {
        gnuplot_cmd(h1,"set terminal postscript eps enhanced dashed color font 'Helvetica,20' dl 4 linewidth 3");
        gnuplot_cmd(h2,"set terminal postscript eps enhanced dashed color font 'Helvetica,20' dl 4 linewidth 3");
    }
    else
    {
        gnuplot_cmd(h1,"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
        gnuplot_cmd(h2,"set terminal postscript eps enhanced solid color font 'Helvetica,20' linewidth 3");
    }

    gnuplot_cmd(h1, ("set output \""+F_PLOT+"PotentialvsW"+".eps\"").c_str());
    gnuplot_cmd(h1, "replot");
    gnuplot_cmd(h2, ("set output \""+F_PLOT+""+"VFvsW"+".eps\"").c_str());
    gnuplot_cmd(h2, "replot");

    gnuplot_close(h1);
    gnuplot_close(h2);

}



/**
 *   \brief Evaluate the different potentials at order m at time n*t: expansion vs real implementation
 **/
void evaluatePotential(cdouble X[],
                       cdouble& Ez,
                       cdouble& Mz,
                       cdouble& Sz,
                       cdouble& Uz,
                       cdouble& Ez2,
                       cdouble& Mz2,
                       cdouble& Sz2,
                       cdouble& Uz2,
                       vector<Oftsc> &En,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Mn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Sn,   //vector of size OFTS_ORDER
                       vector<Oftsc> &Un,   //vector of size OFTS_ORDER
                       vector<Oftsc> &W,   //vector of size OFTS_ORDER
                       vector<Ofsc> &Xe,    //modified position of the Earth
                       vector<Ofsc> &Xm,    //modified position of the Moon
                       vector<Ofsc> &Xs,    //modified position of the Sun
                       double t,
                       QBCP_L& qbcp_l,
                       int m)
{
    double ms = qbcp_l.us.ms;
    double me = qbcp_l.us.me;
    double mm = qbcp_l.us.mm;
    double n  = qbcp_l.us.n;
    double gamma  = qbcp_l.cs.gamma;
    cdouble  xxe, yye, zze, qpe;
    Ofsc xt(OFS_ORDER);
    Ofsc yt(OFS_ORDER);
    Ofsc zt(OFS_ORDER);
    Ofsc temp(OFS_ORDER);
    Ofsc Fz(OFS_ORDER);


    //--------------------------------------------------------
    // From expansions
    //--------------------------------------------------------
    //Earth potential
    Fz.zero();
    for(int i = 0; i <= m; i++) //evaluate all terms up to indix m
    {
        En[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp,  1.0+0.0*I);
    }
    Ez = Fz.evaluate(n*t);


    //Moon potential
    Fz.zero();
    for(int i = 0; i <= m; i++) //evaluate all terms up to indix m
    {
        Mn[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Mz = Fz.evaluate(n*t);


    //Sun potential
    Fz.zero();
    for(int i = 0; i <= m; i++) //evaluate all terms up to indix m
    {
        Sn[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Sz = Fz.evaluate(n*t);



    //Sum potential
    Fz.zero();
    for(int i = 0; i <= m; i++) //evaluate all terms up to indix m
    {
        Un[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Uz = Fz.evaluate(n*t);//Ez+Mz+Sz;//


    //--------------------------------------------------------
    // "True potential"
    //--------------------------------------------------------
    //Evaluate the position
    W[0].evaluate(X, xt, m, OFS_ORDER);
    W[1].evaluate(X, yt, m, OFS_ORDER);
    W[2].evaluate(X, zt, m, OFS_ORDER);

    //-------------------
    //Earth potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xe[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xe[1].evaluate(n*t));
    zze = (zt.evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);
    Ez2 = me/(pow(gamma, 3.0)*qpe);


    //-------------------
    //Moon potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xm[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xm[1].evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);
    Mz2 =  mm/(pow(gamma, 3.0)*qpe);


    //-------------------
    //Sun potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xs[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xs[1].evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);
    Sz2 =  ms/(pow(gamma, 3.0)*qpe);

    //Sum
    Uz2 = Ez2+Mz2+Sz2;


    //--------------------------------------------------------
    // On screen
    //--------------------------------------------------------
    cout << setprecision(2);
    cout << "------------------------------------------------------------------------------" << endl;
    cout << "Evaluate the Potential @ order " << m << "" << endl;
    cout << " |W(1:3)| ~ [" << cabs(xt.evaluate(n*t)) << ", "<< cabs(yt.evaluate(n*t)) << ", "<< cabs(zt.evaluate(n*t)) << "]" << endl;
    cout << " Object  |         True           |       From Expansions    |   Delta (norm) " << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " Earth   | "
         << creal(Ez2) << "  " << cimag(Ez2)
         << "   |   " << creal(Ez) << "  " << cimag(Ez)
         << "   |   " << cabs(Ez2-Ez) << endl;
    cout << " Moon    | "
         << creal(Mz2) << "  " << cimag(Mz2)
         << "   |   " << creal(Mz) << "  " << cimag(Mz)
         << "   |   " << cabs(Mz2-Mz) << endl;
    cout << " Sun     | "
         << creal(Sz2) << "  " << cimag(Sz2)
         << "   |   " << creal(Sz) << "  " << cimag(Sz)
         << "   |   " << cabs(Sz2-Sz) << endl;
    cout << " Sum     | "
         << creal(Uz2) << "  " << cimag(Uz2)
         << "   |   " << creal(Uz) << "  " << cimag(Uz)
         << "   |   " << cabs(Uz2-Uz) << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << setprecision(15);

}


/**
 *   \brief Evaluate the derivatives of the potentials at order m at time n*t: expansion vs real implementation
 **/
void evaluatePotentialDer(cdouble X[],
                          cdouble& Ex, cdouble& Mx, cdouble& Sx, cdouble& Ux,
                          cdouble& Ey, cdouble& My, cdouble& Sy, cdouble& Uy,
                          cdouble& Ez, cdouble& Mz, cdouble& Sz, cdouble& Uz,
                          cdouble& Ex2, cdouble& Mx2, cdouble& Sx2, cdouble& Ux2,
                          cdouble& Ey2, cdouble& My2, cdouble& Sy2, cdouble& Uy2,
                          cdouble& Ez2, cdouble& Mz2, cdouble& Sz2, cdouble& Uz2,
                          vector<Oftsc> &Enx, vector<Oftsc> &Mnx, vector<Oftsc> &Snx, vector<Oftsc> &Unx,
                          vector<Oftsc> &Eny, vector<Oftsc> &Mny, vector<Oftsc> &Sny, vector<Oftsc> &Uny,
                          vector<Oftsc> &Enz, vector<Oftsc> &Mnz, vector<Oftsc> &Snz, vector<Oftsc> &Unz,
                          vector<Oftsc> &W,   vector<Ofsc> &Xe,   vector<Ofsc> &Xm,   vector<Ofsc> &Xs,
                          double t,
                          QBCP_L& qbcp_l,
                          int m)
{
    double ms = qbcp_l.us.ms;
    double me = qbcp_l.us.me;
    double mm = qbcp_l.us.mm;
    double n  = qbcp_l.us.n;
    double gamma  = qbcp_l.cs.gamma;
    cdouble  xxe, yye, zze, qpe;
    Ofsc xt(OFS_ORDER);
    Ofsc yt(OFS_ORDER);
    Ofsc zt(OFS_ORDER);
    Ofsc temp(OFS_ORDER);
    Ofsc Fz(OFS_ORDER);

    //--------------------------------------------------------
    // From expansions
    //--------------------------------------------------------
    //-------------------
    //Earth potential
    //-------------------
    //Enx
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Enx[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Ex = Fz.evaluate(n*t);

    //Eny
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Eny[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Ey = Fz.evaluate(n*t);

    //Enz
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Enz[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp,  1.0+0.0*I);
    }
    Ez = Fz.evaluate(n*t);

    //-------------------
    //Moon potential
    //-------------------
    //Mnx
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Mnx[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Mx = Fz.evaluate(n*t);

    //Mny
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Mny[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    My = Fz.evaluate(n*t);

    //Mnz
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Mnz[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Mz = Fz.evaluate(n*t);


    //-------------------
    //Sun potential
    //-------------------
    //Snx
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Snx[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Sx = Fz.evaluate(n*t);

    //Sny
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Sny[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Sy = Fz.evaluate(n*t);

    //Snz
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Snz[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Sz = Fz.evaluate(n*t);



    //-------------------
    //Sum potential
    //-------------------
    //Unx
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Unx[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Ux = Fz.evaluate(n*t);

    //Uny
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m
    {
        Uny[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Uy = Fz.evaluate(n*t);

    //Unz
    //-------------------
    Fz.zero();
    for(int i = 0; i <= m+1; i++) //evaluate all terms up to indix m+1
    {
        Unz[i].evaluate(X, temp, m, OFS_ORDER); //evaluate up to order m
        Fz.ofs_smult(temp, 1.0+0.0*I);
    }
    Uz = Fz.evaluate(n*t);


    //--------------------------------------------------------
    // "True potential"
    //--------------------------------------------------------
    //Evaluate the position
    W[0].evaluate(X, xt, m, OFS_ORDER);
    W[1].evaluate(X, yt, m, OFS_ORDER);
    W[2].evaluate(X, zt, m, OFS_ORDER);

    //-------------------
    //Earth potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xe[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xe[1].evaluate(n*t));
    zze = (zt.evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);

    Ex2 =  -me/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*xxe;
    Ey2 =  -me/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*yye;
    Ez2 =  -me/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*zze;


    //-------------------
    //Moon potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xm[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xm[1].evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);

    Mx2 =  -mm/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*xxe;
    My2 =  -mm/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*yye;
    Mz2 =  -mm/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*zze;


    //-------------------
    //Sun potential
    //-------------------
    //Le^2 = xe*xe + ye*ye
    xxe = (xt.evaluate(n*t) - Xs[0].evaluate(n*t));
    yye = (yt.evaluate(n*t) - Xs[1].evaluate(n*t));
    qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);

    Sx2 =  -ms/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*xxe;
    Sy2 =  -ms/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*yye;
    Sz2 =  -ms/pow(gamma, 3.0)*cpow(qpe, -3.0+0.0*I)*zze;


    //-------------------
    //Sum potential
    //-------------------
    Ux2 = Ex2 + Mx2 + Sx2;
    Uy2 = Ey2 + My2 + Sy2;
    Uz2 = Ez2 + Mz2 + Sz2;



    //--------------------------------------------------------
    // On screen
    //--------------------------------------------------------
    cout << setprecision(2);
    cout << "------------------------------------------------------------------------------" << endl;
    cout << "Evaluate the potential partial derivatives @ order " << m << "" << endl;
    cout << " Object     |        True            |       From Expansions    |   Delta  " << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " Earth (dx) | "
         << creal(Ex2) << "  " << cimag(Ex2)
         << "   |   " << creal(Ex) << "  " << cimag(Ex)
         << "   |   " << cabs(Ex2-Ex) << endl;
    cout << " Moon  (dx) | "
         << creal(Mx2) << "  " << cimag(Mx2)
         << "   |   " << creal(Mx) << "  " << cimag(Mx)
         << "   |   " << cabs(Mx2-Mx) << endl;
    cout << " Sun   (dx) | "
         << creal(Sx2) << "  " << cimag(Sx2)
         << "   |   " << creal(Sx) << "  " << cimag(Sx)
         << "   |   " << cabs(Sx2-Sx) << endl;
    cout << " Sum   (dx) | "
         << creal(Ux2) << "  " << cimag(Ux2)
         << "   |   " << creal(Ux) << "  " << cimag(Ux)
         << "   |   " << cabs(Ux2-Ux) << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " Earth (dy) | "
         << creal(Ey2) << "  " << cimag(Ey2)
         << "   |   " << creal(Ey) << "  " << cimag(Ey)
         << "   |   " << cabs(Ey2-Ey) << endl;
    cout << " Moon  (dy) | "
         << creal(My2) << "  " << cimag(My2)
         << "   |   " << creal(My) << "  " << cimag(My)
         << "   |   " << cabs(My2-My) << endl;
    cout << " Sun   (dy) | "
         << creal(Sy2) << "  " << cimag(Sy2)
         << "   |   " << creal(Sy) << "  " << cimag(Sy)
         << "   |   " << cabs(Sy2-Sy) << endl;
    cout << " Sum   (dy) | "
         << creal(Uy2) << "  " << cimag(Uy2)
         << "   |   " << creal(Uy) << "  " << cimag(Uy)
         << "   |   " << cabs(Uy2-Uy) << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " Earth (dz) | "
         << creal(Ez2) << "  " << cimag(Ez2)
         << "   |   " << creal(Ez) << "  " << cimag(Ez)
         << "   |   " << cabs(Ez2-Ez) << endl;
    cout << " Moon  (dz) | "
         << creal(Mz2) << "  " << cimag(Mz2)
         << "   |   " << creal(Mz) << "  " << cimag(Mz)
         << "   |   " << cabs(Mz2-Mz) << endl;
    cout << " Sun   (dz) | "
         << creal(Sz2) << "  " << cimag(Sz2)
         << "   |   " << creal(Sz) << "  " << cimag(Sz)
         << "   |   " << cabs(Sz2-Sz) << endl;
    cout << " Sum   (dz) | "
         << creal(Uz2) << "  " << cimag(Uz2)
         << "   |   " << creal(Uz) << "  " << cimag(Uz)
         << "   |   " << cabs(Uz2-Uz) << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << setprecision(15);

}




/**
 *   \brief Evaluate the the vector field at order m at time n*t: expansion vs real implementation
 **/
void evaluateVectorField(cdouble s0[],
                         cdouble VF1[],
                         cdouble VF2[],
                         vector<Oftsc> &W,    //vector of size OFTS_ORDER
                         vector<Oftsc> &FW,   //vector of size OFTS_ORDER
                         double t,
                         QBCP_L& qbcp_l,
                         int m)
{
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(2);
    //------------------------------------------
    // 1. Init
    //------------------------------------------
    double n  = qbcp_l.us.n;
    cdouble z0[6];
    double z0d[6], vfd[6];
    vector<Ofsc> z0t(6);


    //------------------------------------------
    // 2. Evaluate W
    //------------------------------------------
    for(int p = 0; p < 6; p++)
    {
        W[p].evaluate(s0, z0t[p], m, OFS_ORDER);
        z0[p] = z0t[p].evaluate(n*t);
    }


    //------------------------------------------
    // 3. dot(z) = F(z)
    //------------------------------------------
    //Select real part
    for(int p = 0; p < 6; p++) z0d[p] = creal(z0[p]);
    //Compute vector field
    qbfbp_vfn_novar(t, z0d, vfd, &qbcp_l);
    //Store vector field
    for(int p = 0; p < 6; p++) VF1[p] = vfd[p]+0.0*I;

    //------------------------------------------
    // 2. Evaluate FW in VF2
    //------------------------------------------
    for(int p = 0; p < 6; p++)
    {
        FW[p].evaluate(s0, z0t[p], m, OFS_ORDER);
        VF2[p] = z0t[p].evaluate(n*t);
    }

    cout << "------------------------------------------------------------------------------" << endl;
    cout << " Evaluate the Vector Field. @ order " << m  << endl;
    cout << " n°   |       From true VF       |      From expansions     |   Delta " << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    for(int p = 0; p < 6; p++)
    {
        cout << " " << p
             << "   |   " <<creal(VF1[p]) << "  " << cimag(VF1[p])
             << "   |   " << creal(VF2[p]) << "  " << cimag(VF2[p])
             << "   |   " << cabs(VF1[p]-VF2[p]) << endl;
    }
    cout << "------------------------------------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
}

/**
 *  \brief Earth potential = (1-mu)/(gamma^3*||qpe||)
 **/
cdouble EarthPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xe, QBCP_L& qbcp_l)
{
    double me     = qbcp_l.us.me;
    double n      = qbcp_l.us.n;
    double gamma  = qbcp_l.cs.gamma;

    //Le^2 = xe*xe + ye*ye
    cdouble  xxe = (W[0].evaluate(n*t) - Xe[0].evaluate(n*t));
    cdouble  yye = (W[1].evaluate(n*t) - Xe[1].evaluate(n*t));
    cdouble  zze = (W[2].evaluate(n*t));
    cdouble  qpe = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);


    return me/(pow(gamma, 3.0)*qpe);
}


/**
 *  \brief Moon potential = mu/(gamma^3*||qpm||)
 **/
cdouble MoonPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xm, QBCP_L& qbcp_l)
{
    double mm     = qbcp_l.us.mm;
    double n      = qbcp_l.us.n;
    double gamma  = qbcp_l.cs.gamma;

    //Le^2 = xe*xe + ye*ye
    cdouble  xxe = (W[0].evaluate(n*t) - Xm[0].evaluate(n*t));
    cdouble  yye = (W[1].evaluate(n*t) - Xm[1].evaluate(n*t));
    cdouble  zze = (W[2].evaluate(n*t));
    cdouble  qpm = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);

    return mm/(pow(gamma, 3.0)*qpm);
}


/**
 *  \brief Sun potential = ms/(gamma^3*||qps||)
 **/
cdouble SunPotential(double t, vector<Ofsc> &W, vector<Ofsc> &Xs, QBCP_L& qbcp_l)
{
    double ms     = qbcp_l.us.ms;
    double n      = qbcp_l.us.n;
    double gamma  = qbcp_l.cs.gamma;

    //Le^2 = xe*xe + ye*ye
    cdouble  xxe = (W[0].evaluate(n*t) - Xs[0].evaluate(n*t));
    cdouble  yye = (W[1].evaluate(n*t) - Xs[1].evaluate(n*t));
    cdouble  zze = (W[2].evaluate(n*t));
    cdouble  qps = cpow(xxe*xxe + yye*yye + zze*zze, 1.0/2+0.0*I);


    return ms/(pow(gamma, 3.0)*qps);

}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//         Inner routines
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Prints a given vector W of type \c T in a txt files of the form "filename+i.txt", with i = 0, length(W)-1. The type \c T must have an implementation of the operator \c <<.
 **/
template <typename T> void  vector_fprinf(T& W, string filename)
{
    ifstream readStream;
    ofstream myfile;
    string ss1;
    //Loop on all coefficients
    for(int i = 0; i < (int) W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        myfile.open ((filename+"["+ss1+"].txt").c_str());
        myfile << W[i] << endl;
        myfile.close();
    }
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFS version of the recurrence
//
//---------------------------------------------------------------------------------------------------------------------------------------
//Test the recurrence routines
void testLegendreRecurrence_OFS()
{
    //------------------------------------------
    // I/O
    //------------------------------------------
    ostringstream strs;
    strs << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(3);
    cout << setiosflags(ios::right) << std::showpos << setiosflags(ios::scientific) << setprecision(15);

    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(6,6);
    //The matrix Q of the c.o.c. (Floquet part). Q = inv(P)
    matrix<Ofsc> Q(6,6);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(6);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);
    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //Init routine
    initCOC_OFS(P, Q, V, Xe, Xm, Xs, ILe, ILm, ILs, SEML);


    //------------------------------------------
    // Initialisation of the manifold
    //------------------------------------------
    double xMax = 0.05;
    //vector of size 6
    vector<Ofsc> W(6);
    vector<Ofsc> Wc(6);
    for(int j = 0; j < 6;  j++) W[j] = Ofsc(OFS_ORDER);
    //x =  y = z = xMax
    W[0].setCoef(xMax, 0);
    W[1].setCoef(xMax, 0);
    W[2].setCoef(xMax, 0);

    //Update the position
    Wc[0].setCoef(xMax, 0);
    Wc[1].setCoef(xMax, 0);
    Wc[2].setCoef(xMax, 0);
    //Applying the COC
    //applyModifiedCOC_OFS(P, V, Wc, W);


    //Ofs containing W[0]^2+W[1]^2+W[2]^2
    Ofsc Rho2(OFS_ORDER);

    //------------------------------------------
    // Initialisation of the Legendre objects
    //------------------------------------------
    //Potential
    vector<Ofsc> En(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Sn(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Mn(Csts::POTENTIAL_ORDER+1);
    //Potential derivatives
    vector<Ofsc> Enx(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Eny(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Enz(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Mnx(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Mny(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Mnz(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Snx(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Sny(Csts::POTENTIAL_ORDER+1);
    vector<Ofsc> Snz(Csts::POTENTIAL_ORDER+1);

    //Init within the vectors
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) En[j]  = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Sn[j]  = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Mn[j]  = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Enx[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Eny[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Enz[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Mnx[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Mny[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Mnz[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Snx[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Sny[j] = Ofsc(OFS_ORDER);
    for(int j = 0; j < Csts::POTENTIAL_ORDER;  j++) Snz[j] = Ofsc(OFS_ORDER);

    //Update the potential
    updatePotentialExpansion_OFS(W, Rho2, En, Mn, Sn, Enx, Eny, Mnx, Mny, Mnz, Snx, Sny, Snz, Xe, Xm, Xs, ILe, ILm, ILs, SEML);

    //---------------------------------------------------
    //Test
    //---------------------------------------------------
    //---------------------------------------------------
    double t = 1.0;

    //---------------------------------------------------
    //Earth
    //---------------------------------------------------
    cdouble Epn = EarthPotential(t, W, Xe, SEML);
    cdouble Ep  = 0.0+0.0*I;
    for(int i = 0; i<= Csts::POTENTIAL_ORDER; i++) Ep+= En[i].evaluate(SEML.us.n*t);

    cout << setprecision(3);
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Potential at time t = " << t     << endl;
    cout << "At position (" << creal(W[0].ofs_getCoef(0)) << ", " << creal(W[1].ofs_getCoef(0)) << ", " << creal(W[2].ofs_getCoef(0)) << ")"  << endl;
    cout << "In normalized-centered coordinates" << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << setprecision(15);
    cout << "Primary      Numerical result       Analytical result       Error Ratio" << endl;
    cout << "Earth   " << cabs(Ep-Epn)/cabs(Ep) << endl;

    //---------------------------------------------------
    //Moon
    //---------------------------------------------------
    Epn = MoonPotential(t, W, Xm, SEML);
    Ep = 0.0+0.0*I;
    for(int i = 0; i<= Csts::POTENTIAL_ORDER; i++) Ep+= creal(Mn[i].evaluate(SEML.us.n*t));
    cout << "Moon    " << cabs(Ep-Epn)/cabs(Ep) << endl;

    //---------------------------------------------------
    //Sun
    //---------------------------------------------------
    Epn = SunPotential(t, W, Xs, SEML);
    Ep = 0.0+0.0*I;
    for(int i = 0; i<= Csts::POTENTIAL_ORDER; i++) Ep+= creal(Sn[i].evaluate(SEML.us.n*t));
    cout << "Sun     " << cabs(Ep-Epn)/cabs(Ep) << endl;
    //cout << "Sun potential = " << cabs(Ep) << endl;


    //---------------------------------------------------
    //Plot test
    //---------------------------------------------------
    //---------------------------------------------------
    int N =17;
    double REn[N];
    double RMn[N];
    double RSn[N];
    double xMaxv[N];
    double distToOrigin[N];


    //Plotting devices
    //-------------------------------------------------
    double duration;
    char ch; //Used to close the gnuplot windows at the very end of the program
    gnuplot_ctrl  *h1, *h2, *h3;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    string ss, ssd;
    int iter = 1;
    for(int nPot = 10 ; nPot <= Csts::POTENTIAL_ORDER ; nPot+=10)
    {

        //Grid
        xMaxv[0]  = 1e-4;
        xMaxv[1]  = 1e-3;
        xMaxv[2]  = 2.5e-3;
        xMaxv[3]  = 5e-3;
        xMaxv[4]  = 7.5e-3;
        xMaxv[5]  = 1e-2;
        xMaxv[6]  = 2.5e-2;
        xMaxv[7]  = 0.05;
        xMaxv[8]  = 0.075;
        xMaxv[9]  = 0.1;
        xMaxv[10]  = 0.15;
        xMaxv[11] = 0.2;
        xMaxv[12] = 0.25;
        xMaxv[13] = 0.3;
        xMaxv[14] = 0.35;
        xMaxv[15] = 0.4;
        xMaxv[16] = 0.45;
//        xMaxv[17] = 0.5;
//        xMaxv[18] = 0.6;


        //Loop
        for(int k = 0 ; k < N ; k++)
        {
            //Update the distance from the origin
            xMax = xMaxv[k];

            //Update the position
            Wc[0].setCoef(xMax, 0);
            Wc[1].setCoef(xMax, 0);
            Wc[2].setCoef(xMax, 0);

            //Update the position
            W[0].setCoef(xMax, 0);
            W[1].setCoef(xMax, 0);
            W[2].setCoef(xMax, 0);

            //Applying the COC
            applyModifiedCOC_OFS(P, V, Wc, W);

            //Update the analytical expansion of the potential
            tic();
            updatePotentialExpansion_OFS(W, Rho2, En, Mn, Sn, Enx, Eny, Mnx, Mny, Mnz, Snx, Sny, Snz, Xe, Xm, Xs, ILe, ILm, ILs, SEML);
            duration = toc();

            //---------------------------------------------------
            //Earth
            //---------------------------------------------------
            Epn = EarthPotential(t, W, Xe, SEML);
            Ep = 0.0+0.0*I;
            for(int i = 0; i<= nPot; i++) Ep+= En[i].evaluate(SEML.us.n*t);
            REn[k] = cabs(Ep-Epn)/cabs(Ep);

            //---------------------------------------------------
            //Moon
            //---------------------------------------------------
            Epn = MoonPotential(t, W, Xm, SEML);
            Ep = 0.0+0.0*I;
            for(int i = 0; i<= nPot; i++) Ep+= Mn[i].evaluate(SEML.us.n*t);
            RMn[k] = cabs(Ep-Epn)/cabs(Ep);

            //---------------------------------------------------
            //Sun
            //---------------------------------------------------
            Epn = SunPotential(t, W, Xs, SEML);
            Ep = 0.0+0.0*I;
            for(int i = 0; i<= nPot; i++) Ep+= Sn[i].evaluate(SEML.us.n*t);
            RSn[k] = cabs(Ep-Epn)/cabs(Ep);

            distToOrigin[k] = sqrt( pow(cabs(W[0].evaluate(0.0)), 2.0) + pow(cabs(W[1].evaluate(0.0)), 2.0) + pow(cabs(W[2].evaluate(0.0)), 2.0));

        }

        cout << "--------------" << endl;

        //Set grid
        gnuplot_cmd(h1, "set grid");
        gnuplot_cmd(h2, "set grid");
        gnuplot_cmd(h3, "set grid");

        gnuplot_cmd(h1, "set key left top");
        gnuplot_cmd(h2, "set key left top");
        gnuplot_cmd(h3, "set key right bottom");

        gnuplot_cmd(h1, "set title \"Error on the Earth Potential\" ");
        gnuplot_cmd(h2, "set title \"Error on the Moon Potential\" ");
        gnuplot_cmd(h3, "set title \"Error on the Sun Potential\" ");

        gnuplot_cmd(h1, "set logscale y");
        gnuplot_cmd(h2, "set logscale y");
        gnuplot_cmd(h3, "set logscale y");



        gnuplot_cmd(h1, "set format y \"1e\%%L\"");
        gnuplot_cmd(h2, "set format y \"1e\%%L\"");
        gnuplot_cmd(h3, "set format y \"1e\%%L\"");

        gnuplot_set_xlabel(h1, (char*)"Distance to the origin (t = 0) in Normalized-Centered coordinates");
        gnuplot_set_ylabel(h1, (char*)"Error ratio [-]");

        //gnuplot_set_xlabel(h2, (char*)"X  = Y = Z in Normalized-Centered coordinates");
        gnuplot_set_xlabel(h2, (char*)"Distance to the origin (t = 0) in Normalized-Centered coordinates");
        gnuplot_set_ylabel(h2, (char*)"Error ratio [-]");

        gnuplot_set_xlabel(h3, (char*)"Distance to the origin (t = 0) in Normalized-Centered coordinates");
        gnuplot_set_ylabel(h3, (char*)"Error ratio [-]");


        string ss  = static_cast<ostringstream*>( &(ostringstream() << nPot) )->str();
        string ssd = static_cast<ostringstream*>( &(strs << duration) )->str();
        gnuplot_plot_xy(h1, distToOrigin, REn, N, ("N = "+ss+", T/mon. = "+ssd).c_str(), "lines", "1", "3", iter);
        gnuplot_plot_xy(h2, distToOrigin, RMn, N, ("N = "+ss+", T/mon. = "+ssd).c_str(), "lines", "1", "3", iter);
        gnuplot_plot_xy(h3, distToOrigin, RSn, N, ("N = "+ss+", T/mon. = "+ssd).c_str(), "lines", "1", "3", iter);
        strs.str("");
        iter+=1;


        cout << "Last Coef of the Moon potential" << endl;
        cout << "n = " << nPot << "  Mn[0] = " << creal(Mn[nPot].ofs_getCoef(0)) << " " << cimag(Mn[nPot].ofs_getCoef(0)) << endl;
        cout << "-----------------------------------------------------" << endl;

    }


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //Save in EPS format
    gnuplot_cmd(h1, "set terminal postscript eps solid color enhanced");
    gnuplot_cmd(h1, "set output \"plot/EarthPotential_NC.eps\"");
    gnuplot_cmd(h1, "replot");

    gnuplot_cmd(h2, "set terminal postscript eps solid color enhanced");
    gnuplot_cmd(h2, "set output \"plot/MoonPotential_NC.eps\"");
    gnuplot_cmd(h2, "replot");

    gnuplot_cmd(h3, "set terminal postscript eps solid color enhanced");
    gnuplot_cmd(h3, "set output \"plot/SunPotential_NC.eps\"");
    gnuplot_cmd(h3, "replot");

    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);

}

//Init the recurrence on Legendre-derived objects
void initLegendrePoly_OFS( vector<Ofsc> &W,    //vector of size 6
                           vector<Ofsc> &En,   //vector of size OFTS_ORDER
                           vector<Ofsc> &Mn,   //vector of size OFTS_ORDER
                           vector<Ofsc> &Sn,   //vector of size OFTS_ORDER
                           vector<Ofsc> &Enx,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Eny,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Mnx,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Mny,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Snx,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Sny,  //vector of size OFTS_ORDER
                           vector<Ofsc> &Xe,   //modified position of the Earth
                           vector<Ofsc> &Xm,   //modified position of the Moon
                           vector<Ofsc> &Xs,   //modified position of the Sun
                           Ofsc &ILe,          //modified inverse orbit radius of the Earth
                           Ofsc &ILm,          //modified inverse orbit radius of the Moon
                           Ofsc &ILs,          //modified inverse orbit radius of the Sun
                           QBCP_L& qbcp_l)
{
    //------------------------------------------
    // Initialization
    //------------------------------------------
    double gamma =  qbcp_l.cs.gamma;
    double Ke = qbcp_l.us.me/pow(gamma, 3.0);
    double Km = qbcp_l.us.mm/pow(gamma, 3.0);
    double Ks = qbcp_l.us.ms/pow(gamma, 3.0);

    //------------------------------------------
    //Temporary variables
    //------------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);

    //------------------------------------------
    //For the Earth potential
    //------------------------------------------
    //En[0] = Ke*ILe
    //------------------------------------------
    En[0].ofs_mult(ILe,Ke+0.0*I);

    //En[1] = Ke*ILe^3*(Xe*W[0] + Ye*W[1]);
    //------------------------------------------
    //BUX = ILe*En[0] = Ke*ILe^2
    BUX.ofs_prod(En[0], ILe);
    //CUX = ILe*BUX = Ke*ILe^3
    CUX.ofs_prod(BUX, ILe);
    //AUX = CUX*Xe = Ke*ILe^3*Xe
    AUX.ofs_prod(CUX, Xe[0]);
    //BUX = CUX*Ye = Ke*ILe^3*Ye
    BUX.ofs_prod(CUX, Xe[1]);
    //At every order!
    En[1].ofs_sprod(W[0], AUX);
    En[1].ofs_sprod(W[1], BUX);

//    //Test//
//    double n = qbcp_l.qbcp.n;
//    cdouble p1, p2;
//    cout << "--------------------------" << endl;
//    cout << "Test of precision         " << endl;
//    cout << "--------------------------" << endl;
//    p1 = En[1].evaluate(n*1.0);
//    p2 = cpow(ILe.evaluate(n*1.0), 3.0)*Ke*(Xe[0].evaluate(n*1.0)*W[0].evaluate(n*1.0)+Xe[1].evaluate(n*1.0)*W[1].evaluate(n*1.0));
//    cout << "Earth  " << cabs(p1-p2) << endl;

    //------------------------------------------
    //For the derivatives of the Earth potential
    //------------------------------------------
    //Enx[0] = Eny[0] = Enz[0] = 0
    //Enx[1] = Ke*ILe^3*Xe = AUX
    Enx[1].ccopy(AUX);
    //Eny[1] = Ke*ILe^3*Ye = BUX
    Eny[1].ccopy(BUX);
    //Enz[1] = 0


    //------------------------------------------
    //For the Moon potential
    //------------------------------------------
    //Mn[0] = Km*ILm
    //------------------------------------------
    Mn[0].ofs_mult(ILm,Km+0.0*I);

    //Mn[1] = Km*ILm^3*(Xm*W[0] + Ym*W[1]);
    //------------------------------------------
    //BUX = ILm*Mn[0] = Km*ILm^2
    BUX.ofs_prod(Mn[0], ILm);
    //CUX = ILm*BUX = Km*ILm^3
    CUX.ofs_prod(BUX, ILm);
    //AUX = CUX*Xm = Km*ILm^3*Xm
    AUX.ofs_prod(CUX, Xm[0]);
    //BUX = CUX*Ym = Km*ILm^3*Ym
    BUX.ofs_prod(CUX, Xm[1]);
    //At every order!
    Mn[1].ofs_sprod(W[0], AUX);
    Mn[1].ofs_sprod(W[1], BUX);

//    //Test//
//    p1 = Mn[1].evaluate(n*1.0);
//    p2 = cpow(ILm.evaluate(n*1.0), 3.0)*Km*(Xm[0].evaluate(n*1.0)*W[0].evaluate(n*1.0)+Xm[1].evaluate(n*1.0)*W[1].evaluate(n*1.0));
//    cout << "Moon  " << cabs(p1-p2) << endl;


    //------------------------------------------
    //For the derivatives of the Moon potential
    //------------------------------------------
    //Mnx[0] = Mny[0] = Mnz[0] = 0
    //Mnx[1] = Km*ILm^3*Xm = AUX
    Mnx[1].ccopy(AUX);
    //Mny[1] = Km*ILm^3*Ym = BUX
    Mny[1].ccopy(BUX);
    //Mnz[1] = 0


    //------------------------------------------
    //For the Sun potential
    //------------------------------------------
    //Sn[0] = Ks*ILs
    //------------------------------------------
    Sn[0].ofs_mult(ILs,Ks+0.0*I);

    //Sn[1] = Ks*ILs^3*(Xs*W[0] + Ys*W[1]);
    //------------------------------------------
    //BUX = ILs*Sn[0] = Ks*ILs^2
    BUX.ofs_prod(Sn[0], ILs);
    //CUX = ILs*BUX = Ks*ILs^3
    CUX.ofs_prod(BUX, ILs);
    //AUX = CUX*Xs = Ks*ILs^3*Xs
    AUX.ofs_prod(CUX, Xs[0]);
    //BUX = CUX*Ys = Ks*ILs^3*Ys
    BUX.ofs_prod(CUX, Xs[1]);
    //At every order!
    Sn[1].ofs_sprod(W[0], AUX);
    Sn[1].ofs_sprod(W[1], BUX);


//    //Test//
//    p1 = Sn[1].evaluate(n*1.0);
//    p2 = cpow(ILs.evaluate(n*1.0), 3.0)*Ks*(Xs[0].evaluate(n*1.0)*W[0].evaluate(n*1.0)+Xs[1].evaluate(n*1.0)*W[1].evaluate(n*1.0));
//    cout << "Sun  " << cabs(p1-p2) << endl;


    //------------------------------------------
    //For the derivatives of the Sun potential
    //------------------------------------------
    //Snx[0] = Sny[0] = Snz[0] = 0
    //Snx[1] = Ks*ILs^3*Xs = AUX
    Snx[1].ccopy(AUX);
    //Sny[1] = Ks*ILs^3*Ys = BUX
    Sny[1].ccopy(BUX);
    //Snz[1] = 0
}





//Apply recurrence on Legendre-derived objects for m > 1
//Update the order m of the kth coefficient in En, Mn, etc (k < m)
//Provided that the (k-2t)h and (k-1)th coefficient have been updated
void updateLegPoly_OFS(vector<Ofsc> &W,    //vector of size 6
                       Ofsc &Rho2,          //Ofts containing W[0]^2+W[1]^2+W[2]^2
                       vector<Ofsc> &En,   //vector of size OFTS_ORDER
                       vector<Ofsc> &Mn,   //vector of size OFTS_ORDER
                       vector<Ofsc> &Sn,   //vector of size OFTS_ORDER
                       vector<Ofsc> &Xe,    //modified position of the Earth
                       vector<Ofsc> &Xm,    //modified position of the Moon
                       vector<Ofsc> &Xs,    //modified position of the Sun
                       Ofsc &ILe,           //modified inverse orbit radius of the Earth
                       Ofsc &ILm,           //modified inverse orbit radius of the Moon
                       Ofsc &ILs,           //modified inverse orbit radius of the Sun
                       int k)
{
    //At this step, En[0, .., k-1] are at order m
    //We want to update En[k] at order k
    //------------------------------------------
    //Temporary variables
    //------------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);

    //------------------------------------------
    //For the Earth potential
    //------------------------------------------
    //BUX = x*xe + y*ye
    BUX.ofs_sprod(W[0], Xe[0]);
    BUX.ofs_sprod(W[1], Xe[1]);
    //AUX = ILe^2*(2*k-1)/k
    AUX.ofs_mprod(ILe, ILe, (2.0*k - 1.0)/(double)k+0.0*I);
    //CUX = ILe^2*(2*k-1)/k*(x*xe + y*ye)
    CUX.ofs_prod(AUX, BUX);
    //EN[k] += ILe^2*(2*k-1)/k*(x*xe + y*ye)*En[k-1]
    En[k].ofs_sprod(CUX, En[k-1]);
    //AUX = ILe^2*(2*k-1)/k
    AUX.ofs_mprod(ILe, ILe, -(k-1.0)/(double)k+0.0*I);
    //CUX = AUX*Rho2
    CUX.ofs_prod(AUX, Rho2);
    //EN[k] += ILe^2*(2*k-1)/k*Rho2*En[k-2]
    En[k].ofs_sprod(CUX, En[k-2]);


//    double n = qbcp_l.qbcp.n;
//    cdouble p1, p2;
//    //Test//
//    if(k == 2)
//    {
//        p1  = En[k].evaluate(n*1.0);
//        p2  =  (2.0*k - 1.0)/(double)k*ILe.evaluate(n*1.0)*ILe.evaluate(n*1.0)*(W[0].evaluate(n*1.0)*Xe[0].evaluate(n*1.0)
//                                                                                + W[1].evaluate(n*1.0)*Xe[1].evaluate(n*1.0))*En[k-1].evaluate(n*1.0);
//        p2 += -(k-1.0)/(double)k*ILe.evaluate(n*1.0)*ILe.evaluate(n*1.0)*Rho2.evaluate(n*1.0)*En[k-2].evaluate(n*1.0);
//
//        cout << "Earth " << k << " " << cabs(p1-p2) << endl;
//    }


    //------------------------------------------
    //For the Moon potential
    //------------------------------------------
    //BUX = x*xe + y*ye
    BUX.ofs_prod(W[0], Xm[0]);
    BUX.ofs_sprod(W[1], Xm[1]);
    //AUX = ILm^2*(2*k-1)/k
    AUX.ofs_mprod(ILm, ILm, (2.0*k - 1.0)/(double)k+0.0*I);
    //CUX = ILm^2*(2*k-1)/k*(x*xs + y*ys)
    CUX.ofs_prod(AUX, BUX);
    //EN[k] += ILm^2*(2*k-1)/k*(x*xs + y*ys)*Mn[k-1]
    Mn[k].ofs_sprod(CUX, Mn[k-1]);
    //AUX = ILm^2*(2*k-1)/k
    AUX.ofs_mprod(ILm, ILm, -(k-1.0)/(double)k+0.0*I);
    //CUX = AUX*Rho2
    CUX.ofs_prod(AUX, Rho2);
    //Mn[k] += ILm^2*(2*k-1)/k*Rho2*Mn[k-2]
    Mn[k].ofs_sprod(CUX, Mn[k-2]);


//    //Test//
//    if(k == 2)
//    {
//        p1  = Mn[k].evaluate(n*1.0);
//        p2  =  (2.0*k - 1.0)/(double)k*ILm.evaluate(n*1.0)*ILm.evaluate(n*1.0)*(W[0].evaluate(n*1.0)*Xm[0].evaluate(n*1.0)
//                                                                                + W[1].evaluate(n*1.0)*Xm[1].evaluate(n*1.0))*Mn[k-1].evaluate(n*1.0);
//        p2 += -(k-1.0)/(double)k*ILm.evaluate(n*1.0)*ILm.evaluate(n*1.0)*Rho2.evaluate(n*1.0)*Mn[k-2].evaluate(n*1.0);
//
//        cout << "Moon " << k << " " << cabs(p1-p2) << endl;
//    }



    //------------------------------------------
    //For the Sun potential
    //------------------------------------------
    //BUX = x*xe + y*ye
    BUX.ofs_prod(W[0], Xs[0]);
    BUX.ofs_sprod(W[1], Xs[1]);
    //AUX = ILs^2*(2*k-1)/k
    AUX.ofs_mprod(ILs, ILs, (2.0*k - 1.0)/(double)k+0.0*I);
    //CUX = ILs^2*(2*k-1)/k*(x*xs + y*ys)
    CUX.ofs_prod(AUX, BUX);
    //EN[k] += ILs^2*(2*k-1)/k*(x*xs + y*ys)*Sn[k-1]
    Sn[k].ofs_sprod(CUX, Sn[k-1]);
    //AUX = ILs^2*(2*k-1)/k
    AUX.ofs_mprod(ILs, ILs, -(k-1.0)/(double)k+0.0*I);
    //CUX = AUX*Rho2
    CUX.ofs_prod(AUX, Rho2);
    //Sn[k] += ILs^2*(2*k-1)/k*Rho2*Sn[k-2]
    Sn[k].ofs_sprod(CUX, Sn[k-2]);
}

//Update the potential at order m > 1 in s
void updatePotentialExpansion_OFS(vector<Ofsc> &W,    //vector of size 6
                                  Ofsc &Rho2,          //Ofts containing W[0]^2+W[1]^2+W[2]^2
                                  vector<Ofsc> &En,   //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Mn,   //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Sn,   //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Enx,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Eny,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Mnx,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Mny,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Mnz,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Snx,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Sny,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Snz,  //vector of size OFTS_ORDER+1
                                  vector<Ofsc> &Xe,    //modified position of the Earth
                                  vector<Ofsc> &Xm,    //modified position of the Moon
                                  vector<Ofsc> &Xs,    //modified position of the Sun
                                  Ofsc &ILe,           //modified inverse orbit radius of the Earth
                                  Ofsc &ILm,           //modified inverse orbit radius of the Moon
                                  Ofsc &ILs,           //modified inverse orbit radius of the Sun
                                  QBCP_L& qbcp_l)      //QBFPB on a given Li
{

    //------------------------------------------
    //Temporary variables
    //------------------------------------------
    //Ofsc AUX(OFS_ORDER);
    //Ofsc BUX(OFS_ORDER);
    //Ofsc CUX(OFS_ORDER);

    //------------------------------------------
    //Clear
    //------------------------------------------
    Rho2.zero();
    for(int i = 0; i<= Csts::POTENTIAL_ORDER; i++)
    {
        En[i].zero();
        Mn[i].zero();
        Sn[i].zero();

        Mnx[i].zero();
        Mny[i].zero();
        Mnz[i].zero();

        Snx[i].zero();
        Sny[i].zero();
        Snz[i].zero();
    }


    //------------------------------------------
    //Update Rho2
    //------------------------------------------
    Rho2.ofs_sprod(W[0], W[0]);
    Rho2.ofs_sprod(W[1], W[1]);
    Rho2.ofs_sprod(W[2], W[2]);

    //------------------------------------------
    //Initialisation of the recurrence
    //------------------------------------------
    initLegendrePoly_OFS(W, En, Mn, Sn, Enx, Eny, Mnx, Mny, Snx, Sny, Xe, Xm, Xs, ILe, ILm, ILs, qbcp_l);

    //-----------------------------------------------------------------
    //Aplly the recurrence OFTS_ORDER times
    //-----------------------------------------------------------------
    for(int k = 2; k <= Csts::POTENTIAL_ORDER; k++)
    {
        updateLegPoly_OFS(W, Rho2, En, Mn, Sn, Xe, Xm, Xs, ILe, ILm, ILs, k);
    }

}
