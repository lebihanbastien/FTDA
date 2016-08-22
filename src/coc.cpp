#include "coc.h"

/**
 * \file coc.cpp
 * \brief Implements the complete change of coordinates for the Normalized-Centered Hamiltonian of the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  Change of coordinates for the hamiltonian equations of the QBCP
 *  This coc allows to get:
 *  - The dynamical equivalent of Li as a fixed point (get rid of order one of the hamiltonian)
 *  - "normal form" of the order 2 of the hamiltonian
 *
 */


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFTS version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief To TFS format for the whole COC
 **/
void tfs_from_ofs(matrix<Ofsc> &P,
                  matrix<Ofsc> &Q,
                  matrix<Ofsc> &PC,
                  matrix<Ofsc> &PCdot,
                  matrix<Ofsc> &CQ,
                  vector<Ofsc> &Xe,
                  vector<Ofsc> &Xm,
                  vector<Ofsc> &Xs,
                  vector<Ofsc> &V,
                  vector<Ofsc> &Vdot,
                  Ofsc &ILe,
                  Ofsc &ILm,
                  Ofsc &ILs)
{
    tfs_from_ofs_inline(P);
    tfs_from_ofs_inline(Q);
    tfs_from_ofs_inline(PC);
    tfs_from_ofs_inline(PCdot);
    tfs_from_ofs_inline(CQ);
    tfs_from_ofs_inline(V);
    tfs_from_ofs_inline(Xe);
    tfs_from_ofs_inline(Xm);
    tfs_from_ofs_inline(Xs);
    tfs_from_ofs_inline(Vdot);

    Ofsc temp(ILe.getOrder());
    ILe.tfs_from_ofs_inline(temp);
    ILm.tfs_from_ofs_inline(temp);
    ILs.tfs_from_ofs_inline(temp);
}


/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *         This the TFS version of this routine: it retrieves the Fourier series in time domain.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c qbcp_l.F_COC.
 *      * Vdot = dot(V) is computed
 *      * The products P*C and C*Q are computed in PC and QC
 *      * dot(PC) = PCdot is computed
 *
 *  2. Second step
 *  The vectors Xe[3:5], Xm[3:5], and X[3:5] contains the time-dependent positions of the primaries in the xy-plane.
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[3]^2 + Xc[4]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 **/
void tfts_initCOC(matrix<Ofsc> &P,
                  matrix<Ofsc> &Q,
                  matrix<Ofsc> &PC,
                  matrix<Ofsc> &PCdot,
                  matrix<Ofsc> &CQ,
                  vector<Ofsc> &Xe,
                  vector<Ofsc> &Xm,
                  vector<Ofsc> &Xs,
                  vector<Ofsc> &V,
                  vector<Ofsc> &Vdot,
                  Ofsc &ILe,
                  Ofsc &ILm,
                  Ofsc &ILs,
                  QBCP_L& qbcp_l)
{
    //----------------------------
    //Retrieve folder
    //----------------------------
    string F_COC = qbcp_l.cs.F_COC;

    //----------------------------
    //Switch RTBP/QBCP/BCP
    //----------------------------
    if(qbcp_l.model == M_QBCP || qbcp_l.model == M_BCP)
    {
        //Recovering the data: matrix P
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*P.getCA(i,j), F_COC+"P",  i+1, j+1);
        //Recovering the data: matrix Q
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*Q.getCA(i,j), F_COC+"Q",  i+1, j+1);

        //Recovering the data: vector V
        //V = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:V[2] and V[5] are kept null
        readCOC(V[0], F_COC+"G1_",  1, 1);
        readCOC(V[1], F_COC+"G1_",  1, 2);
        readCOC(V[3], F_COC+"G1_",  2, 1);
        readCOC(V[4], F_COC+"G1_",  2, 2);


        //Recovering the data: vector Xm
        readCOC(Xm[0], F_COC+"Xm",  4, 1);
        readCOC(Xm[1], F_COC+"Xm",  5, 1);
        //Note: Xm[2] is kept null

        //Recovering the data: vector Xe
        readCOC(Xe[0], F_COC+"Xe",  4, 1);
        readCOC(Xe[1], F_COC+"Xe",  5, 1);
        //Note: Xe[2] is kept null

        //Recovering the data: vector Xs
        readCOC(Xs[0], F_COC+"Xs",  4, 1);
        readCOC(Xs[1], F_COC+"Xs",  5, 1);
        //Note: Xs[2] is kept null


        //Recovering the data: Ofts ILe, ILm, ILs
        readCOC(ILe, F_COC+"Xe",  6, 1);
        readCOC(ILm, F_COC+"Xm",  6, 1);
        readCOC(ILs, F_COC+"Xs",  6, 1);

    }
    else //RTPB
    {
        //note that V is left untouched (set to zero by default)
        //--------------------

        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, dl1, do1, s1, s2, c2;
        c2 = qbcp_l.cs.c2;
        eta1 = (c2 - 2.0 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2.0 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        dl1 = 2*la1*( (4.0+3*c2)*la1*la1 + 4 + 5*c2 - 6*c2*c2);
        do1 =   om1*( (4.0+3*c2)*om1*om1 - 4 - 5*c2 + 6*c2*c2);
        s1  = sqrt(dl1);
        s2  = sqrt(do1);

        //--------------------
        //Init P
        //--------------------
        P.setCoef(+2*la1/s1,                          0, 1);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 1);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        P.setCoef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        P.setCoef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        P.setCoef(-(om1*om1 - 2*c2 - 1)/s2,           3, 0);

        P.setCoef(-2*la1/s1,                          0, 4);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 4);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        P.setCoef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        P.setCoef(+2*om1/s2,                          0, 3);
        P.setCoef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        P.setCoef(+1.0/sqrt(om2),                     2, 2);
        P.setCoef(+sqrt(om2),                         5, 5);

        //--------------------
        //Check symplectic nature of P
        //--------------------
        /*
        gsl_matrix *Ps = gsl_matrix_calloc(6,6);
        for(int i = 0; i <6; i++)
        {
            for(int j = 0; j <6; j++)
            {
                gsl_matrix_set(Ps, i, j, creal(P.getCoef(i,j).ofs_getCoef(0)));
            }
        }
        realSymplecticMatrixTest(Ps, INVERSE_GSL);
        gsl_matrix_free(Ps);
        //*/

        //--------------------
        //Q = inv(P) (gsl lib)
        //--------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_calloc (NV, NV);
        gsl_matrix* Qc   = gsl_matrix_calloc (NV, NV);
        gsl_permutation * p6 = gsl_permutation_alloc (NV);

        //Init Pc
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) gsl_matrix_set(Pc, i, j, creal(P.getCA(i,j)->ofs_getCoef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);


        //--------------------
        // Init Q
        //--------------------
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) Q.setCoef(gsl_matrix_get(Qc, i, j), i, j);


        double gamma = qbcp_l.cs.gamma;
        switch(qbcp_l.fwrk)
        {
        case F_EM:
            //Sun does not exist
            Xs[0].setCoef(  0.0+0.0*I, 0);
            Xs[1].setCoef(  0.0+0.0*I, 0);
            ILs.setCoef(    0.0+0.0*I, 0);

            //Earth and Moon
            if(qbcp_l.li_EM == 1)
            {
                Xe[0].setCoef((-1.0+gamma)/gamma+0.0*I, 0);
                Xe[1].setCoef( 0.0+0.0*I, 0);
                ILe.setCoef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xm[0].setCoef( +1.0+0.0*I, 0);
                Xm[1].setCoef(  0.0+0.0*I, 0);
                ILm.setCoef(  1.0+0.0*I, 0);
            }
            else
            {
                Xe[0].setCoef((-1.0-gamma)/gamma+0.0*I, 0);
                Xe[1].setCoef( 0.0+0.0*I, 0);
                ILe.setCoef( gamma/(1.0+gamma)+0.0*I, 0);

                Xm[0].setCoef( -1.0+0.0*I, 0);
                Xm[1].setCoef(  0.0+0.0*I, 0);
                ILm.setCoef(    1.0+0.0*I, 0);
            }
            break;

        case F_SEM:
            //Moon does not exist
            Xm[0].setCoef(  0.0+0.0*I, 0);
            Xm[1].setCoef(  0.0+0.0*I, 0);
            ILm.setCoef(    0.0+0.0*I, 0);

            //Sun and Moon
            if(qbcp_l.li_SEM == 1)
            {
                Xs[0].setCoef((-1.0+gamma)/gamma+0.0*I, 0);
                Xs[1].setCoef( 0.0+0.0*I, 0);
                ILs.setCoef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xe[0].setCoef( +1.0+0.0*I, 0);
                Xe[1].setCoef(  0.0+0.0*I, 0);
                ILe.setCoef(  1.0+0.0*I, 0);
            }
            else
            {
                Xs[0].setCoef((-1.0-gamma)/gamma+0.0*I, 0);
                Xs[1].setCoef( 0.0+0.0*I, 0);
                ILs.setCoef( gamma/(1.0+gamma)+0.0*I, 0);

                Xe[0].setCoef( -1.0+0.0*I, 0);
                Xe[1].setCoef(  0.0+0.0*I, 0);
                ILe.setCoef(    1.0+0.0*I, 0);
            }
            break;
        }
    }

    //----------------------------
    // Building PC, CQ
    //----------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize PC and CQ
    vector<int> keyMap(6);
    keyMap[0] = 0;
    keyMap[1] = 1;
    keyMap[2] = 3;
    keyMap[3] = 4;
    keyMap[4] = 2;
    keyMap[5] = 5;


    int ii;
    //Init PC by rows
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,0),   1.0/sqrt(2)+0.0*I, P(ii,3), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 0);
        BUX.ofs_fsum(P(ii,0), I*1.0/sqrt(2), P(ii,3),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 3);
        PC.setCoef(P(ii,1), ii, 1);
        PC.setCoef(P(ii,4), ii, 4);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,2),   1.0/sqrt(2)+0.0*I, P(ii,5), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 2);
        BUX.ofs_fsum(P(ii,2), I*1.0/sqrt(2), P(ii,5),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 5);
    }

    //Init CQ by columns
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(0,ii),    1.0/sqrt(2)+0.0*I, Q(3,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 0, ii);
        BUX.ofs_fsum(Q(0,ii), -1.0/sqrt(2)*I, Q(3,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 3, ii);
        CQ.setCoef(Q(1,ii), 1, ii);
        CQ.setCoef(Q(4,ii), 4, ii);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(2,ii),    1.0/sqrt(2)+0.0*I, Q(5,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 2, ii);
        BUX.ofs_fsum(Q(2,ii), -1.0/sqrt(2)*I, Q(5,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 5, ii);
    }

    //----------------------------
    // Building Vdot
    //----------------------------
    for(int i = 0; i < NV; i++)
    {
        Vdot[i].dot(V[i], qbcp_l.us.n);
    }
    //Init PCdot
    PCdot.dot(PC, qbcp_l.us.n);
}


/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *         This the TFS version of this routine: it retrieves the Fourier series in frequency domain.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c qbcp_l.F_COC.
 *      * Vdot = dot(V) is computed
 *      * The products P*C and C*Q are computed in PC and QC
 *      * dot(PC) = PCdot is computed
 *
 *  2. Second step
 *  The vectors Xe[2], Xm[2], and X[s] contains the time-dependent positions of the primaries in the xy-plane.
 *  Note that this is the TRANSLATED versions of these positions i.e. they are taking into account the motion
 *  of the periodic orbit, dynamical equivalent of the libration point L1,2.
 *
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[0]^2 + Xc[1]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 **/
void initCOC(matrix<Ofsc> &P,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &PCdot,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &Xe,
             vector<Ofsc> &Xm,
             vector<Ofsc> &Xs,
             vector<Ofsc> &V,
             vector<Ofsc> &Vdot,
             Ofsc &ILe,
             Ofsc &ILm,
             Ofsc &ILs,
             QBCP_L& qbcp_l)
{
    //----------------------------
    //Retrieve folder
    //----------------------------
    string F_COC = qbcp_l.cs.F_COC;

    //----------------------------
    //Switch RTBP/QBCP/BCP
    //----------------------------
    if(qbcp_l.model == M_QBCP || qbcp_l.model == M_BCP)
    {
        //Recovering the data: matrix P
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*P.getCA(i,j), F_COC+"P",  i+1, j+1);
        //Recovering the data: matrix Q
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*Q.getCA(i,j), F_COC+"Q",  i+1, j+1);

        //Recovering the data: vector V
        //V = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:V[2] and V[5] are kept null
        readCOC(V[0], F_COC+"G1_",  1, 1);
        readCOC(V[1], F_COC+"G1_",  1, 2);
        readCOC(V[3], F_COC+"G1_",  2, 1);
        readCOC(V[4], F_COC+"G1_",  2, 2);


        //Recovering the data: vector Xm
        readCOC(Xm[0], F_COC+"Xm",  1, 1);
        readCOC(Xm[1], F_COC+"Xm",  2, 1);
        //Note: Xm[2] is kept null

        //Recovering the data: vector Xe
        readCOC(Xe[0], F_COC+"Xe",  1, 1);
        readCOC(Xe[1], F_COC+"Xe",  2, 1);
        //Note: Xe[2] is kept null

        //Recovering the data: vector Xs
        readCOC(Xs[0], F_COC+"Xs",  1, 1);
        readCOC(Xs[1], F_COC+"Xs",  2, 1);
        //Note: Xs[2] is kept null


        //Recovering the data: Ofts ILe, ILm, ILs
        readCOC(ILe, F_COC+"Xe",  3, 1);
        readCOC(ILm, F_COC+"Xm",  3, 1);
        readCOC(ILs, F_COC+"Xs",  3, 1);

    }
    else //RTPB
    {
        //note that V is left untouched (set to zero by default)
        //--------------------

        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, dl1, do1, s1, s2, c2;
        c2 = qbcp_l.cs.c2;
        eta1 = (c2 - 2.0 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2.0 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        dl1 = 2*la1*( (4.0+3*c2)*la1*la1 + 4 + 5*c2 - 6*c2*c2);
        do1 =   om1*( (4.0+3*c2)*om1*om1 - 4 - 5*c2 + 6*c2*c2);
        s1  = sqrt(dl1);
        s2  = sqrt(do1);

        //--------------------
        //Init P
        //--------------------
        P.setCoef(+2*la1/s1,                          0, 1);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 1);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        P.setCoef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        P.setCoef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        P.setCoef(-(om1*om1 - 2*c2 - 1)/s2,           3, 0);

        P.setCoef(-2*la1/s1,                          0, 4);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 4);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        P.setCoef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        P.setCoef(+2*om1/s2,                          0, 3);
        P.setCoef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        P.setCoef(+1.0/sqrt(om2),                     2, 2);
        P.setCoef(+sqrt(om2),                         5, 5);

        //--------------------
        //Check symplectic nature of P
        //--------------------
        /*
        gsl_matrix *Ps = gsl_matrix_calloc(6,6);
        for(int i = 0; i <6; i++)
        {
            for(int j = 0; j <6; j++)
            {
                gsl_matrix_set(Ps, i, j, creal(P.getCoef(i,j).ofs_getCoef(0)));
            }
        }
        realSymplecticMatrixTest(Ps, INVERSE_GSL);
        gsl_matrix_free(Ps);
        //*/

        //--------------------
        //Q = inv(P) (gsl lib)
        //--------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_calloc (NV, NV);
        gsl_matrix* Qc   = gsl_matrix_calloc (NV, NV);
        gsl_permutation * p6 = gsl_permutation_alloc (NV);

        //Init Pc
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) gsl_matrix_set(Pc, i, j, creal(P.getCA(i,j)->ofs_getCoef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);


        //--------------------
        // Init Q
        //--------------------
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) Q.setCoef(gsl_matrix_get(Qc, i, j), i, j);


        double gamma = qbcp_l.cs.gamma;
        switch(qbcp_l.fwrk)
        {
        case F_EM:
            //Sun does not exist
            Xs[0].setCoef(  0.0+0.0*I, 0);
            Xs[1].setCoef(  0.0+0.0*I, 0);
            ILs.setCoef(    0.0+0.0*I, 0);

            //Earth and Moon
            if(qbcp_l.li_EM == 1)
            {
                Xe[0].setCoef((-1.0+gamma)/gamma+0.0*I, 0);
                Xe[1].setCoef( 0.0+0.0*I, 0);
                ILe.setCoef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xm[0].setCoef( +1.0+0.0*I, 0);
                Xm[1].setCoef(  0.0+0.0*I, 0);
                ILm.setCoef(  1.0+0.0*I, 0);
            }
            else
            {
                Xe[0].setCoef((-1.0-gamma)/gamma+0.0*I, 0);
                Xe[1].setCoef( 0.0+0.0*I, 0);
                ILe.setCoef( gamma/(1.0+gamma)+0.0*I, 0);

                Xm[0].setCoef( -1.0+0.0*I, 0);
                Xm[1].setCoef(  0.0+0.0*I, 0);
                ILm.setCoef(    1.0+0.0*I, 0);
            }
            break;

        case F_SEM:
            //Moon does not exist
            Xm[0].setCoef(  0.0+0.0*I, 0);
            Xm[1].setCoef(  0.0+0.0*I, 0);
            ILm.setCoef(    0.0+0.0*I, 0);

            //Sun and Moon
            if(qbcp_l.li_SEM == 1)
            {
                Xs[0].setCoef((-1.0+gamma)/gamma+0.0*I, 0);
                Xs[1].setCoef( 0.0+0.0*I, 0);
                ILs.setCoef( 1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

                Xe[0].setCoef( +1.0+0.0*I, 0);
                Xe[1].setCoef(  0.0+0.0*I, 0);
                ILe.setCoef(  1.0+0.0*I, 0);
            }
            else
            {
                Xs[0].setCoef((-1.0-gamma)/gamma+0.0*I, 0);
                Xs[1].setCoef( 0.0+0.0*I, 0);
                ILs.setCoef( gamma/(1.0+gamma)+0.0*I, 0);

                Xe[0].setCoef( -1.0+0.0*I, 0);
                Xe[1].setCoef(  0.0+0.0*I, 0);
                ILe.setCoef(    1.0+0.0*I, 0);
            }
            break;
        }
    }

    //----------------------------
    // Building PC, CQ
    //----------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize PC and CQ
    vector<int> keyMap(6);
    keyMap[0] = 0;
    keyMap[1] = 1;
    keyMap[2] = 3;
    keyMap[3] = 4;
    keyMap[4] = 2;
    keyMap[5] = 5;


    int ii;
    //Init PC by rows
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,0),   1.0/sqrt(2)+0.0*I, P(ii,3), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 0);
        BUX.ofs_fsum(P(ii,0), I*1.0/sqrt(2), P(ii,3),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 3);
        PC.setCoef(P(ii,1), ii, 1);
        PC.setCoef(P(ii,4), ii, 4);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,2),   1.0/sqrt(2)+0.0*I, P(ii,5), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 2);
        BUX.ofs_fsum(P(ii,2), I*1.0/sqrt(2), P(ii,5),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 5);
    }

    //Init CQ by columns
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(0,ii),    1.0/sqrt(2)+0.0*I, Q(3,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 0, ii);
        BUX.ofs_fsum(Q(0,ii), -1.0/sqrt(2)*I, Q(3,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 3, ii);
        CQ.setCoef(Q(1,ii), 1, ii);
        CQ.setCoef(Q(4,ii), 4, ii);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(2,ii),    1.0/sqrt(2)+0.0*I, Q(5,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 2, ii);
        BUX.ofs_fsum(Q(2,ii), -1.0/sqrt(2)*I, Q(5,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 5, ii);
    }

    //----------------------------
    // Building Vdot
    //----------------------------
    for(int i = 0; i < NV; i++)
    {
        Vdot[i].dot(V[i], qbcp_l.us.n);
    }
    //Init PCdot
    PCdot.dot(PC, qbcp_l.us.n);
}


/**
 *  \brief Apply the change of variables at every order in zIN/zOut.
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Oftsc> &zIn,
              vector<Oftsc> &zOut)
{
    //zeroing the target
    for(unsigned int i = 0; i < zOut.size(); i++) zOut[i].zero();
    //zOut = PC*zIn
    smvprod_u(PC, zIn, zOut);
    //zOut+=V(theta)
    addCoef(V, zOut);
}


/**
 *  \brief Apply the change of variables at order m in zIN/zOut, with or without the zero order V
 *         The change of variables is of the form: zOut = PC * zIN (+ V)
 *         If flag == 1, the zero order is added
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Oftsc> &zIn,
              vector<Oftsc> &zOut,
              int m,
              int flag)
{
    //zOut = PC*zIn
    smvprod_u(PC, zIn, zOut, m);
    //zOut+=V(theta)
    if(m == 0 && flag == 1) addCoef(V, zOut);
}

/**
 *  \brief Apply the change of variables at order m in zIN/zOut, with or without the zero order V
 *         The change of variables is of the form: zOut = PC * zIN (+ V)
 *         If flag == 1, the zero order is added.
 *         This is the TFS version of this routine (time domain).
 **/
void tfts_applyCOC(matrix<Ofsc> &PC,
                   vector<Ofsc> &V,
                   vector<Oftsc> &zIn,
                   vector<Oftsc> &zOut,
                   int m,
                   int flag)
{
    //zOut = PC*zIn
    tfts_smvprod_u(PC, zIn, zOut, m);
    //zOut+=V(theta)
    if(m == 0 && flag == 1) tfts_addCoef(V, zOut);
}


/**
 *  \brief Apply the inverse change of variables at every orders in zIN/zOut.
 *         The change of variables is of the form: zOut = CQ*(zIn - V).
 *
 *   Note: stp is a temporary variable required to store zIn-V
 **/
void applyInvCOC(matrix<Ofsc> &CQ,
                 vector<Ofsc> &V,
                 vector<Oftsc> &zIn,
                 vector<Oftsc> &zOut,
                 vector<Oftsc> &ztp)
{
    //-----------------------------------------
    //ztp = zIn - V
    //-----------------------------------------
    subCoef(V, zIn, ztp);
    //-----------------------------------------
    //zOut = CQ*ztp
    //-----------------------------------------
    smvprod_u(CQ, ztp, zOut);
}


/**
 *  \brief Apply the derivative of the change of variables at every orders in zIN/zOut.
 *         The change of variables is of the form: zdot = PC zh + PC zhdot + Vdot
 **/
void applyDotCOC(matrix<Ofsc> &PC,
                 matrix<Ofsc> &PCdot,
                 vector<Ofsc> &Vdot,
                 vector<Oftsc> &zh,
                 vector<Oftsc> &zhdot,
                 vector<Oftsc> &zdot)
{
    //Zeroing the target
    for(unsigned int i = 0; i < zdot.size(); i++) zdot[i].zero();
    //zdot += PCdot*zh
    smvprod_u(PCdot, zh, zdot);
    //zdot += PC*zhdot
    smvprod_u(PC, zhdot, zdot);
    //zdot +=Vdot(theta)
    addCoef(Vdot, zdot);
}


/**
 *  \brief Apply the derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zdot = PC zh + PC zhdot + Vdot
 **/
void applyDotCOC(matrix<Ofsc> &PC,
                 matrix<Ofsc> &PCdot,
                 vector<Ofsc> &Vdot,
                 vector<Oftsc> &zh,
                 vector<Oftsc> &zhdot,
                 vector<Oftsc> &zdot,
                 int m)
{
    //zdot += PCdot*zh
    smvprod_u(PCdot, zh, zdot, m);
    //zdot += PC*zhdot
    smvprod_u(PC, zhdot, zdot, m);
    //zdot +=Vdot(theta)
    if(m == 0) addCoef(Vdot, zdot);
}


/**
 *  \brief Apply the inverse derivative of the change of variables at every orders m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 *
 *         Warning: zhdot is used as temporary variable within the routine
 *         However, the final result is good.
 **/
void applyInvDotCOC(matrix<Ofsc> &CQ,
                    matrix<Ofsc> &PCdot,
                    vector<Ofsc> &Vdot,
                    vector<Oftsc> &zh,
                    vector<Oftsc> &zdot,
                    vector<Oftsc> &zhdot,
                    vector<Oftsc> &ztp_1,
                    vector<Oftsc> &ztp_2)
{
    //Zeroing the target
    for(unsigned int i = 0; i < zhdot.size(); i++) zhdot[i].zero();
    for(unsigned int i = 0; i < zhdot.size(); i++) ztp_1[i].zero();
    for(unsigned int i = 0; i < zhdot.size(); i++) ztp_2[i].zero();

    //ztp_1 += PCdot*zh @order m
    smvprod_u(PCdot, zh, ztp_1);
    //ztp_2 += zdot - PCdot*zh - Vdot @order m
    for(unsigned int i = 0; i < zh.size(); i++)
    {
        ztp_2[i].ofts_smult_u(zdot[i],  +1.0+0.0*I);
        ztp_2[i].ofts_smult_u(ztp_1[i], -1.0+0.0*I);
    }
    subCoef(Vdot, ztp_2);
    //zhdot = CQ*ztp_2 = CQ*(zdot - PCdot*zh - Vdot) @order m
    smvprod_u(CQ, ztp_2, zhdot);

}


/**
 *  \brief Apply the inverse derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 **/
void applyInvDotCOC(matrix<Ofsc> &CQ,
                    matrix<Ofsc> &PCdot,
                    vector<Ofsc> &Vdot,
                    vector<Oftsc> &zh,
                    vector<Oftsc> &zdot,
                    vector<Oftsc> &zhdot,
                    vector<Oftsc> &ztp_1,
                    vector<Oftsc> &ztp_2,
                    int m)
{
    //ztp_1 += PCdot*zh @order m
    smvprod_u(PCdot, zh, ztp_1, m);
    //ztp_2 += zdot - PCdot*zh - Vdot @order m
    for(unsigned int i = 0; i < zh.size(); i++)
    {
        ztp_2[i].ofts_smult_u(zdot[i],  +1.0+0.0*I, m);
        ztp_2[i].ofts_smult_u(ztp_1[i], -1.0+0.0*I, m);
    }
    if(m == 0) subCoef(Vdot, ztp_2);
    //zhdot = CQ*ztp_2 = CQ*(zdot - PCdot*zh - Vdot) @order m
    smvprod_u(CQ, ztp_2, zhdot, m);
}


/**
 *  \brief Apply the inverse derivative of the change of variables at order m in zIN/zOut.
 *         The change of variables is of the form: zhdot = CQ (zdot - PCdot zh - Vdot)
 *         This is the TFS version of this routine (time domain).
 **/
void tfts_applyInvDotCOC(matrix<Ofsc> &CQ,
                         matrix<Ofsc> &PCdot,
                         vector<Ofsc> &Vdot,
                         vector<Oftsc> &zh,
                         vector<Oftsc> &zdot,
                         vector<Oftsc> &zhdot,
                         vector<Oftsc> &ztp_1,
                         vector<Oftsc> &ztp_2,
                         int m)
{
    //ztp_1 += PCdot*zh @order m
    tfts_smvprod_u(PCdot, zh, ztp_1, m);
    //ztp_2 += zdot - PCdot*zh - Vdot @order m
    for(unsigned int i = 0; i < zh.size(); i++)
    {
        ztp_2[i].tfts_smult_u(zdot[i],  +1.0+0.0*I, m);
        ztp_2[i].tfts_smult_u(ztp_1[i], -1.0+0.0*I, m);
    }
    if(m == 0) tfts_subCoef(Vdot, ztp_2);
    //zhdot = CQ*ztp_2 = CQ*(zdot - PCdot*zh - Vdot) @order m
    tfts_smvprod_u(CQ, ztp_2, zhdot, m);
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          OFS version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Retrieves the Fouriers series within the matrices P and Q, the vectors V, Xe, Xm and Xs and the scalars ILe, ILs and ILm.
 *
 *  - We recall that the change of coordinates between the Normalized-Centered coordinates Z and the Translated-Floquet-Complexified (TFC)
 *  coordinates z is of the form:
 *
 *                  Z = P(t)C z + V(t)  <==> z = Q(t) (z - V(t))  and  Q = inv(P)
 *
 *  - For the vector field, it comes:
 *
 *                  dot(Z) = dot(PC) z + PC dot(z) + dot(V)
 *
 *  In this routine,
 *
 *  1. First step
 *      * P, Q and V are retrieved from txt files set in the folder \c qbcp_l.F_COC.
 *
 *  2. Second step
 *  The vectors Xe[2], Xm[2], and X[s] contains the time-dependent positions of the primaries in the xy-plane.
 *  Note that this is the TRANSLATED versions of these positions i.e. they are taking into account the motion
 *  of the periodic orbit, dynamical equivalent of the libration point L1,2.
 *
 *  The scalars of the form ILc, c = e, m, s are given by: 1/sqrt(Xc[0]^2 + Xc[1]^2)
 *  i.e. inverse of the euclidian distance from the origin, given as Fourier series.
 *
 *  NOTE: This routine is a simplified version of the routine of the same name, used for OFTS computations. It is only used for OFS version of the coc
 *  which itself is only used as a test version of the OFTS coc.
 *
 **/
void initCOC_OFS(matrix<Ofsc> &P,
                 matrix<Ofsc> &Q,
                 vector<Ofsc> &V,
                 vector<Ofsc> &Xe,
                 vector<Ofsc> &Xm,
                 vector<Ofsc> &Xs,
                 Ofsc &ILe,
                 Ofsc &ILm,
                 Ofsc &ILs,
                 QBCP_L& qbcp_l)
{
    //Retrieve folder
    string F_COC = qbcp_l.cs.F_COC;
    //Switch RTBP/QBCP
    if(qbcp_l.model == M_QBCP || qbcp_l.model == M_BCP)
    {
        //Recovering the data: matrix P
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*P.getCA(i,j), F_COC+"P", i+1, j+1);
        //Recovering the data: matrix Q
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*Q.getCA(i,j), F_COC+"Q", i+1, j+1);

        //Recovering the data: vector V
        //V = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:V[2] and V[5] are kept null
        readCOC(V[0], F_COC+"G1_",  1, 1);
        readCOC(V[1], F_COC+"G1_",  1, 2);
        readCOC(V[3], F_COC+"G1_",  2, 1);
        readCOC(V[4], F_COC+"G1_",  2, 2);

        //Recovering the data: vector Xm
        readCOC(Xm[0], F_COC+"Xm",  1, 1);
        readCOC(Xm[1], F_COC+"Xm",  2, 1);
        //Note: Xm[2] is kept null

        //Recovering the data: vector Xe
        readCOC(Xe[0], F_COC+"Xe",  1, 1);
        readCOC(Xe[1], F_COC+"Xe",  2, 1);
        //Note: Xe[2] is kept null

        //Recovering the data: vector Xs
        readCOC(Xs[0], F_COC+"Xs",  1, 1);
        readCOC(Xs[1], F_COC+"Xs",  2, 1);
        //Note: Xs[2] is kept null

        //Recovering the data: Ofts ILe, ILm, ILs
        readCOC(ILe, F_COC+"Xe",  3, 1);
        readCOC(ILm, F_COC+"Xm",  3, 1);
        readCOC(ILs, F_COC+"Xs",  3, 1);

    }
    else //EM RTPB
    {
        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, dl1, do1, s1, s2, c2;
        c2 = qbcp_l.cs.c2;
        eta1 = (c2 - 2 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        dl1 = 2*la1*( (4+3*c2)*la1*la1 + 4 + 5*c2 - 6*c2*c2);
        do1 =   om1*( (4+3*c2)*om1*om1 - 4 - 5*c2 + 6*c2*c2);
        s1  = sqrt(dl1);
        s2  = sqrt(do1);

        //--------------------
        //Init P
        //--------------------
        P.setCoef(+2*la1/s1,                          0, 1);
        P.setCoef(+(la1*la1 - 2*c2 - 1)/s1,           1, 1);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        P.setCoef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        P.setCoef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        P.setCoef(-(om1*om1  - 2*c2 - 1)/s2,          3, 0);

        P.setCoef(-2*la1/s1,                          0, 4);
        P.setCoef(+(la1*la1 - 2*c2 - 1)/s1,           1, 4);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        P.setCoef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        P.setCoef(+2*om1/s2,                          0, 3);
        P.setCoef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        P.setCoef(+1.0/sqrt(om2),                     2, 2);
        P.setCoef(sqrt(om2),                          5, 5);

        //--------------------
        //Q = inv(P) (gsl lib)
        //--------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_calloc (NV, NV);
        gsl_matrix* Qc   = gsl_matrix_calloc (NV, NV);
        gsl_permutation * p6 = gsl_permutation_alloc (NV);

        //Init Pc
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) gsl_matrix_set(Pc, i, j, creal(P.getCA(i,j)->ofs_getCoef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);


        //--------------------
        // Init Q
        //--------------------
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) Q.setCoef(gsl_matrix_get(Qc, i, j), i, j);


        //--------------------
        // Init Xm
        //--------------------
        double gamma = qbcp_l.cs.gamma;

        //Sun does not exist
        Xs[0].setCoef(0.0+0.0*I, 0);
        Xs[1].setCoef(0.0+0.0*I, 0);
        ILs.setCoef(0.0+0.0*I, 0);

        //Earth and Moon
        if(qbcp_l.li_EM == 1)
        {
            Xe[0].setCoef((-1.0+gamma)/gamma+0.0*I, 0);
            Xe[1].setCoef(0.0+0.0*I, 0);
            ILe.setCoef(1.0/sqrt( pow((-1.0+gamma)/gamma, 2.0))+0.0*I, 0);

            Xm[0].setCoef( +1.0+0.0*I, 0);
            Xm[1].setCoef(  0.0+0.0*I, 0);
            ILm.setCoef(    1.0+0.0*I, 0);
        }
        else
        {
            Xe[0].setCoef((-1.0-gamma)/gamma+0.0*I, 0);
            Xe[1].setCoef( 0.0+0.0*I, 0);
            ILe.setCoef( 1.0/sqrt( pow((-1.0-gamma)/gamma, 2.0))+0.0*I, 0);

            Xm[0].setCoef( -1.0+0.0*I, 0);
            Xm[1].setCoef(  0.0+0.0*I, 0);
            ILm.setCoef( 1.0+0.0*I, 0);
        }
    }
}

/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC_OFS(matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  vector<Ofsc> &zIn,
                  vector<Ofsc> &zOut)
{
    //zeroing the target
    for(unsigned int i = 0; i < zOut.size(); i++) zOut[i].zero();
    //zOut = PC*zIn
    smvprod_ofs(PC, zIn, zOut);
    //zOut+=V(theta)
    for(int i = 0; i < (int) zOut.size(); i++) zOut[i] += V[i];
}


///**
// *  \brief Apply the change of variables in zIN/zOut (OFS version).
// *       The change of variables is of the form: zOut = (P*C) zIN + V
// *
// *   This routine shows that using P as it is, without computing PC in the preprocess,
// *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX in the src code).
// **/
//void applyCOC_OFS(matrix<Ofsc> &P, vector<Ofsc> &V, vector<Ofsc> &zIn, vector<Ofsc> &zOut)
//{
//    //-----------------------------------------
//    //Temporary variables
//    //-----------------------------------------
//    Ofsc AUX(OFS_ORDER);
//    Ofsc BUX(OFS_ORDER);
//    Ofsc CUX(OFS_ORDER);
//    Ofsc DUX(OFS_ORDER);
//    //-----------------------------------------
//
//    //Loop on the components of zOut
//    for(int i = 0; i < NV ; i ++)
//    {
//        if(i== 2 || i == 5) //z variables are decoupled from the others (beware of the C++ shift!)
//        {
//            //-----------------------------------------
//            //Z variables
//            //-----------------------------------------
//            //AUX = 1/sqrt(2)*(pi3 + pi6*I)
//            AUX.ofs_fsum(P(i,2), 1.0/sqrt(2)+0.0*I, P(i,5), I*1.0/sqrt(2));
//            //CUX = 1/sqrt(2)*(I*pi3 + pi6)
//            CUX.ofs_fsum(P(i,2), I*1.0/sqrt(2), P(i,5), 1.0/sqrt(2)+0.0*I);
//
//            //zOut[i] =   1/sqrt(2)*(pi3 + pi6*I)*x3
//            //          + 1/sqrt(2)*(I*pi3 + pi6)*y3
//            zOut[i].ofs_prod(AUX, zIn[2]);
//            zOut[i].ofs_sprod(CUX, zIn[5]);
//        }
//        else
//        {
//            //-----------------------------------------
//            //X, Y variables
//            //-----------------------------------------
//            //AUX = 1/sqrt(2)*(pi1 + pi4*I)
//            AUX.ofs_fsum(P(i,0), 1.0/sqrt(2)+0.0*I, P(i,3), I*1.0/sqrt(2));
//            //BUX = pi2
//            BUX.ccopy(P(i,1));
//            //CUX = 1/sqrt(2)*(I*pi1 + pi4)
//            CUX.ofs_fsum(P(i,0), I*1.0/sqrt(2), P(i,3), 1.0/sqrt(2)+0.0*I);
//            //DUX = pi5
//            DUX.ccopy(P(i,4));
//
//            //zOut[i] =   1/sqrt(2)*(p11 + p14*I)*x1
//            //          + p12*x2
//            //          + 1/sqrt(2)*(I*p11 + p14)*y1
//            //          + Vi
//            //Vi is added!!!
//            zOut[i].ofs_prod(zIn[0], AUX);
//            zOut[i].ofs_sprod(zIn[1], BUX);
//            zOut[i].ofs_sprod(zIn[3], CUX);
//            zOut[i].ofs_sprod(zIn[4], DUX);
//            zOut[i] += V[i];
//
//        }
//    }
//}

/**
 *  \brief Apply the change of variables (without the order zero) in zIN/zOut (OFS version).
 *       The change of variables is of the form: zOut = (P*C) zIN
 *
 *   As applyCOC_OFS, this routine shows that using P as it is, without computing PC in the preprocess,
 *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX in the src code).
 **/
void applyModifiedCOC_OFS(matrix<Ofsc> &P, vector<Ofsc> &V, vector<Ofsc> &zIn, vector<Ofsc> &zOut)
{
    //-----------------------------------------
    //Temporary variables
    //-----------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);
    Ofsc DUX(OFS_ORDER);
    Ofsc EUX(OFS_ORDER);

    Ofsc x1(OFS_ORDER);
    Ofsc x2(OFS_ORDER);
    //-----------------------------------------

    for(int i = 0; i < NV ; i ++)
    {
        if(i== 2 || i == 5) //z variables are decoupled from the others (beware of the C++ shift!)
        {
            //-----------------------------------------
            //Z variables
            //-----------------------------------------
            //AUX = 1/sqrt(2)*(pi3 + pi6*I)
            AUX.ofs_fsum(P(i,2), 1.0/sqrt(2)+0.0*I, P(i,5), I*1.0/sqrt(2));
            //CUX = 1/sqrt(2)*(I*pi3 + pi6)
            CUX.ofs_fsum(P(i,2), I*1.0/sqrt(2), P(i,5), 1.0/sqrt(2)+0.0*I);

            //zOut[i] =   1/sqrt(2)*(pi3 + pi6*I)*x3
            //          + 1/sqrt(2)*(I*pi3 + pi6)*y3
            zOut[i].ofs_prod(AUX, zIn[2]);
            zOut[i].ofs_sprod(CUX, zIn[5]);
        }
        else
        {
            //-----------------------------------------
            //X, Y variables
            //-----------------------------------------
            //AUX = 1/sqrt(2)*(pi1 + pi4*I)
            AUX.ofs_fsum(P(i,0), 1.0/sqrt(2)+0.0*I, P(i,3), I*1.0/sqrt(2));
            //BUX = pi2
            BUX.ccopy(P(i,1));
            //CUX = 1/sqrt(2)*(I*pi1 + pi4)
            CUX.ofs_fsum(P(i,0), I*1.0/sqrt(2), P(i,3), 1.0/sqrt(2)+0.0*I);
            //DUX = pi5
            DUX.ccopy(P(i,4));
            //EUX = Vi
            EUX.ccopy(V[i]);

            //zOut[i] =   1/sqrt(2)*(p11 + p14*I)*x1
            //          + p12*x2
            //          + 1/sqrt(2)*(I*p11 + p14)*y1
            //          + p15*y2
            //Vi is not added!!!
            zOut[i].ofs_prod(zIn[0], AUX);
            zOut[i].ofs_sprod(zIn[1], BUX);
            zOut[i].ofs_sprod(zIn[3], CUX);
            zOut[i].ofs_sprod(zIn[4], DUX);
            //zOut[i].addCoef(EUX, 0, 0);

        }
    }
}

/**
 *  \brief Apply the inverse change of variables in zIN/zOut (OFS version).
 *         The change of variables is of the form: zOut = CQ*(zIn - V).
 *
 *   This routine shows that using Q as it is, without computing CQ in the preprocess,
 *   requires to use a lot of side variables (see the use of AUX, BUX, CUX, DUX, EUX in the src code).
 **/
void applyInvCOC_OFS(matrix<Ofsc>& Q, vector<Ofsc>& V, vector<Ofsc>& zIn, vector<Ofsc>& zOut)
{
    //-----------------------------------------
    //Temporary variables
    //-----------------------------------------
    Ofsc AUX(OFS_ORDER);
    Ofsc BUX(OFS_ORDER);
    Ofsc CUX(OFS_ORDER);
    Ofsc DUX(OFS_ORDER);
    Ofsc EUX(OFS_ORDER);
    //-----------------------------------------


    //-----------------------------------------
    //X
    //-----------------------------------------
    AUX.ofs_fsum(Q(0,0), 1.0/sqrt(2)+0.0*I, Q(3,0), -1.0/sqrt(2)*I);
    BUX.ofs_fsum(Q(0,1), 1.0/sqrt(2)+0.0*I, Q(3,1), -1.0/sqrt(2)*I);
    CUX.ofs_fsum(Q(0,3), 1.0/sqrt(2)+0.0*I, Q(3,3), -1.0/sqrt(2)*I);
    DUX.ofs_fsum(Q(0,4), 1.0/sqrt(2)+0.0*I, Q(3,4), -1.0/sqrt(2)*I);
    zOut[0].ofs_prod( zIn[0]-V[0], AUX);
    zOut[0].ofs_sprod(zIn[1]-V[1], BUX);
    zOut[0].ofs_sprod(zIn[3]-V[3], CUX);
    zOut[0].ofs_sprod(zIn[4]-V[4], DUX);


    //-----------------------------------------
    //Y
    //-----------------------------------------
    zOut[1].ofs_prod( zIn[0]-V[0], Q(1,0));
    zOut[1].ofs_sprod(zIn[1]-V[1], Q(1,1));
    zOut[1].ofs_sprod(zIn[3]-V[3], Q(1,3));
    zOut[1].ofs_sprod(zIn[4]-V[4], Q(1,4));


    //-----------------------------------------
    //Z
    //-----------------------------------------
    AUX.ofs_fsum(Q(2,2),    1.0/sqrt(2)+0.0*I, Q(5,2), -1.0/sqrt(2)*I);
    CUX.ofs_fsum(Q(2,5),    1.0/sqrt(2)+0.0*I, Q(5,5), -1.0/sqrt(2)*I);
    zOut[2].ofs_prod(AUX,  zIn[2]);
    zOut[2].ofs_sprod(CUX, zIn[5]);


    //-----------------------------------------
    //PX
    //-----------------------------------------
    AUX.ofs_fsum(Q(0,0), -1.0/sqrt(2)*I, Q(3,0), 1.0/sqrt(2)+0.0*I);
    BUX.ofs_fsum(Q(0,1), -1.0/sqrt(2)*I, Q(3,1), 1.0/sqrt(2)+0.0*I);
    CUX.ofs_fsum(Q(0,3), -1.0/sqrt(2)*I, Q(3,3), 1.0/sqrt(2)+0.0*I);
    DUX.ofs_fsum(Q(0,4), -1.0/sqrt(2)*I, Q(3,4), 1.0/sqrt(2)+0.0*I);
    zOut[3].ofs_prod( zIn[0]-V[0], AUX);
    zOut[3].ofs_sprod(zIn[1]-V[1], BUX);
    zOut[3].ofs_sprod(zIn[3]-V[3], CUX);
    zOut[3].ofs_sprod(zIn[4]-V[4], DUX);


    //-----------------------------------------
    //PY
    //-----------------------------------------
    zOut[4].ofs_prod( zIn[0]-V[0], Q(4,0));
    zOut[4].ofs_sprod(zIn[1]-V[1], Q(4,1));
    zOut[4].ofs_sprod(zIn[3]-V[3], Q(4,3));
    zOut[4].ofs_sprod(zIn[4]-V[4], Q(4,4));


    //-----------------------------------------
    //PZ
    //-----------------------------------------
    AUX.ofs_fsum(Q(2,2),    -1.0/sqrt(2)*I, Q(5,2), 1.0/sqrt(2)+0.0*I);
    CUX.ofs_fsum(Q(2,5),    -1.0/sqrt(2)*I, Q(5,5), 1.0/sqrt(2)+0.0*I);
    zOut[5].ofs_prod(AUX,  zIn[2]);
    zOut[5].ofs_sprod(CUX, zIn[5]);

}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Scalar version of the coc
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Apply the change of variables in zIN/zOut (scalar version).
 *       The change of variables is of the form: zOut = (P*C) zIN + V
 */
void applyCOC(matrix<Ofsc> PC, vector<Ofsc> V, vector<cdouble> &zIn, vector<cdouble> &zOut,  double t)
{
    //zOut = PC*zIn
    mvprod_u(PC, zIn, zOut, t);
    //zOut+=V(theta)
    vvsum_u(V, zOut, t);
}


/**
 *  \brief Apply the derivative change of variables at every order in zIN/zOut (scalar version).
 */
void applyDotCOC(matrix<Ofsc> PC,
                 matrix<Ofsc> PCdot,
                 vector<Ofsc> Vdot,
                 vector<cdouble> &zh,
                 vector<cdouble> &zhdot,
                 vector<cdouble> &zdot,
                 double t)
{
    //zdot = PCdot*zh
    mvprod_u(PCdot, zh, zdot, t);
    //zdot = PC*zhdot
    smvprod_u(PC, zhdot, zdot, t);
    //zdot+=Vdot(theta)
    vvsum_u(Vdot, zdot, t);
}

/**
 *  \brief Apply the inverse of the derivative change of variables in zIN/zOut (scalar version).
 */
void applyInvDotCOC(matrix<Ofsc> CQ,
                    matrix<Ofsc> PCdot,
                    vector<Ofsc> Vdot,
                    vector<cdouble> &zh,
                    vector<cdouble> &zdot,
                    vector<cdouble> &zhdot,
                    vector<cdouble> &ztp,
                    double t)
{
    //ztp = PCdot*zh
    mvprod_u(PCdot, zh, ztp, t);
    //zhdot = zdot - ztp =  zdot - PCdot*zh
    for(unsigned int i = 0; i< zh.size(); i++) zhdot[i] =  zdot[i] - ztp[i];
    //ztp = zhdot - Vdot =  zdot - PCdot*zh - Vdot
    vvsub_u(Vdot, zhdot, ztp, t);
    //zhdot = CQ*ztp = CQ*(zdot - PCdot*zh - Vdot)
    mvprod_u(CQ, ztp, zhdot, t);
}


/**
 *  \brief Apply the inverse of the change of variables at every order in zIN/zOut (scalar version).
 */
void applyInvCOC(matrix<Ofsc> CQ, vector<Ofsc> V, vector<cdouble> &zIn, vector<cdouble> &zOut, vector<cdouble> &ztp, double t)
{
    //-----------------------------------------
    //ztp = zIn - V
    //-----------------------------------------
    vvsub_u(V, zIn, ztp, t);
    //-----------------------------------------
    //zOut = CQ*ztp
    //-----------------------------------------
    mvprod_u(CQ, ztp, zOut, t);
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Tests
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test routine for the change of coordinates: comparison between an object x and COC^{-1}(COC(x)). Tested on OFS and OFTS objects.
 **/
void testCOC()
{
    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);

    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);

    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Test of the OFS routines
    //------------------------------------------
    cout << "-------------------------------------------------" << endl;
    cout << " Test routine for the COC applied on OFS objects " << endl;
    cout << "-------------------------------------------------" << endl;
    cout << std::showpos << setprecision(5) << setiosflags(ios::scientific);

    double maxAbs[2];
    double xMax = 1;          //arbitrary value
    vector<Ofsc> Xo(NV);      //coordinates in NC framework
    vector<Ofsc> Xoc(NV);     //coordinates in diagonalized framework
    vector<Ofsc> Xod(NV);     //coordinates in diagonalized framework
    vector<Ofsc> DeltaXo(NV); //coordinates in diagonalized framework

    //All coordinates to xMax
    for(int i = 0; i < NV; i++) Xoc[i].setCoef(xMax+0.0*I, 0);

    cout << " All coordinates in diagonalized framework are set to " << xMax << endl;
    cout << " We apply the COC and the COC^{-1} " << endl;
    cout << " Comparison with initial coordinates in diagonalized framework" << endl;
    cout << "-------------------------------------------------" << endl;
    cout << " Discrepancy at order zero and maximal value:" << endl;
    cout << " Component      Order zero          Max coef" << endl;

    //Applying the COC
    //applyCOC_OFS(P, V, Xoc, Xo);    //Xo = COC(Xoc)
    applyCOC_OFS(PC, V, Xoc, Xo);         //Xo = COC(Xoc), with the use of PC instead of P
    applyInvCOC_OFS(Q, V, Xo, Xod);   //Xod = COC^(-1)(Xo)

    //Checking the numerical discrepancy
    for(int i = 0; i<NV; i++)
    {

        DeltaXo[i].ofs_fsum(Xoc[i], 1.0+0.0*I, Xod[i], -1.0+0.0*I);
        DeltaXo[i].getCoefMaxNorm(maxAbs);

        cout << "   " << i << "         "  << cabs(DeltaXo[i].ofs_getCoef(0)) << "     "
             <<  maxAbs[0] << "(" <<  resetiosflags( ios::floatfield )
             << maxAbs[1] << setiosflags(ios::scientific) << ")" << endl;
    }

    //------------------------------------------
    // Test of the OFTS routines
    //------------------------------------------
    cout << endl;
    cout << "------------------------------------------------- "  << endl;
    cout << " Test routine for the COC applied on OFTS objects "  << endl;
    cout << "------------------------------------------------- "  << endl;
    cout << " All coordinates in diagonalized framework are set to random OFS" << endl;
    cout << " We apply the COC and the COC^{-1} "                 << endl;
    cout << " Comparison with initial coordinates in diagonalized framework" << endl;
    cout << "-------------------------------------------------"   << endl;
    cout << " Discrepancy at order zero and maximal value:"       << endl;
    cout << " Component      Order zero          Max coef"        << endl;

    vector<Oftsc> W(NV);       //coordinates in NC framework
    vector<Oftsc> Wc(NV);      //coordinates in diagonalized framework
    vector<Oftsc> Wd(NV);      //coordinates in diagonalized framework
    vector<Oftsc> DeltaW(NV);  //coordinates in diagonalized framework
    vector<Oftsc> ztp(NV);     //spare vector<OFTS>

    //All components set to random coefficients
    for(int i = 0; i < NV; i++) Wc[i].setRandomCoefs();

    //Applying the COC
    applyCOC(PC,  V,  Wc, W);          //order zero is added
    //applyCOC(P, V, Wc, W);            //Xo = COC(Xoc)
    applyInvCOC(CQ,  V,  W,  Wd, ztp);  //order zero is added
    //applyInvCOC(Q, V, W, Wd, ztp);     //Xod = COC^(-1)(Xo)

    //Checking the numerical discrepancy
    for(int i = 0; i<NV; i++)
    {
        DeltaW[i].ofts_fsum_u(Wc[i], 1.0+0.0*I, Wd[i], -1.0+0.0*I);
        DeltaW[i].getCoef(0, 0)->getCoefMaxNorm(maxAbs);

        cout << "   " << i << "         "  << cabs(DeltaW[i].getCoef(0, 0)->ofs_getCoef(0))
             << "     " <<  maxAbs[0] << "(" << resetiosflags( ios::floatfield )
             << maxAbs[1] << setiosflags(ios::scientific) << ")" << endl;
    }

}

/**
 *  \brief Test routine for the derivatives of change of coordinates. Tested OFTS objects.
 **/
void testDotCOC()
{
    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Test of the OFTS routines
    //------------------------------------------
    cout << std::showpos << setprecision(5) << setiosflags(ios::scientific);
    cout << "------------------------------------------------- "  << endl;
    cout << " Test routine for the COC applied on OFTS objects "  << endl;
    cout << "------------------------------------------------- "  << endl;
    cout << " All coordinates in diagonalized framework are set to random OFS" << endl;
    cout << " We apply the COC and the COC^{-1} "                 << endl;
    cout << " Comparison with initial coordinates in diagonalized framework" << endl;
    cout << "-------------------------------------------------"   << endl;
    cout << " Discrepancy at order zero and maximal value:"       << endl;
    cout << " Component      Order zero          Max coef"        << endl;

    vector<Oftsc> z(NV);      //NC
    vector<Oftsc> zdot(NV);   //NC
    vector<Oftsc> zdot2(NV);  //NC
    vector<Oftsc> zh(NV);     //TFC
    vector<Oftsc> zhdot(NV);  //TFC
    vector<Oftsc> DeltaW(NV); //Delta in NC
    vector<Oftsc> ztp1(NV);    //spare vector<OFTS>
    vector<Oftsc> ztp2(NV);    //spare vector<OFTS>

    //All components set random coefficients
    for(int i = 0; i < NV; i++) z[i].setRandomCoefs();
    for(int i = 0; i < NV; i++) zdot[i].setRandomCoefs();

    //zh = inCOC(z)
    applyInvCOC(CQ,  V,  z,  zh, ztp1);

    for(unsigned int i = 0; i < zhdot.size(); i++) ztp1[i].zero();
    //zhdot = inCOCdot(zdot)
    //applyInvDotCOC(CQ, PCdot, Vdot, zh, zdot, zhdot, ztp1, ztp2);
    for(int k = 0 ; k<= OFTS_ORDER; k++) applyInvDotCOC(CQ, PCdot, Vdot, zh, zdot, zhdot, ztp1, ztp2, k);
    //zdot2 = COCdot(zhdot)
    applyDotCOC(PC, PCdot, Vdot, zh, zhdot, zdot2);


    double maxAbs[2];
    //Checking the numerical discrepancy
    for(int i = 0; i<NV; i++)
    {
        DeltaW[i].ofts_fsum_u(zdot[i], 1.0+0.0*I, zdot2[i], -1.0+0.0*I);
        DeltaW[i].getCoef(0, 0)->getCoefMaxNorm(maxAbs);

        cout << "   " << i << "         "  << cabs(DeltaW[i].getCoef(0, 0)->ofs_getCoef(0))
             << "     " <<  maxAbs[0] << "(" << resetiosflags( ios::floatfield )
             << maxAbs[1] << setiosflags(ios::scientific) << ")" << endl;
    }

}

/**
 *  \brief Test routine for the change of coordinates on a trajectory integrated on a full period.
 *
 *      The same initial conditions are integrated both in NC and TFC coordinates on a full period T, and the results are compared.
 *      Note that the vector field in TFC makes use of the scalar version of the COC, so this test only proves that the scalar COC is self-consistent
 *      but no more than that.
 **/
void testIntCOC()
{
    double theta1 = SEML.us.T;
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the TFC c.o.c.              " << endl;
    cout << "                                                   " << endl;
    cout << "        Integration both with NC and TFC v.f.      " << endl;
    cout << "                    Comparison.                    " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;

    cout << std::showpos << setprecision(5) << setiosflags(ios::scientific);
    cout << "Arbitrary initial conditions z(0) = z0 are set in NC coordinates." << endl;
    cout << "Then, zh0 = COC(z0) is used to initialize the integration in TFC." << endl;
    cout << "Integration is taken during theta1 = " << theta1/SEML.us.n << " in NC time units (eq. to the system time units, e.g Earth-Moon or Sun-Earth)" << endl;
    cout << "---------------------------------------------------" << endl;

    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //-------------------------------------------------
    //Integration tools
    //-------------------------------------------------
    //For dot(z) = F(z)
    //-------------------
    gsl_odeiv2_system sys;
    sys.function = qbfbp_vfn_novar;
    sys.jacobian = NULL;
    sys.dimension = NV;
    sys.params = &SEML;
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, 1e-6, 1e-14, 1e-14);

    //For dot(zh) = Fh(zh)
    //-------------------
    //COC structure to use in qbfbp_Fh
    vector<cdouble> zh(6);
    vector<cdouble> z(6);
    vector<cdouble> zd(6);
    vector<cdouble> ztp(6);
    vector<cdouble> zhd(6);
    COC coc;
    coc.PC      = &PC;
    coc.PCdot   = &PCdot;
    coc.CQ      = &CQ;
    coc.V       = &V;
    coc.Vdot    = &Vdot;
    coc.qbcp_l  = &SEML;
    coc.zh      = &zh;
    coc.z       = &z;
    coc.zd      = &zd;
    coc.zhd     = &zhd;
    coc.ztp     = &ztp;
    //GSL tools
    gsl_odeiv2_system sys_h;
    sys_h.function = qbfbp_Fh;
    sys_h.jacobian = NULL;
    sys_h.dimension = 12;
    sys_h.params = &coc;
    const gsl_odeiv2_step_type *T_h = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_h = gsl_odeiv2_driver_alloc_y_new (&sys_h, T_h, 1e-6, 1e-14, 1e-14);

    //Time dependencies
    double tinit = 0.0;
    double t;

    //-------------------------------------------------
    //initialization in NC
    //-------------------------------------------------
    double az0[6];
    vector<cdouble> z0(6);
    vector<cdouble> z0_v2(6);
    vector<cdouble> zdot0(6);
    vector<cdouble> zdot0_v2(6);
    for(int i = 0; i< 6; i++)
    {
        az0[i] = 1e-3;
        z0[i] = az0[i]+0.0*I;
        zdot0[i] = az0[i]+0.0*I;  //just for example
    }

    //-------------------------------------------------
    //initialization in TFC
    //-------------------------------------------------
    vector<cdouble> zh0(6);
    vector<cdouble> zhdot0(6);
    vector<cdouble> ztp2(6);
    double ayh0[12];

    //zh0 = COC-1(z0)
    applyInvCOC(CQ, V, z0, zh0, ztp2, SEML.us.n*tinit);
    for(int i =0; i<6; i++) ztp2[i] = 0+0.0*I;

    //Test of COC (uncomment if needed)
    //------------------------
    //applyCOC(PC, V, zh0, z0_v2, n*tinit);
    //cout << "z0_v2 = " << endl;
    //for(int i =0; i<6; i++) cout << creal(z0_v2[i]) << "+" << cimag(z0_v2[i]) << endl;

    //Test of dot(COC) (uncomment if needed)
    //------------------------
    //applyInvDotCOC(CQ, PCdot, Vdot, zh0, zdot0, zhdot0, ztp2, n*tinit);
    //applyDotCOC(PC,PCdot,Vdot, zh0, zhdot0, zdot0_v2 ,n*tinit);
    //cout << "zdot0 = " << endl;
    //for(int i =0; i<6; i++) cout << creal(zdot0[i]) << "+" << cimag(zdot0[i]) << endl;
    //cout << "zdot0_v2 = " << endl;
    //for(int i =0; i<6; i++) cout << creal(zdot0_v2[i]) << "+" << cimag(zdot0_v2[i]) << endl;

    //Update ayh0
    ayh0[0]  = creal(zh0[0]);
    ayh0[1]  = cimag(zh0[0]);
    ayh0[2]  = creal(zh0[1]);
    ayh0[3]  = cimag(zh0[1]);
    ayh0[4]  = creal(zh0[2]);
    ayh0[5]  = cimag(zh0[2]);
    ayh0[6]  = creal(zh0[3]);
    ayh0[7]  = cimag(zh0[3]);
    ayh0[8]  = creal(zh0[4]);
    ayh0[9]  = cimag(zh0[4]);
    ayh0[10] = creal(zh0[5]);
    ayh0[11] = cimag(zh0[5]);

    //-------------------------------------------------
    //Integration in NC
    //-------------------------------------------------
    t = tinit;
    //Integration
    gsl_odeiv2_driver_apply (d, &t, theta1/SEML.us.n, az0); //T = 2pi/n
    //On screen serie
    cout << "---------------------------------" << endl;
    cout << "NC final conditions: " << endl;
    for(int p = 0; p < 6; p++) cout << az0[p] << endl;
    vector<cdouble> z1(6);
    z1[0] = az0[0]+0.0*I;
    z1[1] = az0[1]+0.0*I;
    z1[2] = az0[2]+0.0*I;
    z1[3] = az0[3]+0.0*I;
    z1[4] = az0[4]+0.0*I;
    z1[5] = az0[5]+0.0*I;

    //-------------------------------------------------
    //Integration in TFC
    //-------------------------------------------------
    t = tinit;
    //Integration
    gsl_odeiv2_driver_apply (d_h, &t, theta1/SEML.us.n, ayh0); //T = 2pi/n

    //-------------------------------------------------
    //Back in NC coordinates
    //-------------------------------------------------
    vector<cdouble> zh1(6);
    vector<cdouble> zh1_NC(6);
    vector<cdouble> z1_TFC(6);
    zh1[0] = ayh0[0]  + ayh0[1]*I;
    zh1[1] = ayh0[2]  + ayh0[3]*I;
    zh1[2] = ayh0[4]  + ayh0[5]*I;
    zh1[3] = ayh0[6]  + ayh0[7]*I;
    zh1[4] = ayh0[8]  + ayh0[9]*I;
    zh1[5] = ayh0[10] + ayh0[11]*I;

    //z1_TFC = COC(zh1)
    applyCOC(PC, V, zh1, z1_TFC, SEML.us.n*(theta1/SEML.us.n)); //theta1 = 2pi = n*T = n*(2pi/n)
    //zh1_NC = invCOC(z1)
    applyInvCOC(CQ, V, z1, zh1_NC, ztp2, SEML.us.n*(theta1/SEML.us.n));

    //On screen serie
    cout << "---------------------------------" << endl;
    cout << "TFC final conditions back in NC: " << endl;
    for(int p = 0; p < 6; p++) cout << creal(z1_TFC[p]) << "+" << cimag(z1_TFC[p]) << endl;

    cout << "---------------------------------" << endl;
    cout << "Delta in NC: " << endl;
    for(int p = 0; p < 6; p++) cout << creal(z1_TFC[p]-z1[p]) << "+" << cimag(z1_TFC[p]-z1[p]) << endl;

    cout << "---------------------------------" << endl;
    cout << "Delta in TFC: " << endl;
    for(int p = 0; p < 6; p++) cout << creal(zh1_NC[p]-zh1[p]) << "+" << cimag(zh1_NC[p]-zh1[p]) << endl;
}

/**
 *  \brief Gives the vector field f[12] in TFC coordinates from the state y[12] in TFC coordinates.
 *
 *  The process is the following:
 *  - The state zh is given in the form y[i] = re(zh[i%2]), y[i+1] = Im(zh[i%2]) (real and imag part of the state are stored in y)
 *  - zh is build from y and the scalar version of the COC is applied : z = COC(zh)
 *  - The classical vector field of the QBCP is applied on z: zdot = qbfbp_vfn_novar(z)
 *  - The inverse of the derivatives of the coc is applied on the corresponding vector field: f = invCOCdot(zdot)
 *
 *  This routine is essentially used in testIntCOC to test the COC on a trajectory integrated on a full period.
 *  It makes use of the structure COC, defined in coc.h
 **/
int qbfbp_Fh(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------
    //Initialization
    //-------------------------------
    COC *coc = (COC*) params_void;
    double n = (coc->qbcp_l)->us.n;
    double y2[6];
    double f2[6];

    //-------------------------------
    //Update
    //-------------------------------
    //Reset all temp variables
    fill(coc->zh->begin(),  coc->zh->end(),  0.0);
    fill(coc->z->begin(),   coc->z->end(),   0.0);
    fill(coc->zd->begin(),  coc->zd->end(),  0.0);
    fill(coc->ztp->begin(), coc->ztp->end(), 0.0);
    fill(coc->zhd->begin(), coc->zhd->end(), 0.0);

    //Update zh = current state in TFC variables
    coc->zh->at(0) = y[0]  + I*y[1];
    coc->zh->at(1) = y[2]  + I*y[3];
    coc->zh->at(2) = y[4]  + I*y[5];
    coc->zh->at(3) = y[6]  + I*y[7];
    coc->zh->at(4) = y[8]  + I*y[9];
    coc->zh->at(5) = y[10] + I*y[11];

    //z = COC(zh) = current state in NC variables
    applyCOC(*(coc->PC), *(coc->V), *coc->zh, *coc->z, n*t);
    //Realification: y2 = real(z)
    for(int i = 0; i<6; i++) y2[i] = creal(coc->z->at(i));
    //Realification test.
    for(int i =0; i<6; i++) if(cimag(coc->z->at(i)) != 0) cout << "imag(z) is non null! " <<  cimag(coc->z->at(i)) << endl;

    //-------------------------------
    //Apply the VF in NC variables
    //-------------------------------
    //f2 = F(z2)
    qbfbp_vfn_novar(t, y2, f2, (coc->qbcp_l));

    //Update zd
    coc->zd->at(0) = f2[0]+0.0*I;
    coc->zd->at(1) = f2[1]+0.0*I;
    coc->zd->at(2) = f2[2]+0.0*I;
    coc->zd->at(3) = f2[3]+0.0*I;
    coc->zd->at(4) = f2[4]+0.0*I;
    coc->zd->at(5) = f2[5]+0.0*I;

    //-------------------------------
    //Update in TFC variables
    //-------------------------------
    //Back in TFC coordinates
    applyInvDotCOC(*(coc->CQ), *(coc->PCdot), *(coc->Vdot), *coc->zh, *coc->zd, *coc->zhd, *coc->ztp, n*t);

    //Update f
    f[0]  = creal(coc->zhd->at(0));
    f[1]  = cimag(coc->zhd->at(0));
    f[2]  = creal(coc->zhd->at(1));
    f[3]  = cimag(coc->zhd->at(1));
    f[4]  = creal(coc->zhd->at(2));
    f[5]  = cimag(coc->zhd->at(2));
    f[6]  = creal(coc->zhd->at(3));
    f[7]  = cimag(coc->zhd->at(3));
    f[8]  = creal(coc->zhd->at(4));
    f[9]  = cimag(coc->zhd->at(4));
    f[10] = creal(coc->zhd->at(5));
    f[11] = cimag(coc->zhd->at(5));

    return GSL_SUCCESS;
}

/**
 *  \brief Evaluation of Pij(nt) with *Pij contained in (void) *params.
 *
 *   Note that the use of pointers is required by the C, GSL, root-finding tool used in the routine pij.
 **/
double fpij(double t, void *params)
{
    Ofsc *pij = (Ofsc *) params;
    return creal(pij->evaluate(SEML.us.n*t));
}

/**
 *  \brief Finds the first positive root t0 of the equation Pij(nt) = 0. Note that ii and jj must be given between 0 and 5
 **/
double pij(int ii, int jj)
{
    cout << std::showpos << setprecision(15) << setiosflags(ios::scientific);
    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Bracketing the root
    //------------------------------------------
    int N = 50000;
    int i = 0;
    double ti;
    cdouble p36, pp36;
    p36 = P.getCA(ii,jj)->evaluate(0.0);
    do
    {
        pp36 = p36;
        ti = (double) i*2*M_PI/(N*SEML.us.n);
        p36 = P.getCA(ii,jj)->evaluate(SEML.us.n * ti);
        if(creal(p36)*creal(pp36) < 0)
        {
            cout << "A root of P(3,6) has been found around ti ~" << ti << endl;
            cout << "This time is equal to " << ti/(SEML.us.T/2) << "half periods." << endl;
        }
        i++;
    }
    while((i<N) && (creal(p36)*creal(pp36)>0));
    //At this point, the root is between t[i] and t[i-1]


    //------------------------------------------
    // Isolating the root
    //------------------------------------------
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;

    double x_lo = (double) (i-2)*2*M_PI/(N*SEML.us.n);
    double x_hi = (double) (i-1)*2*M_PI/(N*SEML.us.n);

    //Function for root finding, linked fo fp36 & parameter P
    gsl_function F;
    F.function = &fpij;
    F.params = P.getCA(ii,jj);

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    printf ("Root finding using %s method\n",
            gsl_root_fsolver_name (s));
    printf ("%5s [%9s, %9s] %9s %10s\n",
            "iter", "lower", "upper", "root", "err");
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);

        //Checking convergence
        status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-15);
        //or
        //x_lo = fp36(r, &P);
        //status = gsl_root_test_residual (x_lo , 1e-15);

        if (status == GSL_SUCCESS)
            printf ("Converged:\n");
        printf ("%5d [%.7f, %.7f] %.15f %.7e\n",
                iter, x_lo, x_hi,
                r,
                x_hi - x_lo);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    return r;

}

/**
 *  \brief Plot Pij(nt)
 **/
void pij_plot(gnuplot_ctrl  *h1)
{
    cout << std::showpos << setprecision(15) << setiosflags(ios::scientific);
    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Computing Pij(nt)
    //------------------------------------------
    int N = 5000;
    double ti, pij[36][N], tv[N];
    cdouble pijc;
    for(int i = 0; i < N ; i++)
    {
        ti = (double) 2*M_PI*i/N;
        for(int ii = 0; ii <6 ; ii++)
        {
            for(int jj = 0; jj <6 ; jj++)
            {
                pijc = P.getCA(ii,jj)->evaluate(ti);
                pij[ii*6+jj][i] = creal(pijc);
            }
        }
        tv[i] = ti;
    }


    //------------------------------------------
    // Mean value & Standard deviation
    //------------------------------------------
    double mu[36], sigma[36];

    for(int k = 0; k <36; k++)
    {
        //Mean value
        mu[k] = 0;
        for(int i = 0; i < N ; i++) mu[k] += pij[k][i]/N;

        //Standard deviation
        sigma[k] = 0;
        for(int i = 0; i < N ; i++) sigma[k] += (pij[k][i] - mu[k])*(pij[k][i] - mu[k])/N;
        sigma[k] = sqrt(sigma[k]);

        cout << "pij_plot. mean[" << k << "] = " << mu[k] << endl;
        cout << "pij_plot. sd[" << k << "]   = " << sigma[k] << endl;
    }


    //------------------------------------------
    // Max standard deviation
    //------------------------------------------
    double sdm = sigma[0];
    int argsdm = 0;
    for(int k = 1; k < 36; k++)
    {
        if(sdm < sigma[k])
        {
            sdm = sigma[k];
            argsdm = k;
        }
    }
    cout << "----------------------------" << endl;
    cout << "Maximum sd = " << sdm << endl;
    cout << "Obtained for k = " << argsdm << endl;

    //------------------------------------------
    // Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"t [-]");
    gnuplot_set_ylabel(h1, (char*)"alpha[argsdm](nt) [-]");

    int color = 1;
    gnuplot_plot_xy(h1, tv, pij[argsdm], N, (char*)"", "lines", "1", "1", color++);
}

/**
 *  \brief Plot Pij(nt)
 **/
void pij_plot(int ii, int jj, gnuplot_ctrl  *h1)
{
    cout << std::showpos << setprecision(15) << setiosflags(ios::scientific);
    //------------------------------------------
    // Initialization of the COC
    //------------------------------------------
    //The matrix P of the c.o.c. (Floquet part)
    matrix<Ofsc> P(NV,NV);
    matrix<Ofsc> Q(NV,NV);
    //The vector V of the c.o.c. (Translation part)
    vector<Ofsc> V(NV);
    //The vectors Xe, Xm, Xs, modified position of the primaries used to compute the potential
    vector<Ofsc> Xe(3);
    vector<Ofsc> Xm(3);
    vector<Ofsc> Xs(3);

    //The Ofs ILe, ILm and ILs, modified inverse orbit radii of the primaries used to compute the potential
    Ofsc ILe(OFS_ORDER);
    Ofsc ILm(OFS_ORDER);
    Ofsc ILs(OFS_ORDER);

    //PC=P(theta)*C
    matrix<Ofsc> PC(NV,NV);
    //PCdot=dot(P*C)=dot(P)C
    matrix<Ofsc> PCdot(NV,NV);
    //CQ=Cinv*Pinv=inv(PC)
    matrix<Ofsc> CQ(NV,NV);
    //Vdot=dot(V)
    vector<Ofsc> Vdot(NV);

    //Init routine
    initCOC(P, Q, PC, PCdot, CQ, Xe, Xm, Xs, V, Vdot, ILe, ILm, ILs, SEML);

    //------------------------------------------
    // Computing Pij(nt)
    //------------------------------------------
    int N = 5000;
    double ti, pij[N], tv[N];
    cdouble pijc;

    for(int i = 0; i < N ; i++)
    {
        ti = (double) 2*M_PI*i/N;
        pijc = P.getCA(ii,jj)->evaluate(ti);
        pij[i] = creal(pijc);
        tv[i] = ti;
    }

    //------------------------------------------
    // Mean value & Standard deviation
    //------------------------------------------
    double mu, sigma;

    //Mean value
    mu = 0;
    for(int i = 0; i < N ; i++) mu += pij[i]/N;

    //Standard deviation
    sigma = 0;
    for(int i = 0; i < N ; i++) sigma += (pij[i] - mu)*(pij[i] - mu)/N;
    sigma = sqrt(sigma);

    cout << "pij_plot. mean = " << mu << endl;
    cout << "pij_plot. sd   = " << sigma << endl;

    //------------------------------------------
    // Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"t [-]");
    gnuplot_set_ylabel(h1, (char*)"Pij(nt) [-]");
    gnuplot_plot_xy(h1, tv, pij, N, (char*)"", "lines", "1", "1", 1);
}

/**
 *  \brief Plot the coefficient of qbp
 **/
void coeff_plot(gnuplot_ctrl  *h1, QBCP_L* qbp)
{
    cout << std::showpos << setprecision(15) << setiosflags(ios::scientific);

    //Retrieving the parameters
    double n     = qbp->us.n;

    //-------------------------------------------------------------------------------
    //Plotting parameters
    //-------------------------------------------------------------------------------
    int N = 5000;
    int noc = qbp->numberOfCoefs;
    double ti, tv[N];
    double alpha[noc];
    double alphav[noc][N];

    //------------------------------------------
    // Evaluate the alphas
    //------------------------------------------
    for(int i = 0; i < N ; i++)
    {
        ti = (double) 2*M_PI*i/(n*N);
        evaluateCoef(alpha, ti, n, qbp->nf, qbp->cs.coeffs, noc);
        for(int k = 0; k < noc; k++) alphav[k][i] = alpha[k];
        tv[i] = ti;
    }

    //------------------------------------------
    // Mean value & Standard deviation
    //------------------------------------------
    double mu[noc], sigma[noc];

    for(int k = 0; k <noc; k++)
    {
        //Mean value
        mu[k] = 0;
        for(int i = 0; i < N ; i++) mu[k] += alphav[k][i]/N;

        //Standard deviation
        sigma[k] = 0;
        for(int i = 0; i < N ; i++) sigma[k] += (alphav[k][i] - mu[k])*(alphav[k][i] - mu[k])/N;
        sigma[k] = sqrt(sigma[k]);

        cout << "pij_plot. mean[" << k << "] = " << mu[k] << endl;
        cout << "pij_plot. sd[" << k << "]   = " << sigma[k] << endl;
    }


    //------------------------------------------
    // Max standard deviation
    //------------------------------------------
    double sdm = sigma[0];
    int argsdm = 0;
    for(int k = 1; k < noc; k++)
    {
        if(sdm < sigma[k])
        {
            sdm = sigma[k];
            argsdm = k;
        }
    }
    cout << "----------------------------" << endl;
    cout << "Maximum sd = " << sdm << endl;
    cout << "Obtained for k = " << argsdm << endl;

    //------------------------------------------
    // Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"t [-]");
    gnuplot_set_ylabel(h1, (char*)"alpha[argsdm](nt) [-]");

    int color = 1;
    gnuplot_plot_xy(h1, tv, alphav[argsdm], N, (char*)"", "lines", "1", "1", color++);
}

/**
 *  \brief Plot the potential of the primaries when the s/c is fixed at the Lagrange point. Deprecated.
 **/
void potential_plot(gnuplot_ctrl  *h1, QBCP_L* qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;


    double alpha[8], Ps[3], Pe[3], Pm[3];

    //-------------------------------------------------------------------------------
    //Plotting parameters
    //-------------------------------------------------------------------------------
    int N = 5000;
    double ti, tv[N];
    double pT[3][N];


    //------------------------------------------
    // Set the position at the origine of coordinates
    // i.e. @ the libration point
    //------------------------------------------
    double y[3];
    for(int i = 0; i <3; i++) y[i] = 0.0;

    //------------------------------------------
    // Loop on time
    //------------------------------------------
    for(int i = 0; i < N ; i++)
    {
        ti = (double) 2*M_PI*i/(n*N);
        //-------------------------------------------------------------------------------
        //Evaluate the alphas @t
        //-------------------------------------------------------------------------------
        evaluateCoef(alpha, ti, n, qbp->nf, qbp->cs.coeffs, 8);

        //-------------------------------------------------------------------------------
        //Evaluate the primaries positions @ t
        //-------------------------------------------------------------------------------
        evaluateCoef(Ps, ti, n, qbp->nf, qbp->cs.Ps, 3);
        evaluateCoef(Pe, ti, n, qbp->nf, qbp->cs.Pe, 3);
        evaluateCoef(Pm, ti, n, qbp->nf, qbp->cs.Pm, 3);

        //-------------------------------------------------------------------------------
        // Distances to 2nd power
        //-------------------------------------------------------------------------------
        double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
        double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
        double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

        //-------------------------------------------------------------------------------
        // Potential
        //-------------------------------------------------------------------------------
        pT[0][i] = /*- alpha[5]**/me/pow(qPe2, 1.0/2);
        pT[1][i] = /*- alpha[5]**/mm/pow(qPm2, 1.0/2);
        pT[2][i] = /*- alpha[5]**/ms/pow(qPs2, 1.0/2);
        tv[i] = ti;

    }



    //------------------------------------------
    // Mean value & Standard deviation
    //------------------------------------------
    double mu[3], sigma[3];

    for(int k = 0; k < 3; k++)
    {
        //Mean value
        mu[k] = 0;
        for(int i = 0; i < N ; i++) mu[k] += pT[k][i]/N;

        //Standard deviation
        sigma[k] = 0;
        for(int i = 0; i < N ; i++) sigma[k] += (pT[k][i] - mu[k])*(pT[k][i] - mu[k])/N;
        sigma[k] = sqrt(sigma[k]);

        cout << "pij_plot. mean[" << k << "] = " << mu[k] << endl;
        cout << "pij_plot. sd[" << k << "]   = " << sigma[k] << endl;
    }


    //------------------------------------------
    // Max standard deviation
    //------------------------------------------
    double sdm = sigma[0];
    int argsdm = 0;
    for(int k = 1; k < 3; k++)
    {
        if(sdm < sigma[k])
        {
            sdm = sigma[k];
            argsdm = k;
        }
    }
    cout << "----------------------------" << endl;
    cout << "Maximum sd = " << sdm << endl;
    cout << "Obtained for k = " << argsdm << endl;

    //------------------------------------------
    // Plotting
    //------------------------------------------
    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"t [-]");
    gnuplot_set_ylabel(h1, (char*)"pT[argsdm](nt) [-]");

    int color = 1;
    gnuplot_plot_xy(h1, tv, pT[argsdm], N, (char*)"", "lines", "1", "1", color++);
}


