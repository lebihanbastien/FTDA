/**
 * \file config.cpp
 * \brief Configuration file. It allows to initialize the environment of the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */
#include "config.h"
#include <gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

//---------------------------------------
//QBCP
//---------------------------------------
QBCP SEM;              ///global variable that describes the Sun-Earth-Moon system
QBCP_L SEML;           ///global variable that describes the Sun-Earth-Moon system around Li (EM framework)
QBCP_L SEML_SEM;       ///global variable that describes the Sun-Earth-Moon system around Li (SEM framework)

/**
 *   \brief Initialization of the environnement (Sun, Earth, Moon, Li...).
 *
 *    The global variables SEM and SEML are initialized in order to describe the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem around the Lagrangian point LI,
 *    LI being given as an integer in the file parameters.h. The default coordinates system is the normalized-centered (NC) coordinates of the Earth-Moon system.
 *    Note that, in order for the initialization to be complete - Sun-(Earth+Moon) motion given as Fourier series within SEML -
 *    the routine qbtbp(true) must have been run *once.
 **/
void init_env(int li_EM, int li_SEM, int isNormalized, int model, int fwrk, int pmStyle, int manType_EM, int manType_SEM)
{
    //Init the Sun-Earth-Moon problem
    init_QBCP(&SEM, SUN, EARTH, MOON, model);
    //Init the Sun-Earth-Moon problem focused on one libration point
    init_QBCP_L(&SEML, &SEM, isNormalized, li_EM, li_SEM, false, model, F_EM, pmStyle, manType_EM, manType_SEM);
    init_QBCP_L(&SEML_SEM, &SEM, isNormalized, li_EM, li_SEM, false, model, F_SEM, pmStyle, manType_EM, manType_SEM);
}


//---------------------------------------
//Center manifold
//---------------------------------------
vector<Oftsc>  CM;     ///center manifold in NC coordinates
vector<Oftsc>  CMdot;  ///time derivative of the center manifold in NC coordinates
vector<Oftsc> CMh;     ///center manifold in TFC coordinates
matrix<Oftsc> JCM;     ///jacobian of CM
vector<Oftsc>  Fh;     ///reduced vector field
vector<Oftsc>  DWFh;   ///JCM * Fh

/**
 *   \brief Initialization of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void initCM(QBCP_L &qbcp)
{
    //Memory allocation
    CM     = vector<Oftsc>(NV);
    CMdot  = vector<Oftsc>(NV);
    CMh    = vector<Oftsc>(NV);
    DWFh   = vector<Oftsc>(NV);
    Fh     = vector<Oftsc>(REDUCED_NV);
    JCM    = matrix<Oftsc>(NV,REDUCED_NV);

    //Update
    updateCM(qbcp);
}

/**
 *   \brief Update of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void updateCM(QBCP_L &qbcp)
{
    //Read from bin files
    readVOFTS_bin(CM,    qbcp.cs.F_PMS+"W/W",         OFS_ORDER);
    readVOFTS_bin(CMh,   qbcp.cs.F_PMS+"W/Wh",        OFS_ORDER);
    readVOFTS_bin(Fh,    qbcp.cs.F_PMS+"rvf/fh",      OFS_ORDER);
    readVOFTS_bin(DWFh,  qbcp.cs.F_PMS+"DWf/C_DWf",   OFS_ORDER);
    readVOFTS_bin(CMdot, qbcp.cs.F_PMS+"Wdot/C_Wdot", OFS_ORDER);
    //Building the Jacobian
    for(int i = 0 ; i < NV ; i++) for(int j = 0 ; j < REDUCED_NV; j++) JCM.der(CM[i], j+1, i, j);   //in NC
}


//---------------------------------------
//COC
//---------------------------------------
matrix<Ofsc>  Mcoc;    ///COC matrix
matrix<Ofsc>  Pcoc;    ///COC matrix (Mcoc = Pcoc*Complex matrix)
matrix<Ofsc>  MIcoc;   ///COC matrix = inv(Mcoc)
matrix<Ofsc>  PIcoc;   ///COC matrix = inv(Pcoc)
vector<Ofsc>  Vcoc;    ///COC vector

/**
 *  Main routine for the update of the COC
 **/
void updateCOC(QBCP_L &qbcp)
{
    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, qbcp);
}


/**
 *  Main routine for the initialization of the COC
 **/
void initCOC(QBCP_L &qbcp)
{
    //Memory allocation
    Pcoc  = matrix<Ofsc>(NV,NV);
    Mcoc  = matrix<Ofsc>(NV,NV);
    PIcoc = matrix<Ofsc>(NV,NV);
    MIcoc = matrix<Ofsc>(NV,NV);
    Vcoc  = vector<Ofsc>(NV);

    //Read from files
    updateCOC(qbcp);
}

/**
 *  \brief Update a complex Fourier series given as the component (k,p) of a matrix matrix, from a given txt file.
 *  \param xFFT: the \c ofs<cdouble> to update
 *  \param matrix: the beginning of the name of the source txt file. e.g. "alpha"
 *  \param k the lign indix of the desired component in the matrix \c matrix.
 *  \param p the column indix of the desired component in the matrix \c matrix.
 *
 *  As an example, the call readCOC(xFFT, "alpha", 2, 1), will update xFFT with the file "alpha21.txt".
 **/
void readCOC(Ofsc& xFFT, string matrix, int k, int p)
{
    //Init
    ifstream readStream;
    string ss1, ss2;
    ss1 = static_cast<ostringstream*>( &(ostringstream() << k) )->str();
    ss2 = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    //Reading an OFS from a text file
    readOFS_txt(xFFT, (matrix+ss1+ss2));
}

/**
 *  \brief Lighter version of the init routine, with only P, PC and V retrieved.
 **/
void initCOC(matrix<Ofsc> &P,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &V,
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
    }
    else //EM RTPB
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
    }

    //----------------------------
    // Building PC
    //----------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize PC and CQ
    vector<int> keyMap(4);
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

}

/**
 *   \brief Number to string inner routine
 **/
string numTostring(double num)
{
    string res =  static_cast<ostringstream*>( &(ostringstream() << num) )->str();
    return res;
}
