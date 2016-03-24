#include "pmcoc.h"

/**
 * \file pmcoc.cpp
 * \brief Change of coordinates for the parameterization methods.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  Two types of change of coordinates are implemented:
 *              1. Changes between different manifold coordinates (Real Center to Complex Center...).
 *              2. Evaluations of the parameterization (Real Center to Normalized-Centered and projection the other way around).
 *
 */

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void CCMtoRCM(const cdouble s1[], double si[], int nv)
{
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[2]*I));
    si[2] = creal(1.0/sqrt(2)*(s1[0]*I + s1[2]));
    si[1] = creal(1.0/sqrt(2)*(s1[1]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[1]*I + s1[3]));
    if(nv == 5)
    {
        si[4] = creal(s1[4]);
    }
    else if(nv == 6)
    {
        si[4] = creal(s1[4]);
        si[5] = creal(s1[5]);
    }
}

/**
 *  \brief from RCM to CCM coordinates
 **/
void RCMtoCCM(const double si[], cdouble s1[], int nv)
{
    //From real to complex TFC
    s1[0] = 1.0/sqrt(2)*(si[0] - si[2]*I);
    s1[2] = 1.0/sqrt(2)*(si[2] - si[0]*I);
    s1[1] = 1.0/sqrt(2)*(si[1] - si[3]*I);
    s1[3] = 1.0/sqrt(2)*(si[3] - si[1]*I);
    if(nv == 5)
    {
        s1[4] = si[4]+I*0.0;
    }
    else if(nv == 6)
    {
        s1[4] = si[4]+I*0.0;
        s1[5] = si[5]+I*0.0;
    }
}

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void RCMtoCCM8(const double si[], double s0d[])
{
    //From real to complex TFC
    cdouble s1[REDUCED_NV];
    RCMtoCCM(si, s1, REDUCED_NV);

    //Store real and imag part separately
    CCMtoCCM8(s1, s0d);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void CCM8toRCM(const double s0d[], double si[])
{
    //CCM8 to CCM
    cdouble s1[REDUCED_NV];
    CCM8toCCM(s0d, s1);
    //CCM to RCM
    CCMtoRCM(s1, si, REDUCED_NV);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void CCM8toCCM(const double s0d[], cdouble s1[])
{
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        s1[p]  =   s0d[p2++];
        s1[p] += I*s0d[p2++];
    }
}

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void CCMtoCCM8(const cdouble s1[], double s0d[])
{
    //Store real and imag part separately
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        s0d[p2++] = creal(s1[p]);
        s0d[p2++] = cimag(s1[p]);
    }


//    s0d[0] = creal(s1[0]);
//    s0d[1] = cimag(s1[0]);
//
//    s0d[2] = creal(s1[1]);
//    s0d[3] = cimag(s1[1]);
//
//    s0d[4] = creal(s1[2]);
//    s0d[5] = cimag(s1[2]);
//
//    s0d[6] = creal(s1[3]);
//    s0d[7] = cimag(s1[3]);
}

/**
 *  \brief from TFC to TF coordinates
 **/
void TFCtoTF(const cdouble s1[6], double si[6])
{
    //First center
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[0]*I + s1[3]));
    //Second center
    si[2] = creal(1.0/sqrt(2)*(s1[2]   + s1[5]*I));
    si[5] = creal(1.0/sqrt(2)*(s1[2]*I + s1[5]));
    //Hyperbolic dir
    si[1] = creal(s1[1]);
    si[4] = creal(s1[4]);
}

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluation of the pm
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void RCMtoNC(const double st0[],
             const double t,
             const double n,
             const int order,
             const int ofs_order,
             vector<Oftsc> &W,
             Ofsc &ofs,
             double z1[],
             bool isGS)
{
    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    cdouble s0[REDUCED_NV];
    cdouble z0[6];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, REDUCED_NV);
    //------------------------------------------
    // 2. Update z0
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            if(p == 2 || p == 5)
            {
                //order 1 is sufficient for p = 2,5
                W[p].evaluate(s0, ofs, 1, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
            else
            {
                //For p = 0,1,3,4 normal computation
                W[p].evaluate(s0, ofs, order, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
        }
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            W[p].evaluate(s0, ofs, order, ofs_order);
            z0[p] = ofs.evaluate(n*t, ofs_order);
            z1[p] = creal(z0[p]);
        }
    }

}

/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(s0), t)
 *   \param s0 an array of 4 complex which gives the configuration to input in complex CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoTFC(cdouble s0[],
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        cdouble temp;
        // zIn[0]
        //---------------
        temp = Wh[0].getCoef(1,0)->ofs_getCoef(0);
        zIn[0].setCoef(temp*s0[0],0);
        // zIn[1]
        //---------------
        Wh[1].evaluate(s0, zIn[1], order, ofs_order);
        // zIn[2]
        //---------------
        temp = Wh[2].getCoef(1,1)->ofs_getCoef(0);
        zIn[2].setCoef(temp*s0[1],0);
        // zIn[3]
        //---------------
        temp = Wh[3].getCoef(1,2)->ofs_getCoef(0);
        zIn[3].setCoef(temp*s0[2],0);
        // zIn[4]
        //---------------
        Wh[4].evaluate(s0, zIn[4], order, ofs_order);
        // zIn[5]
        //---------------
        temp = Wh[5].getCoef(1,3)->ofs_getCoef(0);
        zIn[5].setCoef(temp*s0[3],0);
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(s0, zIn[p], order, ofs_order);
        }
    }
}


/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoTFC(const double st0[],
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[REDUCED_NV];

    //------------------------------------------
    // 1. RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, REDUCED_NV);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);
}



/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 8 double which gives the configuration to input in complex CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCM8toTFC(const double st0[],
               const int order,
               const int ofs_order,
               vector<Oftsc> &Wh,
               vector<Ofsc> &zIn,
               bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[REDUCED_NV];

    //------------------------------------------
    // 1. CCM8 to CCM
    //------------------------------------------
    CCM8toCCM(st0, s0);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);
}

/**
 *   \brief Evaluate the configuration zh1 = Wh(g(st0), n*t) in TF coordinates (not complex!)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update in TF coordinates
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void RCMtoTF(const double st0[],
              const double t,
              const double n,
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              Ofsc &ofs,
              double z1[],
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble z0[6];
    //------------------------------------------
    // Update z0
    //------------------------------------------
    RCMtoTFC(st0, t, n, order, ofs_order, Wh, ofs, z0, isGS);
    //------------------------------------------
    // TFC to TF
    //------------------------------------------
    TFCtoTF(z0, z1);
}



/**
 *   \brief Evaluate the configuration zh1 = Wh(g(st0), n*t)
 *   \param s0 an array of 4 cdouble which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z0 the output array to update, in TFC coordinates
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void CCMtoTFC(cdouble s0[],
              const double t,
              const double n,
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              Ofsc &ofs,
              cdouble z0[],
              bool isGS)
{
    //------------------------------------------
    // 1. Update zIn
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            if(p == 1 || p == 4)
            {
                //For p = 1,4 normal computation
                Wh[p].evaluate(s0, ofs, order, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
            }
            else
            {
                //order 1 is sufficient for p = 0,2,3,5
                Wh[p].evaluate(s0, ofs, 1, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
            }
        }
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(s0, ofs, order, ofs_order);
            z0[p] = ofs.evaluate(n*t, ofs_order);
        }
    }
}


/**
 *   \brief Evaluate the configuration zh1 = Wh(g(st0), n*t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z0 the output array to update, in TFC coordinates
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void RCMtoTFC(const double st0[],
              const double t,
              const double n,
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              Ofsc &ofs,
              cdouble z0[],
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[REDUCED_NV];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, REDUCED_NV);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    CCMtoTFC(s0, t, n, order, ofs_order, Wh, ofs, z0, isGS);
}


/**
 *   \brief Evaluate the configuration zh1 = Wh(g(st0), n*t)
 *   \param st0 an array of 8 double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z0 the output array to update, in TFC coordinates
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void CCM8toTFC(const double st0[],
              const double t,
              const double n,
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              Ofsc &ofs,
              cdouble z0[],
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[REDUCED_NV];

    //------------------------------------------
    // CCM8 to CCM
    //------------------------------------------
    CCM8toCCM(st0, s0);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    CCMtoTFC(s0, t, n, order, ofs_order, Wh, ofs, z0, isGS);
}



/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoNCbyTFC(const double st0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  Ofsc &ofs,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // RCM to TFC
    //------------------------------------------
    RCMtoTFC(st0, order, ofs_order, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC_OFS(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}

/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param s0d an array of 8 double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCM8toNCbyTFC(const double s0d[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  Ofsc &ofs,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    cdouble s0[REDUCED_NV];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // CCM8 to TFC
    //------------------------------------------
    CCM8toCCM(s0d, s0);
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC_OFS(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}


/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param s0 an array of 4 complex double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoNCbyTFC(cdouble s0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  Ofsc &ofs,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // CCM to TFC
    //------------------------------------------
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC_OFS(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}


/**
 *   \brief Evaluate the reduced vector field (RVF) dot(s4) = f4(s4, t)
 *   \param s8 an array of 8 doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param f4 the RVF output array of 4 cdouble to update
 **/
void CCM8toRVF(const double s8[],
               const double t,
               const double n,
               const int order,
               const int ofs_order,
               vector<Oftsc> &fh,
               Ofsc &ofs,
               cdouble f4[])
{
    //CCM8 to CCM
    //----------
    cdouble s[REDUCED_NV];
    CCM8toCCM(s8, s);

    //Evaluation of fh
    //----------
    for(int p = 0; p < REDUCED_NV; p++)
    {
//        cout << "CCM8toRVF. s[" << p << "] = " << creal(s[p]) << "   +I*" << cimag(s[p]) << endl;
//        cout << "order = " << order << endl;
//        cout << "ofs_order = " << ofs_order << endl;
        fh[p].evaluate(s, ofs, order, ofs_order);
        f4[p] = ofs.evaluate(n*t, ofs_order);
    }
}

/**
 *   \brief Evaluate the reduced vector field (RVF) dot(s4) = f4(s4, t)
 *   \param s8 an array of 8 doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param f8 the RVF output array of 8 double to update, with separated real and imag parts.
 **/
void CCM8toRVF8(const double s8[],
                const double t,
                const double n,
                const int order,
                const int ofs_order,
                vector<Oftsc> &fh,
                Ofsc &ofs,
                double f8[])
{
    //CCM8 to RVF
    //----------
    cdouble sd[REDUCED_NV];
    CCM8toRVF(s8, t, n, order, ofs_order, fh, ofs, sd);

    //Separation of real and imag parts
    //----------
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        f8[p2++] = creal(sd[p]);
        f8[p2++] = cimag(sd[p]);
    }

//    f8[0] = creal(sd[0]);
//    f8[1] = cimag(sd[0]);
//
//    f8[2] = creal(sd[1]);
//    f8[3] = cimag(sd[1]);
//
//    f8[4] = creal(sd[2]);
//    f8[5] = cimag(sd[2]);
//
//    f8[6] = creal(sd[3]);
//    f8[7] = cimag(sd[3]);
}

/**
 *  \brief Test of the RCMtoNC routines. Requires initCM().
 **/
void testRCMtoNC()
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the RCMtoNC routines        " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);

    //----------------------------------------------
    //init COC
    //----------------------------------------------
    matrix<Ofsc> P(6,6);
    matrix<Ofsc> PC(6,6);
    matrix<Ofsc> Q(6,6);
    matrix<Ofsc> CQ(6,6);
    vector<Ofsc> V(6);
    initCOC(P, PC, Q, CQ, V, SEML);

    //Orders from initCM()
    int order     = CM[0].getOrder();
    int ofs_order = CM[0].getCOrder();

    //Temp OFS variable
    Ofsc ofs(CM[0].getCOrder());

    //Results stored in:
    double z1[6], z1_QBCP[6];
    //Initial RCM conditions
    double st0[REDUCED_NV];

    //----------------------------------------------
    //Computation
    //----------------------------------------------
    //Initial RCM conditions
    for(int i = 0; i < REDUCED_NV; i++) st0[i] = 1.0;

    //RCM to NC with general routine
    RCMtoNC(st0, 1.0, SEML.us.n, order, ofs_order, CM, ofs, z1, false);
    //RCM to NC for QBCP
    RCMtoNCbyTFC(st0, 1.0, SEML.us.n, order, ofs_order, CMh, ofs, PC, V, z1_QBCP, false);

    //----------------------------------------------
    //Results
    //----------------------------------------------
    cout << "General routine             QBCP            Delta" << endl;
    for(int i = 0; i < 6; i++) cout << z1[i] << "   " << z1_QBCP[i] << "   " << z1[i]-z1_QBCP[i] << endl;
    cout << "---------------------------------------------------" << endl;

}

/**
 *  \brief Projection of the current NC state on the central manifold, via CCM coordinates
 **/
void NCprojCCM(const double z[], const double t, const double n, const int ofs_order, matrix<Ofsc> &CQ, vector<Ofsc> &V, double omega1, double omega3, cdouble sc[], int nv)
{
    //-------------------
    //Wh: TFC coordinates
    //-------------------
    cdouble zh[6];

    //-------------------
    //z - V
    //-------------------
    cdouble zd[6];
    for(int p = 0; p < 6; p++) zd[p] = z[p] - V[p].evaluate(n*t, ofs_order);

    //-------------------
    //Update Wh = CQ*(z - V)
    //-------------------
    for(int k = 0; k <6; k++)
    {
        zh[k] = 0.0+0.0*I;
        for(int p = 0; p <6; p++)
        {
            zh[k] += zd[p]* CQ.getCoef(k, p).evaluate(n*t, ofs_order);
        }
    }

    //-------------------
    //Projection on the center manifold
    //-------------------
    TFCprojCCM(zh, omega1, omega3, sc, nv);
}

/**
 *  \brief Projection of the current TFC state on the central manifold, via CCM coordinates
 **/
void TFCprojCCM(const cdouble zh[], double omega1, double omega3, cdouble sc[], int nv)
 {
    //-------------------
    //Projection on the center manifold
    //-------------------
    //Wh1 = i*w1*s1 => s1 = -i/w1*Wh1
    sc[0] = -1.0*I/omega1*zh[0];

    //Wh3 = i*w3*s2 => s2 = -i/w3*Wh3
    sc[1] = -1.0*I/omega3*zh[2];

    //Wh4 = -i*w1*s3 => s3 = +i/w1*Wh4
    sc[2] = +1.0*I/omega1*zh[3];

    //Wh6 = -i*w3*s4 => s4 = +i/w3*Wh6
    sc[3] = +1.0*I/omega3*zh[5];

    if(nv > 4) sc[4] = 0.0;
    if(nv > 5) sc[5] = 0.0;
 }

/**
 *  \brief Test for the projection on the center manifold (both in NC and TFC coordinates)
 **/
void pmProjTest(double si[4])
{
    double z[6], si2[4];
    cdouble sc[4], zh[6];
    Ofs<cdouble> AUX(OFS_ORDER);

    //----------------------------------------
    // First test: NC projected in RCM
    //----------------------------------------
    //RCM to NC
    RCMtoNCbyTFC(si, 1.0, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, CMh, AUX, Mcoc, Vcoc, z, 1);

    //NC proj in CCM
    double omega1 = cimag(Fh[0].getCoef(1,0)->ofs_getCoef(0));
    double omega3 = cimag(Fh[1].getCoef(1,1)->ofs_getCoef(0));
    NCprojCCM(z, 1.0, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, sc, 4);

    //CCM to RCM
    CCMtoRCM(sc, si2, REDUCED_NV);

    cout << "----------------------------------------" << endl;
    cout << "First test: NC projected in RCM         " << endl;
    cout << "----------------------------------------" << endl;
    cout << "Initial si:      After projection:     difference:" << endl;
    for(int i = 0; i < 4; i++) cout << si[i] << "      "  << si2[i] << "      "  << si[i] - si2[i] << endl;


    //----------------------------------------
    // First test: NC projected in RCM
    //----------------------------------------
    //RCM to TFC
    RCMtoTFC(si, 1.0, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, CMh, AUX, zh, 1);
    //Projection: TFC to CCM
    TFCprojCCM(zh, omega1, omega3, sc, 4);
    //CCM to RCM
    CCMtoRCM(sc, si2, REDUCED_NV);

    cout << "----------------------------------------" << endl;
    cout << "Second test: TFC projected in RCM         " << endl;
    cout << "----------------------------------------" << endl;
    cout << "Initial si:      After projection:     difference:" << endl;
    for(int i = 0; i < 4; i++) cout << si[i] << "      "  << si2[i] << "      "  << si[i] - si2[i] << endl;
}
