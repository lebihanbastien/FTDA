/**
 * \file define_env.cpp
 * \brief Define the working environment (the Sun, planets and the moon). All values taken from JPL and Goddard Space Flight Center websites.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

#include "define_env.h"
#include "Constants.h"

//----------------------------------------------------------------------------------------
//            Init routines
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize a QBCP_L structure, i.e. a QBCP focused on two libration points: one the EM system and one of the SEM system.
 *        The libration point must be L1 or L2 both for the EM and SEM systems.
 * \param qbcp_l a pointer on the QBCP_L structure to init.
 * \param qbcp a pointer on the QBCP parent structure.
 * \param isNormalized: are the equations of motion normalized? Probably deprecated, should be always true. Kept for consistency with older code.
 * \param li_EM number of the libration point considered for the EM system
 * \param li_SEM number of the libration point considered for the SEM system
 * \param isNew an integer: equal to 1 if no solution has been previously computed with the routine qbtbp(), 0 otherwise
 * \param model: QBCP, BCP, CRTBP...
 * \param coordsys: default coordinate system for this structure: for example:
 *                      - if coordsys == Csts::EM,  the qbcp_l is focused on the li_EM  point of the EM  system.
 *                      - if coordsys == Csts::SEM, the qbcp_l is focused on the li_SEM point of the SEM system.
 *        The default focus can be change dynamically during computation, via the routines changeCOORDSYS and changeLICOORDSYS.
 * \param pmType: type of parameterization of the manifolds (Csts::GRAPH, Csts::NORMFORM...). Note that the pmType influences the number of coefficients taken into account in the Fourier series! Indeed, for graph method, the reduced vector field is non autonomous, and full Fourier series are used. For normal form, the reduced vector field is quasi autonomous and we can safely reduce the order of the series to 5 (11 coefficients taken into account).
 * \param manType_EM: type of manifold about li_EM (Csts::MAN_CENTER, Csts::MAN_CENTER_S...).
 * \param manType_SEM: type of manifold about li_SEM (Csts::MAN_CENTER, Csts::MAN_CENTER_S...).
 *
 *   Note that the QBCP structure is used only for the initialization of the coordinate systems. More precisely, it contains some parameters
 *   specific to each libration point (gamma), via its CR3BP structures (see init_QBCP, the routine that initializes the QBCP structures).
 **/
void init_QBCP_L(QBCP_L *qbcp_l, QBCP *qbcp, int isNormalized, int li_EM, int li_SEM,
                 int isNew, int model, int coordsys, int pmType, int manType_EM, int manType_SEM)
{
    //====================================================================================
    //      Common to all models
    //      These settings may be needed to initialize the CSYS and USYS structures
    //      For this reason, they are initialized in priority
    //====================================================================================

    //------------------------------------------------------------------------------------
    // - Hard-coded value of the order of the Fourier series. Equals zero for RTBP
    // - Note that for the bicircular models, the order of the Fourier series is still
    // equals to OFS_ORDER.
    // This may need to be changed when some dynamic order will be used, to be able
    // to manipulate Fourier series with an order smaller than the value 30, hard coded
    // in the data files.
    // - Used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    switch(model)
    {
    case Csts::QBCP:
    case Csts::BCP:
    case Csts::ERTBP:
        qbcp_l->nf = OFS_ORDER;
        break;
    case Csts::CRTBP:
        qbcp_l->nf = 0;
        break;
    default:
        cout << "init_QBCP_L. Warning: unknown model." << endl;
    }

    //------------------------------------------------------------------------------------
    //Normalization or not. True is always preferable - may be deprecated
    //Kept for compatibility with old code.
    //------------------------------------------------------------------------------------
    qbcp_l->isNormalized = isNormalized;

    //------------------------------------------------------------------------------------
    //Libration point
    //Note: passing the complete libration point structure as argument may be done in the near future.
    //not necessary for now.
    //------------------------------------------------------------------------------------
    qbcp_l->li_EM  = li_EM;
    qbcp_l->li_SEM = li_SEM;

    //------------------------------------------------------------------------------------
    // Model - used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    qbcp_l->model   = model;

    //------------------------------------------------------------------------------------
    // Parameterization style - used in the CSYS init functions.
    //------------------------------------------------------------------------------------
    qbcp_l->pms   = pmType;

    //------------------------------------------------------------------------------------
    // Effective order of the Fourier series in RVF (reduced vector field)
    //------------------------------------------------------------------------------------
    switch(pmType)
    {
        case Csts::GRAPH:
        case Csts::MIXED:
            qbcp_l->eff_nf = qbcp_l->nf;    //for graph method, the reduced vector field is non autonomous
        break;
        case Csts::NORMFORM:
            qbcp_l->eff_nf = min(qbcp_l->nf, 5);  //for normal form, the reduced vector field is quasi autonomous and we can safely reduce the value to 5
        break;
        default:
        cout << "init_QBCP_L. Warning: unknown pmType." << endl;
    }

    //====================================================================================
    // Unit systems
    //====================================================================================
    init_USYS(&qbcp_l->us_em,  Csts::EM,  model);
    init_USYS(&qbcp_l->us_sem, Csts::SEM, model);

    //====================================================================================
    // Coordinate systems
    //====================================================================================
    qbcp_l->numberOfCoefs = 15; // Number of coefficients in the equations of motion

    //------------------------------------------------------------------------------------
    // Earth-Moon
    //------------------------------------------------------------------------------------
    init_CSYS(&qbcp_l->cs_em_l1, qbcp_l, qbcp,  Csts::EM, 1, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L1
    init_CSYS(&qbcp_l->cs_em_l2, qbcp_l, qbcp,  Csts::EM, 2, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L2
    init_CSYS(&qbcp_l->cs_em_l3, qbcp_l, qbcp,  Csts::EM, 3, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L3

    //------------------------------------------------------------------------------------
    // Sun-Earth
    //------------------------------------------------------------------------------------
    init_CSYS(&qbcp_l->cs_sem_l1, qbcp_l, qbcp,  Csts::SEM, 1, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L1
    init_CSYS(&qbcp_l->cs_sem_l2, qbcp_l, qbcp,  Csts::SEM, 2, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L2
    init_CSYS(&qbcp_l->cs_sem_l3, qbcp_l, qbcp,  Csts::SEM, 3, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L3


    //====================================================================================
    // DEFAULT SETTINGS
    // - The flag li_EM determines the default Earth-Moon libration point.
    // - The flag li_SEM determines the default Sun-(Earth+Moon) libration point.
    // - The flag coordsys determines the default focus; either on the Earth-Moon or
    //   Sun-Earth+Moon framework
    //====================================================================================
    //Coord. syst. for the EM point
    //--------------------
    switch(li_EM)
    {
    case 1:
        qbcp_l->cs_em  = qbcp_l->cs_em_l1;
        break;
    case 2:
        qbcp_l->cs_em  = qbcp_l->cs_em_l2;
        break;
    case 3:
        qbcp_l->cs_em  = qbcp_l->cs_em_l3;
        break;
    }

    //--------------------
    //Coord. syst. for the SE point
    //--------------------
    switch(li_SEM)
    {
    case 1:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l1;
        break;
    case 2:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l2;
        break;
    case 3:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l3;
        break;
    }

    //--------------------
    //Default Li, Unit & Coord. system
    //--------------------
    qbcp_l->coordsys = coordsys;
    switch(coordsys)
    {
    case Csts::EM:
        qbcp_l->us = qbcp_l->us_em;
        qbcp_l->cs = qbcp_l->cs_em;
        qbcp_l->li = qbcp_l->li_EM;
        break;
    case Csts::SEM:
        qbcp_l->us = qbcp_l->us_sem;
        qbcp_l->cs = qbcp_l->cs_sem;
        qbcp_l->li = qbcp_l->li_SEM;
        break;
    default:
        cout << "init_QBCP_L. Warning: unknown coordsys." << endl;
    }

    //------------------------------------------------------------------------------------
    // 6*6 matrix B used for Floquet transformation
    // Do we need it here?
    //------------------------------------------------------------------------------------
    qbcp_l->B  = (double*) calloc(36, sizeof(double));

    //------------------------------------------------------------------------------------
    // When epsilon = 1.0, the Moon is "on". When epsilon = 0.0, the Moon is "off"
    //------------------------------------------------------------------------------------
    //Moon "on" by default. Deprecated value for now.
    qbcp_l->epsilon = 1.0;
}

/**
 * \brief Initializes the Circular Restricted 3-Body Problem
 * \param cr3bp pointer on the CR3BP
 * \param n1 name of the first primary
 * \param n2 name of the second primary
 **/
void init_CR3BP(CR3BP *cr3bp, int n1, int n2)
{
    //Body initialization
    init_body(&(*cr3bp).m1, n1);
    init_body(&(*cr3bp).m2, n2);

    cr3bp->mu = (*cr3bp).m2.M/( (*cr3bp).m1.M + (*cr3bp).m2.M );  // µ = m2/(m1 + m2)
    cr3bp->L  = (*cr3bp).m2.a;                                    // Distance parameter = semi major axis of m2
    cr3bp->T  = (*cr3bp).m2.T;                                    //Time parameter = sidereal period of m2
    cr3bp->R1 = (*cr3bp).m1.Req;
    cr3bp->R2 = (*cr3bp).m2.Req;
    cr3bp->rh = pow((*cr3bp).mu/3,1/3.0);                         //Hill's radius adim formula

    //Li initialization
    init_libp(&cr3bp->l1, *cr3bp, 1);
    init_libp(&cr3bp->l2, *cr3bp, 2);
    init_libp(&cr3bp->l3, *cr3bp, 3);
    init_libp(&cr3bp->l4, *cr3bp, 4);
    init_libp(&cr3bp->l5, *cr3bp, 5);

    //Name
    strcpy(cr3bp->name, cr3bp->m1.name);
    strcat(cr3bp->name, "-");
    strcat(cr3bp->name, cr3bp->m2.name);

    //Distance to manifold approximation
    if(n1 == Csts::EARTH && n2 == Csts::MOON) cr3bp->d_man = 50/cr3bp->L; //50km for the Earth-Moon system
    else cr3bp->d_man = 100/cr3bp->L; //100km otherwise (to be changed for specific systems!)
}

/**
 * \brief Initializes a unit system in the form of a usys structure, such as Earth-Moon or Sun-(Earth+Moon) unit systems.
 * \param usys pointer on the usys structure to init.
 * \param label the type of unit system
 *
 * NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES? Not for now. Consistent with appendices of the PhD manuscript as of this day (30/08/2016).
 **/
void init_USYS(USYS *usys, int label, int model)
{
    //--------------------------
    // These values are used as reference values for all other constants
    //--------------------------
    //EM mass ratio
    usys->mu_EM  = +1.215058162343360e-02;
    //Lunar eccentricity
    usys->lecc   = +0.054900489;

    //---------------------------
    // Physical params in EM units
    // These values are also used
    // as reference values
    // for all other constants
    //--------------------------
    //Pulsation of the system
    double n_EM  = 0.925195985520347;
    //Outer pulsation
    double ns_EM = 1.0-n_EM;
    //Sun radius
    double as_EM = 388.81114;
    //Sun mass
    double ms_EM = ns_EM*ns_EM*pow(as_EM,3.0)-1;
    //Earth mass
    double me_EM = 1.0-usys->mu_EM;
    //Moon mass
    double mm_EM = usys->mu_EM;

    //Rest of the parameters
    switch(model)
    {
    case Csts::QBCP:
    case  Csts::BCP:
    {
        switch(label)
        {
        case Csts::EM :
        {
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
            break;
        }

        case Csts::SEM :
        {
            //Physical params in SEM units
            //--------------------------
            usys->ns = 1.0;
            usys->ni = 1.0/ns_EM;
            usys->n  = usys->ni-usys->ns;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = ms_EM/(1.0+ms_EM);
            usys->me = me_EM/(1.0+ms_EM);
            usys->mm = mm_EM/(1.0+ms_EM);
            break;
        }

        default :  //Earth-Moon
        {
            cout << "init_USYS. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
        }
        }

        //SEM mass ratio
        usys->mu_SEM = (usys->me+usys->mm)/(usys->ms+usys->me+usys->mm);
        //SE mass ratio
        usys->mu_SE  = (usys->me)/(usys->ms+usys->me);
        break;

    }
    case Csts::CRTBP:
    case Csts::ERTBP:
    {
        //SEM mass ratio
        usys->mu_SEM = (me_EM + mm_EM)/(ms_EM + me_EM + mm_EM);
        //SE mass ratio
        usys->mu_SE  = usys->mu_EM*usys->mu_SEM/(1.0 + usys->mu_EM);

        switch(label)
        {
        case Csts::EM:

            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = 0.0;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
            break;

        case Csts::SEM:

            //Physical params in SEM units
            //--------------------------
            usys->ns = 1.0;
            usys->ni = 1.0/ns_EM;
            usys->n  = usys->ni-usys->ns;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = 1.0-usys->mu_SEM;
            usys->me = usys->mu_SEM; //the Earth contains both the Earth's and the Moon's masses.
            usys->mm = 0.0;//mm_EM/(1.0+ms_EM);
            break;

        default:  //Earth-Moon
            cout << "init_USYS. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->n  = n_EM;
            usys->ns = ns_EM;
            usys->ni = 1.0;
            usys->as = as_EM;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;

            usys->ms = 0.0;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
        }
    }
    }


}

/**
 * \brief Initializes a coordinate systems (CSYS structure), with associated vector field coefficients, data folder names, and unit system.
 * \param csys pointer on the CSYS structure to initialize.
 * \param qbcp_l pointer on the QBCP_L structure that contains csys.
 * \param qbcp pointer on the QBCP structure that contains parameters specific to each libration points (namely, gamma)
 * \param coordsys indix of the coordinate system to use (Csts::EM, Csts::SEM).
 * \param li number of the libration point to focus on (L1, L2).
 * \param coefNumber the number of vector field coefficients to initialize. It has been set in the QBCP_init function.
 * \param isNew boolean. if true, the qbtbp has not been computed via the qbtbp() routine, so the vector field coefficients cannot be initialized.
 * \param pmType: type of parameterization of the manifolds (Csts::GRAPH, Csts::NORMFORM...).
 * \param manType: type of manifold about li (Csts::MAN_CENTER, Csts::MAN_CENTER_S...).
 *
 *   Note that the QBCP structure is used only for the initialization of the coordinate systems. More precisely, it contains some parameters
 *   specific to each libration point (gamma), via its CR3BP structures.
 **/
void init_CSYS(CSYS *csys, QBCP_L *qbcp_l, QBCP *qbcp, int coordsys, int li, int coefNumber, int isNew, int pmType, int manType)
{
    int nf = qbcp_l->nf;
    int model = qbcp_l->model;

    //--------------------------------------
    //Complete folders
    //--------------------------------------
    csys->F_GS     = init_F_FOLDER("data/PMFFT",  model, coordsys, li);     //Graph style (PM)
    csys->F_NF     = init_F_FOLDER("data/NF",     model, coordsys, li);     //Normal form style(PM)
    csys->F_MS     = init_F_FOLDER("data/MS",     model, coordsys, li);     //Mixed style (PM)

    csys->F_CS     = init_F_FOLDER("data/CS",     model, coordsys, li);     //Center-stable (PM)
    csys->F_CU     = init_F_FOLDER("data/CU",     model, coordsys, li);     //Center-unstable (PM)
    csys->F_CUS    = init_F_FOLDER("data/CUS",    model, coordsys, li);     //Center-hyperbolic (PM)

    csys->F_COEF   = init_F_FOLDER("data/VF",     model, coordsys, li);     //For integration in a given coord. system
    csys->F_COC    = init_F_FOLDER("data/COC",    model, coordsys, li);     //For the change of coordinates of the PM
    csys->F_PLOT   = init_F_FOLDER("plot",        model, coordsys, li);     //For plot output (gnuplot)
    csys->F_PRINT  = init_F_FOLDER("fprint",      model, coordsys, li);     //For print output (data used postprocessed in R)

    //--------------------------------------
    //Parameterization folders
    //--------------------------------------
    csys->manType = manType;
    switch(manType)
    {
    //If the Center Manifold is selected,
    //we can choose between the styles
    case Csts::MAN_CENTER:

        switch(pmType)
        {
        case Csts::GRAPH:
            csys->F_PMS = csys->F_GS;
            break;
        case Csts::NORMFORM:
            csys->F_PMS = csys->F_NF;
            break;
        case Csts::MIXED:
            csys->F_PMS = csys->F_MS;
            break;

        }
        break;
    //Else, the style is fixed
    case Csts::MAN_CENTER_S:
        csys->F_PMS = csys->F_CS;
        break;

    case Csts::MAN_CENTER_U:
        csys->F_PMS = csys->F_CU;
        break;

    case Csts::MAN_CENTER_US:
        csys->F_PMS = csys->F_CUS;
        break;
    }

    //--------------------------------------
    // Unit system associated with csys
    //--------------------------------------
    switch(coordsys)
    {
    case Csts::EM:
        csys->us = qbcp_l->us_em;
        break;
    case Csts::SEM:
        csys->us = qbcp_l->us_sem;
        break;
    default:
        cout << "init_CSYS. Warning: unknown model." << endl;
    }

    //--------------------------------------
    // c1 and gamma, for normalized computation
    //--------------------------------------
    // First, select the right framework
    // Gives the assosicate CR3BP and mu
    //------------------------------
    CR3BP cr3bp_root;
    switch(coordsys)
    {
    case Csts::EM:
        cr3bp_root  = qbcp->cr3bp1;
        csys->cr3bp = qbcp->cr3bp1;
        csys->mu    = qbcp_l->us_em.mu_EM;
        break;
    case Csts::SEM:
        cr3bp_root  = qbcp->cr3bp2;
        csys->cr3bp = qbcp->cr3bp2;
        csys->mu    = qbcp_l->us_sem.mu_SEM;
        break;
    default:
        cout << "init_CSYS. Warning: unknown model." << endl;
        cr3bp_root  = qbcp->cr3bp1; //cr3bp_root is set here to avoid compilation warning, but the default case should NOT be used!
        csys->cr3bp = qbcp->cr3bp1;
        csys->mu    = qbcp_l->us_em.mu_EM;
    }

    //------------------------------
    //Then set the value of c1 and gamma
    //according to the selected libration point
    //------------------------------
    csys->li = li;
    switch(li)
    {
    case 1:
    {
        csys->gamma = cr3bp_root.l1.gamma_i;
        csys->c1    = (csys->mu-1+csys->gamma)/csys->gamma;
        break;
    }

    case 2:
    {
        csys->gamma =  cr3bp_root.l2.gamma_i;
        csys->c1    = (csys->mu-1-csys->gamma)/csys->gamma;
        break;
    }

    case 3: //same as L1/L2. Different from convention by Jorba & Masdemont 1999
    {
        csys->gamma =  cr3bp_root.l3.gamma_i;
        csys->c1    = (csys->mu + csys->gamma)/csys->gamma;
        break;
    }
    }

    //------------------------------
    //c2 coefficient
    //------------------------------
    csys->c2 = cn(csys->li, csys->gamma, csys->mu, 2);

    //--------------------------------------
    // Creating the arrays of coefficients
    //--------------------------------------
    csys->coeffs = (double*) calloc(coefNumber*(qbcp_l->nf+1), sizeof(double)); //Default set of vector field coefficients
    csys->Ps  = (double*) calloc(3*(nf+1), sizeof(double));  //Sun   position in EM coordinates
    csys->Pe  = (double*) calloc(3*(nf+1), sizeof(double));  //Earth position in EM coordinates
    csys->Pm  = (double*) calloc(3*(nf+1), sizeof(double));  //Moon  position in EM coordinates
    csys->ps  = (double*) calloc(3*(nf+1), sizeof(double));  //Sun   position in NC coordinates
    csys->pe  = (double*) calloc(3*(nf+1), sizeof(double));  //Earth position in NC coordinates
    csys->pm  = (double*) calloc(3*(nf+1), sizeof(double));  //Moon  position in NC coordinates

    //--------------------------------------
    // Retrieving the arrays of coefficients
    //--------------------------------------
    if(isNew) //the QBCP is new, no coefficients exist
    {
        cout << "init_CSYS. The considered QBCP is new:" << endl;
        cout << "Integration coefficient are not retrieved from existing files." << endl;
    }
    else    //the coefficients can be retrieved
    {
        //cout << "init_CSYS. Integration coefficient are retrieved from existing files." << endl;
        //The flag compType is here to tell if we want the coefficients
        //computed from the FFT or from direct computation
        double compType;

        //Switch between the models
        switch(model)
        {
        case Csts::QBCP:
        case Csts::BCP:
        case Csts::ERTBP:
        {
            compType = 1; //from FFTs
            //-------------------------
            // Default set of vector field coefficients
            //-------------------------
            coefRetrieving(csys->F_COEF+"alpha", csys->coeffs, nf, 0, compType, coefNumber);
            //-------------------------
            // primaries position
            //-------------------------
            coefRetrieving(csys->F_COEF+"Ps", csys->Ps, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"Pe", csys->Pe, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"Pm", csys->Pm, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"ps", csys->ps, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"pe", csys->pe, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"pm", csys->pm, nf, 0, compType, 3);
            break;

        }

        case Csts::CRTBP:
        {
            //cout << "init_CSYS. The use of the RTBP has been detected." << endl;
            //----------------------------------------
            // Default set of vector field coefficients
            // A REVOIR: METTRE UN SWITCH GEANT POUR TOUS LES CAS
            // OU: LES METTRE DANS DES FICHIERS TXT via une routine type bcp
            //----------------------------------------
            csys->coeffs[0] = 1.0;
            csys->coeffs[1] = 0.0;
            csys->coeffs[2] = 1.0;
            csys->coeffs[3] = 0.0;
            csys->coeffs[4] = 0.0;
            csys->coeffs[5] = 1.0;

            switch(coordsys)
            {
                case Csts::EM:
                    //Sun position
                    csys->coeffs[6] = 0.0;
                    csys->coeffs[7] = 0.0;
                    //Earth position
                    csys->coeffs[8] = csys->mu;
                    csys->coeffs[9] = 0.0;
                    //Moon position
                    csys->coeffs[10] = csys->mu-1.0;
                    csys->coeffs[11] = 0.0;
                break;

                case Csts::SEM:
                    //Sun position
                    csys->coeffs[6]  = csys->mu;
                    csys->coeffs[7]  = 0.0;
                    //Earth position
                    csys->coeffs[8]  = csys->mu-1.0;
                    csys->coeffs[9]  = 0.0;
                    //Moon position
                    csys->coeffs[10] = 0.0;
                    csys->coeffs[11] = 0.0;
                break;
            }
            //NC coeffs
            csys->coeffs[12] = -csys->c1;     //alpha13 = alpha[12] = -c1
            csys->coeffs[13] = 0.0;           //alpha14 = alpha[13] = 0


            //----------------------------------------
            // Earth, Moon, and Sun €€TODO: a revoir!
            // Mettre dans un fichier txt...
            //----------------------------------------
            switch(coordsys)
            {
                case Csts::EM:
                    csys->Pe[0] = csys->mu;
                    csys->Pe[1] = 0.0;
                    csys->Pe[2] = 0.0;

                    csys->Pm[0] = csys->mu-1.0;
                    csys->Pm[1] = 0.0;
                    csys->Pm[2] = 0.0;

                    csys->Ps[0] = 0.0;
                    csys->Ps[1] = 0.0;
                    csys->Ps[2] = 0.0;
                    break;

                case Csts::SEM:
                    csys->Pe[0] = csys->mu-1.0;
                    csys->Pe[1] = 0.0;
                    csys->Pe[2] = 0.0;

                    csys->Pm[0] = 0.0;
                    csys->Pm[1] = 0.0;
                    csys->Pm[2] = 0.0;

                    csys->Ps[0] = csys->mu;
                    csys->Ps[1] = 0.0;
                    csys->Ps[2] = 0.0;
                    break;
            }

            SYStoNC_prim(csys->Pe, csys->pe, csys->c1, csys->gamma);
            SYStoNC_prim(csys->Pm, csys->pm, csys->c1, csys->gamma);
            SYStoNC_prim(csys->Ps, csys->ps, csys->c1, csys->gamma);

            break;
        }

        default:
            cout << "init_QBCP. Warning: unknown model." << endl;
        }

    }


    //--------------------------------------
    // Storing the solutions of the QBTBP
    //--------------------------------------
    csys->zt = Ofsc(nf);
    csys->Zt = Ofsc(nf);
    csys->ztdot = Ofsc(nf);
    csys->Ztdot = Ofsc(nf);
    switch(model)
    {
    case Csts::CRTBP:
    case Csts::BCP:
    {
        //Perfect circles
        csys->zt.setCoef(1.0+0.0*I, 0);
        csys->Zt.setCoef(1.0+0.0*I, 0);
        break;
    }

    case Csts::QBCP:
    default:
    {
        //Taken from files
        string filename = "data/qbtbp/";
        readOFS_txt(csys->zt, filename+"bjc");
        readOFS_txt(csys->Zt, filename+"cjc");
        //Derivatives
        csys->ztdot.dot(csys->zt, csys->us.n);
        csys->Ztdot.dot(csys->Zt, csys->us.n);
        break;
    }
    }
}

/**
* \brief Initialize the Quasi-Bicircular Four-Body Problem in the form of a QBCP structure.
* \param qbcp pointer on the QBCP structure to init.
* \param BIG name of the first primary (e.g Sun)
* \param MEDIUM name of the second primary (e.g. Earth)
* \param SMALL name of the third primary  (e.g. Moon)
*
* At the end of this routine:
*    - qbcp->cr3bp1 is initialized as the CRTBP of (MEDIUM, SMALL).
*    - qbcp->cr3bp2 is initialized as the CRTBP of (BIG, MEDIUM).
*
* WARNING: in the case of the Sun-Earth-Moon system (in fact, the only interesting case for us...), we need to set (SUN, EARTH_AND_MOON), instead of (SUN, EARTH), in
* qbcp->cr3bp2. So keep in mind that if the configuration is SUN/EARTH/MOON, the result will NOT be (BIG, MEDIUM) in qbcp->cr3bp2.
*
**/
void init_QBCP(QBCP *qbcp, int BIG, int MEDIUM, int SMALL)
{
    //E.g. Earth-Moon init
    init_CR3BP(&qbcp->cr3bp1, MEDIUM, SMALL);
    //E.g. Sun-Earth+Moon init
    //WARNING: in the case of the Sun-Earth-Moon system, we need to set EARTH_AND_MOON, instead of Csts::EARTH alone.
    if(MEDIUM ==  Csts::EARTH && SMALL == Csts::MOON) init_CR3BP(&qbcp->cr3bp2, BIG, Csts::EARTH_AND_MOON);
    else init_CR3BP(&qbcp->cr3bp2, BIG, MEDIUM);
}

/**
 *  \brief Initializes two QBCP_L with two different models for continuation process from one model (epsilon = 0.0) to the other (epsilon = 1.0).
 **/
void init_QBCP_I(QBCP_I *model,
                 QBCP_L *model1,
                 QBCP_L *model2,
                 int n1, int n2, int n3,
                 int isNormalized, int li_EM, int li_SEM,
                 int isNew,
                 int mod1,
                 int mod2,
                 int frwk,
                 int pmStyle)
{
    //Initialize the models
    QBCP qbp1, qbp2;
    init_QBCP(&qbp1, n1, n2, n3);
    init_QBCP(&qbp2, n1, n2, n3);

    //Initialize the models around the given li point
    init_QBCP_L(model1, &qbp1, isNormalized, li_EM, li_SEM, isNew, mod1, frwk, pmStyle, Csts::MAN_CENTER, Csts::MAN_CENTER);
    init_QBCP_L(model2, &qbp2, isNormalized, li_EM, li_SEM, isNew, mod2, frwk, pmStyle, Csts::MAN_CENTER, Csts::MAN_CENTER);

    //Store in model
    model->model1 = *model1;
    model->model2 = *model2;
    model->epsilon = 0.0;
}


/**
 * \brief Initializes a libration point.
 * \param libp a pointer towards the LibrationPoint structure to init.
 * \param cr3bp a CR3BP structure that contains useful coefficients.
 * \param number the indix of the libration point to init.
 **/
void init_libp(LibrationPoint *libp, CR3BP cr3bp, int number)
{

    double gamma_i;

    //Number
    libp->number = number;

    switch(number)
    {
    case 1:
        //Gamma
        gamma_i = cr3bp.rh - 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);                    //initial guess
        gamma_i = rtnewt(polynomialLi, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);   //newton-raphson method
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu - gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l1 = (%5.10f, 0, 0)\n", libp->position[0]);
        break;

    case 2:
        //Gamma
        gamma_i = cr3bp.rh + 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);
        gamma_i = rtnewt(polynomialLi, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu + gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l2 = (%5.10f, 0, 0)\n", libp->position[0]);
        break;

    case 3:
        //Gamma
        gamma_i = 7/12.0*cr3bp.mu + pow(237,2.0)/pow(12,4.0)*pow(cr3bp.mu,3.0);
        gamma_i = rtnewt(polynomialLi, gamma_i, Config::configManager().G_PREC_LIB(), cr3bp.mu, number);
        libp->gamma_i = 1-gamma_i;  //BEWARE: for L3, gamma3 = L3-M1 distance != L3-M2


        //Position
        libp->position[0] = - cr3bp.mu - libp->gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l3 = (%5.10f, 0, 0)\n", libp->position[0]);
        break;

    case 4:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = sqrt(3)/2.0;
        libp->position[2] = 0;
        break;

    case 5:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = -sqrt(3)/2.0;
        libp->position[2] = 0;
        break;
    }

    //Energy & Jacobi constant
    libp->Ei = energy(libp->position, cr3bp.mu);
    libp->Ci = -2*libp->Ei;
}

/**
* \brief Initialize one celestial body
* \param body a pointer on the Body structure to init.
* \param name the name of the body in integer format (consistent with HORIZON numerotation)
**/
void init_body(Body *body, int name)
{

    double days = 86400; //days to seconds

    switch(name)
    {

    case Csts::MERCURY:

        //Physical parameters
        body->Req = 2439.7;        //[km]
        body->Rm = 2439.7;         //[km]
        body->M = 0.330104e24;     //[kg]
        body->GM = 22032;          //[km^3.s^-2]

        //Orbital parameters
        body->a = 57.91e6;         //[kg]
        body->T = 87.9691*days;     //[s]

        strcpy(body->name, "Mercury");
        break;

    case Csts::VENUS:

        //Physical parameters
        body->Req = 6051.8;        //[km]
        body->Rm = 6501.8;         //[km]
        body->M = 4.86732e24;      //[kg]
        body->GM = 324858.63;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 108.21e6;        //[km]
        body->T = 224.701*days;     //[s]

        strcpy(body->name, "Venus");
        break;


    case Csts::EARTH:

        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm  = 6371.00;         //[km]
        body->M   = 5.97219e24;       //[kg]
        body->GM  = 398600.440;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth");
        break;

    case Csts::MOON:

        //Physical parameters
        body->Req = 1737.5;       //[km]
        body->Rm = 1737.5;        //[km]
        body->M =  0.07345814120628661e24;    //[kg] TO BE CONSISTENT WITH HARD CODED VALUE OF mu(Earth-Moon)
        body->GM = 4902.801;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 384400;           //[km]
        body->T = 27.321582*days;    //[s]

        strcpy(body->name, "Moon");
        break;

    case Csts::MARS:

        //Physical parameters
        body->Req = 3396.19;       //[km]
        body->Rm = 3389.50;        //[km]
        body->M = 0.641693e24;     //[kg]
        body->GM = 42828.3;        //[km^3.s^-2]

        //Orbital parameters
        body->a = 227.92e6;       //[kg]
        body->T = 686.98*days;     //[s]

        strcpy(body->name, "Mars");
        break;


    case Csts::JUPITER:

        //Physical parameters
        body->Req =  71492;      //[km]
        body->Rm = 69911;        //[km]
        body->M = 1898.13e24;     //[kg]
        body->GM = 126686511;       //[km^3.s^-2]

        //Orbital parameters
        body->a = 778.57e6;       //[kg]
        body->T = 4332.589*days;     //[s]

        strcpy(body->name, "Jupiter");
        break;

    case Csts::SATURN:

        //Physical parameters
        body->Req =  60268;    //[km]
        body->Rm = 58232;       //[km]
        body->M = 568.319e24;     //[kg]
        body->GM = 37931207.8;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 1433.53e6;       //[kg]
        body->T = 10759.22*days;     //[s]

        strcpy(body->name, "Saturn");
        break;

    case Csts::URANUS:

        //Physical parameters
        body->Req = 25559;      //[km]
        body->Rm = 25362;        //[km]
        body->M =  86.8103e24;    //[kg]
        body->GM =  5793966;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  2872.46e6;      //[kg]
        body->T =  30685.4*days;   //[s]

        strcpy(body->name, "Uranus");
        break;

    case Csts::NEPTUNE:

        //Physical parameters
        body->Req = 24764;      //[km]
        body->Rm = 24622;        //[km]
        body->M = 102.410e24;     //[kg]
        body->GM =  6835107;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  4495.06e6;      //[kg]
        body->T =  60189*days;    //[s]

        strcpy(body->name, "Neptune");
        break;

    case Csts::PLUTO:

        //Physical parameters
        body->Req =  1195;     //[km]
        body->Rm =  1195;       //[km]
        body->M = .01309e24;     //[kg]
        body->GM =  872.4;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  5906.38e6;      //[kg]
        body->T =  90465*days;    //[s]

        strcpy(body->name, "Pluto");
        break;


    case Csts::SUN:

        //Physical parameters
        body->Req = 696342;                //[km]
        body->Rm =  696342;                //[km]
        body->M  = 1988500e24;             //[kg]
        body->GM = 1.3271244004193938e11;  //[km^3.s^-2]

        //Orbital parameters
        body->a = 0;    //[kg]
        body->T = 0;    //[s]

        strcpy(body->name, "Sun");
        break;

    case Csts::EARTH_AND_MOON:
        //Equivalent mass of the Earth+Moon system based at the center of mass
        //additionnal physical properties are those of the Earth for consistency)
        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm = 6371.00;         //[km]
        body->M = 6.04590064229622e+24; //[kg]  TO BE CONSISTENT WITH HARD CODED VALUE OF mu(Sun-(Earth+Moon))
        body->GM = 398600.440+4902.801;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth+Moon");
        break;


    }

}

//-----------------------------------------------------------------------------------------------------------------------------
//            COC: SYSTEM (EM or SEM) to NORMALIZED-CENTERED coordinates, for the primaries
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief From SYSTEM (EM or SEM) to NC coordinates for the primaries. Used in qbtbp_ofs_fft_*
 */
void SYStoNC_prim(double Zc[3], double zc[3], double c1, double gamma)
{
    zc[0] = c1 - Zc[0]/gamma;
    zc[1] =    - Zc[1]/gamma;
    zc[2] =    + Zc[2]/gamma;
}

//-----------------------------------------------------------------------------------------------------------------------------
//            Change Coord. System
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Change the default coordinate system of the QBCP_L structure.
 **/
void changeCOORDSYS(QBCP_L &qbcp_l, int coordsys)
{
    qbcp_l.coordsys = coordsys;
    switch(coordsys)
    {
    case Csts::EM:
        qbcp_l.us = qbcp_l.us_em;
        qbcp_l.cs = qbcp_l.cs_em;
        qbcp_l.li = qbcp_l.li_EM;
        break;
    case Csts::SEM:
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        qbcp_l.li = qbcp_l.li_SEM;
        break;
    default:
        cout << "changeCOORDSYS. Warning: unknown coordsys." << endl;
    }
}

/**
 *  \brief Change the default coordinate system and the libration point for this coordinate system, in the QBCP_L structure.
 **/
void changeLICOORDSYS(QBCP_L &qbcp_l, int coordsys, int li)
{
    //Default settings
    qbcp_l.coordsys = coordsys;  //new default cs
    qbcp_l.li = li;      //new default libration point

    //Change the coord. system approprietly: "coordsys" around "li"
    switch(coordsys)
    {
    case Csts::EM:
        switch(li)
        {
        case 1:
            qbcp_l.cs_em  = qbcp_l.cs_em_l1;
            break;
        case 2:
            qbcp_l.cs_em  = qbcp_l.cs_em_l2;
            break;
        case 3:
            qbcp_l.cs_em  = qbcp_l.cs_em_l3;
            break;
        }
        qbcp_l.us = qbcp_l.us_em;
        qbcp_l.cs = qbcp_l.cs_em;
        break;
    case Csts::SEM:
        switch(li)
        {
        case 1:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l1;
            break;
        case 2:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l2;
            break;
        case 3:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l3;
            break;
        }
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        break;
    default:
        cout << "changeLICOORDSYS. Warning: unknown coordsys." << endl;
    }
}


//-----------------------------------------------------------------------------------------------------------------------------
//            Subroutines
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Return the string corresponding to the libration point number provided (e.g. "L1" if li == 1).
 **/
string init_F_LI(int li)
{
    switch(li)
    {
    case 1:
        return "L1";
    case 2:
        return "L2";
    case 3:
        return "L3";
    case 4:
        return "L4";
    case 5:
        return "L5";
    default:
        cout << "init_F_LI. Warning: supplied libration number is greater than 5." << endl;
    }
    return "L1"; //never here
}

/**
 *  \brief Return the string corresponding to the model indix provided (e.g. "QBCP" if model == Csts::QBCP).
 **/
string init_F_MODEL(int model)
{
    switch(model)
    {
    case Csts::QBCP:
        return "QBCP";
    case Csts::BCP:
        return "BCP";
    case Csts::CRTBP:
        return "RTBP";
    case Csts::ERTBP:
        return "ERTBP";
    default:
        cout << "init_F_MODEL. Warning: unknown model." << endl;
    }
    return "QBCP"; //never here
}

/**
 *  \brief Return the string corresponding to the framework (coord. syst.) indix provided (e.g. "EM" if coordsys == Csts::EM).
 **/
string init_F_COORDSYS(int coordsys)
{
    switch(coordsys)
    {
    case Csts::EM:
        return  "EM";
    case Csts::SEM:
        return  "SEM";
    case Csts::SE:
        return  "SE";
    default:
        cout << "init_F_COORDSYS. Warning: unknown model." << endl;
    }
    return "EM"; //never here
}

/**
 *  \brief Return the folder name corresponding to the prefix/model/framework/libration point number combination provided (e.g. "prefix/QBCP/EM/L1").
 **/
string init_F_FOLDER(string prefix, int model, int coordsys, int li)
{
    return prefix+"/"+init_F_MODEL(model)+"/"+init_F_COORDSYS(coordsys)+"/"+init_F_LI(li)+"/";
}

/**
 * \brief Retrieve a set of coefficients, given as Fourier series from a txt file.
 * \param filename the name of the txt file.
 * \param params a pointer toward the array to update.
 * \param nf the order of the Fourier series.
 * \param shift the indix from which to start the storage of the coefficients in params.
 * \param flag: if flag == 1, the coefficients computed via FFT are used. Otherwise, the expansions obtained through Fourier series algebraic manipulations are used.
 *
 *  Warning: As of now, FFT coefficients must be used for betas and deltas (see QBCP_L structure).
 **/
void coefRetrieving(string filename, double *params, int nf, int shift, int flag, int number)
{
    //Reading tools
    ifstream readStream;
    double cDouble1;
    int alphaNumber = 1;
    string ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    for(int header = shift; header <= (nf+1)*(number-1)+shift; header+=(nf+1))
    {
        if(flag) readStream.open((filename+ss+"c_fft.txt").c_str());
        else readStream.open((filename+ss+"c.txt").c_str());
        for(int i=0; i<= nf; i++)
        {
            readStream >> cDouble1;  //current order
            readStream >> params[i+header];
        }
        readStream.close();
        alphaNumber++;
        ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    }

    if(!flag) cout << "coefRetrieving: the FFT coefficients have not been used." << endl;
}

/**
 * \brief Compute the potential energy for the given state and CR3BP (mu).
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double energy(double y[], double mu)
{
    double r1 = sqrt( pow(y[0]+mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );
    double r2 = sqrt( pow(y[0]- 1 + mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );

    return - ( 1.0/2.0*(pow(y[0],2) + pow(y[1],2)) + (1-mu)/r1 + mu/r2 + 1.0/2.0*mu*(1-mu) );
}

/**
 * \brief Compute the Jacobi constant for the given state and CR3BP (mu)
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double jacobi(double y[], double mu)
{
    return -2.0*energy(y, mu) - (pow(y[3],2.0)+pow(y[4],2.0)+pow(y[5],2.0));
}

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param qbcp_l a reference to the QBCP_L initialized around the selected libration point.
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn(QBCP_L& qbcp_l, int n)
{
    double gamma = qbcp_l.cs.gamma;
    double mu;
    switch(qbcp_l.coordsys)
    {
        case Csts::EM:
            mu   = qbcp_l.us.mu_EM;
        break;
        case Csts::SEM:
            mu   = qbcp_l.us.mu_SEM;
        break;
        case Csts::SE:
            mu   = qbcp_l.us.mu_SE;
        break;
        default: //EM by default
            cout << "WARNING in cn(): unknown framework. EM by default." << endl;
            mu   = qbcp_l.us.mu_EM;
    }


    double res = 0.0;
    switch(qbcp_l.cs.li)
    {
    case 1:
        res =  pow(gamma,-3.0)*(pow(+1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0-gamma), n+1));
        break;
    case 2:
        res =  pow(gamma,-3.0)*(pow(-1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0+gamma), n+1));
        break;
    case 3:
        res =  pow(gamma,-3.0)*(1.0 - mu + mu*pow(gamma/1.0+gamma, n+1)); //convention by Richardson 1980
        break;
    default:
        cout << "cn. Warning: supplied Li number is out of scope. 0.0 is returned." << endl;
    }
    return res;
}

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param li the number of the current libration point (1,2,3)
 * \param gamma the gamma parameter associated to the current libration point
 * \param mu the mass ratio of the current TBP system
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn(int li, double gamma, double mu, int n)
{
    double res = 0.0;
    switch(li)
    {
    case 1:
        res =  pow(gamma,-3.0)*(pow(+1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0-gamma), n+1));
        break;
    case 2:
        res =  pow(gamma,-3.0)*(pow(-1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0+gamma), n+1));
        break;
    case 3:
        res =  pow(gamma,-3.0)*(1.0 - mu + mu*pow(gamma/1.0+gamma, n+1)); //convention by Richardson 1980
        break;
    default:
        cout << "cn. Warning: supplied Li number is out of scope. 0.0 is returned." << endl;
    }
    return res;
}

/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number)
{
    void nrerror(char error_text[]);
    int j;
    double df=0.0, f = 0.0;
    double dx,rtn;

    rtn=x1;   //initial guess

    for (j=1; j<=50; j++)
    {
        (*funcd)(mu, number, rtn,&f,&df);
        dx=f/df;
        rtn -= dx;

        if (fabs(dx) < xacc) return rtn;  //Convergence
    }
    printf("WARNING: Maximum number of iterations exceeded in rtnewt");
    return 0.0;   //Never get here.
}

/**
* \brief Provides the function value and its first derivative for the newton-raphson method.
* f corresponds to the equation satisfied by the Li-m2 distance for the L1/L2 cases
* and by 1-(Li-m1 distance) for the L3 case
**/
void polynomialLi(double mu, int number, double y, double *f, double *df)
{
    switch(number)
    {

    case 1:
        *f =  pow(y,5.0)   - (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) +  2*mu*y - mu;
        *df = 5*pow(y,4.0) - 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    +  2*mu;
        break;

    case 2:
        *f =  pow(y,5.0)   + (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) -  2*mu*y - mu;
        *df = 5*pow(y,4.0) + 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    -  2*mu;
        break;

    case 3:
        //*f =  pow(y,5.0) + (2.0+mu)*pow(y,4.0) + (1+2*mu)*pow(y,3.0) + (1+mu)*pow(y,2.0) +  2*(1-mu)*y + 1-mu;
        *f= pow(y,5.0) + (7+mu)*pow(y,4.0) + (19+6*mu)*pow(y,3.0) -(24+13*mu)*pow(y,2.0) +  (12+14*mu)*y -7*mu;
        *df= 5*pow(y,4.0) + 4*(7+mu)*pow(y,3.0) + 3*(19+6*mu)*pow(y,2.0) -2*(24+13*mu)*pow(y,1.0) +  (12+14*mu);
        //*df = 5*pow(y,4.0) + 4*(2.0+mu)*pow(y,3.0) + 3*(1+2*mu)*pow(y,2.0) + 2*(1+mu)*y +  2*(1-mu);
        break;
    }
}


/**
 *  \brief Prompt "Press Enter to go on"
 **/
void pressEnter(bool isFlag)
{
    if(isFlag)
        {
        char ch;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
    }
}


/**
 *  \brief Prompt msg
 **/
void pressEnter(bool isFlag, string msg)
{
    if(isFlag)
    {
        char ch;
        printf((char*)msg.c_str());
        scanf("%c",&ch);
    }
}
