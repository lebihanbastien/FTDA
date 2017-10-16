#include "Config.h"

using namespace std;

//----------------------------------------------------------------------------------------
//Constructor
//----------------------------------------------------------------------------------------
/**
 *  \brief Initialize all necessary configuration parameters
 **/
Config::Config()
{
    //------------------------------------------------------------------------------------
    // Parameters for ODE structure (see ode.h & cpp)
    //------------------------------------------------------------------------------------
    PREC_ABS    = 1e-15;
    PREC_REL    = 1e-15;
    PREC_ROOT   = 1e-13;
    PREC_LIB    = 1e-16;
    PREC_HSTART = 1e-8;
    DELTA_T     = 1e-6;

    //------------------------------------------------------------------------------------
    // Parameters for differential correction procedures
    //------------------------------------------------------------------------------------
    DC_ITERMAX  = 20;

    //------------------------------------------------------------------------------------
    // Parameters for aesthetics (cout, plotting)
    //------------------------------------------------------------------------------------
    COUT_SMALL_PREC  = 3;
    COUT_MEDIUM_PREC = 5;
    COUT_LARGE_PREC  = 15;//Large cout precision
}


//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Hard precision for integration purposes. It is the default case.
 **/
void Config::C_PREC_HARD()
{
    PREC_ABS = 1e-15;
    PREC_REL = 1e-15;
    cout << "Config manager: the integration precisions (absolute and relative)";
    cout << "have been set to 1e-15" << endl;
}

/**
 *  \brief Soft precision for integration purposes, to fasten correction loops
 **/
void Config::C_PREC_SOFT()
{
    PREC_ABS = 1e-10;
    PREC_REL = 1e-10;
    cout << "Config manager: the integration precisions (absolute and relative)";
    cout << "have been set to 1e-10" << endl;
}

//----------------------------------------------------------------------------------------
//Destructor
//----------------------------------------------------------------------------------------
Config::~Config()
{
    //dtor
}



//----------------------------------------------------------------------------------------
// Aesthetics
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets a small precision in cout.
 **/
void Config::coutsp()
{
    cout <<  std::setprecision(COUT_SMALL_PREC) << resetiosflags(ios::scientific);
}

/**
 *  \brief Sets an average precision in cout.
 **/
void Config::coutmp()
{
    cout <<  std::setprecision(COUT_MEDIUM_PREC) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a big precision in cout.
 **/
void Config::coutlp()
{
    cout <<  std::setprecision(COUT_LARGE_PREC) << std::showpos  <<  setiosflags(ios::scientific);
}


