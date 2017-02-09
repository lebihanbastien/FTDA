#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
using namespace std;

/**
 * \file Config.h
 * \brief List of some parameters that are used throughout numerical compuration,
 *        such as precision in numerical integrators. Defines a Config Manager that allows
 *        to use common values accross the software.
 * \author BLB
 * \date JAN 2017
 */

class Config
{
    private:

        //--------------------------------------------------------------------------------
        //Constants for ODE structure (see ode.h & cpp)
        //--------------------------------------------------------------------------------
        double PREC_HSTART;  //Initial step in numerical integrator
        double PREC_ABS;     //Absolute precision in numerical integrator
        double PREC_REL;     //Relative precision in numerical integrator
        double PREC_ROOT;    //Precision on root (zero) finding
        double PREC_LIB;     //Precision on the position of the libration point (maybe redundant with _ROOT)
        double DELTA_T;      //Very small delta of time necessary to avoid some errors in numerical procedure at t = 0.0.

        //--------------------------------------------------------------------------------
        //Constants for differential correction procedures
        //--------------------------------------------------------------------------------
        int DC_ITERMAX;

    public:
        static Config& configManager()
        {
            static Config instance;
            return instance;
        }

        //--------------------------------------------------------------------------------
        //Getters
        //--------------------------------------------------------------------------------
        double G_PREC_ABS(){return PREC_ABS;}
        double G_PREC_REL(){return PREC_REL;}
        double G_PREC_ROOT(){return PREC_ROOT;}
        double G_PREC_LIB(){return PREC_LIB;}
        double G_PREC_HSTART(){return PREC_HSTART;}
        double G_DELTA_T(){return DELTA_T;}

        int    G_DC_ITERMAX(){return DC_ITERMAX;}

        //--------------------------------------------------------------------------------
        //Setters
        //--------------------------------------------------------------------------------
        /**
         *  \brief Hard precision for integration purposes
         **/
        void C_PREC_HARD();

        /**
         *  \brief Soft precision for integration purposes, to fasten correction loops
         *         It is the default case.
         **/
        void C_PREC_SOFT();

    private:
        Config();
        Config(Config const&);
        ~Config();
        void operator=(Config const&);
};

#endif // CONFIG_H
