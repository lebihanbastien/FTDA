#include "odezero.h"

/**
 * \file odezero.cpp
 * \brief Additional structures for ODE integrationd with event detection, based on GSL implementation.
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

/**
 * \brief return a structure value_output inspired by MATLAB odezero implementation, containing:
 *         - val = y[desired dimension] - desired value
 *         - direction: the direction (increasing, decreasing or both)
 *         - max_events: the maximum events allowed by the user before the end of the integration
 **/
value_output linear_intersection(double t, double yv[], void *params)
{
    value_params *p = (value_params *) params;
    value_output val_par;
    val_par.val        = yv[p->dim] - p->value;
    val_par.direction  = p->direction;
    val_par.max_events = p->max_events;

    return val_par;
}

/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double odezero_event (double t, void *params)
{
    //------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------
    int i;
    struct OdeParams *p = (struct OdeParams *) params;
    double ystart[42];
    for(i=0; i<42; i++) ystart[i] = (p->y0)[i];
    double t0 = p->t0;

    //------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------
    //Starting in the right direction
    p->d->h = (t>t0) ? fabs(p->d->h) : -fabs(p->d->h);
    gsl_odeiv2_driver_apply (p->d, &t0, t, ystart);

    //------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(p->d);

    return p->fvalue->value(t0, ystart, p->fvalue->val_par).val;
}

/**
 * \brief Main routine to stop an integration with a certain condition on the state vector.
 **/
int custom_odezero(double y[], double *tcross, double t1, OdeStruct *ode_s)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int i;                           //loop parameter
    int status;                      //status for GSL function integer output
    double h = 1e-6;                 //first guess for stepper
    double previous_yv[42], yv[42];  //for integration purposes, copy of initial values
    double ys[42];                   //for integration purposes
    double previous_t, t;            //times
    //Root finding
    gsl_function F;                  //called function in the root finding routine
    struct OdeParams params;         //params for F
    double t_low, t_high;            //times for root bracketing
    //Loop parameters
    double r;
    double fy;
    int iter;

    //Copy of initial values for integration purposes
    for(i=0; i<42; i++) yv[i] = y[i];

    //------------------------------------------------------------------------------
    //Reaching y=0
    //------------------------------------------------------------------------------
    double deltat = Config::configManager().G_DELTA_T();
    //Previous value of yv an t are kept
    t = 0.0;
    previous_t = t;
    for(i=0; i<6; i++) previous_yv[i] = yv[i];

    //Stepping until y=0 is crossed
    while (t < deltat || previous_yv[1]*yv[1]>=0)  //NOTE: a condition on time was included to avoid a termination @first step.
    {
        //Updating the previous y and t values
        previous_t = t;
        for(i=0; i<42; i++) previous_yv[i] = yv[i];

        //Evolve one step
        status = gsl_odeiv2_evolve_apply (ode_s->e, ode_s->c, ode_s->s, &ode_s->sys, &t, t1, &h, yv);

        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;

        if(t>=t1-1e-15)
        {
            printf("WARNING: integration time has reached its limit in odezero");
            break;
        }
    }

    //------------------------------------------------------------------------------
    //Root finding: right time to shoot y=0
    //------------------------------------------------------------------------------
    //New start
    for(i=0; i<42; i++) ys[i] = previous_yv[i];

    //Initialization of the parameters for F
    params.t0 = previous_t-1e-15;       //new initial time is previous_t - epsilon
    params.y0 = ys;
    params.d = ode_s->d;

    //Initialization of F
    F.function = &cr3bp_y;
    F.params   = &params;

    //Bracketing the root
    t_low = previous_t;
    t_high = t;

    //Setting the solver
    status = gsl_root_fsolver_set (ode_s->s_root, &F, t_low, t_high);

    //Loop
    fy = ys[1];
    iter = 0;
    do
    {
        status = gsl_root_fsolver_iterate (ode_s->s_root);         //updating the solver
        r = gsl_root_fsolver_root (ode_s->s_root);                 //updating the root
        previous_t = gsl_root_fsolver_x_lower (ode_s->s_root);     //updating t_low
        t = gsl_root_fsolver_x_upper (ode_s->s_root);              //updating t_high

        //Checking convergence
        fy = cr3bp_y (r, &params);
        status = gsl_root_test_residual (fy , ode_s->eps_root);
    }
    while (status == GSL_CONTINUE && (++iter)<100);

    if(iter>=100)
    {
        printf("WARNING: number of iter max exceeded in custom_odezero. Premature ending.\n");
        return GSL_FAILURE;
    }


    //Updating the crossing time
    *tcross = r;

    //Updating the outputs (state + STM)
    t = 0.0;
    status = gsl_odeiv2_driver_apply(ode_s->d, &t, *tcross, y);  //updating y

    //------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------
    reset_ode_structure(ode_s);

    return GSL_SUCCESS;

}

/**
 * \brief return the value y[1] = y for a given integration t from t0 = 0.0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double cr3bp_y (double t, void *params)
{
    //------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------
    int i;
    struct OdeParams *p = (struct OdeParams *) params;
    double ystart[42];
    for(i=0; i<42; i++) ystart[i] = (p->y0)[i];
    double t0 = p->t0;

    //------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------
    gsl_odeiv2_driver_apply (p->d, &t0, t, ystart);

    //------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(p->d);
    return ystart[1];
}
