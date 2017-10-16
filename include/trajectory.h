#ifndef TRAJECTORY_H_INCLUDED
#define TRAJECTORY_H_INCLUDED


/**
 *   \brief Integrates a given trajectory with IC in sti[]
 **/
void trajectory(double sti[], Pmap &pmap, int label);

/**
 *   \brief Integrates a given trajectory up to tf
 **/
int trajectory_integration(Orbit *orbit, double tf);

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(Orbit *orbit,  double t0, double tf, string filename, double **yPlot, double *tplot,  int N, int isResetOn);

/**
 *   \brief Integrates a given trajectory with IC in sti[]
 **/
void trajectory_CU(double sti[], Pmap &pmap, int label, int is3D);

//----------------------------------------------------------------------------------------
//
// Steppers
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(Orbit *orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *ePm,
                   int *nreset,
                   int isResetOn);

/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(Orbit *orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *ePm,
                     int *nreset,
                    int isResetOn);

#endif // TRAJECTORY_H_INCLUDED
