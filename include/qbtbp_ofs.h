#ifndef QBTBP_OFS_H_INCLUDED
#define QBTBP_OFS_H_INCLUDED

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>


#include "ots.h"
#include "ofts.h"
#include "ofs.h"
#include "custom_ode.h"

void qbtbp_ofs (Ofs<double complex> &bj, Ofs<double complex> &cj, Ofts< Ofs<double> > &ofts_z, Ofts< Ofs<double> > &ofts_Z, double n, double ms, double as, double ns);
int qbtbp_derivatives(double t, const double y[], double f[], void *params);
void analyticVsnumeric(double t1, Ofs<double complex> &bjc, Ofs<double complex> &cjc, custom_ode_structure ode_s, double as, double n, double ns);

#endif // QBTBP_OFS_H_INCLUDED
