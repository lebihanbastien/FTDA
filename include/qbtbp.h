#ifndef QBTBP_H_INCLUDED
#define QBTBP_H_INCLUDED

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>


#include "ots.h"
#include "ofts.h"
#include "ofs.h"

void qbtbp(Ofts< Ots<double> > &ofts_z, Ofts< Ots<double> > &ofts_Z, double n, double ms, double as, double ns);


#endif // QBTBP_H_INCLUDED
