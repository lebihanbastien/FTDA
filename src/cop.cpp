#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

/*
    Complex exponentiation
*/
complex cpow(complex x, complex alpha)
{
    double c,d, p, argp;

    c = alpha.real();
    d = alpha.imag();

    p = pow(x.mod(),c)*exp(-d*x.arg());
    argp = c*x.arg() + 0.5*d*log(x.mod()*x.mod());

    return(complex(p*cos(argp), p*sin(argp)));
}
