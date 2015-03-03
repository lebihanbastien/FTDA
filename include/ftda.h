#ifndef FTDA_H_INCLUDED
#define FTDA_H_INCLUDED

class FTDA;


class FTDA
{
    private:

        static int nor;         //static maximum degree
        static int nov;         //static number of variables
        static long  int **psi;  //contains psi(i,j) for i=2..nov.

    public:

        static int init(int nv, int nr);
        static void free();
        static long int nmon(int nv, int nr);
        static void prxkt(int k[], int nv);
};


/*
-----------------------------------------------------------------------
Binomial coefficients
-----------------------------------------------------------------------
*/
unsigned long binomial(unsigned long n, unsigned long k);
unsigned long gcd_ui(unsigned long x, unsigned long y);

#endif // FTDA_H_INCLUDED
