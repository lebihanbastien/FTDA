#include "multimin_test.h"
using namespace std;


//function to evaluate
double funk(double p[], double b[], double c, void * params)
{
    return ((p[0] - b[0])*(p[0] - b[0]) + (p[1] - b[1])*(p[1] - b[1]) + (p[2] - b[2])*(p[2] - b[2]));
}

//gradient of the function funk
void dfunk(double p[], double df[], double b[], double c, void *params )
{
    for(int i = 0; i< 3; i++) df[i] = 2*(p[i] - b[i]);

}

//Test function of the routine frprmn of multimin.c
void test_frprmn()
{
    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the routine frprmn          " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(6);

    //Initialize starting point @ zero
    double *p = dvector(0, 2);
    for(int i = 0; i <3; i++) p[i] = 0.0;

    //Target point (arbitrary)
    double *b = dvector(0, 2);
    b[0] = +1.0;
    b[1] = -1.0;
    b[2] = +0.5;

    //Parameters (needed to compute funk and dfunk)
    double t = 0.0;
    void *params = 0;
    //Parameters of the multimin routine frprmn
    double ftol = 1e-12;
    int iter = 0;
    double fret = 0.0;

    //Multimin routine
    tic();
    frprmn(b, t, params, p, 4, ftol, &iter, &fret, funk, dfunk);
    cout << "End of min. in " << toc() << "s." << endl;

    cout << "Minimum value of the function: " << fret << endl;
    cout << "Target      Location of the minimum: " << endl;
    for(int i =0; i< 3; i++) cout << b[i] << "    " << p[i] << endl;


    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the routine dfrprmn          " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(6);
    //Initialize starting point @ zero
    for(int i = 0; i <3; i++) p[i] = 0.0;

    //Multimin routine
    tic();
    dfrprmn(b, t, params, p, 4, ftol, &iter, &fret, funk, dfunk);
    cout << "End of min. in " << toc() << "s." << endl;

    cout << "Minimum value of the function: " << fret << endl;
    cout << "Target      Location of the minimum: " << endl;
    for(int i =0; i< 3; i++) cout << b[i] << "    " << p[i] << endl;


    cout << "---------------------------------------------------" << endl;
    cout << "                                                   " << endl;
    cout << "               Test of the routine amoeba          " << endl;
    cout << "                                                   " << endl;
    cout << "---------------------------------------------------" << endl;
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(6);
    double **pp   = dmatrix(0, 3, 0, 2);
    double **id   = dmatrix(0, 2, 0, 2);
    double *y     = dvector(0, 3);
    double lambda = 0.1;


    //Identity matrix (ej)
    for(int i= 0; i < 3; i++) for(int j = 0; j < 3; j++) id[i][j] = (i==j)? 1:0;
    //First point is b
    for(int j = 0; j < 3; j++) pp[0][j] = (double) rand()/RAND_MAX;
    //4 others are b+lambda*ej
    for(int i= 1; i <= 3; i++)
    {
        for(int j = 0; j < 3; j++) pp[i][j] = pp[0][j] + lambda*id[i-1][j];
    }
    //Evaluate the 4+1 vertices and store them in y
    for(int i = 0; i <= 3; i++) y[i] = funk(pp[i], b, t, params);


    cout << "Initial vertices: " << endl;
    for(int i =0; i< 3; i++) cout << pp[0][i] << "   "  << pp[1][i] << "   "  << pp[2][i] << "   "  << pp[3][i] << endl;


    int nfunk = 0;
    //multimin routine amoeba
    tic();
    amoeba(b, t, params, pp, y, 3, ftol, funk, &nfunk);
    cout << "End of min. in " << toc() << "s." << endl;


    cout << "Target   |   Location of the minimum: " << endl;
    for(int i =0; i< 3; i++) cout << b[i] << "  |  " << pp[0][i] << "   "  << pp[1][i] << endl;

    //Free memory
    free_dvector(p, 0, 2);
    free_dvector(b, 0, 2);

}


