#include "multimin.h"

//---------------------------------------------------------------------------
#define TINY 1.0e-10 //A small number.
#define EPS  1.0e-10 //A small number.
#define NMAX 5000    //Maximum allowed number of function evalua-
#define ITMAX 200
#define GET_PSUM \
for (j=0;j<ndim;j++) {\
    for (sum=0.0,i=0; i<mpts; i++) sum += p[i][j];\
    psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define FREEALL free_dvector(xi,0,n-1);free_dvector(h,0,n-1);free_dvector(g,0,n-1)
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f)

#define TOL 2.0e-4 //Tolerance passed to dbrent.

#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY_MNBRAK 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))
//---------------------------------------------------------------------------


//Global variables communicate with df1dim.
int ncom;
double *pcom,*xicom,(*nrfunk)(double [], double[], double, void *);
void (*nrdfunk)(double [], double[], double[], double, void *);


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Fletcher-Reeves-Polak-Ribiere minimization
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/*
Given a starting point p[0..n-1] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations.
*/
void frprmn(double z1[], double t, void *params, double p[], int n, double ftol, int *iter, double *fret,
            double (*funk)(double [], double[], double, void *), void (*dfunk)(double [], double[], double[], double, void *))
{
    int j,its;
    double gg,gam,fp,dgg;
    double *g,*h,*xi;

    g =  dvector(0,n-1);
    h =  dvector(0,n-1);
    xi=  dvector(0,n-1);

    //Initializations.
    fp=(*funk)(p, z1, t, params);
    (*dfunk)(p,xi,z1, t, params);

    for (j=0; j<n; j++)
    {
        g[j] = -xi[j];
        xi[j]=h[j]=g[j];
    }

    //Loop over iterations
    for (its=1; its<=ITMAX; its++)
    {
        *iter=its;
        linmin(z1, t, params, p,xi,n,fret,funk); //Next statement is the normal return:
        if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS))
        {
            //FREEALL;
            free_dvector(xi,0,n-1);
            free_dvector(h,0,n-1);
            free_dvector(g,0,n-1);
            return;
        }
        fp= *fret;
        (*dfunk)(p,xi,z1, t, params);
        dgg=gg=0.0;
        for (j=0; j<n; j++)
        {

            gg += g[j]*g[j];
            //dgg += xi[j]*xi[j];       //This statement for Fletcher-Reeves.
            dgg += (xi[j]+g[j])*xi[j];  //This statement for Polak-Ribiere.
        }
        if (gg == 0.0)
        {
            //FREEALL;
            free_dvector(xi,0,n-1);
            free_dvector(h,0,n-1);
            free_dvector(g,0,n-1);
            return;
        }
        gam=dgg/gg;
        for (j=0; j<n; j++)
        {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
    }
    nrerror((char*)"Too many iterations in frprmn");
}

/*
Given a starting point p[0..n-1] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine dlinmin is called to perform line minimizations.
*/
void dfrprmn(double z1[], double t, void *params, double p[], int n, double ftol, int *iter, double *fret,
            double (*funk)(double [], double[], double, void *), void (*dfunk)(double [], double[], double[], double, void *))
{
    int j,its;
    double gg,gam,fp,dgg;
    double *g,*h,*xi;

    g =  dvector(0,n-1);
    h =  dvector(0,n-1);
    xi=  dvector(0,n-1);

    //Initializations.
    fp=(*funk)(p, z1, t, params);
    (*dfunk)(p,xi,z1, t, params);

    for (j=0; j<n; j++)
    {
        g[j] = -xi[j];
        xi[j]=h[j]=g[j];
    }

    //Loop over iterations
    for (its=1; its<=ITMAX; its++)
    {
        *iter=its;
        dlinmin(z1, t, params, p,xi,n,fret,funk, dfunk); //Next statement is the normal return:
        if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS))
        {
            //FREEALL;
            free_dvector(xi,0,n-1);
            free_dvector(h,0,n-1);
            free_dvector(g,0,n-1);
            return;
        }
        fp= *fret;
        (*dfunk)(p,xi,z1, t, params);
        dgg=gg=0.0;
        for (j=0; j<n; j++)
        {

            gg += g[j]*g[j];
            //dgg += xi[j]*xi[j];       //This statement for Fletcher-Reeves.
            dgg += (xi[j]+g[j])*xi[j];  //This statement for Polak-Ribiere.
        }
        if (gg == 0.0)
        {
            //FREEALL;
            free_dvector(xi,0,n-1);
            free_dvector(h,0,n-1);
            free_dvector(g,0,n-1);
            return;
        }
        gam=dgg/gg;
        for (j=0; j<n; j++)
        {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
    }
    nrerror((char*)"Too many iterations in frprmn");
}

/*
Given an n -dimensional point p[0..n-1] and an n -dimensional direction xi[0..n-1] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and dbrent .
*/
void dlinmin(double z1[], double t, void *params, double p[], double xi[], int n, double *fret, double (*funk)(double [], double[], double, void *),
             void (*dfunk)(double [], double[], double[], double, void *))
{
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    ncom=n; //Define the global variables.
    pcom  = dvector(0,n-1);
    xicom = dvector(0,n-1);
    nrfunk=funk;
    nrdfunk=dfunk;
    for (j=0; j<n; j++)
    {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0; //Initial guess for brackets.
    xx=1.0;
    mnbrak(z1, t, params, &ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=dbrent(z1, t, params, ax,xx,bx,f1dim,df1dim,TOL,&xmin);
    for (j=0; j<n; j++)  //Construct the vector results to return.
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_dvector(xicom,0,n-1);
    free_dvector(pcom,0,n-1);
}


/*
Given an n -dimensional point p[0..n-1] and an n -dimensional direction xi[0..n-1] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and brent .
*/
void linmin(double z1[], double t, void *params, double p[], double xi[], int n, double *fret, double (*funk)(double [], double[], double, void *))
{
//    double brent(double[], double, void*, double ax, double bx, double cx,
//    double (*f)(double[], double, void *, double), double tol, double *xmin);
    double f1dim(double z1[], double t, void *params, double x);
    void mnbrak(double[], double, void*, double *ax, double *bx, double *cx, double *fa, double *fb,
    double *fc, double (*func)(double[], double, void *, double));
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    ncom=n; //Define the global variables.
    pcom  = dvector(0,n-1);
    xicom = dvector(0,n-1);
    nrfunk=funk;
    for (j=0; j<n; j++)
    {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0; //Initial guess for brackets.
    xx=1.0;
    mnbrak(z1, t, params,&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=brent(z1, t, params,ax,xx,bx,f1dim,TOL,&xmin);
    for (j=0; j<n; j++) //Construct the vector results to return.
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_dvector(xicom,0,n-1);
    free_dvector(pcom,0,n-1);
}


/*
  Artificial function of one variable along one direction of minimization in the routine linmin
*/
double f1dim(double z1[], double t, void *params, double x)

{
    int j;
    double f,*xt;
    xt=dvector(0,ncom-1);
    for (j=0; j<=ncom-1; j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunk)(xt, z1, t, params);
    free_dvector(xt,0,ncom-1);
    return f;
}

/*
  Gradient of an artificial function of one variable along one direction of minimization in the routine linmin
*/
double df1dim(double z1[], double t, void *params, double x)
{
    int j;
    double df1=0.0;
    double *xt,*df;
    xt = dvector(0,ncom-1);
    df = dvector(0,ncom-1);
    for (j=0; j<ncom; j++) xt[j]=pcom[j]+x*xicom[j];
    (*nrdfunk)(xt,df, z1, t, params);
    for (j=0; j<ncom; j++) df1 += df[j]*xicom[j];
    free_dvector(df,0,ncom-1);
    free_dvector(xt,0,ncom-1);
    return df1;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  1D min
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/*
Given a function func , and given distinct initial points ax and bx , this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax , bx , cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa , fb , and fc .
*/
void mnbrak(double z1[], double t, void *params, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double[], double, void *, double))
{
    double ulim,u,r,q,fu,dum;
    *fa=(*func)(z1, t, params, *ax);
    *fb=(*func)(z1, t, params, *bx);
//Switch roles of a and b so that we can go downhill in the direction from a to b.
    if (*fb > *fa)
    {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax); //First guess for c.
    *fc=(*func)(z1, t, params, *cx);
    while (*fb > *fc)   //Keep returning here until we bracket.
    {

//Compute u by parabolic extrapolation from
//a, b, c. TINY_MNBRAK is used to prevent any pos-
//sible division by zero.
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(DMAX(fabs(q-r),TINY_MNBRAK),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
//We won’t go farther than this. Test various possibilities:
        if ((*bx-u)*(u-*cx) > 0.0)   //Parabolic u is between b and c: try it.
        {
            fu=(*func)(z1, t, params, u);
            if (fu < *fc)   //Got a minimum between b and c.
            {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            }
            else if (fu > *fb)     //Got a minimum between between a and u.
            {
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx); //Parabolic fit was no use. Use default magnification.
            fu=(*func)(z1, t, params, u);

        }
        else if ((*cx-u)*(u-ulim) > 0.0)     //Parabolic fit is between c and its allowed limit.
        {
            fu=(*func)(z1, t, params, u);
            if (fu < *fc)
            {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,(*func)(z1, t, params, u))
            }
        }
        else if ((u-ulim)*(ulim-*cx) >= 0.0)     //Limit parabolic u to maximum allowed value.
        {
            u=ulim;

            fu=(*func)(z1, t, params, u);
        }
        else     //Reject parabolic u, use default magnification.
        {
            u=(*cx)+GOLD*(*cx-*bx);

            fu=(*func)(z1, t, params, u);
        }
        SHFT(*ax,*bx,*cx,u) //Eliminate oldest point and continue.
        SHFT(*fa,*fb,*fc,fu)
    }
}

/*
Brent minimization
Given a function f , and given a bracketing triplet of abscissas ax , bx , cx (such that bx is
between ax and cx , and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value.
*/
double brent(double z1[], double t, void *params, double ax, double bx, double cx, double (*f)(double[], double, void *, double), double tol, double *xmin)

{
    int iter;
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
//This will be the distance moved on
//the step before last.
    double e=0.0;
//a and b must be in ascending order,
//but input abscissas need not be.

//attention from gc++
d = 0.0;
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;  //Initializations...
    fw=fv=fx=(*f)(z1, t, params, x);
    for (iter=1; iter<=ITMAX; iter++) //Main program loop.
    {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))   //Test for done here.
        {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1)   //Construct a trial parabolic fit.
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) d = CGOLD*(e=(x >= xm ? a-x : b-x));
//The above conditions determine the acceptability of the parabolic fit. Here we
//take the golden section step into the larger of the two segments.
            else
            {
                d=p/q;      //Take the parabolic step.
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        }
        else
        {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(z1, t, params, u); //This is the one function evaluation per iteration.
        if (fu <= fx)    //Now decide what to do with our function evaluation.
        {
            if (u >= x) a=x;
            else b=x;

            SHFT(v,w,x,u) //Housekeeping follows:
            SHFT(fv,fw,fx,fu)
        }
        else
        {
            if (u < x) a=u;
            else b=u;
            if (fu <= fw || w == x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v=u;
                fv=fu;
            }
        }  //Done with housekeeping. Back for another iteration.
    }
    nrerror((char*)"Too many iterations in brent");
    *xmin=x;  //Never get here.
    return fx;
}


/*
Brent minimization with gradient
Given a function f and its derivative function df , and given a bracketing triplet of abscissas ax ,
bx , cx [such that bx is between ax and cx , and f(bx) is less than both f(ax) and f(cx) ],
this routine isolates the minimum to a fractional precision of about tol using a modification of
Brent’s method that uses derivatives. The abscissa of the minimum is returned as xmin , and
the minimum function value is returned as dbrent , the returned function value.
*/

double dbrent(double z1[], double t, void *params, double ax, double bx, double cx, double (*f)(double[], double, void *, double),
              double (*df)(double[], double, void *, double), double tol, double *xmin)

{
    int iter,ok1,ok2;  //Will be used as flags for whether proposed steps are acceptable or not.
    double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
    double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

    //attention from gc++
    d = 0.0;

//Comments following will point out only differences from the routine brent. Read that
//routine first.
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(z1, t, params, x);
    dw=dv=dx=(*df)(z1, t, params, x);
    for (iter=1; iter<=ITMAX; iter++)
    {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+ZEPS;
        tol2=2.0*tol1;
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1)
        {
            d1=2.0*(b-a);
            d2=d1;
            if (dw != dx) d1=(w-x)*dx/(dx-dw);
            if (dv != dx) d2=(v-x)*dx/(dx-dv);
            u1=x+d1;
            u2=x+d2;
            ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
            ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
            olde=e;
            e=d;
            if (ok1 || ok2)
            {
                if (ok1 && ok2)
                    d=(fabs(d1) < fabs(d2) ? d1 : d2);
                else if (ok1)
                    d=d1;
                else
                    d=d2;
                if (fabs(d) <= fabs(0.5*olde))
                {
                    u=x+d;
                    if (u-a < tol2 || b-u < tol2)
                        d=SIGN(tol1,xm-x);
                }
                else
                {
                    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
                }
            }
            else
            {
                d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
            }
        }
        else
        {
            d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
        }
        if (fabs(d) >= tol1)
        {
            u=x+d;
            fu=(*f)(z1, t, params, u);
        }
        else
        {
            u=x+SIGN(tol1,d);
            fu=(*f)(z1, t, params, u);
            if (fu > fx)
            {
                *xmin=x;
                return fx;
            }
        }
        du=(*df)(z1, t, params, u);
        if (fu <= fx)
        {
            if (u >= x) a=x;
            else b=x;
            MOV3(v,fv,dv, w,fw,dw);
            MOV3(w,fw,dw, x,fx,dx);
            MOV3(x,fx,dx, u,fu,du);
        }
        else
        {
            if (u < x) a=u;
            else b=u;
            if (fu <= fw || w == x)
            {
                MOV3(v,fv,dv, w,fw,dw);
                MOV3(w,fw,dw, u,fu,du);
            }
            else if (fu < fv || v == x || v == w)
            {
                MOV3(v,fv,dv, u,fu,du);
            }
        }
    }
    nrerror((char*)"Too many iterations in routine dbrent");
    return 0.0;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//
//  Simplex method
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------
/*
Multidimensional minimization of the function funk(x) where x[0..ndim-1] is a vector in ndim
dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[0..ndim]
[0..ndim-1] is input. Its ndim+1 rows are ndim -dimensional vectors which are the vertices of
the starting simplex. Also input is the vector y[0..ndim] , whose components must be pre-
initialized to the values of funk evaluated at the ndim+1 vertices (rows) of p ; and ftol the
fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
nfunk gives the number of function evaluations taken.
*/
void amoeba(double z1[], double t, void *params, double **p, double y[], int ndim, double ftol, double (*funk)(double [], double[], double, void *), int *nfunk)

{
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    double rtol,sum,swap,ysave,ytry,*psum;

    psum= dvector(0,ndim-1);
    *nfunk=0;
    GET_PSUM
    for (;;)
    {
        ilo=0;
        //First we must determine which point is the highest (worst), next-highest, and lowest
        //(best), by looping over the points in the simplex.
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (i=0; i<mpts; i++)
        {
            if (y[i] <= y[ilo]) ilo=i;
            if (y[i] > y[ihi])
            {
                inhi=ihi;
                ihi=i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        //Compute the fractional range from highest to lowest and return if satisfactory.
        if (rtol < ftol)
        {
            //If returning, put best point and value in slot 1.
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[1][i],p[ilo][i])
                break;
        }
        if (*nfunk >= NMAX) nrerror((char*)"NMAX exceeded");
        *nfunk += 2;
        //Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex
        //across from the high point, i.e., reflect the simplex from the high point.
        ytry=amotry(z1, t, params, p,y,psum,ndim,funk,ihi,-1.0);
        if (ytry <= y[ilo])//Gives a result better than the best point, so try an additional extrapolation by a factor 2.
            ytry=amotry(z1, t, params, p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) //The reflected point is worse than the second-highest, so look for an intermediate
        {
            //lower point, i.e., do a one-dimensional contraction.
            ysave=y[ihi];
            ytry=amotry(z1, t, params, p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave)  //Can’t seem to get rid of that high point. Better contract around the lowest (best) point.
            {

                for (i=0; i<mpts; i++)
                {
                    if (i != ilo)
                    {
                        for (j=0; j<ndim; j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=(*funk)(psum, z1, t, params);
                    }
                }
                *nfunk += ndim;  //Keep track of function evaluations.
                //cout << "Evalution n°" << *nfunk << endl;
                GET_PSUM         //Recompute psum.
            }
        }
        else --(*nfunk);    //Correct the evaluation count.
    }    //Go back for the test of doneness and the nextiteration.
    free_dvector(psum,0,ndim-1);
}

/*
Extrapolates by a factor fac through the face of the simplex across from the high point, tries
it, and replaces the high point if the new point is better.
*/
double amotry(double z1[], double t, void *params, double **p, double y[], double psum[], int ndim,
             double (*funk)(double [], double[], double, void *), int ihi, double fac)
{
    int j;
    double fac1,fac2,ytry,*ptry;
    ptry=dvector(0,ndim-1);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry, z1, t, params); //Evaluate the function at the trial point.
    if (ytry < y[ihi])   //If it’s better than the highest, then replace the highest.
    {
        y[ihi]=ytry;
        for (j=0; j<ndim; j++)
        {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    free_dvector(ptry,0,ndim-1);
    return ytry;
}

