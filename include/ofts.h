#ifndef OFTS_H_INCLUDED
#define OFTS_H_INCLUDED


#include <iostream>
#include <iomanip>

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "ftda.h"
#include "oftsh.h"
#include "ofs.h"


//--------------------------------------------------------------------------------//
//  Important note: the delete of Ofts objects is badly handled because
//  I do not konw how to properly delete a recursive object such as every
//  Oftsh in **term.
//  As a consequence, the Ofts object must be handled with care (do not create to many)
//--------------------------------------------------------------------------------//
using namespace std;

template <typename T>
class Ofts;

template<typename T>
std::ostream& operator << (std::ostream& stream, Ofts<T> const& ofts);

template <typename T>
class Ofts
{

private:

    int nv;         //number of variables
    int order;      //order of the fourier_taylor_serie

    int cnv;        //number of variables of the coefficients
    int corder;     //order of the coefficients


    Oftsh<T> **term;   //list of homogeneous terms

public:
    //Create
    Ofts<T>(int newNv, int newOrder, int newCnv, int newCorder);
    Ofts<T>(int newNv, int newOrder, int newCnv, int newCorder, T *coef0);

    //Delete
    ~Ofts<T>();

    //Copy
    Ofts<T>(Ofts<T> const& b);                  //unlinked copy
    Ofts<T>& operator = (Ofts<T> const& b);     //unlinked copy
    Ofts<T>& lcopy(Ofts<T> const& b);           //linked copy (probably useless at this step)
    Ofts<T>& ccopy(Ofts<T> const& b);           //coefficient copy (restricted to same order, same number of variables)
    Ofts<T>& ccopy(Ofts<T> const& b, int const& nrc);           //coefficient copy (restricted to same order, same number of variables)


    //Zeroing
    void zero();

    //Getters
    int getOrder() const;
    int getVariables() const;
    int getCOrder() const;
    int getCVariables() const;
    T* getCoef(int const& mOrder, int const& pos) const;
    Oftsh<T>* getTerm(int const& mOrder) const;


    //Setters
    void setCoefs(T const& m);
    template <typename U> void setCoefs(U const& m);
    template < typename U > void setCoef(U const& m, int const& n);
    Ofts<T>& conjugate();
    Ofts<T>& conjugate(int const& nrc);
    void setRandomCoefs();


    //Functions
    Ofts<T>&                      smprod(Ofts<T> const& a, Ofts<T> const& b, T const& m);
    Ofts<T>&                       mprod(Ofts<T> const& a, Ofts<T> const& b, T const& m);
    Ofts<T>&                       sprod(Ofts<T> const& a, Ofts<T> const& b);
    Ofts<T>&                       sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n);
    Ofts<T>&                        prod(Ofts<T> const& a, Ofts<T> const& b);
    template<typename U> Ofts<T>& smprod(Ofts<T> const& a, Ofts<T> const& b, U const& c);
    template<typename U> Ofts<T>&  mprod(Ofts<T> const& a, Ofts<T> const& b, U const& c);

    Ofts<T>&                             smult(Ofts<T> const& a, T const& m);
    Ofts<T>&                              mult(Ofts<T> const& a, T const& m);
    template<typename U> Ofts< Ots<U> >& smult(Ofts< Ots<U> > const& a, U const& c);
    template<typename U> Ofts< Ots<U> >& smult(Ofts< Ots<U> > const& a, T const& m, U const& c);
    template<typename U> Ofts< Ots<U> >&  mult(Ofts< Ots<U> > const& a, U const& c);
    template<typename U> Ofts< Ots<U> >&  mult(Ofts< Ots<U> > const& a, T const& m, U const& c);
    template<typename U> Ofts< Ots<U> >& smult(Ofts< Ots<U> > const& a, U const& c, int const& k);
    template<typename U> Ofts< Ots<U> >& smult(Ofts< Ots<U> > const& a, T const& m, U const& c, int const& k);

    template<typename U> Ofts< Ofs<U> >& smult(Ofts< Ofs<U> > const& a, U const& c);
    template<typename U> Ofts< Ofs<U> >&  mult(Ofts< Ofs<U> > const& a, U const& c);
    template<typename U> Ofts< Ofs<U> >&  mult(Ofts< Ofs<U> > const& a, T const& m, U const& c);
    template<typename U> Ofts< Ofs<U> >& smult(Ofts< Ofs<U> > const& a, U const& c, int const& k);
    template<typename U> Ofts< Ofs<U> >& smult(Ofts< Ofs<U> > const& a, T const& m, U const& c);
    template<typename U> Ofts< Ofs<U> >& smult(Ofts< Ofs<U> > const& a, T const& m, U const& c, int const& k);

    Ofts<T>& fsum(Ofts<T> const& a,  T const& ma, Ofts<T> const& b, T const& mb);
    Ofts<T>&  divs(Ofts<T> const& a, Ofts<T> const& b);

    template<typename U> Ofts<T>& pows(Ofts<T> const& a,  U const& alpha);
    template<typename U> Ofts<T>& pows(Ofts<T> const& a,  U const& alpha, int const& n);

    void acoef0s(T const& x0);
    T* coef0s(Ofts<T> const& a);

    //Friendly streaming
    friend std::ostream& operator << <>(std::ostream& stream, Ofts<T> const& ofts);
};

//Function
template <typename U>  void fts2fs(Ofs<U> *fs, Ofts< Ots<U> > const& fts_z, double epsilon);
template <typename U>  double complex fts2scalar(Ofts< Ots<U> > const& fts_z, double epsilon, double time);
template <typename U>  double complex fts2scalar(Ofts< Ofs<U> > const& fts_z, double epsilon, double t);

//Include the implementation .tpp
#include "ofts.tpp"



#endif // OFTS_H_INCLUDED
