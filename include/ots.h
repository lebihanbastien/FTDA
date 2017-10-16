#ifndef OTS_H_INCLUDED
#define OTS_H_INCLUDED

#include <iostream>
#include <iomanip>

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "ftda.h"
#include "otsh.h"
#include "parameters.h"
#include "Config.h"

/**
 * \file ots.h
 * \brief Taylor series template class
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */


//--------------------------------------------------------------------------------//
//  Important note: the delete of Ots objects is badly handled because
//  I do not konw how to properly delete a recursive object such as every
//  Otsh in **term.
//  As a consequence, the Ots object must be handled with care (do not create to many)
//--------------------------------------------------------------------------------//
using namespace std;

template <typename T>
class Ots;

typedef Ots<cdouble> Otsc;

template<typename T>
std::ostream& operator << (std::ostream& stream, Ots<T> const& ots);

template <typename T>
class Ots
{

private:

    int nv;           //number of variables
    int order;        //order of the taylor_serie
    Otsh<T> **term;   //list of homogeneous terms

public:
    //Create
    Ots<T>();
    Ots<T>(int newNv, int newOrder);
    Ots<T>(int newNv, int newOrder, T *coef0);
    Ots<T>(int newOrder);

    //Delete
    ~Ots<T>();

    //Copy
    Ots<T>(Ots<T> const& b);                  //unlinked copy
    Ots<T>& operator = (Ots<T> const& b);     //unlinked copy
    Ots<T>& lcopy(Ots<T> const& b);           //linked copy (probably useless at this step)
    Ots<T>& ccopy(Ots<T> const& b);           //coefficient copy (restricted to same order, same number of variables)

    //Zeroing
    void zero();

    //Getters
    int getOrder() const;
    int getNV() const;
    T getCoef(int const& mOrder, int const& pos) const;
    Otsh<T>* getTerm(int const& mOrder) const;

    //Setters
    void setAllCoefs(T const& m);
    void setCoef (T const& m, int const& pos);
    void setCoef (T const& m, int pos, int i);
    void conjugate();
    void setRandomCoefs();

    void conjugateforOFS();

    //Functions
    Ots<T>& smprod(Ots<T> const& a, Ots<T> const& b, T const& m, int const& n);
    Ots<T>& smprod(Ots<T> const& a, Ots<T> const& b, T const& m);
    Ots<T>&  mprod(Ots<T> const& a, Ots<T> const& b, T const& m);
    Ots<T>&  sprod(Ots<T> const& a, Ots<T> const& b);
    Ots<T>&  sprod(Ots<T> const& a, Ots<T> const& b, int const& n);
    Ots<T>&   prod(Ots<T> const& a, Ots<T> const& b);

    Ots<T>& smult(Ots<T> const& a, T const& m);
    Ots<T>&  mult(Ots<T> const& a, T const& m);
    Ots<T>&  mult(T const& m);

    Ots<T>& fsum(Ots<T> const& a,  T const& ma, Ots<T> const& b, T const& mb);

    Ots<T>&  divs(Ots<T> const& a, Ots<T> const& b);
    Ots<T>& mdivs(Ots<T> const& a, Ots<T> const& b, T const& m);

    Ots<T>& pows(Ots<T> const& a,  T const& alpha);

    void acoef0s(T const& x0);
    T coef0s(Ots<T> const& a);

    //Derivative
    Ots<T>&  der(Ots< T > const& a, int ni);
    Ots<T>& sder(Ots< T > const& a, int ni);
    Ots<T>&  der(Ots< T > const& a, int ni, int m);

    //Friendly streaming
    friend std::ostream& operator << <>(std::ostream& stream, Ots<T> const& ots);

    //Evaluate
    template<typename U> void evaluate(U X[], T& z);
    template<typename U> void evaluate_conjugate(U X[], T& z);
};

//Function
template<typename T> T cpow(T const& x, T const& alpha);
template<typename T>   Ots<T> operator / (Ots<T> const& a, Ots<T> const& b);
template <typename T>  Ots<T> operator * (T const& c, Ots<T> const& a);

//Include the implementation .tpp
#include "ots.tpp"


#endif // OTS_H_INCLUDED


