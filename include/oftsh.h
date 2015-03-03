#ifndef OFTSH_H_INCLUDED
#define OFTSH_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ftda.h"
#include "ots.h"
#include "ofs.h"

using namespace std;


template <typename T>
class Oftsh;

template<typename T>
std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh);


template <typename T>
class Oftsh
{

private:

    int nv;             //number of variables
    int order;          //degree
    T *coef;            //address of the first coefficient
    Oftsh<T> *term;      //recursivity

    static int globalOrder;
    static int globalVariable;

public:

    //Create
    Oftsh<T>();
    Oftsh<T>(int newNv, int newOrder);
    Oftsh<T>(int newNv, int newOrder, T *coef0);

    //Delete
    ~Oftsh<T>();


    //Copy
    Oftsh<T>(Oftsh<T> const& b);                //unlinked copy
    Oftsh<T>& operator = (Oftsh<T> const& b);   //unlinked copy
    Oftsh<T>& lcopy(Oftsh<T> const& b);         //linked copy (for recursivity)
    Oftsh<T>& ccopy(Oftsh<T> const& b);         //coefficient copy (restricted to same order, same number of variables)


    //Linking
    void setCoefs(T *coef0);

    //Zeroing
    void zero();

    //Setters
    void setCoef(T value, int pos);
    void addCoef(T value, int pos);

    template<typename U> void setTCoef(U value, int pos);
    template<typename U> void setT0Coef(U value, int pos);

    void setRandomCoefs();

    Oftsh<T>& conjugate();

    //Getters
    Oftsh<T> getTerm();
    Oftsh<T> getTerm(int i) const;
    T getCoef(int i);
    int getOrder();
    int getVariables();
    T* getCoefAddress();

    //Operators
    bool isEqual(Oftsh<T> const& b) const;
    /*Oftsh<T>& operator += (Oftsh<T> const& b);
    Oftsh<T>& operator -= (Oftsh<T> const& b);
    Oftsh<T>& operator *= (T const& c);
    Oftsh<T>& operator /= (T const& c);
*/

    Oftsh<T>& smult(Oftsh<T> const& a, T const& c);

    template<typename U>
    Oftsh< Ots<U> >& smult(Oftsh< Ots<U> > const& a, U const& m);
    template<typename U>
    Oftsh< Ofs<U> >& smult(Oftsh< Ofs<U> > const& a, U const& c);

    template<typename U>
    Oftsh< Ots<U> >& smult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, U const& c);
    template<typename U>
    Oftsh< Ofs<U> >& smult(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c);

    template<typename U>
    Oftsh< Ots<U> >& mult(Oftsh< Ots<U> > const& a, U const& c);
    template<typename U>
    Oftsh< Ofs<U> >& mult(Oftsh< Ofs<U> > const& a, U const& c);


    template<typename U>
    Oftsh< Ots<U> >& mult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, U const& c);
    template<typename U>
    Oftsh< Ofs<U> >& mult(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c);


    template<typename U>
    Oftsh< Ots<U> >& smult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, Ots<U> const& r, U const& c,  Ots<U> & temp);



    Oftsh<T>&  mult(Oftsh<T> const& a, T const& c);
    Oftsh<T>&  mult(T const& c);

    Oftsh<T>&  sprod(Oftsh<T> const& a, Oftsh<T> const& b);

    template<typename U>
    Oftsh< Ots<U> >& smprod(Oftsh< Ots<U> > const& a, Oftsh< Ots<U> > const& b, U const& m);
    template<typename U>
    Oftsh< Ofs<U> >& smprod(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m);


    template<typename U>
    Oftsh< Ots<U> >& smprod(Oftsh< Ots<U> > const& a, Oftsh< Ots<U> > const& b, Ots<U> const& r, U const& m, Ots<U> &temp);

    //Friendly streaming
    friend std::ostream& operator << <>(std::ostream& stream, Oftsh<T> const& oftsh);
};

// Functions
//---------------------------------------------------------------------------
template<typename T> bool operator==(Oftsh<T> const& a, Oftsh<T> const& b);
template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b);
template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b);


//Include the implementation .tpp
#include "oftsh.tpp"

#endif // OFTSH_H_INCLUDED
