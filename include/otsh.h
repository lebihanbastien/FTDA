#ifndef OTSH_H_INCLUDED
#define OTSH_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <vector>


#include "ftda.h"

using namespace std;


template <typename T>
class Otsh;

template<typename T>
std::ostream& operator << (std::ostream& stream, Otsh<T> const& otsh);


template <typename T>
class Otsh
{

private:

    int nv;             //number of variables
    int order;          //degree
    T   *coef;          //address of the first coefficient
    Otsh<T> *term;      //recursivity

    static int globalOrder;
    static int globalVariable;

public:

    //Create
    Otsh<T>();
    Otsh<T>(int newNv, int newOrder);
    Otsh<T>(int newNv, int newOrder, T *coef0);

    //Delete
    ~Otsh<T>();


    //Copy
    Otsh<T>(Otsh<T> const& b);                //unlinked copy
    Otsh<T>& operator = (Otsh<T> const& b);   //unlinked copy
    Otsh<T>& lcopy(Otsh<T> const& b);         //linked copy (for recursivity)
    Otsh<T>& ccopy(Otsh<T> const& b);         //coefficient copy (restricted to same order, same number of variables)


    //Linking
    void setCoefs(T *coef0);

    //Zeroing
    void zero();

    //Setters
    void setCoef(T value, int pos);
    void setCoefTot(T value, int pos);
    void addCoef(T value, int pos);
    void conjugate();
    void setRandomCoefs();

    //Getters
    Otsh<T> getTerm();
    Otsh<T> getTerm(int i);
    T getCoef(int i);
    int getOrder();
    int getVariables();
    T* getCoefAddress();

    //Operators
    bool isEqual(Otsh<T> const& b) const;
    Otsh<T>& operator += (Otsh<T> const& b);
    Otsh<T>& operator -= (Otsh<T> const& b);
    Otsh<T>& operator *= (T const& c);
    Otsh<T>& operator /= (T const& c);


    Otsh<T>& smprod(Otsh<T> const& a, Otsh<T> const& b, T const& m);
    Otsh<T>&  sprod(Otsh<T> const& a, Otsh<T> const& b);

    Otsh<T>& smult(Otsh<T> const& a, T const& c);
    Otsh<T>&  mult(Otsh<T> const& a, T const& c);
    Otsh<T>&  mult(T const& c);


    //Friendly streaming
    friend std::ostream& operator << <>(std::ostream& stream, Otsh<T> const& otsh);
};

// Functions
//---------------------------------------------------------------------------
template<typename T> bool operator==(Otsh<T> const& a, Otsh<T> const& b);
template<typename T> Otsh<T> operator + (Otsh<T> const& a, Otsh<T> const& b);
template<typename T> Otsh<T> operator - (Otsh<T> const& a, Otsh<T> const& b);


//Include the implementation .tpp
#include "otsh.tpp"


#endif // OTSH_H_INCLUDED
