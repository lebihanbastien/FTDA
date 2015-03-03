#ifndef OFS_H_INCLUDED
#define OFS_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ots.h"

using namespace std;

template<typename T>
class Ofs;

template<typename T>
std::ostream& operator << (std::ostream& stream, Ofs<T> const& ofs);

template<typename T>
class Ofs
{


private:
    int order;
    T *coef;
    static int globalOrder;

public:
    //Create
    Ofs<T>();
    Ofs<T>(const int newOrder);
    Ofs<T>(const int newOrder, T coef0);
    Ofs<T>(const int newOrder, T coefn, int n);


    //Copy
    Ofs<T>(Ofs const& ofs_);
    Ofs<T>& ccopy(Ofs<T> const& b);

    //Delete
    ~Ofs<T>();

    //Setters
    void setCoef(T value, int pos);
    void addCoef(T value, int pos);
    void setCoefs(T value);
    void conjugate();

    //Getters
    int getOrder() const;
    Ofs* getAddress() const;
    T getCoef(int pos) const;

    //Zeroing
    void zero();

    //Operators
    Ofs<T>& operator=(Ofs<T> const& b);
    Ofs<T>& operator  = (T const& coef0);
    Ofs<T>& operator += (Ofs<T> const& ofsToCopy);
    Ofs<T>& operator -= (Ofs<T> const& ofsToCopy);

    Ofs<T>& operator *= (T const& c);
    Ofs<T>& operator /= (T const& c);


    Ofs<T>& sprod(Ofs<T> const& a, Ofs<T> const& b);
    Ofs<T>&  prod(Ofs<T> const& a, Ofs<T> const& b);
    Ofs<T>& sprodn(Ofs<T> const& a, Ofs<T> const& b, int n);
    Ofs<T>& smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);
    Ofs<T>&  mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m);

    Ofs<T>& smult(Ofs<T> const& a, T const& c);
    Ofs<T>&  mult(Ofs<T> const& a, T const& c);

    bool isEqual(Ofs<T> const& b) const;

    //Ots 2 ofs
    Ofs<T>& ts2fs(Ots<T> const& ts);

    //Evaluate
    double complex evaluate(double const& t);


    //Ofstream friendly
    friend std::ostream& operator << <>(std::ostream& stream, Ofs<T> const& ofs);
};

//Functions

template<typename T> bool   operator == (Ofs<T> const& a, Ofs<T> const& b);
template<typename T> Ofs<T> operator - (Ofs<T> const& ofsToCopy);
template<typename T> Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b);
template<typename T> Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b);
template<typename T> Ofs<T> operator * (Ofs<T> const& a, T const& c);
template<typename T> Ofs<T> operator / (Ofs<T> const& a, T const& c);
template<typename T> Ofs<T> operator * (T const& c, Ofs<T> const& a);
template<typename T> Ofs<T> operator *  (Ofs<T> const& a, Ofs<T> const& b);

//Ofs to Ots
template<typename T> Ots<T>& fs2ts(Ots<T> & ts, Ofs<T> const& fs);

//Include the implementation .tpp
#include "ofs.tpp"

#endif // OFS_H_INCLUDED
