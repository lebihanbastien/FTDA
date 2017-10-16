#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ofts.h"


/**
 * \file matrix.h
 * \brief An extension of the vector class of C++ for matrix manipulation. See the .tpp file for comments.
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//--------------------------------------------------------------------------
// matrix<T> class
//--------------------------------------------------------------------------
template<typename T>
class matrix;



template <typename T>
class matrix
{

private:

public:

    int size1;      //rows
    int size2;      //columns
    vector<T> coef; //coefficients

    //--------------------------------------------
    //Create
    //--------------------------------------------
    matrix<T>();
    matrix<T>(const int size1_, const int size2_);

    //--------------------------------------------
    //Copy
    //--------------------------------------------
    matrix<T>(matrix const& ofs_);
    matrix<T>& ccopy(matrix<T> const& b);
    matrix<T>& lcopy(matrix<T> const& b);

    //--------------------------------------------
    //Delete
    //--------------------------------------------
    ~matrix<T>();

    //--------------------------------------------
    //Setters
    //--------------------------------------------
    void setCoef(T const & value, int i, int j);
    void addCoef(T const & value, int i, int j);
    template <typename U> void setCoef(U const & value, int i, int j);
    void setCoef(cdouble const & value, int i, int j);

    //--------------------------------------------
    //Getters
    //--------------------------------------------
    T getCoef(int i, int j) const;
    T* getCA(int i, int j) const;
    int getSize(int num) const;

    //--------------------------------------------
    //Operators
    //--------------------------------------------
    matrix<T>& operator  = (matrix<T> const& b);

    //--------------------------------------------
    //Operations
    //--------------------------------------------
    void der(T const &a, int ni, int i, int j);
    void der(T const &a, int ni, int i, int j, int k);
    void dot(T const &a, double n, int i, int j);
    void dot(matrix<T> const &a, double n);
    void zero();
    void tfts_der(T const &a, int ni, int i, int j, int k);


    T operator [](int i) const    {return coef[i];}
    T & operator [](int i) {return coef[i];}
    T operator ()(int i, int j) const    {return coef[size2*i+j];}
    T & operator ()(int i, int j) {return coef[size2*i+j];}

};

//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >
//---------------------------------------------------------------------------
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut);
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const& m);

template <typename T , typename U> void addCoef(vector<U> const& a, vector<T>& vOut);
template <typename T , typename U> void subCoef(vector<U> const& a, vector<T>& vOut);

template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut);
template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m);

//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, TFS format
//---------------------------------------------------------------------------
template <typename T , typename U> void tfts_smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m);
template <typename T>  void tfts_smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m);
template <typename T , typename U> void tfts_subCoef(vector<U> const& a, vector<T>& vOut);

//---------------------------------------------------------------------------
//TFS <--> OFS format
//---------------------------------------------------------------------------
inline void tfs_from_ofs_inline(matrix< Ofsc >& a);
inline void tfs_to_ofs_inline(matrix< Ofsc >& a);
inline void tfs_from_ofs_inline(vector< Ofts< Ofsc > >& a, int m);
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a, int m);
inline void tfs_to_ofs_inline(matrix< Ofts< Ofsc > >& a);

//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, to use only in tests
//---------------------------------------------------------------------------
template <typename T , typename U> void addCoef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut);
template <typename T , typename U> void subCoef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut);


//---------------------------------------------------------------------------
//Functions used with T = Ofs<U>
//---------------------------------------------------------------------------
void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut);

//---------------------------------------------------------------------------
//Functions used with U = cdouble
//---------------------------------------------------------------------------
template <typename U> void mvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void smvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void vvsum_u(vector< Ofs<U> > & a, vector<U>& vOut, double const& t);
template <typename U> void vvsub_u(vector< Ofs<U> > & a, vector<U>& vOut, double const& t);
template <typename U> void vvsub_u(vector< Ofs<U> > & a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void smmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut);
template <typename U> void smtmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut);

//---------------------------------------------------------------------------
// Read & Write
//---------------------------------------------------------------------------
/**
 * \brief Writes a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, size1(W)-1
 *        and j = 0, size2(W)-1.
 **/
void writeMOFTS_bin(matrix<Ofts<Ofsc > > &W, string filename);

/**
 * \brief Reads a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, size1(W)-1
 *        and j = 0, size2(W)-1.
 **/
inline void readMOFTS_bin(matrix<Ofts<Ofsc > > &W, string filename);

//DEPRECATED
void readMOFTS(matrix<Ofts<Ofsc > > &W, string filename, int fftN);
void writeMOFTS_txt(matrix<Ofts<Ofsc > > &W, string filename, int fftN);


//---------------------------------------------------------------------------
//Include the implementation .tpp
//---------------------------------------------------------------------------
#include "matrix.tpp"

//--------------------------------------------------------------------------
// end of matrix<T> class
//--------------------------------------------------------------------------

#endif // MATRIX_H_INCLUDED
