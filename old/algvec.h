#ifndef ALVECTOR_H_INCLUDED
#define ALVECTOR_H_INCLUDED

/**
 * \file vector.h
 * \brief Some extension of the vector class of C++ for algebraic manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ofts.h"
#include "parameters.h"

using namespace std;

//--------------------------------------------------------------------------
// vector<T> class
//--------------------------------------------------------------------------
template<typename T>
class vector;


/** \class vector
 *  \brief Some extension of the vector class of C++ for algebraic manipulation
 *
 *  The vector class is a (very) small extension of the default C++ vector class in order to include some new operations, such as derivation wrt to a given variable
 *  defined within the coefficients of the vector (time for the vector<Ofs>, another variable for the vector<Ofts>).
 *
 */
template <typename T>
class vector   // : public vector<T> heritage from vector: a good idea? a useful one?
{

private:

    int size1;      ///size of the vector
    vector<T> coef; ///array coefficients

public:
    //------------------
    //Create
    //------------------
    /**
     *  \brief Default constructor of the class vector<T>.
     */
    vector<T>();
    /**
     *  \brief Constructor with a given size
     */
    vector<T>(const int size1_);
    /**
     *  \brief Constructor from a given vector object (without any link).
     */
    vector<T>(vector const& ofs_);

    //------------------
    //Copy
    //------------------
    /**
     *  \brief  Copy from a given vector object (only the coefficients).
     */
    vector<T>& ccopy(vector<T> const& b);
    /**
     *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
     */
    vector<T>& lcopy(vector<T> const& b);
    /**
     *  \brief  An operator. Constructor from a given vector object (only the coefficients).
     */
    vector<T>& operator  = (vector<T> const& b);

    //------------------
    //Delete
    //------------------
    /**
     *  \brief Default destructor of the class vector<T>.
     */
    ~vector<T>();

    //------------------
    //Setters
    //------------------
    /**
     *  \brief Sets a coefficient at a given position in the vector
     */
    void setCoef(T const & value, int i);
    /**
     *  \brief Adds a coefficient at a given position in the serie.
     */
    void addCoef(T const & value, int i);

    //------------------
    //Getters
    //------------------
    /**
     *  \brief  Gets the coefficient at a given position.
     *
     *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
     */
    T getCoef(int i) const;
    /**
     *  \brief  Gets the address of the ith coefficient.
     */
    T* getCA(int i) const;
    /**
     *  \brief  Gets the size of the vector.
     */
    int getSize() const;


    //------------------
    //Operations
    //------------------
    /**
     *  \brief Derivation of the ith coefficient wrt to the variable ni.
     */
    void der(T const &a, int ni, int i);
    /**
     *  \brief Derivation of a the ith coefficient wrt to time.
     */
    void dot(T const &a, double n, int i);
    /**
     *  \brief Derivation wrt to time.
     */
    void dot(vector<T> const &a, double n);
};


//Include the implementation .tpp
#include "vector.tpp"

//--------------------------------------------------------------------------
// end of vector<T> class
//--------------------------------------------------------------------------

#endif // ALVECTOR_H_INCLUDED
