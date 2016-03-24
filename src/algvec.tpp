//############################################################################
// Implementation of the vector template class
//############################################################################

/**
 * \file vector.tpp
 * \brief Some extension of the vector class of C++ for algebraic manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class vector<T>.
 */
template <typename T> vector<T>::vector()
{
    size1 = NV;
    coef = vector<T>(size1);
}

/**
 *  \brief Constructor with a given size
 */
template <typename T> vector<T>::vector(const int size1_)
{
    size1 = size1_;
    coef = vector<T>(size1);
}

/**
 *  \brief Constructor from a given vector object (without any link).
 */
template <typename T> vector<T>::vector(vector const& b)
{
    size1 = b.size1;
    coef  = vector<T>(size1);
    for(int i = 0 ; i< size1; i++) coef[i].ccopy(b.coef[i]);
}


//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
//Do we need to implement something here?
/**
 *  \brief Default destructor of the class vector<T>.
 */
template <typename T> vector<T>::~vector<T>()
{
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief  An operator. Constructor from a given vector object (only the coefficients).
 */
template <typename T> vector<T>& vector<T>::operator = (vector<T> const& b)
{
    if(this != &b)
    {
        size1 = b.size1;
        coef  = vector<T>(size1);
        for(int i = 0 ; i < size1; i++) coef[i].ccopy(b.coef[i]);
    }
    return *this; //same object if returned
}

/**
 *  \brief  Copy from a given vector object (only the coefficients).
 */
template <typename T> vector<T>& vector<T>::ccopy(vector<T> const& b)
{
    if(size1 != b.size1)
    {
        cout << "Erreur in ccopy for vector: sizes do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {

    for(int i = 0 ; i < size1; i++) coef[i].ccopy(b.coef[i]);
    return *this; //same object if returned
    }
}

/**
 *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
 */
template <typename T> vector<T>& vector<T>::lcopy(vector<T> const& b)
{
    size1 = b.size1;
    for(int i = 0 ; i < size1; i++) coef[i].lcopy(b.coef[i]);
    return *this;
}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the vector
 */
template <typename T> void vector<T>::setCoef(T const & value, int i)
{
    this->getCA(i)->ccopy(value);
}

/**
 *  \brief Adds a coefficient at a given position in the serie.
 */
template <typename T> void vector<T>::addCoef(T const & value, int i)
{
   this->getCA(i)->smult(value, 1.0);
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <typename T> T vector<T>::getCoef(int i) const
{
    if( i >= size1)
    {
        cout << "Error in vector<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0];
    }
    else return coef[i];
}

/**
 *  \brief  Gets the address of the ith coefficient.
 */
template <typename T> T* vector<T>::getCA(int i) const
{
    if( i >= size1)
    {
        cout << "Error in vector<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0].getAddress();
    }
    else return coef[i].getAddress();
}

/**
 *  \brief  Gets the size of the vector.
 */
template <typename T> int vector<T>::getSize() const
{
    return size1;
}


//---------------------------------------------------------------------------
//Operations
//---------------------------------------------------------------------------
/**
 *  \brief Derivation wrt to time.
 */
template <typename T> void vector<T>::dot(vector<T> const &a, double n)
{
    if(a.size1 != size1)
    {
        cout << "Error in vector<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else
    {
        for(int i = 0; i< size1; i++) this->dot(a.getCoef(i), n, i);
    }
}

/**
 *  \brief Derivation of a the ith coefficient wrt to time.
 *
 *   The class T of the coefficient must contain the derivation function dot(U const &a, double n), as in the Ofs<T> class.
 */
template <typename T> void vector<T>::dot(T const &a, double n, int i)
{
    if(i >= size1)
    {
        cout << "Error in vector<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else getCA(i)->dot(a, n);
}

/**
 *  \brief Derivation of the ith coefficient wrt to the variable ni.
 *
 *   The class T of the coefficient must contain the derivation function der(U const &a, double ni), as in the Ofts<T> class.
 */
template <typename T> void vector<T>::der(T const &a, int ni, int i)
{
    if(i > size1)
    {
        cout << "Error in vector<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else coef[i].der(a, ni);
}
