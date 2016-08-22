//############################################################################
// Implementation of the matrix template class
//############################################################################

/**
 * \file matrix.tpp
 * \brief Some extension of the vector class of C++ for matrix manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default creator for matrix class. Never use.
 **/
template <typename T> matrix<T>::matrix()
{
    size1 = 0;//NV;
    size2 = 0;//REDUCED_NV;
    coef = *(new vector<T>());
}

/**
 *  \brief Default creator for matrix class, with size size1 x size2.
 **/
template <typename T> matrix<T>::matrix(const int size1_, const int size2_)
{
    size1 = size1_;
    size2 = size2_;
    coef = *(new vector<T>(size1*size2));
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief Create a matrix equal to the matrix b. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>::matrix(matrix const& b)
{
    size1 = b.size1;
    size2 = b.size2;
    coef  = *(new vector<T>(size1*size2));
    for(int i = 0 ; i< size1*size2; i++) coef[i].ccopy(b.coef[i]);
}

/**
 *  \brief Create a matrix equal to the matrix b. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::operator = (matrix<T> const& b)
{
    if(this != &b)
    {
        size1 = b.size1;
        size2 = b.size2;
        coef  = *(new vector<T>(size1*size2));
        for(int i = 0 ; i < size1*size2; i++) coef[i].ccopy(b.coef[i]);
    }
    return *this; //same object if returned
}

/**
 *  \brief Copy the matrix b in this. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::ccopy(matrix<T> const& b)
{
    if(size1 != b.size1 || size2 != b.size2)
    {
        cout << "Erreur in ccopy for matrix: sizes do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {

    for(int i = 0 ; i < size1*size2; i++) coef[i].ccopy(b.coef[i]);
    return *this; //same object if returned
    }
}

/**
 *  \brief Link the matrix b with this. Needs the routine:
 *          - lcopy
 **/
template <typename T> matrix<T>& matrix<T>::lcopy(matrix<T> const& b)
{
    size1 = b.size1;
    size2 = b.size2;
    for(int i = 0 ; i < size1*size2; i++) coef[i].lcopy(b.coef[i]);
    return *this;
}

//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Delete function, empty.
 *         Do we need to implement something here? Probably not, not pointer used.
 **/
template <typename T> matrix<T>::~matrix<T>()
{
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief Gets the coefficient at the position (i,j)
 **/
template <typename T> T matrix<T>::getCoef(int i, int j) const
{
    if( i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0];
    }
    else return coef[i*size2 + j];
}

/**
 *  \brief Gets the address of the coefficient at the position (i,j)
 **/
template <typename T> T* matrix<T>::getCA(int i, int j) const
{
    if( i >= size1 || j >=size2)
    {
        cout << "Error in matrix<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0].getAddress();
    }
    else return coef[i*size2 + j].getAddress();
}

/**
 *  \brief Gets the size (either size1 or size2) of the matrix.
 **/
template <typename T> int matrix<T>::getSize(int num) const
{
    if(num == 1) return size1;
    else if(num == 2) return size2;
    else cout << "Error in matrix<T>::getSize: required number must be 1 or 2." << endl; return 0;

}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets the coefficient T at the position (i,j). Requires ccopy.
 **/
template <typename T> void matrix<T>::setCoef(T const & value, int i, int j)
{
    this->getCA(i,j)->ccopy(value);
}

/**
 *  \brief Sets the subcoefficient value at the order zero of the coefficient (i,j). Requires setCoef.
 **/
template <typename T> template <typename U> void matrix<T>::setCoef(U const & value, int i, int j)
{
    this->getCA(i,j)->setCoef(value, (int const) 0);
}

/**
 *  \brief Adds the coefficient T at the position (i,j). Requires smult.
 **/
template <typename T> void matrix<T>::addCoef(T const & value, int i, int j)
{
   this->getCoef(i,j).smult(value, 1.0);
}

/**
 *  \brief Zeroing of the matrix. Requires zero().
 **/
template <typename T> void matrix<T>::zero()
{
    for(int i = 0 ; i < size1*size2; i++) coef[i].zero();
}

//---------------------------------------------------------------------------
//Operations
//---------------------------------------------------------------------------
/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix. Requires der.
 **/
template <typename T> void matrix<T>::der(T const &a, int ni, int i, int j)
{
    if(i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else coef[i*size2 +j].der(a, ni);
}

/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix, at order k of the coefficients. Requires der (at order k).
 **/
template <typename T> void matrix<T>::der(T const &a, int ni, int i, int j, int k)
{
    if(i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else coef[i*size2 +j].der(a, ni, k);
}

/**
 *  \brief Derivation of a wrt to time, set at position (i,j) in the matrix. Frequency is n. Requires dot.
 **/
template <typename T> void matrix<T>::dot(T const &a, double n, int i, int j)
{
    if(i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else getCA(i,j)->dot(a, n);
}

/**
 *  \brief Derivation of a wrt to time of the whole matrix. Frequency is n. Requires dot.
 **/
template <typename T> void matrix<T>::dot(matrix<T> const &a, double n)
{
    if(a.size1 != size1 || a.size2 != size2)
    {
        cout << "Error in matrix<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else
    {
        for(int i = 0; i< size1; i++)
            for(int j = 0; j< size2; j++) this->dot(a.getCoef(i,j), n, i, j);
    }
}

//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >
//---------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with T = Ofts< Ofsc >. Requires ofts_sprod.
 **/
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.getSize(2) != vIn.size() || a.getSize(1) != vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (T*)&vOut[i])->ofts_sprod(a.getCoef(i,j), vIn[j]);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with T = Ofts< Ofsc >. Requires ofts_sprod.
 **/
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m)
{
    if( a.getSize(2) != (int) vIn.size() || a.getSize(1) != vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                vOut[i].ofts_sprod(*a.getCA(i,j), vIn[j], m);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >. Requires ofts_smult_t.
 **/
template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if( a.getSize(2) != (int) vIn.size() || a.getSize(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (T*)&vOut[i])->ofts_smult_t(vIn[j], (U&) a.coef[i*a.size2 + j]);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >. Requires ofts_smult_t.
 **/
template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m)
{
    if( a.getSize(2) != (int) vIn.size() || a.getSize(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (T*)&vOut[i])->ofts_smult_t(vIn[j], (U&) a.coef[i*a.size2 + j], m);
            }
        }
    }
}

/**
 *  \brief Add the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
 **/
template <typename T , typename U> void addCoef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->addCoef(a[i], 0, 0);
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
 **/
template <typename T , typename U> void subCoef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->addCoef(-a[i], 0, 0);
        }
    }
}


//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, TFS format
//---------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >, TFS format. Requires ofts_smult_t.
 **/
template <typename T , typename U> void tfts_smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m)
{
    if( a.getSize(2) != (int) vIn.size() || a.getSize(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (T*)&vOut[i])->tfts_smult_t(vIn[j], (U&) a.coef[i*a.size2 + j], m);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with T = Ofts< Ofsc >, TFS format. Requires ofts_sprod.
 **/
template <typename T>  void tfts_smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const &m)
{
    if( a.getSize(2) != (int) vIn.size() || a.getSize(1) != vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                vOut[i].tfts_sprod(*a.getCA(i,j), vIn[j], m);
            }
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
 *         Note the routine version for TFS format is identical to the one for OFS format because the underlying routines are used the same way.
 **/
template <typename T , typename U> void tfts_subCoef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->addCoef(-a[i], 0, 0);
        }
    }
}

/**
 *  \brief Add the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
 *         Note the routine version for TFS format is identical to the one for OFS format because the underlying routines are used the same way.
 **/
template <typename T , typename U> void tfts_addCoef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->addCoef(a[i], 0, 0);
        }
    }
}

/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix, at order k of the coefficients. Requires tfts_der (at order k).
 **/
template <typename T> void matrix<T>::tfts_der(T const &a, int ni, int i, int j, int k)
{
    if(i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else coef[i*size2 +j].tfts_der(a, ni, k);
}


//---------------------------------------------------------------------------
//TFS <--> OFS format
//---------------------------------------------------------------------------
/**
 *  \brief  Inline from frequency domain to time domain, for matrix< Ofsc > object.
 */
inline void tfs_from_ofs_inline(matrix< Ofsc >& a)
{
    Ofsc temp(a.getCA(0,0)->getOrder());
    for(int i =0; i < a.getSize(1) ; i++)
    {
        for(int j =0; j < a.getSize(2); j++)
        {
            a.getCA(i,j)->tfs_from_ofs_inline(temp);
        }
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for matrix< Ofsc > object.
 */
inline void tfs_to_ofs_inline(matrix< Ofsc >& a)
{
    for(int i =0; i < a.getSize(1) ; i++)
    {
        for(int j =0; j < a.getSize(2); j++)
        {
            a.getCA(i,j)->tfs_to_ofs_inline();
        }
    }
}

/**
 *  \brief  Inline from frequency domain to time domain, for vector< Ofsc > object.
 */
inline void tfs_from_ofs_inline(vector< Ofsc >& a)
{
    Ofsc temp(a[0].getOrder());
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_from_ofs_inline(temp);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofsc > object.
 */
inline void tfs_to_ofs_inline(vector< Ofsc >& a)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_to_ofs_inline();
    }
}

/**
 *  \brief  Inline from frequency domain to time domain, for vector< Ofts< Ofsc >  > object.
 */
inline void tfs_from_ofs_inline(vector< Ofts< Ofsc > >& a, int m)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_from_ofs_inline(m);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a, int m)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_to_ofs_inline(m);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a)
{
    for(int m = 0; m <= a[0].getOrder(); m++)
    {
        for(int i =0; i < (int) a.size() ; i++)
        {
            a[i].tfs_to_ofs_inline(m);
        }
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for matrix< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(matrix< Ofts< Ofsc > >& a)
{
    for(int i =0; i < a.getSize(1) ; i++)
    {
        for(int j =0; j < a.getSize(2); j++)
        {
            a.getCA(i,j)->tfs_to_ofs_inline();
        }
    }
}


//---------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, to use only in tests
//---------------------------------------------------------------------------
/**
 *  \brief Add the algebraic vector a to the order zero of vIn and set the result in vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
 *         Careful: this routine leads to the copy of an entire Ofts structure. Very heavy! Use only in tests.
 **/
template <typename T , typename U> void addCoef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.getSize() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < vOut.size() ; i++)
        {
                ((T*)&vOut[i])->ccopy(vIn[i]);
                ((T*)&vOut[i])->addCoef(a.getCoef(i), 0, 0);
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vIn and set the result in vOut. Used with T = Ofts< Ofsc >. Requires addCoef from vector class.
*         Careful: this routine leads to the copy of an entire Ofts structure. Very heavy! Use only in tests.
 **/
template <typename T , typename U> void subCoef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in addCoef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ((T*)&vOut[i])->ccopy(vIn[i]);
                ((T*)&vOut[i])->addCoef(-a[i], 0, 0);
        }
    }
}


//---------------------------------------------------------------------------
//Functions used with T = Ofs<U>
//---------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with T = Ofsc. Requires sprod.
 **/
inline void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut)
{
    if((unsigned int)  a.getSize(2) != vIn.size() || (unsigned int) a.getSize(1) != vOut.size() )
    {
        cout << "Error in smvprod_ofs (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (Ofsc*)&vOut[i])->ofs_sprod(a.getCoef(i,j), vIn[j]);
            }
        }
    }
}

//---------------------------------------------------------------------------
//Functions used with U = cdouble
//---------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product. Used with T = Ofs<U>, U = cdouble. Requires evaluate.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void mvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    //Zeroing the target
    for(int i = 0; i< (int) vOut.size(); i++) vOut[i] = 0+0.0*I;
    //Loop on coefficients
    for(int i =0; i < a.getSize(1) ; i++)
    {
        for(int j =0; j < a.getSize(2); j++)
        {
            vOut[i] += vIn[j]*a.getCA(i,j)->evaluate(t);
        }
    }
}

/**
 *  \brief Matrix-vector product (with sum). Used with T = Ofs<U>, U = cdouble. Requires evaluate.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    //Loop on coefficients
    for(int i =0; i < a.getSize(1) ; i++)
    {
        for(int j =0; j < a.getSize(2); j++)
        {
            vOut[i] += vIn[j]*a.getCA(i,j)->evaluate(t);
        }
    }
}

/**
 *  \brief Matrix-matrix product a x b. Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut)
{
    //Loop on coefficients
    for(int i =0; i < mOut.getSize(1) ; i++)
    {
        for(int j =0; j < mOut.getSize(2); j++)
        {
            for(int k =0; k < a.getSize(2); k++) mOut.getCA(i,j)->ofs_sprod(*a.getCA(i,k), *b.getCA(k,j));
        }
    }
}

/**
 *  \brief Matrix-matrix product a^T x b. Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smtmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut)
{
    //Loop on coefficients
    for(int i =0; i < mOut.getSize(1) ; i++)
    {
        for(int j =0; j < mOut.getSize(2); j++)
        {
            for(int k =0; k < a.getSize(2); k++) mOut.getCA(i,j)->ofs_sprod(*a.getCA(k,i), *b.getCA(k,j));
        }
    }
}

/**
 *  \brief vOut += a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsum_u(vector< Ofs<U> >& a, vector<U>& vOut, double const& t)
{
    for(int i =0; i < (int) vOut.size() ; i++)
    {
        vOut[i] += a[i].evaluate(t);
    }
}

/**
 *  \brief vOut -= a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsub_u(vector< Ofs<U> >& a, vector<U>& vOut, double const& t)
{
    for(int i =0; i < vOut.size() ; i++) vOut[i] -= a[i].evaluate(t);
}

/**
 *  \brief vOut = vIn - a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsub_u(vector< Ofs<U> >& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    for(int i =0; i < (int) vOut.size() ; i++) vOut[i] = vIn[i]-a[i].evaluate(t);
}

//---------------------------------------------------------------------------
// Read & Write
//---------------------------------------------------------------------------
/**
 * \brief Writes a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, size1(W)-1
 *        and j = 0, size2(W)-1.
 **/
inline void writeMOFTS_bin(matrix<Ofts<Ofsc > > &W, string filename)
{
    string ss1, ss2;
    //Loop on all coefficients
    for(int i = 0; i < W.getSize(1); i++)
    {
        for(int j = 0; j < W.getSize(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            writeOFTS_bin(*W.getCA(i,j), (filename+"["+ss1+"]"+"["+ss2+"].bin"));
        }
    }
}

/**
 * \brief Reads a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, size1(W)-1
 *        and j = 0, size2(W)-1.
 **/
inline void readMOFTS_bin(matrix<Ofts<Ofsc > > &W, string filename, int fftN)
{
    string ss1, ss2;
    //Loop on all coefficients
    for(int i = 0; i < W.getSize(1); i++)
    {
        for(int j = 0; j < W.getSize(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            readOFTS_bin(*W.getCA(i,j), (filename+"["+ss1+"]["+ss2+"].bin"), fftN);
        }
    }
}

/**
 *  \brief writes a matrix of Ofts objects from files. DEPRECATED.
 **/
inline void writeMOFTS_txt(matrix<Ofts<Ofsc > > &W, string filename)
{
    ifstream readStream;
    ofstream myfile;
    string ss1, ss2;
    //Loop on all coefficients

    for(int i = 0; i < W.getSize(1); i++)
    {
        for(int j = 0; j < W.getSize(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            writeOFTS_txt(*W.getCA(i,j), (filename+"["+ss1+"]["+ss2+"].txt"));
        }
    }
}

/**
 *  \brief read a matrix of Ofts objects from files. DEPRECATED.
 **/
inline void readMOFTS(matrix<Ofts<Ofsc > > &W, string filename, int fftN)
{
    ifstream readStream;
    ofstream myfile;
    string ss1, ss2;
    //Loop on all coefficients

    for(int i = 0; i < W.getSize(1); i++)
    {
        for(int j = 0; j < W.getSize(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            readOFTS_txt(*W.getCA(i,j), (filename+"["+ss1+"]["+ss2+"].txt"), fftN);
        }
    }
}



