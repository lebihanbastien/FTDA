//############################################################################
// Implementation of the Ofts template class
//############################################################################

/**
 * \file ofts.tpp
 * \brief Fourier-Taylor series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofts<T>.
 */
template<typename T> Ofts<T>::Ofts()
{
    int i, index;
    nv     = REDUCED_NV;
    order  = OFTS_ORDER;
    cnv    = OFS_NV;
    corder = OFS_ORDER;

    //New array of coefficients
    coefs = (Ofs<complex double>*) calloc(binomial(nv+order,nv), sizeof(Ofs<complex double>)); //new T[binomial(nv+order, nv)]();
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k] = Ofs<complex double>(corder);

    //Allocation of the homogeneous polynomials
    term = new Oftsh<T>*[order+1];
    //term = (Oftsh<T>**) calloc(binomial(nv+order,nv), sizeof(Oftsh<T>*));//new Oftsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh<T>(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->linkCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
 */
template<typename T> Ofts<T>::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    coefs = (T*) calloc(binomial(nv+order,nv), sizeof(T)); //new T[binomial(nv+order, nv)]();

    //Allocation of the homogeneous polynomials
    term = new Oftsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->linkCoefs(coefs+index);
        //term[i]->linkCoefs(new T[FTDA::nmon(nv,i)]);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients, Ofsc case.
 */
template<> inline Ofts< Ofs<complex double> >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    coefs = (Ofs<complex double>*) calloc(binomial(nv+order,nv), sizeof(Ofs<complex double>)); //new T[binomial(nv+order, nv)]();
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ofs<complex double>(corder);

    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ofs<complex double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh< Ofs<complex double> >(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->linkCoefs(coefs+index);
        //term[i]->linkCoefs(new T[FTDA::nmon(nv,i)]);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients. Ofts< Ofsd > case.
 */
template<> inline Ofts< Ofsd >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    coefs = (Ofsd*) calloc(binomial(nv+order,nv), sizeof(Ofsd));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ofsd(corder);
    //-------------------------------------------------------

    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ofsd >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh< Ofsd >(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->linkCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor from a given Ofts object (without any link).
 */
template<typename T> Ofts<T>::Ofts(Ofts<T> const& b)
{
    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    //coefs = new T[binomial(nv+b.order, b.nv)]();
    coefs = (Ofs<complex double>*) calloc(binomial(nv+order,nv), sizeof(Ofs<complex double>)); //new T[binomial(nv+order, nv)]();
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ofs<complex double>(corder);

    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coefs[index+i] = b.term[nrc]->getCoef(i);
        index+=  FTDA::nmon(nv,nrc);
    }

    //Allocation of the homogeneous polynomials
    term = new Oftsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=order; nrc++)
    {
        //Allocation of each hp
        term[nrc] = new Oftsh<T>(nv, nrc);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[nrc]->linkCoefs(coefs+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

/**
 *  \brief  An operator. Constructor from a given Ofts object (only the coefficients).
 */
template<typename T> Ofts<T>& Ofts<T>::operator = (Ofts<T> const& b)
{
    if(this != &b)
    {
        int nrc, index;

        //Same nv/order
        nv = b.nv;
        order = b.order;
        cnv = b.cnv;
        corder = b.corder;

        //Copy of all the coefficients at every order in new array
        T *coef0 = new T[binomial(nv+order, nv)]();
        index = 0;
        for(int nrc=0; nrc<= order; nrc++)
        {
            for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i] = b.term[nrc]->getCoef(i);
            index+=  FTDA::nmon(nv,nrc);
        }

        delete term;

        //Allocation of the homogeneous polynomials
        term = new Oftsh<T>*[order+1];

        //Allocation of the coefficients
        index = 0;
        for(nrc=0; nrc<=order; nrc++)
        {
            //Allocation of each hp
            term[nrc] = new Oftsh<T>(nv, nrc);//allocate_homog(nv, i);
            //Link h to coefs at each level of the tree
            term[nrc]->linkCoefs(coef0+index);
            index+=  FTDA::nmon(nv,nrc);
        }

    }

    return *this;
}


//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofts<T>. WARNING: memory leak here, through the terms of type Oftsh.
 *
 * Certainly a problem here: the line thats delete the coefficient is commented, because it leads to an error when programs end.
 * May lead to memory leak if the objects are created "on the fly", which may be the case in some inner routines like smprod_t.
 */
template<typename T> Ofts<T>::~Ofts<T>()
{
    //if(coefs != NULL) delete coefs;
    if(term != NULL)
    {
        //Certainly a problem at this point: since the delete routine of Oftsh is empty, only the first leaf of the Oftsh tree is deleted...
        for(int i =0; i<= order ; i++) delete term[i];
    }
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
 */
template<typename T> Ofts<T>& Ofts<T>::lcopy (Ofts<T> const& b)
{
    order = b.order;
    nv = b.nv;
    term = b.term;
    cnv = b.cnv;
    corder = b.corder;
    return *this;
}

/**
 *  \brief  Copy from a given Ofts object (only the coefficients).
 *
 *  Restricted to same order, same number of variables
 */
template<typename T> Ofts<T>& Ofts<T>::ccopy (Ofts<T> const& b)
{
    if(order != b.order || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(b.term[nrc]->getCoef(i), i);
        return *this;
    }
}

/**
 *  \brief  Copy from a given Ofts object (only the coefficients) at order nrc.
 *
 *  Restricted to same order, same number of variables
 */
template<typename T> Ofts<T>& Ofts<T>::ccopy (Ofts<T> const& b, int const& nrc)
{
    if(nrc > min(order, b.order) || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(b.term[nrc]->getCoef(i), i);
        return *this;
    }
}


//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given order \c ord and a given position \c i at this order in the serie.
 *
 * Set of given coefficient at term ord and position i in this term
 *   - ord gives the order of the homogeneous polynomial (hp)
 *   - i gives the positions wihtin this hp
 * Ex: setCoef(m, 2, 0) sets the x^2 coef
 * Ex: setCoef(m, 2, 1) sets the x*y coef
 */
template<typename T> void Ofts<T>::setCoef(T const& m, int ord, int i)
{
    term[ord]->setCoef(m, i);
}

/**
 *  \brief Adds a coefficient at a given order \c ord and a given position \c i at this order in the serie.
 *
 * Set of given coefficient at term ord and position i in this term
 *   - ord gives the order of the homogeneous polynomial (hp)
 *   - i gives the positions wihtin this hp
 * Ex: addCoef(m, 2, 0) adds to the x^2 coef
 * Ex: addCoef(m, 2, 1) adds to the x*y coef
 */
template<typename T> void Ofts<T>::addCoef(T const& m, int ord, int i)
{
    term[ord]->addCoef(m, i);
}

/**
 *  \brief Sets of a given double/complex (typename U) subcoefficient everywhere (at each order and in each coefficient).
 */
template<typename T> template < typename U > void Ofts<T>::setAllCoefs(U const& m)
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setSubCoef(m, i);
}

/**
 *  \brief Sets random coefficients to all positions in the serie.
 */
template<typename T> void Ofts<T>::setRandomCoefs()
{
    for(int nrc=0; nrc<= order/2; nrc++)
    {
        term[nrc]->setRandomCoefs();
    }
}

/**
 *  \brief Sets of a given U subcoefficient at order zero of the coefficient at position \c i of term of order \c n.
 */
template<typename T> template < typename U > void Ofts<T>::setCoef0(U const& m, int const& ord, int const& i)
{
    //term[ord]->setT0Coef(m, i);
    term[ord]->setSubCoef(m, i, 0);
}


//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the serie.
 */
template<typename T> int Ofts<T>::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the order of the coefficients.
 */
template<typename T> int Ofts<T>::getCOrder() const
{
    return corder;
}

/**
 *  \brief  Gets the number of variables of serie.
 */
template<typename T> int Ofts<T>::getNV() const
{
    return nv;
}

/**
 *  \brief  Gets the number of variables of the coefficients.
 */
template<typename T> int Ofts<T>::getCVariables() const
{
    return cnv;
}

/**
 *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
 */
template<typename T>  T* Ofts<T>::getCoef(int const& ord, int const& pos) const
{
    if(ord > order || pos >= FTDA::nmon(nv, ord))
    {
        cout << "Error in getCoef: out of range. First term is returned" << endl;
        cout << "Requested order: " << ord << ", Maximum allowed: " <<  order << endl;
        cout << "Requested pos: " << pos << ", Maximum allowed: " <<  FTDA::nmon(nv, ord) << endl;
        return this->term[0]->getCA();
    }
    else return this->term[ord]->getCA()+pos;
}

/**
 *  \brief  Gets the adress of the term at order \c ord
 */
template<typename T> Oftsh<T>* Ofts<T>::getTerm(int const& ord) const
{
    if(ord > order)
    {
        cout << "Error in getTerm: out of range. First term is returned" << endl;
        cout << "Requested order: " << ord << ", Maximum allowed: " <<  order << endl;
        return this->term[0];
    }
    else return this->term[ord];
}

/**
 *  \brief  Gets the adress of the Ofts object
 */
template <typename T> Ofts<T>* Ofts<T>::getAddress() const
{
    return (Ofts<T>*) this;
}



//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template<typename T> void Ofts<T>::zero()
{
    for(int nrc=0; nrc<= order; nrc++) term[nrc]->zero();
}


//--------------------------------------------------------------------------------
//Operations
//--------------------------------------------------------------------------------
//------------------
// Conjugate
//------------------
/**
 *  \brief Conjugates  all terms (Oftsh object), and only them! To be used with evaluate_conjugate to have the true conjugate.
 */
template<typename T> Ofts<T>& Ofts<T>::conjugate()
{
    for(int nrc=0; nrc<= order; nrc++) term[nrc]->conjugate();
    return *this;
}

/**
 *  \brief Conjugates the order \c nrc. To be used with evaluate_conjugate to have the true conjugate.
 */
template<typename T> Ofts<T>& Ofts<T>::conjugate(int const& nrc)
{
    term[nrc]->conjugate();
    return *this;
}

//------------------
// smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient.
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= order; k++) term[k]->oftsh_smult_t(*a.term[k], m);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->oftsh_smult_t(*a.term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}


/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_mult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->oftsh_mult_t(*a.term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}

//---------------------------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= order; k++) term[k]->oftsh_smult_u(*a.term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k)
{
    if(k > order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->oftsh_smult_u(*a.term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

//---------------------------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->oftsh_smult_tu(*a.term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient, at order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c, int const& k)
{
    if(k > order  || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->oftsh_smult_tu(*a.term[k], m, c);  //ps[k] += m*a[k]
        return *this;
    }
}


//------------------
// mult
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  += m a \f$ with m a coefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_mult_t(Ofts<T> const& a, T const& m)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= order; k++)
        {
            term[k]->mult(*a.term[k], m);  //ps[k] = m*a[k], contains the zeroing
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient at order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k)
{
    if(k > order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->oftsh_mult_u(*a.term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->oftsh_mult_u(*a.term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c m a \f$ with m a coefficient and c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->mult(*a.term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}


//------------------
// sfsum
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= order; k++)
        {
            term[k]->oftsh_smult_t(*a.term[k], ma);   //ps[k]  += ma*a[k]
            term[k]->oftsh_smult_t(*b.term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->oftsh_smult_t(*a.term[n], ma);   //ps[k]  += ma*a[k]
        term[n]->oftsh_smult_t(*b.term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  c: a reference to an Ofts object
 *  \param  mc: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_tt(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv|| order != c.order || nv != c.nv)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->oftsh_smult_t(*a.term[n], ma);   //ps[k]  += ma*a[k]
        term[n]->oftsh_smult_t(*b.term[n], mb);   //ps[k]  += mb*b[k]
        term[n]->oftsh_smult_t(*c.term[n], mc);   //ps[k]  += mb*b[k]
        return *this;
    }

}

//------------------
// fsum
//------------------

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_fsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->oftsh_mult_t(*a.term[k], ma);   //ps[k]   = ma*a[k], contains the zeroing
            term[k]->oftsh_smult_t(*b.term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_fsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->oftsh_mult_t(*a.term[n], ma);    //ps[k]   = ma*a[k], contains the zeroing
        term[n]->oftsh_smult_t(*b.term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= order; k++)
        {
            term[k]->oftsh_mult_u(*a.term[k], ca);        //ps[k]   = ca*a[k], contains the zeroing
            term[k]->oftsh_smult_u(*b.term[k], cb);       //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 *  \param  m: a reference to the order to update
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb, int const& m)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[m]->oftsh_mult_u(*a.term[m], ca);        //ps[k]   = ca*a[k], contains the zeroing
        term[m]->oftsh_smult_u(*b.term[m], cb);       //ps[k]  += mb*b[k]
        return *this;
    }

}


//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$.
 *
 *  Handle the case for which n >= max(a.order, b.order)
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b)
{
    int k, i, i0, i1;
    //Product
    for(k=0; k<=order; k++)
    {
        i0 = min(b.order, k);
        i1 = min(a.order, k);
        for(i= k-i0; i<=i1; i++) term[k]->oftsh_sprod(*a.term[i], *b.term[k-i]);
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.order, b.order)
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i, i0, i1;
    i0 = min(b.order, n);
    i1 = min(a.order, n);
    //Product
    for(i= n-i0; i<=i1; i++) term[n]->oftsh_sprod(*a.term[i], *b.term[n-i]);
    return *this;
}

//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, T& temp)
{
//    //temp = 0 sizeOf(this)
//    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
//    //temp = a*b
//    temp.ofts_sprod(a,b);
//    //this = m*temp
//    this->ofts_smult_t(temp, m);
//    return *this;

//-----------------------------------------------------------------------
//Other possibility: use a temporary variable in Ofs, and not Ofts
    int k, i, i0, i1;
    //Product
    for(k=0; k<=order; k++)
    {
        i0 = min(b.order, k);
        i1 = min(a.order, k);
        for(i= k-i0; i<=i1; i++) term[k]->oftsh_smprod_t(*a.term[i], *b.term[k-i], m, temp);
    }
    return *this;
//-----------------------------------------------------------------------
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, int const& n, T& temp)
{
//    //temp = 0 sizeOf(this)
//    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
//    //temp = a*b at order n
//    temp.ofts_sprod(a,b,n);
//    //this = m*temp at order n
//    this->ofts_smult_t(temp,m, n);


    //-----------------------------------------------------------------------
    //Other possibility: use a temporary variable in Ofs, and not Ofts
    int i, i0, i1;

        i0 = min(b.order, n);
        i1 = min(a.order, n);
        for(i= n-i0; i<=i1; i++) term[n]->oftsh_smprod_t(*a.term[i], *b.term[n-i], m, temp);

    return *this;
    //-----------------------------------------------------------------------


    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c*a*b \f$ with c a subcoefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
    //temp = a*b
    temp.ofts_sprod(a,b);
    //this = c*temp
    this->ofts_smult_u(temp,c);
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 *
 * NOTE: we test smprod at all order. If the test is good, then we change smprod at order n
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
    //temp = a*b at order n
    temp.ofts_sprod(a,b,n);
    //this = c*temp at order n
    this->ofts_smult_u(temp,c, n);
    return *this;
}


//------------------
// prod
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with m a coefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  m: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, T const& m)
{
    this->zero();
    this->smprod(a,b,m);
    return *this;
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with c a subcoefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    this->zero();
    this->smprod(a,b,c);
    return *this;
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = a*b \f$.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 */
template<typename T> Ofts<T>& Ofts<T>::prod(Ofts<T> const& a, Ofts<T> const& b)
{
    this->zero();
    this->sprod(a,b);
    return *this;
}


//------------------
// pows
//------------------
/**
 *   \brief Power function: p = a^alpha at order n. This generic routine is not directly used, only specialized versions are (in particular for the Ofts<Ofs> form).
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_pows(Ofts<T> const& a,  U const& alpha, int const& n)
{
    T x0;
    x0 = coef0s(a);
    if(n==0)
    {
        this->acoef0s(cpow(x0, alpha));
    }
    else
    {
        //Sets every coefficients to zero for order n
        this->term[n]->zero();
        for(int j=0; j<= n-1; j++) this->term[n]->smprod(*a.term[n-j], *this->term[j], alpha*(n-j)-j);// smprodh(ps->term[k], s->term[k-j], ps->term[j], alpha*(k-j)-j);
        this->term[n]->mult(1.0/(n*x0)); //multh(ps->term[k], 1.0/(x0*k));
    }

    return *this;
}

/**
 *   \brief Power function: p = a^alpha at order n. Ofsc case: ONLY FOR OFS_ORDER = 0
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::ofts_pows(Ofts< Ofsc > const& a,  U const& alpha, int const& n)
{
    if(n==0)
    {
        Ofsc temp(corder);
        Ofsc temp2(corder);
        //Initialization of order zero
        temp.ofs_pows(*coef0s(a), alpha);
        this->acoef0s(temp);
    }
    else
    {
        Ofsc temp(corder);
        Ofsc temp2(corder);
        Ofsc temp3(corder);
        //temp2 = 1/x0
        temp2.ofs_pows(*coef0s(a), -1.0+0.0*I);
        //Sets every coefficients to zero for order n
        this->term[n]->zero();
        //Recurrence scheme @order n
        int i0 = min(a.order, n);
        for(int j= n-i0; j< n; j++)
        {
            temp.ofs_mult(temp2, (alpha*(n-j)-j)/n);
            this->term[n]->oftsh_smprod_t(*a.term[n-j], *this->term[j], temp, temp3);
        }
    }

    return *this;
}

//------------------
// Order 0 routines
//------------------
/**
 * \brief Returns the address of the first coefficient of order 0 of the taylor serie s
 */
template<typename T> T* Ofts<T>::coef0s(Ofts<T> const& a)
{
    return a.term[0][0].getCA();
}

/**
 * \brief Sets the coefficient of order 0 of the taylor serie s equal to x0
 */
template<typename T> void Ofts<T>::acoef0s(T const& x0)
{
    this->term[0][0].setCoef(x0,0);
}

//--------------------------------------------------------------------------------
//Operations with TFS coefficients - pure operations
//--------------------------------------------------------------------------------
//----------------
// Frequency to Time domain
//---------------
/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs(Ofts<T> const& a)
{
    Ofs< cdouble > tfs(corder);
    for(int nrc=0; nrc<= order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)
        {
            tfs.tfs_from_ofs(a.term[nrc]->getCoef(i));
            this->term[nrc]->setCoef(tfs, i);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs_inline()
{
    Ofs< cdouble > tfs(corder);
    for(int nrc=0; nrc<= order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)
        {
            (this->term[nrc]->getCA()+i)->tfs_from_ofs_inline(tfs);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs_inline(int nrc)
{
    Ofs< cdouble > tfs(corder);
    //Current homogeneous polynomial
    for (int i=0; i< FTDA::nmon(nv, nrc); i++)
    {
        (this->term[nrc]->getCA()+i)->tfs_from_ofs_inline(tfs);
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this.
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs(Ofts<T> const& a)
{
    Ofs< cdouble > ofs(corder);
    for(int nrc=0; nrc<= order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)
        {
            ofs.tfs_to_ofs(a.term[nrc]->getCoef(i));
            this->term[nrc]->setCoef(ofs, i);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs_inline()
{
    for(int nrc=0; nrc<= order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)
        {
            (this->term[nrc]->getCA()+i)->tfs_to_ofs_inline();
        }
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs_inline(int nrc)
{
    //Current homogeneous polynomial
    for (int i=0; i< FTDA::nmon(nv, nrc); i++)
    {
        (this->term[nrc]->getCA()+i)->tfs_to_ofs_inline();
    }
    return *this;
}


//----------------
// Pows
//----------------
/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_pows(Ofts<T> const& a,  U const& alpha)
{
    cout << "tfts_pows is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_pows(Ofts< Ofsc > const& a,   U const& alpha)
{
    Ofsc temp(corder);

    //Initialization of order zero
    coef0s(*this)->tfs_pows(*coef0s(a), alpha);

    //temp2 = 1/x0
    temp.tfs_pows(*coef0s(a), -1.0+0.0*I);

    //Recurrence scheme
    for(int k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        //Loop on all previously computed homogeneous terms
        for(int j=0; j<= k-1; j++)
        {
            this->term[k]->tftsh_smprod_tu(*a.term[k-j], *this->term[j], temp, (alpha*(k-j)-j)/k);
        }
    }
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients at order n
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_pows(Ofts<T> const& a,  U const& alpha, int const& n)
{
    cout << "tfts_pows is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients at order n. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_pows(Ofts< Ofsc > const& a,   U const& alpha, int const& n)
{
    if(n == 0)
    {
        //Initialization of order zero
        coef0s(*this)->tfs_pows(*coef0s(a), alpha);
    }
    else
    {
        Ofsc temp(corder);
        //temp2 = 1/x0
        temp.tfs_pows(*coef0s(a), -1.0+0.0*I);
        //Sets every coefficients to zero for order n
        this->term[n]->zero();
        //Loop on previously computed coefficient
        for(int j=0; j<= n-1; j++)
        {
            this->term[n]->tftsh_smprod_tu(*a.term[n-j], *this->term[j], temp, (alpha*(n-j)-j)/n);
        }

    }

    return *this;
}

/**
 *   \brief Compute the order zero of the power function: p = a^alpha, with Tfs coefficients at order n
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_spows_zero(Ofts<T> const& a,  U const& alpha, int const& n)
{
    cout << "tfts_pows_zero is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Compute the order zero of the power function: p = a^alpha, with Tfs coefficients at order n. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_spows_zero(Ofts< Ofsc > const& a,   U const& alpha, int const& n)
{
    if(n == 0)
    {
        //Initialization of order zero
        coef0s(*this)->tfs_pows(*coef0s(a), alpha);
    }
    else
    {
        Ofsc temp(corder);
        //temp2 = 1/x0
        temp.tfs_pows(*coef0s(a), -1.0+0.0*I);
        //Only order 0
        int j = 0;
        this->term[n]->tftsh_smprod_tu(*a.term[n-j], *this->term[j], temp, (alpha*(n-j)-j)/n);
    }

    return *this;
}



//----------------
// smult
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_smult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->tftsh_smult_t(*a.term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order n.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::tfts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& n)
{
    if(n > order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->tftsh_smult_u(*a.term[n], c);  //ps[k] += m*a[k]
        return *this;
    }
}

//----------------
// sprod
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.order, b.order)
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i, i0, i1;
    i0 = min(b.order, n);
    i1 = min(a.order, n);
    //Product
    for(i= n-i0; i<=i1; i++) term[n]->tftsh_sprod(*a.term[i], *b.term[n-i]);
    return *this;
}

/**
 *  \brief  An operation. Adds the order zero of the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.order, b.order)
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sprod_zero(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i0 = min(b.order, n);
    int i1 = min(a.order, n);
    //Product
    term[n]->tftsh_sprod(*a.term[n-i0], *b.term[i0]);
    term[n]->tftsh_sprod(*a.term[i1], *b.term[n-i1]);
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //-----------------------------------------------------------------------
    //Other possibility: use a temporary variable in Ofs, and not Ofts
    int i, i0, i1;
    i0 = min(b.order, n);
    i1 = min(a.order, n);
    for(i= n-i0; i<=i1; i++) term[n]->tftsh_smprod_u(*a.term[i], *b.term[n-i], c);

    return *this;
    //-----------------------------------------------------------------------
}


/**
 *  \brief  An operation. Adds the order zero of the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_smprod_u_zero(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //-----------------------------------------------------------------------
    //Other possibility: use a temporary variable in Ofs, and not Ofts
    int i0 = min(b.order, n);
    int i1 = min(a.order, n);
    term[n]->tftsh_smprod_u(*a.term[n-i0], *b.term[i0], c);
    term[n]->tftsh_smprod_u(*a.term[i1], *b.term[n-i1], c);
    return *this;
    //-----------------------------------------------------------------------
}

//----------------
// sfsum
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->tftsh_smult_t(*a.term[n], ma);   //ps[k]  += ma*a[k]
        term[n]->tftsh_smult_t(*b.term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

};

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  c: a reference to an Ofts object
 *  \param  mc: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sfsum_tt(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv|| order != c.order || nv != c.nv)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[n]->tftsh_smult_t(*a.term[n], ma);   //ps[k]  += ma*a[k]
        term[n]->tftsh_smult_t(*b.term[n], mb);   //ps[k]  += mb*b[k]
        term[n]->tftsh_smult_t(*c.term[n], mc);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 *  \param  m: a reference to the order to update
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb, int const& m)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[m]->tftsh_mult_u(*a.term[m], ca);        //ps[k]   = ca*a[k], contains the zeroing
        term[m]->tftsh_smult_u(*b.term[m], cb);       //ps[k]  += mb*b[k]
        return *this;
    }

}

//----------------
// der
//----------------
/**
 *   \brief Partial derivative at order m: works for order = a.order. TFS format.
 *          WARNING: order m is derived and set in order m-1 of this!
 *          If m==0, nothing is done.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfts_der(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in tfts_der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else term[m-1]->tfts_derh(*a.term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *   \brief Same as der but the result is added to the current Ofts instance.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfts_sder(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in tfts_sder @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else term[m-1]->tfts_sderh(*a.term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}


//--------------------------------------------------------------------------------
//Operations with specific conditions on the Fourier-Taylor series
//--------------------------------------------------------------------------------
/**
 *  \fn    template<typename U> Ofts< Ofsd >& Ofts< Ofsd >::ofts_pows(Ofts< Ofsd > const& a,   U const& alpha, int const& n)
 *
 *  \brief Power function: p = a^alpha at order n, for Ofts<Ofsd> objects with Fourier coefficients in Frequency domain.
 *         Moreover, the order zero of the FT series must be unitary.
 **/
template<> template<typename U> Ofts< Ofsd >& Ofts< Ofsd >::ofts_pows(Ofts< Ofsd > const& a,   U const& alpha, int const& n)
{
    if(n==0)
    {
        //Initialization of order zero
        this->acoef0s(*coef0s(a));  //a[0]^alpha = 1.0^alpha = 1;
        return *this;
    }
    else
    {
        //Sets every coefficients to zero for order k
        this->term[n]->zero();
        int i0 = min(a.order, n);
        for(int i= n-i0; i< n; i++) this->term[n]->oftsh_smprod_u(*a.term[n-i], *this->term[i], (double) (alpha*(n-i)-i)/n);
        return *this;
    }
}


/**
 *  \brief Power function: p = a^alpha at order n, for Ofts<Ofsd> objects with Fourier coefficients in Frequency domain.
 *         Moreover, the order zero of the FT series a[0] must satisfy: a[0] >> a[i], for all i > 0
 *         In this routine, the coefficient a0inv = 1/(a[0]) and a0palpha = a[0]^alpha
 **/
template<> template<typename U> Ofts< Ofsd >& Ofts< Ofsd >::pows(Ofts< Ofsd > const& a,  //intitial Ofts
                                                                 Ofsd a0inv,             //inverse of order 0
                                                                 Ofsd a0palpha,          //order 0 ^alpha
                                                                 U const& alpha)         //power coef
{
    //Initialization of order zero
    this->acoef0s(a0palpha);  //a[0]^alpha = a0palpha provided
    //Recurrence scheme
    for(int k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        for(int j=0; j<= k-1; j++)
        {
            this->term[k]->oftsh_smprod_u(*a.term[k-j], *this->term[j], (double) (alpha*(k-j)-j)/k);
            this->term[k]->oftsh_mult_t(*this->term[k], a0inv);
        }

    }
    return *this;
}


//--------------------------------------------------------------------------------
//Stream
//--------------------------------------------------------------------------------
/**
 *   \brief Stream operator << for Ofts objects.
 **/
template<typename T> std::ostream& operator << (std::ostream& stream, Ofts<T> const& ofts)
{
    int i,j, nrc;
    int k[ofts.nv];

    for(nrc=0; nrc<= ofts.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ofts.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ofts.nv, nrc); i++)
        {
            for(j=0; j<ofts.nv; j++) stream <<   setiosflags(ios::right) <<  k[j] << " ";
            stream << endl;
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(16)  <<  ofts.term[nrc]->getCoef(i) << std::noshowpos << endl;

            if(i< FTDA::nmon(ofts.nv, nrc)-1)  FTDA::prxkt(k, ofts.nv);
        }
    }
    return stream;
}

/**
 *   \brief Stream operator for Ofts objects. Print only the order zero of each coefficients.
 **/
template<typename T> void Ofts<T>::fprint_0(ofstream& stream)
{
    int i,j, nrc;
    int k[nv];
    for(nrc=0; nrc<= order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(nv, nrc); i++)
        {
            //for(j=0; j<nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            for(j=0; j<nv; j++) stream <<   setiosflags(ios::right) <<  k[j] << " ";
            stream << endl;
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
            term[nrc]->getCoef(i).fprint_0(stream);
            stream << std::noshowpos << endl;
            if(i< FTDA::nmon(nv, nrc)-1)  FTDA::prxkt(k, nv);
        }
    }
}



//--------------------------------------------------------------------------------
//Evaluate
//--------------------------------------------------------------------------------
/**
 *  \brief Generic routine for the evaluation of an Ofts object at the state X: z = this(X).
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate(U X[], T& z)
{
    z.zero();
    for(int k = order; k >= 0 ; k--)
    {
        term[k]->sevaluate(X, z);
    }
}

/**
 *  \brief Generic routine for the conjugate evaluation of an Ofts object at the state X: z = conj(this(X)).
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate_conjugate(U X[], T& z)
{
    z.zero();
    for(int k = order; k >= 0 ; k--)
    {
        term[k]->sevaluate_conjugate(X, z);
    }
}

/**
 *  \brief Toutine for the evaluation of an Ofts object at the state X, at order m: z = [this(X)]_m.
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate(U X[], T& z, int const& m, int const& ofs_order)
{
    z.zero();
    for(int k = m; k >= 0 ; k--)
    {
        term[k]->sevaluate(X, z, ofs_order);
    }
}

/**
 *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
 **/
template<typename T> template<typename U> cdouble Ofts<T>::fevaluate(U X[], double const& t, int const& m, int const& ofs_order)
{
    //Initialize the state
    cdouble z = 0.0+0.0*I;

    //Initialize the cosinus/sinus arrays
    double cR[ofs_order];
    double sR[ofs_order];

    cR[0] = cos(t);
    sR[0] = sin(t);
    for(int i = 1; i< ofs_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Initialize the exponents array
    int kv[nv];

    //Loop on all desired orders
    for(int k = m; k >= 0 ; k--)
    {
        term[k]->fevaluate(X, z, kv, cR, sR, ofs_order);
    }

    return z;
}

/**
 *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
 *
 *         Contrary to the routine fevaluate(U X[], double const& t, int const& m, int const& ofs_order), the cosinus/sinus arrays are given as inputs:
 *              - cR[] = [cos(t), ..., cos(ofs_order*t)]
 *              - sR[] = [sin(t), ..., sin(ofs_order*t)]
 **/
template<typename T> template<typename U> cdouble Ofts<T>::fevaluate(U X[], double cR[], double sR[], int const& m, int const& ofs_order)
{
    //Initialize the state
    cdouble z = 0.0+0.0*I;

    //Initialize the exponents array
    int kv[nv];

    //Loop on all desired orders
    for(int k = m; k >= 0 ; k--)
    {
        term[k]->fevaluate(X, z, kv, cR, sR, ofs_order);
    }

    return z;
}

/**
 *  \brief Contribution of the order m of this to the evaluation: z = [this(X)]_=m).
 **/
template<typename T> template<typename U> void Ofts<T>::contribution(U X[], T& z, int const& m)
{
    z.zero();
    term[m]->sevaluate(X, z);
}


//---------------------------------------------------------------------------
//Derivation
//---------------------------------------------------------------------------
/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n: \c this \f$ = \frac{\partial a}{\partial z_{n_i}} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::der(Ofts<T> const& a, int ni)
{
    for(int k = 0; k < order ; k++) //Careful here: the sum goes up to order-1!
    {
        term[k]->derh(*a.term[k+1], ni);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n: \c this \f$ += \frac{\partial a}{\partial z_{n_i}} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::sder(Ofts< T > const& a, int ni)
{
    for(int k = 0; k < order ; k++) //Careful here: the sum goes up to order-1!
    {
        term[k]->sderh(*a.term[k+1], ni);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n, at order m of the expansions:
 *          \f$ this_{m-1}  = \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
 *
 *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
 **/
template<typename T> Ofts<T>& Ofts<T>::der(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else term[m-1]->derh(*a.term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n, at order m of the expansions:
 *          \f$ this_{m-1}  += \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
 *
 *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
 **/
template<typename T> Ofts<T>& Ofts<T>::sder(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else term[m-1]->sderh(*a.term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *  \brief Partial derivative wrt to time:
 *          \f$ this  += \left[ \frac{\partial a}{\partial t} \right] \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::dot(Ofts<T> const& a, double const&  n)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
    }
    else
    {
        //Loop
        for(int k = 0; k <= order ; k++) term[k]->dot(*a.term[k], n);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the time, at order k of the expansions:
 *          \f$ this_{k-1}  += \left[ \frac{\partial a}{\partial t} \right]_k \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::dot(Ofts<T> const& a, double const&  n, int const& k)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
    }
    else
    {
        term[k]->dot(*a.term[k], n);
    }
    return *this;
}

//---------------------------------------------------------------------------
//Integral
//---------------------------------------------------------------------------
/**
 *  \brief Primitive wrt to the variable z[ni], with ni = 1,...n: \c this \f$ = \int a dz_{n_i} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::sprim(Ofts< T > const& a, int ni)
{
    for(int k = 0; k < order ; k++) //Careful here: the sum goes up to order-1!
    {
        term[k+1]->sprimh(*a.term[k], ni);
    }
    return *this;
}

//---------------------------------------------------------------------------
//Norms
//---------------------------------------------------------------------------
/**
 *  \brief L1 norm of the term of order m: returns \f$ L_1 \left( [this]_m \right) \f$
 **/
template<typename T> double Ofts<T>::l1norm(int const& m)
{
    return term[m]->l1norm();
}

/**
 *  \brief Infinity norm of the term of order m: returns \f$ L_\infty \left( [this]_m \right) \f$
 **/
template<typename T> double Ofts<T>::linfnorm(int const& m)
{
    return term[m]->linfnorm();
}

/**
 *  \brief  Number of small divisors under a certain value sdmax in the term of order m
 */
template<typename T> int Ofts<T>::nsd(int const& m, int odmax, double sdmax)
{
    return term[m]->nsd(odmax, sdmax);
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Reading & writing
//
//---------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------
// Text format, write
//----------------------------------------------
/**
 * \brief Writes a given object W of type \c Ofts<Ofsc >  in a txt files of the form "filename".
 **/
inline void  writeOFTS_txt(Ofts<Ofsc > &W, string filename)
{
    ofstream myfile;
    myfile.open ((filename).c_str(), ios::out);
    myfile << W << endl;
    myfile.close();
}


/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void  writeVOFTS_txt(vector<Ofts<Ofsc > > &W, string filename)
{
    ofstream myfile;
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        myfile.open ((filename+"["+ss1+"].txt").c_str(), ios::out);
        myfile << W[i] << endl;
        myfile.close();
    }
}

//----------------------------------------------
// Text format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in txt format.
 **/
inline void readOFS_txt(Ofsc &xFFT, ifstream &readStream, int fftN)
{
    //Init
    double ct, cr, ci;
    //Reading
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.setCoef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in txt format.
 **/
inline int readOFTS_txt(Ofts<Ofsc > &x, string filename, int fftN)
{
    //Init
    ifstream readStream;
    string ct;
    //Reading
    readStream.open((filename).c_str(), ios::in);

    //Check that the opening went well
    if (!readStream.is_open())
    {
        cout << "readOFTS_txt. Cannot open file " << filename << endl;
        cout << "Check the text data exist." << endl;
        return GSL_FAILURE;
    }
    else
    {
        for(int k = 0 ; k <= x.getOrder() ; k++)
        {
            for(int p = 0; p < FTDA::nmon(x.getNV(), k); p++)
            {
                //Current kv
                getline(readStream, ct);
                //Reading the coefficient
                readOFS_txt(*x.getCoef(k,p), readStream, fftN);
                getline(readStream, ct);
                getline(readStream, ct);
            }
        }
        readStream.close();
        return GSL_SUCCESS;
    }
}

/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void readVOFTS_txt(vector<Ofts<Ofsc > >  &W, string filename, int fftN)
{
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        readOFTS_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
    }
}


//----------------------------------------------
// Binary format, write
//----------------------------------------------
/**
 * \brief Writes a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void  writeOFS_bin(Ofsc  &xFFT, fstream &myfile)
{
    //Init
    int fftN = xFFT.getOrder();
    double res;

    //Writing
    for(int i = -fftN; i<=fftN; i++)
    {
        //Real part
        res = creal(xFFT.ofs_getCoef(i));
        myfile.write((char*) &res, sizeof(double));

        //Imag part
        res = cimag(xFFT.ofs_getCoef(i));
        myfile.write((char*) &res, sizeof(double));
    }
}

/**
 * \brief Writes a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline void  writeOFTS_bin(Ofts<Ofsc > &W, string filename)
{
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::out);
    //Loop on order
    for(int nrc=0; nrc<= W.getOrder(); nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< FTDA::nmon(W.getNV(), nrc); i++)
        {
            //Write each Ofs coefficient
            writeOFS_bin(*W.getCoef(nrc,i), myfile);
        }
    }
    myfile.close();
}

/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void  writeVOFTS_bin(vector<Ofts<Ofsc > > &W, string filename)
{
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        writeOFTS_bin(W[i], (filename+"["+ss1+"].bin"));
    }
}

//----------------------------------------------
// Binary format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void readOFS_bin(Ofsc  &xFFT, fstream &myfile)
{
    int fftN = xFFT.getOrder();
    double cr, ci;
    //Writing
    for(int i = -fftN; i<=fftN; i++)
    {
        //Real part
        myfile.read((char*)&cr, sizeof(double));
        //Imag part
        myfile.read((char*)&ci, sizeof(double));
        //Put in current position
        xFFT.setCoef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline int readOFTS_bin(Ofts<Ofsc > &W, string filename, int fftN)
{
    //Init
    fstream myfile;

    //Open the stream
    myfile.open((filename).c_str(), ios::binary | ios::in);

    //Check that the opening went well
    if (!myfile.is_open())
    {
        cout << "readOFTS_bin. Cannot open file " << filename << endl;
        cout << "readOFTS_bin. Check the binary data exist." << endl;
        return GSL_FAILURE;
    }
    else
    {
        //Loop on order
        for(int nrc=0; nrc<= W.getOrder(); nrc++)
        {
            //Current homogeneous polynomial
            for (int i=0; i< FTDA::nmon(W.getNV(), nrc); i++)
            {
                //Read each Ofs coefficient
                readOFS_bin(*W.getCoef(nrc,i), myfile);
            }
        }
        myfile.close();
        return GSL_SUCCESS;
    }
}

/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void readVOFTS_bin(vector<Ofts<Ofsc > >  &W, string filename, int fftN)
{
    string ss1;
    int status, global_status;
    global_status = 0;

    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();

        //Read binary format
        status = readOFTS_bin(W[i], (filename+"["+ss1+"].bin"), fftN);

        //Try txt format if failure
        if(status != GSL_SUCCESS)
        {
            cout << "readVOFTS_bin. Last reading went wrong. Trying to find data in txt format..." << endl;
            status = readOFTS_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
            if(status != GSL_SUCCESS)
            {
                cout << "readVOFTS_bin. Txt format also went wrong. Check data manually." << endl;
                global_status--;
            }
            else
            {
                cout << "readVOFTS_bin. Success with txt format." << endl;
                global_status++;
            }
        }
    }

    //If global success, we can save in binary format
    if(global_status == (int) W.size())
    {
        int ch;
        cout << "readVOFTS_bin. Success with txt format for all components." << endl;
        cout << "Please enter 1 if you want to save the vector in binary files:" << endl;
        scanf("%d",&ch);
        if(ch == 1) writeVOFTS_bin(W, filename);
    }

}

//----------------------------------------------
// Text 2 binary
//----------------------------------------------
/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 *        Writes it again in binary form, in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void writeVOFTS_txt2bin(vector<Ofts<Ofsc > > &W, string filename, int fftN)
{
    //Read from txt files
    readVOFTS_txt(W, filename, fftN);
    //Write in binary format
    writeVOFTS_bin(W, filename);
}

//--------------------------------------------------------------------------------
// Ofts 2 Ofs
//--------------------------------------------------------------------------------
/**
 *  \brief Transform a one-variable Ofts< Ots<U> > fts_z into an Ofs<U> object fs object such that fts_z(epsilon) = fs.
 **/
template <typename U>  void fts2fs(Ofs<U> *fs, Ofts< Ots<U> > const& fts_z, double epsilon)
{

    int k, l;
    double tempc;
    int nf = fts_z.getCOrder()/2;
    Ofs<U> tempfs(nf);
    Oftsh< Ots<U> > *tempfth;
    //Cleaning of fs
    fs->zero();
    //Orders >= 0
    for(k = -nf; k <= nf; k++)
    {
        tempc = 0;
        for(l = 0 ; l <= fts_z.getOrder(); l++)
        {
            //Trigonometric coefficient of order l
            tempfth = fts_z.getTerm(l);
            tempfs.ts2fs(tempfth->getCoef(0));
            //tempc += ulj*epsilon^l
            tempc += tempfs.ofs_getCoef(k)*pow(epsilon, (double) l);
        }
        fs->setCoef(tempc, k);
    }
}

/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into an Ofs<U> object fs object such that fts_z(epsilon) = fs.
 **/
template <typename U>  void fts2fs(Ofs<U> *fs, Ofts< Ofs<U> > const& fts_z, double epsilon)
{

    int k, l;
    double tempc;
    int nf = fts_z.getCOrder();
    Ofs<U> tempfs(nf);
    //Cleaning of fs
    fs->zero();
    //Orders >= 0
    for(k = -nf; k <= nf; k++)
    {
        tempc = 0;
        for(l = 0 ; l <= fts_z.getOrder(); l++)
        {
            //Trigonometric coefficient of order l
            tempfs = fts_z.getTerm(l)->getCoef(0);
            //tempc += ulj*epsilon^l
            tempc += tempfs.ofs_getCoef(k)*pow(epsilon, (double) l);
        }
        fs->setCoef(tempc, k);
    }
}

/**
 *  \brief Transform a one-variable Ofts< Ots<U> > fts_z into a complex number fts_z(epsilon, t).
 **/
template <typename U>  cdouble fts2scalar(Ofts< Ots<U> > const& fts_z, double epsilon, double t)
{
    Ofs<U> fs(fts_z.getCOrder()/2);
    fts2fs(&fs, fts_z, epsilon);
    return fs.evaluate(t);
}

/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into a complex number fts_z(epsilon, t).
 **/
template <typename U>  cdouble fts2scalar(Ofts< Ofs<U> > const& fts_z, double epsilon, double t)
{
    Ofs<U> fs(fts_z.getCOrder());
    fts2fs(&fs, fts_z, epsilon);
    return fs.evaluate(t);
}
