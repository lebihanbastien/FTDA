//############################################################################
// Implementation of the Oftsh template class
//############################################################################
/**
 * \file oftsh.tpp
 * \brief Homogeneous Fourier-Taylor series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Oftsh<T>.
 */
template<typename T> Oftsh<T>::Oftsh()
{
    order = 0;
    nv = 0;
    coef = new T(0);
}

/**
 *  \brief Constructor with given order and number of variables.
 *
 * Allocates memory  without any link to a coefficient array (requires the inmediate used of linkCoefs afterwards).
 */
template<typename T> Oftsh<T>::Oftsh(const int newNv, const int newOrder)
{
    order = newOrder;
    nv = newNv;
    coef = 0;

    //Tree
    if(nv==1)  //1-variable polynomial
    {
        term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>)); //calloc is necessary if we are initializing a fourier-taylor serie, with taylor-serie coefficients
        if (term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        //Again, calloc is necessary if we are initializing a fourier-taylor serie, with taylor-serie coefficients
        //At this step, all sons of the current Oftsh<T> have nv = 0, order = 0 and coef points @ one T.
        term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>));
        if (term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=order; i++)
        {
            //term[i].lcopy(Oftsh<T>(newNv-1, newOrder-i));     //WHY is it not working here but in the other allocation routines ??
            term[i].lcopy(*(new Oftsh<T>(newNv-1, newOrder-i)));
            //At this step, all sons of the current Oftsh<T> have nv = nv-1, order = order-i and coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
        }
    }
}

/**
 *  \brief Constructor from a given Oftsh object (without any link).
 */
template<typename T> Oftsh<T>::Oftsh(Oftsh<T> const& b)
{
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Oftsh<T>(b.nv, b.order));
// Goal: avoid the creation of a temporary Oftsh
//-----------------------------------------------------------
    order = b.order;
    nv = b.nv;
    coef = new T(0);//[ Manip::nmon(nv, order)]();
    //Tree
    if(nv==1)  //1-variable polynomial
    {
        term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>));
        if (term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>));  //At this step, all sons of the current Oftsh<T> have nv = 0, order = 0 and coef points @ one T.
        if (term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=order; i++)
        {
            term[i].lcopy(Oftsh<T>(nv-1, order-i));
            //At this step, all sons of the current Oftsh<T> have nv = nv-1, order = order-i and coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
        }
    }
//-----------------------------------------------------------
// End of replacement code
//-----------------------------------------------------------

    //Copy in linking of the coefficients
    T *coef0 = new T[ Manip::nmon(b.nv, b.order)]();
    for(int i = 0 ; i <  Manip::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
    this->linkCoefs(coef0);
}

/**
 *  \brief  An operator. Constructor from a given Oftsh object (only the coefficients).
 */
template<typename T> Oftsh<T>& Oftsh<T>::operator = (Oftsh<T> const& b)
{
    if(this != &b)
    {
        delete term;
        delete coef;
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Oftsh<T>(b.nv, b.order));
// Goal: avoid the creation of a temporary Oftsh
//-----------------------------------------------------------
        order = b.order;
        nv = b.nv;
        coef = new T(0);//[ Manip::nmon(nv, order)]();
        //Tree
        if(nv==1)  //1-variable polynomial
        {
            term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>));
            if (term == NULL)
            {
                puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
                exit(1);
            }
        }
        else
        {
            int i;
            term = (Oftsh<T>*) calloc(order+1, sizeof(Oftsh<T>));  //At this step, all sons of the current Oftsh<T> have nv = 0, order = 0 and coef points @ one T.
            if (term == NULL)
            {
                puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
                exit(1);
            }

            for(i=0; i<=order; i++)
            {
                term[i].lcopy(Oftsh<T>(nv-1, order-i));
                //At this step, all sons of the current Oftsh<T> have nv = nv-1, order = order-i and coef points @ one T.
                //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
            }
        }
//-----------------------------------------------------------
// End of replacement code
//-----------------------------------------------------------
        //Copy in linking of the coefficients
        T *coef0 = new T[ Manip::nmon(b.nv, b.order)]();
        for(int i = 0 ; i <  Manip::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
        this->linkCoefs(coef0);
    }
    return *this; //same object if returned
}


//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Oftsh<T>. WARNING: memory leak here.
 *
 *  TO BE DETERMINED: how properly delete with recursivity? Seems to work fine like this, but memory leak?
 */
template<typename T> Oftsh<T>::~Oftsh<T>()
{
    delete term;
//    if(coef != 0 && coef != NULL)
//    {
//        delete coef;
//        coef = 0;
//    }
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief  Linked copy from a given Oftsh object (exact same object is obtained).
 *  \param  b: a reference to the Oftsh object to copy
 *  \return a reference to the current object
 *
 *  Note: restricted to same order, same number of variables.
 */
template<typename T> Oftsh<T>& Oftsh<T>::lcopy (Oftsh<T> const& b)
{
    order = b.order;
    nv = b.nv;
    term = b.term;
    coef = b.coef;
    return *this;
}

/**
 *  \brief  Copy from a given Oftsh object (only the coefficients).
 *  \param  b: a reference to the Oftsh object to copy
 *  \return a reference to the current object
 */
template<typename T> Oftsh<T>& Oftsh<T>::ccopy (Oftsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int i = 0 ; i< Manip::nmon(nv, order) ; i++) coef[i] = b.coef[i];
        return *this;
    }
}


//---------------------------------------------------------------------------
//Linking
//---------------------------------------------------------------------------
/**
 *  \brief  Performs linking between the Oftsh<T> and a coefficient array.
 */
template<typename T> void Oftsh<T>::linkCoefs(T *coef0)
{
    if(nv==1) //1-variable homogeneous polynomial
    {
        coef = coef0;
    }
    else
    {
        T *coefc;
        int i;
        //allocate the position of the first coefficient
        coef = coef0;
        //perform recursice allocation
        for(i=0, coefc=coef0; i<=order; coefc+= Manip::nmon(nv-1,order-i), i++)
        {
            term[i].linkCoefs(coefc);
        }
    }
}


//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the polynomial.
 */
template<typename T> void Oftsh<T>::setCoef(T const& value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(nv, order)) coef[pos] = value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Adds a coefficient at a given position in the polynomial.
 */
template<typename T> void Oftsh<T>::addCoef(T const& value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(nv, order)) coef[pos] += value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets all subcoefficient of the coefficient \c pos to \c value.
 */
template<typename T> template<typename U> void Oftsh<T>::setSubCoef(U value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(nv, order)) coef[pos].setAllCoefs(value);
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets subcoefficient \c i of the coefficient \c pos to \c value.
 */
template<typename T> template<typename U> void Oftsh<T>::setSubCoef(U value, int pos, int i)
{
    if(pos >=0 && pos <  Manip::nmon(nv, order)) coef[pos].setCoef(value, i);
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets random coefficients to all positions in the polynomial.
 */
template<typename T> void Oftsh<T>::setRandomCoefs()
{
    for(int pos = 0; pos< Manip::nmon(nv, order); pos++)
    {
            coef[pos].setRandomCoefs();
            //coef[pos]/=(double)(order+1);
    }
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the first child in the tree.
 */
template<typename T> Oftsh<T> Oftsh<T>::getTerm() const
{
    return term[0];
}

/**
 *  \brief  Gets the child \c i in the tree.
 *
 *  If position is out of scope, the first term is returned, as in Oftsh<T>::getTerm()
 */
template<typename T> Oftsh<T> Oftsh<T>::getTerm(int i) const
{
    if(i >=0 && i <= order) return term[i];
    else
    {
        cout << "Error in Oftsh<T>::getTerm(int i): position is out of scope." << endl;
        cout << "First term is returned" << endl;
        return term[0];
    }
}

/**
 *  \brief  Gets the address of the first coefficient
 */
template<typename T> T* Oftsh<T>::getCA() const
{
    return coef;
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *  If position is out of scope, the first coefficient is returned.
 */
template<typename T>  T Oftsh<T>::getCoef(int i) const
{
    if(i >=0 && i <  Manip::nmon(nv, order)) return coef[i];
    else
    {
        cout << "Error in Oftsh<T>::getCoef(int i): position is out of scope." << endl;
        cout << "First coefficient is returned" << endl;
        return coef[0];
    }
}

/**
 *  \brief  Gets the order of the polynomial.
 */
template<typename T> int Oftsh<T>::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the number of variables of the polynomial.
 */
template<typename T> int Oftsh<T>::getNV() const
{
    return nv;
}


//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template<typename T> void Oftsh<T>::zero()
{
    for(int i = 0 ; i<  Manip::nmon(nv, order); i++) this->setCoef(0.0, i);
}

/**
 *  \brief  Sets all coefficients to zero. Ofsd case
 */
template<> inline void Oftsh< Ofsd >::zero()
{
    for(int i = 0 ; i<  Manip::nmon(nv, order); i++) this->coef[i].zero();
}

/**
 *  \brief  Sets all coefficients to zero. Ofsc case
 */
template<> inline void Oftsh< Ofsc >::zero()
{
    for(int i = 0 ; i<  Manip::nmon(nv, order); i++) this->coef[i].zero();
}


//--------------------------------------------------------------------------------
//Operations
//--------------------------------------------------------------------------------
//------------------
// Conjugate
//------------------
/**
 *  \brief Conjugates the coefficients the Oftsh object (and only them!). To be used with evaluate_conjugate to have the true conjugate.
 *   The conjugate of the serie \f$ T_n = \sum \limits_{|r| = n} c_r x^r \f$ is
 *  \f$ \bar{T}_n = \sum \sum \limits_{|r| = n} \bar{c}_{r} x^r\f$
 */
template<typename T> Oftsh<T>& Oftsh<T>::conjugate()
{
    for(int pos=0; pos< Manip::nmon(nv, order); pos++) coef[pos].conjugate();
    return *this;
}

//------------------
// Smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a coefficient. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_sprod(a.coef[i], c); //coef[i] += c*a.coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$ with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2, T& temp)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_smprod_t(a.coef[i], c1, c2, temp); //coef[i] += c*a.coef[i]; //
        return *this;
    }
}

//------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_smult(a.coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m * r_a * a \f$ with m a subcoefficient and ra a coefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_smprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}




//------------------
// mult
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a coefficient. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_prod(a.coef[i], c);
        return *this;
    }
}

//------------------

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_mult(a.coef[i], c);
        return *this;
    }
}


//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with a and b Oftsh objects
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_sprod(*aa , *bb);
                    pp->oftsh_smult_t(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->oftsh_smult_t(*bb, *(aa->coef)), bb++ , pp++);
                *(pp->coef)+= *(aa->coef) * *(bb->coef);
            }
            else this->oftsh_smult_t(a, *(b.coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; *pp+= *aa * *bb , bb++ , pp++ ) ;
    }
    else
    {
        *(this->coef)+= *(a.coef ) * *(b.coef); //1-variate homogeneous polynomial
    }
    return *this;
}

//------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with a and b Oftsh objects. Oftsh< Ofsd > double case.
 * Specialization of the routine: template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
 */
template<> inline Oftsh< Ofsd >& Oftsh< Ofsd >::oftsh_sprod(Oftsh< Ofsd > const& a, Oftsh< Ofsd > const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                for(int i = 0; i < a.order ; i++)
                {
                    for(int j = 0; j < b.order; j++) term[i+j].oftsh_sprod(a.term[i], b.term[j]);
                    term[i+b.order].oftsh_smult_t(a.term[i], b.term[b.order].coef[0]);
                }
                for(int j = 0; j < b.order; j++) term[a.order+j].oftsh_smult_t(b.term[j], a.term[a.order].coef[0]);
                term[a.order+b.order].coef[0].ofs_sprod(a.term[a.order].coef[0], b.term[b.order].coef[0]);
            }
            else this->oftsh_smult_t(a, *(b.coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        for(int i = 0; i<= a.order; i++) for(int j = 0; j <= b.order; j++) coef[i+j].ofs_sprod(a.coef[i], b.coef[j]);
    }
    else
    {
        //1-variate homogeneous polynomial
        coef[0].ofs_sprod(a.coef[0], b.coef[0]);
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with a and b Oftsh objects. Oftsh< Ofsc > cdouble case.
 * Specialization of the routine: template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
 */
template<> inline Oftsh< Ofsc >& Oftsh< Ofsc >::oftsh_sprod(Oftsh< Ofsc > const& a, Oftsh< Ofsc > const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                for(int i = 0; i < a.order ; i++)
                {
                    for(int j = 0; j < b.order; j++) term[i+j].oftsh_sprod(a.term[i], b.term[j]);
                    term[i+b.order].oftsh_smult_t(a.term[i], b.term[b.order].coef[0]);
                }
                for(int j = 0; j < b.order; j++) term[a.order+j].oftsh_smult_t(b.term[j], a.term[a.order].coef[0]);
                term[a.order+b.order].coef[0].ofs_sprod(a.term[a.order].coef[0], b.term[b.order].coef[0]);
            }
            else this->oftsh_smult_t(a, *(b.coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        for(int i = 0; i<= a.order; i++) for(int j = 0; j <= b.order; j++) coef[i+j].ofs_sprod(a.coef[i], b.coef[j]);
    }
    else
    {
        //1-variate homogeneous polynomial
        coef[0].ofs_sprod(a.coef[0], b.coef[0]);
    }
    return *this;
}


//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$ with a and b Oftsh objects and c a coefficient.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, T& temp)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_smprod_t(*aa , *bb, c, temp);
                    pp->oftsh_smult_tt(*aa, *(bb->coef), c, temp);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->oftsh_smult_tt(*bb, *(aa->coef), c, temp), bb++ , pp++);
                pp->coef->ofs_smprod_t(*(aa->coef), *(bb->coef), c, temp);
            }
            else this->oftsh_smult_tt(a, *(b.coef), c, temp); //b is scalar
        }
        else this->oftsh_smult_tt(b, *(a.coef), c, temp);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->ofs_smprod_t(*aa, *bb, c, temp);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->ofs_smprod_t(*(a.coef ), *(b.coef), c, temp); //1-variate homogeneous polynomial
    }
    return *this;

//Other implementation less pointers
/*    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {

                for(int i = 0; i < a.order ; i++)
                {
                    for(int j = 0; j < b.order; j++) term[i+j].oftsh_smprod_t(a.term[i], b.term[j], c, temp);
                    term[i+b.order].oftsh_smult_tt(a.term[i], b.term[b.order].coef[0], c, temp);
                }

                for(int j = 0; j < b.order; j++) term[a.order+j].oftsh_smult_tt(b.term[j], a.term[a.order].coef[0], c, temp);
                term[a.order+b.order].coef[0].ofs_smprod_t(a.term[a.order].coef[0], b.term[b.order].coef[0], c, temp);
            }
            else this->oftsh_smult_tt(a, *(b.coef), c, temp); //b is scalar
        }
        else this->oftsh_smult_tt(b, *(a.coef), c, temp);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
          for(int i = 0; i<= a.order; i++) for(int j = 0; j <= b.order; j++) coef[i+j].ofs_smprod_t(a.coef[i], b.coef[j], c, temp);
    }
    else
    {
        coef[0].ofs_sprod(a.coef[0], b.coef[0]);
    }
    return *this;
*/
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += m a \times b \f$ with a and b Oftsh objects. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ofs<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_smprod_u(*aa , *bb, m);
                    pp->oftsh_smult_tu(*aa, *(bb->coef), m);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->oftsh_smult_tu(*bb, *(aa->coef), m), bb++ , pp++);
                pp->coef->ofs_smprod(*(aa->coef), *(bb->coef), m);
            }
            else this->oftsh_smult_tu(a, *(b.coef),m); //b is scalar
        }
        else this->oftsh_smult_tu(b, *(a.coef), m);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        Ofs<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->ofs_smprod(*aa, *bb, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->ofs_smprod(*(a.coef ), *(b.coef), m); //1-variate homogeneous polynomial
    }
    return *this;
}



//---------------------------------------------------------------------------
// TFS operations
//---------------------------------------------------------------------------
//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with a and b Oftsh objects
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_sprod(*aa , *bb);
                    pp->tftsh_smult_t(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->tftsh_smult_t(*bb, *(aa->coef)), bb++ , pp++);
                pp->coef->tfs_sprod(*(aa->coef), *(bb->coef));
            }
            else this->tftsh_smult_t(a, *(b.coef)); //b is scalar
        }
        else this->tftsh_smult_t(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_sprod(*aa, *bb);
    }
    else
    {
        this->coef->tfs_sprod(*(a.coef ), *(b.coef));//1-variate homogeneous polynomial
    }
    return *this;
}

//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$ with a and b Oftsh objects and c a coefficient.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_t(*aa , *bb, c);
                    pp->tftsh_smult_tt(*aa, *(bb->coef), c);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->tftsh_smult_tt(*bb, *(aa->coef), c), bb++ , pp++);
                pp->coef->tfs_smprod_t(*(aa->coef), *(bb->coef), c);
            }
            else this->tftsh_smult_tt(a, *(b.coef), c); //b is scalar
        }
        else this->tftsh_smult_tt(b, *(a.coef), c);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod_t(*aa, *bb, c);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->tfs_smprod_t(*(a.coef ), *(b.coef), c); //1-variate homogeneous polynomial
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$ with a and b Oftsh objects and c a coefficient.
 *
 *  oftsh_smult_tt, ofs_smprod_t
 */
template<typename T> template<typename U> Oftsh<T>& Oftsh<T>::tftsh_smprod_tu(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, U const& m)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_tu(*aa , *bb, c, m);
                    pp->tftsh_smult_ttu(*aa, *(bb->coef), c, m);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->tftsh_smult_ttu(*bb, *(aa->coef), c, m), bb++ , pp++);
                pp->coef->tfs_smprod_tu(*(aa->coef), *(bb->coef), c, m);
            }
            else this->tftsh_smult_ttu(a, *(b.coef), c, m); //b is scalar
        }
        else this->tftsh_smult_ttu(b, *(a.coef), c, m);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod_tu(*aa, *bb, c, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->tfs_smprod_tu(*(a.coef ), *(b.coef), c, m); //1-variate homogeneous polynomial
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += m a \times b \f$ with a and b Oftsh objects. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ofs<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_u(*aa , *bb, m);
                    pp->tftsh_smult_tu(*aa, *(bb->coef), m);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->tftsh_smult_tu(*bb, *(aa->coef), m), bb++ , pp++);
                pp->coef->tfs_smprod(*(aa->coef), *(bb->coef), m);
            }
            else this->tftsh_smult_tu(a, *(b.coef),m); //b is scalar
        }
        else this->tftsh_smult_tu(b, *(a.coef), m);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        Ofs<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod(*aa, *bb, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->tfs_smprod(*(a.coef ), *(b.coef), m); //1-variate homogeneous polynomial
    }
    return *this;
}

//------------------
// smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum: note ofs_smult can be used because it does the same as tfs_smult would do.
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_smult(a.coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum: note ofs_mult can be used because it does the same as tfs_smult would do.
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].ofs_mult(a.coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a coefficient. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smult_t(Oftsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].tfs_sprod(a.coef[i], c); //coef[i] += c*a.coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$ with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].tfs_smprod_t(a.coef[i], c1, c2); //coef[i] += c*a.coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$ with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> template<typename U> Oftsh<T>& Oftsh<T>::tftsh_smult_ttu(Oftsh<T> const& a, T const& c1, T const& c2, U const& m)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].tfs_smprod_tu(a.coef[i], c1, c2, m); //coef[i] += c*a.coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m * r_a * a \f$ with m a subcoefficient and ra a coefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].tfs_smprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}

//------------------
// derh
//------------------
/**
 *  \brief  An operation. Applies the partial derivative with respect to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tfts_derh(Oftsh<T> const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->tfts_derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->setCoef( ((double)a.order+0.0*I)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->tftsh_mult_u(*pp, (double)(pp - a.term)+0.0*I);
            }
        }
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect to the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tfts_sderh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->tfts_sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->addCoef( (a.order+0.0*I)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->tftsh_smult_u(*pp, (pp - a.term)+0.0*I);
            }
        }
    }

    return *this;
}


//---------------------------------------------------------------------------
// Derivation
//---------------------------------------------------------------------------
/**
 *  \brief  An operation. Applies the partial derivative with respect to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::derh(Oftsh<T> const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->setCoef( ((double)a.order+0.0*I)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_mult_u(*pp, (double)(pp - a.term)+0.0*I);
            }
        }
    }
    return *this;
}

/**
 *  \brief  An operation. Applies the partial derivative with respect to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *   Specialization of  Oftsh<T>::derh(Oftsh< T > const& a, int ni) to Ofsd coefficients.
 */
template<> inline Oftsh<Ofsd >& Oftsh<Ofsd >::derh(Oftsh<Ofsd  > const& a, int ni)
{
    Oftsh<Ofsd > *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->setCoef( ((double)a.order)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_mult_u(*pp, (double) (pp - a.term));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect to the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->addCoef( (a.order+0.0*I)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, (pp - a.term)+0.0*I);
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial primitive with respect to the variable \c ni
 *
 *  Notes:
 *  1. If a is of order n, this is of order n+1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::sprimh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {
        for(pp = a.term, dd=this->term; pp <= pf; pp++,dd++)
        {
            dd->sprimh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->addCoef( 1.0/(a.order+1.0+0.0*I)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term+1, pp=a.term; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, 1.0/(pp - a.term+0.0*I+1.0));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect to the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *   Specialization of  Oftsh<T>::sderh(Oftsh< T > const& a, int ni) to Ofsd coefficients.
 */
template<> inline Oftsh<Ofsd >& Oftsh<Ofsd >::sderh(Oftsh< Ofsd > const& a, int ni)
{
    Oftsh<Ofsd > *dd, *pp, *pf;
    pf = a.term+a.order;
    if(this->nv > ni)  //the number of variables is greater than ni
    {

        for(pp = a.term, dd=this->term; pp < pf; pp++,dd++)
        {
            dd->sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->addCoef( ((double)a.order)*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, (double) (pp - a.term));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
 */
template<typename T> Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Derivation of all the coefficients
        for (int i=0 ; i <  Manip::nmon(nv, order); i++) coef[i].dot(a.coef[i], n);
        return *this;
    }
}

//---------------------------------------------------------------------------
// Evaluate
//---------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ = T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::evaluate(U X[], T& z)
{
    //Zeroing z
    z.zero();
    int *kv = (int*) calloc(order, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (U) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
    free(kv);
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ = T_n(X) \f$. double version.
 */
template<typename T>  void Oftsh<T>::evaluate(double X[], T& z)
{
    //Zeroing z
    z.zero();
    int *kv = (int*) calloc(order, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (double) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate(U X[], T& z)
{
    int *kv = (int*) calloc(order, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < nv; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z.ofs_smult(coef[i], (U) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
    free(kv);
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$ at order ofs_order. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate(U X[], T& z, int const& ofs_order)
{
    int *kv = (int*) calloc(nv, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //Evaluate one coefficient
        if(!coef[i].isnull(1))
        {
            //z += X[ii]^kv[ii]*coef(i)
            bux = (U) (1.0+0.0*I);
            for(int ii = 0; ii < nv; ii++)
            {
                aux = (U) (1.0+0.0*I);
                if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
                bux *= aux;
            }
            z.ofs_smult(coef[i], (U) bux, ofs_order);
       }
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv); //update the exponents
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and time t and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluatedc(U X[], U& z, double const& t, int const& ofs_order)
{
    int *kv = (int*) calloc(nv, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < nv; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z += coef[i].evaluate(t, ofs_order)*bux;
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and time t and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::fevaluate(U X[], U& z, int kv[], double cR[], double sR[], int const& ofs_order)
{
    //cout << "oftsh_fevaluate. order = " << order << endl;
    U aux, bux;
    //Initialize kv
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    //Loop on all monomials
    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < nv; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z += coef[i].fevaluate(cR, sR, ofs_order)*bux;
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$. double version
 */
template<typename T> void Oftsh<T>::sevaluate(double X[], T& z)
{
    int *kv = (int*) calloc(order, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    double aux, bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }

        z.ofs_smult(coef[i], bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ = T_n(\bar{X}) \f$.
 */
template<typename T> template<typename U> void Oftsh<T>::evaluate_conjugate(U X[], T& z)
{
    //Zeroing z
    z.zero();
    int kv[nv];
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= conj(X[ii]);
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (U) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
}

/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ = T_n(\bar{X}) \f$.
 */
template<typename T> void Oftsh<T>::evaluate_conjugate(double X[], T& z)
{
    //Zeroing z
    z.zero();
    int kv[nv];
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (double) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
}


/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate_conjugate(U X[], T& z)
{
    int kv[nv];
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= conj(X[ii]);
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (U) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
}

/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
 */
template<typename T> void Oftsh<T>::sevaluate_conjugate(double X[], T& z)
{
    int kv[nv];
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(nv, order); i++)
    {
        //z += X[ii]^kv[ii]*coef(i)
        bux = 1.0;
        for(int ii = 0; ii < nv; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(coef[i], (double) bux);
        if(i< Manip::nmon(nv, order)-1)  Manip::prxkt(kv, nv);
    }
}


/**
 *  \brief  Evaluates the \f$ L_1 \f$ norm of the current Oftsh object.
 */
template<typename T> double Oftsh<T>::l1norm()
{
    double l1n = 0.0;
    if(MODEL_TYPE == Csts::QBCP) for(int k =0; k < Manip::nmon(nv, order); k++) l1n += cabs(this->getCoef(k).l1norm());
    else if(MODEL_TYPE == Csts::CRTBP) for(int k =0; k < Manip::nmon(nv, order); k++) l1n += cabs(this->getCoef(k).ofs_getCoef(0));
    return l1n;
}



/**
 *  \brief  Number of small divisors under a certain value
 */
template<typename T> int Oftsh<T>::nsd(int odmax, double sdmax)
{
    int res = 0;
    for(int k =0; k < Manip::nmon(nv, order); k++) res += this->getCoef(k).nsd(odmax, sdmax);
    return res;
}


/**
 *  \brief  Evaluates the \f$ L_{\infty} \f$ norm of the current Oftsh object.
 */
template<typename T> double Oftsh<T>::linfnorm()
{
    double lin = 0.0;
    if(MODEL_TYPE == Csts::QBCP)
    {
        lin = cabs(this->getCoef(0).l1norm()+0.0*I);
        for(int k =1; k < Manip::nmon(nv, order); k++)
        {
            if(lin < fabs(this->getCoef(k).l1norm())) lin = fabs(this->getCoef(k).l1norm());
        }
    }
    else if(MODEL_TYPE == Csts::CRTBP)
    {
        lin = cabs(this->getCoef(0).ofs_getCoef(0));
        for(int k =1; k < Manip::nmon(nv, order); k++)
        {
           if(lin < cabs(this->getCoef(k).ofs_getCoef(0))) lin = cabs(this->getCoef(k).ofs_getCoef(0));
        }
    }
    return lin;
}

//---------------------------------------------------------------------------
// Functions
//---------------------------------------------------------------------------
/**
 * \fn template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sum a+b
 */
template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop+=b;
    return cop;
}

/**
 * \fn template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sub a-b
 */
template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop-=b;
    return cop;
}

/**
 *  \fn template <typename T>  cdouble smult_expectedError(Oftsh<T> const& a, T const& c, U X[], double const& t)
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t. Works only when a.order = b.order which is the default case.
 */
template <typename U> cdouble smult_expectedError(Oftsh<Ofs<U> > const& a, Ofs<U> const& c, U X[], double const& t)
{
    //Initialisation
    cdouble result = 0.0+0.0*I;
    int *kv = (int*) calloc(a.getOrder(), sizeof(int));
    kv[0] = a.getOrder();
    for(int i=1; i< a.getNV(); i++) kv[i] = 0;
    U aux, bux;

    //Loop on all the monomials in a
    for (int i=0; i< Manip::nmon(a.getNV(), a.getOrder()); i++)
    {
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < a.getNV(); ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }
        //Updating the result
        result += bux*sprod_expectedError(a.getCoef(i), c, t);

        //updating the exponents in kv
        if(i< Manip::nmon(a.getNV(), a.getOrder())-1)  Manip::prxkt(kv, a.getNV());
    }
    free(kv);
    return result;
}

/**
 *  \fn template <typename T>  cdouble smult_expectedError(Oftsh<T> const& a, T const& c, U X[], double const& t)
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t. Works only when a.order = b.order which is the default case. double version
 */
inline cdouble smult_expectedError(Oftsh<Ofsd > const& a, Ofsd const& c, double X[], double const& t)
{
    //Initialisation
    cdouble result = 0.0+0.0*I;
    int *kv = (int*) calloc(a.getOrder(), sizeof(int));
    kv[0] = a.getOrder();
    for(int i=1; i< a.getNV(); i++) kv[i] = 0;
    double aux, bux;

    //Loop on all the monomials in a
    for (int i=0; i< Manip::nmon(a.getNV(), a.getOrder()); i++)
    {
        bux = 1.0;
        for(int ii = 0; ii < a.getNV(); ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }
        //Updating the result
        result += bux*sprod_expectedError(a.getCoef(i), c, t);

        //updating the exponents in kv
        if(i< Manip::nmon(a.getNV(), a.getOrder())-1)  Manip::prxkt(kv, a.getNV());
    }
    free(kv);
    return result;
}



//---------------------------------------------------------------------------
//Stream
//---------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
template<typename T> std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh)
{
    int i,j;
    int k[oftsh.nv];
    k[0] = oftsh.order;
    for(i=1; i<oftsh.nv; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Degree: " << oftsh.order    << endl;
    stream << "#Variables: " << oftsh.nv    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< Manip::nmon(oftsh.nv, oftsh.order); i++)
    {
        for(j=0; j<oftsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  oftsh.coef[i] << std::noshowpos << endl;

        if(i< Manip::nmon(oftsh.nv, oftsh.order)-1)  Manip::prxkt(k, oftsh.nv);
    }
    return stream;
}

/**
 *  \brief  A stream operator. Oftsh<complex double>. Specialization of  template<typename T> std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh)
 */
template<> inline std::ostream& operator << (std::ostream& stream, Oftsh<complex double> const& oftsh)
{
    int i,j;
    int k[oftsh.nv];
    k[0] = oftsh.order;
    for(i=1; i<oftsh.nv; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Degree: " << oftsh.order    << endl;
    stream << "#Variables: " << oftsh.nv    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< Manip::nmon(oftsh.nv, oftsh.order); i++)
    {
        for(j=0; j<oftsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  creal(oftsh.coef[i]) << " " <<  cimag(oftsh.coef[i]) << std::noshowpos << endl;

        if(i< Manip::nmon(oftsh.nv, oftsh.order)-1)  Manip::prxkt(k, oftsh.nv);
    }
    return stream;
}
