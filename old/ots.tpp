//############################################################################
// Implementation of the Ots template class
//############################################################################

/**
 * \file ots.tpp
 * \brief Taylor series template class
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */



//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ots<T>.
 */
template<typename T> Ots<T>::Ots()
{
    int i, index;

    nv     = REDUCED_NV;
    order  = OTS_ORDER;
    T *coefs = new T[binomial(nv+order, nv)]();


    //Allocation of the homogeneous polynomials
    term = new Otsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Otsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setAllCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor with given order of for the Taylor series. The number of variables is set to REDUCED_NV by default.
 */
template<typename T> Ots<T>::Ots(int newOrder)
{
    int i, index;

    nv       = REDUCED_NV;
    order    = newOrder;
    T *coefs = new T[binomial(nv+order, nv)]();


    //Allocation of the homogeneous polynomials
    term = new Otsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Otsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setAllCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables for the Taylor series.
 */
template<typename T> Ots<T>::Ots(int newNv, int newOrder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    T *coefs = new T[binomial(nv+order, nv)]();


    //Allocation of the homogeneous polynomials
    term = new Otsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Otsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setAllCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}


//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/*
    Construction without any link (entirely new object!)
*/
template<typename T> Ots<T>::Ots(Ots<T> const& b)
{

    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;

    //Copy of all the coefficients at every order in new array
    T *coef0 = new T[binomial(b.nv+b.order, b.nv)]();
    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i] = b.term[nrc]->getCoef(i);
        index+=  FTDA::nmon(nv,nrc);
    }

    //Allocation of the homogeneous polynomials
    term = new Otsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=order; nrc++)
    {
        //Allocation of each hp
        term[nrc] = new Otsh<T>(nv, nrc);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[nrc]->setAllCoefs(coef0+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

/*
    Copy without any link (entirely new object!)
*/
template<typename T> Ots<T>& Ots<T>::operator = (Ots<T> const& b)
{
    if(this != &b)
    {
        int nrc, index;

        //Same nv/order
        nv = b.nv;
        order = b.order;


        //Copy of all the coefficients at every order in new array
        T *coef0 = new T[binomial(b.nv+b.order, b.nv)]();
        index = 0;
        for(int nrc=0; nrc<= order; nrc++)
        {
            for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i] = b.term[nrc]->getCoef(i);
            index+=  FTDA::nmon(nv,nrc);
        }

        delete term;

        //Allocation of the homogeneous polynomials
        term = new Otsh<T>*[order+1];

        //Allocation of the coefficients
        index = 0;
        for(nrc=0; nrc<=order; nrc++)
        {
            //Allocation of each hp
            term[nrc] = new Otsh<T>(nv, nrc);//allocate_homog(nv, i);
            //Link h to coefs at each level of the tree
            term[nrc]->setAllCoefs(coef0+index);
            index+=  FTDA::nmon(nv,nrc);
        }

    }

    return *this;
}

/*
    Linked copy
*/
template<typename T> Ots<T>& Ots<T>::lcopy (Ots<T> const& b)
{
    order = b.order;
    nv = b.nv;
    term = b.term;
    return *this;
}

/*
    Coefficient copy (restricted to same order, same number of variables)
*/
template<typename T> Ots<T>& Ots<T>::ccopy (Ots<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ots<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(b.term[nrc]->getCoef(i), i);
        return *this;
    }
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
//Order of the expansion
template<typename T> int Ots<T>::getOrder() const
{
    return order;
}

//Coefficient at position pos
template<typename T> T Ots<T>::getCoef(int const& mOrder, int const& pos) const
{
    if(mOrder > order || pos >= FTDA::nmon(nv, mOrder))
    {
        cout << "Error in getCoef: out of range. First term is returned" << endl;
        cout << "Requested order: " << mOrder << "Maximum allowed: " <<  order << endl;
        cout << "Requested pos: " << pos << "Maximum allowed: " <<  FTDA::nmon(nv, mOrder) << endl;
        return this->term[0]->getCoef(0);
    }
    else return this->term[mOrder]->getCoef(pos);
}

//Number of variables
template<typename T> int Ots<T>::getNV() const
{
    return nv;
}

//Get Term of order mOrder (OTSH object)
template<typename T> Otsh<T>* Ots<T>::getTerm(int const& mOrder) const
{
    if(mOrder > order)
    {
        cout << "Error in getTerm: out of range. First term is returned" << endl;
        cout << "Requested order: " << mOrder << ", Maximum allowed: " <<  order << endl;
        return this->term[0];
    }
    else return this->term[mOrder];
}

//---------------------------------------------------------------------------
//Delete
//  TO BE DETERMINED: how properly delete with recursivity?
// Seems to work fine like this, but memory leak?
//---------------------------------------------------------------------------
template<typename T> Ots<T>::~Ots<T>()
{
    if(term != NULL)
    {
        //Certainly a problem at this point: since the delete routine of Otsh is empty, only the first leaf of the Otsh tree is deleted...
        for(int i =0; i<= order ; i++) delete term[i];
        term = 0;
    }
}

//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
template<typename T> void Ots<T>::zero()
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(0.0, i);
}

template<> inline void Ots<cdouble>::zero()
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(0.0+0.0*I, i);
}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
template<typename T> void Ots<T>::setAllCoefs(T const& m)
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(m, i);
}

template<typename T> void Ots<T>::setCoef(T const& m, int const& pos)
{
    if(pos >= binomial(nv+order, nv))
    {
        cout << "Error in setCoef: out of range. Nothing is done" << endl;
        cout << "Requested pos: " << pos << "Maximum allowed: " <<  FTDA::nmon(nv, order) << endl;
    }
    else this->term[0]->setCoefTot(m, pos);
}

template<typename T> void Ots<T>::setCoef(T const& m, int pos, int i)
{
    term[pos]->setCoef(m, i);
}

template <typename T> void Ots<T>::conjugate()
{
    for(int i=0; i<= order; i++) term[i]->conjugate();
}

template <typename T> void Ots<T>::conjugateforOFS()
{
    for(int i=0; i<= order; i++) term[i]->conjugateforOFS();
}

template<typename T> void Ots<T>::setRandomCoefs()
{
    for(int nrc=0; nrc<= order/2; nrc++) term[nrc]->setRandomCoefs();
}

//---------------------------------------------------------------------------
// Functions
//---------------------------------------------------------------------------
/*
    Product: p = p + a*b at order p->order
*/
template<typename T> Ots<T>& Ots<T>::sprod(Ots<T> const& a, Ots<T> const& b)
{
    int k, i, i0, i1;
    //Product
    for(k=0; k<=order; k++)
    {
        i0 = min(b.order, k);
        i1 = min(a.order, k);
        for(i= k-i0; i<=i1; i++) term[k]->sprod(*a.term[i], *b.term[k-i]);
    }
    return *this;
}

/*
    Product: p = p + a*b only at order n
    Handle the case for which n >= max(a.order, b.order)
*/
template<typename T> Ots<T>& Ots<T>::sprod(Ots<T> const& a, Ots<T> const& b, int const& n)
{
    int i, i0, i1;

    i0 = min(b.order, n);
    i1 = min(a.order, n);

    //Product
    for(i= n-i0; i<=i1; i++) term[n]->sprod(*a.term[i], *b.term[n-i]);

    return *this;
}

/*
    Product: p = p + a*b at order p->order
*/
template<typename T> Ots<T>& Ots<T>::prod(Ots<T> const& a, Ots<T> const& b)
{
    this->zero();
    this->sprod(a,b);

    return *this;
}

/*
    Product: p = p + m*a*b at order p->order
*/
template<typename T> Ots<T>& Ots<T>::smprod(Ots<T> const& a, Ots<T> const& b, T const& m)
{
    int k, i;
    //Product
    for(k=0; k<=order; k++)
    {
        for(i=0; i<=k; i++)
        {
            term[k]->smprod(*a.term[i], *b.term[k-i], m);
        }
    }

    return *this;
}

/*
    Product: p = p + m*a*b at order n
*/
template<typename T> Ots<T>& Ots<T>::smprod(Ots<T> const& a, Ots<T> const& b, T const& m, int const& n)
{
    int i, i0, i1;

    i0 = min(b.order, n);
    i1 = min(a.order, n);
    //Product
    for(i= n-i0; i<=i1; i++) term[n]->smprod(*a.term[i], *b.term[n-i], m);
    return *this;
}

/*
    Product: p = m*a*b at order p->order
*/
template<typename T> Ots<T>& Ots<T>::mprod(Ots<T> const& a, Ots<T> const& b, T const& m)
{
    this->zero();
    this->smprod(a,b, m);
    return *this;
}

/*
    S-Mult: ps = ps + r*a at order ps->order
*/
template<typename T> Ots<T>& Ots<T>::smult(Ots<T> const& a, T const& m)
{

    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ots<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->smult(*a.term[k], m);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

/*
    Mult: ps = r*a at order ps->order
*/
template<typename T> Ots<T>& Ots<T>::mult(Ots<T> const& a, T const& m)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Ots<T> is returned" << endl;
        return *this;
    }
    else
    {
        int k;
        for(k=0; k <= order; k++)
        {
            term[k]->mult(*a.term[k], m);  //ps[k] = m*a[k], contains the zeroing
        }

        return *this;
    }
}

/*
    Mult: ps = r*ps at order ps->order
*/
template<typename T> Ots<T>& Ots<T>::mult(T const& m)
{
        for(int k=0; k <= order; k++) term[k]->mult(m);  //ps[k] *= m
        return *this;
}

/*
    Factored Sum: p = ma*a+mb*b at order p->order
*/
template<typename T> Ots<T>& Ots<T>::fsum(Ots<T> const& a, T const& ma, Ots<T> const& b, T const& mb)
{
    if(order != a.order || nv != a.nv || order != b.order || nv != b.nv)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ots<T> is returned" << endl;
        return *this;
    }
    else
    {
        int k;
        for(k=0; k <= order; k++)
        {
            term[k]-> mult(*a->term[k], ma);   //ps[k]   = ma*a[k], contains the zeroing
            term[k]->smult(*b->term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/*
    Division of taylor functions: p = a/b
*/
template<typename T> Ots<T>& Ots<T>::divs(Ots<T> const& a, Ots<T> const& b)
{
    int k, j;

    //Order 0
    this->acoef0s(coef0s(a) / coef0s(b));

    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //sets every coefficients to zero for order k
        this->term[k]->zero();
        //Copy of a[k] in p[k]
        this->term[k]->smult(*a.term[k], 1.0);
        for(j=0; j<= k-1; j++)  this->term[k]->smprod(*b.term[k-j], *this->term[j], -1.0);
        this->term[k]->mult(1.0/coef0s(b));
    }

    return *this;
}

/*
    cdouble spec
*/
template<> inline Ots<cdouble>& Ots<cdouble>::divs(Ots<cdouble> const& a, Ots<cdouble> const& b)
{
    int k, j;

    //Order 0
    this->acoef0s(coef0s(a) / coef0s(b));

    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //sets every coefficients to zero for order k
        this->term[k]->zero();
        //Copy of a[k] in p[k]
        this->term[k]->smult(*a.term[k], 1.0+0.0*I);
        for(j=0; j<= k-1; j++)  this->term[k]->smprod(*b.term[k-j], *this->term[j], -1.0+0.0*I);
        this->term[k]->mult(1.0/coef0s(b));
    }

    return *this;
}

/*
    Division of taylor functions: p = a/(m*b)
*/
template<typename T> Ots<T>& Ots<T>::mdivs(Ots<T> const& a, Ots<T> const& b, T const& m)
{
    int k, j;

    //Order 0
    this->acoef0s(coef0s(a) / (m*coef0s(b)));

    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        //Copy of a[k] in p[k]
        this->term[k]->smult(*a.term[k], 1.0);
        for(j=0; j<= k-1; j++)  this->term[k]->smprod(*b.term[k-j], *this->term[j], -m);
        this->term[k]->mult(1.0/(m*coef0s(b)));
    }

    return *this;
}

/**
 *   \brief Power function: p = a^alpha
 **/
template<typename T> Ots<T>& Ots<T>::pows(Ots<T> const& a,  T const& alpha)
{
    int k, j;
    T x0;

    //Initialization of order zero
    x0 = coef0s(a);

    this->acoef0s(cpow(x0, alpha));
    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        for(j=0; j<= k-1; j++) this->term[k]->smprod(*a.term[k-j], *this->term[j], alpha*(k-j)-j);// smprodh(ps->term[k], s->term[k-j], ps->term[j], alpha*(k-j)-j);
        this->term[k]->mult(1.0/(k*x0)); //multh(ps->term[k], 1.0/(x0*k));
    }
    return *this;
}

/*
    Returns the address of the first coefficient of order 0 of the taylor serie s
*/
template<typename T> T Ots<T>::coef0s(Ots<T> const& a)
{
    return a.term[0]->getCoef(0);
}

/*
   Sets the coefficient of order 0 of the taylor serie s equal to x0
*/
template<typename T> void Ots<T>::acoef0s(T const& x0)
{
    this->term[0]->setCoef(x0, 0);
}

/*
    Complex exponentiation
*/
template<typename T> T cpow(T const& x, T const& alpha)
{
    //Generic version does nothing!
    return cpow(x, alpha);
}

template<> inline complex double cpow< complex double>(complex double const& x, complex double const& alpha)
{
    //Complex version
    return cpow(x, alpha);
}

template<> inline double cpow<double>(double const& x, double const& alpha)
{
    //Check that we dont have simultaneously: x < 0 and alpha not integer:
    if( x <=0 && floor(alpha) != alpha)
    {
        cout << "error in  cpow<double>: result is complex. 1.0 is returned" << endl;
        return 1.0;
    }

    //Double version is simply power
    if(alpha > 0) return pow(x, alpha);
    else return 1.0/pow(x, -alpha);
}


/*
    Returns a/c at order a.order
*/
template <typename T>  Ots<T> operator / (Ots<T> const& a, Ots<T> const& b)
{
    Ots<T> cop(a);
    cop.divs(a,b);
    return cop;
}


/*
    Returns c*a at order a.order
*/
template <typename T>  Ots<T> operator * (T const& c, Ots<T> const& a)
{
    Ots<T> cop(a);
    cop.mult(c);
    return cop;
}


//Derivation
//---------------------------------------------------------------------------
//Partial derivation; Works for order = a.order and order = a.order-1
template<typename T> Ots<T>& Ots<T>::der(Ots< T > const& a, int ni)
{
    int kmax;
    if(order == a.order) kmax = order-1;
    else if(order == a.order-1) kmax = order;
    else
    {
        cout << "Error in der routines (ots object). Orders do not match." << endl;
        return *this;
    }


    for(int k = 0; k <= kmax ; k++)
    {
        term[k]->derh(*a.term[k+1], ni);
    }

    return *this;
}

//Partial derivation; Works for order = a.order and order = a.order-1
template<typename T> Ots<T>& Ots<T>::sder(Ots< T > const& a, int ni)
{
    int kmax;
    if(order == a.order) kmax = order-1;
    else if(order == a.order-1) kmax = order;
    else
    {
        cout << "Error in der routines (ots object). Orders do not match." << endl;
        return *this;
    }


    for(int k = 0; k <= kmax ; k++)
    {
        term[k]->sderh(*a.term[k+1], ni);
    }

    return *this;
}

//Partial derivation @order m;
//WARNING: order m of a is derived and put order m-1 in this
template<typename T> Ots<T>& Ots<T>::der(Ots< T > const& a, int ni, int m)
{
    if(m > 0) term[m-1]->derh(*a.term[m], ni);
    return *this;
}


//---------------------------------------------------------------------------
//Stream
//---------------------------------------------------------------------------
//Double version
template<typename T> std::ostream& operator << (std::ostream& stream, Ots<T> const& ots)
{
    int i,j, nrc;
    int k[ots.nv];

    stream << "#Taylor serie"             << endl;
    stream << "#Degree: " << ots.order    << endl;
    stream << "#Variables: " << ots.nv    << endl;
    stream << "--------------------------"<< endl;

    for(nrc=0; nrc<= ots.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ots.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ots.nv, nrc); i++)
        {
            for(j=0; j<ots.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ots.term[nrc]->getCoef(i) << std::noshowpos << endl;

            if(i< FTDA::nmon(ots.nv, nrc)-1)  FTDA::prxkt(k, ots.nv);
        }
    }
    return stream;
}

//Complex version
template<> inline std::ostream& operator << (std::ostream& stream, Ots< complex double > const& ots)
{
    int i,j, nrc;
    int k[ots.nv];

    stream << "#Taylor serie"             << endl;
    stream << "#Degree: " << ots.order    << endl;
    stream << "#Variables: " << ots.nv    << endl;
    stream << "--------------------------"<< endl;

    for(nrc=0; nrc<= ots.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ots.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ots.nv, nrc); i++)
        {
            for(j=0; j<ots.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  creal(ots.term[nrc]->getCoef(i)) << " " << cimag(ots.term[nrc]->getCoef(i)) << std::noshowpos << endl;

            if(i< FTDA::nmon(ots.nv, nrc)-1)  FTDA::prxkt(k, ots.nv);
        }
    }
    return stream;
}


//---------------------------------------------------------------------------
//Evaluate
//---------------------------------------------------------------------------
//Evaluate
template<typename T> template<typename U> void Ots<T>::evaluate(U X[], T& z)
{
    z = 0;
    for(int k = order; k >= 0 ; k--)
    {
        term[k]->sevaluate(X, z);
    }
}

//Evaluate
template<typename T> template<typename U> void Ots<T>::evaluate_conjugate(U X[], T& z)
{
    z = 0;
    for(int k = order; k >= 0 ; k--)
    {
        term[k]->sevaluate_conjugate(X, z);
    }
}
