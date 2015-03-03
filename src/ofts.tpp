//############################################################################
// Implementation of the Ofts template class
//############################################################################

//Create
//____________________________________________________________________________________________________________________

template<typename T> Ofts<T>::Ofts(int newNv, int nr, int newCnv, int newCorder, T *coefs)
{
    int i, index;

    nv = newNv;
    order = nr;
    cnv = newCnv;
    corder = newCorder;

    //Allocation of the homogeneous polynomials
    term = new Oftsh<T>*[nr+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=nr; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

//Generic (works for double and double complex coefficients)
template<typename T> Ofts<T>::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    T *coefs = new T[binomial(nv+order, nv)]();


    //Allocation of the homogeneous polynomials
    term = new Oftsh<T>*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh<T>(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}

// Ofs<double> case
template<> inline Ofts< Ofs<double> >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;
    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    Ofs<double> *coefs = (Ofs<double>*) calloc(binomial(nv+order,nv), sizeof(Ofs<double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ofs<double>(corder);
    //-------------------------------------------------------


    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ofs<double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh< Ofs<double> >(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}


//Ots<double> case
template<> inline Ofts< Ots<double> >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;

    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    Ots<double> *coefs = (Ots<double>*) calloc(binomial(nv+order,nv), sizeof(Ots<double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ots<double>(cnv, corder);
    //-------------------------------------------------------


    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ots<double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh< Ots<double> >(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}


// Ots<complex double> case
template<> inline Ofts< Ots<complex double> >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    nv = newNv;
    order = newOrder;

    cnv = newCnv;
    corder = newCorder;

    //New array of coefficients
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();

    Ots<complex double> *coefs = (Ots<complex double>*) calloc(binomial(nv+order,nv), sizeof(Ots<complex double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coefs[k]  = Ots<complex double>(cnv, corder);

    //-------------------------------------------------------
    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ots<complex double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=order; i++)
    {
        //Allocation of each hp
        term[i] = new Oftsh< Ots<complex double> >(nv, i);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[i]->setCoefs(coefs+index);
        index+=  FTDA::nmon(nv,i);
    }
}


//Copy
//____________________________________________________________________________________________________________________

// Construction without any link (entirely new object!)
template<typename T> Ofts<T>::Ofts(Ofts<T> const& b)
{

    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    T *coef0 = new T[binomial(nv+b.order, b.nv)]();
    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i] = b.term[nrc]->getCoef(i);
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
        term[nrc]->setCoefs(coef0+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

// Ofs<double> case
template<> inline Ofts<Ofs< double> >::Ofts(Ofts< Ofs<double> > const& b)
{

    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    Ofs<double> *coef0 = (Ofs<double>*) calloc(binomial(nv+order,nv), sizeof(Ofs<double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coef0[k]  = Ofs<double>(corder);
    //-------------------------------------------------------

    //Copy of the coefficients
    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i] = b.term[nrc]->getCoef(i);
        index+=  FTDA::nmon(nv,nrc);
    }

    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ofs<double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=order; nrc++)
    {
        //Allocation of each hp
        term[nrc] = new Oftsh< Ofs<double> >(nv, nrc);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[nrc]->setCoefs(coef0+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

// Ots<double> case
template<> inline Ofts<Ots< double> >::Ofts(Ofts< Ots<double> > const& b)
{

    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    Ots<double> *coef0 = (Ots<double>*) calloc(binomial(nv+order,nv), sizeof(Ots<double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coef0[k]  = Ots<double>(cnv, corder);
    //-------------------------------------------------------

    //Copy of the coefficients
    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i].ccopy(b.term[nrc]->getCoef(i));
        index+=  FTDA::nmon(nv,nrc);
    }

    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ots<double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=order; nrc++)
    {
        //Allocation of each hp
        term[nrc] = new Oftsh< Ots<double> >(nv, nrc);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[nrc]->setCoefs(coef0+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

// Ots<complex double> case
template<> inline Ofts<Ots< complex double> >::Ofts(Ofts< Ots<complex double> > const& b)
{

    int nrc, index;

    //Same nv/order
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    //-------------------------------------------------------
    //Replacing T *coefs = new T[binomial(nv+order, nv)]();
    Ots<complex double> *coef0 = (Ots<complex double>*) calloc(binomial(nv+order,nv), sizeof(Ots<complex double>));
    for(unsigned k= 0; k< binomial(nv + order, nv); k++) coef0[k]  = Ots<complex double>(cnv, corder);
    //-------------------------------------------------------

    //Copy of the coefficients
    index = 0;
    for(int nrc=0; nrc<= order; nrc++)
    {
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i].ccopy(b.term[nrc]->getCoef(i));
        index+=  FTDA::nmon(nv,nrc);
    }

    //Allocation of the homogeneous polynomials
    term = new Oftsh< Ots<complex double> >*[order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=order; nrc++)
    {
        //Allocation of each hp
        term[nrc] = new Oftsh< Ots<complex double> >(nv, nrc);//allocate_homog(nv, i);
        //Link h to coefs at each level of the tree
        term[nrc]->setCoefs(coef0+index);
        index+=  FTDA::nmon(nv,nrc);
    }
}

//---------------------------------------


// Copy without any link (entirely new object!)
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
            term[nrc]->setCoefs(coef0+index);
            index+=  FTDA::nmon(nv,nrc);
        }

    }

    return *this;
}

//Ots<double> case
template<> inline Ofts< Ots< double> >& Ofts< Ots< double> >::operator = (Ofts< Ots< double> > const& b)
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
        //-------------------------------------------------------
        //Replacing T *coefs = new T[binomial(nv+order, nv)]();
        Ots<double> *coef0 = (Ots<double>*) calloc(binomial(nv+order,nv), sizeof(Ots<double>));
        for(unsigned k= 0; k< binomial(nv + order, nv); k++) coef0[k]  = Ots<double>(cnv, corder);
        //-------------------------------------------------------

        //Copy of the coefficients
        index = 0;
        for(int nrc=0; nrc<= order; nrc++)
        {
            for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i].ccopy(b.term[nrc]->getCoef(i));
            index+=  FTDA::nmon(nv,nrc);
        }


        delete term;
        //Allocation of the homogeneous polynomials
        term = new Oftsh< Ots<double> >*[order+1];



        //Allocation of the coefficients
        index = 0;
        for(nrc=0; nrc<=order; nrc++)
        {
            //Allocation of each hp
            term[nrc] = new Oftsh< Ots<double> >(nv, nrc);//allocate_homog(nv, i);
            //Link h to coefs at each level of the tree
            term[nrc]->setCoefs(coef0+index);
            index+=  FTDA::nmon(nv,nrc);
        }

    }

    return *this;
}

// Ots<complex double> case
template<> inline Ofts< Ots< double complex> >& Ofts< Ots< double complex> >::operator = (Ofts< Ots< double complex> > const& b)
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
        //-------------------------------------------------------
        //Replacing T *coefs = new T[binomial(nv+order, nv)]();
        Ots<double complex> *coef0 = (Ots<double complex>*) calloc(binomial(nv+order,nv), sizeof(Ots<double complex>));
        for(unsigned k= 0; k< binomial(nv + order, nv); k++) coef0[k]  = Ots<double complex>(cnv, corder);
        //-------------------------------------------------------

        //Copy of the coefficients
        index = 0;
        for(int nrc=0; nrc<= order; nrc++)
        {


            for (int i=0; i< FTDA::nmon(nv, nrc); i++)  coef0[index+i].ccopy(b.term[nrc]->getCoef(i));
            index+=  FTDA::nmon(nv,nrc);
        }

        delete term;
        //Allocation of the homogeneous polynomials
        term = new Oftsh< Ots<double complex> >*[order+1];
        //Allocation of the coefficients
        index = 0;
        for(nrc=0; nrc<=order; nrc++)
        {
            //Allocation of each hp
            term[nrc] = new Oftsh< Ots<double complex> >(nv, nrc);//allocate_homog(nv, i);
            //Link h to coefs at each level of the tree
            term[nrc]->setCoefs(coef0+index);
            index+=  FTDA::nmon(nv,nrc);
        }

    }

    return *this;
}

//---------------------------------------

// Linked copy
template<typename T> Ofts<T>& Ofts<T>::lcopy (Ofts<T> const& b)
{
    order = b.order;
    nv = b.nv;
    term = b.term;
    cnv = b.cnv;
    corder = b.corder;
    return *this;
}

//---------------------------------------

// Coefficient copy (restricted to same order, same number of variables)
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

// Ots<double> case
template<> inline Ofts< Ots<double> >& Ofts< Ots<double> >::ccopy (Ofts< Ots<double> > const& b)
{
    if(order != b.order || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  (term[nrc]->getCoefAddress()+i)->ccopy(b.term[nrc]->getCoef(i));
        return *this;
    }
}

// Ots<complex double> case
template<> inline Ofts< Ots<complex double> >& Ofts< Ots<complex double> >::ccopy (Ofts< Ots<complex double> > const& b)
{
    if(order != b.order || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  (term[nrc]->getCoefAddress()+i)->ccopy(b.term[nrc]->getCoef(i));
        return *this;
    }
}


// Coefficient copy (restricted to same order, same number of variables)
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

// Ots<double> case
template<> inline Ofts< Ots<double> >& Ofts< Ots<double> >::ccopy (Ofts< Ots<double> > const& b, int const& nrc)
{
    if(nrc > min(order, b.order) || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  (term[nrc]->getCoefAddress()+i)->ccopy(b.term[nrc]->getCoef(i));
        return *this;
    }
}

// Ots<complex double> case
template<> inline Ofts< Ots<complex double> >& Ofts< Ots<complex double> >::ccopy (Ofts< Ots<complex double> > const& b, int const& nrc)
{
    if(nrc > min(order, b.order) || nv != b.nv || corder != b.corder || cnv != b.cnv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for (int i=0; i< FTDA::nmon(nv, nrc); i++)  (term[nrc]->getCoefAddress()+i)->ccopy(b.term[nrc]->getCoef(i));
        return *this;
    }
}



//Getters
//____________________________________________________________________________________________________________________

template<typename T> int Ofts<T>::getOrder() const
{
    return order;
}

template<typename T> int Ofts<T>::getCOrder() const
{
    return corder;
}

template<typename T> int Ofts<T>::getVariables() const
{
    return nv;
}

template<typename T> int Ofts<T>::getCVariables() const
{
    return cnv;
}

template<typename T> T* Ofts<T>::getCoef(int const& mOrder, int const& pos) const
{
    if(mOrder > order || pos >= FTDA::nmon(nv, mOrder))
    {
        cout << "Error in getCoef: out of range. First term is returned" << endl;
        cout << "Requested order: " << mOrder << ", Maximum allowed: " <<  order << endl;
        cout << "Requested pos: " << pos << ", Maximum allowed: " <<  FTDA::nmon(nv, mOrder) << endl;
        return this->term[0]->getCoefAddress();
    }
    else return this->term[mOrder]->getCoefAddress()+pos;
}

template<typename T> Oftsh<T>* Ofts<T>::getTerm(int const& mOrder) const
{
    if(mOrder > order)
    {
        cout << "Error in getTerm: out of range. First term is returned" << endl;
        cout << "Requested order: " << mOrder << ", Maximum allowed: " <<  order << endl;
        return this->term[0];
    }
    else return this->term[mOrder];
}


//Delete
//____________________________________________________________________________________________________________________

template<typename T> Ofts<T>::~Ofts<T>()
{
    if(term != NULL)
    {
        //Certainly a problem at this point: since the delete routine of Oftsh is empty, only the first leaf of the Oftsh tree is deleted...
        for(int i =0; i<= order ; i++) delete term[i];
        term = 0;
    }
}

//Zeroing
//____________________________________________________________________________________________________________________

template<typename T> void Ofts<T>::zero()
{
    for(int nrc=0; nrc<= order; nrc++) term[nrc]->zero();
}


//Setters
//____________________________________________________________________________________________________________________

// Set of given coefficient, whatever his type is
template<typename T> void Ofts<T>::setCoefs(T const& m)
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setCoef(m, i);
}
//Set of a given double/complex (typename U) coefficient everywhere
template<typename T> template < typename U > void Ofts<T>::setCoefs(U const& m)
{
    for(int nrc=0; nrc<= order; nrc++) for (int i=0; i< FTDA::nmon(nv, nrc); i++)  term[nrc]->setTCoef(m, i);
}

template<typename T> void Ofts<T>::setRandomCoefs()
{
    //srand(time(NULL)); //random seed
    for(int nrc=0; nrc<= order; nrc++) term[nrc]->setRandomCoefs();
}

template<typename T> template < typename U > void Ofts<T>::setCoef(U const& m, int const& n)
{
    term[n]->setT0Coef(m, 0);
}


//Conjugate
template<typename T> Ofts<T>& Ofts<T>::conjugate()
{
     for(int nrc=0; nrc<= order; nrc++) term[nrc]->conjugate();
     return *this;
}

template<typename T> Ofts<T>& Ofts<T>::conjugate(int const& nrc)
{
     term[nrc]->conjugate();
     return *this;
}


// Functions
//____________________________________________________________________________________________________________________

// S-Mult: p = p + m*a at order p->order
template<typename T> Ofts<T>& Ofts<T>::smult(Ofts<T> const& a, T const& m)
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
            term[k]->smult(*a.term[k], m);  //ps[k] += m*a[k]
        }

        return *this;
    }

}

//---------------------------------------

// Mult: p = p + m*a at order p->order
template<typename T> Ofts<T>& Ofts<T>::mult(Ofts<T> const& a, T const& m)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
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

//---------------------------------------

// S-Mult: p = p + c*a at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::smult(Ofts< Ots<U> > const& a, U const& c)
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
            term[k]->smult(*a.term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}
// S-Mult: p = p + c*a at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::smult(Ofts< Ofs<U> > const& a, U const& c)
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
            term[k]->smult(*a.term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}



// S-Mult: p = p + c*a at order k with c of type <typename U>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::smult(Ofts< Ots<U> > const& a, U const& c, int const& k)
{
    if(k > order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->smult(*a.term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}
// S-Mult: p = p + c*a at order k with c of type <typename U>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::smult(Ofts< Ofs<U> > const& a, U const& c, int const& k)
{
    if(k > order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->smult(*a.term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}



// S-Mult: p = c*a at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::mult(Ofts< Ots<U> > const& a, U const& c)
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
            term[k]->mult(*a.term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}
// S-Mult: p = c*a at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::mult(Ofts< Ofs<U> > const& a, U const& c)
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
            term[k]->mult(*a.term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}



// S-Mult: p = p + c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::smult(Ofts< Ots<U> > const& a, T const& m, U const& c)
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
            term[k]->smult(*a.term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}
// S-Mult: p = p + c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::smult(Ofts< Ofs<U> > const& a, T const& m, U const& c)
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
            term[k]->smult(*a.term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

// S-Mult: p = p + c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::smult(Ofts< Ots<U> > const& a, T const& m, U const& c, int const& k)
{
    if(k > order  || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->smult(*a.term[k], m, c);  //ps[k] += m*a[k]
        return *this;
    }
}
// S-Mult: p = p + c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::smult(Ofts< Ofs<U> > const& a, T const& m, U const& c, int const& k)
{
    if(k > order  || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        term[k]->smult(*a.term[k], m, c);  //ps[k] += m*a[k]
        return *this;
    }
}

// S-Mult: p = c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ots<U> >& Ofts<T>::mult(Ofts< Ots<U> > const& a, T const& m, U const& c)
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
// S-Mult: p = c*m*a at order p->order with c of type <typename U> and m of type <typename T>
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::mult(Ofts< Ofs<U> > const& a, T const& m, U const& c)
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


//---------------------------------------

// Product: p = p + a*b at order p->order
template<typename T> Ofts<T>& Ofts<T>::sprod(Ofts<T> const& a, Ofts<T> const& b)
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


// Product: p = p + a*b only at order n
// Handle the case for which n >= max(a.order, b.order)
template<typename T> Ofts<T>& Ofts<T>::sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i, i0, i1;

    i0 = min(b.order, n);
    i1 = min(a.order, n);

    //Product
    for(i= n-i0; i<=i1; i++) term[n]->sprod(*a.term[i], *b.term[n-i]);

    return *this;
}


// Product: p = p + a*b at order p->order
template<typename T> Ofts<T>& Ofts<T>::prod(Ofts<T> const& a, Ofts<T> const& b)
{
    this->zero();
    this->sprod(a,b);

    return *this;
}


// Product: p = p + m*a*b at order p->order
template<typename T> Ofts<T>& Ofts<T>::smprod(Ofts<T> const& a, Ofts<T> const& b, T const& m)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
    //temp = a*b
    temp.sprod(a,b);
    //this = m*temp
    this->smult(temp,m);
    return *this;
}

// Product: p = m*a*b at order p->order
template<typename T> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, T const& m)
{
    this->zero();
    this->smprod(a,b,m);
    return *this;
}


// Product: p = p + c*a*b at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts<T>& Ofts<T>::smprod(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->nv, this->order, this->cnv, this->corder);
    //temp = a*b
    temp.sprod(a,b);
    //this = m*temp
    this->smult(temp,c);
    return *this;
}

// Product: p = p + c*a*b at order p->order with c of type <typename U>
template<typename T> template<typename U> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    this->zero();
    this->smprod(a,b,c);
    return *this;
}


//---------------------------------------


// Factored Sum: p = ma*a+mb*b at order p->order
template<typename T> Ofts<T>& Ofts<T>::fsum(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)
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
            term[k]-> mult(*a.term[k], ma);   //ps[k]   = ma*a[k], contains the zeroing
            term[k]->smult(*b.term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}


//____________________________________________________________________________________________________________________
//
// END OF VALIDITY OF Ofs<double> implementation (no division, no power...)
//____________________________________________________________________________________________________________________

// Division of taylor functions: p = a/b
template<typename T>  Ofts<T>& Ofts<T>::divs(Ofts<T> const& a, Ofts<T> const& b)
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

// Double case
template<>  inline Ofts< Ots<double> >& Ofts< Ots<double> >::divs(Ofts< Ots<double> > const& a, Ofts< Ots<double> > const& b)
{
    int k, j;

    //Order 0
    this->acoef0s(*coef0s(a) / *coef0s(b));

    //Temporary Ots
    Ots< double > ots_temp(this->cnv, this->corder);
    //Temporary Ofts
    Ofts< Ots<double> > ofts_temp(this->nv, this->order, this->cnv, this->corder);

    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //sets every coefficients to zero for order k
        this->term[k]->zero();
        //Copy of a[k] in p[k]
        this->term[k]->smult(*a.term[k], 1.0);
        for(j=0; j<= k-1; j++)  this->term[k]->smprod(*b.term[k-j], *this->term[j],  -1.0);
        //ots_temp = 1/coef0s(b);
        ots_temp.pows(*coef0s(b), -1.0);
        //ofts_temp.term[k] += ots_temp*this->term[k]
        ofts_temp.term[k]->smult(*this->term[k], ots_temp);
        //copy ofts_temp.term[k] in this->term[k]
        this->term[k]->ccopy(*ofts_temp.term[k]);
    }

    return *this;
}

// Double complex case
template<>  inline Ofts< Ots<double complex> >& Ofts< Ots<double complex> >::divs(Ofts< Ots<double complex> > const& a, Ofts< Ots<double complex> > const& b)
{
    int k, j;

    //Order 0
    this->acoef0s(*coef0s(a) / *coef0s(b));

    //Temporary Ots
    Ots< double complex > ots_temp(this->cnv, this->corder);
    //Temporary Ofts
    Ofts< Ots<double complex> > ofts_temp(this->nv, this->order, this->cnv, this->corder);

    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //sets every coefficients to zero for order k
        this->term[k]->zero();
        //Copy of a[k] in p[k]
        this->term[k]->smult(*a.term[k], (double complex) 1.0);
        for(j=0; j<= k-1; j++)  this->term[k]->smprod(*b.term[k-j], *this->term[j],  (double complex) -1.0);
        //ots_temp = 1/coef0s(b);
        ots_temp.pows(*coef0s(b), -1.0);
        //ofts_temp.term[k] += ots_temp*this->term[k]
        ofts_temp.term[k]->smult(*this->term[k], ots_temp);
        //copy ofts_temp.term[k] in this->term[k]
        this->term[k]->ccopy(*ofts_temp.term[k]);
    }

    return *this;
}



//---------------------------------------

// Power function: p = a^alpha
template<typename T> template<typename U> Ofts<T>& Ofts<T>::pows(Ofts<T> const& a,  U const& alpha)
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

// Double case
template<> template<typename U> Ofts< Ots<double> >& Ofts< Ots<double> >::pows(Ofts< Ots<double> > const& a,   U const& alpha)
{
    int k, j;

    //Temporary Ots
    Ots< double > ots_temp(this->cnv, this->corder);

    //Initialization of order zero
    coef0s(*this)->pows(*coef0s(a), alpha);
    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        for(j=0; j<= k-1; j++) this->term[k]->smprod(*a.term[k-j], *this->term[j], *coef0s(a), (alpha*(k-j)-j)/k, ots_temp);
    }
    return *this;
}

// Double complex case
template<> template<typename U> Ofts< Ots<double complex> >& Ofts< Ots<double complex> >::pows(Ofts< Ots<double complex> > const& a,   U const& alpha)
{
    int k, j;

    //Temporary Ots
    Ots< double complex > ots_temp(this->cnv, this->corder);

    //Initialization of order zero
    coef0s(*this)->pows(*coef0s(a), alpha);
    //Recurrence scheme
    for(k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        for(j=0; j<= k-1; j++) this->term[k]->smprod(*a.term[k-j], *this->term[j], *coef0s(a), (complex double) (alpha*(k-j)-j)/k, ots_temp);
    }
    return *this;
}


// Ofs<double> case ONLY FOR a[0] = 1.0 (unitary Fourier coefficient)
template<> template<typename U> Ofts< Ofs<double> >& Ofts< Ofs<double> >::pows(Ofts< Ofs<double> > const& a,   U const& alpha)
{
    //Initialization of order zero
    this->acoef0s(*coef0s(a));  //a[0]^alpha = 1.0^alpha = 1;
    //Recurrence scheme
    for(int k=1; k <= order; k++)
    {
        //Sets every coefficients to zero for order k
        this->term[k]->zero();
        for(int j=0; j<= k-1; j++) this->term[k]->smprod(*a.term[k-j], *this->term[j], (double) (alpha*(k-j)-j)/k);

    }
    return *this;
}


// Power function: p = a^alpha @order n
template<typename T> template<typename U> Ofts<T>& Ofts<T>::pows(Ofts<T> const& a,  U const& alpha, int const& n)
{
    T x0;
    x0 = coef0s(a);
    if(n==0)
    {
        this->acoef0s(cpow(x0, alpha));
        return *this;

    }
    else
    {
        //Sets every coefficients to zero for order n
        this->term[n]->zero();
        //Recurrence scheme @order n
        int i0 = min(a.order, n);
        for(int i= n-i0; i< n; i++) this->term[n]->smprod(*a.term[n-i], *this->term[i], alpha*(n-i)-i);// smprodh(ps->term[k], s->term[k-j], ps->term[j], alpha*(k-j)-j);
        this->term[n]->mult(1.0/(n*x0)); //multh(ps->term[k], 1.0/(x0*k));
        return *this;
    }
}

// Double case
template<> template<typename U> Ofts< Ots<double> >& Ofts< Ots<double> >::pows(Ofts< Ots<double> > const& a,   U const& alpha, int const& n)
{
    if(n==0)
    {
        //Initialization of order zero
        coef0s(*this)->pows(*coef0s(a), alpha);
        return *this;
    }
    else
    {
        //Temporary Ots
        Ots< double > ots_temp(this->cnv, this->corder);
        //Sets every coefficients to zero for order k
        this->term[n]->zero();
        int i0 = min(a.order, n);
        for(int i= n-i0; i< n; i++) this->term[n]->smprod(*a.term[n-i], *this->term[i], *coef0s(a), (double) (alpha*(n-i)-i)/n, ots_temp);
        return *this;
    }
}


// Double complex case
template<> template<typename U> Ofts< Ots<double complex> >& Ofts< Ots<double complex> >::pows(Ofts< Ots<double complex> > const& a,   U const& alpha, int const& n)
{

    if(n==0)
    {
        //Initialization of order zero
        coef0s(*this)->pows(*coef0s(a), alpha);
        return *this;

    }

    else
    {
        //Temporary Ots
        Ots< double complex > ots_temp(this->cnv, this->corder);
        //Initialization of order zero
        coef0s(*this)->pows(*coef0s(a), alpha);
        //Sets every coefficients to zero for order k
        this->term[n]->zero();
        int i0 = min(a.order, n);
        for(int i= n-i0; i< n; i++) this->term[n]->smprod(*a.term[n-i], *this->term[i], *coef0s(a), (double complex) (alpha*(n-i)-i)/n, ots_temp);
        return *this;
    }
}



// Ofs<double> case ONLY FOR a[0] = 1.0 (unitary Fourier coefficient)
template<> template<typename U> Ofts< Ofs<double> >& Ofts< Ofs<double> >::pows(Ofts< Ofs<double> > const& a,   U const& alpha, int const& n)
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
        for(int i= n-i0; i< n; i++) this->term[n]->smprod(*a.term[n-i], *this->term[i], (double) (alpha*(n-i)-i)/n);
        return *this;
    }
}

//---------------------------------------

/*
    Returns the address of the first coefficient of order 0 of the taylor serie s
*/
template<typename T> T* Ofts<T>::coef0s(Ofts<T> const& a)
{
    return a.term[0][0].getCoefAddress();
}

/*
   Sets the coefficient of order 0 of the taylor serie s equal to x0
*/
template<typename T> void Ofts<T>::acoef0s(T const& x0)
{
    this->term[0][0].setCoef(x0,0);
}



//Stream
//____________________________________________________________________________________________________________________


template<typename T> std::ostream& operator << (std::ostream& stream, Ofts<T> const& ofts)
{
    int i,j, nrc;
    int k[ofts.nv];

    stream << "#Taylor serie"             << endl;
    stream << "#Degree: " << ofts.order    << endl;
    stream << "#Variables: " << ofts.nv    << endl;
    stream << "--------------------------"<< endl;

    for(nrc=0; nrc<= ofts.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ofts.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ofts.nv, nrc); i++)
        {
            for(j=0; j<ofts.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofts.term[nrc]->getCoef(i) << std::noshowpos << endl;

            if(i< FTDA::nmon(ofts.nv, nrc)-1)  FTDA::prxkt(k, ofts.nv);
        }
    }
    return stream;
}


template<> inline std::ostream& operator << (std::ostream& stream, Ofts< Ots <complex double> > const& ofts)
{
    int i,j, nrc;
    int k[ofts.nv];

    stream << "#Taylor serie"             << endl;
    stream << "#Degree: " << ofts.order    << endl;
    stream << "#Variables: " << ofts.nv    << endl;
    stream << "--------------------------"<< endl;

    Ofs<complex double> ofs(ofts.corder/2);

    for(nrc=0; nrc<= ofts.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ofts.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ofts.nv, nrc); i++)
        {
            ofs.ts2fs(ofts.term[nrc]->getCoef(i));
            for(j=0; j<ofts.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofs << std::noshowpos << endl;
            //stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofts.term[nrc]->getCoef(i) << std::noshowpos << endl;

            if(i< FTDA::nmon(ofts.nv, nrc)-1)  FTDA::prxkt(k, ofts.nv);
        }
    }
    return stream;
}


template<> inline std::ostream& operator << (std::ostream& stream, Ofts< Ots <double> > const& ofts)
{
    int i,j, nrc;
    int k[ofts.nv];

    stream << "#Taylor serie"             << endl;
    stream << "#Degree: " << ofts.order    << endl;
    stream << "#Variables: " << ofts.nv    << endl;
    stream << "--------------------------"<< endl;

    Ofs<double> ofs(ofts.corder/2);

    for(nrc=0; nrc<= ofts.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ofts.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(ofts.nv, nrc); i++)
        {
            ofs.ts2fs(ofts.term[nrc]->getCoef(i));
            for(j=0; j<ofts.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofs << std::noshowpos << endl;
            //stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofts.term[nrc]->getCoef(i) << std::noshowpos << endl;

            if(i< FTDA::nmon(ofts.nv, nrc)-1)  FTDA::prxkt(k, ofts.nv);
        }
    }
    return stream;
}


// fts 2 fs
//____________________________________________________________________________________________________________________
/*
    fts 2 fs for fts of one variable only
*/
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
            tempc += tempfs.getCoef(k)*pow(epsilon, (double) l);
        }
        fs->setCoef(tempc, k);
    }
}

/*
    fts 2 fs for fts of one variable only
*/
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
            tempc += tempfs.getCoef(k)*pow(epsilon, (double) l);
        }
        fs->setCoef(tempc, k);
    }
}

// fts 2 scalar
//____________________________________________________________________________________________________________________
template <typename U>  double complex fts2scalar(Ofts< Ots<U> > const& fts_z, double epsilon, double t)
{
    Ofs<U> fs(fts_z.getCOrder()/2);
    fts2fs(&fs, fts_z, epsilon);
    return fs.evaluate(t);
}

template <typename U>  double complex fts2scalar(Ofts< Ofs<U> > const& fts_z, double epsilon, double t)
{
    Ofs<U> fs(fts_z.getCOrder());
    fts2fs(&fs, fts_z, epsilon);
    return fs.evaluate(t);
}


