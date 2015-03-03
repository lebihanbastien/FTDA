//############################################################################
// Implementation of the Oftsh template class
//############################################################################

template<typename T> int Oftsh<T>::globalOrder = 20;
template<typename T> int Oftsh<T>::globalVariable = 3;


//Create
//---------------------------------------------------------------------------

/*
    Default constructor (should not be used anywhere else than in this file)
*/
template<typename T> Oftsh<T>::Oftsh()
{
    order = 0;
    nv = 0;
    coef = new T(0);
}

/*
   Allocates memory  without any link to a coefficient array (requires the inmediate used of setCoefs afterwards).
*/
template<typename T> Oftsh<T>::Oftsh(int newNv, int newOrder)
{
    order = newOrder;
    nv = newNv;
    coef = 0;//new T(0);//[ FTDA::nmon(nv, order)]();

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
            term[i].lcopy(Oftsh<T>(newNv-1, newOrder-i));
            //At this step, all sons of the current Oftsh<T> have nv = nv-1, order = order-i and coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
        }
    }
}

/*
   Performs linking between the Oftsh<T> and a coefficient array

   parameters:
   h: the polynomial to link;
   coef0: the array of coefficients;
*/
template<typename T> void Oftsh<T>::setCoefs(T *coef0)
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
        for(i=0, coefc=coef0; i<=order; coefc+= FTDA::nmon(nv-1,order-i), i++)
        {
            term[i].setCoefs(coefc);
        }
    }
}



//Copy
//---------------------------------------------------------------------------
/*
    Construction without any link
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
    coef = new T(0);//[ FTDA::nmon(nv, order)]();
    //Tree
    if(nv==1)  //1-variable polynomial
    {
        term = new Oftsh<T>[order+1];
        if (term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        term = new Oftsh<T>[order+1];  //At this step, all sons of the current Oftsh<T> have nv = 0, order = 0 and coef points @ one T.
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

    //Copy in linking of the coefficients
    T *coef0 = new T[ FTDA::nmon(b.nv, b.order)]();
    for(int i = 0 ; i <  FTDA::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
    this->setCoefs(coef0);
}

/*
    Copy without any link (entirely new object!)
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
        coef = new T(0);//[ FTDA::nmon(nv, order)]();
        //Tree
        if(nv==1)  //1-variable polynomial
        {
            term = new Oftsh<T>[order+1];
            if (term == NULL)
            {
                puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
                exit(1);
            }
        }
        else
        {
            int i;
            term = new Oftsh<T>[order+1];  //At this step, all sons of the current Oftsh<T> have nv = 0, order = 0 and coef points @ one T.
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
        //Copy in linking of the coefficients
        T *coef0 = new T[ FTDA::nmon(b.nv, b.order)]();
        for(int i = 0 ; i <  FTDA::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
        this->setCoefs(coef0);
    }
    return *this; //same object if returned
}

/*
    Linked copy
*/
template<typename T> Oftsh<T>& Oftsh<T>::lcopy (Oftsh<T> const& b)
{
    order = b.order;
    nv = b.nv;
    term = b.term;
    coef = b.coef;
    return *this;
}

/*
    Coefficient copy (restricted to same order, same number of variables)
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
        for(int i = 0 ; i< FTDA::nmon(nv, order) ; i++) coef[i] = b.coef[i];
        return *this;
    }
}


//Delete       TO BE DETERMINED: how properly delete with recursivity?
//              Seems to work fine like this, but memory leak?
//---------------------------------------------------------------------------
template<typename T> Oftsh<T>::~Oftsh<T>()
{
    /*if(order == 0)
    {
        delete term;
    }*/
}

//Zeroing
//---------------------------------------------------------------------------
template<typename T> void Oftsh<T>::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->setCoef(0.0, i);
}


// Ofs<double> case
template<> inline void Oftsh< Ofs<double> >::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->coef[i].zero();
}


// Ots<double> case
template<> inline void Oftsh< Ots<double> >::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->coef[i].zero();
}


/*
    Ots<double complex> case
*/
template<> inline void Oftsh< Ots< complex double> >::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->coef[i].zero();
}

//Setters
//---------------------------------------------------------------------------
template<typename T> void Oftsh<T>::setCoef(T value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos] = value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

template<typename T> template<typename U> void Oftsh< T >::setTCoef(U value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos].setCoefs(value);
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

template<typename T> template<typename U> void Oftsh< T >::setT0Coef(U value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos].setCoef(value, 0);
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

template<typename T> void Oftsh<T>::addCoef(T value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos] += value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

// Ots<double> case
template<> inline void Oftsh< Ots<double> >::setCoefs( Ots<double> *coef0)
{
    if(nv==1) //1-variable homogeneous polynomial
    {
        coef = coef0;
    }
    else
    {
         Ots<double>  *coefc;
        int i;
        //allocate the position of the first coefficient
        coef = coef0;
        //perform recursice allocation
        for(i=0, coefc=coef0; i<=order; coefc+= FTDA::nmon(nv-1,order-i), i++)
        {
            term[i].setCoefs(coefc);
        }
    }
}


template<typename T> void Oftsh<T>::setRandomCoefs()
{
    for(int pos = 0; pos< FTDA::nmon(nv, order); pos++) coef[pos].setRandomCoefs();
}

//Conjugate
template<typename T> Oftsh<T>& Oftsh<T>::conjugate()
{
     for(int pos=0; pos< FTDA::nmon(nv, order); pos++) coef[pos] = conj(coef[pos]);

     return *this;
}

// Ots<double complex> case
template<> inline Oftsh< Ots<double complex> >& Oftsh< Ots<double complex> >::conjugate()
{
     for(int pos=0; pos< FTDA::nmon(nv, order); pos++)
     {
        coef[pos].conjugate();
     }
     return *this;
}

// Ots<double> case
template<> inline Oftsh< Ots<double> >& Oftsh< Ots<double> >::conjugate()
{
     for(int pos=0; pos< FTDA::nmon(nv, order); pos++)
     {
        coef[pos].conjugate();
     }
     return *this;
}


// Ofs<double> case
template<> inline Oftsh< Ofs<double> >& Oftsh< Ofs<double> >::conjugate()
{
     for(int pos=0; pos< FTDA::nmon(nv, order); pos++)
     {
        coef[pos].conjugate();
     }
     return *this;
}

//Getters
//---------------------------------------------------------------------------
template<typename T> Oftsh<T> Oftsh<T>::getTerm()
{
    return term[0];
}

template<typename T> T* Oftsh<T>::getCoefAddress()
{
    return coef;
}

template<typename T> Oftsh<T> Oftsh<T>::getTerm(int i) const
{
    if(i >=0 && i <= order) return term[i];
    else
    {
        cout << "Error in getTerm(int i): position is out of scope." << endl;
        cout << "First term is returned" << endl;
        return term[0];
    }
}

/*
    Double && Complex case
    If position is out of scope, 0.0 is returned by default
*/
template<typename T>  T Oftsh<T>::getCoef(int i)
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else return 0.0;
}


//Ots<double> case
template<>  inline Ots<double> Oftsh< Ots<double> >::getCoef(int i)
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else
    {
        cout << "Warning in Oftsh< Ots<double> >::getCoef: required position is out of scope, null coefficient is returned." << endl;
        int nvc, nrc;
        nvc = coef[0].getVariables();
        nrc = coef[0].getOrder();
        Ots<double> ots(nvc, nrc);
        return ots;
    }
}

//Ots<complex double> case
template<> inline Ots<complex double> Oftsh< Ots<complex double> >::getCoef(int i)
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else
    {
        cout << "Warning in Oftsh< Ots<complex> >::getCoef: required position is out of scope, null coefficient is returned." << endl;
        int nvc, nrc;
        nvc = coef[0].getVariables();
        nrc = coef[0].getOrder();
        Ots<complex double> ots(nvc, nrc);
        return ots;
    }
}


template<typename T> int Oftsh<T>::getOrder()
{
    return order;
}

template<typename T> int Oftsh<T>::getVariables()
{
    return nv;
}



//Operators
//---------------------------------------------------------------------------
/*
template<typename T> bool Oftsh<T>::isEqual(Oftsh<T> const& b) const
{
    if(order != b.order || nv != b.nv) return false;
    else
    {
        bool result = true;
        for(int i = 0 ; i<=  FTDA::nmon(nv, order); i++) result = result&&(coef[i] == b.coef[i]);
        return result;
    }
}

template<typename T> Oftsh<T>& Oftsh<T>::operator += (Oftsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using += : the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] += b.coef[i];
        return *this;
    }
}

template<typename T> Oftsh<T>& Oftsh<T>::operator -= (Oftsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using += : the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] -= b.coef[i];
        return *this;
    }
}

template<typename T> Oftsh<T>& Oftsh<T>::operator *= (T const& c)
{
    for(int i=0; i <  FTDA::nmon(nv, order); i++) coef[i] *= c;
    return *this;
}

template<typename T> Oftsh<T>& Oftsh<T>::operator /= (T const& c)
{
    for(int i=0; i <  FTDA::nmon(nv, order); i++) coef[i] /= c;
    return *this;
}

*/


/*
   Performs the scalar-polynomial multiplication + sum h = h + r*a

   parameters:
   h: the polynomial output
   r: the scalar
   a: the polynomial input
*/
template<typename T> Oftsh<T>& Oftsh<T>::smult(Oftsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] += c*a.coef[i];
        return *this;
    }
}

//Double case
template<> inline Oftsh< Ots<double> >& Oftsh< Ots<double> >::smult(Oftsh< Ots<double> > const& a,  Ots<double>  const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << order << " " << a.order << endl;
        cout << "Error using smult: the order and/or number of variables does not match (here). Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].sprod(a.coef[i], c);
        return *this;
    }
}
//Double complex case
template<> inline Oftsh< Ots<double complex> >& Oftsh< Ots<double complex> >::smult(Oftsh< Ots<double complex> > const& a,  Ots<double complex>  const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << order << " " << a.order << endl;
        cout << "Error using smult: the order and/or number of variables does not match (here). Initial Oftsh< Ots<double complex> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].sprod(a.coef[i], c);
        return *this;
    }
}


// Case with c is <typename U>
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::smult(Oftsh< Ots<U> > const& a, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].smult(a.coef[i], c);
        return *this;
    }
}
// Case with c is <typename U>
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::smult(Oftsh< Ofs<U> > const& a, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].smult(a.coef[i], c);
        return *this;
    }
}



// Case with c is <typename U>
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::mult(Oftsh< Ots<U> > const& a, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].mult(a.coef[i], c);
        return *this;
    }
}
// Case with c is <typename U>
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::mult(Oftsh< Ofs<U> > const& a, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].mult(a.coef[i], c);
        return *this;
    }
}


// Case with c is <typename U> && Ots<U>
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::smult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].smprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}
// Case with c is <typename U> && Ots<U>
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::smult(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].smprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}





template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::mult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].mprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::mult(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
        if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].mprod(a.coef[i], ra, c);  //sprod(a.coef[i], c*ra)
        return *this;
    }
}


// Case with c is <typename U> && Ots<U>
//void smulthdiv_tsc(fth *h, fth *a, taylor_serie r1, taylor_serie r2, complex m, taylor_serie *temp)
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::smult(Oftsh< Ots<U> > const& a, Ots<U> const& ra, Ots<U> const& r, U const& c,  Ots<U> & temp)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++)
        {
            //temp = a.coef[i] / r
            temp.mdivs(a.coef[i], r, 1.0);
            //coef[i]+= m*ra*temp
            coef[i].smprod(temp, ra, c);
        }
        return *this;
    }
}


/*
   Performs the scalar-polynomial multiplication  h = r*a

   parameters:
   h: the polynomial output
   r: the scalar
   a: the polynomial input
*/
template<typename T> Oftsh<T>& Oftsh<T>::mult(Oftsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] = c*a.coef[i];
        return *this;
    }
}

//Double case
template<> inline Oftsh< Ots<double> >& Oftsh< Ots<double> >::mult(Oftsh< Ots<double> > const& a,  Ots<double>  const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].prod(a.coef[i], c);
        return *this;
    }
}

//Double complex case
template<> inline Oftsh< Ots<double complex> >& Oftsh< Ots<double complex> >::mult(Oftsh< Ots<double complex> > const& a,  Ots<double complex>  const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ots<double complex> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i].prod(a.coef[i], c);
        return *this;
    }
}


/*
   Performs the homogeneous polynomial-polynomial multiplication + sum: p = p + a*b

   parameters:
   p: the polynomial output
   a: the polynomial input
   b: the polynomial input
*/
template<typename T> Oftsh<T>& Oftsh<T>::sprod(Oftsh<T> const& a, Oftsh<T> const& b)
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
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->sprod(*aa , *bb);
                    pp->smult(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef)), bb++ , pp++);
                *(pp->coef)+= *(aa->coef) * *(bb->coef);
            }
            else this->smult(a, *(b.coef)); //b is scalar
        }
        else this->smult(b, *(a.coef));   //a is scalar
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


// Ofs<double> case
template<> inline Oftsh< Ofs<double> >& Oftsh< Ofs<double> >::sprod(Oftsh< Ofs<double> > const& a, Oftsh< Ofs<double> > const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ofs<double> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->sprod(*aa , *bb);
                    pp->smult(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef)), bb++ , pp++);
                pp->coef->sprod(*(aa->coef), *(bb->coef));
            }
            else this->smult(a, *(b.coef)); //b is scalar
        }
        else this->smult(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
         Ofs<double>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->sprod(*aa, *bb);
    }
    else
    {
        this->coef->sprod(*(a.coef ), *(b.coef));//1-variate homogeneous polynomial
    }
    return *this;
}

// Ots<double> case
template<> inline Oftsh< Ots<double> >& Oftsh< Ots<double> >::sprod(Oftsh< Ots<double> > const& a, Oftsh< Ots<double> > const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ots<double> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->sprod(*aa , *bb);
                    pp->smult(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef)), bb++ , pp++);
                pp->coef->sprod(*(aa->coef), *(bb->coef));
                //*(pp->coef)+= *(aa->coef) * *(bb->coef);
            }
            else this->smult(a, *(b.coef)); //b is scalar
        }
        else this->smult(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
         Ots<double>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->sprod(*aa, *bb);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->sprod(*(a.coef ), *(b.coef));//*(this->coef)+= *(a.coef ) * *(b.coef); //1-variate homogeneous polynomial
    }
    return *this;
}

// Ots<double complex> case
template<> inline Oftsh< Ots<double complex> >& Oftsh< Ots<double complex> >::sprod(Oftsh< Ots<double complex> > const& a, Oftsh< Ots<double complex> > const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ots<double complex> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->sprod(*aa , *bb);
                    pp->smult(*aa, *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef)), bb++ , pp++);
                pp->coef->sprod(*(aa->coef), *(bb->coef));
                //*(pp->coef)+= *(aa->coef) * *(bb->coef);
            }
            else this->smult(a, *(b.coef)); //b is scalar
        }
        else this->smult(b, *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
         Ots<double complex>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->sprod(*aa, *bb);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->sprod(*(a.coef ), *(b.coef));//*(this->coef)+= *(a.coef ) * *(b.coef); //1-variate homogeneous polynomial
    }
    return *this;
}


// p = p + m*a*b
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::smprod(Oftsh< Ots<U> > const& a, Oftsh< Ots<U> > const& b, U const& m)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ots<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->smprod(*aa , *bb, m);
                    pp->smult(*aa, *(bb->coef), m);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef), m), bb++ , pp++);
                pp->coef->smprod(*(aa->coef), *(bb->coef),m);
            }
            else this->smult(a, *(b.coef),m); //b is scalar
        }
        else this->smult(b, *(a.coef),m);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        Ots<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->smprod(*aa, *bb,m);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->smprod(*(a.coef ), *(b.coef),m); //1-variate homogeneous polynomial
    }
    return *this;
}


// Ofs<double> case
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::smprod(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)
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
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->smprod(*aa , *bb, m);
                    pp->smult(*aa, *(bb->coef), m);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef), m), bb++ , pp++);
                pp->coef->smprod(*(aa->coef), *(bb->coef),m);
            }
            else this->smult(a, *(b.coef),m); //b is scalar
        }
        else this->smult(b, *(a.coef),m);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        Ofs<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->smprod(*aa, *bb,m);//*pp+= *aa * *bb;
    }
    else
    {
        this->coef->smprod(*(a.coef ), *(b.coef),m); //1-variate homogeneous polynomial
    }
    return *this;
}

/*
    p = p + m*a*b/r
*/
template<typename T> template<typename U> Oftsh< Ots<U> >& Oftsh<T>::smprod(Oftsh< Ots<U> > const& a, Oftsh< Ots<U> > const& b, Ots<U> const& r, U const& m, Ots<U> & temp)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Oftsh< Ots<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->smprod(*aa , *bb, r, m, temp);
                    pp->smult(*aa, *(bb->coef), r, m, temp);
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, *(aa->coef), r, m, temp), bb++ , pp++)
                //temp = aa->coef/r
                temp.mdivs(*(aa->coef), r, 1.0);
                //pp->coef += m*bb->coef*temp
                pp->coef->smprod(temp, *(bb->coef), m);
            }
            else this->smult(a, *(b.coef), r, m, temp); //b is scalar
        }
        else this->smult(b, *(a.coef), r, m, temp);   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        Ots<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )
            {
                //temp = aa/r
                temp.mdivs(*aa, r, 1.0);
                //pp += m*bb*temp
                pp->smprod(temp, *bb, m);
            }
    }
    else
    {
        //temp = a/r
        temp.mdivs(*a.coef, r, 1.0);
        //this-Coef += m*bb*temp
        this->coef->smprod(temp, *b.coef, m);
    }
    return *this;
}


// Functions
//---------------------------------------------------------------------------
template<typename T> bool operator==(Oftsh<T> const& a, Oftsh<T> const& b)
{
    return a.isEqual(b);
}

/*
    Returns a+b
*/
template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop+=b;
    return cop;
}

/*
    Returns a-b
*/
template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop-=b;
    return cop;
}

//Stream
//---------------------------------------------------------------------------

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

    for (i=0; i< FTDA::nmon(oftsh.nv, oftsh.order); i++)
    {
        for(j=0; j<oftsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  oftsh.coef[i] << std::noshowpos << endl;

        if(i< FTDA::nmon(oftsh.nv, oftsh.order)-1)  FTDA::prxkt(k, oftsh.nv);
    }
    return stream;
}

//Complex version
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

    for (i=0; i< FTDA::nmon(oftsh.nv, oftsh.order); i++)
    {
        for(j=0; j<oftsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  creal(oftsh.coef[i]) << " " <<  cimag(oftsh.coef[i]) << std::noshowpos << endl;

        if(i< FTDA::nmon(oftsh.nv, oftsh.order)-1)  FTDA::prxkt(k, oftsh.nv);
    }
    return stream;
}

template<> inline std::ostream& operator << (std::ostream& stream, Oftsh< Ots <double> > const& oftsh)
{
    int i,j, nrc;
    int k[oftsh.nv];

    stream << "#Homogeneous polynomial"             << endl;
    stream << "#Degree: " << oftsh.order    << endl;
    stream << "#Variables: " << oftsh.nv    << endl;
    stream << "--------------------------"<< endl;

    Ofs<double> ofs(oftsh.coef[0].getVariables(), oftsh.coef[0].getOrder()/2);

    cout << "oftsh.coef[0].getOrder()/2: " << oftsh.coef[0].getOrder()/2 << endl;

    nrc = oftsh.order;
        //Update k
        k[0] = nrc;
        for(i=1; i<oftsh.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< FTDA::nmon(oftsh.nv, nrc); i++)
        {
            ofs.ts2fs(oftsh.coef[i]);
            for(j=0; j<oftsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  ofs << std::noshowpos << endl;
            //stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  oftsh.coef[i] << std::noshowpos << endl;

            if(i< FTDA::nmon(oftsh.nv, nrc)-1)  FTDA::prxkt(k, oftsh.nv);
        }
    return stream;
}
