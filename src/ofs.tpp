//############################################################################
// Implementation of the Ofs template class
//############################################################################

template <typename T>
int Ofs<T>::globalOrder = 20;

/*
    The Ofs class computes the operations of fourier series. Note that the array containing the coefficients
    is indexed from 0 to 2*order of the Fourier series, with the positions 0 to order-1 containing the terms
    of order -order to -1, and the positions order to 2*order containing the terms of order 0 to order.
    The getters and the setters should be the ONLY routines to take into account this shift.
*/
//Create
//---------------------------------------------------------------------------
template <typename T> Ofs<T>::Ofs()
{
    cout << "An Ofs object was created without any selected order:" << endl;
    cout << "Order " << Ofs::globalOrder << " is used by default." << endl;
    order = Ofs::globalOrder;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
}

template <typename T>
Ofs<T>::Ofs(const int newOrder)
{
    order = newOrder;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
}

template <typename T> Ofs<T>::Ofs(const int newOrder, T coef0)
{
    order = newOrder;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
    this->setCoef(coef0,0);
}

template <typename T> Ofs<T>::Ofs(const int newOrder, T coefn, int n)
{
    order = newOrder;
    coef = new T[2*order+1];
    for(int i = 0 ; i< 2*order + 1; i++) coef[i] = 0.0;
    if(fabs(n) <= order) this->setCoef(coefn, n);
    else cout << "Error in constructor: position is out of scope\n No coefficient is set (all zeros).\n" << endl;
}

//Copy
//---------------------------------------------------------------------------
template <typename T> Ofs<T>::Ofs(Ofs const& b)
{
    order = b.order;
    coef = new T[2*order+1];
    for(int i = -order ; i<= order; i++) this->setCoef(b.getCoef(i), i);
}

template <typename T> Ofs<T>& Ofs<T>::operator = (Ofs<T> const& b)
{
    if(this != &b)
    {
        order = b.order;
        delete coef;
        coef = new T[2*order+1];
        for(int i = -order ; i<= order; i++) this->setCoef(b.getCoef(i), i);
    }
    return *this; //same object if returned
}

template <typename T> Ofs<T>& Ofs<T>::ccopy(Ofs<T> const& b)
{
    if(order != b.order)
    {
        cout << "Erreur in ccopy for Ofs: orders do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {

    for(int i = -order ; i<= order; i++) this->setCoef(b.getCoef(i), i);
    return *this; //same object if returned
    }
}

//Single coefficients
template <typename T> Ofs<T>& Ofs<T>::operator  = (T const& coef0)
{
    order = 0;//Ofs::globalOrder;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
    coef[order] = coef0;
    return *this;
}


//Delete
//---------------------------------------------------------------------------
template <typename T> Ofs<T>::~Ofs<T>()
{
    delete coef;
    coef = 0;
}

//Zeroing
//---------------------------------------------------------------------------
template <typename T> void Ofs<T>::zero()
{
    for(int i = -order ; i<= order; i++) this->setCoef(0.0, i);
}

//Setters
//---------------------------------------------------------------------------
template <typename T> void Ofs<T>::setCoef(T value, int pos)
{
    if(fabs(pos) <= order) coef[pos+order] = value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

template <typename T> void Ofs<T>::setCoefs(T value)
{
    for(int pos = -order; pos <= order; pos++) this->setCoef(value, pos);
}

template <typename T> void Ofs<T>::addCoef(T value, int pos)
{
    if(fabs(pos) <= order) coef[pos+order] += value;
    else cout << "Error in addCoef: position is out of scope\n" << endl;
}

template <typename T> void Ofs<T>::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -order; k <=order; k++) this->setCoef(ofs_temp.getCoef(-k), k);
}

template <> inline void Ofs<double complex >::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -order; k <=order; k++) this->setCoef(conj(ofs_temp.getCoef(-k)), k);
}

//Getters
//---------------------------------------------------------------------------

/*
    If position is out of scope, 0.0 is returned by default
*/
template <typename T> T Ofs<T>::getCoef(int pos) const
{
    if(fabs(pos) <= order)  return coef[pos+order];
    else return 0.0;
}

template <typename T> int Ofs<T>::getOrder() const
{
    return order;
}

template <typename T> Ofs<T>* Ofs<T>::getAddress() const
{
    return (Ofs<T>*) this;
}


//Operators
//---------------------------------------------------------------------------
template <typename T> bool Ofs<T>::isEqual(Ofs<T> const& b) const
{
    if(order != b.order) return false;
    else
    {
        bool result = true;
        for(int i = 0 ; i< 2*order + 1; i++) result = result&&(coef[i] == b.coef[i]);
        return result;
    }

}

template <typename T>  Ofs<T>& Ofs<T>::operator += (Ofs<T> const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        T temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new T[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] += b.coef[i+b.order];
    return *this;
}

template <typename T>  Ofs<T>& Ofs<T>::operator -= (Ofs<T> const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        T temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new T[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] -= b.coef[i+b.order];
    return *this;
}

template <typename T>  Ofs<T>& Ofs<T>::operator *= (T const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] *= c;
    return *this;
}

template <typename T>  Ofs<T>& Ofs<T>::operator /= (T const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] /= c;
    return *this;
}


/*
    Returns a*b at order max(a.order,b.order).

    WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in:
                    temp.addCoef(a.getCoef(p)*b.getCoef(n-p), n);

    The getCoef function set these coefficients to zero, which guarantees the good result.
    However, unecessary product are made. psup and pinf must be redefined.
*/
template <typename T>  Ofs<T>  operator * (Ofs<T> const& a, Ofs<T> const& b)
{
    int n, p, tO, psup, pinf;
    tO = max(a.getOrder(),b.getOrder());

    //Virgin ofs of same order
    Ofs<T> temp(tO);

    //Product
    for(n=-tO ; n<= tO; n++)
    {
        psup = min(n+a.getOrder(),  a.getOrder());
        pinf = max(n-b.getOrder(), -b.getOrder());
        for(p=pinf; p<= psup; p++)
        {
            //indix
            temp.addCoef(a.getCoef(p)*b.getCoef(n-p), n);
        }
    }

    return temp;
}

/*
    This = this + a*b with new order = max(order, a.order, b.order): su is truncated at max(a.getOrder(),b.getOrder())

    WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in:
                    temp.addCoef(a.getCoef(p)*b.getCoef(n-p), n);

    The getCoef function set these coefficients to zero, which guarantees the good result.
    However, unecessary product are made. psup and pinf must be redefined.
*/
template <typename T>  Ofs<T>& Ofs<T>::sprod(Ofs<T> const& a, Ofs<T> const& b)
{
    int n, p, tO, psup, pinf;
    tO = max(a.getOrder(),b.getOrder());

    //Product
    for(n=-tO ; n<= tO; n++)
    {
        psup = min(n+a.getOrder(),  a.getOrder());
        pinf = max(n-b.getOrder(), -b.getOrder());
        for(p=pinf; p<= psup; p++)
        {
            //indix
            this->addCoef(a.getCoef(p)*b.getCoef(n-p), n);
        }
    }
    return *this;
}

/*
    This = this + a*b with new order = max(order, a.order, b.order): su is truncated at max(a.getOrder(), b.getOrder())

    WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in:
                    temp.addCoef(a.getCoef(p)*b.getCoef(n-p), n);

    The getCoef function set these coefficients to zero, which guarantees the good result.
    However, unecessary product are made. psup and pinf must be redefined.
*/
template <typename T>  Ofs<T>& Ofs<T>::smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    int n, p, tO, psup, pinf;
    tO = max(a.getOrder(),b.getOrder());

    //Product
    for(n=-tO ; n<= tO; n++)
    {
        psup = min(n+a.getOrder(),  a.getOrder());
        pinf = max(n-b.getOrder(), -b.getOrder());
        for(p=pinf; p<= psup; p++)
        {
            //indix
            this->addCoef(m*a.getCoef(p)*b.getCoef(n-p), n);
        }
    }
    return *this;
}

/*
    This = a*b with new order = max(a.order, b.order).
*/
template <typename T>  Ofs<T>& Ofs<T>::prod(Ofs<T> const& a, Ofs<T> const& b)
{
    this->zero();
    this->sprod(a,b);
    return *this;
}
/*
    This = m*a*b with new order = max(a.order, b.order).
*/
template <typename T>  Ofs<T>& Ofs<T>::mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    this->zero();
    this->smprod(a,b,m);
    return *this;
}

/*
    this|n += [a*b]|n. With new order = max(order, n).
*/
template <typename T>  Ofs<T>& Ofs<T>::sprodn(Ofs<T> const& a, Ofs<T> const& b, int n)
{
    int p, psup, pinf;

    if(n > order) //if n > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        T temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new T[2*n+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+n] = temp[i+order];
        order = n;
    }

    //Product
    psup = min(n+a.getOrder(),  a.getOrder());
    pinf = max(n-b.getOrder(), -b.getOrder());
    for(p=pinf; p<= psup; p++)
    {
        this->addCoef(a.getCoef(p)*b.getCoef(n-p), n);
    }


    return *this;
}

/*
   Performs the scalar-polynomial multiplication + sum h = h + r*a

   parameters:
   h: the polynomial output
   r: the scalar
   a: the polynomial input
*/
template<typename T> Ofs<T>& Ofs<T>::smult(Ofs<T> const& a, T const& c)
{
    if(order != a.order)
    {
        cout << "Error using smult: the order of variables does not match. Initial Ofs<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for(int i=0; i<2*order+1; i++) coef[i] += c*a.coef[i];
        return *this;
    }
}

/*
   Performs the scalar-polynomial multiplication h = r*a

   parameters:
   h: the polynomial output
   r: the scalar
   a: the polynomial input
*/
template<typename T> Ofs<T>& Ofs<T>::mult(Ofs<T> const& a, T const& c)
{
    if(order != a.order)
    {
        cout << "Error using smult: the order of variables does not match. Initial Ofs<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for(int i=0; i<2*order+1; i++) coef[i] = c*a.coef[i];
        return *this;
    }
}


/*
    Puts a taylor serie of TWO variables back into a fourier serie assuming x0 = exp[-it], x1 = exp[+it].
    Both series have to be properly initialized
    Note that:
    fs is of order J, of size 2J+1
    ts is of order 2*J
    ts coeffs is of size binomial(2*J+2, 2)
*/
template <typename T>  Ofs<T>& Ofs<T>::ts2fs(Ots<T> const& ts)
{
    int k, l, fspow;
    //Cleaning of fs
    this->zero();
    //Order 0
    this->coef[this->order] = (T) ts.getCoef(0,0);
    //Orders > 0
    for(k = 1; k<= ts.getOrder() ; k++)
    {
        for(l=0; l<= k; l++)
        {
            fspow = 2*l-k;
            if(fspow <= this->order && fspow >= -this->order)
            {
                //printf("k: %i, l: %i, fspower: %i, cc: %i, indix of ts: %i\n", k, l, fspow, cc, cc+l);
                this->coef[this->order+fspow] += ts.getCoef(k, l);
            }

        }
    }

    return *this;
}

/*
    Puts a fourier serie into a taylor serie of TWO variables, with x0 = exp[-it], x1 = exp[+it].
    Both series have to be properly initialized
    Note that:
    fs is of order J, of size 2J+1
    ts is of order 2*J
    ts coeffs is of size binomial(2*J+2, 2)
*/
template<typename T> Ots<T>& fs2ts(Ots<T> *ts, Ofs<T> const& fs)
{
    int j, cc;
    //Order 0
    ts->setCoef(fs.getCoef(0), 0);

    //Orders > 0
    cc = 1;
    for(j = 1; j<= fs.getOrder(); j++)
    {
        ts->setCoef(fs.getCoef(-j), cc);
        cc+=FTDA::nmon(2,j);
        ts->setCoef(fs.getCoef(j), cc-1);
    }

    return *ts;
}


// Functions
//---------------------------------------------------------------------------
template <typename T>  bool operator==(Ofs<T> const& a, Ofs<T> const& b)
{
    return a.isEqual(b);
}

/*
    Returns a+b at order max(a.order, b.order)
*/
template <typename T>  Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop+=b;
    return cop;
}

/*
    Returns a-b at order max(a.order, b.order)
*/
template <typename T>  Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop-=b;
    return cop;
}

/*
    Returns -b at order b.order
*/
template <typename T>  Ofs<T> operator - (Ofs<T> const& b)
{
    Ofs<T> ofs(b.getOrder());
    for(int i = -b.getOrder(); i<= b.getOrder() ; i++) ofs.setCoef(-b.getCoef(i),i);
    return ofs;
}

/*
    Returns a*c at order a.order
*/
template <typename T>  Ofs<T> operator * (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/*
    Returns c*a at order a.order
*/
template <typename T>  Ofs<T> operator * (T const& c, Ofs<T> const& a)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/*
    Returns a/c at order a.order
*/
template <typename T>  Ofs<T> operator / (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop/=c;
    return cop;
}

//Stream
//---------------------------------------------------------------------------
template <typename T>  std::ostream& operator << (std::ostream& stream, Ofs<T> const& ofs)
{
    stream << "Fourier serie" << endl;
    //Order
    stream << "Order : " << ofs.order << endl;
    //Coefficients
    for(int i = 0 ; i< 2*ofs.order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.order << "   " <<  setiosflags(ios::scientific) << setprecision(15) << ofs.coef[i] << endl;

    }
    return stream;
}

// Complex case
template <>  inline std::ostream& operator << (std::ostream& stream, Ofs< double complex > const& ofs)
{
    stream << "Fourier serie" << endl;
    //Order
    stream << "Order : " << ofs.order << endl;
    //Coefficients
    for(int i = 0 ; i< 2*ofs.order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.order << "   " <<  setiosflags(ios::scientific) << setprecision(15) << creal(ofs.coef[i]) << "  " << cimag(ofs.coef[i]) << endl;

    }
    return stream;
}


//Evaluate
template <typename T>  double complex Ofs<T>::evaluate(double const& t)
{
    double complex result = 0;
    //result += coef[n]*(cos(nt)+i*sin(nt));
    for(int n=-order; n<=order; n++) result+= getCoef(n)*(cos(n*t) + I*sin(n*t));
    return result;
}
