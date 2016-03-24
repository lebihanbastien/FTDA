//############################################################################
// Implementation of the Ofs template class
//############################################################################

/**
 * \file ofs.tpp
 * \brief Fourier series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofs<T>.
 */
template <typename T> Ofs<T>::Ofs()
{
    order = OFS_ORDER;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor with a given order.
 */
template <typename T>
Ofs<T>::Ofs(const int newOrder)
{
    order = newOrder;
    coef = new T[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor from a given Ofs object.
 */
template <typename T> Ofs<T>::Ofs(Ofs const& b)
{
    order = b.order;
    coef = new T[2*order+1];
    for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order]; //this->setCoef(b.getCoef(i), i);
}

//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofs<T>.
 */
template <typename T> Ofs<T>::~Ofs<T>()
{
    if(coef != NULL) delete coef;
    coef = 0;
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief  Copy from a given Ofs object (only the coefficients).
 */
template <typename T> Ofs<T>& Ofs<T>::ccopy(Ofs<T> const& b)
{
    if(order != b.order)
    {
        cout << "Erreur in Ofs<T>::ccopy: orders do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {
        for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order];//this->setCoef(b.getCoef(i), i);
        return *this;
    }
}

/**
 *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
 */
template <typename T> Ofs<T>& Ofs<T>::lcopy(Ofs<T> const& b)
{
    order = b.order;
    coef = b.coef;
    return *this;
}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the serie.
 */
template <typename T> void Ofs<T>::setCoef(T const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] = value;
    else cout << "Error in Ofs<T>::setCoef: position is out of scope. No coefficient is set." << endl;
}

/**
 *  \brief Sets a coefficient at a given position in the serie (cdouble version)
 */
template <typename T> template <typename U> void Ofs<T>::setCoef(U const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] = value+0.0*I;
    else cout << "Error in Ofs<T>::setCoef: position is out of scope. No coefficient is set." << endl;
}

/**
 *  \brief Adds a coefficient at a given position in the serie.
 */
template <typename T> void Ofs<T>::addCoef(T const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] += value;
    else cout << "Error in Ofs<T>::addCoef: position is out of scope\n" << endl;
}

/**
 *  \brief Sets a coefficient to all positions in the serie.
 */
template <typename T> void Ofs<T>::setAllCoefs(T const& value)
{
    for(int pos = -order; pos <= order; pos++) this->setCoef(value, pos);
}

/**
 *  \brief Sets random coefficients to all positions in the serie.
 */
template <typename T> void Ofs<T>::setRandomCoefs()
{
    //all coefficients between -1 and 1
    //for(int pos = -order; pos <= order; pos++)  this->setCoef(2.0*((double) rand()/RAND_MAX - 0.5), pos);
    for(int pos = -order; pos <= order; pos++)  this->setCoef((double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX), pos);
}

/**
 *  \brief Sets random cdouble coefficients to all positions in the serie.
 */
template <> inline void Ofsc::setRandomCoefs()
{
    //all coefficients between expect order zero.
    for(int pos = -order; pos <= order; pos++)  this->setCoef((double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX) + I*(double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX), pos);
    //warning: order zero is set real! (important for precision).
    double rd = (double) rand();
    cdouble fac = (cdouble) (rd/(10.0*(fabs(pow(0,7.0))+1)*RAND_MAX) + 0.0I);
    this->setCoef(fac, 0);
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the serie.
 */
template <typename T> int Ofs<T>::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the pointer address of the Ofs object
 */
template <typename T> Ofs<T>* Ofs<T>::getAddress() const
{
    return (Ofs<T>*) this;
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <typename T> T Ofs<T>::tfs_getCoef(int pos) const
{
    if(pos < 2*order+1) return coef[pos];
    else
    {
        cout << "Warning in Ofs<T>::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return (T) 0.0;
    }
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <> inline cdouble Ofsc::tfs_getCoef(int pos) const
{
    if(pos < 2*order+1)  return coef[pos];
    else
    {
        cout << "Warning in Ofs<T>::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return 0.0+0.0*I;
    }
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <typename T> T Ofs<T>::ofs_getCoef(int pos) const
{
    if(fabs(pos) <= order)  return coef[pos+order];
    else
    {
        cout << "Warning in Ofs<T>::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return (T) 0.0;
    }
}

/**
 *  \brief  Gets the coefficient at a given position. cdouble case
 */
template <> inline cdouble Ofsc::ofs_getCoef(int pos) const
{
    if(fabs(pos) <= order)  return coef[pos+order];
    else
    {
        cout << "Warning in Ofs<T>::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return  0.0+0.0*I;
    }
}

/**
 *  \brief  Computes the maximum coefficient in norm.
 */
template <typename T> void Ofs<T>::getCoefMaxNorm(double maxAbs[]) const
{
    int n = this->getOrder();
    maxAbs[0] = cabs(coef[-n+order]);//cabs(this->getCoef(-n));
    maxAbs[1] = -n;

    //Loop on all the coefficient but -n
    for(int i = -n+1; i <= n; i++)
    {
        if(cabs(coef[i+order]) > maxAbs[0])
        {
            maxAbs[0] = cabs(coef[i+order]);
            maxAbs[1] = i;
        }
    }
}

/**
 *  \brief  Computes the maximum coefficient in norm.
 */
template <typename T> double Ofs<T>::getCoefMaxNorm() const
{
    int n = this->getOrder();
    double maxAbs = cabs(coef[-n+order]);//cabs(this->getCoef(-n));

    //Loop on all the coefficient but -n
    for(int i = -n+1; i <= n; i++)
    {
        if(cabs(coef[i+order]) > maxAbs)
        {
            maxAbs =cabs(coef[i+order]);
        }
    }

    return maxAbs;
}

/**
 *  \brief  Is the Ofs object equal to zero, at order ofs_order?
 */
template <typename T> bool Ofs<T>::isnull(const int ofs_order) const
{
    for(int i = -min(ofs_order, order); i <= min(ofs_order, order); i++)
    {
        if(cabs(ofs_getCoef(i)) != 0.0) return false;
    }
    return true;
}

//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template <typename T> void Ofs<T>::zero()
{
    for(int i = -order ; i<= order; i++) coef[i+order] = 0.0;
}

/**
 *  \brief  Sets all coefficients to zero. cdouble case
 */
template <> inline void Ofsc::zero()
{
    for(int i = -order ; i<= order; i++) coef[i+order] = 0.0+0.0*I;
}


//---------------------------------------------------------------------------
// Frequency domain <--> Time domain
//---------------------------------------------------------------------------
/**
 *  \brief  From Frequency domain to time domain.
 */
template <typename T> void Ofs<T>::tfs_from_ofs(Ofs<T> const& a)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(order < a.getOrder())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    int N = 2*order+1;
    //---------------------
    //Copy the coefficients in coef
    //---------------------
    for(int i = 0; i < N; i++) coef[i] = ((Ofs<T>)a).evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  Inline from Frequency domain to time domain.
 */
template <typename T> void Ofs<T>::tfs_from_ofs_inline(Ofs<T>& temp)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(order != temp.getOrder())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    //Copy in temp
    //---------------------
    temp.ccopy(*this);

    int N = 2*order+1;
    //---------------------
    //Copy the coefficients in coef
    //---------------------
    for(int i = 0; i < N; i++) coef[i] = temp.evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  From Time domain to Frequency domain.
 */
template <typename T> void Ofs<T>::tfs_to_ofs(Ofs<T> const& a)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(order > a.getOrder())
    {
        cout << "tfs_to_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    // FFT structures
    //---------------------
    int N = 2*a.getOrder()+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_alloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(a.coef[i]), cimag(a.coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without setCoef
    //this->coef[order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without setCoef
        //this->coef[order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->coef[order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

/**
 *  \brief  Inline from Time domain to Frequency domain.
 */
template <typename T> void Ofs<T>::tfs_to_ofs_inline()
{
    //---------------------
    // FFT structures
    //---------------------
    int N = 2*order+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_alloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(this->coef[i]), cimag(this->coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without setCoef
    //this->coef[order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without setCoef
        //this->coef[order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->coef[order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

/**
 *  \brief  From Time domain to Frequency domain. Alternative version with Fortran code, for real values only
 */
template <typename T> void Ofs<T>::tfs_to_ofs_F(Ofs<T>& a)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(order > a.getOrder())
    {
        cout << "tfs_to_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    //FFT
    //---------------------
    int N = 2*order+1;
    int M = order;
    double CSF[M+1], SIF[M+1], F[N];
    for(int i = 0; i < N; i++) F[i]  = creal(a.coef[i]);
    foun_(F, &N, &M, CSF, SIF);


    //---------------------
    //Order 0
    //---------------------
    this->setCoef(CSF[0],  0);

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(0.5*(CSF[i] + I*SIF[i]), -i);
        //Positive frequencies
        this->setCoef(0.5*(CSF[i] - I*SIF[i]),  i);
    }

}

/**
 *  \brief  From Time domain to Frequency domain. Version with externalized GSL structures.
 */
template <typename T> void Ofs<T>::tfs_to_ofs(Ofs<T>& a, gsl_vector_complex *data_fft, gsl_fft_complex_wavetable *wavetable, gsl_fft_complex_workspace *workspace)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(order != a.getOrder())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    // FFT structures
    //---------------------
    int N = 2*order+1;

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(a.coef[i]), cimag(a.coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without setCoef
    //this->coef[order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without setCoef
        //this->coef[order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->coef[order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }
}


//---------------------------------------------------------------------------
//Operators (+=, -=, ...)
//---------------------------------------------------------------------------
/**
 *  \brief  An operator. Constructor from a given Ofs object (only the coefficients).
 */
template <typename T> Ofs<T>& Ofs<T>::operator = (Ofs<T> const& b)
{
    if(this != &b)
    {
        order = b.order;
        if(coef != NULL) delete coef;
        coef = new T[2*order+1];
        for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order]; //this->setCoef(b.getCoef(i), i);
    }
    return *this; //same object if returned
}

/**
 *  \brief  An operator. Constructor from a given coefficient. The returned object is of order zero.
 */
template <typename T> Ofs<T>& Ofs<T>::operator  = (T const& coef0)
{
    order = 0;
    coef = new T[2*order+1];
    coef[order] = coef0;
    return *this;
}

/**
 * \brief An operator. Compares two Ofs objects
 */
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

/**
 *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
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

/**
 *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
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

/**
 *  \brief  An operator. Multiplies all coefficients by a given \c T coefficient.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator *= (T const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] *= c;
    return *this;
}

/**
 *  \brief  An operator. Divides all coefficients by a given \c T coefficient.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator /= (T const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] /= c;
    return *this;
}


//---------------------------------------------------------------------------
//Operators (sprod, smult, ...)
//---------------------------------------------------------------------------
/**
 *  \brief Conjugates the Ofs object.
 */
template <typename T> void Ofs<T>::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -order; k <=order; k++) coef[k+order] = ofs_temp.coef[-k+order];//this->setCoef(ofs_temp.getCoef(-k), k);
}

/**
 *  \brief Conjugates the Ofsc object (the conjugate of each term is taken).
 */
template <> inline void Ofs<cdouble >::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -order; k <=order; k++) coef[k+order] = conj(ofs_temp.coef[-k+order]);//this->setCoef(conj(ofs_temp.getCoef(-k)), k);
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.

    Notes:
    1. \c this \f$ += a \times b \f$. with new order = max(order, a.order, b.order): the sum is truncated at max(a.getOrder(),b.getOrder()).
    2. WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in: this->addCoef(a.getCoef(p)*b.getCoef(n-p), n). The getCoef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.order = b.order which is the default case.
 */
template <typename T>  void Ofs<T>::ofs_sprod(Ofs<T> const& a, Ofs<T> const& b)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        //psup = n>0? order: n+order;//min(n+order,  order);
        //pinf = n>0? n-order: -order;//max(n-order, -order);
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        for(int p=pinf; p<= psup; p++) coef[n+order] += a.coef[p+order]*b.coef[n-p+order];
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$.
 */
template <typename T>  void Ofs<T>::ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, Ofs<T>& temp)
{
    //temp = a*b
    temp.ofs_prod(a,b);
    //this = temp*c = a*b*c
    this->ofs_sprod(temp, c);
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \times b \f$.

    Notes:
    1. \c this \f$ += a \times b \f$. with new order = max(order, a.order, b.order): the sum is truncated at max(a.getOrder(),b.getOrder()).
    2. WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in: this->addCoef(a.getCoef(p)*b.getCoef(n-p), n). The getCoef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.order = b.order which is the default case.
 */
template <typename T>  void Ofs<T>::ofs_smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            coef[n+order] += m*a.coef[p+order]*b.coef[n-p+order];
        }
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = a \times b \f$.
 *
 *  Notes: see smprod.
 */
template <typename T> void Ofs<T>::ofs_prod(Ofs<T> const& a, Ofs<T> const& b)
{
    int psup, pinf;
    //Product

    for(int n=-order ; n<= order; n++)
    {
        //psup = n>0? order: n+order;//min(n+order,  order);
        //pinf = n>0? n-order: -order;//max(n-order, -order);
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        coef[n+order] = 0.0; //reset

       for(int p=pinf; p<= psup; p++) coef[n+order] += a.coef[p+order]*b.coef[n-p+order];
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = a \times b \f$. cdouble case
 *
 *  Notes: see smprod.
 */
template <> inline void Ofsc::ofs_prod(Ofsc const& a, Ofsc const& b)
{
    int psup, pinf;
    //Product

    for(int n=-order ; n<= order; n++)
    {
        //psup = n>0? order: n+order;//min(n+order,  order);
        //pinf = n>0? n-order: -order;//max(n-order, -order);
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        coef[n+order] = 0.0+0.0*I; //reset

       for(int p=pinf; p<= psup; p++) coef[n+order] += a.coef[p+order]*b.coef[n-p+order];
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = m a \times b \f$.
 *
 *  Notes: see smprod.
 */
template <typename T>  void Ofs<T>::ofs_mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        coef[n+order] = 0.0; //reset
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            coef[n+order] += m*a.coef[p+order]*b.coef[n-p+order];
        }
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = m a \times b \f$. cdouble case
 *
 *  Notes: see smprod.
 */
template <>  inline  void Ofsc::ofs_mprod(Ofsc const& a, Ofsc const& b, cdouble const& m)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        coef[n+order] = 0.0+0.0*I; //reset
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            coef[n+order] += m*a.coef[p+order]*b.coef[n-p+order];
        }
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_smult(Ofs<T> const& a, T const& c)
{
    if(order != a.order)
    {
        cout << "Error using smult: the order of variables does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        //Sum
        for(int i = -order; i <= order; i++) coef[i+order] += c*a.coef[i+order];//addCoef(c*a.getCoef(i), i);
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at a certain order eff_order.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_smult(Ofs<T> const& a, T const& c, int eff_order)
{
    //Sum
    for(int i = -eff_order; i <= eff_order; i++)
    {
        addCoef(c*a.ofs_getCoef(i), i);
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_mult(Ofs<T> const& a, T const& c)
{
    if(order != a.order)
    {
        cout << "Error using smult: the order of variables does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        //Sum
        for(int i=0; i<2*order+1; i++) coef[i] = c*a.coef[i];
    }
}

/**
 *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb)
{
    if(order != a.order ||  order != b.order)
    {
        cout << "Error using fsum: the order does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        for(int i=-order; i<=order; i++)
        {
            coef[i+order] = ma*a.coef[i+order] + mb*b.coef[i+order];//setCoef(ma*a.getCoef(i)+mb*b.getCoef(i), i);
        }
    }
}

/**
 *  \brief  An operation. Performs the expansion this \f$ = a^{\alpha} \f$ when \c a is of the form \c this \f$ a = 1 + d \f$ with \f$ |d|_1 << 1 \f$.
 *
 *  Note: works well with at least a factor \f$ 1e-4 \f$ between the order zero and the order one.
 */
template<typename T> void Ofs<T>::ofs_epspow(Ofs<T> const& a, T const& alpha)
{
    //set to zero
    this->zero();
    //order 0
    this->setCoef(1.0, 0);
    //order 1
    Ofs<T> a0(a);
    a0.setCoef(0.0,0); //a0 = a without the first order
    Ofs<T> ak(a0);
    Ofs<T> akc(a0);
    int nf = a.getOrder();
    double facinv = 1;
    T coef = alpha;
    //this = this + alpha*a
    this->ofs_smult(ak, facinv*coef);
    for(int k = 2; k<= nf; k++)
    {
        //1/fac(k)
        facinv*=1.0/k;
        //alpha*(alpha-1)...*(alpha-k+1)
        coef*= alpha - k + 1;
        //akc = ak;
        akc.ccopy(ak);
        //ak = akc*a0
        ak.ofs_prod(akc, a0);
        //this += facinv*coef*ak
        this->ofs_smult(ak, facinv*coef);
    }
}

/**
 *  \brief An operation. Performs the expansion this \f$ = a^{\alpha} \f$ using inverse FFT, the power function in time domain, and finally direct FFT.
 **/
template<typename T> void Ofs<T>::ofs_pows(Ofs<T> const& a, T const& alpha)
{
    this->tfs_from_ofs(a);
    this->tfs_pows(alpha);
    this->tfs_to_ofs_inline();
}

//---------------------------------------------------------------------------
// TFS operations
//---------------------------------------------------------------------------
//-----------------
// pows
//-----------------
/**
 *  \brief  An operation. Performs the power this = a^alpha in time domain.
 */
template<typename T> void Ofs<T>::tfs_pows(Ofs<T> const& a, T const& alpha)
{
    for(int i = 0; i < 2*order+1; i++) coef[i] = cpow(a.coef[i], alpha);
}

/**
 *  \brief  An operation. Performs the power this = this^alpha in time domain.
 */
template<typename T> void Ofs<T>::tfs_pows(T const& alpha)
{
    for(int i = 0; i < 2*order+1; i++) coef[i] = cpow(coef[i], alpha);
}

//-----------------
// sprod
//-----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ in time domain.

    Notes:
    1. \c this \f$ += a \times b \f$. with new order = max(order, a.order, b.order): the sum is truncated at max(a.getOrder(),b.getOrder()).
    2. WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in: this->addCoef(a.getCoef(p)*b.getCoef(n-p), n). The getCoef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.order = b.order which is the default case.
 */
template <typename T>  void Ofs<T>::tfs_sprod(Ofs<T> const& a, Ofs<T> const& b)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += a.coef[k]*b.coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T>  void Ofs<T>::tfs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += a.coef[k]*b.coef[k]*c.coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T> template<typename U> void Ofs<T>::tfs_smprod_tu(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, U const& m)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += m*a.coef[k]*b.coef[k]*c.coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T> void Ofs<T>::tfs_smprod(Ofs<T> const& a, Ofs<T> const& b,  T const& m)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += m*a.coef[k]*b.coef[k];
}

//---------------------------------------------------------------------------
//Derivation
//---------------------------------------------------------------------------
/**
 *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$ in frequency domain.
 */
template<typename T> void Ofs<T>::dot(Ofs<T> const& a, double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-order; k<= order; k++) coef[k+order] = k*n*I*a.coef[k+order];
}

/**
 *  \brief  An operation. Applies the time derivative with pulsation \f$ \omega = n \f$ directly to this in frequency domain.
 */
template<typename T> void Ofs<T>::dot(double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-order; k<= order; k++) coef[k+order] = k*n*I*coef[k+order];
}

//---------------------------------------------------------------------------
//Change of format
//---------------------------------------------------------------------------
/**
 *  \brief  Performs the transformation from an Ots object (Object Taylor Serie).
 *
 *  ts is a Taylor serie of two variables, namely \f$ \sigma_1 = e^(-it) \f$ and \f$ \sigma_2 = e^(-it) \f$.
 *  Notes:
 *  1. Both series have to be properly initialized
 *  2. If this is of order J (size 2J+1), ts is of order 2J (of size binomial(2J+2,2)
 */
template <typename T>  void Ofs<T>::ts2fs(Ots<T> const& ts)
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
}

/**
 * \fn template<typename T> Ots<T>& fs2ts(Ots<T> & ts, Ofs<T> const& fs)
 * \brief Performs the transformation from an Ofs object (fs) to an Ots object (ts).
 *
 *  ts is a Taylor serie of two variables, namely \f$ \sigma_1 = e^(-it) \f$ and \f$ \sigma_2 = e^(-it) \f$.
 *  Notes:
 *  1. Both series have to be properly initialized
 *  2. If this is of order J (size 2J+1), ts is of order 2J (of size binomial(2J+2,2)
 */
template<typename T> void fs2ts(Ots<T> *ts, Ofs<T> const& fs)
{
    int j, cc;
    //Order 0
    ts->setCoef(fs.ofs_getCoef(0), 0);

    //Orders > 0
    cc = 1;
    for(j = 1; j<= fs.getOrder(); j++)
    {
        ts->setCoef(fs.ofs_getCoef(-j), cc);
        cc+=FTDA::nmon(2,j);
        ts->setCoef(fs.ofs_getCoef(j), cc-1);
    }
}

//---------------------------------------------------------------------------
// Functions (+, -,...)
//---------------------------------------------------------------------------
/**
 * \fn template<typename T> bool   operator == (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Compares two Ofs objects
 */
template <typename T>  bool operator == (Ofs<T> const& a, Ofs<T> const& b)
{
    return a.isEqual(b);
}

/**
 * \fn template<typename T> Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a+b
 */
template <typename T>  Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop+=b;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a-b
 */
template <typename T>  Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop-=b;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& b)
 * \brief Returns -b at order b.order
 */
template <typename T>  Ofs<T> operator - (Ofs<T> const& b)
{
    Ofs<T> ofs(b.getOrder());
    for(int i = -b.getOrder(); i<= b.getOrder() ; i++) ofs.setCoef(-b.ofs_getCoef(i),i);
    return ofs;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& b). cdouble case
 * \brief Returns -b at order b.order
 */
template <>  inline  Ofsc operator - (Ofsc const& b)
{
    Ofsc ofs(b.getOrder());
    for(int i = -b.getOrder(); i<= b.getOrder() ; i++) ofs.setCoef(0.0*I-b.ofs_getCoef(i),i);
    return ofs;
}

/**
 * \fn template<typename T> Ofs<T> operator * (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the product c*a
 */
template <typename T>  Ofs<T> operator * (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator * (T const& c, Ofs<T> const& a)
 * \brief An operator. Makes the product c*a
 */
template <typename T>  Ofs<T> operator * (T const& c, Ofs<T> const& a)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator / (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the division a/c
 */
template <typename T>  Ofs<T> operator / (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop/=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator *  (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the product a*b.
 *
 *  Was supposed to be used instead of sprod/prod. In fact, makes the code a bit unclear
 *  and adds a hidden temporary variable to the implementation.
 *  As a consequence, use sprod/prod in priority.
 */
template <typename T>  Ofs<T>  operator * (Ofs<T> const& a, Ofs<T> const& b)
{
    int n, p, tO, psup, pinf;
    tO = max(a.getOrder(), b.getOrder());

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
            temp.addCoef(a.ofs_getCoef(p)*b.ofs_getCoef(n-p), n);
        }
    }

    return temp;
}


//---------------------------------------------------------------------------
// Functions (change of format)
//---------------------------------------------------------------------------
/**
 * \fn void inline doubleToComplex(Ofsd& xr, Ofsc& xc)
 * \brief Copy from Ofsd to Ofsc.
 */
void inline doubleToComplex(Ofsd const& xr, Ofsc& xc)
{
    int nf = xr.getOrder(); //order of the expansion
    if(nf != xc.getOrder()) //checking that the orders match
    {
        cout << "doubleToComplex: orders do not match." << endl;
    }
    else for(int l = -nf; l<=nf; l++) xc.setCoef((cdouble) (xr.ofs_getCoef(l)+I*0.0), l); //copy from double to cdouble
}


//---------------------------------------------------------------------------
// Functions (Real and Imaginary part)
//---------------------------------------------------------------------------
/**
 * \fn void inline realPart(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xr = real(x).
 */
void inline realPart(Ofsc const& x, Ofsc& xr)
{
    int nf = x.getOrder();
    //xc = conj(xc)
    Ofsc xc(x);
    xc.conjugate();
    //Storing the real part
    for(int l = -nf; l<=nf; l++) xr.setCoef(0.5*(x.ofs_getCoef(l) + xc.ofs_getCoef(l)), l);
}

/**
 * \fn void inline realPart(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xi = imag(x).
 */
void inline imagPart(Ofsc const& x, Ofsc& xi)
{
    int nf = x.getOrder();
    //xc = conj(xc)
    Ofsc xc(x);
    xc.conjugate();
    //Storing the imag part
    for(int l = -nf; l<=nf; l++) xi.setCoef(-0.5*I*(x.ofs_getCoef(l) - xc.ofs_getCoef(l)), l);
}


//---------------------------------------------------------------------------
// Functions (evaluate)
//---------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::fevaluate(double cR[], double sR[], int eff_order)
{
    cdouble result = 0+0.0*I;
    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];
    //Return result
    return result;
}

/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluate(double const& theta, int eff_order)
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    //Initialisation of the cosinus/sinus arrays
    initcRsR(theta, cR, sR, eff_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);


    return result;
    //Obsolete: for(int k=-order; k<=order; k++) result+= getCoef(k)*(cos((double)k*theta) + I*sin((double)k*theta));
}

/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluate(double const& theta)
{
    cdouble result = 0+0.0*I;
    double cR[order];
    double sR[order];

    //Initialisation of the cosinus/sinus arrays
    initcRsR(theta, cR, sR, order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);

    return result;
    //Obsolete: for(int k=-order; k<=order; k++) result+= getCoef(k)*(cos((double)k*theta) + I*sin((double)k*theta));
}

/**
 *  \brief  Evaluates the derivatives of the Ofs object at angle theta, with pulsation n (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluatedot(double const& theta, double const& n)
{
    cdouble result = 0+0.0*I;
    double cR[order];
    double sR[order];

    //Initialisation of the cosinus/sinus arrays
    initcRsR(theta, cR, sR, order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=order; k>=1; k--) result+= -n*k*sR[k-1]*(coef[k+order] + coef[-k+order]) + I*n*k*cR[k-1]*(coef[k+order] - coef[-k+order]);

    return result;
    //Obsolete: for(int k=-order; k<=order; k++) result+= getCoef(k)*(cos((double)k*t) + I*sin((double)k*t));
}

/**
 *  \brief  Evaluates the derivatives of the Ofs object at angle theta, with pulsation n (theta = nt), at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluatedot(double const& theta, double const& n, int eff_order)
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    //Initialisation of the cosinus/sinus arrays
    initcRsR(theta, cR, sR, eff_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= -n*k*sR[k-1]*(coef[k+order] + coef[-k+order]) + I*n*k*cR[k-1]*(coef[k+order] - coef[-k+order]);

    return result;
    //Obsolete: for(int k=-order; k<=order; k++) result+= getCoef(k)*(cos((double)k*t) + I*sin((double)k*t));
}

/**
 *  \brief  Expected error on the product: \c this \f$ += a \times b \f$ at times t. Works only when a.order = b.order which is the default case.
 */
template <typename T>  cdouble sprod_expectedError(Ofs<T> const& a, Ofs<T> const& b, double const& t)
{
    cdouble result = 0.0+0.0*I;
    for(int n= a.getOrder()+1 ; n<= 2*a.getOrder(); n++)
    {
        for(int p=n-a.getOrder(); p<= a.getOrder(); p++)
        {
            //indix
            result += a.ofs_getCoef(p)*b.ofs_getCoef(n-p)*(cos(n*t) + I*sin(n*t));
            result += a.ofs_getCoef(-p)*b.ofs_getCoef(-n+p)*(cos(-n*t) + I*sin(-n*t));
        }
    }
    return result;
}

/**
 *  \brief  Expected error on the product: \c this \f$ += c*a \times b \f$ at times t. Works only when a.order = b.order which is the default case.
 */
template <typename T>  cdouble smprod_expectedError(Ofs<T> const& a, Ofs<T> const& b, T const& c, double const& t)
{
    cdouble result = 0.0+0.0*I;
    for(int n= a.getOrder()+1 ; n<= 2*a.getOrder(); n++)
    {
        for(int p=n-a.getOrder(); p<= a.getOrder(); p++)
        {
            //indix
            result += c*a.ofs_getCoef(p)*b.ofs_getCoef(n-p)*(cos(n*t) + I*sin(n*t));
            result += c*a.ofs_getCoef(-p)*b.ofs_getCoef(-n+p)*(cos(-n*t) + I*sin(-n*t));
        }
    }
    return result;
}

//Stream
//---------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
template <typename T>  std::ostream& operator << (std::ostream& stream, Ofs<T> const& ofs)
{
    //stream << "Fourier serie" << endl;
    //Order
    //stream << "Order : " << ofs.order << endl;
    //Coefficients
    for(int i = 0 ; i< 2*ofs.order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.order << "   " <<  setiosflags(ios::scientific) << setprecision(15) << ofs.coef[i] << endl;

    }
    return stream;
}

/**
 *  \brief  A stream operator in the \c double \c complex case.
 */
template <>  inline std::ostream& operator << (std::ostream& stream, Ofs< cdouble > const& ofs)
{
    //Coefficients
    for(int i = 0 ; i< 2*ofs.order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.order << "   " <<  setiosflags(ios::scientific) << setprecision(16) << creal(ofs.coef[i]) << "  " << cimag(ofs.coef[i]) << endl;

    }
    return stream;
}

/**
 *  \brief  Evaluates an upper bound for the \f$ L_1 \f$ norm of the current Ofs object.
 */
template<typename T> double Ofs<T>::l1norm()
{
//    double theta, temp;
//    double l1n = cabs(evaluate(0.0));
//    int N = 1000;
//    //Loop on a N+1 point grid from theta = 0 to 2pi
//    for(int idt = 1; idt <= N ; idt++)
//    {
//        theta = ((double)idt)/N*2*M_PI;
//        temp = cabs(evaluate(theta));
//        if(temp > l1n) l1n = temp;
//    }
//    return l1n;
    double l1n = 0.0;
    for(int i = 0 ; i< 2*order + 1; i++) l1n += cabs(coef[i]);
    return l1n;
}

/**
 *  \brief  Number of small divisors under a certain value
*/
template<typename T> int Ofs<T>::nsd(int odmax, double sdmax)
{
    int res = 0;
    for(int i = -odmax; i <= odmax; i++) if(cabs(this->ofs_getCoef(i)) < sdmax)
        {
            //cout << i << " sd = " << creal(this->ofs_getCoef(i)) << " " << cimag(this->ofs_getCoef(i)) << endl;
            res++;
        }
    return res;
}

/**
 *  \brief  A stream operator. Print only the autonomous term (term of order zero).
 */
template<typename T> void Ofs<T>::fprint_0(ofstream& stream)
{
    stream << setw(3) << setiosflags(ios::right) << std::showpos << 0 << "   " <<  setiosflags(ios::scientific) << setprecision(15) << creal(this->ofs_getCoef(0)) << "  " << cimag(this->ofs_getCoef(0)) << endl;
}


//-------------------------------------------------------------------------------------------------------
// Reading an OFS from a text file
//-------------------------------------------------------------------------------------------------------
/**
 * \fn void inline readOFS_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 */
void inline readOFS_txt(Ofsc& xFFT, string filename)
{
    //Init
    ifstream readStream;
    double ct, cr, ci;
    int fftN = xFFT.getOrder();

    //Reading
    readStream.open((filename+".txt").c_str());
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.setCoef(cr+I*ci, i);
    }
    readStream.close();
}

//-------------------------------------------------------------------------------------------------------
// Array cos/sin computation
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Initializes the arrays cR[] and sR[], containing the numbers cos(t), cos(nt), ... cos(ofs_order*n*t) and sin(t), sin(nt), ... sin(ofs_order*n*t), respectively.
 **/
void inline initcRsR(double t, double cR[], double sR[], int ofs_order)
{
    cR[0] = cos(t);
    sR[0] = sin(t);
    for(int i = 1; i< ofs_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }
}

//--------------------------------------------------------------------------
// Single storage
//--------------------------------------------------------------------------
/**
 *  \brief Single storage of a QBTBP Ofs object in txt files + its cosinus/sinus version.
 */
inline void ofs_sst(Ofsc &xOFS, string filename, int flag, string suffix)
{
    int nf = xOFS.getOrder();
    ofstream curentStream;

    //Storage in txt file
    curentStream.open((filename+suffix+".txt").c_str());
    curentStream <<  xOFS << endl;
    curentStream.close();

    curentStream.open((filename+"c"+suffix+".txt").c_str());
    if(flag) //even case
    {
        //Cosinus expansion version
        curentStream << 0 << " " << creal(xOFS.ofs_getCoef(0)) << endl;
        for(int l = 1; l<=nf; l++) curentStream << l << " " << creal(xOFS.ofs_getCoef(-l) + xOFS.ofs_getCoef(l))  <<  endl;

    }
    else //odd case
    {
        //Sinus expansion version
        for(int l = 0; l<=nf; l++)
            curentStream << l << " " <<  cimag(xOFS.ofs_getCoef(-l) - xOFS.ofs_getCoef(l)) << endl;
        curentStream.close();
    }
    curentStream.close();
}
