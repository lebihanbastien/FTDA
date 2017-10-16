//############################################################################
// Implementation of the Otsh template class
//############################################################################

/**
 * \file otsh.tpp
 * \brief Homogeneous Taylor series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


template<typename T> int Otsh<T>::globalOrder = 20;
template<typename T> int Otsh<T>::globalVariable = 3;


//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------

/*
    Default constructor (should not be used anywhere else than in this file)
*/
template<typename T> Otsh<T>::Otsh()
{
    order = 0;
    nv = 0;
    coef = new T();
}

/*
   Allocates memory  without any link to a coefficient array (requires the inmediate used of setAllCoefs afterwards).
*/
template<typename T> Otsh<T>::Otsh(int newNv, int newOrder)
{
    order = newOrder;
    nv = newNv;
    coef = new T();//[ FTDA::nmon(nv, order)]();

    //Tree
    if(nv==1)  //1-variable polynomial
    {
        term = new Otsh<T>[order+1];
        if (term == NULL)
        {
            puts("Otsh<T>::Otsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        term = new Otsh<T>[order+1];  //At this step, all sons of the current Otsh<T> have nv = 0, order = 0 and coef points @ one T.
        if (term == NULL)
        {
            puts("Otsh<T>::Otsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=order; i++)
        {
            term[i].lcopy(Otsh<T>(newNv-1, newOrder-i));
            //At this step, all sons of the current Otsh<T> have nv = nv-1, order = order-i and coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
        }
    }
}

/*
   Performs linking between the Otsh<T> and a coefficient array

   parameters:
   h: the polynomial to link;
   coef0: the array of coefficients;
*/
template<typename T> void Otsh<T>::setAllCoefs(T *coef0)
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
            term[i].setAllCoefs(coefc);
        }
    }
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/*
    Construction without any link
*/
template<typename T> Otsh<T>::Otsh(Otsh<T> const& b)
{
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Otsh<T>(b.nv, b.order));
// Goal: avoid the creation of a temporary Otsh
//-----------------------------------------------------------
    order = b.order;
    nv = b.nv;
    coef = new T();//[ FTDA::nmon(nv, order)]();
    //Tree
    if(nv==1)  //1-variable polynomial
    {
        term = new Otsh<T>[order+1];
        if (term == NULL)
        {
            puts("Otsh<T>::Otsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        term = new Otsh<T>[order+1];  //At this step, all sons of the current Otsh<T> have nv = 0, order = 0 and coef points @ one T.
        if (term == NULL)
        {
            puts("Otsh<T>::Otsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=order; i++)
        {
            term[i].lcopy(Otsh<T>(nv-1, order-i));
            //At this step, all sons of the current Otsh<T> have nv = nv-1, order = order-i and coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
        }
    }
//-----------------------------------------------------------

    //Copy in linking of the coefficients
    T *coef0 = new T[ FTDA::nmon(b.nv, b.order)]();
    for(int i = 0 ; i <  FTDA::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
    this->setAllCoefs(coef0);
}

/*
    Copy without any link (entirely new object!)
*/
template<typename T> Otsh<T>& Otsh<T>::operator = (Otsh<T> const& b)
{
    if(this != &b)
    {
        delete term;
        delete coef;
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Otsh<T>(b.nv, b.order));
// Goal: avoid the creation of a temporary Otsh
//-----------------------------------------------------------
        order = b.order;
        nv = b.nv;
        coef = new T(0);//[ FTDA::nmon(nv, order)]();
        //Tree
        if(nv==1)  //1-variable polynomial
        {
            term = new Otsh<T>[order+1];
            if (term == NULL)
            {
                puts("Otsh<T>::Otsh<T>: out of memory (2)");
                exit(1);
            }
        }
        else
        {
            int i;
            term = new Otsh<T>[order+1];  //At this step, all sons of the current Otsh<T> have nv = 0, order = 0 and coef points @ one T.
            if (term == NULL)
            {
                puts("Otsh<T>::Otsh<T>: out of memory (2)");
                exit(1);
            }

            for(i=0; i<=order; i++)
            {
                term[i].lcopy(Otsh<T>(nv-1, order-i));
                //At this step, all sons of the current Otsh<T> have nv = nv-1, order = order-i and coef points @ one T.
                //As a consequence, a proper linkage is required to ensured that coef points at a complete array of coefficients
            }
        }
//-----------------------------------------------------------
        //Copy in linking of the coefficients
        T *coef0 = new T[ FTDA::nmon(b.nv, b.order)]();
        for(int i = 0 ; i <  FTDA::nmon(b.nv, b.order) ; i++) coef0[i] = b.coef[i];
        this->setAllCoefs(coef0);
    }
    return *this; //same object if returned
}

/*
    Linked copy
*/
template<typename T> Otsh<T>& Otsh<T>::lcopy (Otsh<T> const& b)
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
template<typename T> Otsh<T>& Otsh<T>::ccopy (Otsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Otsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int i = 0 ; i< FTDA::nmon(nv, order) ; i++) coef[i] = b.coef[i];
        return *this;
    }
}


//---------------------------------------------------------------------------
//Delete
// TO BE DETERMINED: how properly delete with recursivity?
// Seems to work fine like this, but memory leak?
//---------------------------------------------------------------------------
template<typename T> Otsh<T>::~Otsh<T>()
{
    /*if(coef != NULL)
    {
        delete coef;
        coef = 0;
    }*/
}

//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
template<typename T> void Otsh<T>::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->setCoef(0.0, i);
}

template<> inline void Otsh<cdouble>::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->setCoef(0.0+0.0*I, i);
}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
template<typename T> void Otsh<T>::setCoef(T value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos] = value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

/**
 * \brief Used to get a coef from term[0] of an Ots. As a consequence,
 * Positions bigger than FTDA::nmon(nv, order) can be set, but no
 * bigger than binomial(nv + ots.order, nv).
 * Use with care!!
 */
template<typename T> void Otsh<T>::setCoefTot(T value, int pos)
{
    coef[pos] = value;
}

template<typename T> void Otsh<T>::addCoef(T value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos] += value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

template<typename T> void Otsh<T>::setRandomCoefs()
{
    for(int pos = 0; pos< FTDA::nmon(nv, order); pos++)
    {
        coef[pos] = (double) rand()/RAND_MAX; //between 0 and 1
        coef[pos] = coef[pos]/((double)(order+1)*(order+1));
    }
}

template<> inline void Otsh<cdouble>::setRandomCoefs()
{
    for(int pos = 0; pos< FTDA::nmon(nv, order); pos++)
    {
        coef[pos] =(double) rand()/RAND_MAX*I+(double) rand()/RAND_MAX; //between 0 and 1
        coef[pos] = coef[pos]/((double)(order+1)*(order+1));
    }
}

template <typename T> void Otsh<T>::conjugate()
{
    for(int pos = 0; pos< FTDA::nmon(nv, order); pos++)
    {
        coef[pos] = (sizeof(T) == sizeof(double))? (T) coef[pos]:conj(coef[pos]);
    }
}

template <> inline void Otsh<double>::conjugate()
{
    //nothing is done here.
}


// Conjugate (2 variables only !)
template <typename T> void Otsh<T>::conjugateforOFS()
{
    if(nv != 2) cout << "Error in conjugate() of Otsh: number of variables is different from 2." << endl;
    else
    {
        int nc = order+1;
        int i;
        //Temporary copy of the coeffcients at order this.order
        T temp[nc];
        for(i =0; i< nc; i++) temp[i] = coef[i];
        //Conjugate
        for(i =0; i< nc; i++)
        {
            coef[i] = conj(temp[nc-1-i]);
        }
    }
}

// Double case
template <> inline void Otsh<double>::conjugateforOFS()
{
    if(nv != 2) cout << "Error in conjugate() of Otsh: number of variables is different from 2." << endl;
    else
    {
        int nc = order+1;
        int i;
        //Temporary copy of the coeffcients at order this.order
        double temp[nc];
        for(i =0; i< nc; i++) temp[i] = coef[i];
        //Conjugate
        for(i =0; i< nc; i++)
        {
            coef[i] = temp[nc-1-i];
        }
    }
}

//Getters
//---------------------------------------------------------------------------
template<typename T> Otsh<T> Otsh<T>::getTerm()
{
    return term[0];
}

template<typename T> Otsh<T> Otsh<T>::getTerm(int i)
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
    If position is out of scope, 0.0 is returned by default
*/
template<typename T>  T Otsh<T>::getCoef(int i) const
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else return 0.0;
}

template<>  inline  cdouble Otsh<cdouble>::getCoef(int i) const
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else return 0.0+0.0*I;
}

template<typename T> T* Otsh<T>::getCA()
{
    return coef;
}


template<typename T> int Otsh<T>::getOrder()
{
    return order;
}

template<typename T> int Otsh<T>::getNV()
{
    return nv;
}



//Operators
//---------------------------------------------------------------------------
template<typename T> bool Otsh<T>::isEqual(Otsh<T> const& b) const
{
    if(order != b.order || nv != b.nv) return false;
    else
    {
        bool result = true;
        for(int i = 0 ; i<=  FTDA::nmon(nv, order); i++) result = result&&(coef[i] == b.coef[i]);
        return result;
    }
}

template<typename T> Otsh<T>& Otsh<T>::operator += (Otsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using += : the order and/or number of variables does not match. Initial Otsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] += b.coef[i];
        return *this;
    }
}

template<typename T> Otsh<T>& Otsh<T>::operator -= (Otsh<T> const& b)
{
    if(order != b.order || nv != b.nv)
    {
        cout << "Error using += : the order and/or number of variables does not match. Initial Otsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] -= b.coef[i];
        return *this;
    }
}

template<typename T> Otsh<T>& Otsh<T>::operator *= (T const& c)
{
    for(int i=0; i <  FTDA::nmon(nv, order); i++) coef[i] *= c;
    return *this;
}

template<typename T> Otsh<T>& Otsh<T>::operator /= (T const& c)
{
    for(int i=0; i <  FTDA::nmon(nv, order); i++) coef[i] /= c;
    return *this;
}

/*
   Performs the homogeneous polynomial-polynomial multiplication + sum: p = p + a*b

   parameters:
   p: the polynomial output
   a: the polynomial input
   b: the polynomial input
*/
template<typename T> Otsh<T>& Otsh<T>::sprod(Otsh<T> const& a, Otsh<T> const& b)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Otsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
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

/*
   Performs the homogeneous polynomial-complex/polynomial T multiplication + sum: p = p + m*a*b

   parameters:
   p: the polynomial output
   a: the polynomial input
   b: the polynomial input
   m: the factor input
*/
template<typename T> Otsh<T>& Otsh<T>::smprod(Otsh<T> const& a, Otsh<T> const& b, T const& m)
{
    if (nv > 2)  //3+-variate polynomial
    {
        if (a.order)
        {
            if (b.order)
            {
                Otsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.term+a.order;
                bf = b.term+b.order;
                for ( aa= a.term , pp0= this->term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->smprod(*aa , *bb, m);
                    pp->smult(*aa, m * *(bb->coef));
                }
                for ( bb= b.term , pp= pp0 ; bb< bf ; pp->smult(*bb, m * *(aa->coef)), bb++ , pp++);
                *(pp->coef)+= m * *(aa->coef) * *(bb->coef);
            }
            else this->smult(a, m * *(b.coef)); //b is scalar
        }
        else this->smult(b, m * *(a.coef));   //a is scalar
    }
    else if (nv == 2)  //2-variate homogeneous polynomial
    {
        T *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.coef +a.order;
        bf = b.coef +b.order;
        for ( aa= a.coef , pp0= this->coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.coef , pp= pp0 ; bb<= bf ; *pp+= m * *aa * *bb , bb++ , pp++ ) ;
    }
    else
    {
        *(this->coef)+= m * *(a.coef ) * *(b.coef); //1-variate homogeneous polynomial
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
template<typename T> Otsh<T>& Otsh<T>::smult(Otsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Otsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] += c*a.coef[i];
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
template<typename T> Otsh<T>& Otsh<T>::mult(Otsh<T> const& a, T const& c)
{
    if(order != a.order || nv != a.nv)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Otsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] = c*a.coef[i];
        return *this;
    }
}

template<typename T> Otsh<T>& Otsh<T>::mult(T const& c)
{
    for (int i=0 ; i <  FTDA::nmon(nv, order); i++) coef[i] = c*coef[i];
    return *this;
}


// Functions
//---------------------------------------------------------------------------
template<typename T> bool operator==(Otsh<T> const& a, Otsh<T> const& b)
{
    return a.isEqual(b);
}

/*
    Returns a+b
*/
template<typename T> Otsh<T> operator + (Otsh<T> const& a, Otsh<T> const& b)
{
    Otsh<T> cop(a);
    cop+=b;
    return cop;
}

/*
    Returns a-b
*/
template<typename T> Otsh<T> operator - (Otsh<T> const& a, Otsh<T> const& b)
{
    Otsh<T> cop(a);
    cop-=b;
    return cop;
}

// Derivative
//---------------------------------------------------------------------------
template<typename T> Otsh<T>& Otsh<T>::derh(Otsh< T > const& a, int ni)
{
    Otsh<T> *dd, *pp, *pf;
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
            this->setCoef((T) a.order*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->mult(*pp, (T) (pp - a.term));
            }
        }
    }

    return *this;
}


template<> inline Otsh<cdouble>& Otsh<cdouble>::derh(Otsh< cdouble > const& a, int ni)
{
    Otsh<cdouble> *dd, *pp, *pf;
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
            this->setCoef(a.order*a.getCoef(0)+0.0*I, 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->mult(*pp, (pp - a.term)+0.0*I);
            }
        }
    }

    return *this;
}


template<typename T> Otsh<T>& Otsh<T>::sderh(Otsh< T > const& a, int ni)
{
    Otsh<T> *dd, *pp, *pf;
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
            this->addCoef((T) a.order*a.getCoef(0), 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->smult(*pp, (T) (pp - a.term));
            }
        }
    }

    return *this;
}


template<> inline Otsh<cdouble>& Otsh<cdouble>::sderh(Otsh< cdouble > const& a, int ni)
{
    Otsh<cdouble> *dd, *pp, *pf;
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
            this->addCoef(a.order*a.getCoef(0)+0.0*I, 0);
        }
        else
        {
            for(dd=this->term, pp=a.term+1; pp<=pf; pp++,dd++)
            {
                dd->smult(*pp, (pp - a.term)+0.0*I);
            }
        }
    }

    return *this;
}


//Stream
//---------------------------------------------------------------------------
template<typename T> std::ostream& operator << (std::ostream& stream, Otsh<T> const& otsh)
{
    int i,j;
    int k[otsh.nv];
    k[0] = otsh.order;
    for(i=1; i<otsh.nv; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Degree: " << otsh.order    << endl;
    stream << "#Variables: " << otsh.nv    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< FTDA::nmon(otsh.nv, otsh.order); i++)
    {
        for(j=0; j<otsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  otsh.coef[i] << std::noshowpos << endl;

        if(i< FTDA::nmon(otsh.nv, otsh.order)-1)  FTDA::prxkt(k, otsh.nv);
    }
    return stream;
}


//Complex version
template<> inline std::ostream& operator << (std::ostream& stream, Otsh< complex double > const& otsh)
{
    int i,j;
    int k[otsh.nv];
    k[0] = otsh.order;
    for(i=1; i<otsh.nv; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Degree: " << otsh.order    << endl;
    stream << "#Variables: " << otsh.nv    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< FTDA::nmon(otsh.nv, otsh.order); i++)
    {
        for(j=0; j<otsh.nv; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  creal(otsh.coef[i]) << " " <<  cimag(otsh.coef[i]) << std::noshowpos << endl;

        if(i< FTDA::nmon(otsh.nv, otsh.order)-1)  FTDA::prxkt(k, otsh.nv);
    }
    return stream;
}



//---------------------------------------------------------------------------
// Evaluate
//---------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Otsh<T>::sevaluate(U X[], T& z)
{
    int *kv = (int*) calloc(order, sizeof(int));
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< FTDA::nmon(nv, order); i++)
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
        z += coef[i]*bux;
        if(i< FTDA::nmon(nv, order)-1)  FTDA::prxkt(kv, nv);
    }
    free(kv);
}



/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
 */
template<typename T> template<typename U> void Otsh<T>::sevaluate_conjugate(U X[], T& z)
{
    int kv[nv];
    kv[0] = order;
    for(int i=1; i<nv; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< FTDA::nmon(nv, order); i++)
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

        z += coef[i]*bux;
        if(i< FTDA::nmon(nv, order)-1)  FTDA::prxkt(kv, nv);
    }
}
