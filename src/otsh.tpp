//############################################################################
// Implementation of the Otsh template class
//############################################################################

template<typename T> int Otsh<T>::globalOrder = 20;
template<typename T> int Otsh<T>::globalVariable = 3;


//Create
//---------------------------------------------------------------------------

/*
    Default constructor (should not be used anywhere else than in this file)
*/
template<typename T> Otsh<T>::Otsh()
{
    order = 0;
    nv = 0;
    coef = new T(0);
}

/*
   Allocates memory  without any link to a coefficient array (requires the inmediate used of setCoefs afterwards).
*/
template<typename T> Otsh<T>::Otsh(int newNv, int newOrder)
{
    order = newOrder;
    nv = newNv;
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
template<typename T> void Otsh<T>::setCoefs(T *coef0)
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
template<typename T> Otsh<T>::Otsh(Otsh<T> const& b)
{
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
    this->setCoefs(coef0);
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
        this->setCoefs(coef0);
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


//Delete       TO BE DETERMINED: how properly delete with recursivity?
//              Seems to work fine like this, but memory leak?
//---------------------------------------------------------------------------
template<typename T> Otsh<T>::~Otsh<T>()
{
    /*if(coef != NULL)
    {
        delete coef;
        coef = 0;
    }*/
}

//Zeroing
//---------------------------------------------------------------------------
template<typename T> void Otsh<T>::zero()
{
    for(int i = 0 ; i<  FTDA::nmon(nv, order); i++) this->setCoef(0.0, i);
}

//Setters
//---------------------------------------------------------------------------
template<typename T> void Otsh<T>::setCoef(T value, int pos)
{
    if(pos >=0 && pos <  FTDA::nmon(nv, order)) coef[pos] = value;
    else cout << "Error in setCoef: position is out of scope\n No coefficient is set." << endl;
}

//Used to get a coef from term[0] of an Ots. As a consequence,
//Positions bigger than FTDA::nmon(nv, order) can be set, but no
//bigger than binomial(nv + ots.order, nv).
//Use with care!!
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
        coef[pos] = 2.0*((double) rand()/RAND_MAX - 0.5); //between -1 and 1
    }
}


// Conjugate (2 variables only !)
template <typename T> void Otsh<T>::conjugate()
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
template <> inline void Otsh<double>::conjugate()
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
template<typename T>  T Otsh<T>::getCoef(int i)
{
    if(i >=0 && i <  FTDA::nmon(nv, order)) return coef[i];
    else return 0.0;
}

template<typename T> T* Otsh<T>::getCoefAddress()
{
    return coef;
}


template<typename T> int Otsh<T>::getOrder()
{
    return order;
}

template<typename T> int Otsh<T>::getVariables()
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
