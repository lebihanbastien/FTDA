#include "gslc.h"

/**
 * \file gslc.cpp
 * \brief Additional operations on GSL objects.
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//---------------------------------------------------------------------------------------
// Matrix and vectors
//---------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);
        }
}

/**
 * \brief Transform a vector into a complex matrix with a given shift in the initial vector
 **/
void gslc_vectorToComplexMatrix(gsl_matrix_complex *m, const double y[], int rows, int columns, int shift)
{
    int i,j,k;
    gsl_complex zc;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            GSL_SET_COMPLEX(&zc, y[shift+k-1], 0.0);
            gsl_matrix_complex_set(m, i-1, j-1, zc);
        }
}

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            y[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
        }
}

/**
 * \brief Transform a complex matrix into a vector with a given shift in the final vector
 **/
void gslc_complexMatrixToVector(double y[], const gsl_matrix_complex  *m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            y[shift+k-1] = GSL_REAL(gsl_matrix_complex_get(m, i-1, j-1));
        }
}

//---------------------------------------------------------------------------------------
// Complex numbers
//---------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from 2 doubles
 **/
gsl_complex gslc_complex(double real, double imag)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, real, imag);
    return one_c;
}

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, creal(x), cimag(x));
    return one_c;
}

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c)
{
    return GSL_REAL(c) + I*GSL_IMAG(c);
}

//---------------------------------------------------------------------------------------
// Printing a real matrix
//---------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix *M)
{
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.15e ", gsl_matrix_get(M, i, j));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix *M)
{
    int im = M->size1;
    int jm = M->size2;
    double c;
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++){
                c = gsl_matrix_get(M, i, j);
                if(fabs(c) > 1e-5) printf("%+1.3e ", c);
                else printf("%+1.3e ", 0.0);
         }
        printf("\n");
    }
}

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix *M, char* fileName)
{

    FILE *f;
    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", gsl_matrix_get(M, i, j));
        fprintf(f, "\n");
    }

    fclose(f);
}

//---------------------------------------------------------------------------------------
// Printing a complex matrix
//---------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            /*if( GSL_IMAG(gsl_matrix_complex_get(M, i, j)) != 0.0)*/ printf("%+1.0e%+1.0ei ", GSL_REAL(gsl_matrix_complex_get(M, i, j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
            //else printf("%+1.0e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));

        }
        printf("\n");
    }
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex *M)
{
    int im = M->size1;
    int jm = M->size2;
    double c;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
            {
                c = GSL_REAL(gsl_matrix_complex_get(M, i, j));
                if(fabs(c) > 1e-5) printf("%+1.3e ", c);
                else printf("%+1.3e ", 0.0);
            }
        printf("\n");
    }
}

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;

    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e %+5.16e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;
    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    fprintf(f,"real part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fprintf(f,"imag part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }
    fclose(f);
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;

    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+2.0e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}

//---------------------------------------------------------------------------------------
// Printing a complex vector
//---------------------------------------------------------------------------------------
/**
 * \brief Print a complex vector into a txt file.
 **/
void gslc_vector_complex_fprintf(const gsl_vector_complex *V, char* fileName)
{
    FILE *f;

    f = fopen(fileName, "w");
    int im = V->size;

    for(int i = 0; i < im ; i++) fprintf(f, "%+5.15e  %+5.15e i\n", GSL_REAL(gsl_vector_complex_get(V, i)), GSL_IMAG(gsl_vector_complex_get(V, i)));
    fclose(f);
}
/**
 * \brief Print a complex vector.
 **/
void gslc_vector_complex_printf(const gsl_vector_complex *V)
{
    int im = V->size;

    for(int i = 0; i < im ; i++) printf("%+5.15e  %+5.15e i\n", GSL_REAL(gsl_vector_complex_get(V, i)), GSL_IMAG(gsl_vector_complex_get(V, i)));
}

//---------------------------------------------------------------------------------------
// Printing a vector
//---------------------------------------------------------------------------------------
/**
 * \brief Print a real vector.
 **/
void gslc_vector_printf(const gsl_vector *V)
{
    int im = V->size;
    for(int i = 0; i < im ; i++) printf("%+5.15e \n", gsl_vector_get(V, i));
}

//---------------------------------------------------------------------------------------
// Printing an eigensystem
//---------------------------------------------------------------------------------------
/**
 * \brief Print an eigensystem with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_printf(gsl_vector_complex *eval, gsl_matrix_complex *evec, int M)
{
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < M; i++)
    {
        eval_i = gsl_vector_complex_get (eval, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            printf("%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
}

/**
 * \brief Print an eigensystem with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_printf(gsl_matrix_complex *eval, gsl_matrix_complex *evec, int M)
{
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < M; i++)
    {
        eval_i = gsl_matrix_complex_get (eval, i, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            printf("%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
}

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_fprintf(gsl_vector_complex *eval, gsl_matrix_complex *evec, int M, char* fileName)
{
    FILE *f;
    f = fopen(fileName, "w");
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < M; i++)
    {
        eval_i = gsl_vector_complex_get (eval, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        fprintf (f,"eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        fprintf (f,"||eigenvalue|| = %+5.15e\n",
                sqrt(GSL_REAL(eval_i)*GSL_REAL(eval_i) + GSL_IMAG(eval_i)*GSL_IMAG(eval_i)));
        fprintf (f,"eigenvector = \n");
        for (int j = 0; j < M; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            fprintf(f,"%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }

     fclose(f);
}

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_fprintf(gsl_matrix_complex *eval, gsl_matrix_complex *evec, int M, char* fileName)
{
    FILE *f;
    f = fopen(fileName, "w");
    gsl_vector_complex_view evec_i;
    gsl_complex eval_i;
    for (int i = 0; i < M; i++)
    {
        eval_i = gsl_matrix_complex_get (eval, i, i);
        evec_i = gsl_matrix_complex_column (evec, i);
        fprintf (f,"eigenvalue = %+5.15e + %+5.15ei\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        fprintf (f,"eigenvector = \n");
        for (int j = 0; j < 6; ++j)
        {
            gsl_complex z =
                gsl_vector_complex_get(&evec_i.vector, j);
            fprintf(f,"%+5.15e + %+5.15ei\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
     fclose(f);
}

//---------------------------------------------------------------------------------------
// Reading GSL objects
//--------------------------------------------------------------------------------------
/**
 * \brief Read a complex vector from a txt file, obtained with the routine gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName).
 **/
void glsc_matrix_complex_read(gsl_matrix_complex *M, string filename)
{
    int n1 = M->size1;
    int n2 = M->size2;
        //Init
    ifstream readStream;
    double cr, ci;

    //Reading
    cout << "glsc_matrix_complex_read. cr = " << filename+".txt" << endl;
    readStream.open((filename+".txt").c_str());
    for(int i = 0; i < n1 ; i++)
    {
        for(int j = 0; j< n2 ; j++)
        {
            readStream >> cr;  //real part
            readStream >> ci;  //imag part

            gsl_matrix_complex_set(M, i, j, gslc_complex(cr, ci));
        }
    }
    readStream.close();
}


//---------------------------------------------------------------------------------------
// Views of a matrix
//--------------------------------------------------------------------------------------
/**
 * \brief Set the kth column of the complex matrix M in the complex vector V using the GSL "views" structures.
 **/
void gslc_matrix_complex_column(gsl_vector_complex *V, gsl_matrix_complex *M, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_column(M, k);     //Select the kth column of M
    gsl_vector_complex_memcpy(V, &vecview.vector); //Copy in V
}

/**
 * \brief Set the kth row of the complex matrix M in the complex vector V using the GSL "views" structures.
 **/
void gslc_matrix_complex_row(gsl_vector_complex *V, gsl_matrix_complex *M, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_row(M, k);        //Select the kth column of M
    gsl_vector_complex_memcpy(V, &vecview.vector); //Copy in V
}

/**
 * \brief Set the complex vector V in the kth column of the complex matrix M using the GSL "views" structures.
 **/
void gslc_matrix_complex_column_V(gsl_matrix_complex *M, gsl_vector_complex *V, int k)
{
    gsl_vector_complex_view vecview;               //Vector views
    vecview = gsl_matrix_complex_column(M, k);     //Select the kth column of M
    gsl_vector_complex_memcpy(&vecview.vector, V); //Copy of V in the kth column of M
}


//---------------------------------------------------------------------------------------
// Misc manipulations
//---------------------------------------------------------------------------------------
/**
 * \brief Gives the infinity norm of a complex matrix M
 **/
double gslc_matrix_complex_infinity_norm(gsl_matrix_complex *M)
{
    int irow = M->size1;
    int icol = M->size2;
    //Mr = |M|
    gsl_matrix *Mr = gsl_matrix_calloc(irow, icol);
    for(int i = 0; i < irow; i++)
    for(int j = 0; j < icol; j++) gsl_matrix_set(Mr, i, j, gsl_complex_abs(gsl_matrix_complex_get(M, i, j)));
    //maxr = max (Mr)
    double maxr = gsl_matrix_max(Mr);
    //Memory release
    gsl_matrix_free(Mr);
    //Return maxr
    return maxr;
}

/**
 *  \brief Normalization: xm = xm/norm(xm)
 **/
void gslc_vector_complex_normalize(gsl_vector_complex *xm)
{
    gsl_complex VN1invC;
    double VN1 = gsl_blas_dznrm2 (xm);
    GSL_SET_COMPLEX(&VN1invC, 1/VN1, 0.0);
    gsl_vector_complex_scale (xm, VN1invC);
}

/**
 * \brief Isolate the real part of a complex matrix M: Mr = real(M)
 **/
void gslc_matrix_complex_real(gsl_matrix *Mr, gsl_matrix_complex *M)
{
    int irow = M->size1;
    int icol = M->size2;
    for(int i=0; i<irow; i++) for(int j=0; j<icol; j++) gsl_matrix_set(Mr, i, j,  GSL_REAL(gsl_matrix_complex_get(M, i, j)));
}

/**
 * \brief Copy the conjugate of a complex vector xm into xc.
 **/
void gslc_vector_complex_conjugate_memcpy(gsl_vector_complex *xc, gsl_vector_complex *xm)
{
    for(int i=0; i< (int)xm->size; i++) gsl_vector_complex_set(xc, i, gsl_complex_conjugate(gsl_vector_complex_get(xm, i)));
}


//---------------------------------------------------------------------------------------
// Specific routines for Monodromy and STM matrices manipulations
//---------------------------------------------------------------------------------------
/**
 * \brief Delete one raw and one column of a given square GSL matrix
 **/
gsl_matrix_complex * gslc_matrix_complex_deleteRC(gsl_matrix_complex *M, int k)
{
    if(k >= (int)M->size1)
    {
        cout << "Error in gslc_matrix_complex_deleteRC: the provided indix is out of scope." << endl;
        return NULL;
    }
    else if((int)M->size1 != (int)M->size2)
    {
        cout << "Error in gslc_matrix_complex_deleteRC: the matrix is not square" << endl;
        return NULL;
    }
    else
    {
        int im = (int)M->size1;
        gsl_matrix_complex *result = gsl_matrix_complex_calloc(im-1, im-1);

        if(k == 0)
        {
            gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (M , 1 , 1 , im-1 , im-1);
            gsl_matrix_complex_memcpy(result, &M11.matrix);
        }
        else if(k==im-1)
        {
            gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (M , 0 , 0 , im-1 , im-1);
            gsl_matrix_complex_memcpy(result, &M11.matrix);

        }
        else
        {
                gsl_matrix_complex_view M11 = gsl_matrix_complex_submatrix (M , 0   , 0   , k   , k  );
                gsl_matrix_complex_view M12 = gsl_matrix_complex_submatrix (M , 0   , k+1 , k   , im-1-k);
                gsl_matrix_complex_view M21 = gsl_matrix_complex_submatrix (M , k+1 , 0   , im-1-k , k  );
                gsl_matrix_complex_view M22 = gsl_matrix_complex_submatrix (M , k+1 , k+1 , im-1-k , im-1-k);


                gsl_matrix_complex_view R11 = gsl_matrix_complex_submatrix (result , 0 , 0 , k    , k   );
                gsl_matrix_complex_view R12 = gsl_matrix_complex_submatrix (result , 0 , k , k    , im-1-k);
                gsl_matrix_complex_view R21 = gsl_matrix_complex_submatrix (result , k , 0 , im-1-k , k   );
                gsl_matrix_complex_view R22 = gsl_matrix_complex_submatrix (result , k , k , im-1-k , im-1-k);

                gsl_matrix_complex_memcpy(&R11.matrix, &M11.matrix);
                gsl_matrix_complex_memcpy(&R12.matrix, &M12.matrix);
                gsl_matrix_complex_memcpy(&R21.matrix, &M21.matrix);
                gsl_matrix_complex_memcpy(&R22.matrix, &M22.matrix);

        }

    return result;
    }

}

/**
 * \brief Inverse transformation of a vector during Wielandt deflation: w = 1/vw*(w + vx/(vw-vx)*(z.w)*x)
 **/
void gslc_wielandt_inv_trans(gsl_vector_complex *w, gsl_complex vw, gsl_vector_complex const *x, gsl_complex vx, gsl_vector_complex *z)
{
    gsl_complex lm;
    gsl_complex lm2;
    gsl_vector_complex *xc = gsl_vector_complex_calloc(x->size);

    //lm2 = vx/(vw-vx)
    lm2 = gsl_complex_sub(vw, vx);
    lm2 = gsl_complex_div(vx, lm2);
    //lm = vx/(vw-vx)*(z^T*w) = vx/(vw-vx)*(z.w)
    gsl_blas_zdotu(z , w , &lm);  //zdotu or zdotc??
    lm = gsl_complex_mul(lm, lm2);
    //xc = lm*x
    gsl_vector_complex_memcpy(xc, x);
    gsl_vector_complex_scale(xc, lm);
    //w = w - xc;
    gsl_vector_complex_add(w, xc);
    //w = 1/vw*w
    gsl_vector_complex_scale(w, gsl_complex_inverse(vw));
    //Normalisation
    gsl_vector_complex_scale(w, gslc_complex(1/gsl_blas_dznrm2 (w),0.0));
    //memory release
    gsl_vector_complex_free(xc);
}

/**
 *  \brief Inverse a symplectic complex matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_complex_symplectic_inverse(const gsl_matrix_complex *S0, gsl_matrix_complex *Sinv)
{
    int n = S0->size1/2;

    gsl_complex minus_one_c = gslc_complex(-1.0, 0.0);

    gsl_matrix_complex *S = gsl_matrix_complex_calloc (2*n, 2*n);
    gsl_matrix_complex_memcpy(S, S0);

    //--------------------------------------------------------------------------------------------------
    //Views of S
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_complex_view S11 = gsl_matrix_complex_submatrix (S , 0 , 0 , n , n );   //S11
    gsl_matrix_complex_view S12 = gsl_matrix_complex_submatrix (S , 0 , n , n , n );   //S12
    gsl_matrix_complex_view S21 = gsl_matrix_complex_submatrix (S , n , 0 , n , n );   //S21
    gsl_matrix_complex_view S22 = gsl_matrix_complex_submatrix (S , n , n , n , n );   //S22

    //--------------------------------------------------------------------------------------------------
    //Views of Sinv
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_complex_view Sinv11 = gsl_matrix_complex_submatrix (Sinv , 0 , 0 , n , n ); //Sinv11
    gsl_matrix_complex_view Sinv12 = gsl_matrix_complex_submatrix (Sinv , 0 , n , n , n ); //Sinv12
    gsl_matrix_complex_view Sinv21 = gsl_matrix_complex_submatrix (Sinv , n , 0 , n , n ); //Sinv21
    gsl_matrix_complex_view Sinv22 = gsl_matrix_complex_submatrix (Sinv , n , n , n , n ); //Sinv22


    //--------------------------------------------------------------------------------------------------
    // Storage
    //--------------------------------------------------------------------------------------------------
    //Sij_m = -S12
    gsl_matrix_complex *Sij_m = gsl_matrix_complex_calloc (n, n);
    gsl_matrix_complex_memcpy(Sij_m, &S12.matrix);
    gsl_matrix_complex_scale(Sij_m, minus_one_c);
    //Sinv12 = -S12^T = Sij_m^T
    gsl_matrix_complex_transpose_memcpy(&Sinv12.matrix, Sij_m);
    //Sij_m = -S21
    gsl_matrix_complex_memcpy(Sij_m, &S21.matrix);
    gsl_matrix_complex_scale(Sij_m, minus_one_c);
    //Sinv12 = -S21^T = Sij_m^T
    gsl_matrix_complex_transpose_memcpy(&Sinv21.matrix, Sij_m);
    //Sinv11 = S22^T
    gsl_matrix_complex_transpose_memcpy(&Sinv11.matrix, &S22.matrix);
    //Sinv22 = S11^T
    gsl_matrix_complex_transpose_memcpy(&Sinv22.matrix, &S11.matrix);

    //Memory release
    gsl_matrix_complex_free(Sij_m);
}

/**
 *  \brief Inverse a symplectic real matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_symplectic_inverse(const gsl_matrix *S0, gsl_matrix *Sinv)
{
    int n = S0->size1/2;
    gsl_matrix *S = gsl_matrix_calloc (2*n, 2*n);
    gsl_matrix_memcpy(S, S0);

    //--------------------------------------------------------------------------------------------------
    //Views of S
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_view S11 = gsl_matrix_submatrix (S , 0 , 0 , n , n );   //S11
    gsl_matrix_view S12 = gsl_matrix_submatrix (S , 0 , n , n , n );   //S12
    gsl_matrix_view S21 = gsl_matrix_submatrix (S , n , 0 , n , n );   //S21
    gsl_matrix_view S22 = gsl_matrix_submatrix (S , n , n , n , n );   //S22

    //--------------------------------------------------------------------------------------------------
    //Views of Sinv
    //--------------------------------------------------------------------------------------------------
    gsl_matrix_view Sinv11 = gsl_matrix_submatrix (Sinv , 0 , 0 , n , n ); //Sinv11
    gsl_matrix_view Sinv12 = gsl_matrix_submatrix (Sinv , 0 , n , n , n ); //Sinv12
    gsl_matrix_view Sinv21 = gsl_matrix_submatrix (Sinv , n , 0 , n , n ); //Sinv21
    gsl_matrix_view Sinv22 = gsl_matrix_submatrix (Sinv , n , n , n , n ); //Sinv22


    //--------------------------------------------------------------------------------------------------
    // Storage
    //--------------------------------------------------------------------------------------------------
    //Sij_m = -S12
    gsl_matrix *Sij_m = gsl_matrix_calloc (n, n);
    gsl_matrix_memcpy(Sij_m, &S12.matrix);
    gsl_matrix_scale(Sij_m, -1.0);
    //Sinv12 = -S12^T = Sij_m^T
    gsl_matrix_transpose_memcpy(&Sinv12.matrix, Sij_m);
    //Sij_m = -S21
    gsl_matrix_memcpy(Sij_m, &S21.matrix);
    gsl_matrix_scale(Sij_m, -1.0);
    //Sinv12 = -S21^T = Sij_m^T
    gsl_matrix_transpose_memcpy(&Sinv21.matrix, Sij_m);
    //for(int i =0; i<n; i++) for(int j=0; j<n; j++) gsl_matrix_set(&Sinv21.matrix, i, j, gsl_complex_conjugate(gsl_matrix_get(&Sinv21.matrix, i, j)));
    //Sinv11 = S22^T
    gsl_matrix_transpose_memcpy(&Sinv11.matrix, &S22.matrix);
    //for(int i =0; i<n; i++) for(int j=0; j<n; j++) gsl_matrix_set(&Sinv11.matrix, i, j, gsl_complex_conjugate(gsl_matrix_get(&Sinv11.matrix, i, j)));
    //Sinv22 = S11^T
    gsl_matrix_transpose_memcpy(&Sinv22.matrix, &S11.matrix);
    //for(int i =0; i<n; i++) for(int j=0; j<n; j++) gsl_matrix_set(&Sinv22.matrix, i, j, gsl_complex_conjugate(gsl_matrix_get(&Sinv22.matrix, i, j)));


    //Memory release
    gsl_matrix_free(Sij_m);
}

/**
 *  \brief Transpose + conjugate the complex matrix S intro SH: SH = S^H != S^T
 **/
void gslc_matrix_complex_H_memcpy(gsl_matrix_complex *SH, const gsl_matrix_complex *S)
{
    int n = S->size1;
    //gsl_matrix_complex_transpose_memcpy(SH, S);
    for(int i =0; i<n; i++) for(int j=0; j<n; j++) gsl_matrix_complex_set(SH, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(S, j, i)));
}

/**
 *  \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  real(J) = | 0  In |   and imag(J) = 0
 *            |-In 0  |
 **/
void glsc_matrix_complex_set_J(gsl_matrix_complex *J)
{
    if(J->size1%2!=0)
    {
        cout << "glsc_matrix_complex_set_J. J is not of the form 2n*2n." << endl;
        return;
    }
    else
    {
        gsl_complex one_c  = gslc_complex(1.0, 0.0);
        gsl_complex minus_one_c = gslc_complex(-1.0, 0.0);
        int n = J->size1/2;

        gsl_matrix_complex_set_zero(J);
        for(int j = 0; j< n; j++)
        {
            gsl_matrix_complex_set(J, j, j+n, one_c);
            gsl_matrix_complex_set(J, j+n, j, minus_one_c);
        }
    }

}

/**
 * \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  J = | 0  In |
 *      |-In 0  |
 **/
void glsc_matrix_set_J(gsl_matrix *J)
{
    if(J->size1%2!=0)
    {
        cout << "glsc_matrix_complex_set_J. J is not of the form 2n*2n." << endl;
        return;
    }
    else
    {
        int n = J->size1/2;

        gsl_matrix_set_zero(J);
        for(int j = 0; j< n; j++)
        {
            gsl_matrix_set(J, j, j+n,  1.0);
            gsl_matrix_set(J, j+n, j, -1.0);
        }
    }

}

/**
 * \brief Use the symmetry of the QBCP to compute the stable (resp. unstable) from the unstable (resp. stabl) eigenvector
 **/
void gslc_vector_complex_symcpy(gsl_vector_complex *VEP2, gsl_vector_complex *VEP1)
{
    //Opposite for i = 0, 2, 4
    gsl_vector_complex_set(VEP2, 0, gsl_complex_negative(gsl_vector_complex_get(VEP1, 0)));
    gsl_vector_complex_set(VEP2, 2, gsl_complex_negative(gsl_vector_complex_get(VEP1, 2)));
    gsl_vector_complex_set(VEP2, 4, gsl_complex_negative(gsl_vector_complex_get(VEP1, 4)));
    //Equal for i = 1, 3, 5
    gsl_vector_complex_set(VEP2, 1, gsl_vector_complex_get(VEP1, 1));
    gsl_vector_complex_set(VEP2, 3, gsl_vector_complex_get(VEP1, 3));
    gsl_vector_complex_set(VEP2, 5, gsl_vector_complex_get(VEP1, 5));
}

/**
 * \brief Matrix-vector product when the matrix is given as a product of matrices: ym = DAT[1]...DAT[M] * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_vector_product(gsl_matrix_complex **DAT, const gsl_vector_complex *xm, gsl_vector_complex *ym, int M)
{
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_vector_complex *ym1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex *ym2 = gsl_vector_complex_calloc(6);

    //ym1 = xm
    gsl_vector_complex_memcpy(ym1, xm);
    for(int i =1; i<=M;i++)
    {
        //ym2 = DAT[i]*ym1
        gsl_blas_zgemv (CblasNoTrans, one_c , DAT[i], ym1 , zero_c , ym2);
        //ym1 = ym2
        gsl_vector_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_vector_complex_memcpy(ym, ym2);

    gsl_vector_complex_free(ym1);
    gsl_vector_complex_free(ym2);
}

/**
 * \brief Matrix-matrix product when the matrix is given as a product of matrices: ym = DAT[1]...DAT[M] * xm, with
 *        ym  a 6x6 matrix
 *        xm  a 6x6 matrix
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_matrix_product(gsl_matrix_complex **DAT, const gsl_matrix_complex *xm, gsl_matrix_complex *ym, int M)
{
    gsl_complex one_c  = gslc_complex(1.0, 0.0);
    gsl_complex zero_c = gslc_complex(0.0, 0.0);

    gsl_matrix_complex *ym1 = gsl_matrix_complex_calloc(6,6);
    gsl_matrix_complex *ym2 = gsl_matrix_complex_calloc(6,6);

    //ym1 = xm
    gsl_matrix_complex_memcpy(ym1, xm);
    for(int i =1; i<=M;i++)
    {
        //ym2 = DAT[i]*ym1
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one_c , DAT[i], ym1 , zero_c , ym2);
        //ym1 = ym2
        gsl_matrix_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_matrix_complex_memcpy(ym, ym2);

    gsl_matrix_complex_free(ym1);
    gsl_matrix_complex_free(ym2);
}

/**
 * \brief Matrix inverse-vector product when the matrix is given as a product of matrices: ym = DAT[M]^(-1)...DAT[1]^(-1) * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_vector_invproduct(gsl_matrix_complex **DAT,  const gsl_vector_complex *xm, gsl_vector_complex *ym, int M)
{
    int s;
    gsl_permutation * p = gsl_permutation_alloc (6);
    gsl_matrix_complex *AUX = gsl_matrix_complex_calloc(6,6);
    gsl_vector_complex *ym1 = gsl_vector_complex_calloc(6);
    gsl_vector_complex *ym2 = gsl_vector_complex_calloc(6);

    //ym1 = xm
    gsl_vector_complex_memcpy(ym1, xm);
    for(int i =M; i>=1;i--)
    {
        //AUX = DAT(K). Here, memcpy is mandatory, otherwise the original DAT would be modified by the LU decomposition
        gsl_matrix_complex_memcpy(AUX, DAT[i]);
        //LU decomposition of AUX
        gsl_linalg_complex_LU_decomp (AUX, p, &s);
        //BUX = AUX^(-1)*VEP
        gsl_linalg_complex_LU_solve (AUX , p , ym1 , ym2);
        //ym1 = ym2;
        gsl_vector_complex_memcpy(ym1, ym2);
    }
    //ym = ym2
    gsl_vector_complex_memcpy(ym, ym2);

    gsl_vector_complex_free(ym1);
    gsl_vector_complex_free(ym2);
    gsl_matrix_complex_free(AUX);
    gsl_permutation_free(p);
}

/**
 * \brief Initialize and return a matrix as a product of complex matrices (DAT)
 *        CAREFUL: the array of matrices DAT is shifted of one & so that the storage is easier in other routines (e.g. vepro)
 *        ==> DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 **/
gsl_matrix_complex ** gslc_matrix_complex_product_alloc(int size1, int size2, int M)
{
    gsl_matrix_complex **DAT    = (gsl_matrix_complex**) malloc((M+1)*sizeof(gsl_matrix_complex*));

    for(int i = 0; i<= M; i++)
    {
        DAT[i] = gsl_matrix_complex_calloc(size1,size2);
        gsl_matrix_complex_set_zero(DAT[i]);
    }

    return DAT;
}

/**
 * \brief Free a matrix as a product of complex matrices (DAT)
 *        CAREFUL: the array of matrices DAT is shifted of one & so that the storage is easier in other routines (e.g. vepro)
 *        ==> DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 **/
void gslc_matrix_complex_product_free(gsl_matrix_complex ** P, int M)
{
    for(int i = 0; i<= M; i++) gsl_matrix_complex_free(P[i]);
}


/**
 * \brief Initialize and return a matrix as a product of real (DAT)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
gsl_matrix ** gslc_matrix_array_alloc(int size1, int size2, int M)
{
    gsl_matrix **DAT = (gsl_matrix**) malloc((M)*sizeof(gsl_matrix*));

    for(int i = 0; i< M; i++)
    {
        DAT[i] = gsl_matrix_calloc(size1,size2);
        gsl_matrix_set_zero(DAT[i]);
    }

    return DAT;
}

/**
 * \brief Free a matrix as a product of real (DAT)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
void gslc_matrix_array_free(gsl_matrix ** P, int M)
{
    for(int i = 0; i< M; i++) gsl_matrix_free(P[i]);
}





