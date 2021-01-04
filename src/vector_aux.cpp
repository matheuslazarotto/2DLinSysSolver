#include "vector_aux.hpp"

double remainder(double quot, double divid)
{
    /** Beware of divisions between opposite 
     * sign values, it will return 0 **/
    return quot - floor(quot / divid) * divid;
}

void lin_sys_solv_LU(vector<vector<double>> A, vector<double> &x, vector<double> b) 
{
    /** Solve linear system:  Ax = b   by LU decomposition 
      * based on GLS Linear Algebra library 
      **/
    gsl_matrix *gslA;
    gsl_vector *gslx;
    gsl_vector *gslb;

    unsigned int size_b = b.size();
    unsigned int size_mA = A.size();
    unsigned int size_nA = A[0].size();
    
    // Allocation
    gslA = gsl_matrix_alloc(size_mA, size_nA);
    gslx = gsl_vector_alloc(size_nA);
    gslb = gsl_vector_alloc(size_b);

    // Set GSL matrices
    for (unsigned int i = 0; i < size_mA; i++) 
    {
        gsl_vector_set(gslb, i, b[i]);
        for (unsigned int j = 0; j < size_nA; j++)
        {
            gsl_matrix_set(gslA, i, j, A[i][j]);
        }
    }

    /*** Solve system ***/
    int gsl_signum;
    gsl_permutation *gsl_permA;
    gsl_permA = gsl_permutation_alloc(size_mA);

    // matrix-vector product: \vec{y} = alfa * A * vec{x} + beta * \vec{y}
    // gsl_blas_dgemv(CblasNoTrans, 1.0, gslA, gslx, 0.0, gslb);
    
    gsl_linalg_LU_decomp(gslA, gsl_permA, &gsl_signum);
    gsl_linalg_LU_solve(gslA, gsl_permA, gslb, gslx);

    // Assign 'gslx' to 'vector<double> *x'
    for (unsigned int i = 0; i < size_nA; i++) 
    { 
        x[i] = gsl_vector_get(gslx, i); 
    }

    gsl_permutation_free(gsl_permA);
    gsl_matrix_free(gslA);
    gsl_vector_free(gslb);
    gsl_vector_free(gslx);
}

void lin_sys_solv_SV(vector<vector<double>> A, vector<double> &x, vector<double> b) 
{
    /** Solve linear system:  Ax = b  by SVD 
     *  decomposition:
     * 
     *  A = U S V^T
      *
      *  based on GLS Linear Algebra library 
      **/
    gsl_matrix *gslA;
    gsl_vector *gslx;
    gsl_vector *gslb;

    unsigned int size_b = b.size();
    unsigned int size_mA = A.size();
    unsigned int size_nA = A[0].size();
    
    // Allocation
    gslA = gsl_matrix_alloc(size_mA, size_nA);
    gslx = gsl_vector_alloc(size_nA);
    gslb = gsl_vector_alloc(size_b);

    // Set GSL matrices
    for (unsigned int i = 0; i < size_mA; i++)
    {
        gsl_vector_set(gslb, i, b[i]);
        for (unsigned int j = 0; j < size_nA; j++)
        {
            gsl_matrix_set(gslA, i, j, A[i][j]);
        }
    }

    /*** Solve system ***/
    gsl_matrix *gslV;
	gsl_vector *gslS;
	gsl_vector *workspace;

	gslV = gsl_matrix_alloc(size_nA, size_nA);
	gslS = gsl_vector_alloc(size_nA);
	workspace = gsl_vector_alloc(size_nA);

    gsl_linalg_SV_decomp(gslA, gslV, gslS, workspace);
    // gslA matrix is changed to U under SVD decomposition
    gsl_linalg_SV_solve(gslA, gslV, gslS, gslb, gslx);

    // Assign 'gslx' to 'vector<double> *x'
    for (unsigned int i = 0; i < size_nA; i++) 
    { 
        x[i] = gsl_vector_get(gslx, i); 
    }

    gsl_matrix_free(gslA);
    gsl_matrix_free(gslV);
    gsl_vector_free(gslS);
    gsl_vector_free(gslb);
    gsl_vector_free(gslx);
    gsl_vector_free(workspace);
}

void print_vector_i(vector<int> v)
{
    cout << "[";
    for (unsigned int i = 0; i < v.size(); i++) 
    { 
        if (i == (v.size() - 1)) { cout << v[i]; }
        else { cout << v[i] << ", "; }
    }
    cout << "]\n" << endl;    
}

void print_vector_f(vector<double> v) 
{
    cout << "[";
    for (unsigned int i = 0; i < v.size(); i++) 
    { 
        if (i == (v.size() - 1)) { cout << v[i]; }
        else { cout << v[i] << ", "; }
    }
    cout << "]\n" << endl;
}

void print_vector_c(vector<complex<double>> v)
{
    cout << "[";
    for (unsigned int i = 0; i < v.size(); i++) 
    {
        if (i == (v.size() - 1)) { cout << v[i]; }
        else { cout << v[i] << ", "; }
    }
    cout << "]\n" << endl;
}

void print_matrix_i(vector<vector<double>> M)
{
    for (unsigned int i = 0; i < M.size(); i++) 
    {
        cout << "[";
        for (unsigned int j = 0; j < M[0].size(); j++) 
        {
            cout << M[i][j] << "   ";
        }
        cout << "\b\b" << "]\n" << endl;
    }    
}

void print_matrix_f(vector<vector<double>> M)
{
    for (unsigned int i = 0; i < M.size(); i++) 
    {
        cout << "[";
        for (unsigned int j = 0; j < M[0].size(); j++) 
        {
            cout << M[i][j] << "   ";
        }
        cout << "\b\b" << "]\n" << endl;
    }    
}

void print_matrix_c(vector<vector<complex<double>>> M)
{
    for (unsigned int i = 0; i < M.size(); i++) 
    {
        cout << "[";
        for (unsigned int j = 0; j < M[0].size(); j++) 
        {
            cout << M[i][j] << "   ";
        }
        cout << "\b\b" << "]\n" << endl;
    }    
}

int kron_delta_i(int i, int j)
{
    return ((i == j) ? 1 : 0); 
}

double kron_delta_f(int i, int j)
{
    return ((i == j) ? 1.0 : 0.0);
}

double norm(vector<double> v)
{
    double n = 0.0;
    
    for (auto vi: v) { n += vi * vi; }
    
    return sqrt(n);
}

void normalize(vector<double> &v)
{
    /* Normalize a vector */
    double n = norm(v);
    
    assert(n != 0.0);

    for (unsigned int k = 0; k < v.size(); k++) { v[k] = v[k] / n; }
}

int matrix_rank_f(vector<vector<double>> M, double tol)
{
	/** Performs a SDV (Singular Value Decomposition) 
	 *  on input matrix 'M'. The elements of the singular 
	 *  matrix S are counted, and the final rank corresponds 
	 *  to the number of non-zero elements, i.e. elements 
	 *  above 'tol'.
	 **/
	unsigned int rows = M.size();
	unsigned int cols = M[0].size();
	
	assert(rows > 0);
	assert(cols > 0);
	
	/** A general rectangular [m, n] matrix A has a singular 
	 * value decomposition into the product of an [m, n] 
	 * orthogonal matrix U, a [n, n] diagonal matrix of sin- 
	 * gular values S and the transpose of a [n, n] orthogo-
	 * nal square matrix V:
	 * 
	 * 	A = U S V^T
	 * 
	 * The singular values S[i][i] are all non-negative and 
	 * are generally chosen to form a non-increasing sequence.
	 **/
	
	gsl_matrix *gslM;
	gsl_matrix *gslV;
	gsl_vector *gslS;
	gsl_vector *workspace;
	
	gslM = gsl_matrix_alloc(rows, cols);
	gslV = gsl_matrix_alloc(cols, cols);
	gslS = gsl_vector_alloc(cols);
	workspace = gsl_vector_alloc(cols);
	
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			gsl_matrix_set(gslM, i, j, M[i][j]);
		}
	}
	
	/* GSL restraints the calculation for [m >= n] */
	assert(rows >= cols);
	/* After decomposed, matrix 'gslM' is altered to 'U' */
	gsl_linalg_SV_decomp(gslM, gslV, gslS, workspace);
	
	/* Rank = (# of singular values above 'tol') */
	unsigned int rank = 0;
	for (unsigned int i = 0; i < rows; i++)
	{
		if (gsl_vector_get(gslS, i) > tol) { rank += 1; }
	}
	
	gsl_matrix_free(gslM);
	gsl_matrix_free(gslV);
	gsl_vector_free(gslS);
	gsl_vector_free(workspace);
	
	return rank;
}

double trace_matrix_f(vector<vector<double>> M)
{
    double Tr = M[0][0];
    unsigned int rows = M.size();    
    unsigned int cols = M[0].size();

    assert(rows == cols); // !Trace of non square Matrix can not be taken!
    assert(rows != 0);     // !Trace of null dimension Matrix can not be taken!
    assert(cols != 0);     // !Trace of null dimension Matrix can not be taken!
    
    for (unsigned int i = 1; i < rows; i++) 
    {
        Tr += M[i][i];
    }

    return Tr;
}

double deter_matrix_f(vector<vector<double>> M)
{
    double det;
    unsigned int rows = M.size();
    unsigned int cols = M[0].size();

    assert(rows == cols); // !Determinant of non square Matrix can not be taken!

    // Allocation 
    gsl_matrix *gslM;
    gslM = gsl_matrix_alloc(rows, cols);

    // Set GSL matrix
    for (unsigned int i = 0; i < rows; i++) 
    {
        for (unsigned int j = 0; j < cols; j++) 
        {
            gsl_matrix_set(gslM, i, j, M[i][j]);
        }
    }

    // Evaluate determinant
    int gsl_signum;
    gsl_permutation *gslM_perm;
    gslM_perm = gsl_permutation_alloc(rows);

    gsl_linalg_LU_decomp(gslM, gslM_perm, &gsl_signum);
    det = gsl_linalg_LU_det(gslM, gsl_signum);

    gsl_permutation_free(gslM_perm);
    gsl_matrix_free(gslM);

    return det;
}

double prod_vector_f(vector<double> v1, vector<double> v2)
{
    // !Vector length incompatible with scalar product!
    assert(v1.size() == v2.size());
    
    double c = 0.0;
    for (unsigned int i = 0; i < v1.size(); i++) { c += v1[i] * v2[i]; }
    return c;    
}

void eigen_stuff_f(vector<vector<double>> M, 
                   vector<complex<double>> &eig_vals,
                   vector<vector<complex<double>>> &eig_vecs)
{
    gsl_matrix *gslM;
    gsl_matrix_complex *gslEig_vecs;
    gsl_vector_complex *gslEig_vals;
    gsl_eigen_nonsymmv_workspace *w;

    unsigned int rows = M.size();
    unsigned int cols = M[0].size();
    
    assert(rows == cols);  // !Eigenstuff can not be calculated for a non square Matrix!

    // Allocation
    w = gsl_eigen_nonsymmv_alloc(rows);
    gslM = gsl_matrix_alloc(rows, cols);
    gslEig_vals = gsl_vector_complex_alloc(rows);
    gslEig_vecs = gsl_matrix_complex_alloc(rows, cols);
    
    // Set GSL matrix
    for (unsigned int i = 0; i < rows; i++) 
    {
        for (unsigned int j = 0; j < cols; j++) 
        {
            gsl_matrix_set(gslM, i, j, M[i][j]);
        }
    }

    // Calculate eigenvalues and eigenvectors
    // (normalized eigenvectors)
    gsl_eigen_nonsymmv(gslM, gslEig_vals, gslEig_vecs, w);
    gsl_eigen_nonsymmv_sort(gslEig_vals, gslEig_vecs, GSL_EIGEN_SORT_ABS_DESC);

    for (unsigned int i = 0; i < rows; i++) 
    {
        complex<double> z(GSL_REAL(gsl_vector_complex_get(gslEig_vals, i)),
                          GSL_IMAG(gsl_vector_complex_get(gslEig_vals, i)));
        eig_vals[i] = z;
        for (unsigned int j = 0; j < rows; j++) 
        {
            complex<double> vec_z(GSL_REAL(gsl_matrix_complex_get(gslEig_vecs, j, i)),
                                  GSL_IMAG(gsl_matrix_complex_get(gslEig_vecs, j, i)));
            eig_vecs[j][i] = vec_z;
        }
    }
    
    gsl_eigen_nonsymmv_free(w);
    gsl_matrix_free(gslM);
    gsl_matrix_complex_free(gslEig_vecs);
    gsl_vector_complex_free(gslEig_vals);
}

vector<int> linspace_i(int xi, int xf, unsigned int skip) 
{
    vector<int> v;
    assert(xi <= xf);

    if (abs(xf - xi) < (int) skip) { v.push_back(xi); }
    else
    {
        for (int n = 0; n < abs(xf - xi) / (int) skip; n++)
        {
            v.push_back(xi + n * skip);
        }
    }

    return v;
}

vector<double> linspace_f(double xi, double xf, unsigned int N) 
{
    vector<double> v(N+1);

    for (int n = 0; n <= (int) N; n++)
    {
        v[n] = xi + ((double) n / (double) N) * (xf - xi);
    }

    return v;
}

vector<double> arr_to_vec(double arr[], unsigned int n)
{
    vector<double> vec(n);
    
    for (unsigned int k = 0; k < n; k++) { vec[k] = arr[k]; }

    return vec;
}

vector<double> vec_max(vector<double> v) 
{
    /** Returns a vector with the maximum element of a vector 
     * and its index position (in float units!) **/
    
    assert(v.size() > 0); // !Vector must have elements!

    double v_n = v[0];
    int i_max = 0;
    
    vector<double> c = {v_n, (double) i_max};

    for (unsigned int i = 1; i < v.size(); i++) 
    {
        v_n = v[i];
        
        if (v_n > c[0]) 
        { 
            c[0] = v_n;
            c[1] = (double) i;
        } 
    }

    return c;
}

vector<double> vec_min(vector<double> v)
{
    /** Returns a vector with the minimum element of a vector 
     * and its index position (in float units!) **/
    
    assert(v.size() > 0); // !Vector must have elements!
    
    double v_n = v[0];
    double min = v_n;
    int i_min = 0;

    vector<double> c = {v_n, (double) i_min};

    for (unsigned int i = 1; i < v.size(); i++) 
    {
        v_n = v[i];
        
        if (v_n < min) 
        { 
            c[0] = v_n;
            c[1] = (double) i;
        } 
    }

    return c;
}

vector<double> vec_minmax(vector<double> v) 
{
    /** Returns a vector with the maximum and minimum elements 
     * of a vector and its index position (in FLOAT units!) **/

    assert(v.size() > 0);   // !Vector must have elements!

    double v_n = v[0];
    double max = v_n;
    double min = v_n;
    int i_max = 0;
    int i_min = 0;

    vector<double> c = {v_n, (double) i_max, v_n, (double) i_min};

    for (unsigned int i = 0; i < v.size(); i++)
    {
        v_n = v[i]; 

        if (v_n > max) 
        {
            c[0] = v_n;
            c[1] = (double) i;
        } 
        else if (v_n < min) 
        {
            c[2] = v_n;
            c[3] = (double) i;
        }
    }

    /* c = {max, i_max, min, i_min} */
    return c;
}

vector<double> complex_to_float(vector<complex<double>> v)
{
    /** Returns the real part of a complex vector **/
    assert(v.size() > 0);   // !Vector must have non zero size!
    
    vector<double> v_real(v.size());
    
    for (unsigned int i = 0; i < v.size(); i++) { v_real[i] = v[i].real(); }

    return v_real;    
}

vector<complex<double>> float_to_complex(vector<double> v)
{
    /** Returns the a complex vector with the real part of a real vector **/
    assert(v.size() > 0);   // !Vector must have non zero size!

    vector<complex<double>> v_complex(v.size());

    for (unsigned int i = 0; i < v.size(); i++)
    {
        complex<double> val(v[i], 0.0);
        v_complex[i] = val;
    }

    return v_complex;
}

vector<double> scalar_vector_f(double c, vector<double> v)
{
    assert(v.size() > 0);       // !Vector must have non zero size!
    
    vector<double> r(v.size());
    for (unsigned int i = 0; i < v.size(); i++) { r[i] = c * v[i]; }
    
    return r;
}

vector<double> prod_matrix_vector_f(vector<vector<double>> M, vector<double> v)
{
    unsigned int size_rM = M.size();     /* rows */
    unsigned int size_cM = M[0].size();  /* columns */
    unsigned int size_v = v.size();      /* length */
    /* Output vector */
    vector<double> x(size_rM, 0.0);

    // !Matrix / vector dimension incompatible for product!
    assert(size_cM == size_v);
    
    for (unsigned int i = 0; i < size_rM; i++) 
    {
        for (unsigned int j = 0; j < size_cM; j++) 
        {
            x[i] += M[i][j] * v[j];
        }
    }
    return x;
}

vector<double> linear_comb_vector_f(double a1, vector<double> v1, double a2, vector<double> v2)
{
    // !Vector size incompatible for sum!
    assert(v1.size() == v2.size());
    
    vector<double> c(v1.size());
    for (unsigned int i = 0; i < v1.size(); i++) 
    { 
        c[i] = a1 * v1[i] + a2 * v2[i]; 
    }
    
    return c;
}

vector<vector<double>> diagonal_matrix(vector<vector<double>> M)
{
    assert(M.size() == M[0].size()); // !Matrix must be square shape!

    unsigned int size_M = M.size();     /* Rows|Columns */
    vector<vector<double>> D(size_M, vector<double>(size_M));

    for (unsigned int i = 0; i < size_M; i++)
    {
        for (unsigned int j = 0; j < size_M; j++)
        {
            if (i == j)
            {
                D[i][i] = M[i][i];
            }
            else
            {
                D[i][j] = 0.0;
            }
        }
    }

    return D;
}

vector<vector<double>> transpose_matrix(vector<vector<double>> M)
{
    unsigned int size_rM = M.size();     /* Rows */
    unsigned int size_cM = M[0].size();  /* Columns */
    vector<vector<double>> Mt(size_cM, vector<double>(size_rM));
    
    for (unsigned int i = 0; i < size_cM; i++)
    {
        for (unsigned int j = 0; j < size_rM; j++)
        {
            Mt[i][j] = M[j][i];
        }
    }
    
    return Mt;
}

vector<vector<double>> scalar_matrix_f(double c, vector<vector<double>> M)
{
    vector<vector<double>> R(M.size(), vector<double>(M[0].size()));

    for (unsigned int i = 0; i < M.size(); i++) 
    {
        for (unsigned int j = 0; j < M[0].size(); j++) 
        {
            R[i][j] = c * M[i][j];
        }
    }
    return R;
}

vector<vector<double>> prod_matrix_f(vector<vector<double>> A, 
                                     vector<vector<double>> B)
{
    unsigned int size_mA = A.size();     /* rows */
    unsigned int size_nA = A[0].size();  /* columns */
    unsigned int size_mB = B.size();     /* rows */
    unsigned int size_nB = B[0].size();  /* columns */
    /* output matrix */
    vector<vector<double>> C(size_mA, vector<double> (size_nB));

    // !Matrix dimension incompatible for product!
    assert(size_nA == size_mB);

    for (unsigned int i = 0; i < size_mA; i++) 
    {
        for (unsigned int j = 0; j < size_nB; j++) 
        {
            C[i][j] = 0.0;
            for (unsigned int k = 0; k < size_nA; k++) 
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<vector<double>> prod_dyadic_f(vector<double> a, vector<double> b)
{
    assert(a.size() > 0);   // !Vector must not have zero size!
    assert(b.size() > 0);   // !Vector must not have zero size!

    vector<vector<double>> M(a.size(), vector<double>(b.size()));
    for (unsigned int i = 0; i < a.size(); i++) 
    {
        for (unsigned int j = 0; j < b.size(); j++) 
        {
            M[i][j] = a[i] * b[j];
        }
    }

    return M;
}

vector<vector<double>> linear_comb_matrix_f(double a1, vector<vector<double>> M1, 
                                            double a2, vector<vector<double>> M2)
{
    //  !Matrix size incompatible for sum!
    assert(M1.size() == M2.size() && M1[0].size() == M2[0].size());
    
    vector<vector<double>> C(M1.size(), vector<double>(M1[0].size()));
    for (unsigned int i = 0; i < M1.size(); i++) 
    {
        for (unsigned int j = 0; j < M1[0].size(); j++) 
        {
            C[i][j] = a1 * M1[i][j] + a2 * M2[i][j];
        }
    }

    return C;
}
