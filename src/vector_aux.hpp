#ifndef _VEC_AUX_H_
#define _VEC_AUX_H_

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>

using namespace std;

double remainder(double quot, double divid);

void lin_sys_solv_LU(vector<vector<double>> A, vector<double> &x, vector<double> b);
void lin_sys_solv_SV(vector<vector<double>> A, vector<double> &x, vector<double> b);

void print_vector_i(vector<int>);
void print_vector_f(vector<double>);
void print_vector_c(vector<complex<double>>);
void print_matrix_i(vector<vector<int>>);
void print_matrix_f(vector<vector<double>>);
void print_matrix_c(vector<vector<complex<double>>>);

int kron_delta_i(int i, int j);
double kron_delta_f(int i, int j);
double norm(vector<double> v);
void normalize(vector<double> &v);

int matrix_rank_f(vector<vector<double>> M, double tol);
double trace_matrix_f(vector<vector<double>> A);
double deter_matrix_f(vector<vector<double>> A);
double prod_vector_f(vector<double> v1, vector<double> v2);
void eigen_stuff_f(vector<vector<double>> A, 
                   vector<complex<double>> &eig_vals,
                   vector<vector<complex<double>>> &eig_vecs);
vector<int> linspace_i(int xi, int xf = 1, unsigned int skip = 1);
vector<double> linspace_f(double xi, double xf, unsigned int N);
vector<double> arr_to_vec(double arr[], unsigned int n);
vector<double> vec_max(vector<double> v);
vector<double> vec_min(vector<double> v);
vector<double> vec_minmax(vector<double> v);
vector<double> complex_to_float(vector<complex<double>> v);
vector<complex<double>> float_to_complex(vector<double> v);
vector<double> scalar_vector_f(double c, vector<double> v);
vector<double> prod_matrix_vector_f(vector<vector<double>> M, vector<double> v);
vector<double> linear_comb_vector_f(double a1, vector<double> v1, double a2, vector<double> v2);
vector<vector<double>> diagonal_matrix(vector<vector<double>> M);
vector<vector<double>> transpose_matrix(vector<vector<double>> M);
vector<vector<double>> scalar_matrix_f(double c, vector<vector<double>> M);
vector<vector<double>> prod_matrix_f(vector<vector<double>> A, vector<vector<double>> B);
vector<vector<double>> prod_dyadic_f(vector<double> a, vector<double> b);
vector<vector<double>> linear_comb_matrix_f(double a1, vector<vector<double>> M1, 
                                            double a2, vector<vector<double>> M2);

#endif
