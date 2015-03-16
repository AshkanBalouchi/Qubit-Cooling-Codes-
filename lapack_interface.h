
#include <cstdlib>
#include <cmath> 
#include <vector>
#include <complex> 

extern "C" {
	// -------- LAPACK routines 
	// diagonalize symmetric real 
	void dsyev_(const char* choose, const char* uplow, long int* N, double* A, long int* rownum, double* Evals, double* WORK, long int* LWORK, long int* INFO );
	// diagonalize Hermitian 
	void zheev_(const char* choose, const char* uplow, long int* N, double* A, long int* rownum, double* Evals, double* WORK, long int* LWORK, double* RWORK, long int* INFO ); 
	// -------- BLAS routines 
	// matrix operating on a vector 
	void zgemv_(const char*, long int*, long int*, double*, double*, long int*, double*, long int*, double*, double*, long int*); 
	// matrix x matrix
	void zgemm_(const char*, const char*, long int*, long int*, long int*, double*, double*, long int*, double*, long int*, double*, double*, long int*); 
}

long int diagonalize(long int Dim, std::complex<double>* Mat, double* Evals);
long int tensor_prod(long int DimA, std::complex<double>* A, long int DimB, std::complex<double>* B, long int& DimC, std::complex<double>* C);
long int tprod_addto(long int DimA, std::complex<double>* A, long int DimB, std::complex<double>* B, std::complex<double>* C);
long int diagonalize_real(long int Dim, double* Mat, double* Evals); 
long int tensor_prod_real(long int DimA, double* A, long int DimB, double* B, long int& DimC, double* C);
long int tprod_addto_real(long int DimA, double* A, long int DimB, double* B, double* C);
long int A_on_v(    long int N, std::complex<double>* v, std::complex<double>* A); 
long int AT_on_v(   long int N, std::complex<double>* v, std::complex<double>* A); 
long int A_times_B( long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B); 
long int AT_times_B(long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B); 
long int A_times_BT(long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B); 
long int Acols_timesby_v( long int N, std::complex<double>* v, std::complex<double>* A); 
long int Arows_timesby_v( long int N, std::complex<double>* v, std::complex<double>* A);
long int copy_mat(long int N, std::complex<double>* A, std::complex<double>* B); 
std::complex<double> trace_mat(long int N, std::complex<double>* A); 
long int trace_bath(long int Ntot, long int Nsys, std::complex<double>* psi, std::complex<double>* rho); 
long int trace_out_sys( std::complex<double>* rho_in, std::complex<double>* rho_out, long int sys, std::vector<unsigned int>& dims, long int& new_size); 

