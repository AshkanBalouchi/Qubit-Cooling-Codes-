
//  This file implements the classes for quantum states
//  and quantum operators. This includes the member functions
//  that allow the operators to act on the states, and the 
//  functions that give binary operations.
//
//  The classes are:      qstate
//                        qmatrix
//                        qoperator
//
//  The operators are:    qoperator: tensor,*,*(scalar),+,- 
//                        qstate:    tensor,*(scalar),+,>>,<<
//                        qmatrix:   tensor,*,*(scalar),+
//
//                        general: qoperator*qstate, 
//                                 qoperator*qmatrix, 
//                                 qmatrix*qoperator
//

#include "qobjects.h"

// Note: This class library includes code to use one of five external fft libraries.
//       This code is currently commented out, so that the user can use qobjects without 
//       having to install an fft library. In this case the meamber functions which 
//       invoke the fft's simply do nothing. Comments are included which direct 
//       the user as to which segments to uncomment so as to use one of the fft 
//       libraries. The current supported libraries are: GNU's GSL, netlibs Fftpack, 
//       SCS's Library (zzfft), FFTW and an fft written in f90 by Daniel Steck. This
//       last one can be requested from the author.  
//       -KJ

// To use GNU's GSL uncomment:
// #include <gsl/gsl_fft_complex.h>
// int fft(std::vector< std::complex<double> >& x, gsl_fft_direction);

// To use netlibs Fftpack uncomment:
// ------------------------------------------
// extern "C" {
// 
//   void zffti_(int* n, double* wsave); 
//   
//   void zfftf_(int* n, double* x, double* wsave); 
//   
//   void zfftb_(int* n, double* x, double* wsave);
// 
// }
//-------------------------------------------

// To use the SCS FFT uncomment:
// -------------------------------------------
//#include <scsl_fft.h>
// -------------------------------------------

// To use Dan Steck's fft1d uncomment:
// ------------------------------------------
//extern "C" {  void ffts1d_( int* fftsign, double* data, int* n ); }
// ------------------------------------------

// To use FFTW (Fastest FT in the West) uncomment the indicated lines in qobjects.h
// and the lines indicated in this file (below). To find them search for "FFTW"

// Diagonalisation using netlibs LAPACK 
// This not yet implemented here, but for future use, the calling interface is as follows 
// extern "C" {
//	 
//	 void ZHEEVD (char* JOBZ,  char* UPLO, int* N, double* A, int* LDA, 
//				  double* W, double* WORK, int* LWORK, double* RWORK, 
//				  int* LRWORK, int* IWORK, int* LIWORK, int* INFO)
// 
//     // note: DOUBLE COMPLEX A(LDAB,*), WORK(*)   
// }

//  ****************************
//  definitions for class qstate
//  ****************************

qstate::qstate()
{
	qcoeff.push_back(1);
	dims.push_back(1);
}

qstate::qstate(unsigned int dim)
{
	if (dim>0) {
		for (unsigned int i=0;i<dim;++i) {
			qcoeff.push_back(0);
		}
		dims.push_back(dim);
	}
}

qstate::qstate(std::string sname, unsigned int dim )
{
	if (dim>0) {
		bool match = false;
		if (sname == "down") {
			match = true;
			for (unsigned int i=1;i<=dim;++i) {
				if (i==dim) { qcoeff.push_back(1); }
				else { qcoeff.push_back(0); }
			}
			dims.push_back(dim);
		}
		if (sname == "ground" || sname == "up") {
			match = true;
			for (unsigned int i=1;i<=dim;++i) {
				if (i==1) { qcoeff.push_back(1); }
				else { qcoeff.push_back(0); }
			}
			dims.push_back(dim);
		}
		if (sname == "+") {
			match = true;
			for (unsigned int i=1;i<=dim;++i) {
				qcoeff.push_back(1/std::sqrt(1.0*dim));
			}
			dims.push_back(dim);
		}
		if (sname == "-") {
			match = true;
			int sign = 1;
			for (unsigned int i=1;i<=dim;++i) {
				qcoeff.push_back(sign/std::sqrt(1.0*dim));
				sign *= -1;
			}
			dims.push_back(dim);
		}
		if (sname == "y+") {
			match = true;
			dim = 2;
			std::complex<double> ei(0,1);
			qcoeff.push_back(-ei*std::sqrt(0.5));
			qcoeff.push_back(std::sqrt(0.5));
			dims.push_back(dim);	 
		}
		if (sname == "y-") {
			match = true;
			dim = 2;
			std::complex<double> ei(0,1);
			qcoeff.push_back(ei*std::sqrt(0.5));
			qcoeff.push_back(std::sqrt(0.5));
			dims.push_back(dim);	 
		}
		if (!match) { std::cout << "Error in qstate contructor: requested state not found " << std::endl; }
	}
}

qstate::qstate(std::string sname, std::complex<double> alpha, unsigned int dim)
{
	if (dim>0) {
		if (sname == "coherent") {
			std::complex<double> f_n(1.0,0.0);
			for (unsigned int i=0;i<dim;++i) {
				qcoeff.push_back(f_n);
				f_n = f_n*alpha/std::sqrt(i+1.0);
			}
			dims.push_back(dim);
			this->normalise();
		}
		if (sname == "number") {
			unsigned int number = static_cast<unsigned int>(std::floor(alpha.real()));
			for (unsigned int i=0;i<dim;++i) {
				if (i==number) { qcoeff.push_back(1.0); } 
				else { qcoeff.push_back(0.0); }
			}
			dims.push_back(dim);
		}
	}
	
	if ((dim>0)&&(dims.size()==0)) {
		std::cout << "Could not create state: requested state unknown" << std::endl;
	}
}

// std::vector< std::complex<double> > hermites(std::complex<double> x, int n) {
//	std::vector< std::complex<double> > h(n+1);
//   	h[0] = 1.0; if (n==0) {return h;}
//	h[1] = x;   if (n==1) {return h;} 
//	if (n>1) { for (unsigned int m=2;m<=n;m++) { h[m] = x*h[m-1] - 2.0*(m-1.0)*h[m-2]; std::cout << h[m] << std::endl; } }
//    return h;
// }
// 
// std::vector<double> factorial(int n) {
//  std::vector<double> fac(n+1);
//  fac[0] = 1; if (n==0) {return fac;}
//  	if (n>0) { for (unsigned int m=1;m<=n;m++) { fac[m] = m*fac[m-1];  std::cout << fac[m] << std::endl;} }
//    return fac;
// }

std::vector< std::complex<double> > HOSFSP(std::complex<double> x, std::complex<double> p, int n) {  
    // First n Hermites * Sqrt( P^n / n! )  
	std::vector< std::complex<double> > hof(n+1);
   	hof[0] = 1.0; if (n==0) {return hof;} 
	hof[1] = x*std::sqrt(p);   if (n==1) {return hof;} 
	if (n>1) { 
		for (unsigned int m=2;m<=n;m++) { 
			hof[m] = x*std::sqrt(p/(1.0*m))*hof[m-1] - 2.0*p*(m-1.0)*hof[m-2]/( std::sqrt(1.0*m*(m-1)) ); 
		} 
	}
    return hof; 
}

qstate::qstate(std::string sname, std::complex<double> alpha, std::complex<double> squeeze, unsigned int dim)
{
	if (dim>0) {                                           // Formula is in AJP 58, 1003 (1990) 
		if (sname == "squeezed" && std::abs(squeeze)>0.0) {
			std::complex<double> ei(0,1);
			double r = std::abs(squeeze);
			double theta = std::arg(squeeze);
			std::complex<double> gamma = alpha*std::cosh(r) + std::conj(alpha)*std::exp(ei*theta)*std::sinh(r);
			std::complex<double> x = gamma*std::exp(-ei*theta/2.0)/std::sqrt(std::sinh(2*r));				
			std::vector< std::complex<double> > hof = HOSFSP(x, std::exp(ei*theta)*std::tanh(r), dim-1); 
			for (int m=0;m<dim;m++) { qcoeff.push_back(hof[m]); }
			double norm = 0; for (int i=0;i<qcoeff.size();i++) {norm += std::abs(qcoeff[i]*std::conj(qcoeff[i]));}
			for (int i=0;i<qcoeff.size();i++) {qcoeff[i] /= std::sqrt(norm);}
			dims.push_back(dim);
		}
	}
	if ((dim>0)&&(dims.size()==0)) { 
		std::cout << "Could not create state: requested state unknown or parameter out of bounds" << std::endl;
	}
}

qstate::qstate( std::vector< std::complex<double> >& vec )
{
	dims.push_back(vec.size());
	qcoeff = vec;
}

qstate::qstate(std::string sname, unsigned int nx, 
			   const double& L,      const double& hbar,
			   const double& Vx0,    const double& Vp0,
			   const double& meanx0, const double& meanp0 )
{
	if (nx>0) {
		if (sname == "Gaussian") {
			double dx = L/nx;
			double pi = 3.1415927;
			//double e  = 2.7182818;
			std::complex<double> ei(0,1);
			double xi = 0;
			double norm0 = std::pow((2*pi*Vx0),-0.25)*std::sqrt(dx);
			double zeta0 = std::sqrt(std::fabs(Vp0 - hbar*hbar/(4*Vx0)));
			
			for (unsigned int i=0;i<nx;++i) {
				xi = -L/2 + i*dx;
				qcoeff.push_back( norm0*std::exp(- (xi - meanx0)*(xi - meanx0)/(4*Vx0)
												 + ei*meanp0*(xi + meanx0)/hbar
												 + ei*zeta0*(xi*xi)/(2*hbar*sqrt(Vx0))
												 )               );  
			}
			dims.push_back(nx);
		}
	}
	if ((nx>0)&&(dims.size()==0)) {
		std::cout << "Could not create state: requested state unknown" << std::endl;
	}
}

int qstate::clear()
{
	for (unsigned int i=0;i<qcoeff.size();++i) {
		qcoeff[i] = 0;
	}
	return 0;
}

int qstate::assign( std::vector< std::complex<double> >& vec)
{
	dims.assign(1,vec.size());
	qcoeff = vec;
	return 0;
}

int qstate::assign( std::complex<double>* start, std::complex<double>* finish)
{
	qcoeff.assign(start,finish);
	dims.assign(1,qcoeff.size());
	return 0;
}

int qstate::writeblock(int block) const
{
	if (dims.size()>1) {
		for (unsigned int i=0;i<dims[0];++i) {
			for (unsigned int j=0;j<dims[1];++j) {
				std::cout << qcoeff[i + dims[0]*j + (block-1)*dims[0]*dims[1]] << "  ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	else {
		for (unsigned int i=0;i<dims[0];++i) {
			std::cout << qcoeff[i] << std::endl;
		}
		std::cout << std::endl;
	}
	return 0;
}

double qstate::norm() const
{
	return ((*this).inner_prod(*this)).real();
}

qmatrix qstate::trace_out( int sys ) const 
{
	
	std::vector<unsigned int> new_dims(dims.size()-1); 
	for (int i=0;i<dims.size();i++) { 
		if (i < sys-1) { new_dims[i]=dims[i]; }
		if (i > sys-1) { new_dims[i-1] = dims[i]; }
	}
	
	int new_size(1);
	for (int i=0;i<new_dims.size();i++) { new_size *= new_dims[i]; }
	
	qmatrix result(new_size,new_dims); 
	
	int chunk(1); for (int i=0;i<sys-1;i++) { chunk *= dims[i]; } 
	int dim(dims[sys-1]); 
	int num(1);  for (int i=sys;i<dims.size();i++) { num *= dims[i]; } 
	
	for (int i=0;i<dim;i++) {
		
		for (int j1=0;j1<num;j1++) {
			for (int j2=0;j2<num;j2++) {
				
				for (int k1=0;k1<chunk;k1++) {
					for (int k2=0;k2<chunk;k2++) {
						
						result.addto_elem(k1 + j1*chunk, k2 + j2*chunk, 
										  qcoeff[k1 + i*chunk + j1*chunk*dim]*std::conj(qcoeff[k2 + i*chunk + j2*chunk*dim]) );
						
					}}
				
			}}
		
	}
	
	return result;
}

qmatrix qstate::trace_all_but( int sys ) const 
{
	
	std::vector<unsigned int> new_dims(1); 
	new_dims[0] = dims[sys-1]; 
		
	int new_size(new_dims[0]);
		
	qmatrix result(new_size,new_dims); 
	
	int chunk(1); for (int i=0;i<sys-1;i++) { chunk *= dims[i]; }  
	int num(1);  for (int i=sys;i<dims.size();i++) { num *= dims[i]; } 
	
	for (int i=0;i<chunk;i++) { 
		for (int j=0;j<num;j++) { 
			
			for (int k1=0;k1<new_size;k1++) { 
				for (int k2=0;k2<new_size;k2++) { 
					
					std::complex<double> val =             qcoeff[i + k1*chunk + j*chunk*new_size] 
											   *std::conj( qcoeff[i + k2*chunk + j*chunk*new_size] );
					result.addto_elem(k1,k2, val );
										
				}
			}
						
		}
	}
	
	return result;	
}


std::complex<double> qstate::get_elem( unsigned int i ) const
{
	return qcoeff[i];
}

int qstate::addto_elem( unsigned int i, std::complex<double> val )
{
    qcoeff[i] += val;
	return 0;
}

void* qstate::get_array() 
{
	//return reinterpret_cast<void*>(&(qcoeff[0])); 
	return &(qcoeff[0]);
}

unsigned int qstate::len() const
{
	return qcoeff.size(); 
}

std::vector<unsigned int> qstate::get_dims() const
{
	return dims;
}

int qstate::fftshift()    // thie fftshifts subsystem 1  
{
	unsigned int n = qcoeff.size();
	unsigned int half = n/2;
	std::vector< std::complex<double> > temp(qcoeff);  
	for (int i=0;i<half;i++) { qcoeff[i] = temp[i + half]; }
	for (int i=0;i<half;i++) { qcoeff[i + half] = temp[i]; }
	return 0;
}

int qstate::normalise() 
{
	double scale = std::sqrt((*this).norm());
	// std::cout << "Cannot normalise state: norm is zero" << std::endl;
	for (unsigned int i=0;i<qcoeff.size();++i) {
		qcoeff[i] = qcoeff[i]/scale;
	}
	return 0;
}

std::complex<double> qstate::inner_prod(const qstate& psi) const
{
	std::complex<double> result = 0;   
	if (qcoeff.size() == psi.qcoeff.size()) {
		for (unsigned int i=0;i<qcoeff.size();++i) {
			result += std::conj(qcoeff[i])*psi.qcoeff[i];
		}
	}
	else {
		std::cout << "Cannot take inner product: state size mismatch\n";
		return 0;
	}
	return result;
}

std::complex<double> inner_prod(const qstate& psi, const qstate& chi)
{
	return psi.inner_prod(chi);
}

double qstate::parity() const
{
	unsigned int nx = qcoeff.size();
	qstate psi_temp(*this), psi_refl(nx); 
	psi_temp.qcoeff[0] = 0.0;    // xtra point  (x = -L)
	psi_temp.qcoeff[nx/2] = 0.0; // grid center (x = 0 ) 
	psi_refl = psi_temp; 
	for (unsigned int i=1;i<nx;i++) { psi_refl.qcoeff[i] = psi_temp.qcoeff[nx-i]; } 
	return 2*(psi_temp.inner_prod(psi_refl)).real() - 1; 
}

int qstate::apply(const qoperator& op) 
{
	qstate temp(*this);
	clear();
	if (op.matdim()==temp.qcoeff.size()) {
		for (unsigned int i=0;i<op.size();++i) {
			qcoeff[op.getrow(i)] += op.getval(i)*temp.qcoeff[op.getcol(i)];
		}
	}
	else {
		std::cout << "Cannot apply operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qstate::apply(const qoperator& op, std::complex<double> cx) 
{
	qstate temp(*this);
	clear();
	if (op.matdim()==temp.qcoeff.size()) {
		for (unsigned int i=0;i<op.size();++i) {
			qcoeff[op.getrow(i)] += op.getval(i)*cx*temp.qcoeff[op.getcol(i)];
		}
	}
	else {
		std::cout << "Cannot apply operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qstate::apply_Hconj(const qoperator& op) 
{
	qstate temp(*this);
	clear();
	if (op.matdim()==temp.qcoeff.size()) {
		for (unsigned int i=0;i<op.size();++i) {
			qcoeff[op.getcol(i)] += std::conj(op.getval(i))*temp.qcoeff[op.getrow(i)];
		}
	}
	else {
		std::cout << "Cannot apply operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qstate::apply_diag( std::complex<double>* vec, int dim )
{
	int state_dim = 1; 
	for (unsigned int i=0;i<dims.size();i++) { state_dim *= dims[i]; } 
	if (dim == state_dim) { 
		std::complex<double> *s_ptr = &(qcoeff[0]), *fin = &(qcoeff[dim-1]); fin++; 
		std::complex<double>* v_ptr = vec; 
		while (s_ptr != fin) { *(s_ptr++) *= *(v_ptr++); }
	} else { std::cout << "Cannot apply diag vector: size mismatch\n" << std::endl; } 
	return 0; 
}

int	qstate::apply_FORT_cplx_mat( std::complex<double>* mat, int dim )   // requires lapack 
{
	int state_dim = 1; 
	for (unsigned int i=0;i<dims.size();i++) { state_dim *= dims[i]; } 
	if (dim == state_dim) { 
		std::complex<double>* state_array = &(qcoeff[0]); 
		long int Dim = dim; 
		A_on_v(Dim, state_array, mat); // calls lapack 
	} else { std::cout << "Cannot apply FORT matrix: size mismatch\n" << std::endl; } 
	return 0; 
}

int	qstate::applyT_FORT_cplx_mat( std::complex<double>* mat, int dim )  // requires lapack
{
	int state_dim = 1; 
	for (unsigned int i=0;i<dims.size();i++) { state_dim *= dims[i]; } 
	if (dim == state_dim) { 
		std::complex<double>* state_array = &(qcoeff[0]); 
		long int Dim = dim; 
		AT_on_v(Dim, state_array, mat); // calls lapack 
	} else { std::cout << "Cannot apply FORT matrix: size mismatch\n" << std::endl; } 
	return 0; 
}

const qstate qstate::operator+(const qstate& psi) const
{
	qstate result(*this);
	if (result.qcoeff.size()==psi.qcoeff.size()) {
		for (unsigned int i=0;i<result.qcoeff.size();++i) {
			result.qcoeff[i] += psi.qcoeff[i];
		}
	}
	else {
		std::cout << "Cannot add states: size mismatch\n";
	}
	return result;
}

const qstate qstate::operator-(const qstate& psi) const
{
	qstate result(*this);
	if (result.qcoeff.size()==psi.qcoeff.size()) {
		for (unsigned int i=0;i<result.qcoeff.size();++i) {
			result.qcoeff[i] -= psi.qcoeff[i];
		}
	}
	else {
		std::cout << "Cannot add states: size mismatch\n";
	}
	return result;
}

const qstate& qstate::operator=(const qstate& psi) 
{
	if (this->qcoeff.size()==psi.qcoeff.size()) {
		for (unsigned int i=0;i<psi.qcoeff.size();++i) {
			this->qcoeff[i] = psi.qcoeff[i];
		}
		for (unsigned int i=0;i<psi.dims.size();++i) {
			this->dims[i] = psi.dims[i];
		}
	}
	else {
		this->qcoeff = psi.qcoeff;
		this->dims = psi.dims;
	}
	return *this;
}

const qstate& qstate::operator+=(const qstate& psi)
{
	if (this->qcoeff.size()==psi.qcoeff.size()) {
		for (unsigned int i=0;i<psi.qcoeff.size();++i) {
			this->qcoeff[i] += psi.qcoeff[i];
		}
	}
	else {
		std::cout << "Cannot add to (+=) the state: size mismatch\n";
	}
	return *this;
}

const qstate& qstate::operator*=(std::complex<double> z)
{
	for (unsigned int i=0;i<qcoeff.size();++i) { qcoeff[i] *= z; }
	return *this;
}

const qstate& qstate::operator*=(double c)
{
	for (unsigned int i=0;i<qcoeff.size();++i) { qcoeff[i] *= c; }
	return *this;
}

const qstate qstate::operator*(const std::complex<double>& cx) const
{
	qstate result = *this;
	for (unsigned int i=0;i<result.qcoeff.size();++i) {
		result.qcoeff[i] *= cx;
	}
	return result;
}

const qstate operator*(const std::complex<double>& cx, const qstate& psi)
{
	return psi*cx;
}

const qstate qstate::operator*(const double& cx) const
{
	qstate result = *this;
	for (unsigned int i=0;i<result.qcoeff.size();++i) {
		result.qcoeff[i] *= cx;
	}
	return result;
}

const qstate operator*(const double& cx, const qstate& psi)
{
	return psi*cx;
}

int qstate::tensor(const qstate& chi)
{
	qstate psi = *this; 
	qstate result(0);
	if ((chi.qcoeff.size()>0)&&(psi.qcoeff.size()>0)){
		for (unsigned int i=0;i<chi.qcoeff.size();++i) {
			for (unsigned int j=0;j<psi.qcoeff.size();++j) {
				result.qcoeff.push_back(psi.qcoeff[j]*chi.qcoeff[i]);
			}
		}
		for (unsigned int i=0;i<psi.dims.size();++i) {
			if (psi.dims[i]>1) {result.dims.push_back(psi.dims[i]);}     // dims[0] is the fastest index 
		}
		for (unsigned int i=0;i<chi.dims.size();++i) {
			if (chi.dims[i]>1) {result.dims.push_back(chi.dims[i]);}
		}
	}
	*this = result; 
	return 0; 
}

qstate tensor(const qstate& psi, const qstate& chi)
{
	qstate result = psi;
	result.tensor(chi);
	return result;
}

qstate tensor(const qstate& psi1, const qstate& psi2, const qstate& psi3)
{
	qstate result = tensor(psi1,psi2);
	result.tensor(psi3);
	return result;
}

int qstate::shift( const int nn ) 
{
	if (nn!=0) {
		unsigned int nx = qcoeff.size();
		if (nn>0) {
			unsigned int n = static_cast<unsigned int>(nn);
			if (n>nx) {std::cout << " shift too big! " << std::endl;}
			if (n==nx) {std::cout << " warning: shifted the full length -> all elements are now zero" << std::endl;}
			for (unsigned int j=0;j<(nx-n);++j) {qcoeff[j] = qcoeff[j+n];}
			for (unsigned int j=(nx-n);j<nx;++j) {qcoeff[j] = 0;}
		}
		else {
			unsigned int n = static_cast<unsigned int>(-nn);
			if (n>nx) {std::cout << " shift too big! " << std::endl;}
			if (n==nx) {std::cout << " warning: shifted the full length -> result is the zero vector" << std::endl;}
			for (unsigned int j=nx-1;j>=n;j--) { qcoeff[j] = qcoeff[j-n];}
			for (unsigned int j=n;j>=1;j--) { qcoeff[j-1] = 0;}           // cant condition on <zero!
		}
	}
	return 0;
}

int qstate::rotate( const int nn )
{
	if (nn!=0) {
		unsigned int nx = qcoeff.size();
		if (nn>0) {
			unsigned int n = static_cast<unsigned int>(nn);
			if (n>nx) {std::cout << " rotation too big! " << std::endl;}
			std::vector< std::complex<double> > temp;
			for (unsigned int j=0;j<n;++j) {temp.push_back(qcoeff[j]);}
			for (unsigned int j=0;j<(nx-n);++j) {qcoeff[j] = qcoeff[j+n];}
			for (unsigned int j=nx-n;j<nx;++j) {qcoeff[j] = temp[j-nx+n];} 
		}
		else {
			unsigned int n = static_cast<unsigned int>(-nn);
			if (n>nx) {std::cout << " rotation too big! " << std::endl;}
			std::vector< std::complex<double> > temp;
			for (unsigned int j=nx-n;j<nx;++j) {temp.push_back(qcoeff[j]);}
			for (unsigned int j=nx-1;j>=n;--j) {qcoeff[j] = qcoeff[j-n];}
			for (unsigned int j=0;j<n;++j) {qcoeff[j] = temp[j];}
		}
	}
	return 0;
}

// note: this assumes the system consists of qubits only
int qstate::cnot( const int control, const int target )
{
	std::complex<double> ztemp;
	
	int t_stride = static_cast<int>(std::ceil(std::pow(2.0,target)));
	int t_blocks = qcoeff.size()/(2*t_stride);
	std::vector<int> targetbits(qcoeff.size());
	for (int i=0;i<t_blocks;i++) { 
		for (int j=2*i*t_stride;j<2*i*t_stride+t_stride;j++) { // loop over target 1's
			targetbits[j] = 1;
		}
	}
	
	int c_stride = static_cast<int>(std::ceil(std::pow(2.0,control)));
	int c_blocks = qcoeff.size()/(2*c_stride);
	for (int i=0;i<c_blocks;i++) { 
		for (int j=2*i*c_stride;j<2*i*c_stride+c_stride;j++) {  // loop over control 1's
			if (targetbits[j]==1) { 
				ztemp = qcoeff[j];
				qcoeff[j] = qcoeff[j + t_stride];
				qcoeff[j + t_stride] = ztemp;
			}
		}
	}
	
	return 0;
}

// note: this assumes the system consists of qubits only
int qstate::flip( const int target )
{
	std::complex<double> ztemp;
	
	int t_stride = static_cast<int>(std::ceil(std::pow(2.0,target)));
	int t_blocks = qcoeff.size()/(2*t_stride);
	std::vector<int> targetbits(qcoeff.size());
	for (int i=0;i<t_blocks;i++) { 
		for (int j=2*i*t_stride;j<2*i*t_stride+t_stride;j++) { // loop over target 1's
			targetbits[j] = 1;
		}
	}
	
	for (int j=0;j<qcoeff.size();j++) {  // loop over control 1's
		if (targetbits[j]==1) { 
			ztemp = qcoeff[j];
			qcoeff[j] = qcoeff[j + t_stride];
			qcoeff[j + t_stride] = ztemp;
		}
	}
	
	return 0;
}

// ----- to use FFTW uncomment: --------------------------------------
//  int qstate::fftw_init( fftw_plan& plan1, fftw_plan& plan2 )
// { 
//    int chunknum = qcoeff.size()/dims[0];
//    int n = static_cast<int>(dims[0]);
//    void* itt = &qcoeff[0];
//    
//    std::vector< std::complex<double> > temp = qcoeff;
//                                        
//    plan1 = fftw_plan_many_dft(1, &n, chunknum, reinterpret_cast<fftw_complex*>(itt), &n, 1, n, 
//                                                reinterpret_cast<fftw_complex*>(itt), &n, 1, n, 
//                                                -1, FFTW_MEASURE);
//    plan2 = fftw_plan_many_dft(1, &n, chunknum, reinterpret_cast<fftw_complex*>(itt), &n, 1, n, 
//                                                reinterpret_cast<fftw_complex*>(itt), &n, 1, n, 
//                                                 1, FFTW_MEASURE);
//                                                                                  
//    for (int i=0;i<qcoeff.size();i++) { qcoeff[i] = temp[i]; }
//    return 0;
// }
// 
// int qstate::fftw_fin( fftw_plan& plan1, fftw_plan& plan2 )
// {
//    fftw_destroy_plan(plan1);
//    fftw_destroy_plan(plan2);
//    return 0;
// }
// 
// int qstate::fftw( fftw_plan& plan )
// {
//    fftw_execute(plan);
//    return 0;
// }
// ----------------------------------------------------------------------

int qstate::fft_init( double** init_array )
{
    std::vector< std::complex<double> >::iterator itt;
    int chunknum = qcoeff.size()/dims[0];
    int n = static_cast<int>(dims[0]);
	
	// -----to use netlibs Fftpack uncoment: ----------------------------------
	//       *init_array = new double[4*n + 15];
	//       zffti_(&n, *init_array);
	//-------------------------------------------------------------------------
	
	// -----to use SCS uncomment: -----------------------------------------------------------
	//       itt = qcoeff.begin();
	//       double scale = sqrt(1.0/n);
	//       int isys[2]; isys[0] = 1;
	//       *init_array = new double[2*n + 256];
	//       double* work = new double[2*n];
	//       zzfft(0, n, scale, reinterpret_cast<scsl_zomplex*>(itt), 
	//                          reinterpret_cast<scsl_zomplex*>(itt), *init_array, work, isys);
	// 
	//       delete[] work;
	// --------------------------------------------------------------------------------------  
	
	// -----Note: to use GNU's GSL nothing needs to be uncommented here -------
	
	// -----Note: to use Steck's FFT nothing needs to be uncommented here -----
	
    return 0;
}

int qstate::fft_ss1b(int sign, double* init_array)
{
    std::vector< std::complex<double> >::iterator itt;
    int chunknum = qcoeff.size()/dims[0];
    int n = static_cast<int>(dims[0]);
    for (int i=0;i<chunknum;++i) {
		itt = qcoeff.begin() + i*dims[0];
		
		// -----to use netlibs Fftpack uncoment: ----------------------------------
		//       if (sign==1) {
		//         zfftf_(&n, reinterpret_cast<double*>(itt), init_array);
		//       } else {
		//         zfftb_(&n, reinterpret_cast<double*>(itt), init_array);
		//       }
		//-------------------------------------------------------------------------
		
		// -----to use SCS uncomment: -----------------------------------------------------------
		//       double scale = sqrt(1.0/n);
		//       int isys[2]; isys[0] = 1;
		//       double* work = new double[2*n];
		//       zzfft(-sign, n, scale, reinterpret_cast<scsl_zomplex*>(itt), 
		//                              reinterpret_cast<scsl_zomplex*>(itt), init_array, work, isys);
		//       delete[] work;
		// --------------------------------------------------------------------------------------  
		
		// -----to use GNU's GSL uncomment: -----------------------------------------------------
		//      gsl_fft_complex_radix2_transform(reinterpret_cast<gsl_complex_packed_array>(itt),
		//                                       1,dims[0],static_cast<gsl_fft_direction>(sign));
		// --------------------------------------------------------------------------------------  
		
		// -----to use Dan Steck's fft1d uncomment:----------------------------------------------
		//       sign = -sign;
		//       ffts1d_(&sign, reinterpret_cast<double*>(itt), &n);
		// --------------------------------------------------------------------------------------   
		
    }
    return 0;
}

int qstate::fft_subsys1(int sign)
{
    std::vector< std::complex<double> >::iterator itt;
    int chunknum = qcoeff.size()/dims[0];
    int n = static_cast<int>(dims[0]);
    for (int i=0;i<chunknum;++i) {
		itt = qcoeff.begin() + i*dims[0];
		
		// -----to use netlibs Fftpack uncoment: ----------------------------------
		//       double* wsave = new double[4*n + 15];
		//       zffti_(&n, wsave);
		//       if (sign==1) {
		//         zfftf_(&n, reinterpret_cast<double*>(itt), wsave);
		//       } else {
		//         zfftb_(&n, reinterpret_cast<double*>(itt), wsave);
		//       }
		//       delete[] wsave;
		//-------------------------------------------------------------------------
		
		// -----to use SCS uncomment: -----------------------------------------------------------
		//       double scale = sqrt(1.0/n);
		//       int isys[2]; isys[0] = 1;
		//       double* table = new double[2*n + 256];
		//       double* work = new double[2*n];
		//       zzfft(0, n, scale, reinterpret_cast<scsl_zomplex*>(itt), 
		//                          reinterpret_cast<scsl_zomplex*>(itt), table, work, isys);
		//       zzfft(-sign, n, scale, reinterpret_cast<scsl_zomplex*>(itt), 
		//                              reinterpret_cast<scsl_zomplex*>(itt), table, work, isys);
		//       delete[] work;
		//       delete[] table;
		// --------------------------------------------------------------------------------------  
		
		// ----to use GNU's GSL uncomment: ------------------------------------------------------
		//      gsl_fft_complex_radix2_transform(reinterpret_cast<gsl_complex_packed_array>(itt),
		//                                       1,dims[0],static_cast<gsl_fft_direction>(sign));
		// --------------------------------------------------------------------------------------  
		
		// -----to use Dan Steck's fft1d uncomment:----------------------------------------------
		//       sign = -sign;
		//       ffts1d_(&sign, reinterpret_cast<double*>(itt), &n);
		// --------------------------------------------------------------------------------------
		
    }
    return 0;
}

std::ostream& qstate::output(std::ostream& s) const 
{
	s << dims.size() << std::endl;
	for (unsigned int i=0;i<dims.size();i++) {s << dims[i] << " ";}
	s << std::endl;
	s << qcoeff.size() << std::endl;
	unsigned int elems_per_line = 4;
	div_t qr = std::div(static_cast<int>(qcoeff.size()),static_cast<int>(elems_per_line));
	unsigned int fullrows = static_cast<unsigned int>(qr.quot);
	unsigned int lastrown = static_cast<unsigned int>(qr.rem);
	for (unsigned int i=0;i<fullrows;i++) {
		for (unsigned int j=0;j<elems_per_line;j++) {
			s << qcoeff[i*elems_per_line + j] << " ";
		}
		s << std::endl;
	}
	for (unsigned int i=0;i<lastrown;i++) {s << qcoeff[fullrows*elems_per_line + i] << " ";}
	s << std::endl;
	return s;
}

std::ostream& operator<<(std::ostream& s, const qstate& psi) 
{
	psi.output(s);
	return s;
}

std::istream& qstate::input(std::istream& s)
{
	qcoeff.assign(1,0); qcoeff.pop_back();     // clear the vectors
	dims.assign(1,0); dims.pop_back();
	unsigned int dims_size, dims_elem;
	s >> dims_size;
	for (unsigned int i=0;i<dims_size;i++) {
		s >> dims_elem;
		dims.push_back(dims_elem);
	}
	unsigned int qcoeff_size;
	std::complex<double> qcoeff_elem;
	s >> qcoeff_size;
	for (unsigned int i=0;i<qcoeff_size;i++) {
		s >> qcoeff_elem;
		qcoeff.push_back(qcoeff_elem);
	}
	return s;
}

std::istream& operator>>(std::istream& s, qstate& psi)
{
	psi.input(s);
	return s;
}

int qstate::output(std::ofstream& out_real, std::ofstream& out_imag) const 
{
	for (int i=0;i<qcoeff.size();i++) { out_real << qcoeff[i].real() << " ";}
	out_real << std::endl;
	for (int i=0;i<qcoeff.size();i++) { out_imag << qcoeff[i].imag() << " ";}
	out_imag << std::endl;   
	return 0;
}

int qstate::output_dist(std::ofstream& outfile, unsigned int subsys)
{
	if ((subsys<=dims.size())&&(subsys>0)) {
		int chunksize = 1;
		for (unsigned int i=0;i<(subsys-1);++i) {
			chunksize = chunksize*dims[i];
		}
		int chunknum = 1;
		for (unsigned int i=subsys;i<dims.size();++i) {
			chunknum = chunknum*dims[i];
		}
		int chunksep = chunksize*dims[subsys-1];
		std::complex<double> zval;
		double val;
		for (unsigned int i=0;i<dims[subsys-1];++i) {
			val = 0;
			for (int j=0;j<chunknum;++j) {
				for (int k=0;k<chunksize;++k) {
					zval = qcoeff[k + i*chunksize + j*chunksep];
					val += zval.real()*zval.real() + zval.imag()*zval.imag();
				}
			}
			outfile << val << " ";
		}
		outfile << std::endl;
	}
	else {
		std::cout << dims.size() << std::endl;
		std::cout << "Cannot calculate distribution: subsystem does not exist" << std::endl;
		return 1;
	}
	return 0;
}

int qstate::display( std::string name ) const
{
	std::cout << name << std::endl;
	if (dims.size()>2) {
		int blocktotal = 1;
		for (unsigned int i=2;i<dims.size();++i) { 
			blocktotal *= dims[i];
		}
		for (int block=1;block<=blocktotal;++block) {
			writeblock(block);
		}
	} 
	else {
		writeblock(1);
	}
	return 0;
}

int qstate::display1D( std::string name ) const
{
	std::cout << name << std::endl;
	for (int i=0;i<qcoeff.size();i++) {
		std::cout << qcoeff[i] << std::endl;
	}
	return 0;
}

int showstate(const qstate& psi, std::string name )
{
	psi.display(name);
	return 0;
}


//  *****************************
//  definitions for class qmatrix
//  ***************************** 

qmatrix::qmatrix(unsigned int dim)
{
	if (dim>0) {
		std::vector< std::complex<double> > temp;
		for (unsigned int i=0;i<dim;++i) {temp.push_back(0);}
		for (unsigned int j=0;j<dim;++j) {elements.push_back(temp);}
		dims.push_back(dim);
	}
}

qmatrix::qmatrix( unsigned int dim, std::vector<unsigned int> dims_vec )
{
	if (dim>0) {
		std::vector< std::complex<double> > temp;
		for (unsigned int i=0;i<dim;++i) {temp.push_back(0);}
		for (unsigned int j=0;j<dim;++j) {elements.push_back(temp);}
		dims = dims_vec;
	}
}

qmatrix::qmatrix(double p, double phi)
{
	std::vector< std::complex<double> > temp;
	temp.push_back(p);
	temp.push_back(std::sqrt(p*(1-p))*std::polar(1.0,phi));
	elements.push_back(temp);
	temp[0] = std::conj(temp[1]);
	temp[1] = (1-p);
	elements.push_back(temp);
	dims.push_back(2);
}

qmatrix::qmatrix(double p, double q, double phi)
{
	std::vector< std::complex<double> > temp;
	std::complex<double> corner = std::sqrt(p*q)*std::polar(1.0,phi); 
	temp.push_back(p);
	temp.push_back(corner);
	elements.push_back(temp);
	temp[0] = std::conj(corner);
	temp[1] = q;
	elements.push_back(temp);
	dims.push_back(2);
}

qmatrix::qmatrix(std::string sname, unsigned int dim )
{
	if (dim>0) {
		std::vector< std::complex<double> > temp;     // *** should call qmatrix(dim) here
		for (unsigned int i=0;i<dim;++i) {temp.push_back(0);}
		for (unsigned int j=0;j<dim;++j) {elements.push_back(temp);}
		dims.push_back(dim);
		
		if (sname == "up")   {elements[dim-1][dim-1] = 1;}
		if (sname == "down") {elements[0][0] = 1;}
		if (sname == "maxmix") { for (unsigned int i=0;i<dim;++i) {elements[i][i] = 1.0/dim;}}
	}
}

qmatrix::qmatrix(std::string sname, double p, unsigned int dim )
{
	if (dim>0) {
		std::vector< std::complex<double> > temp;     // *** should call qmatrix(dim) here
		for (unsigned int i=0;i<dim;++i) {temp.push_back(0);}
		for (unsigned int j=0;j<dim;++j) {elements.push_back(temp);}
		dims.push_back(dim);
		
		if (sname == "thermal") { 
			double avn = p;  
			double s = avn/(1 + avn);
			elements[0][0] = 1.0;
			double norm = elements[0][0].real();
			for (unsigned int i=1;i<dim;++i) { 
				elements[i][i] = s*elements[i-1][i-1]; 
				norm += elements[i][i].real(); 
			} 
			for (unsigned int i=0;i<dim;++i) { elements[i][i] /= norm; }
		}
	}
}


qmatrix::qmatrix( const qstate& psi ) 
{
	std::vector<unsigned int> psi_dims = psi.get_dims();
	unsigned int size = 1;
	for (unsigned int i=0;i<psi_dims.size();i++) {size *= psi_dims[i];}
    
	std::vector< std::complex<double> > temp;    // *** should call qmatrix(dim) here
	for (unsigned int i=0;i<size;++i) {temp.push_back(0);}
	for (unsigned int j=0;j<size;++j) {elements.push_back(temp);}
	dims = psi_dims;
    
	for (unsigned int i=0;i<size;++i) {
		for (unsigned int j=0;j<size;++j) {
			elements[i][j] = std::conj(psi.get_elem(i))*psi.get_elem(j);
		}
	}
}

qmatrix::qmatrix( const qoperator& op ) 
{
	if (op.matdim()>0) {
		std::vector< std::complex<double> > temp;     // *** should call qmatrix(dim) here
		for (unsigned int i=0;i<op.matdim();++i) {temp.push_back(0);}
		for (unsigned int i=0;i<op.matdim();++i) {elements.push_back(temp);}
		dims.push_back(op.matdim());
		
		for (unsigned int i=0;i<op.size();++i) {
			elements[op.getcol(i)][op.getrow(i)] = op.getval(i);
		}
	}
}

qmatrix::qmatrix( const qoperator& op, std::vector<unsigned int> dims_vec) 
{
	unsigned int totdim = 1; 
	for (unsigned int i=0;i<dims_vec.size();++i) {totdim *= dims_vec[i];}
	if (op.matdim()==totdim) {
		std::vector< std::complex<double> > temp;     // *** should ideally call qmatrix(dim,dims_vec) here
		for (unsigned int i=0;i<op.matdim();++i) {temp.push_back(0);}
		for (unsigned int i=0;i<op.matdim();++i) {elements.push_back(temp);}
		dims = dims_vec;
		
		for (unsigned int i=0;i<op.size();++i) {
			elements[op.getcol(i)][op.getrow(i)] = op.getval(i);
		}
	} 
	else {
		std::cout << "error constructing qmatrix from qoperator: dimension mismatch" << std::endl; 
	}
}

qmatrix::qmatrix( std::complex<double>* mat, int N, std::vector<unsigned int> dims_vec)
{
	std::vector< std::complex<double> > temp(N); 
	for (unsigned int j=0;j<N;j++) { 
		for (unsigned int i=0;i<N;i++) {temp[i] = mat[i + j*N]; } 
		elements.push_back(temp);
	}
	dims = dims_vec; 
}

int qmatrix::clear()
{
	for (unsigned int i=0;i<elements.size();++i) {
		for (unsigned int j=0;j<elements.size();++j) { 
			elements[i][j] = 0;
		}
	}
	return 0;
}

std::complex<double> qmatrix::get_elem( unsigned int row, unsigned int col ) const 
{
	return elements[col][row];
}

int qmatrix::addto_elem( unsigned int row, unsigned int col, std::complex<double> val ) 
{
	elements[col][row] += val;
	return 0;
}     

std::vector<unsigned int> qmatrix::get_dims( ) const 
{
	return dims;
}

int qmatrix::display(std::string name, unsigned int prec ) const
{
	std::cout << name << std::endl;
	for (unsigned int j=0;j<elements.size()/dims[0];++j) {
		for (unsigned int j2=0;j2<dims[0];++j2) {
			
			//------------------------------------- write a line
			std::cout << " ";                                      
			for (unsigned int i=0;i<elements.size()/dims[0];++i) { 
				for (unsigned int i2=0;i2<dims[0];++i2) {
					
					double zr=elements[ i2+dims[0]*i ][ j2+dims[0]*j ].real();
					double zi=elements[ i2+dims[0]*i ][ j2+dims[0]*j ].imag();
					
					if (std::fabs(zr)*std::fabs(zr) < 1e-20) {zr = 0.0;}
					if (zi*zi < 1e-20) {zi = 0.0;}
					
					if ((std::fabs(zr) > std::fabs(zi)*1000) || 
						((std::fabs(zr)==0.0) && (std::fabs(zi)==0.0))) {
						std::cout << std::setw(prec) << std::setprecision(prec) << zr << "  ";
					}
					if (std::fabs(zi) > std::fabs(zr)*1000) {
						std::cout << std::setw(prec) << std::setprecision(prec) << zi << "i ";
					}
					else {
						if (std::fabs(zr) < std::fabs(zi)*1000) {
							std::cout << std::setw(prec) << std::setprecision(prec) 
							<< zr << "+" << zi << "i ";
						}
					}
					
				}
				std::cout << " ";
			}
			std::cout << std::endl;
			
		}
		std::cout << std::endl;
	}
	return 0;
}

int qmatrix::display_abs(std::string name, unsigned int prec ) const
{
	std::cout << name << std::endl;
	for (unsigned int j=0;j<elements.size()/dims[0];++j) {
		for (unsigned int j2=0;j2<dims[0];++j2) {
			
			//------------------------------------- write a line
			std::cout << " ";                                      
			for (unsigned int i=0;i<elements.size()/dims[0];++i) { 
				for (unsigned int i2=0;i2<dims[0];++i2) {
					
					double za=std::abs(elements[ i2+dims[0]*i ][ j2+dims[0]*j ]);
					
					if (za < 1e-5) {za = 0.0;}
					std::cout << std::setw(prec) << std::setprecision(prec) 
					<< za; 
				}
				std::cout << " ";
			}
			std::cout << std::endl;
			
		}
		std::cout << std::endl;
	}
	return 0;
}


std::complex<double> qmatrix::trace() const
{
	std::complex<double> temp(0);
	for (unsigned int i=0;i<elements.size();++i) {temp += elements[i][i];}
	return temp;
}

qmatrix qmatrix::trace_out( int sys ) const 
{
	
	std::vector<unsigned int> new_dims(dims.size()-1); 
	for (int i=0;i<dims.size();i++) { 
		if (i < sys-1) { new_dims[i]=dims[i]; }
		if (i > sys-1) { new_dims[i-1] = dims[i]; }
	}
	
	int new_size(1);
	for (int i=0;i<new_dims.size();i++) { new_size *= new_dims[i]; }
	qmatrix result(new_size); 
	result.dims = new_dims;
	
	int chunk(1); for (int i=0;i<sys-1;i++) { chunk *= dims[i]; } 
	int dim(dims[sys-1]); 
	int num(1);  for (int i=sys;i<dims.size();i++) { num *= dims[i]; } 
	
	for (int i=0;i<dim;i++) {
		
		for (int j1=0;j1<num;j1++) {
			for (int j2=0;j2<num;j2++) {
				
				for (int k1=0;k1<chunk;k1++) {
					for (int k2=0;k2<chunk;k2++) {
						
						result.elements[k2 + j2*chunk][k1 + j1*chunk]  
						+= elements[k2 + i*chunk + j2*chunk*dim][k1 + i*chunk + j1*chunk*dim];
						
					}}
				
			}}
		
	}
	
	return result;
}

qmatrix qmatrix::trace_all_but( int sys ) const 
{
	
	std::vector<unsigned int> new_dims(1); 
	new_dims[0] = dims[sys-1]; 
		
	int new_size(dims[sys-1]);
	
	qmatrix result(new_size); 
	result.dims = new_dims;
		
	int chunk(1); for (int i=0;i<sys-1;i++) { chunk *= dims[i]; }  
	int num(1);  for (int i=sys;i<dims.size();i++) { num *= dims[i]; } 
	
	for (int i=0;i<chunk;i++) {
		for (int j=0;j<num;j++) {
			
			for (int k1=0;k1<new_size;k1++) {
				for (int k2=0;k2<new_size;k2++) {
					
					result.elements[k2][k1] += elements[i + k2*chunk + j*chunk*new_size][i + k1*chunk + j*chunk*new_size]; 
					
				}
			}
			
		}
	}
	
	return result;	
}

std::complex<double> trace( const qmatrix& rho) {
	return rho.trace();
}

int qmatrix::normalise() 
{
	*this = (1.0/this->trace())**this;
	return 0;
}

std::complex<double> qmatrix::inner_prod(const qmatrix& sigma) const
{
	std::complex<double> result = 0;   
	if (elements.size() == sigma.elements.size()) {
		for (unsigned int i=0;i<elements.size();++i) {
			std::complex<double> temp(0);
			for (unsigned int k=0;k<elements.size();++k) {
				temp += elements[k][i]*sigma.elements[i][k];
			}
			result += temp;
		}
	}
	else {
		std::cout << "Cannot take inner product: state size mismatch\n";
		return 0;
	}
	return result;
}

std::complex<double> inner_prod(const qmatrix& rho, const qmatrix& sigma)
{
	return rho.inner_prod(sigma);
}

double qmatrix::fidelity( const qstate& psi) const {
	qoperator op(*this);   // note: not efficient 
	return std::sqrt(op.expect(psi).real()); 
} 

int qmatrix::apply(const qoperator& op)       // A \rho A^\dagger
{
	qmatrix temp(elements.size());
	if (op.matdim()==temp.elements.size()) {
		for (unsigned int i=0;i<op.size();++i) {
			for (unsigned int k=0;k<elements.size();++k) {
				temp.elements[k][op.getrow(i)] += op.getval(i)*elements[k][op.getcol(i)];
			}
		}
		clear();
		for (unsigned int j=0;j<op.size();++j) {
			for (unsigned int i=0;i<elements.size();++i) {
				elements[op.getrow(j)][i] += std::conj(op.getval(i))*temp.elements[op.getcol(j)][i];
			}
		} 
	}
	else {
		std::cout << "Cannot apply operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qmatrix::apply_L(const qoperator& op)       // A \rho 
{
	qmatrix temp = *this;
	clear();
	if (op.matdim()==temp.elements.size()) {
		for (unsigned int i=0;i<op.size();++i) {
			for (unsigned int k=0;k<temp.elements.size();++k) {
				elements[k][op.getrow(i)] += op.getval(i)*temp.elements[k][op.getcol(i)];
			}
		}
	}
	else {
		std::cout << "Cannot apply_L operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qmatrix::apply_R(const qoperator& op)       // \rho A 
{
	qmatrix temp = *this;
	if (op.matdim()==temp.elements.size()) {
		clear();
		for (unsigned int j=0;j<op.size();++j) {
			for (unsigned int i=0;i<elements.size();++i) {
				elements[op.getcol(j)][i] += op.getval(j)*temp.elements[op.getrow(j)][i];
			}
		} 
	}
	else {
		std::cout << "Cannot apply_R operator: size mismatch\n";
		return 1;
	}
	return 0;
}

int qmatrix::apply(const qoperator& op, std::complex<double> cx) 
{
	apply(op);
	*this = cx**this;
	return 0;
}

const qmatrix qmatrix::operator+(const qmatrix& rho) const
{
	qmatrix result(*this);
	unsigned int size = result.elements.size();
	if (size==rho.elements.size()) {
		for (unsigned int i=0;i<size;++i) {
			for (unsigned int j=0;j<size;++j) {
				result.elements[i][j] += rho.elements[i][j];
			}
		}
	}
	else {
		std::cout << "Cannot add states: size mismatch\n";
	}
	return result;
}

const qmatrix qmatrix::operator-(const qmatrix& rho) const
{
	qmatrix result(*this);
	unsigned int size = result.elements.size();
	if (size==rho.elements.size()) {
		for (unsigned int i=0;i<size;++i) {
			for (unsigned int j=0;j<size;++j) {
				result.elements[i][j] -= rho.elements[i][j];
			}
		}
	}
	else {
		std::cout << "Cannot add states: size mismatch\n";
	}
	return result;
}

const qmatrix& qmatrix::operator+=(const qmatrix& rho)
{
	unsigned int size = elements.size();
	if (size==rho.elements.size()) {
		for (unsigned int i=0;i<size;++i) {
			for (unsigned int j=0;j<size;++j) {
				elements[i][j] += rho.elements[i][j];
			}
		}
	}
	else {
		std::cout << "Cannot add-to (+=) the state: size mismatch\n";
	}
	return *this;
}

const qmatrix qmatrix::operator*(const qmatrix& sig) const
{
	qmatrix rho(this->elements.size());
	unsigned int size = this->elements.size();
	unsigned int sigsize = sig.elements.size();
	if (size == sigsize) {
		for (unsigned int i=0;i<size;++i) {
			for (unsigned int k=0;k<size;++k) {
				for (unsigned int j=0;j<size;++j) {
					rho.elements[k][i] +=  this->elements[j][i]*sig.elements[k][j];
				}
			}
		}
	}
	else {     
		std::cout << "Cannot multiply density matricies: size missmatch" << std::endl;
	}
	return rho;
}

const qmatrix qmatrix::operator*(const std::complex<double>& cx) const
{
	qmatrix result = *this;
	unsigned int size = elements.size();
	for (unsigned int i=0;i<size;++i) {
		for (unsigned int j=0;j<size;++j) {
			result.elements[i][j] *= cx;
		}
	}
	return result;
}

const qmatrix qmatrix::operator*(const double& c) const
{
	qmatrix result = *this;
	unsigned int size = elements.size();
	for (unsigned int i=0;i<size;++i) {
		for (unsigned int j=0;j<size;++j) {
			result.elements[i][j] *= c;
		}
	}
	return result;
}

const qmatrix operator*(const std::complex<double>& cx, const qmatrix& rho)
{
	return rho*cx;
}

const qmatrix operator*(const double& c, const qmatrix& rho)
{
	return rho*c;
}

int qmatrix::tensor(const qmatrix& sig)
{
	unsigned int rhosize = this->elements.size(), 
	sigsize = sig.elements.size();
	
	qmatrix result(rhosize*sigsize); 
	
	if ((rhosize>0)&&(sigsize>0)) {
		
		for (unsigned int coli=0;coli<sigsize;coli++) {
			for (unsigned int colj=0;colj<rhosize;colj++) {
				
				for (unsigned int i=0;i<sigsize;i++) {
					for (unsigned int j=0;j<rhosize;j++) {
						result.elements[colj + coli*rhosize][j + i*rhosize] 
						= this->elements[colj][j]*sig.elements[coli][i];
					}
				}
				
			}
		}
		
		result.dims = this->dims;
		for (unsigned int i=0;i<sig.dims.size();++i) {result.dims.push_back(sig.dims[i]);}
		
		*this = result;
	}
	return 0;
}

qmatrix tensor(const qmatrix& rho, const qmatrix& sigma)
{
	qmatrix result = rho; 
	result.tensor(sigma); 
	return result; 
} 

qmatrix tensor(const qmatrix& rho1, const qmatrix& rho2, const qmatrix& rho3)
{
	qmatrix result = tensor(rho1,rho2); 
	result.tensor(rho3); 
	return result; 
} 

qoperator qmatrix::eigvecs_2by2( ) const 
{
    qoperator U;
	U.eigvecs_2by2( *this );
    return U;
}

// ----- to use FFTW uncomment: --------------------------------------

// int qmatrix::fftw_init( fftw_plan* plan1r, fftw_plan& plan1c, 
//                         fftw_plan* plan2r, fftw_plan& plan2c, std::complex<double>* z_ptr )
// {
//    if (dims.size()==2) { 
//      if (dims[1]!=2) { 
//        std::cout << "Error (matrix::fftw): matrix::fftw only works if the second subsystem has dimension 2" << std::endl;
//      }
//    }
//    if (dims.size()>2) {
//        std::cout << "Error (matrix::fftw): matrix::fftw only works for density matrices with 2 subsystems" << std::endl;
//    }
//    else {
//      void *v_ptr; v_ptr = z_ptr;
//      fftw_complex *ptr = reinterpret_cast<fftw_complex*>(v_ptr);
//      
//      int dim1 = dims[0];
//      int dim2;
//      if (dims.size()==2) { dim2 = dims[1]; } else { dim2 = 1; } 
//      int howmany_cols = dim1*dim2*dim2;
//            
//      // collumn transform, forward and back
//      plan1c = fftw_plan_many_dft(1, &dim1, howmany_cols, ptr, &dim1, 1, dim1, ptr, &dim1, 1, dim1, -1, FFTW_MEASURE);
//      plan2c = fftw_plan_many_dft(1, &dim1, howmany_cols, ptr, &dim1, 1, dim1, ptr, &dim1, 1, dim1,  1, FFTW_MEASURE);
//                                                              
//      // row transform, forward and back
//      fftw_iodim row_layout;
//      row_layout.n = dim1;
//      row_layout.is = row_layout.os = dim1*dim2; 
//      
//      fftw_iodim howmany_rows[1];
//      howmany_rows[0].n = dim1*dim2;
//      howmany_rows[0].is = howmany_rows[0].os = 1;
//      //howmany_rows[1].n = dim2;
//      //howmany_rows[1].is = howmany_rows[0].os = dim1*dim1*dim2;
//
//      plan1r[0] = fftw_plan_guru_dft(1, &row_layout, 1, howmany_rows, ptr, ptr,  1, FFTW_MEASURE);
//      plan1r[1] = fftw_plan_guru_dft(1, &row_layout, 1, howmany_rows, ptr + dim1*dim1*dim2, ptr + dim1*dim1*dim2,  1, FFTW_MEASURE);      
//      plan2r[0] = fftw_plan_guru_dft(1, &row_layout, 1, howmany_rows, ptr, ptr, -1, FFTW_MEASURE);
//      plan2r[1] = fftw_plan_guru_dft(1, &row_layout, 1, howmany_rows, ptr + dim1*dim1*dim2, ptr + dim1*dim1*dim2, -1, FFTW_MEASURE);
//    }
//    return 0;
// }
// 
// int qmatrix::fftw_fin( fftw_plan* plan1, fftw_plan& plan2, fftw_plan* plan3, fftw_plan& plan4 )
// {
//    fftw_destroy_plan(plan1[0]);
//    fftw_destroy_plan(plan1[1]);
//    fftw_destroy_plan(plan2);
//    fftw_destroy_plan(plan3[0]);
//    fftw_destroy_plan(plan3[1]);
//    fftw_destroy_plan(plan4);
//    return 0;
// }
// 
// int qmatrix::fftw( fftw_plan* plan1, fftw_plan& plan2 )
// {
//    fftw_execute(plan1[0]);
//    fftw_execute(plan1[1]);
//    fftw_execute(plan2);
//    return 0;
// }
// ----------------------------------------------------------------------

int	qmatrix::to_FORT_cplx_mat( std::complex<double>* mat, int dim ) const 
{
	int dim_tot(1); 
	for (unsigned int j=0;j<dims.size();j++) { dim_tot *= dims[j]; } 
	if (dim == dim_tot) {
		for (unsigned int j=0;j<dim;j++) { 
			for (unsigned int i=0;i<dim;i++) { 
				mat[i + j*dim] = elements[j][i]; 
			} 
		}
	} else { std::cout << "Error: cannot create FORT matrix from qmatrix: dimension mismatch" << std::endl; }
	return 0; 
}

int	qmatrix::get_FORT_cplx_mat( std::complex<double>* mat, int dim, std::vector<unsigned int> new_dims )
{
	int dim_tot(1); 
	for (unsigned int j=0;j<new_dims.size();j++) { dim_tot *= new_dims[j]; } 
	if (dim == dim_tot) {
		for (unsigned int j=0;j<dim;j++) { 
			for (unsigned int i=0;i<dim;i++) { 
				elements[j][i] = mat[i + j*dim]; 
			} 
		}
		dims = new_dims; 
	} else { std::cout << "Error: cannot get FORT matrix: dimension mismatch" << std::endl; }
	return 0; 
}

int qmatrix::copy_to_array_col( std::complex<double> * ptr, unsigned int len)
{
    if (len != elements.size()*elements[0].size()) {
		std::cout << "Error (qmatrix::copy_to_array): length mismatch" << std::endl;
    }
    else {
		unsigned int side = elements.size();
		for (unsigned int j=0;j<side;j++) {
			for (unsigned int k=0;k<side;k++) {
				ptr[k + side*j] = elements[j][k];                    // collumns go fast
			}
		}
    }
    return 0;
}

int qmatrix::copy_from_array_col( std::complex<double> * ptr, unsigned int len)
{
    if (len != elements.size()*elements[0].size()) {
		std::cout << "Error (qmatrix::copy_to_array): length mismatch" << std::endl;
    }
    else {
		unsigned int side = elements.size();
		for (unsigned int j=0;j<side;j++) {
			for (unsigned int k=0;k<side;k++) {
				elements[j][k] = ptr[k + side*j];             // collumns go fast
			}
		}
    }
    return 0;
}

int qmatrix::copy_to_array_row( std::complex<double> * ptr, unsigned int len)
{
    if (len != elements.size()*elements[0].size()) {
		std::cout << "Error (qmatrix::copy_to_array): length mismatch" << std::endl;
    }
    else {
		unsigned int side = elements.size();
		for (unsigned int j=0;j<side;j++) {
			for (unsigned int k=0;k<side;k++) {
				ptr[k + side*j] = elements[k][j];                   // rows go fast
			}
		}
    }
    return 0;
}

int qmatrix::copy_from_array_row( std::complex<double> * ptr, unsigned int len)
{
    if (len != elements.size()*elements[0].size()) {
		std::cout << "Error (qmatrix::copy_to_array): length mismatch" << std::endl;
    }
    else {
		unsigned int side = elements.size();
		for (unsigned int j=0;j<side;j++) {
			for (unsigned int k=0;k<side;k++) {
				elements[k][j] = ptr[k + side*j];                 // rows go fast
			}
		}
    }
    return 0;
}

int qmatrix::output_dist(std::ofstream& outfile, unsigned int subsys) const 
{
	if ((subsys<=dims.size())&&(subsys>0)) {
		int chunksize = 1;
		for (unsigned int i=0;i<(subsys-1);++i) {
			chunksize = chunksize*dims[i];
		}
		int chunknum = 1;
		for (unsigned int i=subsys;i<dims.size();++i) {
			chunknum = chunknum*dims[i];
		}
		int chunksep = chunksize*dims[subsys-1];
		double val;
		for (unsigned int i=0;i<dims[subsys-1];++i) {
			val = 0;
			for (int j=0;j<chunknum;++j) {
				for (int k=0;k<chunksize;++k) {
					val += elements[k + i*chunksize + j*chunksep][k + i*chunksize + j*chunksep].real();
				}
			}
			outfile << val << " ";
		}
		outfile << std::endl;
	}
	else {
		std::cout << dims.size() << std::endl;
		std::cout << "Cannot calculate distribution: subsystem does not exist" << std::endl;
		return 1;
	}
	return 0;
}

int qmatrix::output_1line(std::ofstream& out_real, std::ofstream& out_imag) const 
{
	double r_val(0), i_val(0);
	
	for (int j=0;j<elements.size();++j) {
		for (int k=0;k<elements.size();++k) {
			out_real << elements[k][j].real() << " ";
		}
	}
	out_real << std::endl;
	
	for (int j=0;j<elements.size();++j) {
		for (int k=0;k<elements.size();++k) { 
			out_imag << elements[k][j].imag() << " ";
		}
	}
	out_imag << std::endl;
	
	return 0;
}

int qmatrix::output(std::ofstream& out_real, std::ofstream& out_imag) const 
{
	double r_val(0), i_val(0);
	
	for (int j=0;j<elements.size();++j) {
		for (int k=0;k<elements.size();++k) {
			out_real << elements[k][j].real() << " ";
		}
		out_real << std::endl;
	}
	
	for (int j=0;j<elements.size();++j) {
		for (int k=0;k<elements.size();++k) { 
			out_imag << elements[k][j].imag() << " ";
		}
		out_imag << std::endl;
	}
	return 0;
}

int showstate(const qmatrix& rho, std::string name )
{
	rho.display(name);
	return 0;
}

//  *******************************
//  definitions for class qoperator
//  *******************************

qoperator::qoperator()
{
	dimension = 1;
	row.push_back(0);
	col.push_back(0);
	val.push_back(1);
}

qoperator::qoperator( unsigned int size )
{
	dimension = size;
	for (unsigned int i=0;i<size;++i) {
		row.push_back(i);
		col.push_back(i);
		val.push_back(1);
	}
}

qoperator::qoperator(const qmatrix& rho)
{
	std::vector<unsigned int> dim_vec = rho.get_dims();
	dimension = 1;
	for (unsigned int i=0;i<dim_vec.size();++i) { dimension *= dim_vec[i]; }
	std::complex<double> elem;
	double tol(1e-7/dimension);
	
	for (unsigned int i=0;i<dimension;++i) {
		for (unsigned int j=0;j<dimension;++j) {
			elem = rho.get_elem(i,j);
			if (std::abs(elem) > tol) {
				row.push_back(i);
				col.push_back(j);
				val.push_back(elem);
			}
		}
	}
}

qoperator::qoperator(std::string opname, unsigned int size )
{
	dimension = size;
	if (size>0) {
		if (opname=="identity") {
			for (unsigned int i=0;i<size;++i) {
				row.push_back(i);
				col.push_back(i);
				val.push_back(1);
			}
		}
		if ((opname=="spinX")&&(size==2)) {
			row.push_back(0); col.push_back(1); val.push_back(1.0);
			row.push_back(1); col.push_back(0); val.push_back(1.0);
		}
		if ((opname=="spinY")&&(size==2)) {
			std::complex<double> ei(0,1);
			row.push_back(0); col.push_back(1); val.push_back(-ei);
			row.push_back(1); col.push_back(0); val.push_back(ei);
		}
		if ((opname=="spinZ")&&(size==2)) {
			row.push_back(0); col.push_back(0); val.push_back(1.0);
			row.push_back(1); col.push_back(1); val.push_back(-1.0);
		}
		if ((opname=="proj_up")&&(size==2)) {
			row.push_back(0); col.push_back(0); val.push_back(1.0);
		}
		if ((opname=="proj_down")&&(size==2)) {
			row.push_back(1); col.push_back(1); val.push_back(1.0);
		}		
		if ((opname=="destroy")&&(size>1)) {
			for (unsigned int i=0;i<(size-1);++i) {
				row.push_back(i);
				col.push_back(i+1);
				val.push_back(sqrt(1.0*(i+1)));
			}
		}
		if ((opname=="create")&&(size>1)) {
			for (unsigned int i=0;i<(size-1);++i) {
				row.push_back(i+1);
				col.push_back(i);
				val.push_back(sqrt(1.0*(i+1)));
			}
		}
		if ((opname=="proj_ground")&&(size>1)) {
			row.push_back(0); col.push_back(0); val.push_back(1.0);
		}		
		if ((opname=="Hadamard")&&(size==2)) {
			for (unsigned int i=0;i<size;++i) {
				for (unsigned int j=0;j<size;++j) {
					row.push_back(i);
					col.push_back(j);
					val.push_back(sqrt(0.5));
				}
			}
			val[3] = -val[3];
		}
		if ((opname=="rotate:z->x")&&(size==2)) { 
			row.push_back(0); col.push_back(0); val.push_back(std::sqrt(0.5));
			row.push_back(0); col.push_back(1); val.push_back(std::sqrt(0.5));
			row.push_back(1); col.push_back(0); val.push_back(std::sqrt(0.5));
			row.push_back(1); col.push_back(1); val.push_back(-std::sqrt(0.5));
		}	 
		if ((opname=="rotate:z->y")&&(size==2)) { 
			std::complex<double> ei(0,1);
			row.push_back(0); col.push_back(0); val.push_back(-ei*std::sqrt(0.5));
			row.push_back(0); col.push_back(1); val.push_back( ei*std::sqrt(0.5));
			row.push_back(1); col.push_back(0); val.push_back(std::sqrt(0.5));
			row.push_back(1); col.push_back(1); val.push_back(std::sqrt(0.5));
		}	     
		if ((opname=="rotate:z->x+z")&&(size==2)) {
			row.push_back(0); col.push_back(0); val.push_back(-0.92387953);
			row.push_back(0); col.push_back(1); val.push_back( 0.38268343);
			row.push_back(1); col.push_back(0); val.push_back(-0.38268343);
			row.push_back(1); col.push_back(1); val.push_back(-0.92387953);         
		}	     
		if ((opname=="rotate:z->x-z")&&(size==2)) { 
			row.push_back(0); col.push_back(0); val.push_back( 0.38268343);
			row.push_back(0); col.push_back(1); val.push_back(-0.92387953);
			row.push_back(1); col.push_back(0); val.push_back( 0.92387953);
			row.push_back(1); col.push_back(1); val.push_back( 0.38268343);
		}	 	 
		if ((opname=="rotate:z->x+y")&&(size==2)) { 
			std::complex<double> ei(0,1);
			row.push_back(0); col.push_back(0); val.push_back(0.5 - 0.5*ei);
			row.push_back(0); col.push_back(1); val.push_back(0.5 - 0.5*ei);
			row.push_back(1); col.push_back(0); val.push_back( std::sqrt(0.5));
			row.push_back(1); col.push_back(1); val.push_back(-std::sqrt(0.5));
		}	 	
		if ((opname=="rotate:z->x-y")&&(size==2)) { 
			std::complex<double> ei(0,1);
			row.push_back(0); col.push_back(0); val.push_back(0.5 + 0.5*ei);
			row.push_back(0); col.push_back(1); val.push_back(0.5 + 0.5*ei);
			row.push_back(1); col.push_back(0); val.push_back( std::sqrt(0.5));
			row.push_back(1); col.push_back(1); val.push_back(-std::sqrt(0.5)); 
		}	 		 
		if ((opname=="lineig")&&(size>1)) {
			for (unsigned int i=0;i<size;++i) {
				row.push_back(i);
				col.push_back(i);
				val.push_back(i/(size-1.0) - 0.5);
			}
		}
	}
	if ((size>1)&&(val.size()==0)) {
		std::cout << "Error creating operator. You requested: " << opname << std::endl;
	}
}

qoperator::qoperator(std::string opname, double param, unsigned int size )   // elem has range: 0 ... D^2-1 
{
	std::complex<double> ei(0,1); 
	if (size>0) {
		if (opname=="basis") {
			int elem = static_cast<int>(std::floor(param)); 
			div_t qr = std::div(elem,static_cast<int>(size));
			int l = qr.quot;
			int k = qr.rem;
			
            qoperator X, Z, R;
			X.dimension = size; Z.dimension = size; 
			X.row.clear(); X.col.clear(); X.val.clear();
			Z.row.clear(); Z.col.clear(); Z.val.clear();
			
			for (unsigned int i=0;i<size;++i) { 
				X.row.push_back(static_cast<unsigned int>(std::div(i+k,static_cast<int>(size)).rem));
				X.col.push_back(i);
				X.val.push_back(1);
			}
			for (unsigned int i=0;i<size;++i) { 
				Z.row.push_back(i);
				Z.col.push_back(i);
				Z.val.push_back(std::exp(ei*(6.283185307*i*l)/(1.0*size)));
			}
		    *this = X*Z;
		}
		if ((opname=="rotate_about:y")&&(size==2)) { 
			dimension = size;
			double angle = param; 
			row.push_back(0); col.push_back(0); val.push_back(std::cos(angle/2));
			row.push_back(0); col.push_back(1); val.push_back(-std::sin(angle/2));
			row.push_back(1); col.push_back(0); val.push_back(std::sin(angle/2));
			row.push_back(1); col.push_back(1); val.push_back(std::cos(angle/2));
		}
	}
	if ((size>1)&&(val.size()==0)) {
		std::cout << "Error creating operator. You requested: " << opname << std::endl;
	}
}

qoperator::qoperator(std::string opname, std::vector<double> & x, double beta, unsigned int n )
{
	dimension = x.size();
	if (n>dimension) {
		std::cout << "Warning: new basis cannot be larger than current dimension" << std::endl;
		n = dimension;
	}
	if (((opname=="transform:x->Fock")||(opname=="transform:p->Fock"))&&(dimension>1)&&(n>1)) {
		
		// generate the first n Hermite polynomials 
		std::vector<double> one(x.size(),1.0);
		std::vector< std::vector<double> > HO_eigs(0), HO_imag(0);
		
		HO_eigs.push_back(one);
		std::vector<double> poly(0);
		for (unsigned int i=0;i<dimension;++i) { poly.push_back(2*beta*x[i]); }
		HO_eigs.push_back(poly);  
		
		for (unsigned int m=2;m<n;++m) { 
			for (unsigned int i=0;i<dimension;++i) { 
				poly[i] = 2*beta*x[i]*HO_eigs[m-1][i] - 2*(m-1)*HO_eigs[m-2][i]; 
			}
			HO_eigs.push_back(poly);   
		}
		
		// turn these into the Harmonic Oscillator eigenfunctions
		std::vector<double> base_func(dimension);
		for (unsigned int i=0;i<dimension;++i) { base_func[i] = std::sqrt(beta/std::sqrt(3.1415927))
		*std::exp(-0.5*beta*beta*x[i]*x[i]); }
		std::vector<double> n_func(n);
		n_func[0] = 1;
		for (unsigned int m=1;m<n;++m) { 
			n_func[m] = n_func[m-1]/std::sqrt(2.0*m);
		}
		
		for (unsigned int m=0;m<n;++m) { 
			for (unsigned int i=0;i<dimension;++i) { HO_eigs[m][i] *= n_func[m]*base_func[i]; }  
		}
		
		// output the eigenfunctions to a file (as cols)
		//     std::ofstream oeigs("HO_eigs.dat",std::ios_base::trunc);
		//     for (unsigned int i=0;i<dimension;++i) {
		//       for (unsigned int m=0;m<n;++m) {oeigs << HO_eigs[m][i] << " ";} 
		//       oeigs << std::endl;
		//     }
		//     oeigs.flush(); oeigs.close();
		
		// now use Gramm-Schmidt to make all the HO eigenfunctions exactly orthogonal
		double norm,inprod;
		for (unsigned int m=1;m<n;++m) { 
			
			norm = 0.0;    // normalise orthogonalised vectors (first is trivially orthogonalised)
			for (unsigned int i=0;i<dimension;++i) { norm += HO_eigs[m-1][i]*HO_eigs[m-1][i]; } 
			norm = std::sqrt(norm); 
			for (unsigned int i=0;i<dimension;++i) { HO_eigs[m-1][i] /= norm; }
			
			for (unsigned int l=m;l>0;--l) {    // orthogonalise
				inprod = 0.0; 
				for (unsigned int i=0;i<dimension;++i) { inprod += HO_eigs[m][i]*HO_eigs[l-1][i]; }
				for (unsigned int i=0;i<dimension;++i) { HO_eigs[m][i] -= inprod*HO_eigs[l-1][i]; }
			}  
		}
		
		// normalise last (nth) vector
		norm = 0.0;    // normalise orthogonalised vectors (first is trivially orthogonalised)
		for (unsigned int i=0;i<dimension;++i) { norm += HO_eigs[n-1][i]*HO_eigs[n-1][i]; } 
		norm = std::sqrt(norm); 
		for (unsigned int i=0;i<dimension;++i) { HO_eigs[n-1][i] /= norm; }    
		
		// transform these to the p basis if requested p->Fock 
		qstate psi;  
		std::vector< std::complex<double> > temp(dimension); 
		HO_imag = HO_eigs;
		if (opname=="transform:p->Fock") {
			for (unsigned int m = 0;m<n;++m) {
				for (unsigned int i=0;i<dimension;++i) { temp[i] = HO_eigs[m][i]; };
				psi.assign(temp);
				psi.fft_subsys1(1); 
				for (unsigned int i=0;i<dimension;++i) { 
					HO_eigs[m][i] = psi.get_elem(i).real(); 
					HO_imag[m][i] = psi.get_elem(i).imag();
				}
			}
			// normalise again in case the fft changed the normalisation
			double fftnorm;
			for (unsigned int i=0;i<dimension;++i) { 
				fftnorm += HO_eigs[0][i]*HO_eigs[0][i] + HO_imag[0][i]*HO_imag[0][i]; 
			}
			fftnorm = 1.0/std::sqrt(fftnorm);
			for (unsigned int m = 0;m<n;++m) {
				for (unsigned int i=0;i<dimension;++i) { 
					HO_eigs[m][i] *= fftnorm; 
					HO_imag[m][i] *= fftnorm;
				}
			}
		}
		
		std::complex<double> iu(0,1);
		
		// the HO-eigenfunctions are the rows of the transform operator
		for (unsigned int i=0;i<dimension;++i) {
			for (unsigned int m=0;m<n;++m) {
				row.push_back(m);
				col.push_back(i);
				val.push_back(HO_eigs[m][i] + iu*HO_imag[m][i]);
			}
		}
		
		// also imaginary parts of the HO-eigenfunctions if p->Fock (currently done above)
		
		//     if (opname=="transform:p->Fock") {
		//       for (unsigned int i=0;i<dimension;++i) {
		//         for (unsigned int m=0;m<n;++m) {
		// 	  val[m + i*n] += iu*HO_imag[m][i];
		//         }
		//       }
		//     }
		
	}
	
	if ((n>0)&&(dimension>0)&&(val.size()==0)) {
		std::cout << "Error creating transform. You requested: " << opname << std::endl;
	}
}

qoperator::qoperator( std::string opname, const qstate& psi)
{
	if (opname=="swap") { 
		std::vector<unsigned int> sys_dims = psi.get_dims();
		if ((sys_dims.size() == 2) && (sys_dims[0]==sys_dims[1])) {
			for (unsigned int i=0;i<sys_dims[0];++i) { 
				for (unsigned int j=0;j<sys_dims[0];++j) {  
					row.push_back(j + i*sys_dims[0]);
					col.push_back(i + j*sys_dims[0]);
					val.push_back(1);
				}
			}
		} else {std::cout << "Error in creating swap operator: systems have different size, or more than two systems" << std::endl;}
	}
}

qoperator::qoperator(const std::vector< std::complex<double> >& Fx)
{
	dimension = Fx.size();
	for (unsigned int i=0;i<dimension;++i) {
		row.push_back(i);
		col.push_back(i);
		val.push_back(Fx[i]);
	}
}

qoperator::qoperator(const std::vector<double>& Fx)
{
	dimension = Fx.size();
	for (unsigned int i=0;i<dimension;++i) {
		row.push_back(i);
		col.push_back(i);
		val.push_back(Fx[i]);
	}
}

int qoperator::addto_elem( int i, int j, std::complex<double> z) 
{
    row.push_back(i);
	col.push_back(j);
	val.push_back(z); 
	return 0; 
}

qstate qoperator::apply(const qstate& psi) const
{
	qstate result(psi);
	result.apply(*this);
	return result;
}

qstate qoperator::apply(const qstate& psi, std::complex<double> cx) const
{
	qstate result(psi);
	result.apply(*this,cx);
	return result;
}

std::complex<double> qoperator::expect(const qstate& psi) const
{
	qstate temp(psi);
	temp.apply(*this);
	return inner_prod(psi,temp);
} 

void qoperator::normalize()
{
	double opnorm = 0;
	for (unsigned int i=0;i<val.size();++i) { opnorm += std::norm(val[i]); }
	opnorm = std::sqrt(opnorm);
	for (unsigned int i=0;i<val.size();++i) { val[i] /= opnorm; }
}

void qoperator::Hconj()
{
	std::vector< int > temp(row);
	row = col;
	col = temp;
	for (unsigned int i=0;i<val.size();++i) { val[i] = conj(val[i]); }
}

void qoperator::exp_diag( std::complex<double> t )
{
	std::vector<int> d_empty(dimension);  
	for (unsigned int n=0;n<dimension;n++) { d_empty[n] = n; } 
	for (unsigned int i=0;i<val.size();i++) { 
		if (row[i]==col[i]) { val[i] = std::exp(t*val[i]); d_empty[row[i]] = -1; } 
	}
	for (unsigned int n=0;n<dimension;n++) { 
		if (d_empty[n] >= 0) { 
			row.push_back(n); 
			col.push_back(n); 
			val.push_back(1.0); 
		} 
	} 
}

void qoperator::exp_it_diag( double t ) 
{
	std::complex<double> ei(0,1);
	std::vector<int> d_empty(dimension);  
	for (unsigned int n=0;n<dimension;n++) { d_empty[n] = n; } 
	for (unsigned int i=0;i<val.size();i++) { 
		if (row[i]==col[i]) { val[i] = std::exp(ei*t*val[i]); d_empty[row[i]] = -1; } 
	}
	for (unsigned int n=0;n<dimension;n++) { 
		if (d_empty[n] >= 0) { 
			row.push_back(n); 
			col.push_back(n); 
			val.push_back(1.0); 
		} 
	} 
}

std::complex<double> qoperator::op_inner_prod( const qoperator& qopin) const 
{
	std::complex<double> result = 0; 
	for (unsigned int i=0;i<val.size();i++) { 
		for (unsigned int j=0;j<qopin.val.size();j++) { 
			if ( (qopin.row[j]==row[i]) && ((qopin.col[j]==col[i])) ) { result += qopin.val[j]*val[i]; }
		}
	}
	return result; 
}

int	qoperator::to_FORT_real_mat( double* Mat, int dim ) const  // fortan format: columns are the fast index 
{
	if (dim == dimension) {
		int tot_len = dim*dim; 
		for (unsigned int i=0;i<tot_len;++i) { Mat[i] = 0; }
		for (unsigned int i=0;i<val.size();++i) { Mat[ row[i] + dim*col[i] ] = val[i].real(); }
	} else { std::cout << "Error in to_FORT_real_mat: dimension mismatch" << std::endl; }
	return 0; 
}

int	qoperator::to_FORT_cplx_mat( std::complex<double>* Mat, int dim ) const  // fortan format: columns fast  
{
	if (dim == dimension) {
		int tot_len = dim*dim; 
		for (unsigned int i=0;i<tot_len;++i) { Mat[i] = 0; }
		for (unsigned int i=0;i<val.size();++i) { Mat[ row[i] + dim*col[i] ] = val[i]; }
	} else { std::cout << "Error in to_FORT_cplx_mat: dimension mismatch" << std::endl; }
	return 0; 
}

const qoperator qoperator::operator+(const qoperator& op2) const
{
	qoperator op1 = *this;
	qoperator op = op2;
	if (op1.dimension==op2.dimension) {
		op.dimension = op1.dimension;
	}
	else {
		std::cout << "cannot add operators: different dimensions\n";
	}
	bool exists;
	for (unsigned int i=0;i<op1.size();++i) {
		exists = false;
		for (unsigned int j=0;j<op2.size();++j) {
			if ((op2.row[j] == op1.row[i])&&(op2.col[j] == op1.col[i])) {
				op.val[j] += op1.val[i];
				exists = true;
			}
		}
		if (exists==false) {
			op.row.push_back(op1.row[i]);
			op.col.push_back(op1.col[i]);
			op.val.push_back(op1.val[i]);
		}
	}
	return op;
}

const qoperator qoperator::operator-(const qoperator& op1) const
{
	qoperator op2 = *this;
	qoperator op = op2;
	if (op1.dimension==op2.dimension) {
		op.dimension = op1.dimension;
	}
	else {
		std::cout << "cannot add operators: different dimensions\n";
	}
	bool exists;
	for (unsigned int i=0;i<op1.size();++i) {
		exists = false;
		for (unsigned int j=0;j<op2.size();++j) {
			if ((op2.row[j] == op1.row[i])&&(op2.col[j] == op1.col[i])) {
				op.val[j] -= op1.val[i];
				exists = true;
			}
		}
		if (exists==false) {
			op.row.push_back(op1.row[i]);
			op.col.push_back(op1.col[i]);
			op.val.push_back(-op1.val[i]);
		}
	}
	return op;
}

const qoperator qoperator::operator*(const qoperator& op2) const
{
	qoperator op1 = *this;
	qoperator op("identity",0);
	if (op1.dimension==op2.dimension) {
		op.dimension = op1.dimension;
	}
	else {
		std::cout << "cannot multiply operators: different dimensions\n";
	}
	bool exists;
	for (unsigned int i=0;i<op1.size();++i) {
		for (unsigned int j=0;j<op2.size();++j) {
			if (op1.col[i] == op2.row[j]) {
				exists = false;
				for (unsigned int k=0;k<op.size();++k) {
					if ((op.row[k]==op1.row[i])&&(op.col[k]==op2.col[j])) {
						op.val[k] += op1.val[i]*op2.val[j];
						exists = true;
					}
				}
				if (exists==false) {
					op.row.push_back(op1.row[i]);
					op.col.push_back(op2.col[j]);
					op.val.push_back(op1.val[i]*op2.val[j]);
				}
			}
		}
	}
	return op;
}

const qoperator qoperator::operator*(const std::complex<double>& cx) const
{
	qoperator op = *this;
	for (unsigned int i=0;i<op.val.size();++i) {op.val[i] *= cx;}
	return op;
}

const qoperator qoperator::operator*(const double& c) const
{
	qoperator op = *this;
	for (unsigned int i=0;i<op.val.size();++i) {op.val[i] *= c;}
	return op;
}

const qoperator operator*(const std::complex<double>& cx, const qoperator& op)
{
	return op*cx;
}

const qoperator operator*(const double& c, const qoperator& op)
{
	return op*c;
}

const qoperator& qoperator::operator+=(const qoperator& op1) 
{
	qoperator op2 = *this;
	if (op1.dimension!=op2.dimension) {
		std::cout << "cannot add operators: different dimensions\n";
	}
	bool exists;
	for (unsigned int i=0;i<op1.size();++i) {
		exists = false;
		for (unsigned int j=0;j<op2.size();++j) {
			if ((op2.row[j] == op1.row[i])&&(op2.col[j] == op1.col[i])) {
				(*this).val[j] += op1.val[i];
				exists = true;
			}
		}
		if (exists==false) {
			(*this).row.push_back(op1.row[i]);
			(*this).col.push_back(op1.col[i]);
			(*this).val.push_back(op1.val[i]);
		}
	}
	return *this;
}

const qoperator& qoperator::operator*=(const std::complex<double>& z) 
{
	for (unsigned int i=0;i<this->val.size();i++) {
		this->val[i] *= z; 
	}
	return *this;
}

int qoperator::tensor(const qoperator& op2)
{
    qoperator op1 = *this;
    qoperator op(0);
    op.dimension = op1.dimension*op2.dimension;
    for (unsigned int j=0;j<op2.size();++j) {
		for (unsigned int i=0;i<op1.size();++i) {
			op.row.push_back(op1.row[i] + op1.dimension*op2.row[j]);
			op.col.push_back(op1.col[i] + op1.dimension*op2.col[j]);
			op.val.push_back(op1.val[i]*op2.val[j]);
		}
    }
    *this = op;
    return 0;
}

qoperator tensor(const qoperator& op1, const qoperator& op2)
{
	qoperator result = op1;
	result.tensor(op2);
	return result; 
}

qoperator tensor(const qoperator& op1, const qoperator& op2, const qoperator& op3)
{
	qoperator result = tensor(op1,op2);
	result.tensor(op3); 
	return result;
}

qoperator Hconj(const qoperator& op) { 
	qoperator result = op;
	result.Hconj();
	return result;
}

void qoperator::eigvecs_2by2( const qmatrix& rho ) 
{
    std::complex<double> a,b,c,d;
	a = rho.get_elem(0,0);
	b = rho.get_elem(0,1);
	c = rho.get_elem(1,0);
	d = rho.get_elem(1,1);
	
	std::complex<double> lam1 = (a+d)/2.0 + std::sqrt( (a+d)*(a+d)/4.0 + b*c - a*d );
	std::complex<double> lam2 = (a+d)/2.0 - std::sqrt( (a+d)*(a+d)/4.0 + b*c - a*d );
	
	std::complex<double> y1 = std::sqrt( (lam1-a)*std::conj(lam1-a) / ( (lam1-a)*std::conj(lam1-a) + b*std::conj(b)  )   );
	std::complex<double> x1 = y1*b/(lam1-a);
	
	std::complex<double> y2 = std::sqrt( x1*std::conj(x1) / ( x1*std::conj(x1) + y1*std::conj(y1) ) );
	std::complex<double> x2 = -y2*std::conj(y1)/std::conj(x1);
	
	val[0] = x1;
	row.push_back(1); col.push_back(0); val.push_back( y1 );
	row.push_back(0); col.push_back(1); val.push_back( x2 );
	row.push_back(1); col.push_back(1); val.push_back( y2 );         
	
	dimension = 2; 
}

int qoperator::display( std::string name, unsigned int prec ) const
{
	qmatrix rho(*this);
	rho.display(name, prec);
	return 0;
}

int qoperator::displaylist() const
{
	for (unsigned int i=0;i<val.size();++i) {
		std::cout << row[i] << " " << col[i] << ": " << val[i] << std::endl;
	}
	std::cout << std::endl;
	return 0;
}

// **********************
// inter-object functions
// **********************

const qstate operator*(const qoperator& op, const qstate& psi) {
	qstate temp = psi;
	temp.apply(op);
	return temp;
}

const qmatrix operator*(const qoperator& op, const qmatrix& rho) {
	qmatrix temp = rho;
	temp.apply_L(op);
	return temp;
}

const qmatrix operator*(const qmatrix& rho, const qoperator& op) {
	qmatrix temp = rho;
	temp.apply_R(op);
	return temp;
}

int RK4_stepper(const qoperator& dU, qstate& psi) {
	static qstate K1(psi), K2(psi), K3(psi), K4(psi);
   	K1 = dU*psi;
	K2 = dU*(psi + 0.5*K1); 
	K3 = dU*(psi + 0.5*K2); 
	K4 = dU*(psi + K3);
	psi += 0.166666*(K1 + 2.0*K2 + 2.0*K3 + K4);
	psi.normalise(); 
	return 0; 
}

int RK4_stepper(const qoperator& dU, qmatrix& rho) {
	static qmatrix K1(rho), K2(rho), K3(rho), K4(rho), rho_temp(rho);
	K1 = dU*rho - rho*dU;
	rho_temp = rho + 0.5*K1; 
	K2 = dU*rho_temp - rho_temp*dU;
	rho_temp = rho + 0.5*K2;
	K3 = dU*rho_temp - rho_temp*dU;
	rho_temp = rho + K3;
	K4 = dU*rho_temp - rho_temp*dU;
	rho += 0.166666*(K1 + 2.0*K2 + 2.0*K3 + K4);
	rho.normalise(); 
	return 0; 
}

int RK4_stepper( qmatrix (*calc_drho) (qmatrix&), qmatrix& rho) {
	static qmatrix K1(rho), K2(rho), K3(rho), K4(rho), rho_temp(rho); 
	K1 = (*calc_drho) (rho); 
	rho_temp = rho + 0.5*K1; 
	K2 = (*calc_drho) (rho_temp); 
	rho_temp = rho + 0.5*K2; 
	K3 = (*calc_drho) (rho_temp); 
	rho_temp = rho + K3; 
	K4 = (*calc_drho) (rho_temp); 
	rho += 0.166666*(K1 + 2.0*K2 + 2.0*K3 + K4); 
	rho.normalise(); 
	return 0; 
}

//  *******************
//  extra functions 
//  *******************

int make_noise( std::vector<double>& noise)
{
	int ns = noise.size()/2;
	double tn1,tn2,rn;
	for (unsigned int i=0;i<ns;++i) {
		rn = 2;
		while (rn >= 1 || rn == 0) {
            tn1 = 2*static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1;
            tn2 = 2*static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1;
            rn = tn1*tn1 + tn2*tn2;
		}
		rn = std::sqrt((-2*std::log(rn)/rn));
		noise[2*i] = tn1*rn;
		noise[2*i+1] = tn2*rn;
	}
	return 0;
}

double make_gauss( )
{
	double tn1,tn2,rn;
	rn = 2;
	while (rn >= 1 || rn == 0) {
		tn1 = 2*static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1;
		tn2 = 2*static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1;
		rn = tn1*tn1 + tn2*tn2;
	}
	rn = std::sqrt((-2*std::log(rn)/rn));
	return tn1*rn;
}

double make_cauchy( )
{
	double r1;
	do { 
		r1 = static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 0.5; 
	} 
	while ((r1 == 0.5) || (r1 == -0.5));
	r1 = 3.14159*r1;   // number in (-pi/2 and pi/2)
	return std::tan(r1);
}

double make_sym_stable( double dt, double alpha, double sigma)
{
	double r,W;
	
	if (alpha>0.0) { 
		
		if (alpha == 1) { 
			
			return dt*sigma*make_cauchy();
			
		} else {
			
			double s = 1/alpha;
			
			do { r = static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 0.5; } 
			while ((r == 0.5) || (r == -0.5));
			r = 3.14159*r;
			
			do { W = static_cast<double>(rand())/static_cast<double>(RAND_MAX); }
			while ((W == 0) || (W == 1));
			W = -std::log(W); 
			
			double k = std::cos((1-alpha)*r)/W;
			
			return sigma*std::pow(dt,s)*std::pow(k,(1-alpha)*s)*std::sin(alpha*r)/std::pow(std::cos(r),s); 
			
		}
		
	} else {
		std::cout << "error: Levy parameter alpha must be greater than zero." << std::endl;
		return 0;
	}
	
} 

double make_poisson( double dt, double lambda, double jump )
{
	double W = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
	
	if ( W < (1 - std::exp(-lambda*dt)) ) {return jump;} 
	else {return 0.0;}
}

double make_stable_sub_poisson( double dt, double alpha, double sigma, double lambda , double jump)
{
	double dN = make_poisson(dt,lambda,jump);
	return make_sym_stable(dN,alpha,sigma); 
} 

double make_stable_sub_wiener( double dt, double alpha, double sigma, double g_width )
{
	double dW = make_gauss()*std::sqrt(dt)*g_width;
	return make_sym_stable(dW,alpha,sigma); 
} 

// ---------- requires fftw: ---------------------------------------------------------------
// int stable_density( double alpha, double sigma, double dx, std::vector<double>& density, 
//                                                            std::vector<double>& xgrid )  // n should be even
// {
//   double pi = 3.1415927;
//   int nx = density.size();  
//   double Lx = dx*nx, 
//          ds = 2*pi/Lx,
//		  Ls = nx*ds;  
//		  
//   std::vector<double> sgrid(nx);
//   std::vector< std::complex<double> > chi(nx); 
//   
//   void* itt = &chi[0];
//   fftw_plan plan;
//   
//   plan = fftw_plan_many_dft(1, &nx, 1, reinterpret_cast<fftw_complex*>(itt), &nx, 1, nx, 
//									    reinterpret_cast<fftw_complex*>(itt), &nx, 1, nx, -1, FFTW_MEASURE);
//										
//   for (int i=0;i<nx/2;i++) {sgrid[i] = i*ds;} 
//   for (int i=0;i<nx/2;i++) {sgrid[i+nx/2] = i*ds - Ls/2;} 
//   
//   //for (int i=0;i<nx;i++) {sgrid[i] = i*ds - Ls/2;} 
//   for (int i=0;i<nx;i++) {chi[i] = std::exp(-std::pow(sigma*sigma*sgrid[i]*sgrid[i],alpha/2));} 
//
//   fftw_execute(plan);
//   fftw_destroy_plan(plan);
//   
//   for (int i=0;i<nx/2;i++) {density[i] = chi[i+nx/2].real();}   // center is at i = nx/2 
//   for (int i=0;i<nx/2;i++) {density[i+nx/2] = chi[i].real();}
//   for (int i=0;i<nx;i++) {xgrid[i] = i*dx - Lx/2;}
// 
//   double norm(0);
//   for (int i=0;i<nx;i++) {norm += density[i]*dx;}
//   for (int i=0;i<nx;i++) {density[i] /= norm;}
//   
//   return 0; 
// }
// -------------------------------------------------------------------------

int sample( const std::vector<double>& dist )  
{
	double z = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	for (int i=0;i<dist.size();i++) { if ( z < dist[i]) {return i;} } 
	return dist.size()-1; 
}

//  *******************
//  auxillary functions 
//  *******************

// the following requires GNU's gsl
// int fft(std::vector< std::complex<double> >& x, gsl_fft_direction sign)
// {
//  std::vector< std::complex<double> >::iterator itt = x.begin();
//  return gsl_fft_complex_radix2_transform(reinterpret_cast<gsl_complex_packed_array>(itt),
//                                          1,x.size(),sign);
// }

