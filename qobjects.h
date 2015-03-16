
//  This is the header file containing the classes for quantum states
//  and quantum operators. These include the member functions that 
//  allow the operators to act on the states, and the functions that give 
//  binary operations.
//
//  The classes are:      qstate
//                        qmatrix
//                        qoperator
//
//  The operators are:    qoperator: tensor,*,*(scalar),+,-
//                        qstate:    tensor,*(scalar),+,-,>>,<<
//                        qmatrix:   tensor,*,*(scalar),+,-
//
//                        general: qoperator*qstate, 
//                                 qoperator*qmatrix, 
//                                 qmatrix*qoperator
//

#ifndef QSTATE_H
#define QSTATE_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 
#include <vector>
#include <complex>
#include <string>

// ----to use LAPACK uncomment: ---------------------------------
#include "lapack_interface.h"
// --------------------------------------------------------------

// ----to use FFTW uncomment: -----------------------------------
// #include "fftw3.h"
// --------------------------------------------------------------

 class qoperator;
 class qmatrix;

 class qstate {
       std::vector< std::complex<double> > qcoeff;
       std::vector<unsigned int> dims;
       
       int writeblock( int ) const;
 public:
       qstate();
       qstate( std::vector< std::complex<double> >& );
       qstate( unsigned int dim );
       qstate( std::string, unsigned int dim = 2 );
       qstate( std::string, std::complex<double>, unsigned int );
       qstate( std::string, std::complex<double>,  std::complex<double>, unsigned int );
       qstate( std::string, unsigned int, const double&, const double&, const double&, const double&, const double&, const double& );

       int                       clear();
       int                       normalise();
       int                       assign( std::vector< std::complex<double> >& );
       int                       assign( std::complex<double>* start, std::complex<double>* finish);
       std::complex<double>      get_elem( unsigned int ) const;
       int                       addto_elem( unsigned int, std::complex<double> ); 
       void*                     get_array();
	   unsigned int				 len() const;
       std::vector<unsigned int> get_dims() const; 
       int                       shift( const int n ); 
       int                       rotate( const int n ); 
	   int						 fftshift(); 

       int display1D( std::string name = "qstate: " ) const; 
       int display( std::string name = "qstate: " ) const;

       double               norm() const;
       qmatrix              trace_out( int sys ) const;  // trace out a single system
	   qmatrix              trace_all_but( int sys ) const;  // trace over all but one system 
       std::complex<double> inner_prod( const qstate& ) const;
       double				parity() const; 
	   int                  apply( const qoperator& );
       int                  apply( const qoperator&, std::complex<double> );
       int                  apply_Hconj( const qoperator& );
	   int					apply_diag( std::complex<double>*, int dim ); 
       int                  tensor( const qstate& );
	 
	   int apply_FORT_cplx_mat( std::complex<double>* mat, int dim );  // requires lapack 
	   int applyT_FORT_cplx_mat( std::complex<double>* mat, int dim ); // requires lapack
       
       // ------ to use fftw uncomment ---------
       // int fftw_init( fftw_plan&, fftw_plan& );   // requires library fftw3
       // int fftw_fin( fftw_plan&, fftw_plan& );    // requires library fftw3
       // int fftw( fftw_plan& );                    // requires library fftw3
       // --------------------------------------
         
       int fft_init( double** );           // requires an fft library
       int fft_ss1b( int sign, double* );  // requires an fft library
       int fft_subsys1(int sign);          // requires an fft library
       
       // these assume the system contains only qubits
       int cnot( const int control, const int target ); // cnot operation
       int flip( const int target ); // not operation
       
       const qstate operator+( const qstate& ) const;
       const qstate operator-( const qstate& ) const;
       const qstate operator*( const std::complex<double>& ) const;
       const qstate operator*( const double& ) const;
	   const qstate& operator=( const qstate& );
       const qstate& operator+=( const qstate& );
	   const qstate& operator*=( std::complex<double> );
	   const qstate& operator*=( double );

       std::ostream& output( std::ostream& ) const; 
       std::istream& input( std::istream& );
       
       int output(std::ofstream& out_real, std::ofstream& out_imag) const;
       int output_dist( std::ofstream&, unsigned int );
 };
 
 // qstate associated functions
 
 qstate               tensor(const qstate&, const qstate&);
 qstate               tensor(const qstate&, const qstate&, const qstate&);
 const qstate         operator*(const std::complex<double>&, const qstate&);
 const qstate         operator*(const double&, const qstate&);
 std::complex<double> inner_prod(const qstate&, const qstate&);
 int                  showstate(const qstate&, std::string name = "qstate: ");
 std::ostream&        operator<<(std::ostream&, const qstate& );
 std::istream&        operator>>(std::istream&, qstate& );

 
 class qmatrix {
       std::vector< std::vector< std::complex<double> > > elements;    // format: elements is a vector of the columns of rho
       std::vector<unsigned int> dims; 
 
 public:
       qmatrix( unsigned int dim = 2 );
       qmatrix( unsigned int dim, std::vector<unsigned int> dims_vec );
       qmatrix( std::string, unsigned int dim = 2 );
	   qmatrix( std::string, double, unsigned int dim );
       qmatrix( double, double );
       qmatrix( double, double, double );
       qmatrix( const qstate& ); 
       qmatrix( const qoperator& );
	   qmatrix( const qoperator&, std::vector<unsigned int> ); 
	   qmatrix( std::complex<double>*, int, std::vector<unsigned int> );   // create from FORTRAN matrix 
       
       int                       clear();
       int                       normalise();
       std::complex<double>      get_elem( unsigned int row, unsigned int col ) const; 
       int                       addto_elem( unsigned int row, unsigned int col, std::complex<double> val ); 
       std::vector<unsigned int> get_dims() const;
	   

       int display( std::string name = "qmatrix: ", unsigned int prec = 2 ) const;
	   int display_abs( std::string name = "qmatrix: ", unsigned int prec = 2 ) const;
       
       std::complex<double> trace() const;
       qmatrix              trace_out( int sys ) const;   // trace_out a single system 
	   qmatrix              trace_all_but( int sys ) const;  // trace over all but one system 
       std::complex<double> inner_prod( const qmatrix& ) const;
       double               fidelity( const qstate& ) const;
       int                  apply( const qoperator& );
       int                  apply_L( const qoperator& );
       int                  apply_R( const qoperator& );
       int                  apply( const qoperator&, std::complex<double> );
       int                  tensor( const qmatrix& );
	   qoperator            eigvecs_2by2() const; 

       // ------- to use fftw uncomment --------------------------------
       // int fftw_init( fftw_plan*, fftw_plan&, 
       //                fftw_plan*, fftw_plan&, std::complex<double> * );   // requires library fftw3
       // int fftw_fin( fftw_plan*, fftw_plan&, fftw_plan*, fftw_plan& );    // requires library fftw3
       // int fftw( fftw_plan* , fftw_plan& );                               // requires library fftw3 
       // --------------------------------------------------------------
         
	   int to_FORT_cplx_mat( std::complex<double>* mat, int dim ) const; 
	   int get_FORT_cplx_mat( std::complex<double>* mat, int dim, std::vector<unsigned int> new_dims ); 

       int copy_to_array_col( std::complex<double> *, unsigned int );  
       int copy_from_array_col( std::complex<double> *, unsigned int );
       int copy_to_array_row( std::complex<double> *, unsigned int );  
       int copy_from_array_row( std::complex<double> *, unsigned int );
       
       const qmatrix operator+( const qmatrix& ) const;
       const qmatrix operator-( const qmatrix& ) const;
       const qmatrix operator*( const qmatrix& ) const;
       const qmatrix operator*( const std::complex<double>& ) const;
       const qmatrix operator*( const double& ) const;
       const qmatrix& operator+=( const qmatrix& );

       int output_1line(std::ofstream& out_real, std::ofstream& out_imag) const; 
       int output(std::ofstream& out_real, std::ofstream& out_imag) const;
       int output_dist( std::ofstream&, unsigned int ) const;
 };
 
 // qmatrix associated functions

 std::complex<double> trace(const qmatrix&);
 qmatrix              tensor(const qmatrix&, const qmatrix&);
 qmatrix              tensor(const qmatrix&, const qmatrix&, const qmatrix&);
 const qmatrix        operator*(const std::complex<double>&, const qmatrix&);
 const qmatrix        operator*(const double&, const qmatrix&);
 int                  showstate(const qmatrix&, std::string name = "qmatrix: " );
 
 
 class qoperator {
       std::vector< int > row;
       std::vector< int > col;
       std::vector< std::complex<double> > val;
       unsigned int dimension;
 public:
       qoperator();
       qoperator( unsigned int size );
       qoperator( const qmatrix& );
       qoperator( std::string opname, unsigned int size = 2 );
	   qoperator( std::string opname, double, unsigned int );
       qoperator( std::string, std::vector<double> &, double, unsigned int );
       qoperator( std::string, const qstate& );

       qoperator( const std::vector<double>& );
       qoperator( const std::vector< std::complex<double> >& );
       
       int                  getrow( int i) const { return row[i]; }
       int                  getcol( int i) const { return col[i]; }
       std::complex<double> getval( int i) const { return val[i]; }
	   int					addto_elem( int, int, std::complex<double> );
       unsigned int         matdim() const { return dimension; }
       unsigned int         size() const { return val.size(); }  // size and dim are different
       
       int tensor( const qoperator& );
       
       void                 normalize();
       void                 Hconj();
       std::complex<double> expect( const qstate& ) const; 
       void                 eigvecs_2by2( const qmatrix& ); 
       void                 exp_diag( std::complex<double> );   // formerly exp_t_elem()
       void                 exp_it_diag( double );              // formerly exp_it_elem() 
       std::complex<double> op_inner_prod( const qoperator& ) const; 
       qstate               apply( const qstate& ) const;
       qstate               apply( const qstate&, std::complex<double> ) const;
	 
	   int					to_FORT_real_mat( double* mat, int dim ) const; 
	   int					to_FORT_cplx_mat( std::complex<double>* mat, int dim ) const; 
	   
       const qoperator operator+( const qoperator& ) const;
       const qoperator operator-( const qoperator& ) const;
       const qoperator operator*( const qoperator& ) const;
       const qoperator operator*( const std::complex<double>& ) const;
       const qoperator operator*( const double& ) const;
       const qoperator& operator+=( const qoperator& );
	   const qoperator& operator*=( const std::complex<double>& );

       int display( std::string name = "qoperator: ", unsigned int prec = 2 ) const;
       int displaylist() const; 
 };
 
 // qoperator associated functions
 
 qoperator       tensor(const qoperator&, const qoperator&); 
 qoperator       tensor(const qoperator&, const qoperator&, const qoperator&);
 qoperator       Hconj(const qoperator&);
 const qoperator operator*(const std::complex<double>&, const qoperator&);
 const qoperator operator*(const double&, const qoperator&);
 
 // inter-object associated functions
 
 const qstate  operator*(const qoperator&, const qstate&);
 const qmatrix operator*(const qoperator&, const qmatrix&);
 const qmatrix operator*(const qmatrix&, const qoperator&);
 
 int RK4_stepper(const qoperator& dU, qstate& psi); 
 int RK4_stepper(const qoperator& dU, qmatrix& rho); 
 int RK4_stepper(qmatrix (*calc_drho) (qmatrix&), qmatrix& rho); 
 
 // functions for generating random variables 
 
 int make_noise( std::vector<double>& );    // note: size of vector should be even
 
 double make_cauchy();
 double make_gauss();
 double make_sym_stable( double dt, double alpha, double sigma );
 double make_poisson( double dt, double lambda, double jump);
 double make_stable_sub_poisson( double dt, double alpha, double sigma, double lambda, double jump );
 double make_stable_sub_wiener( double dt, double alpha, double sigma, double g_width );
 
 // requires fftw: 
 //int stable_density( double alpha, double sigma, double dx, std::vector<double>& density, std::vector<double>& xgrid ); 
 
 int sample( const std::vector<double>& x_distribution );

#endif
 
 
