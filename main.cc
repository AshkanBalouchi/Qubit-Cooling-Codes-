
#include "main.h"
#include "qobjects.h"
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
using namespace std;



double sq2( double x) { return x*x; }

double my_f (const gsl_vector *x, void *params); 
void  my_df (const gsl_vector *x, void *params, gsl_vector *df); 
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df); 

// ------------- global values ----------------
int	Nspace = 2; 		// # of parameters in the search space 
int	calls = 0; 
// --------------------------------------------

int main(int argc, char* argv[]){ 
	
	const double pi = 3.14159; 
	
	const int seed = 123456;    // random seed
	std::srand(seed); 
	
	int max_iter = 1000; 
						
	// ----------------- optimization log output file
	std::ofstream optlogfl("opt_calc.log",std::ios_base::trunc);
    std::ofstream optdatfl("opt_calc.dat",std::ios_base::trunc);
	std::ofstream optpulfl("opt_protocol.dat",std::ios_base::trunc);

	// ----------------- initialize timing 
	time_t time0, *time_ptr = 0;   
	time0 = time(time_ptr);

	// ----------------- optimiser set-up 
	int status;
	double val; 
		
	gsl_multimin_function_fdf my_func;
	double par[1] = { 0 }; // function parameters if desired 
	my_func.n = static_cast<size_t>(Nspace);
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = par;
	
	double linemin_tol = 0.01; 
	double step_size = 0.1; 
	double grad_bound = 1e-4; 
		
	unsigned int iter = 0;
	
	// initial point 
	gsl_vector *x;
	x = gsl_vector_alloc(static_cast<size_t>(Nspace));
	for (int i=0;i<Nspace;i++) { gsl_vector_set(x, i, 0.5); }
	
	// initialize optimizer 
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc(T, static_cast<size_t>(Nspace));
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, linemin_tol);
							
	// -------- optimize 
	do {
		
		status = gsl_multimin_fdfminimizer_iterate(s); 
		if (status) break;  // minimum has been found during set-up, so end the loop
		
		status = gsl_multimin_test_gradient(s->gradient, grad_bound);
		 
		val = s->f;  // put the current value into val 
		
		// output optimization progress to log file 
		if (status == GSL_SUCCESS) {optlogfl << "Minimum found at: " << std::endl;}
		optlogfl << "iter: " << iter << " val: " << val << std::endl << std::flush; 
		optdatfl << iter << " " << val << std::endl << std::flush;
		if (std::div(iter,500).rem==0) {
			optlogfl << "iter: " << iter << " val: " << val << std::endl << std::flush; 
			optlogfl << "ctrl_pars: "; 
			for (int i=0;i<Nspace;i++) { optlogfl << gsl_vector_get(s->x, i) << " "; } 
			optlogfl << std::endl << std::flush; 
		}
		
	} while (status == GSL_CONTINUE && iter++ < max_iter);
	
	// output final value to log file 
	if (status != GSL_SUCCESS) { 
		optlogfl << "Minimum not yet found: " << std::endl; 
		optlogfl << "final iter: " << iter << " val: " << val << std::endl << std::flush; 
	}  
	
	double fin_x[Nspace];  
	for (int i=0;i<Nspace;i++) { fin_x[i] = gsl_vector_get(s->x, i); }
				
	// tidy up minimizer 
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
		
	// -------------------- outputting results to standard out 
	if (status == GSL_SUCCESS) { std::cout << "Minimum found: " << std::endl; } 
	else { std::cout << "Minimum not yet found: " << std::endl; }
	std::cout << "val = " << val << std::endl; 
	std::cout << "ctrl_pars: "; 
	for (int i=0;i<Nspace;i++) { std::cout << fin_x[i] << " "; } 
	std::cout << std::endl; 
	std::cout << "total function calls: " << calls << std::endl; 
	optpulfl << std::endl;
	
	// output simulation time
	double took = difftime(time(time_ptr),time0); 
	optlogfl << std::endl << "Simulation took " << took    << " secs";
	optlogfl << std::endl << "             or " << took/60 << " mins"; 
	optlogfl << std::endl << "             or " << took/(60*60) << " hours" << std::endl;
//	}

	optpulfl.close(); 
	optlogfl.close();
	optdatfl.close();
//	MPI_Finalize();
					 
	return 0;
}
