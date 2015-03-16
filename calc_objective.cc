
#include "main.h"
#include "lapack_interface.h" 
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



// -------------- global values (defined in main.cc) --------------- 
extern int Nspace;	// # of parameters in the search space
extern int calls; 

// simulation routine (defined below)   
double sim(double* x);   

// function to minimize
double my_f(const gsl_vector* x_vec, void* params) { 
	double val, x[Nspace];  	
    for (int i=0;i<Nspace;i++) { x[i] = gsl_vector_get(x_vec, i); }
	return sim(x); 
}

// first-order finite-diff gradient 
void  my_df(const gsl_vector* x_vec, void* params, gsl_vector* df) { 
	static double h = 1e-6; 
	double x[Nspace], y[Nspace]; 
    for (int i=0;i<Nspace;i++) { y[i] = x[i] = gsl_vector_get(x_vec, i); }
	double val = sim(x); 
	for (int j=0;j<Nspace;j++) {
		y[j] += h;  
		gsl_vector_set(df, j, (sim(y) - val)/h );  // first-order approx. of the grad 
		y[j] = x[j]; 
	} 
	cout << "my_df:  "<<gsl_vector_get(df, 0)<<endl;
}

// first-order finite-diff gradient 
void my_fdf(const gsl_vector *x_vec, void *params, double *f, gsl_vector *df) {  
	static double h = 1e-6; 
	double x[Nspace], y[Nspace]; 
    for (int i=0;i<Nspace;i++) { y[i] = x[i] = gsl_vector_get(x_vec, i); }
	double val = sim(x); 
	for (int j=0;j<Nspace;j++) {
		y[j] += h;  
		gsl_vector_set(df, j, (sim(y) - val)/h );  // first-order approx. of the grad 
		y[j] = x[j]; 
	} 
	*f = val;
    cout << "my_fdf:  "<<gsl_vector_get(df, 0)<<endl;
	
}

// --------------- Simulation --------------------
/*
 double coss(double x){
	return cos(x);
}


double sim(double* x) {  
	
	static bool first_call = true;  // in case you want to do something only the first time 
	static qoperator I; // define operators etc only once if they dont change 
	double val(0); 
	
	
	// just a dummy function to minimize (simultation would go here)
	val = x[0]*x[0];  
	
	
	calls++;  // record how many times the simulation is called 
	if (first_call) { first_call = false; }  
	cout << "calls:  "<<calls<<"  val:  "<<val<<"  x[0]:  "<<x[0]<<endl;

	return val;
}
 */
double sim(double* x){
	
	//---------------------------------------initialize parameters
	
	double Gamma = 0.01;
	double PT=0.1;
	double k=1;
	double w = 2*k;
	double alpha = x[0];

	
	const int seed = 12345;               // random seed
	const std::complex<double> ei(0,1);    // imaginary unit 
	const std::complex<double> one(1,0);   // real unit 
	const unsigned int N = 2,              // size of quantum system
	nt = 250;           // number of trajectories
	const double T = 1;                  // total time
	//dt = T/ns,
	//ns = 1000,         // number of time steps
	const double   dt = 0.001;     //time step
	const int ns = floor(T/dt);  //number of time intervals	const double rdt = std::sqrt(dt);
	const double rdt = std::sqrt(dt);

	double theta,theta_rec,dtheta;
	std::complex<double> axi,ayi,azi;
	double ax,ay,az,a;
	double avPground,RavTheta;
	double val(0); 
	
	
	static bool first_call = true;
	
	//--------------------------------initialise states
	
	qmatrix Up("up",N); 
	
	//--------------------------------initialize operators
	
	
	qoperator SigmaX("spinX",N);
	qoperator SigmaZ("spinZ",N);
	qmatrix I("maxmix",N); 
	
	qoperator Sigma = sin(alpha)*SigmaX+cos(alpha)*SigmaZ;
	qoperator Sigma2 = (sin(alpha)*SigmaX+cos(alpha)*SigmaZ)*(sin(alpha)*SigmaX+cos(alpha)*SigmaZ);
	
	
	
	
	//---------------------------------generate Wiener noise (dW's & dV's)
	static std::vector<double> noise(ns*nt);
	if (first_call) {
	make_noise(noise); 
	}
	
	//------------------------------------------------------------
	//                        Simulation
	//------------------------------------------------------------
	
	double t = 0,                      // initialise variables
	dW = 0;
	//  dV = 0;
	double n = PT/(1-2*PT);
	
	std::complex<double> d_recX = 0;
	
	std::vector<double> avTheta(ns+1);
	std::vector<double> Pground(ns+1);
	
	
	qmatrix rho;
	
	
	for (unsigned int i=0;i<ns+1;i++) {avTheta[i] = 0; Pground[i]=0;} 
	
	for (unsigned int j=0;j<nt;j++) {               // average multiple trajectories
		
		//rho = Up;
		ax = 0;
		az = 0;   // initialize ground state
		rho = I+(0.5)*ax*SigmaX+(0.5)*az*SigmaZ;
		avTheta[0] = 0;
		Pground[0] = nt/2;
		
		for (unsigned int i=0;i<ns;i++) {             // integrate a single trajectory
			
			
			ax += -(Gamma/2)*ax*dt;               
			az += -(Gamma)*(az + 1 - 2*n)*dt;                 // thermal noise evolution
			
		//	theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
		//	a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
		
			
			rho = I+(0.5)*ax*SigmaX+(0.5)*az*SigmaZ;	 // constructing density matrix
			
			dW = rdt*noise[i+j*ns];
			
			d_recX = 4*k*dt*trace(Sigma*rho) + std::sqrt(2*k)*dW;
			rho  += (-k*dt)*(Sigma2*rho + rho*Sigma2 + (-2.0)*Sigma*rho*Sigma)      // Sigma measurement
			+ d_recX*(Sigma*rho + rho*Sigma)
			+ k*(dW*dW - dt)*(Sigma2*rho + rho*Sigma2 + 2.0*Sigma*rho*Sigma);     // Milstien
			
			rho.normalise();
			
			
			axi = trace(rho*SigmaX);
			azi = trace(rho*SigmaZ);
			ax = axi.real();
			az = azi.real();
			theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
			a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
			
			
			if (theta > 0) {
				dtheta = -w*dt;
			} 
			if (theta < 0) {
				dtheta = w*dt;
			}
			
			if (std::abs(dtheta) > std::abs(theta) ) {
				theta = 0;
				ax = 0;
				az = -a;
			}   else {
				theta += dtheta;
				ax = std::sin(theta)*a;
				az = -std::cos(theta)*a;
			}          
		/*	//--------------------
			theta = 0;
			ax = 0;
			az = -a;
			//--------------------
		*/	
			
			avTheta[i+1] +=theta;
			Pground[i+1] += 0.5*(1-az);   // finds the probability of being in ground state
			
			
		}  //single trajectory integration
		
		
	}  //average mulitple trajectories
	
	
	for (unsigned int i=0;i<ns+1;i++) {avTheta[i] /= nt;
		Pground[i]/=nt;} // averages over trajectories
	
	
	
	
	//---------------
	// output results
	//---------------
	/*RavTheta = 0;
	avPground = 0;
	for (unsigned int i=20; i<ns+1; i++) {
		avPground += Pground[i];
        // RavTheta += abs(avTheta[i]);
	}
	avPground /= (ns-20+1);*/
   // RavTheta /= (ns-4000);
	//-----------------------
	// output simulation time
	//-----------------------
//	std::cout << "k=" << k << "  theta=" << alpha << " is done"  << std::endl;
    calls++;  	
    if (first_call) { first_call = false; } 
	val =1-Pground[ns];
	cout << "calls:  "<<calls<<"  avPground:  "<<val<<"  k:  "<<k<<"  alpha:  "<<alpha<<endl;
	return (val);
	
}

