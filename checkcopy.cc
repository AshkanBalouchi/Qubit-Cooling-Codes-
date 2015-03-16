/*
 *  check.cc
 *  
 *
 *  Created by Ashkan Balouchi on 2/20/13.
 *  Copyright 2013 LSU. All rights reserved.
 *
 */

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

//#include "mpi.h"




double  SteadyStateProbGround(double x){
	
	//---------------------------------------initialize parameters
	
	
	double Gamma = 0.1;
	double PT=0.1;
	double w =60;
	double k=1;
	double stop=2.5;
	
	
	double alpha ;

	
	const int seed = 123;    // random seed
	std::srand(seed);
	const std::complex<double> ei(0,1);    // imaginary unit 
	const std::complex<double> one(1,0);   // real unit 
	const unsigned int N = 2,              // size of quantum system
	nt = 1000;                             // number of trajectories
	const double T=(1/k)*stop;             // total time
	const double dt = 0.0001;               // time step
	const int ns = static_cast<int>(floor(T/dt));  // number of time intervals  // constant dt
	
	const double rdt = std::sqrt(dt);
	double theta,theta_rec,dtheta;
	std::complex<double> axi,ayi,azi;
	double ax,ay,az,a;
	double avPground,RavTheta;
	double val; 
	
	
	static bool first_call = true;
	
	//--------------------------------initialise states
	
	qmatrix Up("up",N); 
	
	//--------------------------------initialize operators
	
	
	qoperator SigmaX("spinX",N);
	qoperator SigmaZ("spinZ",N);
	qmatrix I("maxmix",N); 
	
	
	//---------------------------------generate Wiener noise (dW's & dV's)
	static std::vector<double> noise(ns*nt);
	if (first_call) {
		make_noise(noise); 
	}
	
	//------------------------------------------------------------
	//                        Simulation
	//------------------------------------------------------------
	
			double b;        
	
			int c=-30;  
			//b=-0.682456+0.189997 * exp( -0.518407*w );        
			b=0;		
			double t = 0,                      // initialise variables
			dW = 0;
			// double n = PT/(1-2*PT);
			double n=0.1;
			
			
			std::complex<double> d_recX = 0;
			
			std::vector<double> avTheta(ns+1);
			std::vector<double> Pground(ns+1) , Pground_all(ns+1), axa(ns+1) , axa_all(ns+1);
			double axa_final;
			
			qmatrix rho;
			
			
			for (unsigned int i=0;i<ns+1;i++) {avTheta[i] = 0; Pground[i]=0; axa[i]=0;} 
			
			for (unsigned int j=0;j<nt;j++) {               // average multiple trajectories
				
				//rho = Up;
				ax = 0;
				az =-0.83334;   // initialize ground state
				rho = I+(0.5)*ax*SigmaX+(0.5)*az*SigmaZ;
				//avTheta[0] = 0;
				//Pground[0] = nt/2;
				
				for (unsigned int i=0;i<ns;i++) {             // integrate a single trajectory
					
					
					ax += -(Gamma/2)*ax*dt;               
					az += -(Gamma)*(az + 1 - 2*n)*dt;                 // thermal noise evolution
					
					//	theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
					theta = std::atan2(ax,-az);
                	a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
					
					
					//	alpha=-0.7*theta;
					//             alpha=1.570796;
					//                alpha=0;
					//              alpha =-b*theta;
					alpha =(c/60.00)*M_PI-b*theta;
					
					/*	for (int h=0; h<Nspace; h++) {
					 alpha+=x[h]*pow(theta,h);
					 }	
					 */	
					
					qoperator Sigma = sin(alpha)*SigmaX+cos(alpha)*SigmaZ;
					qoperator Sigma2 = (sin(alpha)*SigmaX+cos(alpha)*SigmaZ)*(sin(alpha)*SigmaX+cos(alpha)*SigmaZ);
					
					
					
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
					//			theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
					theta = std::atan2(ax,-az);
					
					a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
				
				/*	
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
					}      */   
				 
					//--------------------
						theta = 0;
					 ax = 0;
					 az = -a;
					//--------------------
					
					axa[i+1] +=std::abs(ax);	
					Pground[i+1] += 0.5*(1-az);   // finds the probability of being in ground state
					
					
				}  //single trajectory integration
				
				
				cout << "traj no: "<<j<<std::endl;
			}  //average mulitple trajectories
			
		
    val =Pground[ns]/nt;
	if (first_call) { first_call = false; } 
	return val;
	
}



int main(int argc, char* argv[])
{ 	time_t time0, *time_ptr = 0;       // timing initialisation`
	time0 = time(time_ptr);
	
	
	double P = SteadyStateProbGround(1);
	
		std::cout << "val= "<< P << std::endl;
		std::cout << std::endl << "Simulation took " << difftime(time(time_ptr),time0) 
		<< " seconds" << std::endl;
		return 0;
}
