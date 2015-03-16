a
/*
 *  check.cc
 *  
 *
 *  Created by Ashkan Balouchi on 10/17/12.
 *  Copyright 2012 LSU. All rights reserved.
 *
 */

//#include "check.h"

/*
 *  FTN-function.cc
 *  
 *
 *  Created by Ashkan Balouchi on 8/6/12.
 *  Copyright 2012 LSU. All rights reserved.
 *
 */

//#include "FTN-function.h"
#include "qobjects.h"
//#include "FTN-opt.h"
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

double SteadyStateProbGround (double stop,double Gamma,double PT,double w,double k,double alpha)
{
	//---------------------------------------initialize parameters
	
	const int seed = 12345;               // random seed
	const std::complex<double> ei(0,1);    // imaginary unit 
	const std::complex<double> one(1,0);   // real unit 
	const unsigned int N = 2,              // size of quantum system
	nt = 50;           // number of trajectories
	const double T=(1/k)*stop;  //total time
	
	 const double   dt = 0.001;     //time step
	const int ns = floor(T/dt);  //number of time intervals   // constant dt
	/*
	 const int ns=1000;           //number of time intervals
	const double dt=T/ns;         
	  */                                                    // constant number of intervals
	
	const double rdt = std::sqrt(dt);

	 //	double k;                                 // measurement rate
	//	int w;                          // feedback rate
	//	double Gamma = 0.2;                   // damping rate
	//	double PT = 0.1;                      // population of excited state
	double theta,theta_rec,dtheta;
	std::complex<double> axi,ayi,azi;
	double ax,ay,az,a;
	double avPground;
	
	
	static bool x = true;
	//--------------------------------initialise files
	
	
	
	
	//--------------------------------initialise states
	
	qmatrix Up("up",N); 
	
	//--------------------------------initialize operators
	
	
	qoperator SigmaX("spinX",N);
	//qoperator SigmaX2 = SigmaX*SigmaX;
	qoperator SigmaZ("spinZ",N);
	//qoperator SigmaZ2 = SigmaZ*SigmaZ;
	qoperator Sigma = sin(alpha)*SigmaX+cos(alpha)*SigmaZ;
	qoperator Sigma2 = (sin(alpha)*SigmaX+cos(alpha)*SigmaZ)*(sin(alpha)*SigmaX+cos(alpha)*SigmaZ);
	qmatrix I("maxmix",N); 
	
	
	//---------------------------------generate Wiener noise (dW's & dV's)
	
	static std::vector<double> noise(ns*nt);
	if (x) {
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
		
		ax = 0;
		az = 0;   // initialize ground state
		rho = I+(0.5)*ax*SigmaX+(0.5)*az*SigmaZ;
		avTheta[0] = 0;
		Pground[0] = nt/2;
		
		for (unsigned int i=0;i<ns;i++) {             // integrate a single trajectory
			
			
			ax += -(Gamma/2)*ax;               
			az += -(Gamma)*(az + 1 - 2*n);                 // thermal noise evolution
			
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
		/*		//------------------------------
			theta = 0;
			ax = 0;
			az = -a;
			//------------------------------
						 */

			//	avTheta[i+1] +=theta;
			Pground[i+1] += 0.5*(1-az);   // finds the probability of being in ground state
			
			
		}  //single trajectory integration
		
		//	std::cout <<"for k="<< k <<"  finished traj no. " << j+1 <<"  theta_rec:"<<theta<<"   ax:"<<ax<<"   az:"<<az<< std::endl;
		
	}  //average mulitple trajectories
	
	
	for (unsigned int i=0;i<ns+1;i++) {avTheta[i] /= nt;
		Pground[i]/=nt;} // averages over trajectories
	
	
	
	
	//---------------
	// output results
	//---------------
	
	/* avPground = 0;
	 for (unsigned int i=ns-19; i<ns+1; i++) {
	 avPground += Pground[i];
	 }
	 avPground /= 20;*/
	
	/*
	std::ofstream result("Check-Gamma0.PT0.w0.AlphaPiHalf.G0.01.PT0.1.k50.txt");
	//std::ofstream result("Check-Gamma0.PT0.w0.AlphaZero.G0.01.PT0.1.k100.txt");

	
	for (int i=0; i<ns+1; i++) {
		result << i <<"  "<< Pground[i] << std::endl;
		
	}
	
	result.close();*/
	
	
	
	//-----------------------
	// output simulation time
	//-----------------------
//std::cout << "k=" << k << "  theta=" << alpha << " is done"  << std::endl;
	if (x) {
		x = false;
	}
	avPground=Pground[ns];
	return (1000*avPground);
	
}
//-----------------------------------
//main simulation
//-----------------------------------



int main ()
{ 
	double PT=0.1,
	G=0.01,
	
	//ang=1.5707963,
	//ang=0,
	ang,
	k=1;
	double w=2*k;
	double stop = 0.2;
	
	time_t time0, *time_ptr = 0;       // timing initialisation`
	time0 = time(time_ptr);
	
	//double P = SteadyStateProbGround(G,PT,w,k,ang);
	
	std::ofstream result("Gamma0.01.PT0.1k1stop0.2dt.txt");
	for (ang=0; ang<6.4; ang+=0.1)
	 {
	 double P = SteadyStateProbGround(stop,G,PT,w,k,ang);
	 //cout << "k = "<< k<<" alpha= "<<ang<< "  avPground = " << P << std::endl;
	 result << ang  << "  " << P << std::endl;
		 cout << ang<<"  done"<<endl;
	 }
			
	result.close();
	std::cout << std::endl << "Simulation took " << difftime(time(time_ptr),time0) 
	<< " seconds" << std::endl;
	
	return 0;
}
