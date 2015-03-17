/*
 *  sim.h
 *  
 *
 *  Created by Ashkan Balouchi on 10/8/12.
 *  Copyright 2012 LSU. All rights reserved.
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

//***************************************************************************************
//desctritizor function
//***************************************************************************************
/*
void func(double Bloch_len, double Bloch_ang,double mapv[25][2], double &k,double &alpha){
	
	double angr = (Bloch_ang+3.14159)/1.256637;
	double lenr = Bloch_len/0.2;
	int count;
	count=floor(angr)*5+floor(lenr)+1;
	k = mapv[count][1];
	alpha = mapv[count][2];
	
	
}



//***************************************************************************************
//main simulation function
//***************************************************************************************

double SteadyStateProbGround (double Gamma,double PT,double w,double mapv[25][2],double &k,double &alpha)
{
	//---------------------------------------initialize parameters
	
	const int seed = 123456;               // random seed
	const std::complex<double> ei(0,1);    // imaginary unit 
	const std::complex<double> one(1,0);   // real unit 
	const unsigned int N = 2,              // size of quantum system
	ns = 6000,         // number of time steps
	nt = 10;           // number of trajectories
	const double T = 6,                   // total time
	dt = T/ns,
	rdt = std::sqrt(dt);
	double theta,theta_rec,dtheta;
	std::complex<double> axi,ayi,azi;
	double ax,ay,az,a;
	double avPground;
	
	
	static bool firstcall = true;
	
	//--------------------------------initialise states
	
	qmatrix Up("up",N); 
	
	//--------------------------------initialize operators
	
	
	qoperator SigmaX("spinX",N);
	qoperator SigmaZ("spinZ",N);
	qmatrix I("maxmix",N); 
	
	
	//---------------------------------generate Wiener noise (dW's & dV's)
	
	std::vector<double> noise(ns*nt);
	if (firstcall = true) {
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
		
		rho = Up;
		ax = 0;
		az = -1;   // initialize ground state
		avTheta[0] = 0;
		Pground[0] = nt;
		
		for (unsigned int i=0;i<ns;i++) {             // integrate a single trajectory
			
			
			ax += -(Gamma/2)*ax;               
			az += -(Gamma)*(az + 1 - 2*n);                 // thermal noise evolution
			
			theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
			a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
			
			func(a,theta,mapv,k,alpha);
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
			
			
			
			avTheta[i+1] +=theta;
			Pground[i+1] += 0.5*(1-az);   // finds the probability of being in ground state
			
			
		}  //single trajectory integration
		
		
	}  //average mulitple trajectories
	
	
	for (unsigned int i=0;i<ns+1;i++) {avTheta[i] /= nt;
		Pground[i]/=nt;} // averages over trajectories
	
	
	
	
	//---------------
	// output results
	//---------------
	
	avPground = 0;
	for (unsigned int i=2000; i<ns+1; i++) {
		avPground += Pground[i];
	}
	avPground /= (ns-2000);
	//-----------------------
	// output simulation time
	//-----------------------
	std::cout << "k=" << k << "  theta=" << alpha << " is done"  << std::endl;
	firstcall = false;
	return (avPground);
	
}*/

//double sim(double* x)
//{
	//---------------------------------------initialize parameters
int main()
{
	double Gamma = 0.1;
	double PT=0.2;
	double w = 50;
	double alpha = 0;
	double k=5;
	
	const int seed = 123456;               // random seed
	const std::complex<double> ei(0,1);    // imaginary unit 
	const std::complex<double> one(1,0);   // real unit 
	const unsigned int N = 2,              // size of quantum system
	ns = 6000,         // number of time steps
	nt = 10;           // number of trajectories
	const double T = 6,                   // total time
	dt = T/ns,
	rdt = std::sqrt(dt);
	double theta,theta_rec,dtheta;
	std::complex<double> axi,ayi,azi;
	double ax,ay,az,a;
	double avPground;
	
	
	static bool firstcall = true;
	
	//--------------------------------initialise states
	
	qmatrix Up("up",N); 
	
	//--------------------------------initialize operators
	
	
	qoperator SigmaX("spinX",N);
	qoperator SigmaZ("spinZ",N);
	qmatrix I("maxmix",N); 
	
	
	//---------------------------------generate Wiener noise (dW's & dV's)
	
	std::vector<double> noise(ns*nt);
	if (firstcall = true) {
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
		
		rho = Up;
		ax = 0;
		az = -1;   // initialize ground state
		avTheta[0] = 0;
		Pground[0] = nt;
		
		for (unsigned int i=0;i<ns;i++) {             // integrate a single trajectory
			
			
			ax += -(Gamma/2)*ax;               
			az += -(Gamma)*(az + 1 - 2*n);                 // thermal noise evolution
			
			theta = std::atan(ax/-az);       //the angle of the vector with +z axes			
			a = std::sqrt(az*az + ax*ax);    //length of the Bloch vector
			
			//	func(a,theta,mapv,k,alpha);
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
			
			
			
			avTheta[i+1] +=theta;
			Pground[i+1] += 0.5*(1-az);   // finds the probability of being in ground state
			
			
		}  //single trajectory integration
		
		
	}  //average mulitple trajectories
	
	
	for (unsigned int i=0;i<ns+1;i++) {avTheta[i] /= nt;
		Pground[i]/=nt;} // averages over trajectories
	
	
	
	
	//---------------
	// output results
	//---------------
	
	avPground = 0;
	for (unsigned int i=2000; i<ns+1; i++) {
		avPground += Pground[i];
	}
	avPground /= (ns-2000);
	//-----------------------
	// output simulation time
	//-----------------------
	//	std::cout << "k=" << k << "  theta=" << alpha << " is done"  << std::endl;
    //calls++;  	
    firstcall = false;
//	return (-avPground);
	cout << "avPground=  "<<avPground<<endl;
return (0);
}

