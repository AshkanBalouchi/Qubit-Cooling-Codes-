
#include "lapack_interface.h"

long int diagonalize(long int Dim, std::complex<double>* Mat, double* Evals) { 
	
	// The unitary U comes back as Mat
	// If the input matrix is Mat = H, then 
	//        H = U diag(Evals) U^T
	//   and  H U = U diag(eig)     so that the columns of U are the eigenvectors
	
	long int info = 0;
	const char* do_vecs = "V", *up = "U";
	double work_dummy[1];
	double* rwork = new double[Dim*3-2]; 
	double* Mat_ptr = reinterpret_cast<double*>(Mat); 
	double* Ev_ptr = reinterpret_cast<double*>(Evals); 
	
	// query to see what the optimal work size is
	long int lwork = -1; 
	zheev_(do_vecs, up, &Dim, Mat_ptr, &Dim, Ev_ptr, work_dummy, &lwork, rwork, &info);
	
	lwork = static_cast<long int>(work_dummy[0]);
	double* work = new double[2*lwork]; 
	
	// calculate the eigenvalues and eigenvectors of A
	zheev_(do_vecs, up, &Dim, Mat_ptr, &Dim, Ev_ptr, work, &lwork, rwork, &info);
	
	delete[] work; 
	delete[] rwork; 
	
	return info;	
}

long int diagonalize_real(long int Dim, double* Mat, double* Evals) { 
	
	// The unitary U comes back as Mat
	// If the input matrix is Mat = H, then 
	//        H = U diag(Evals) U^T
	//   and  H U = U diag(eig)     so that the columns of U are the eigenvectors
	
	long int info = 0;
	const char* do_vecs = "V", *up = "U";
	double work_dummy[1];
	
	// query to see what the optimal work size is
	long int lwork = -1; 
	dsyev_(do_vecs, up, &Dim, Mat, &Dim, Evals, work_dummy, &lwork, &info);
	
	lwork = static_cast<long int>(work_dummy[0]);
	double* work = new double[lwork]; 
	
	// calculate the eigenvalues and eigenvectors of A
	dsyev_(do_vecs, up, &Dim, Mat, &Dim, Evals, work, &lwork, &info);
	
	delete[] work; 
	
	return info;	
}

long int tensor_prod_real(long int DimA, double* A, long int DimB, double* B, long int& DimC, double* C) { 
	// Note: B is the fast index 
	long int Dimtot = DimA*DimB; 
	DimC = Dimtot; 
	double *ptr, *fin; 
	for (long int k=0;k<DimA;k++) { 
		for (long int n=0;n<DimA;n++) {
			double A_nk = A[n + k*DimA];
			for (long int l=0;l<DimB;l++) { 
				ptr = &C[n*DimB + Dimtot*(l + k*DimB)]; fin = ptr + DimB;
				double* B_ptr = &B[l*DimB]; 
				while (ptr!=fin) { *(ptr++) = A_nk*(*(B_ptr++));  }
				//				for (long int m=0;m<DimB;m++) {
				//						C[m + n*DimB + Dimtot*(l + k*DimB)] += A[n + k*DimA]*B[m + l*DimB]; 
				//				} 	
			}
		}
	} 	
	return 0; 
}

long int tensor_prod(long int DimA, std::complex<double>* A, long int DimB, std::complex<double>* B, long int& DimC, std::complex<double>* C) { 
	// Note: B is the fast index 
	long int Dimtot = DimA*DimB; 
	DimC = Dimtot; 
	std::complex<double> *ptr, *fin; 
	for (long int k=0;k<DimA;k++) { 
		for (long int n=0;n<DimA;n++) {
			std::complex<double> A_nk = A[n + k*DimA];
			for (long int l=0;l<DimB;l++) { 
				ptr = &C[n*DimB + Dimtot*(l + k*DimB)]; fin = ptr + DimB;
				std::complex<double>* B_ptr = &B[l*DimB]; 
				while (ptr!=fin) { *(ptr++) = A_nk*(*(B_ptr++));  }
			}
		}
	} 	
	return 0; 
}

long int tprod_addto_real(long int DimA, double* A, long int DimB, double* B, double* C) { 
	// Note: B is the fast index 
	double *ptr, *fin; 
	long int Dimtot = DimA*DimB; 
	for (long int k=0;k<DimA;k++) { 
		for (long int n=0;n<DimA;n++) {
			double A_nk = A[n + k*DimA];
			for (long int l=0;l<DimB;l++) { 
				ptr = &C[n*DimB + Dimtot*(l + k*DimB)]; fin = ptr + DimB;
				double* B_ptr = &B[l*DimB]; 
				while (ptr!=fin) { *(ptr++) += A_nk*(*(B_ptr++));  }   // add_to 
			}
		}
	} 	
	return 0; 
}

long int tprod_addto(long int DimA, std::complex<double>* A, long int DimB, std::complex<double>* B, std::complex<double>* C) { 
	// Note: B is the fast index 
	std::complex<double> *ptr, *fin; 
	long int Dimtot = DimA*DimB; 
	for (long int k=0;k<DimA;k++) { 
		for (long int n=0;n<DimA;n++) {
			std::complex<double> A_nk = A[n + k*DimA];
			for (long int l=0;l<DimB;l++) { 
				ptr = &C[n*DimB + Dimtot*(l + k*DimB)]; fin = ptr + DimB;
				std::complex<double>* B_ptr = &B[l*DimB]; 
				while (ptr!=fin) { *(ptr++) += A_nk*(*(B_ptr++));  }   // add_to 
			}
		}
	} 	
	return 0; 
}


long int A_on_v(long int N, std::complex<double>* v, std::complex<double>* A) { 
	
	const char* Mat_op = "N"; 
	std::complex<double> alpha[1] = {1.0}, beta[1] = {0.0}; 
	double* alphav = reinterpret_cast<double*> (alpha); 
	double* betav  = reinterpret_cast<double*> (beta); 
	long int incx[1] = {1}, incy[1] = {1}; 
	std::complex<double> x[N]; 
	// copy v to x 
	std::complex<double> *p_x = x, *fin_x = &x[N-1], *p_v = v; fin_x++; 
	while (p_x!=fin_x) { *(p_x++) = *(p_v++); } 
	
	double* A_ptr = reinterpret_cast<double*> (A);
	double* v_ptr = reinterpret_cast<double*> (v); 
	double* x_ptr = reinterpret_cast<double*> (x);
	
	zgemv_(Mat_op, &N, &N, alphav, A_ptr, &N, x_ptr, incx, betav, v_ptr, incy); 
	return 0; 
}

long int AT_on_v(long int N, std::complex<double>* v, std::complex<double>* A) { 
	
	const char* Mat_op = "C"; 
	std::complex<double> alpha[1] = {1.0}, beta[1] = {0.0}; 
	double* alphav = reinterpret_cast<double*> (alpha); 
	double* betav  = reinterpret_cast<double*> (beta); 
	long int incx[1] = {1}, incy[1] = {1}; 
	std::complex<double> x[N]; 
	// copy v to x 
	std::complex<double> *p_x = x, *fin_x = &x[N-1], *p_v = v; fin_x++; 
	while (p_x!=fin_x) { *(p_x++) = *(p_v++); } 
	
	double* A_ptr = reinterpret_cast<double*> (A);
	double* v_ptr = reinterpret_cast<double*> (v); 
	double* x_ptr = reinterpret_cast<double*> (x);
	
	zgemv_(Mat_op, &N, &N, alphav, A_ptr, &N, x_ptr, incx, betav, v_ptr, incy); 
	return 0; 
}


long int A_times_B(long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B) { 
	
	const char *A_Mat_op = "N", *B_Mat_op = "N"; 
	std::complex<double> alpha[1] = {1.0}, beta[1] = {0.0}; 
	double* alphav = reinterpret_cast<double*> (alpha); 
	double* betav  = reinterpret_cast<double*> (beta); 
	
	double* A_ptr = reinterpret_cast<double*> (A);
	double* B_ptr = reinterpret_cast<double*> (B); 
	double* C_ptr = reinterpret_cast<double*> (C);
	
	zgemm_(A_Mat_op, B_Mat_op, &N, &N, &N, alphav, A_ptr, &N, B_ptr, &N, betav, C_ptr, &N); 
	return 0; 
}

long int A_times_BT(long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B) { 
	
	const char *A_Mat_op = "N", *B_Mat_op = "C"; 
	std::complex<double> alpha[1] = {1.0}, beta[1] = {0.0}; 
	double* alphav = reinterpret_cast<double*> (alpha); 
	double* betav  = reinterpret_cast<double*> (beta); 
	
	double* A_ptr = reinterpret_cast<double*> (A);
	double* B_ptr = reinterpret_cast<double*> (B); 
	double* C_ptr = reinterpret_cast<double*> (C);
	
	zgemm_(A_Mat_op, B_Mat_op, &N, &N, &N, alphav, A_ptr, &N, B_ptr, &N, betav, C_ptr, &N); 
	return 0; 
}

long int AT_times_B(long int N, std::complex<double>* C, std::complex<double>* A, std::complex<double>* B) { 
	
	const char *A_Mat_op = "C", *B_Mat_op = "N"; 
	std::complex<double> alpha[1] = {1.0}, beta[1] = {0.0}; 
	double* alphav = reinterpret_cast<double*> (alpha); 
	double* betav  = reinterpret_cast<double*> (beta); 
	
	double* A_ptr = reinterpret_cast<double*> (A);
	double* B_ptr = reinterpret_cast<double*> (B); 
	double* C_ptr = reinterpret_cast<double*> (C);
	
	zgemm_(A_Mat_op, B_Mat_op, &N, &N, &N, alphav, A_ptr, &N, B_ptr, &N, betav, C_ptr, &N); 
	return 0; 
}

long int Acols_timesby_v( long int N, std::complex<double>* v, std::complex<double>* A) {
	
	std::complex<double> *elem, *end_col, *v_ptr = v; 
	for (long int n=0;n<N;n++) {  
		elem = &(A[n*N]); end_col = &(A[(n+1)*N-1]); end_col++;  
		while (elem!=end_col) { *(elem++) *= *v_ptr; }
		v_ptr++; 
	}
	return 0; 
}

long int Arows_timesby_v( long int N, std::complex<double>* v, std::complex<double>* A){
	
	std::complex<double> *elem, *end_col, *v_ptr = v; 
	for (long int n=0;n<N;n++) {  
		elem = &(A[n*N]); end_col = &(A[(n+1)*N-1]); end_col++; 
		v_ptr = v; 
		while (elem!=end_col) { *(elem++) *= *(v_ptr++); }
	}
	return 0; 
}

long int copy_mat(long int N, std::complex<double>* A, std::complex<double>* B) {
	long int n_max = N*N; 
	for (long int n=0;n<n_max;n++) { A[n] = B[n]; } 
	return 0; 
}

std::complex<double> trace_mat(long int N, std::complex<double>* A) {
	std::complex<double> sum(0); 
	long int n_max = N*N; 
	for (long int n=0;n<n_max;n+=N+1) { sum += A[n]; } 
	return sum; 
}

long int trace_bath(long int Ntot, long int Nsys, std::complex<double>* psi, std::complex<double>* rho) { 
	
	long int Nbath = Ntot/Nsys; 
	std::complex<double> sum; 
	for (long int n=0;n<Nsys;n++) { 
		for (long int m=0;m<Nsys;m++) { 
			std::complex<double> *n_ztr = &psi[n*Nbath], *n_zin = &psi[(n+1)*Nbath-1]; n_zin++; 
			std::complex<double> *m_ztr = &psi[m*Nbath]; 
			sum = 0.0; 
			while (n_ztr!=n_zin) { sum += (*(m_ztr++))*std::conj(*(n_ztr++)); }
			rho[m + n*Nsys] = sum;
		}
	}
	
	return 0;  
}

long int trace_out_sys( std::complex<double>* rho_in, std::complex<double>* rho_out, long int sys, std::vector<unsigned int>& dims, long int& new_size){
	
	std::vector<unsigned int> new_dims(dims.size()-1); 
	for (int i=0;i<dims.size();i++) { 
		if (i < sys-1) { new_dims[i] = dims[i]; }
		if (i > sys-1) { new_dims[i-1] = dims[i]; }
	}

	new_size = 1;
	for (int i=0;i<new_dims.size();i++) { new_size *= new_dims[i]; }
  
	int chunk(1); for (int i=0;i<sys-1;i++) { chunk *= dims[i]; } 
	int dim(dims[sys-1]); 
	int num(1);  for (int i=sys;i<dims.size();i++) { num *= dims[i]; } 
	int old_size = chunk*dim*num; 
	
	// clear rho_out 
	for (int i=0;i<new_size*new_size;i++) { rho_out[i] = 0; } 
	
	for (int i=0;i<dim;i++) {
		
		for (int j1=0;j1<num;j1++) {
			for (int j2=0;j2<num;j2++) {
				
				for (int k1=0;k1<chunk;k1++) {
					for (int k2=0;k2<chunk;k2++) {
						
						rho_out[k2 + j2*chunk + new_size*(k1 + j1*chunk)]  
						+= rho_in[k2 + i*chunk + j2*chunk*dim + old_size*(k1 + i*chunk + j1*chunk*dim)];
						
					}}
				
			}}
		
	}
	
	dims = new_dims; 
	
	return 0; 
}



