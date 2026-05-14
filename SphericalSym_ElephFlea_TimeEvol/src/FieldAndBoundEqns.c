#include "1+1_Free_Boundaries.h"
// -------------------------------------------------------------------------------
void LinearOperator_Matrix(parameters par, int idom,  int i0, int i1, double **M){
	int II_Re_psi       = Index(par, idom, 0, i1, i0),
		// II_Im_psi       = Index(par, idom, 1, i1, i0),
		II_Re_dpsi_dtau = Index(par, idom, 1, i1, i0),
		// II_Im_dpsi_dtau = Index(par, idom, 3, i1, i0),	
	    j0, j1, N0 = par.N0, N1 = par.N1[idom];
	
	for(j0=0;j0<=N0; j0++){
		for(j1=0;j1<=N1; j1++){
			int JJ_Re_psi       = Index(par, idom, 0, j1, j0),
				// JJ_Im_psi       = Index(par, idom, 1, j1, j0),
				JJ_Re_dpsi_dtau = Index(par, idom, 1, j1, j0);
				// JJ_Im_dpsi_dtau = Index(par, idom, 3, j1, j0);

				//Field = 0 --- Re(phi) couples with Re(dphi_dtau) (Field 2) via first order reduction
				M[II_Re_psi][JJ_Re_psi]  = - par.Dtau_2D[II_Re_psi][JJ_Re_psi];
				M[II_Re_psi][JJ_Re_dpsi_dtau] += par.Id_2D[II_Re_dpsi_dtau][JJ_Re_dpsi_dtau];

				// //Field = 1 --- Im(phi) couples with Im(dphi_dtau) (Field 3) via first order reduction
				// M[II_Im_psi][JJ_Im_psi]  =  - par.Dtau_2D[II_Im_psi][JJ_Im_psi];
				// M[II_Im_psi][JJ_Im_dpsi_dtau] +=  par.Id_2D[II_Im_dpsi_dtau][JJ_Im_dpsi_dtau];


				//Field = 2 --- Re(dphi_dtau) couples with Re(phi) (Field 0) via wave equation
				M[II_Re_dpsi_dtau][JJ_Re_dpsi_dtau]   =  -par.W[II_Re_dpsi_dtau]*par.Dtau_2D[II_Re_dpsi_dtau][JJ_Re_dpsi_dtau] 
				       								    + par.beta_1[II_Re_dpsi_dtau]*par.Dsig_2D[II_Re_dpsi_dtau][JJ_Re_dpsi_dtau]
													    + par.beta_0[II_Re_dpsi_dtau]*par.Id_2D[II_Re_dpsi_dtau][JJ_Re_dpsi_dtau];
				M[II_Re_dpsi_dtau][JJ_Re_psi]  		 +=  par.alpha_2[II_Re_dpsi_dtau]*par.Dsigsig_2D[II_Re_psi][JJ_Re_psi]
														+ par.alpha_1[II_Re_dpsi_dtau]*par.Dsig_2D[II_Re_psi][JJ_Re_psi]
														+ par.alpha_0[II_Re_dpsi_dtau]*par.Id_2D[II_Re_psi][JJ_Re_psi];

				// //Field = 3 --- Im(dphi_dtau) couples with Im(phi) (Field 1) via wave equation
				// M[II_Im_dpsi_dtau][JJ_Im_dpsi_dtau]   = -par.W[II_Im_dpsi_dtau]*par.Dtau_2D[II_Im_dpsi_dtau][JJ_Im_dpsi_dtau] 
				//        									+ par.beta_1[II_Im_dpsi_dtau]*par.Dsig_2D[II_Im_dpsi_dtau][JJ_Im_dpsi_dtau]
				// 									    + par.beta_0[II_Im_dpsi_dtau]*par.Id_2D[II_Im_dpsi_dtau][JJ_Im_dpsi_dtau];
				// M[II_Im_dpsi_dtau][JJ_Im_psi]  		 +=   par.alpha_2[II_Im_dpsi_dtau]*par.Dsigsig_2D[II_Im_psi][JJ_Im_psi]
				// 										+ par.alpha_1[II_Im_dpsi_dtau]*par.Dsig_2D[II_Im_psi][JJ_Im_psi]
				// 										+ par.alpha_0[II_Im_dpsi_dtau]*par.Id_2D[II_Im_psi][JJ_Im_psi];		
		}
	}
	
	
 	return;
}
// -------------------------------------------------------------------------------
void LinearBoundaryConditions_Matrix(parameters par, int idom,  int i0, int i1, double **M){
	// LinearOperator_Matrix(par, idom,  i0, i1, M);
	int iF, j0, j1, N0 = par.N0, N1 = par.N1[idom];

	if(i1==N1 && idom ==0) LinearOperator_Matrix(par, idom,  i0, i1, M);
	else if(i1==0 && idom==nDoms-1) LinearOperator_Matrix(par, idom,  i0, i1, M);


	else{
		for(iF=0; iF<nFields; iF++){
			int II = Index(par, idom, iF, i1, i0);
			
			if(i1==N1 && idom!=0){
				int JJ = Index(par, idom-1, iF, 0, i0);
				M[II][II] = par.Id_2D[II][II];
				M[II][JJ] = -par.Id_2D[JJ][JJ];;

			}
			else if(i1==0 && idom!=nDoms-1){
				int II_nextDom = Index(par, idom+1, iF, par.N1[idom+1], i0);
				for(j0=0;j0<=N0; j0++){
					for(j1=0;j1<=N1; j1++){
						int JJ = Index(par, idom, iF, j1, j0), JJ_nextDom = Index(par, idom+1, iF, j1, j0);
						M[II][JJ] = -par.Dsig_2D[II][JJ];
						M[II][JJ_nextDom] = par.Dsig_2D[II_nextDom][JJ_nextDom];
					}
				}
			}

		}

	}
	

		
	
	
 	return;
}

double window_function(parameters par, double t){
	double window =1.;

	return window;
}
//----------------------------------------		
double jump_field(parameters par, int iF, double x0){
	double jump = 0.;	

	jump = iF%2 == 0? 0.: 0.;
	return jump;
}
//----------------------------------------		
double jump_field_dtau(parameters par, int iF, double x0){
	double jump = 0.;	

	return jump;
}	
//----------------------------------------			    	
double jump_derivative(parameters par, int iF, double x0){
	double jump = 0.;
	
	return jump;
}
//----------------------------------------	
double jump_derivative_dtau(parameters par, int iF, double x0){
	double jump = 0.;
	
	return jump;
}
//----------------------------------------	
double complex Get_Source(parameters par, int iDom,  double x0, int i1){
	double t = get_t_from_x0(par, x0);
	double chi_1 = par.grid_points_x1[iDom][i1];
	
	double complex Source;
	




	if(par.HomoEq == 0){
		Source = 0;
	}
	else if( par.HomoEq == 0){
		double complex ss=par.ss, S; //ss_conj = conj(ss), S_conj; 

		S = Clenshaw_Chebyshev(par.cheb_real_Source[iDom], par.N_input_Source[iDom], chi_1) 
			+I*Clenshaw_Chebyshev(par.cheb_imag_Source[iDom], par.N_input_Source[iDom], chi_1);
		
		//S_conj = conj(S);
		

		double sigma, dx_dsigma;
		get_sigma(par, iDom, i1, &sigma, &dx_dsigma);
		Source = cexp(ss * t)*S;// + cexp(ss_conj * t)*S_conj; 		
	
		
	}
	else{
		fprintf(stderr, "Error in function SecondOrderSoruce: PertOrder = %d not implemented\n", par.HomoEq);
		exit(-1);
	}
	
	return Source;
}
// -------------------------------------------------------------------------------
void F_of_XY(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, double *F)
{
	derivs W_tau_sigma;

	double sigma, dx;
	get_sigma(par, iDom, i1, &sigma, &dx);

	allocate_derivs(&W_tau_sigma, 0, nFields-1);
	scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, W, W_tau_sigma, nFields);
	
	double  w, a2, a1, a0, b1, b0,
			Re_u, Re_u_dsigma, Re_u_d2sigma, Re_u_dt,
			Re_du, Re_du_dt, Re_du_dsigma;
			// Im_u, Im_u_dsigma, Im_u_d2sigma, Im_u_dt,
			// Im_du, Im_du_dt, Im_du_dsigma;
	double complex Source = Get_Source(par, iDom, x0, i1);
	

	a2 = func_alpha2(par, sigma);
	a1 = func_alpha1(par, sigma);
	a0 = func_alpha0(par, sigma);
	b1 = func_beta1(par, sigma);
	b0 = func_beta0(par, sigma);
	w  = func_w(par, sigma);

	Re_u = W_tau_sigma.X[0]; Re_u_dsigma = W_tau_sigma.X1[0]; Re_u_d2sigma = W_tau_sigma.X11[0]; Re_u_dt = W_tau_sigma.Y[0]; 
	// Im_u = W_tau_sigma.X[1]; Im_u_dsigma = W_tau_sigma.X1[1]; Im_u_d2sigma = W_tau_sigma.X11[1]; Im_u_dt = W_tau_sigma.Y[1]; 
	Re_du = W_tau_sigma.X[1]; Re_du_dt = W_tau_sigma.Y[1]; Re_du_dsigma = W_tau_sigma.X1[1];	
	// Im_du = W_tau_sigma.X[3]; Im_du_dt = W_tau_sigma.Y[3]; Im_du_dsigma = W_tau_sigma.X1[3];

	F[0] = Re_du - Re_u_dt;
	// F[1] = Im_du - Im_u_dt;
	F[1] = -w*Re_du_dt + (a2*Re_u_d2sigma + a1*Re_u_dsigma + a0*Re_u) + (b1*Re_du_dsigma + b0*Re_du) - creal(Source);
	// F[3] = -w*Im_du_dt + (a2*Im_u_d2sigma + a1*Im_u_dsigma + a0*Im_u) + (b1*Im_du_dsigma + b0*Im_du) - cimag(Source);


	free_derivs(&W_tau_sigma, 0, nFields-1);
	return;
}
// -------------------------------------------------------------------------------
void JDX_of_XY(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *JDX)
{
	derivs DW_tau_sigma;
	double sigma, dx;
	get_sigma(par, iDom, i1, &sigma, &dx);
	
	allocate_derivs(&DW_tau_sigma, 0, nFields-1);
	scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, DW, DW_tau_sigma, nFields);
	
	double  w, a2, a1, a0, b1, b0,
			D_Re_u, D_Re_u_dsigma, D_Re_u_d2sigma, D_Re_u_dt,
			D_Re_du, D_Re_du_dt, D_Re_du_dsigma;
			// D_Im_u, D_Im_u_dsigma, D_Im_u_d2sigma, D_Im_u_dt,
			// D_Im_du, D_Im_du_dt, D_Im_du_dsigma;

	a2 = func_alpha2(par, sigma);
	a1 = func_alpha1(par, sigma);
	a0 = func_alpha0(par, sigma);
	b1 = func_beta1(par, sigma);
	b0 = func_beta0(par, sigma);
	w  = func_w(par, sigma);

	D_Re_u = DW_tau_sigma.X[0]; D_Re_u_dsigma = DW_tau_sigma.X1[0]; D_Re_u_d2sigma = DW_tau_sigma.X11[0]; D_Re_u_dt = DW_tau_sigma.Y[0]; 
	// D_Im_u = DW_tau_sigma.X[1]; D_Im_u_dsigma = DW_tau_sigma.X1[1]; D_Im_u_d2sigma = DW_tau_sigma.X11[1]; D_Im_u_dt = DW_tau_sigma.Y[1]; 

	D_Re_du = DW_tau_sigma.X[2]; D_Re_du_dt = DW_tau_sigma.Y[2]; D_Re_du_dsigma = DW_tau_sigma.X1[2];	
	// D_Im_du = DW_tau_sigma.X[3]; D_Im_du_dt = DW_tau_sigma.Y[3]; D_Im_du_dsigma = DW_tau_sigma.X1[3];

	JDX[0] =  D_Re_du - D_Re_u_dt;
	// JDX[1] =  D_Im_du - D_Im_u_dt;
	JDX[1] =  -w*D_Re_du_dt + (a2*D_Re_u_d2sigma + a1*D_Re_u_dsigma + a0*D_Re_u) + (b1*D_Re_du_dsigma + b0*D_Re_du);	
	// JDX[3] =  -w*D_Im_du_dt + (a2*D_Im_u_d2sigma + a1*D_Im_u_dsigma + a0*D_Im_u) + (b1*D_Im_du_dsigma + b0*D_Im_du);
	
	free_derivs(&DW_tau_sigma,      0, nFields-1);
	
	return;
	                                     
}
// -------------------------------------------------------------------------------
// void Jump_Field(parameters par, int iDom, double x0, int i1, derivs V, derivs ID, double *F){
// 	derivs W_grid_Dom, W_grid_prevDom;
// 	allocate_derivs(&W_grid_Dom, 0, NFields);
// 	allocate_derivs(&W_grid_prevDom, 0, NFields);
	
// 	Get_Fields2D_atGrid( par, iDom,  i0, i1, ID,  V, W_grid_Dom);
// 	Get_Fields2D_atGrid( par, iDom-1, i0, 0, ID,  V, W_grid_prevDom);
					
// 	for(iField=0; iField < nFields; iField++){
// 	 F_2D_Equations[iField] = W_grid.X[iField] - W_grid_prevDom.X[iField] - jump_field(par, iField, x0);
// 	}

// 	free_derivs(&W_grid_Dom, 0, NFields);
// 	free_derivs(&W_grid_prevDom, 0, NFields);
// 	return;
// }
// -------------------------------------------------------------------------------
void Boundary_Data_Left(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, double *F){
	F_of_XY(par, iDom, x0, i1, W, W1D, F);
	return;
}
// -------------------------------------------------------------------------------
void LinearBoundary_Data_Left(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *JDX){
	JDX_of_XY(par, iDom, x0, i1, W, W1D, DW, DW1D, JDX);
	return;
}
// -------------------------------------------------------------------------------
void Boundary_Data_Right(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, double *F){
	F_of_XY(par, iDom, x0, i1, W, W1D, F);
	return;
}
// -------------------------------------------------------------------------------
void LinearBoundary_Data_Right(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *JDX){
	JDX_of_XY(par, iDom, x0, i1, W, W1D, DW, DW1D, JDX);
	return;
}
// -------------------------------------------------------------------------------
void FBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, double *F)
{
	 
	return;
}
// -------------------------------------------------------------------------------
void JDXBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, derivs DWslice, derivs DW1D, double *JDX)
{
	return;
}
// -------------------------------------------------------------------------------
