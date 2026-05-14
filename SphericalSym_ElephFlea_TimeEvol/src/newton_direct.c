#include "1+1_Free_Boundaries.h"

//-----------------------------------------------------------------------------
void Jacobian(parameters par, derivs V, derivs ID, double **J)
{
	int j, l, ntotal = par.ntotal, Ntotal=ntotal-1;
	double *DX, *JDX;

	DX  = dvector(0, Ntotal);
	JDX = dvector(0, Ntotal);
	
	fill0_dvector(DX,   0, Ntotal);
	fill0_dvector(JDX,  0, Ntotal);

	for(j=0; j <= Ntotal; j++){
		DX[j] = 1.;
		J_times_DX(par, DX, V, ID, JDX);
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDX[l];
		DX[j] = 0.;
	}
	
	free_dvector(DX,  0, Ntotal);
	free_dvector(JDX, 0, Ntotal);
}
// -------------------------------------------------------------------------------
void DF_of_X(parameters par, double *X, derivs V, derivs ID, double **J)
{
	int j, k, ntotal = par.ntotal, Ntotal=ntotal-1;
	double *Xp, *Fp, *Xm, *Fm, eps = 5.e-06;
	
	Xp = dvector(0, Ntotal);
	Fp = dvector(0, Ntotal);
	Xm = dvector(0, Ntotal);
	Fm = dvector(0, Ntotal);
	
	for(j=0; j<= Ntotal; j++)
		Xp[j] = Xm[j] = X[j];

	for(j=0; j<= Ntotal; j++){
		Xp[j] += eps; 
		Xm[j] -= eps; 
		F_of_X(par, Xp, V, ID, Fp);
		F_of_X(par, Xm, V, ID, Fm);
		for(k=0; k<= Ntotal; k++)
			J[k][j] = 0.5*(Fp[k]-Fm[k])/eps;

		Xp[j] = Xm[j] = X[j];
	}
	
	free_dvector(Xp, 0, Ntotal);
	free_dvector(Xm, 0, Ntotal);
	free_dvector(Fm, 0, Ntotal);
	free_dvector(Fp, 0, Ntotal);
}
// -------------------------------------------------------------------------------
int newton_direct(parameters par, derivs ID, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int ntotal = par.ntotal, Ntotal=ntotal-1, *indx, iter=0, j, FLAG=-1;
	double *F, *DX, **J, d, norm;
	derivs V;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	indx  = ivector(0, Ntotal);
	J     = dmatrix(0, Ntotal, 0, Ntotal);
	allocate_derivs(&V, 0, Ntotal);
	
	F_of_X(par, X, V, ID, F);
	// PrintVector(F, 0, Ntotal);
	// exit(-1);

	norm = norm2(F, 0, Ntotal);
	copy_dvector(DX, F, 0, Ntotal); 
	
	if (Newton_verb == 1){
		fprintf(par.fout," Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm);
		fprintf(par.fout," ------------------------------------------------------------------\n");
		// exit(1);
	}

	if( par.EquationType_FLAG == 0 && par.ncount == 0){
		// Jacobian(par, V, ID, par.J_LU);
		Get_LinearOperator(par, par.J_LU);
		ludcmp(par.J_LU, Ntotal, par.indx_LU, &d, 0);
	}
	else if (par.EquationType_FLAG == 1){
		// Jacobian(par, V, ID, par.J_LU);
		Get_LinearOperator(par, par.J_LU);
		ludcmp(par.J_LU, Ntotal, par.indx_LU, &d, 0);
	}
	// Jacobian(par, V, ID, J);		
	// PrintMatrix("J.txt", J, 0, Ntotal, 0, Ntotal);
		
	// DF_of_X( par, X, V, ID, J);
	// PrintMatrix("J_DF.txt", J, 0, Ntotal, 0, Ntotal);
	// exit(-1);
	// ludcmp(J, Ntotal, indx, &d, 0);

	double tol=Newton_tol;//norm*Newton_tol;
//  	while(iter < Newton_itmin || (norm > tol  && iter < Newton_itmax)){
	while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
		if(norm < tol) FLAG++;
		iter += 1;

		if(par.EquationType_FLAG == 2){
			Jacobian(par, V, ID, par.J_LU);
			ludcmp(par.J_LU, Ntotal, par.indx_LU, &d, 0);
		}
		
		// Jacobian(par, V, ID, J);		
		// PrintMatrix("J.txt", J, 0, Ntotal, 0, Ntotal);
		
		// DF_of_X( par, X, V, ID, J);
		// PrintMatrix("J_DF.txt", J, 0, Ntotal, 0, Ntotal);
		// exit(-1);

		// ludcmp(J, Ntotal, indx, &d, 0);
		
		// lubksb(J, Ntotal, indx, DX, 0);
		
		lubksb(par.J_LU, Ntotal, par.indx_LU, DX, 0);

		for (j = 0; j < ntotal; j++) 
			X[j] -= DX[j];

		F_of_X(par, X, V, ID, F);
		norm = norm2(F, 0, Ntotal);
		copy_dvector(DX, F, 0, Ntotal); 
		
		if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
			fprintf(par.fout,"\n No Convergence of the Newton Raphson Method (norm =%3.15e). Now exiting to system.\n\n", norm);
			exit(1);
		}
		if (Newton_verb == 1){

			fprintf(par.fout, " Newton: iter = %3d \t |F| = %e (%3.2e), %d \n", iter, norm, tol, FLAG);
		}	
	}
	if(norm > tol)
		fprintf(par.fout,
		" Newton Raphson Method failed to converge to prescribed tolerance (thread = %d). Now move on to the next sequence element.\n", par.i_omp);

	free_dvector(F,    0, Ntotal);
	free_dvector(DX,   0, Ntotal);
	free_ivector(indx, 0, Ntotal);
	free_dmatrix(J,    0, Ntotal, 0, Ntotal);
	free_derivs(&V,    0, Ntotal);

	return iter;
}
//-----------------------------------------------------------------------------
void Allocate_LinearStatic_Jacobian(parameters *par){
	
	(*par).J_LU  = dmatrix(0, (*par).Ntotal, 0, (*par).Ntotal);
	(*par).indx_LU  = ivector(0,(*par).Ntotal);

	return;
}
//-----------------------------------------------------------------------------
void Free_LinearStatic_Jacobian(parameters *par){
	
	free_dmatrix((*par).J_LU, 0, (*par).Ntotal, 0, (*par).Ntotal);
	free_ivector((*par).indx_LU , 0, (*par).Ntotal);

	return;
}
