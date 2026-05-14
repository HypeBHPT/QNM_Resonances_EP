#include "1+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
int Newton_SDIRK_IG(parameters par, double x0, double *Mj, double *kj)
{
	int ntotal = par.ntotal_slice, Ntotal = ntotal-1, *indx, iter=0, j;
	double *F, *Dkj, **J, d, norm;
	
	F     = dvector(0, Ntotal);
	Dkj   = dvector(0, Ntotal);
	indx  = ivector(0, Ntotal);
	J     = dmatrix(0, Ntotal, 0, Ntotal);	
	
	fill0_dvector(kj, 0, Ntotal);
	
	F_of_kj_IG(par, x0, Mj, kj, F);
	
	norm = norm2(F, 0, Ntotal);
	copy_dvector(Dkj, F, 0, Ntotal); 
	
	if (Newton_SDIRK_verb == 1){
		printf("\t Initial Residual: \t |F| = %e \n", norm);
		printf("\t -------------------------------------------------------------------------\n");
	}
 	while(iter < Newton_SDIRK_itmin || (norm > Newton_SDIRK_tol  && iter < Newton_SDIRK_itmax)){
		iter += 1;

		/*************ATTENTION: DF_of_X_SDIRK_IG NOT WORKING***********/
		// DF_of_X_SDIRK_IG( par, x0, Mj, kj, /*Wslice, W1D,*/ J); 


		Jacobian_SDIRK_IG(par, x0, Mj, kj, J);
		
		//******************************************
		// ELLIPTIC LU SOLVER FOR INITIAL GUESS
		ludcmp(J, Ntotal, indx, &d, 0);
		lubksb(J, Ntotal, indx, Dkj, 0);
		//******************************************


		for (j = 0; j < ntotal; j++) 
			kj[j] -= Dkj[j];

		F_of_kj_IG(par, x0, Mj, kj, /*Wslice, W1D,*/ F);
		norm = norm2(F, 0, Ntotal);
		copy_dvector(Dkj, F, 0, Ntotal); 
		if (Newton_SDIRK_verb == 1){
// 			printf(" Newton_SDIRK: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e \n",
			if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				printf("\n\t No Convergence of the Newton_SDIRK Raphson Method. Now exiting to system.\n\n");
// 				exit(1);
			}
			printf("\t Newton_SDIRK: iter = %3d \t |F| = %e (%e)\n", iter, norm, Newton_SDIRK_tol);
		}	
	}
	if(norm > Newton_SDIRK_tol)
		printf(
		"\t Newton_SDIRK Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,      0, Ntotal);
	free_dvector(Dkj,    0, Ntotal);
	free_ivector(indx,   0, Ntotal);
	free_dmatrix(J,      0, Ntotal, 0, Ntotal);


	return iter;
}//-----------------------------------------------------------------------------
void Jacobian_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J){
	int j, l, ntotal = par.ntotal_slice, Ntotal = ntotal - 1;
	double *Dkj, *JDkj;
	
	Dkj  = dvector(0, Ntotal);
	JDkj = dvector(0, Ntotal);
	
	fill0_dvector(Dkj,   0, Ntotal);
	fill0_dvector(JDkj,  0, Ntotal);

	for(j=0; j <= Ntotal; j++){
		Dkj[j] = 1.;
		J_times_Dkj_IG(par, x0, Mj, kj, Dkj, JDkj);
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDkj[l];
		Dkj[j] = 0.;
	}
	
	free_dvector(Dkj,  0, Ntotal);
	free_dvector(JDkj, 0, Ntotal);
}
// -------------------------------------------------------------------------------
void DF_of_X_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J){

	/*************ATTENTION: DF_of_X_SDIRK_IG NOT WORKING***********/
	int j, k, ntotal = par.ntotal_slice, Ntotal = ntotal - 1;
	double *kjp, *Fp, *kjm, *Fm, eps = 5.e-06;
	
	kjp = dvector(0, Ntotal);
	Fp  = dvector(0, Ntotal);
	kjm = dvector(0, Ntotal);
	Fm  = dvector(0, Ntotal);
	
	for(j=0; j<= Ntotal; j++)
		kjp[j] = kjm[j] = kj[j];

	for(j=0; j<= Ntotal; j++){
		kjp[j] += eps; 
		kjm[j] -= eps; 
		F_of_kj_IG(par, x0, Mj, kjp, Fp);
		F_of_kj_IG(par, x0, Mj, kjm, Fp);
		for(k=0; k<= Ntotal; k++)
			J[k][j] = 0.5*(Fp[k]-Fm[k])/eps;
		kjp[j] = kjm[j] = kj[j];
	}
	
	free_dvector(kjp, 0, Ntotal);
	free_dvector(kjm, 0, Ntotal);
	free_dvector(Fm,  0, Ntotal);
	free_dvector(Fp,  0, Ntotal);

	return;
}
