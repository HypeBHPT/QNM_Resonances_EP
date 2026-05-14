#include "1+1_Free_Boundaries.h"

//   -----------------------------------------------
int main() 
{
	
		parameters par;
		par.fout = stdout;

		
		double *X;
		derivs ID, Sol;
		int n;

		
		int spin = 2, l=2, m = 0;
		
		double eps = 204./100000., 
			   a0_over_lambda = 11.5722,//50,//10,//120000./10000., 
			   ID_r0_over_M = 100, 
			   ID_Gauss_width_lambda = 1;

		int N_gaussian = 150, N_bump = 150;
		
		
		set_parameters(&par, spin, l, m, eps, a0_over_lambda, ID_r0_over_M, ID_Gauss_width_lambda, N_gaussian, N_bump);

		fprintf(par.fout,"**** Solving for spin =%d, ell = %d ****\n\n", spin, par.ell); fflush(par.fout);

		for(n=0; n<par.nmax; n++){
		  	par.ncount = n;
		  	par.t0     = par.Tmin + n*par.DT;
			par.t1     = par.Tmin + (n+1)*par.DT;		     
			par.Del_t  = par.t1 - par.t0;			
			
			allocate_derivs(&ID, 0, par.ntotal_slice - 1);
			if(n==0){
				get_Initial_Data(par, ID);	
				Get_PhysDiffMatrices_2D(&par);
				Get_Diff_Operator_Functions(&par);
			}	     
			else{
				get_Initial_Data_newTimeStep(par, Sol, ID);
		      	free_derivs(&Sol, 0, par.ntotal - 1);
		    }
			
			

		    X = dvector(0, par.Ntotal);
			get_initial_guess(par, ID, X);			
			
			allocate_derivs(&Sol, 0, par.Ntotal);
						
			solve_equations(par, ID, X, Sol);  
			
				     
			// output_cheb(par, Sol);
			output_Sol(par, Sol);
			// output_TimeEvolution(par, ID, Sol);
			// output_cheb_TimeEvolution(par, ID, Sol);
			// output_Source(par);
			// output_Sol_derv_sigma(par, Sol);
						

			free_derivs(&ID, 0, par.ntotal_slice - 1);
			free_dvector(X,  0, par.Ntotal);
		}
		free_derivs(&Sol, 0, par.Ntotal);
		free_parameters(&par);

	

	
	return 0;
}

