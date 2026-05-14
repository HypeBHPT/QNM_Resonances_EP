#include "1+1_Free_Boundaries.h"

//   -----------------------------------------------
int main() 
{
	

	int n_omp = omp_get_max_threads();
	double start_time, final_time;
	start_time = omp_get_wtime();	

	int l, l_min=0, l_max=20;

	// #pragma omp parallel for 
	for(l=l_min; l<=l_max; l++){
		parameters par;

		par.i_omp = omp_get_thread_num();
		printf("Solving ell = %d on thread = %d (%d)\n", l, par.i_omp,n_omp );fflush(stdout);


		double *X;
		derivs ID, Sol;
		int n;

		char fn[500];
		

		sprintf(fn,"output_thread_%d.out", par.i_omp);
		par.fout = fopen(fn, "w");

		int m, m_min = l%2, m_max = l;
		for(m=m_min; m<=m_max; m+=2){

			set_parameters(&par, l, m);	
			fprintf(par.fout,"**** Solving ell = %d, m = %d on thread = %d (%d) ****\n", par.ell, par.m, par.i_omp,n_omp ); fflush(par.fout);

			  for(n=0; n<par.nmax; n++){

			  	par.ncount = n;
			  	par.t0     = par.Tmin + n*par.DT;
				par.t1     = par.Tmin + (n+1)*par.DT;		     
				par.Del_t  = par.t1 - par.t0;			
				
				fprintf(par.fout,"Solving for t0=%lf, t1=%lf - n=%d (%d)\n", par.t0, par.t1, n, par.nmax);

				allocate_derivs(&ID, 0, par.ntotal_slice - 1);
				if(n==0){
					get_Initial_Data(par, ID);	
				}	     
				else{
					get_Initial_Data_newTimeStep(par, Sol, ID);
			      	free_derivs(&Sol, 0, par.ntotal - 1);
			    }

			    X = dvector(0, par.Ntotal);
				get_initial_guess(par, ID, X);
				
				allocate_derivs(&Sol, 0, par.Ntotal);
				
				// solve_equations(par, ID, X, Sol);  		    
					     
				// output_cheb(par, Sol);
				// output_Sol(par, Sol);
				// output_Sol_derv_sigma(par, Sol);
				// output_Fluxes_and_Ft(par, Sol);
				

				free_derivs(&ID, 0, par.ntotal_slice - 1);
				free_dvector(X,  0, par.Ntotal);
			}
			free_derivs(&Sol, 0, par.Ntotal);
			free_parameters(&par);
		}
	}
	final_time = omp_get_wtime();
	printf("number of threads = %d, Total time = %3.15e\n", n_omp, final_time-start_time);
	
	return 0;
}

