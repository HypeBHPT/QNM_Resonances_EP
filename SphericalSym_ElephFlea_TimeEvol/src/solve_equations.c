#include "1+1_Free_Boundaries.h"

void solve_equations(parameters par, derivs ID, double *X, derivs Sol){
	static int n=0;
	
	int newt_iter;
	double start_time, final_time;
	
	fprintf(par.fout, "\nSolving equations for t0 = %lf, t1 = %lf --- n=%d (%d) \n", par.t0, par.t1, n, par.nmax-1);

	start_time = clock(); //Start measuring time
	if(par.Solver_FLAG == 0){

		newt_iter = newton_direct(par, ID, X);
	}
	else if(par.Solver_FLAG == 1)
		newt_iter = newton(par, ID, X);
	else{
		fprintf(stderr, "Error in solve_equations: Method = %d not known\n", par.Solver_FLAG );
		exit(-1);
	}
	final_time = clock(); //Stop measuring time	
	fprintf(par.fout,"Newton exit after %d iterations (time = %3.3e s)\n\n", newt_iter, (final_time-start_time)/CLOCKS_PER_SEC);

	get_solution(par, ID, X, Sol);
	n++;
	return;
}
//--------------------------------------------------------------------
void get_solution(parameters par, derivs ID, double *X, derivs W){
  
  int iDom, iField, i0, i1,
      N0=par.N0 ,Ntotal=par.ntotal-1;
  derivs V;
  allocate_derivs(&V, 0, Ntotal);
  Get_V_From_X(par, X, V);

  	for(i0=0; i0<=N0; i0++) {
  		double x0 = par.grid_points_x0[i0], x0p1 = x0+1.;
  		for(iField=0; iField < nFields1D; iField++){
  			int indx1D       = Index_field1D(par, iField, i0),
						indx1D_slice = Index_field1D_slice(par, iField);

				W.X[indx1D]    = ID.X[indx1D_slice]   + x0p1*V.X[indx1D];
				W.Y[indx1D]    = V.X[indx1D]          + x0p1*V.Y[indx1D];
			}
			for(iDom=0; iDom<nDoms; iDom++){
				int N1=par.N1[iDom];
				for(i1=0; i1<=N1; i1++){
					for(iField=0; iField < nFields; iField++){
						int indx       = Index(par, iDom, iField, i1, i0),
								indx_slice = Index_slice(par, iDom, iField, i1);

						W.X[indx]    = ID.X[indx_slice]   + x0p1*V.X[indx];
						W.X1[indx]   = ID.X1[indx_slice]  + x0p1*V.X1[indx];
						W.X11[indx]  = ID.X11[indx_slice] + x0p1*V.X11[indx];
						W.Y[indx]    = V.X[indx]          + x0p1*V.Y[indx];
						W.Y1[indx]   = V.X1[indx]         + x0p1*V.Y1[indx];
						W.Y11[indx]  = V.X11[indx]        + x0p1*V.Y11[indx];	


						// double tau = get_t_from_x0(par, x0), x1=par.grid_points_x1[iDom][i1], sigma=get_sigma_from_x1(par, iDom, x1);
						// if(sigma==1. && iField==0){
						// 	printf("tau = %lf, sigma = %lf, phi(%d) = %3.15e, dphi_dtau(%d)=%3.15e\n", tau,  sigma, iField, W.X[indx], iField, W.Y[indx]);
						// 	// pause();
						// }

				}
				
			}
		}
	}
  
  free_derivs(&V,  0, Ntotal);
  return;
}