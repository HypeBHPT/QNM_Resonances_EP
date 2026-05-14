#include "1+1_Free_Boundaries.h"

#define jmax 500
#define imax 500

// -------------------------------------------------------------------------------
void test(parameters par, derivs ID, double *X)
{
	int ntotal = par.ntotal, Ntotal=ntotal-1, NB = par.NB, NA = par.NA[0],
		NVslice_total = par.nVslice + nFields1D - 1, Nmax=par.Nmax,
    iField, i, j;
	double *F, *p, *cheb_p;
	cheb_derivs chebV;
	derivs V;
	FILE *stream;

	F        = dvector(0, Ntotal);
	p        = dvector(0, Nmax);  
	cheb_p   = dvector(0, Nmax);   
	allocate_derivs(&V, 0, Ntotal);
	allocate_cheb_derivs(&chebV, 0, NVslice_total, 0, NB);
	
	Get_V_From_X(par, X, V);
	get_chebV(par, V, chebV);
	

	for(iField=0; iField < nFields; iField++){
		
		if(iField == 0) stream = fopen("U.dat",   "w");
		if(iField == 1) stream = fopen("w.dat",   "w");
		if(iField == 2) stream = fopen("rho.dat", "w");
		
		for(j=0;j<=jmax;j++){
			double B = -1.+2.*j/jmax, 
				Bp1 = B+1., tau = 0.5*Bp1, t = par.t0 + par.Del_t*tau;
			fprintf(stream, "\n# Time=%e\n", t);
			for(i=0; i<=NA; i++){
				int iV_slice = Index_slice(par, 0, iField, i);
				p[i] = ID.X[iV_slice] + Bp1*Clenshaw_Chebyshev(chebV.X[iV_slice], NB, B);
			}
			Chebyshev_Coefficients_Lobatto(p, cheb_p, NA);
			for(i=0;i<=imax;i++){
				double A = -1.+2.*i/imax, s = 0.5*(1.+A);
				fprintf(stream, "%e %e \n", s, Clenshaw_Chebyshev(cheb_p, NA, A));
			}
		}
		fclose(stream);

		if(iField == 0) stream = fopen("cheb_U.dat",   "w");
		if(iField == 1) stream = fopen("cheb_w.dat",   "w");
		if(iField == 2) stream = fopen("cheb_rho.dat", "w");
		
		for(i=0; i<=NA; i++)
			fprintf(stream, "%d %e \n", i, cheb_p[i]);
		fclose(stream);
	}
	stream = fopen("R.dat",   "w");
	for(j=0;j<=jmax;j++){
		double B = -1.+2.*j/jmax, 
			Bp1 = B+1., tau = 0.5*Bp1, t = par.t0 + par.Del_t*tau;
		int iV_slice = par.nVslice + 0;
		fprintf(stream, "%e %e \n", t, ID.X[iV_slice] + Bp1*Clenshaw_Chebyshev(chebV.X[iV_slice], NB, B));
	}
	fclose(stream);

	stream = fopen("cheb_R.dat",   "w");
	for(j=0;j<=NB;j++){
		int iV_slice = par.nVslice + 0;
		fprintf(stream, "%d %e \n", j, chebV.X[iV_slice][j]);
	}
	fclose(stream);

	free_dvector(F,      0, Ntotal);
	free_dvector(p,      0, Nmax);  
	free_dvector(cheb_p, 0, Nmax);
	free_derivs(&V,  0, Ntotal);
	free_cheb_derivs(&chebV, 0, NVslice_total, 0, NB);
}
