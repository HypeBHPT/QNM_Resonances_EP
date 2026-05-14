#include "1+1_Free_Boundaries.h"

//-----------------------------------------------------------------------------
int newton(parameters par, derivs ID, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int ntotal = par.ntotal, Ntotal=ntotal-1, iter=0, j, bicgstab_iter, FLAG=-1;
	double *F, *DX, norm, normres;
	derivs V;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	allocate_derivs(&V, 0, Ntotal);



	F_of_X(par, X, V, ID, F);
	norm = norm2(F, 0, Ntotal);
	// PrintVector(F, 0, Ntotal);
	// exit(-1);

	double tol=Newton_tol;//norm*Newton_tol;
	if (Newton_verb == 1){
		printf("\n\n");
		printf(" Newton Raphson Method: Initial Residual: \t |F| = %6.1e (%6.1e) \n", norm, tol);
		printf(" ------------------------------------------------------------------\n");
	}
	
//  	while(iter < Newton_itmin || (norm > tol  && iter < Newton_itmax)){
  	while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
		if(norm < tol) FLAG++;
		iter += 1;
		fill0_dvector(DX, 0, Ntotal);

		par.bicgstab_tol = bicgstab_decr*norm;
		bicgstab_iter    = bicgstab(par, V, ID, F, DX, &normres);

		for (j = 0; j < ntotal; j++) X[j] -= DX[j];

		F_of_X(par, X, V, ID, F);
		norm = norm2(F, 0, Ntotal);
		
		if (isinf(norm) || isnan(norm) || norm > 1.0e+10){
			printf("\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
			exit(1);
		}
		
		if (Newton_verb == 1){
			printf("\t\t\t\t\t Newton: iter = %1d \t Number of bicgstab-steps = %2d \t |F| = %6.1e (%6.1e) %d \n", 
						 iter, bicgstab_iter, norm, tol, FLAG);			
		}
	}
	
	if(norm > tol)
		printf(" Newton Raphson Method failed to converge to prescribed tolerance.\n");
	

	free_dvector(F,  0, Ntotal);
	free_dvector(DX, 0, Ntotal);
	free_derivs(&V,  0, Ntotal);

	if (Newton_verb == 1) printf("\n\n");
	
	return iter;
}
// -------------------------------------------------------------------------------
int bicgstab(parameters par, derivs V, derivs ID, double *F, double *DX, double *normres)
{	// Biconjugate Gradient Stabilized Method
	int ntotal = par.ntotal, Ntotal=ntotal-1, N0 = par.N0, 
		Ntotal_slice = par.ntotal_slice, iter, n;
	double alpha = 0., beta = 0., rho = 0., rho1 = 1., 
		tol = par.bicgstab_tol, rhotol = 1e-50, omega = 0., omegatol = 1e-50,
		*p, *ph, *q, *r, *rt, *s, *sh, *t;
	cheb_derivs chebV;
	
	p  = dvector(0, Ntotal);    ph = dvector(0, Ntotal);
	q  = dvector(0, Ntotal);    r  = dvector(0, Ntotal);
	s  = dvector(0, Ntotal);    sh = dvector(0, Ntotal);
	rt = dvector(0, Ntotal);    t  = dvector(0, Ntotal);
	allocate_cheb_derivs(&chebV, 0, Ntotal_slice, 0, N0);
	
	get_chebV(par, V, chebV);

	// 	compute initial residual rt = r = p = F - J*DX        
	J_times_DX(par, DX, V, ID, r);

	for (n = 0; n < ntotal; n++) rt[n] = r[n] = p[n] = F[n] - r[n];

	
	*normres = norm2(r, 0, Ntotal);
	if (*normres <= tol) return 0;

	for (iter = 0; iter < bicgstab_itmax; iter++) {
		rho = scalarproduct(rt, r, 0, Ntotal);
		if (fabs(rho) < rhotol) break;
		if (iter > 0){ // compute direction vector p
			beta = (rho/rho1)*(alpha/omega);
			for (n = 0; n < ntotal; n++) 
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		// compute direction adjusting vector ph and scalar alpha :

		PreCond(par, chebV, p, ID, ph);  // solves J*ph=p approximately for ph
		J_times_DX(par, ph, V, ID, q);    // q = J*ph

		alpha = rho/scalarproduct(rt, q, 0, Ntotal);
		for (n = 0; n < ntotal; n++) s[n] = r[n] - alpha * q[n];
		// early check of tolerance:
		*normres = norm2(s, 0, Ntotal);
		if (*normres <= tol) {
			for (n = 0; n < ntotal; n++)
				DX[n] += alpha * ph[n];
			break;
		}
		// compute stabilizer vector sh and scalar omega:

		PreCond(par, chebV, s, ID, sh);  // solves J*sh=s approximately for sh
		J_times_DX(par, sh, V, ID, t);    // t = J*sh
		
		omega = scalarproduct(t, s, 0, Ntotal) / scalarproduct(t, t, 0, Ntotal);
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++) {
			DX[n] += alpha * ph[n] + omega * sh[n];
			r[n]   = s[n] - omega * t[n];
		}
		// are we done? 
		*normres = norm2(r, 0, Ntotal);
		if (bicgstab_verb == 1) {printf("BiCGStab: iter = %2d  norm = %6.1e\r", iter+1, *normres);fflush(0); }
		if (*normres <= tol) break;
		rho1 = rho;
		if (fabs(omega) < omegatol) break;
	}

	free_dvector(p,      0, Ntotal);    free_dvector(ph, 0, Ntotal); 
	free_dvector(q,      0, Ntotal);    free_dvector(r,  0, Ntotal);
	free_dvector(s,      0, Ntotal);    free_dvector(sh, 0, Ntotal); 
	free_dvector(rt,     0, Ntotal);    free_dvector(t,  0, Ntotal);
	free_cheb_derivs(&chebV, 0, Ntotal_slice, 0, N0);

	/* iteration failed */
	if (iter > bicgstab_itmax) return -1;
	
	/* breakdown */
	if (fabs(rho)   < rhotol)   return -10;
	if (fabs(omega) < omegatol) return -11;
	
	/* success! */
	if (bicgstab_verb == 1) {printf("BiCGStab: iter = %2d  norm = %6.1e\n", iter+1, *normres);fflush(0); }
	return iter+1;
}
// // -------------------------------------------------------------------------------
// void get_cheb(parameters par, int i_slice, double *p, double *V, double **chebV)
// {
// 	int i0, N0=par.N0, ntotal_2Dslice = par.ntotal_2Dslice;
	
// 	for(i0=0; i0<=N0; i0++){
// 		int iV = i_slice +(NB-j)*ntotal_2Dslice;
// 		p[j] = V[iV];
// 	}
// 	Chebyshev_Coefficients_Radau_RHS(p, chebV[iV_slice], NB);
// }
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
void get_Chebyshev_Coefficients_derivs(derivs psi, derivs c, int N, char *grid){
	Chebyshev_Coefficients(psi.X,   c.X,   N, grid);
	Chebyshev_Coefficients(psi.X1,  c.X1,  N, grid);
	Chebyshev_Coefficients(psi.X11, c.X11, N, grid);



	Chebyshev_Coefficients(psi.Y,   c.Y,   N, grid);
	Chebyshev_Coefficients(psi.Y1,  c.Y1,  N, grid);
	Chebyshev_Coefficients(psi.Y11, c.Y11, N, grid);
}
// -------------------------------------------------------------------------------
void get_chebV(parameters par, derivs V, cheb_derivs chebV)
{
	derivs p, cp;
	int iDom, iField, i1, i0, N0=par.N0, N1;

	allocate_derivs(&p, 0, N0);
	allocate_derivs(&cp, 0, N0);

	for(iDom=0; iDom<nDoms; iDom++){
		N1 = par.N1[iDom];
		for(iField=0; iField<nFields; iField++){
			for(i1=0; i1<=N1; i1++){
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(par, iDom, iField, i1, i0);
					p.X[i0] = V.X[Idx];
					p.X1[i0] = V.X1[Idx];
					p.X11[i0] = V.X11[Idx];

					p.Y[i0] = V.Y[Idx];
					p.Y1[i0] = V.Y1[Idx];
					p.Y11[i0] = V.Y11[Idx];
				}
				get_Chebyshev_Coefficients_derivs(p, cp, N0, par.grid_x0);

				int Idx_slice=Index_slice(par, iDom, iField, i1);
				for(i0=0; i0<=N0; i0++){
					chebV.X[Idx_slice][i0]   = cp.X[i0];
					chebV.X1[Idx_slice][i0]  = cp.X1[i0];
					chebV.X11[Idx_slice][i0] = cp.X11[i0];

					chebV.Y[Idx_slice][i0]	 = cp.Y[i0];
					chebV.Y1[Idx_slice][i0]	 = cp.Y1[i0];
					chebV.Y11[Idx_slice][i0] = cp.Y11[i0];	
				}
			}
		}
	}

	for(iField=0; iField<nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_field1D(par, iField, i0);
			p.X[i0] = V.X[Idx];
			p.Y[i0] = V.Y[Idx];
		}
		int Idx_slice= Index_field1D(par, iField, i0);;
		get_Chebyshev_Coefficients_derivs(p, cp, N0, par.grid_x0);
		for(i0=0; i0<=N0; i0++){
			chebV.X[Idx_slice][i0]   = cp.X[i0];
			chebV.Y[Idx_slice][i0]	 = cp.Y[i0];
		}
	}



	free_derivs(&p, 0, N0);
	free_derivs(&cp, 0, N0);
	return; 
}

