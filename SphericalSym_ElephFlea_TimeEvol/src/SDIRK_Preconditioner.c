#include "1+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
void get_chebF(parameters par, double *F, double **chebF)
{
	
	int i0, N0=par.N0 , i1, N1, iField, iDom;
	double *p=dvector(0,N0), *cp=dvector(0,N0);

	for(iDom=0; iDom<nDoms; iDom++){
		N1 = par.N1[iDom];
		for(iField=0; iField<nFields; iField++){
			for(i1=0; i1<=N1; i1++){
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(par, iDom, iField, i1, i0);
					p[i0] = F[Idx];
				}				
				Chebyshev_Coefficients(p, cp, N0, par.grid_x0);

				int Idx_slice=Index_slice(par, iDom, iField, i1);
				copy_dvector(chebF[Idx_slice], cp, 0, N0);
			}
		}
	}

	for(iField=0; iField<nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_field1D(par, iField, i0);
			p[i0] = F[Idx];
		}
		Chebyshev_Coefficients(p, cp, N0, par.grid_x0);

		int Idx_slice=Index_field1D_slice(par, iField);
		copy_dvector(chebF[Idx_slice], cp, 0, N0);
	}


	free_dvector(p, 0, N0);
	free_dvector(cp, 0, N0); 

	return;
}
// -------------------------------------------------------------------------------
void PreCond(parameters par, cheb_derivs chebV, double *F, derivs ID, double *X)
{
	
	int i0, ntotal_slice = par.ntotal_slice, Ntotal_slice=ntotal_slice-1, 
		N0=par.N0, iDom, iField, i1, N1;

	double *X_previous = dvector(0,Ntotal_slice), *X_next = dvector(0,Ntotal_slice),
			**chebF = dmatrix(0, Ntotal_slice, 0, N0);
	
		
	fill0_dvector(X_previous,   0, Ntotal_slice);
	fill0_dvector(X_next,   0, Ntotal_slice);
	get_chebF(par, F, chebF);

	for(i0=N0; i0>=0; i0--){// Indetx order of vector in time direction is contrary to time flow
		double x0_previous = i0 == N0? -1.: par.grid_points_x0[i0+1],
			   x0_next=par.grid_points_x0[i0],
			   h = x0_next - x0_previous;
		SDIRK_Step(par, ID, chebV, chebF, x0_previous, h, X_previous, X_next);

		for(iDom=0; iDom<nDoms; iDom++){
			N1 = par.N1[iDom];
			for(iField=0; iField<nFields; iField++){
				for(i1=0; i1<=N1; i1++){
					int i_slice = Index_slice(par, iDom, iField, i1),
						idx = Index(par, iDom, iField, i1, i0);
					X[idx] = X_next[i_slice]/(1.+x0_next);
				}
			}
		}

	   	copy_dvector(X_previous, X_next, 0, Ntotal_slice);
	}

	free_dvector(X_previous,0,Ntotal_slice); free_dvector(X_next,0,Ntotal_slice);
	free_dmatrix(chebF, 0, Ntotal_slice, 0, N0);
	return;
}//-----------------------------------------------------------------------------
void SDIRK_Step(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double h, double *X_previous, double *X_next){
	int i_RK, s = par.SDIRK_s, i_slice, Ntotal_slice=par.ntotal_slice-1, m;
	double **k  = dmatrix(1, s, 0, Ntotal_slice), *Mj=dvector(0, Ntotal_slice);

	par.gah = par.gamma*h;

	fill0_dvector(X_next, 0, Ntotal_slice);	
	copy_dvector(Mj, X_previous, 0, Ntotal_slice);

	if (Newton_LinSolve_SDIRK_verb == 1)  printf("\n\t\t\t\t Time t=%lf, (x0 = %lf) \n", get_t_from_x0(par, x0), x0);
	for(i_RK=1; i_RK<=s;i_RK++){	
		double x0_RK = x0 + h*par.SDIRK_C[i_RK];
		if (Newton_LinSolve_SDIRK_verb == 1) printf("\n\t\t\t\t Runge Kutta Step = %d\n", i_RK);

		LinSolve_SDIRK(par, ID, chebV, chebF, x0_RK, Mj, k[i_RK]);

		fill0_dvector(Mj, 0, Ntotal_slice);
		for(i_slice=0; i_slice<= Ntotal_slice; i_slice++){
			for(m=1; m<=i_RK; m++){
				if(i_RK<s) Mj[i_slice]     += par.SDIRK_A[i_RK+1][m]*k[m][i_slice];
				else       X_next[i_slice] += par.SDIRK_B[m]     *k[m][i_slice];
			}
			if(i_RK<s)   Mj[i_slice]      = X_previous[i_slice] + h*Mj[i_slice];
			else      X_next[i_slice]  = X_previous[i_slice] + h*X_next[i_slice];
		}
	}
	if (Newton_LinSolve_SDIRK_verb == 1) pause();

	free_dvector(Mj, 0, Ntotal_slice);
	free_dmatrix(k,  1, s, 0, Ntotal_slice);
	return;
}
//-------------------------------------------------------------------------
void LinSolve_SDIRK(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj){
	int	nslice_total=par.ntotal_slice, Nslice_total = nslice_total-1, 
		*indx, j, l;
	double *F0, **K, d, *F0_m_Kkj, norm=0., tol=0.;

	F0       = dvector(0, Nslice_total);             // F(qj) = F0 - K * qj
	F0_m_Kkj = dvector(0, Nslice_total);
	indx     = ivector(0, Nslice_total);
	K        = dmatrix(0, Nslice_total, 0, Nslice_total); 


	fill0_dvector(kj, 0, Nslice_total);
	F_of_kj(par, ID, chebV, chebF, x0, Mj, kj, F0);

	if (Newton_LinSolve_SDIRK_verb == 1){
		norm = norm2(F0, 0, Nslice_total);
		tol=norm*Newton_LinSolve_SDIRK_tol;
		printf("\t\t\t\t Initial Residual: \t |F| = %e (%e)\n", norm,tol );
		printf("\t\t\t\t  ------------------------------------------------------------------\n");
	}


	for(j=0; j <= Nslice_total; j++){
		kj[j] = 1.;
		F_of_kj(par, ID, chebV, chebF, x0, Mj, kj, F0_m_Kkj);
		for(l=0; l <= Nslice_total; l++)
			K[l][j] = F0[l] - F0_m_Kkj[l]; 
		kj[j] = 0.;
	}
	
	//********ELIPTICAL LU Solver*******
	ludcmp(K, Nslice_total, indx, &d, 0);
	copy_dvector(kj, F0, 0, Nslice_total);
	lubksb(K, Nslice_total, indx, kj, 0);
	//---***************************----
	

	if (Newton_LinSolve_SDIRK_verb == 1){
		F_of_kj(par, ID, chebV, chebF, x0, Mj, kj, F0);
		norm = norm2(F0, 0, Nslice_total);
		printf("\t\t\t\t LU linear solver: |F| = %e (%3.2e) \n", norm, tol);
	}

	free_dvector(F0,     0, Nslice_total);
	free_ivector(indx,   0, Nslice_total);
	free_dmatrix(K,      0, Nslice_total, 0, Nslice_total);


	return;
}
// -------------------------------------------------------------------------------
void Get_Wslice_W1D_Fslice(parameters par, double x0, derivs ID, double **chebF, cheb_derivs chebV, derivs Wslice, derivs W1D, double *Fslice)
{
int i_slice, iField, ntotal_2Dslice=par.ntotal_2Dslice, N0=par.N0;
	double x0p1 = x0+1.;
	
	for(i_slice=0; i_slice < ntotal_2Dslice; i_slice ++){
		double 
			VX    = Clenshaw_Chebyshev(chebV.X[i_slice],   N0, x0),
			VX1   = Clenshaw_Chebyshev(chebV.X1[i_slice],  N0, x0),
			VX11  = Clenshaw_Chebyshev(chebV.X11[i_slice], N0, x0),
			VY    = Clenshaw_Chebyshev(chebV.Y[i_slice],   N0, x0),
			VY1   = Clenshaw_Chebyshev(chebV.Y1[i_slice],  N0, x0),
			VY11  = Clenshaw_Chebyshev(chebV.Y11[i_slice], N0, x0);
			
		Wslice.X[i_slice]    = ID.X[i_slice]   + x0p1*VX;
		Wslice.X1[i_slice]   = ID.X1[i_slice]  + x0p1*VX1;
		Wslice.X11[i_slice]  = ID.X11[i_slice] + x0p1*VX11;
		
		Wslice.Y[i_slice]    = VX                  + x0p1*VY;
		Wslice.Y1[i_slice]   = VX1                 + x0p1*VY1;
		Wslice.Y11[i_slice]  = VX11                + x0p1*VY11;


		
		Fslice[i_slice] = Clenshaw_Chebyshev(chebF[i_slice], N0, x0);
		
	}
	for(iField=0; iField < nFields1D; iField++){
		i_slice = Index_field1D_slice(par, iField);
		double 
			VX   = Clenshaw_Chebyshev(chebV.X[i_slice],   N0, x0),
			VY   = Clenshaw_Chebyshev(chebV.Y[i_slice],   N0, x0);
			
		W1D.X[iField]    = ID.X[i_slice] + x0p1*VX;
		W1D.Y[iField]    = VX            + x0p1*VY;

		Fslice[i_slice] = Clenshaw_Chebyshev(chebF[i_slice], N0, x0);
	}

	return;
}
// -------------------------------------------------------------------------------
void F_of_kj(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F){
	int iDom, Ntotal_slice = par.ntotal_slice -1, Ntotal_2Dslice = par.ntotal_2Dslice -1, iField, i1, N1, NFields = nFields - 1, NFields1D = nFields1D - 1;
	double *F_2D_Equations, *F_1D_Equations, *Fslice;
	derivs W_grid, Wslice, W1D, DW_grid, DWslice, DW1D;
	
	allocate_derivs(&Wslice, 0, Ntotal_2Dslice); 
	allocate_derivs(&DWslice, 0, Ntotal_2Dslice);
	allocate_derivs(&W_grid, 0, NFields); 
	allocate_derivs(&DW_grid, 0, NFields);
	
	F_2D_Equations = dvector(0, NFields);
	Fslice = dvector(0, Ntotal_slice);

	if(nFields1D!=0){
		F_1D_Equations = dvector(0, NFields1D);
		allocate_derivs(&W1D,    0, NFields1D);	allocate_derivs(&DW1D,    0, NFields1D);		
	}
	else{
		F_1D_Equations = NULL;
		allocate_derivs_NULL(&W1D); allocate_derivs_NULL(&DW1D);
	}

	Get_Wslice_W1D_Fslice(par, x0, ID, chebF, chebV, Wslice, W1D, Fslice);

	get_Ws_from_M_and_k(par, Mj, kj, DWslice, DW1D);

	for(iDom=0; iDom<nDoms; iDom++){
		N1=par.N1[iDom];
		for(i1=0; i1<=N1; i1++){
			Get_Fields2D_Slice_atGrid(par, iDom, i1, Wslice, W_grid);
			Get_Fields2D_Slice_atGrid(par, iDom, i1, DWslice, DW_grid);

			if(i1==N1 && iDom!=0){ // Transition_Condition_Left: Continuity on Function
				derivs W_grid_prevDom, DW_grid_prevDom;
				allocate_derivs(&W_grid_prevDom, 0, NFields); 
				allocate_derivs(&DW_grid_prevDom, 0, NFields);
				Get_Fields2D_Slice_atGrid( par, iDom-1, 0, Wslice, W_grid_prevDom);
				Get_Fields2D_Slice_atGrid( par, iDom-1, 0, DWslice, DW_grid_prevDom);
				for(iField=0; iField < nFields; iField++){
					F_2D_Equations[iField] = DW_grid.X[iField] - DW_grid_prevDom.X[iField];
				}
				free_derivs(&W_grid_prevDom, 0, NFields); 
				free_derivs(&DW_grid_prevDom, 0, NFields);
			}
			else if(i1==0 && iDom!=nDoms-1){// Transition_Condition_Right: Continuity on Function's derivative
				derivs W_grid_NextDom, W_grid_sigma, W_grid_NextDom_sigma,
					  DW_grid_NextDom, DW_grid_sigma, DW_grid_NextDom_sigma;
		
				allocate_derivs(&W_grid_sigma, 0, NFields); 
				allocate_derivs(&DW_grid_sigma, 0, NFields);

				scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, W_grid, W_grid_sigma, nFields);
				scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, DW_grid, DW_grid_sigma, nFields);

				allocate_derivs(&W_grid_NextDom, 0, NFields); 
				allocate_derivs(&DW_grid_NextDom, 0, NFields);
				allocate_derivs(&W_grid_NextDom_sigma, 0, NFields);	
				allocate_derivs(&DW_grid_NextDom_sigma, 0, NFields);					

				Get_Fields2D_Slice_atGrid( par, iDom+1, par.N1[iDom+1], Wslice, W_grid_NextDom);
				Get_Fields2D_Slice_atGrid( par, iDom+1, par.N1[iDom+1], DWslice, DW_grid_NextDom);			
				scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], W_grid_NextDom, W_grid_NextDom_sigma, nFields);
				scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], DW_grid_NextDom, DW_grid_NextDom_sigma, nFields);
				
				for(iField=0; iField < nFields; iField++){
					F_2D_Equations[iField] = DW_grid_NextDom_sigma.X1[iField] - DW_grid_sigma.X1[iField];
				}

				free_derivs(&W_grid_sigma, 0, NFields); free_derivs(&DW_grid_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom_sigma, 0, NFields); free_derivs(&DW_grid_NextDom_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom, 0, NFields); free_derivs(&DW_grid_NextDom, 0, NFields);				
			}
			else if(i1==N1 && iDom==0)
				LinearBoundary_Data_Left(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, F_2D_Equations);
			else if(i1==0 && iDom==nDoms-1)
				LinearBoundary_Data_Right(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, F_2D_Equations);
			else
				JDX_of_XY(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, F_2D_Equations);

			for(iField=0; iField < nFields; iField++){
				int indx = Index_slice(par, iDom, iField, i1);
				F[indx] = F_2D_Equations[iField] - Fslice[indx];
			}
		}
	}

	if(nFields1D!=0){
		JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, F_1D_Equations);
		for(iField=0; iField < nFields1D; iField++){
			int indx = Index_field1D_slice(par, iField);
			F[indx] = F_1D_Equations[iField] - Fslice[indx];
		}
	}

	free_dvector(F_2D_Equations,   0, NFields);
	free_dvector(Fslice,   0, Ntotal_slice);
	free_derivs(&W_grid,      0, NFields); free_derivs(&DW_grid,      0, NFields);
	free_derivs(&Wslice, 0, Ntotal_2Dslice); free_derivs(&DWslice, 0, Ntotal_2Dslice);
	if(nFields1D!=0){
		free_dvector(F_1D_Equations,   0, NFields1D);
		free_derivs(&W1D,    0, NFields1D); free_derivs(&DW1D,    0, NFields1D);
	}
	return;
}
// -------------------------------------------------------------------------------
