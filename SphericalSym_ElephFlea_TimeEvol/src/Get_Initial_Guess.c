#include "1+1_Free_Boundaries.h"

void get_initial_guess(parameters par, derivs ID, double *X){

	if(par.IG_FLAG == 0)
		fill0_dvector(X,0, par.ntotal-1);
	else if(par.IG_FLAG == 1){
			get_Initial_Guess_SDIRK(par, ID, X);
		}
	else{
		fprintf(stderr, "Error in get_initial_guess: IG_FLAG = %d not implemented\n", par.IG_FLAG);
		exit(-1);
	}

	return;
}
//---------------------------------------------------------------------------------
void get_Initial_Guess_SDIRK(parameters par, derivs ID, double *X)
{
	
	int i0, i_slice, ntotal_slice = par.ntotal_slice, Ntotal_slice=ntotal_slice-1, 
		N0=par.N0, iDom, iField, i1, N1;

	double *X_previous = dvector(0,Ntotal_slice), *X_next = dvector(0,Ntotal_slice);

	for(i_slice=0; i_slice < ntotal_slice; i_slice++) X_previous[i_slice] = ID.X[i_slice];

	for(i0=N0; i0>=0; i0--){// Indetx order of vector in time direction is contrary to time flow
		double x0_previous = i0 == N0? -1.: par.grid_points_x0[i0+1],
			   x0_next=par.grid_points_x0[i0],
			   h = x0_next - x0_previous;

		SDIRK_Step_IG(par, x0_previous, h, X_previous, X_next);
		
		for(iDom=0; iDom<nDoms; iDom++){
			N1 = par.N1[iDom];
			for(iField=0; iField<nFields; iField++){
				for(i1=0; i1<=N1; i1++){
					int i_slice = Index_slice(par, iDom, iField, i1),
						idx = Index(par, iDom, iField, i1, i0);
					X[idx] = (X_next[i_slice] - ID.X[i_slice])/(1.+x0_next);
				}
			}
		}

		copy_dvector(X_previous, X_next, 0, Ntotal_slice);

	}
	free_dvector(X_previous,0,Ntotal_slice); free_dvector(X_next,0,Ntotal_slice);
	return;
}
// -------------------------------------------------------------------------------
void SDIRK_Step_IG(parameters par, double x0, double h, double *X_previous, double *X_next)
{
	int i_RK, s = par.SDIRK_s, i_slice, Ntotal_slice=par.ntotal_slice-1, m;
	double **k  = dmatrix(1, s, 0, Ntotal_slice), *Mj=dvector(0, Ntotal_slice);

	par.gah = par.gamma*h;

	fill0_dvector(X_next, 0, Ntotal_slice);	
	copy_dvector(Mj, X_previous, 0, Ntotal_slice);

	
	for(i_RK=1; i_RK<=s;i_RK++){
		double x0_RK = x0 + h*par.SDIRK_C[i_RK];
		if (Newton_SDIRK_verb == 1) printf("\n\t Runge Kutta Step = %d\n", i_RK);

		Newton_SDIRK_IG(par, x0_RK, Mj, k[i_RK]);

		fill0_dvector(Mj, 0, Ntotal_slice);

		for(i_slice=0; i_slice<= Ntotal_slice; i_slice++){
			for(m=1; m<=i_RK; m++){
				if(i_RK<s) Mj[i_slice]     += par.SDIRK_A[i_RK+1][m]*k[m][i_slice];
				else       X_next[i_slice] += par.SDIRK_B[m]     *k[m][i_slice];
			}
			if(i_RK<s)   Mj[i_slice]   = X_previous[i_slice] + h*Mj[i_slice];
			else      X_next[i_slice]  = X_previous[i_slice] + h*X_next[i_slice];
		}
	 }

	free_dvector(Mj, 0, Ntotal_slice);
	free_dmatrix(k,  1, s, 0, Ntotal_slice);

	return;
}
// -------------------------------------------------------------------------------
void get_Ws_from_M_and_k(parameters par, double *Mj, double *kj, derivs Wslice, derivs W1D){
	int islice, nslice=par.ntotal_slice;

	for(islice=0; islice<nslice; islice ++){
		Wslice.Y[islice] = kj[islice];
		Wslice.X[islice] = Mj[islice] + par.gah*kj[islice];

	}
	Derivatives_x1(par, Wslice);

	return;
}
// -------------------------------------------------------------------------------
void Get_Fields2D_Slice_atGrid(parameters par, int iDom, int i1, derivs Wslice, derivs W){
	int iField;

	for(iField=0; iField < nFields; iField++){
		int indx_slice = Index_slice(par, iDom, iField, i1);
		
		W.X[iField]    = Wslice.X[indx_slice];  
		W.X1[iField]   = Wslice.X1[indx_slice]; 
		W.X11[iField]  = Wslice.X11[indx_slice];

		W.Y[iField]    = Wslice.Y[indx_slice];  
		W.Y1[iField]   = Wslice.Y1[indx_slice]; 
		W.Y11[iField]  = Wslice.Y11[indx_slice];
	}

	return;
}
// -------------------------------------------------------------------------------
void F_of_kj_IG(parameters par, double x0, double *Mj, double *kj, double *F){
	int iDom, Ntotal_2Dslice = par.ntotal_2Dslice -1, iField, i1, N1, NFields = nFields - 1, NFields1D = nFields1D - 1;

	double *F_2D_Equations, *F_1D_Equations;
	F_2D_Equations = dvector(0, NFields);

	derivs W_grid, Wslice, W1D;
	allocate_derivs(&Wslice, 0, Ntotal_2Dslice);
	allocate_derivs(&W_grid, 0, NFields);

	if(nFields1D!=0){
		F_1D_Equations = dvector(0, NFields1D);
		allocate_derivs(&W1D,    0, NFields1D);		
	}
	else{
		F_1D_Equations = NULL;
		allocate_derivs_NULL(&W1D);
	}

	get_Ws_from_M_and_k(par, Mj, kj, Wslice, W1D); 

	for(iDom=0; iDom<nDoms; iDom++){
		N1=par.N1[iDom];
		for(i1=0; i1<=N1; i1++){
			Get_Fields2D_Slice_atGrid(par, iDom, i1, Wslice, W_grid);

			if(i1==N1 && iDom!=0){ // Transition_Condition_Left: Continuity on Function
				derivs W_grid_prevDom;
				allocate_derivs(&W_grid_prevDom, 0, NFields);
				Get_Fields2D_Slice_atGrid( par, iDom-1, 0, Wslice, W_grid_prevDom);
				for(iField=0; iField < nFields; iField++){
					F_2D_Equations[iField] = W_grid.X[iField] - W_grid_prevDom.X[iField] - jump_field(par, iField, x0);
				}
				free_derivs(&W_grid_prevDom, 0, NFields);
			}
			else if(i1==0 && iDom!=nDoms-1){// Transition_Condition_Right: Continuity on Function's derivative
				derivs W_grid_NextDom, W_grid_sigma, W_grid_NextDom_sigma;
		
				allocate_derivs(&W_grid_sigma, 0, NFields);
				scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, W_grid, W_grid_sigma, nFields);

				allocate_derivs(&W_grid_NextDom, 0, NFields);
				allocate_derivs(&W_grid_NextDom_sigma, 0, NFields);					

				Get_Fields2D_Slice_atGrid( par, iDom+1, par.N1[iDom+1], Wslice, W_grid_NextDom);			
				scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], W_grid_NextDom, W_grid_NextDom_sigma, nFields);
				
				for(iField=0; iField < nFields; iField++){
					F_2D_Equations[iField] = W_grid_NextDom_sigma.X1[iField] - W_grid_sigma.X1[iField] - jump_derivative(par, iField, x0);
				}

				free_derivs(&W_grid_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom, 0, NFields);
					
			}
			else if(i1==N1 && iDom==0)    	Boundary_Data_Left(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);
			else if(i1==0 && iDom==nDoms-1)  Boundary_Data_Right(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);
			else 							F_of_XY(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);

			for(iField=0; iField < nFields; iField++){
				int indx = Index_slice(par, iDom, iField, i1);	
				F[indx] = F_2D_Equations[iField];
			}
		}
			
	}
	if(nFields1D!=0){
		FBound_of_XY(par, x0, Wslice, W1D, F_1D_Equations);
		for(iField=0; iField < nFields1D; iField++){
			int indx = Index_field1D_slice(par, iField);
			F[indx] = F_1D_Equations[iField];
		}
	}

	free_dvector(F_2D_Equations,   0, NFields);
	free_derivs(&W_grid,      0, NFields);
	free_derivs(&Wslice, 0, Ntotal_2Dslice);
	if(nFields1D!=0){
		free_dvector(F_1D_Equations,   0, NFields1D);
		free_derivs(&W1D,    0, NFields1D);
	}

	return;
}
// -------------------------------------------------------------------------------
void J_times_Dkj_IG(parameters par, double x0, double *Mj, double *kj, double *Dkj, double *JDqj){
	int iDom, Ntotal_2Dslice = par.ntotal_2Dslice -1, Ntotal_slice = par.ntotal_slice - 1, iField, i1, N1, NFields = nFields - 1, NFields1D = nFields1D - 1;

	double *JDX_2D_Equations, *JDX_1D_Equations, *Null_M;
	JDX_2D_Equations = dvector(0, NFields);

	derivs W_grid, DW_grid, 
		   Wslice, DWslice, 
		   W1D, DW1D;
	allocate_derivs(&Wslice, 0, Ntotal_2Dslice); 
	allocate_derivs(&DWslice, 0, Ntotal_2Dslice);
	allocate_derivs(&W_grid, 0, NFields); 
	allocate_derivs(&DW_grid, 0, NFields);

	if(nFields1D!=0){
		JDX_1D_Equations = dvector(0, NFields1D);
		allocate_derivs(&W1D,    0, NFields1D);	
		allocate_derivs(&DW1D,    0, NFields1D);		
	}
	else{
		JDX_1D_Equations = NULL;
		allocate_derivs_NULL(&W1D);
		allocate_derivs_NULL(&DW1D);
	}

	get_Ws_from_M_and_k(par, Mj, kj, Wslice, W1D); 

	Null_M = dvector(0, Ntotal_slice);
	fill0_dvector(Null_M, 0, Ntotal_slice);

	get_Ws_from_M_and_k(par, Null_M, Dkj, DWslice, DW1D); 

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
					JDX_2D_Equations[iField] = DW_grid.X[iField] - DW_grid_prevDom.X[iField];
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
					JDX_2D_Equations[iField] = DW_grid_NextDom_sigma.X1[iField] - DW_grid_sigma.X1[iField];
				}

				free_derivs(&W_grid_sigma, 0, NFields); free_derivs(&DW_grid_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom_sigma, 0, NFields); free_derivs(&DW_grid_NextDom_sigma, 0, NFields);
				free_derivs(&W_grid_NextDom, 0, NFields); free_derivs(&DW_grid_NextDom, 0, NFields);
					
			}
			else if(i1==N1 && iDom==0)    	LinearBoundary_Data_Left(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);
			else if(i1==0 && iDom==nDoms-1)  LinearBoundary_Data_Right(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);
			else 							JDX_of_XY(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);

			for(iField=0; iField < nFields; iField++){
				int indx = Index_slice(par, iDom, iField, i1);	
				JDqj[indx] = JDX_2D_Equations[iField];
			}
		}
			
	}
	if(nFields1D!=0){
		JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, JDX_1D_Equations);
		for(iField=0; iField < nFields1D; iField++){
			int indx = Index_field1D_slice(par, iField);
			JDqj[indx] = JDX_1D_Equations[iField];
		}
	}

	free_dvector(JDX_2D_Equations,   0, NFields);
	free_derivs(&W_grid,      0, NFields); free_derivs(&DW_grid,      0, NFields);
	free_derivs(&Wslice, 0, Ntotal_2Dslice); free_derivs(&DWslice, 0, Ntotal_2Dslice);
	free_dvector(Null_M,   0, Ntotal_slice);
	if(nFields1D!=0){
		free_dvector(JDX_1D_Equations,   0, NFields1D);
		free_derivs(&W1D,    0, NFields1D); free_derivs(&DW1D,    0, NFields1D);
	}

	return;
}
// // -------------------------------------------------------------------------------
