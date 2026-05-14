 #include "1+1_Free_Boundaries.h"

//-------------------------------------------------------------------------------
void Get_V_From_X(parameters par, double *X, derivs V){
  int i, Ntotal=par.Ntotal;
  
//   #pragma omp parallel for
  for(i=0; i<=Ntotal; i++){
    V.X[i]=X[i];
  }

  Derivatives_x0x1(par, V);
  return;
}
// -------------------------------------------------------------------------------
void Get_Fields1D(parameters par, int i0, derivs ID, derivs V, derivs W1D){
	int iField;	
	double x0 = par.grid_points_x0[i0], x0p1=x0+1.;

	for(iField=0; iField < nFields1D; iField++){
		int 
			indx       = Index_field1D(par, iField, i0),
			indx_slice = Index_field1D_slice(par, iField);
			
		W1D.X[iField]    = ID.X[indx_slice]   + x0p1*V.X[indx];
		W1D.Y[iField]    = V.X[indx]          + x0p1*V.Y[indx];

	}

	return;
}
// -------------------------------------------------------------------------------
void Get_Fields2D_atGrid(parameters par, int iDom, int i0, int i1, derivs ID, derivs V, derivs W_grid){
	int iField;
	double x0 = par.grid_points_x0[i0], x0p1=x0+1.;

	for(iField=0; iField < nFields; iField++){
		int indx       = Index(par, iDom, iField, i1, i0), 
			indx_slice = Index_slice(par, iDom, iField, i1);

		W_grid.X[iField]    = ID.X[indx_slice]   + x0p1*V.X[indx];
		W_grid.X1[iField]   = ID.X1[indx_slice]  + x0p1*V.X1[indx];
		W_grid.X11[iField]  = ID.X11[indx_slice] + x0p1*V.X11[indx];
		W_grid.Y[iField]    = V.X[indx]          + x0p1*V.Y[indx];
		W_grid.Y1[iField]   = V.X1[indx]         + x0p1*V.Y1[indx];
		W_grid.Y11[iField]  = V.X11[indx]        + x0p1*V.Y11[indx];

	}

	return;
}
// -------------------------------------------------------------------------------
void Get_Fields2D_atSlice(parameters par, int i0, int i1, derivs ID, derivs V, derivs Wslice){
	int iDom, iField;
	double x0 = par.grid_points_x0[i0], x0p1=x0+1.;

	for(iDom=0; iDom<nDoms; iDom++){		
		for(iField=0; iField < nFields; iField++){
			int indx       = Index(par, iDom, iField, i1, i0), 
				indx_slice = Index_slice(par, iDom, iField, i1);

			Wslice.X[indx_slice]    = ID.X[indx_slice]   + x0p1*V.X[indx];
			Wslice.X1[indx_slice]   = ID.X1[indx_slice]  + x0p1*V.X1[indx];
			Wslice.X11[indx_slice]  = ID.X11[indx_slice] + x0p1*V.X11[indx];
			
			Wslice.Y[indx_slice]    = V.X[indx]          + x0p1*V.Y[indx];
			Wslice.Y1[indx_slice]   = V.X1[indx]         + x0p1*V.Y1[indx];
			Wslice.Y11[indx_slice]  = V.X11[indx]        + x0p1*V.Y11[indx];
		}		
	}
	
	return;
}
// -------------------------------------------------------------------------------
void F_of_X(parameters par, double *X, derivs V, derivs ID, double *F){
	int i1, N1, iDom, i0, N0 = par.N0, iField, NFields1D = nFields1D-1, NFields = nFields-1, Ntotal_2Dslice = par.ntotal_2Dslice -1;
	double *F_2D_Equations, *F_1D_Equations;

	derivs W_grid, W1D, Wslice;
	F_2D_Equations = dvector(0, NFields);
	
	allocate_derivs(&W_grid, 0, NFields);
	allocate_derivs(&Wslice, 0, Ntotal_2Dslice);	 

	Get_V_From_X(par, X, V);

	if(nFields1D!=0){
		F_1D_Equations = dvector(0, NFields1D);
		allocate_derivs(&W1D,    0, NFields1D);		
	}
	else{
		F_1D_Equations = NULL;
		allocate_derivs_NULL(&W1D);
	} 

	for(i0=0; i0<=N0; i0++){
		Get_Fields1D(par, i0, ID, V, W1D);
		double x0 = par.grid_points_x0[i0];

		for(iDom=0; iDom<nDoms; iDom++){
			N1=par.N1[iDom];
			for(i1=0; i1<=N1; i1++){
				Get_Fields2D_atSlice(par, i0, i1, ID,  V, Wslice);
				Get_Fields2D_atGrid( par, iDom, i0, i1, ID,  V, W_grid);

				if(i1==N1 && iDom!=0){ // Transition_Condition_Left: Continuity on Function
					derivs W_grid_prevDom;
					allocate_derivs(&W_grid_prevDom, 0, NFields);
					Get_Fields2D_atGrid( par, iDom-1, i0, 0, ID,  V, W_grid_prevDom);
					for(iField=0; iField < nFields; iField++){
						F_2D_Equations[iField] = W_grid.X[iField] - W_grid_prevDom.X[iField] - jump_field(par, iField, x0);
					// 	//-----
					// double tau = get_t_from_x0(par, x0), x1=par.grid_points_x1[iDom][i1], sigma=get_sigma_from_x1(par, iDom, x1);
					// printf("iF = %d idom = %d, tau=%3.15e\t sigma=%3.15e\t \nphi_dom1=%3.15e\t phi_dom2=%3.15e\t dphi = %3.15e, jump = %3.15e\t F[%d]=%3.5e\n", iField, iDom, tau, sigma, W_grid_prevDom.X[iField], W_grid.X[iField], W_grid.X[iField] - W_grid_prevDom.X[iField], jump_field(par, iField, x0),
					// 	iField,F_2D_Equations[iField]);
					// //--------
					}

			

					free_derivs(&W_grid_prevDom, 0, NFields);
				}
				else if(i1==0 && iDom!=nDoms-1){// Transition_Condition_Right: Continuity on Function's derivative
					derivs W_grid_NextDom, W_grid_sigma, W_grid_NextDom_sigma;
		
					allocate_derivs(&W_grid_sigma, 0, NFields);
					scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, W_grid, W_grid_sigma, nFields);

					allocate_derivs(&W_grid_NextDom, 0, NFields);
					allocate_derivs(&W_grid_NextDom_sigma, 0, NFields);					

					Get_Fields2D_atGrid( par, iDom+1, i0, par.N1[iDom+1], ID,  V, W_grid_NextDom);					
					scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], W_grid_NextDom, W_grid_NextDom_sigma, nFields);
					
					for(iField=0; iField < nFields; iField++){
						F_2D_Equations[iField] = W_grid_NextDom_sigma.X1[iField] - W_grid_sigma.X1[iField] - jump_derivative(par, iField, x0);
					
					// //-----
					// double tau = get_t_from_x0(par, x0), x1=par.grid_points_x1[iDom][i1], sigma=get_sigma_from_x1(par, iDom, x1);
					// printf("iF = %d idom = %d, tau=%3.15e\t sigma=%3.15e\t \nphi_dom1,sigma=%3.15e\t phi_dom2,sigma=%3.15e\t dphi,sigma = %3.15e, jump = %3.15e\t F[%d]=%3.5e\n", iField, iDom, tau, sigma, W_grid_sigma.X1[iField], W_grid_NextDom_sigma.X1[iField], W_grid_NextDom_sigma.X1[iField] - W_grid_sigma.X1[iField], jump_derivative(par, iField, x0),
					// 	iField,F_2D_Equations[iField]);
					// //--------
					}
										

					free_derivs(&W_grid_sigma, 0, NFields);
					free_derivs(&W_grid_NextDom_sigma, 0, NFields);
					free_derivs(&W_grid_NextDom, 0, NFields);
					
				}  
				else if(i1==N1 && iDom==0)    		Boundary_Data_Left(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);
				else if(i1==0 && iDom==nDoms-1)  Boundary_Data_Right(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);
				else 														F_of_XY(par, iDom, x0, i1, W_grid, W1D, F_2D_Equations);

				for(iField=0; iField < nFields; iField++){
					int indx = Index(par, iDom, iField, i1, i0);	
					F[indx] = F_2D_Equations[iField];
				}
			}


			
		}

		if(nFields1D!=0){
			FBound_of_XY(par, i0, Wslice, W1D, F_1D_Equations);
			for(iField=0; iField < nFields1D; iField++){
				int indx = Index_field1D(par, iField, i0);
				F[indx] = F_1D_Equations[iField];
			}
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
void J_times_DX(parameters par, double *DX, derivs V, derivs ID, double *JDX){
	int i1, N1, iDom, i0, N0 = par.N0, iField, Ntotal = par.Ntotal, NFields1D = nFields1D-1, NFields = nFields-1, Ntotal_slice = par.ntotal_slice -1, Ntotal_2Dslice = par.ntotal_2Dslice -1;
	double *JDX_2D_Equations, *JDX_1D_Equations=NULL;

	derivs W_grid, W1D, Wslice, Null_ID,
		   DV, DW_grid, DW1D, DWslice;

	JDX_2D_Equations = dvector(0, NFields);
	
	allocate_derivs(&DV, 0, Ntotal);
	allocate_derivs(&W_grid, 0, NFields); allocate_derivs(&DW_grid, 0, NFields);
	allocate_derivs(&Wslice, 0, Ntotal_2Dslice); allocate_derivs(&DWslice, 0, Ntotal_2Dslice);
	allocate_derivs_NULL(&W1D); allocate_derivs_NULL(&DW1D);	 

	Get_V_From_X(par, DX, DV);
	if(nFields1D!=0){
		JDX_1D_Equations = dvector(0, NFields1D);
		allocate_derivs(&W1D,    0, NFields1D); allocate_derivs(&DW1D,    0, NFields1D);
		
	}

	allocate_derivs(&Null_ID, 0, Ntotal_slice);
	fill0_derivs(Null_ID, 0, Ntotal_slice);


	for(i0=0; i0<=N0; i0++){
		Get_Fields1D(par, i0, ID, V, W1D); Get_Fields1D(par, i0, Null_ID, DV, DW1D);
		double x0 = par.grid_points_x0[i0];

		for(iDom=0; iDom<nDoms; iDom++){
			N1=par.N1[iDom];
			for(i1=0; i1<=N1; i1++){
				Get_Fields2D_atSlice(par, i0, i1, ID,  V, Wslice); 
				Get_Fields2D_atSlice(par, i0, i1, Null_ID,  DV, DWslice);
				Get_Fields2D_atGrid( par, iDom, i0, i1, ID,  V, W_grid); 
				Get_Fields2D_atGrid( par, iDom, i0, i1, Null_ID,  DV, DW_grid);

				
				if(i1==N1 && iDom!=0){ // Transition_Condition_Left: Continuity on Function
					derivs W_grid_prevDom, DW_grid_prevDom;
					allocate_derivs(&W_grid_prevDom, 0, NFields); allocate_derivs(&DW_grid_prevDom, 0, NFields);
					Get_Fields2D_atGrid( par, iDom-1, i0, 0, ID,  V, W_grid_prevDom);
					Get_Fields2D_atGrid( par, iDom-1, i0, 0, Null_ID,  DV, DW_grid_prevDom);

					for(iField=0; iField < nFields; iField++){
						JDX_2D_Equations[iField] = DW_grid.X[iField] - DW_grid_prevDom.X[iField];
					}
					free_derivs(&W_grid_prevDom, 0, NFields); free_derivs(&DW_grid_prevDom, 0, NFields);
				}
				else if(i1==0 && iDom!=nDoms-1){// Transition_Condition_Right: Continuity on Function's derivative
					derivs W_grid_NextDom, DW_grid_NextDom, 
								 W_grid_sigma, DW_grid_sigma, 
								 W_grid_NextDom_sigma, DW_grid_NextDom_sigma;
		
					allocate_derivs(&W_grid_sigma, 0, NFields); allocate_derivs(&DW_grid_sigma, 0, NFields);
					scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, W_grid, W_grid_sigma, nFields);
					scale_derivs_x0x1_to_tau_sigma(par, iDom, i1, DW_grid, DW_grid_sigma, nFields);

					allocate_derivs(&W_grid_NextDom, 0, NFields); allocate_derivs(&DW_grid_NextDom, 0, NFields);
					allocate_derivs(&W_grid_NextDom_sigma, 0, NFields); allocate_derivs(&DW_grid_NextDom_sigma, 0, NFields);					

					Get_Fields2D_atGrid( par, iDom+1, i0, par.N1[iDom+1], ID,  V, W_grid_NextDom);
					Get_Fields2D_atGrid( par, iDom+1, i0, par.N1[iDom+1], Null_ID,  DV, DW_grid_NextDom);					
					scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], W_grid_NextDom, W_grid_NextDom_sigma, nFields);
					scale_derivs_x0x1_to_tau_sigma(par, iDom+1, par.N1[iDom+1], DW_grid_NextDom, DW_grid_NextDom_sigma, nFields);
					
					for(iField=0; iField < nFields; iField++){
						JDX_2D_Equations[iField] = DW_grid_NextDom_sigma.X1[iField] - DW_grid_sigma.X1[iField];
					}

					free_derivs(&W_grid_sigma, 0, NFields); free_derivs(&DW_grid_sigma, 0, NFields);
					free_derivs(&W_grid_NextDom_sigma, 0, NFields); free_derivs(&DW_grid_NextDom_sigma, 0, NFields);
					free_derivs(&W_grid_NextDom, 0, NFields); free_derivs(&DW_grid_NextDom, 0, NFields);
					
				} 
				else if(i1==N1 && iDom==0)    				 LinearBoundary_Data_Left(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);
				else if(i1==0 && iDom==nDoms-1)  LinearBoundary_Data_Right(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);
				else 							
				JDX_of_XY(par, iDom, x0, i1, W_grid, W1D, DW_grid, DW1D, JDX_2D_Equations);
				

				for(iField=0; iField < nFields; iField++){
					int indx = Index(par, iDom, iField, i1, i0);
					JDX[indx] = JDX_2D_Equations[iField];
				}
			}

			
		}

		if(nFields1D!=0){
			JDXBound_of_XY(par, i0, Wslice, W1D, DWslice, DW1D, JDX_1D_Equations);
			for(iField=0; iField < nFields1D; iField++){
				int indx = Index_field1D(par, iField, i0);
				JDX[indx] = JDX_1D_Equations[iField];
			}
		}
	}

	free_dvector(JDX_2D_Equations,   0, NFields);
	free_derivs(&W_grid,      0, NFields); free_derivs(&DW_grid,      0, NFields);
	free_derivs(&Wslice, 0, Ntotal_2Dslice); free_derivs(&DWslice, 0, Ntotal_2Dslice);
	free_derivs(&DV, 0, Ntotal);
	if(nFields1D!=0){
		free_dvector(JDX_1D_Equations,   0, NFields1D);
		free_derivs(&W1D,    0, NFields1D); free_derivs(&DW1D,    0, NFields1D);
	}

	return;
}
// -------------------------------------------------------------------------------	
void Get_LinearOperator(parameters par, double **M){
	int idom, i0, i1, N0, N1, N=par.Ntotal;
	fill0_dmatrix(M, 0, N, 0, N);

	N0 = par.N0;
	for(idom=0;idom<nDoms; idom++){
		N1 = par.N1[idom];
		for(i0=0; i0<=N0; i0++){
			for(i1=0; i1<=N1; i1++){
				if(i1*( i1-N1 )==0 ){//BOUNDARY CONDITION
            		LinearBoundaryConditions_Matrix(par, idom, i0, i1, M);
          		}
				else{//BULK EQUATIONS
					LinearOperator_Matrix(par, idom, i0, i1, M);
				}				
			}
		}
	}


	return;
}
// -------------------------------------------------------------------------------