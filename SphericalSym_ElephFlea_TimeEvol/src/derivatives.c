#include "1+1_Free_Boundaries.h"
//--------------------------------------------------------------
void Get_DiffMatrices_2D(parameters *par){
	//ATTENTION: Time Domain is perfomed in auxiliary field V defined via
	// W = InitalData(x1) + (1+x0) * V(x0,x1)
	// Thus, matrices pick a factor (1+x0), and x0_derivative pick extra terms with chain rule
	int N = (*par).Ntotal,  i0, j0, i1, j1, idom, iF, N0, N1;;
	double **D0,  **D1, **D11;

	
	(*par).Id_2D = dmatrix(0,N, 0,N);
	fill0_dmatrix((*par).Id_2D, 0, N, 0, N);
	(*par).D0_2D = dmatrix(0,N, 0,N);
	fill0_dmatrix((*par).D0_2D, 0, N, 0, N);
	(*par).D1_2D = dmatrix(0,N, 0,N);
	fill0_dmatrix((*par).D1_2D, 0, N, 0, N);
	// (*par).D00_2D = dmatrix(0,N, 0,N);
	// fill0_dmatrix((*par).D00_2D, 0, n, 0, n);
	(*par).D11_2D = dmatrix(0,N, 0,N);
	fill0_dmatrix((*par).D11_2D, 0, N, 0, N);
	(*par).D01_2D = dmatrix(0,N, 0,N);
	fill0_dmatrix((*par).D01_2D, 0, N, 0, N);



	N0 = (*par).N0;
	for(idom=0; idom<nDoms; idom++){
		
		D0 = dmatrix(0,N0,0,N0);
		Chebyshev_Differentiation_Matrix(N0, D0, (*par).grid_x0);
		
		N1 = (*par).N1[idom];
		D1 = dmatrix(0,N1,0,N1); D11 = dmatrix(0,N1,0,N1);		
		Chebyshev_First_and_Second_Differentiation_Matrix(N1, D1, D11, (*par).grid_x1[idom]);


		for(iF=0; iF<nFields; iF++){
			for(i0=0; i0<=N0; i0++){
				double x0 = (*par).grid_points_x0[i0], x0p1 = x0+1.;
				for(i1=0; i1<=N1; i1++){
					int II = Index(*par, idom, iF, i1, i0);
					(*par).Id_2D[II][II] = 1.;

					
					
					for(j0=0; j0<=N0;j0++){
						int JJ = Index(*par, idom, iF, i1, j0);

						
						(*par).D0_2D[II][JJ] = (*par).Id_2D[II][JJ] +  x0p1*D0[i0][j0];
						// (*par).D00_2D[II][JJ] = D00[i0][j0];
						
						
						for(j1=0; j1<=N1;j1++){
							int JJ = Index(*par, idom, iF, j1, j0);
							(*par).D01_2D[II][JJ] =  x0p1*D0[i0][j0]*D1[i1][j1];
						}
						
					}
					for(j1=0; j1<=N1;j1++){
						int JJ = Index(*par, idom, iF, j1, i0);
						(*par).D1_2D[II][JJ] = x0p1* D1[i1][j1];
						(*par).D01_2D[II][JJ] +=  D1[i1][j1];
						(*par).D11_2D[II][JJ] = x0p1* D11[i1][j1];
					}
					
					(*par).Id_2D[II][II] *= x0p1;

				}
			}
		}
		free_dmatrix(D0, 0,N0,0,N0);
		free_dmatrix(D1, 0,N1,0,N1); free_dmatrix(D11, 0,N1,0,N1);
	}

	return;
}
//-----------------------------------------------------------
void Allocate_PhysDiffMatrices_2D(parameters *par){
	int N = (*par).Ntotal;
	
	(*par).Dtau_2D = dmatrix(0,N, 0,N);
	(*par).Dsig_2D = dmatrix(0,N, 0,N);
	// (*par).Dtautau_2D = dmatrix(0,N, 0,N);
	(*par).Dsigsig_2D = dmatrix(0,N, 0,N);
	(*par).Dtausig_2D = dmatrix(0,N, 0,N);
	
	return;
}
//-----------------------------------------------------------
void Get_PhysDiffMatrices_2D(parameters *par){
	int N=(*par).Ntotal, iF, i0, j0, i1, j1, idom, N0, N1;
	
	derivs SpecMatrices, PhysMatrices;
	allocate_derivs(&SpecMatrices, 0, 0);
	allocate_derivs(&PhysMatrices, 0, 0);
	

	 	
	fill0_dmatrix((*par).Dtau_2D, 0, N, 0, N);	
	fill0_dmatrix((*par).Dsig_2D, 0, N, 0, N);	
	// fill0_dmatrix((*par).Dtautau_2D, 0, n, 0, n);	
	fill0_dmatrix((*par).Dsigsig_2D, 0, N, 0, N);	
	fill0_dmatrix((*par).Dtausig_2D, 0, N, 0, N);



	N0 = (*par).N0;
	
	for(idom=0; idom<nDoms; idom++){
		N1 = (*par).N1[idom];
	
		for(iF=0; iF<nFields; iF++){
			for(i0=0; i0<=N0; i0++){
				

				for(i1=0; i1<=N1; i1++){
					int II = Index(*par, idom, iF, i1, i0);
					
						for(j0=0; j0<=N0;j0++){
							for(j1=0; j1<=N1;j1++){
								int JJ = Index(*par, idom, iF, j1, j0);
								SpecMatrices.Y[0] = (*par).D0_2D[II][JJ];
								// SpecMatrices.Z = (*par).D00_2D[II][JJ];
								SpecMatrices.Y1[0] = (*par).D01_2D[II][JJ];
								SpecMatrices.X1[0] = (*par).D1_2D[II][JJ];
								SpecMatrices.X11[0] = (*par).D11_2D[II][JJ];
								

								
								scale_derivs_x0x1_to_tau_sigma( (*par), idom, i1, SpecMatrices, PhysMatrices, 1);
								
								(*par).Dtau_2D[II][JJ] = PhysMatrices.Y[0];
								// (*par).Dtautau_2D[II][JJ] = PhysMatrices.d2sigma;
								(*par).Dtausig_2D[II][JJ] = PhysMatrices.Y1[0];
								(*par).Dsig_2D[II][JJ] = PhysMatrices.X1[0];
								(*par).Dsigsig_2D[II][JJ] = PhysMatrices.X11[0];
													

							}
						}
				}
			}
		}
	}


	free_derivs(&SpecMatrices, 0, 0);
	free_derivs(&PhysMatrices, 0, 0);
	return;
}
//---------------------------------------------------------------------
// -------------------------------------------------------------------------------
void free_DiffMatrices_2D(parameters *par){
	int N = (*par).Ntotal;
	free_dmatrix((*par).D0_2D, 0,N,0,N);
	free_dmatrix((*par).D1_2D, 0,N,0,N);
	// free_dmatrix((*par).D00_2D, 0,N,0,N);
	free_dmatrix((*par).D11_2D, 0,N,0,N);
	free_dmatrix((*par).D01_2D, 0,N,0,N);

	free_dmatrix((*par).Id_2D, 0,N,0,N);
	free_dmatrix((*par).Dtau_2D, 0,N,0,N);
	// free_dmatrix((*par).Dtautau_2D, 0,N,0,N);
	free_dmatrix((*par).Dsig_2D, 0,N,0,N);
	free_dmatrix((*par).Dsigsig_2D, 0,N,0,N);
	free_dmatrix((*par).Dtausig_2D, 0,N,0,N);



	return;
}
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
void Derivatives_x1(parameters par, derivs Wslice)
{
	double *p, *d1p, *d2p, *q, *d1q, *d2q;
	int iDom, i1, iField, Nmax=par.Nmax;
		
	p   = dvector(0, Nmax);  q   = dvector(0, Nmax);   
	d1p = dvector(0, Nmax);  d1q = dvector(0, Nmax); 
	d2p = dvector(0, Nmax);  d2q = dvector(0, Nmax);

	for(iDom=0; iDom<nDoms; iDom++){
		int N1=par.N1[iDom];
		for(iField=0; iField < nFields; iField++){
			for(i1=0; i1<=N1; i1++){
					int indx = Index_slice(par, iDom, iField, i1);
					p[i1] = Wslice.X[indx];
					q[i1] = Wslice.Y[indx];
			}
			Chebyshev_Collocations_Derivatives(  p, d1p, par.Dx1[iDom], N1);
			Chebyshev_Collocations_Derivatives(d1p, d2p, par.Dx1[iDom], N1);
			Chebyshev_Collocations_Derivatives(  q, d1q, par.Dx1[iDom], N1);
			Chebyshev_Collocations_Derivatives(d1q, d2q, par.Dx1[iDom], N1);
			for(i1=0; i1<=N1; i1++){
				int indx = Index_slice(par, iDom, iField, i1);
				Wslice.X1[indx]  = d1p[i1];
				Wslice.X11[indx] = d2p[i1];
				Wslice.Y1[indx]  = d1q[i1];
				Wslice.Y11[indx] = d2q[i1];
			}
		}
	}
	free_dvector(  p, 0, Nmax);  free_dvector(  q, 0, Nmax);
	free_dvector(d1p, 0, Nmax);  free_dvector(d1q, 0, Nmax);
	free_dvector(d2p, 0, Nmax);  free_dvector(d2q, 0, Nmax);

	return;
}
//---------------------------------------------------------------
void Derivatives_x0x1(parameters par, derivs V)
{
	double *p, *d1p, *d2p, *q, *d1q, *d2q;
	int iDom, i0, i1, iField, N0=par.N0,Nmax=par.Nmax;
		
	p   = dvector(0, Nmax);  q   = dvector(0, Nmax);   
	d1p = dvector(0, Nmax);  d1q = dvector(0, Nmax);   
	d2p = dvector(0, Nmax);  d2q = dvector(0, Nmax);
	
	// #pragma omp parallel for
	for(iDom=0; iDom<nDoms; iDom++){
		int N1=par.N1[iDom];
		for(iField=0; iField < nFields; iField++){
			
			
			for(i1=0; i1<=N1; i1++){
				for(i0=0; i0<=N0; i0++){
					int indx = Index(par, iDom, iField, i1, i0);
					p[i0] = V.X[indx];
				}
				Chebyshev_Collocations_Derivatives(  p, d1p, par.Dx0, N0);

				for(i0=0; i0<=N0; i0++){
					int indx = Index(par, iDom, iField, i1, i0);
					V.Y[indx] = d1p[i0];					
				}
			}

			for(i0=0; i0<=N0; i0++){
				for(i1=0; i1<=N1; i1++){
					int indx = Index(par, iDom, iField, i1, i0);
					p[i1] = V.X[indx];
					q[i1] = V.Y[indx];
				
				}
				Chebyshev_Collocations_Derivatives(  p, d1p, par.Dx1[iDom], N1);
				Chebyshev_Collocations_Derivatives(d1p, d2p, par.Dx1[iDom], N1);
				Chebyshev_Collocations_Derivatives(  q, d1q, par.Dx1[iDom], N1);
				Chebyshev_Collocations_Derivatives(d1q, d2q, par.Dx1[iDom], N1);
				
				for(i1=0; i1<=N1; i1++){
					int indx = Index(par, iDom, iField, i1, i0);
					V.X1[indx]  = d1p[i1];
					V.X11[indx] = d2p[i1];
					V.Y1[indx]  = d1q[i1];
					V.Y11[indx] = d2q[i1];				
				}
			}
		}
	}

	for(iField=0; iField < nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){ 
			int indx = Index_field1D(par, iField, i0);
			p[i0] = V.X[indx];
		}
		Chebyshev_Collocations_Derivatives(  p, d1p, par.Dx0, N0);
		
		for(i0=0; i0<=N0; i0++){
			int indx = Index_field1D(par, iField, i0);
			V.Y[indx] = d1p[i0];			
		}
	}
	free_dvector(  p, 0, Nmax);  free_dvector(  q, 0, Nmax);  
	free_dvector(d1p, 0, Nmax);  free_dvector(d1q, 0, Nmax);  
	free_dvector(d2p, 0, Nmax);  free_dvector(d2q, 0, Nmax);

	return;
}