#include "1+1_Free_Boundaries.h"

double delta(double a, double b){
  
  double out = a==b? 1:0;
  return out;
  
}

double ID_Function(parameters par, double sigma){
	double f, 
	lambda_rStar = par.lambda_rStar, rp = par.ro_over_M,
	delta_rp, x, dx_dsigma, rStar_over_M, rStar_P, x_P, dxP_dsigma, sigma_p, dsigmap_dr;

	func_tortoise_x(par,  sigma, &x, &dx_dsigma);
  	rStar_over_M = par.lambda_over_M * x;


  
   func_sigma_of_r(par, rp, &sigma_p, &dsigmap_dr);
   func_tortoise_x(par,  sigma_p, &x_P, &dxP_dsigma);
   rStar_P = par.lambda_over_M * x_P;

  delta_rp = rStar_over_M - rStar_P;
 if(sigma==0.){
	f=0.;
 }
 else if(sigma==1.){
	f=0.;
 }
 else{
	f = 1./ ( sqrt(2 * Pi  ) * lambda_rStar) * exp(  -0.5 * sqr(  delta_rp/lambda_rStar       )    );
 }

	
// f = 1.;

	return f;
}


// -------------------------------------------------------------------------------
void get_Initial_Data(parameters par, derivs ID)
{ 
	int i1,idom, N1;
		
	
	for(idom = 0; idom<nDoms; idom++){
	  N1=par.N1[idom];

	  for(i1=0; i1<=N1; i1++){
				int Indx_phi_Real  = Index_slice(par, idom, 0, i1),
				    // Indx_phi_Imag  = Index_slice(par, idom, 1,  i1),
				    
					Indx_psi_Real  = Index_slice(par, idom, 1, i1);
				    // Indx_psi_Imag  = Index_slice(par, idom, 3,  i1);

				double sigma, dx_dsigma;
				get_sigma(par, idom, i1, &sigma, &dx_dsigma);
				

				if(par.HomoEq==0){
					
					ID.X[Indx_phi_Real] = ID_Function(par, sigma);
					// ID.X[Indx_phi_Imag] = 0.;
					ID.X[Indx_psi_Real] = 0.;
					// ID.X[Indx_psi_Imag] = 0;

				}
				else if(par.HomoEq==1){
					ID.X[Indx_phi_Real] = 0.;
					// ID.X[Indx_phi_Imag] = 0.;
					ID.X[Indx_psi_Real] = 0.;
					// ID.X[Indx_psi_Imag] = 0;
				}				
		}
	}
	
	Derivatives_x1(par, ID);	
	Output_Initial_Data(par, ID);
	Output_Cheb_Initial_Data(par,ID);
	// exit(-1);
	return;
}
//--------------------------------------------------------------------------------
void get_Initial_Data_newTimeStep(parameters par, derivs W, derivs ID){
  
  int iDom, iField, i0, i1;
  
  
	i0=0;	//Final Time step
	for(iField=0; iField < nFields1D; iField++){
	  
		int indx       = Index_field1D(par, iField, i0),
				indx_slice = Index_field1D_slice(par, iField);
			
		ID.X[indx_slice] = W.X[indx];
	}
	for(iDom=0; iDom<nDoms; iDom++){
		int N1=par.N1[iDom];
		for(i1=0; i1<=N1; i1++){
			for(iField=0; iField < nFields; iField++){
				int 
					indx       = Index(par, iDom, iField, i1, i0), 
					indx_slice = Index_slice(par, iDom, iField, i1);	
				
				ID.X[indx_slice] = W.X[indx];
			}
		}
	}

	Derivatives_x1(par, ID);
	
}
//-----------------
// -------------------------------------------------------------------------------
void Output_Initial_Data(parameters par, derivs ID)
{ 
	int i1,idom, iF, N1;
	char fn_ID[500];
	FILE *fp;
		
	
	for(idom = 0; idom<nDoms; idom++){
	  N1=par.N1[idom];
	  for(iF=0; iF<nFields; iF++){
		sprintf(fn_ID, "data/%s/ID_dom_%d_field_%d.dat", par.SimName, idom, iF);
		fp   =fopen(fn_ID, "w");  

		for(i1=0; i1<=N1; i1++){
			int Indx = Index_slice(par, idom, iF, i1);			

			double  sigma, dx_dsigma;
			get_sigma(par, idom, i1, &sigma, &dx_dsigma);

			double r, dr;
			func_r_of_sigma(par, sigma, &r, &dr);
			
			fprintf(fp, "%3.15e\t %3.15e\t %3.15e\n", sigma, r, ID.X[Indx]);
		}
		fclose(fp);
	  }
	}
		
	return;
}
// -------------------------------------------------------------------------------
void Output_Cheb_Initial_Data(parameters par,derivs ID)
{ 
	int i1, iF, idom, N1;
	char fn_ID[500];
	FILE *fp;
	double *f, *cf;
		
	
	for(idom = 0; idom<nDoms; idom++){
	  N1=par.N1[idom];
	  f=dvector(0,N1);
	  cf=dvector(0,N1);
	  for(iF=0; iF<nFields; iF++){

		for(i1=0; i1<=N1; i1++){
			int Indx = Index_slice(par, idom, iF, i1);	
			f[i1] = ID.X[Indx];	
		}
		Chebyshev_Coefficients(f, cf, N1, par.grid_x1[idom]);

		sprintf(fn_ID, "data/%s/cheb_ID_dom_%d_field_%d.dat", par.SimName, idom, iF);
	  	fp   =fopen(fn_ID, "w"); 
		for(i1=0;i1<=N1;i1++) fprintf(fp, "%d\t %3.15e\n", i1, cf[i1]);
		fclose(fp);
	  }
	}
		
	return;
}
