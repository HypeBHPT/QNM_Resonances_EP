#include "1+1_Free_Boundaries.h"
#define imax 20
#define jmax 5

//-------------------------------------------------------------------------------
void output_cheb(parameters par, derivs Sol){
  
  FILE *fp;
  char fn[500];

  
  int iDom , iField,
      i0, i1, N0=par.N0, N1;
  
  for(iDom=0; iDom<nDoms; iDom++){
    N1=par.N1[iDom];
    for(iField=0; iField<nFields; iField++){
      double *p=dvector(0, N1), *cp=dvector(0, N1), **c=dmatrix(0,N0, 0,N1);
      for(i0=N0; i0>=0; i0--){
      	for(i1=0; i1<=N1; i1++){
      		int indx = Index(par, iDom,iField,i1,i0);
      		p[i1]=Sol.X[indx];
      	}
      	Chebyshev_Coefficients(p, cp, N1, par.grid_x1[iDom]);

      	for(i1=0; i1<=N1; i1++){
      		c[i0][i1] = cp[i1];
      	}
      }
      free_dvector(p,0,N1); free_dvector(cp, 0, N1);
       
      sprintf(fn, "data/%s/cheb_field_%d_dom_%d_s_%03d.dat", par.SimName,iField,iDom, par.ncount);
      fp=fopen(fn,"w");
           
      //Header
      fprintf(fp,"#1: i1\t");
      for(i0=N0; i0>=0; i0--){
      	double t=get_t(par, i0);
      	fprintf(fp,"#%d:c_i(t=%3.5e)\t", N0-i0+2, t);
      }
      fprintf(fp, "\n");
      //end Header
      
      //Data
      for(i1=0; i1<=N1; i1++){
      	fprintf(fp,"%d\t", i1);
      	for(i0=N0; i0>=0; i0--){
      		fprintf(fp,"%3.15e\t", c[i0][i1] );
      	}
      	fprintf(fp,"\n");
      }
      fclose(fp);
      //end Data

      p=dvector(0,N0); cp=dvector(0,N0);

      for(i1=0; i1<=N1; i1++){
      	for(i0=0; i0<=N0; i0++){
      		int indx = Index(par, iDom,iField,i1,i0);
      		p[i0]=Sol.X[indx];
      	}
      	Chebyshev_Coefficients(p, cp, N0, par.grid_x0);

      	for(i0=0; i0<=N0; i0++){
      		c[i0][i1] = cp[i0];
      	}
      }
      free_dvector(p,0,N0); free_dvector(cp, 0, N0);
      
      sprintf(fn, "data/%s/cheb_field_%d_dom_%d_t_%03d.dat", par.SimName,iField,iDom, par.ncount);
      fp=fopen(fn,"w");
           
      //Header
      fprintf(fp,"#1: i0\t");
      for(i1=N1; i1>=0; i1--){
      	double sigma, dx;
      	get_sigma(par, iDom, i1, &sigma, &dx);
      	fprintf(fp,"#%d:c_i0(s=%3.5e)\t", N1-i1+2, sigma);
      }
      fprintf(fp, "\n");
      //end Header
      
      //Data
      for(i0=0; i0<=N0; i0++){
      	fprintf(fp,"%d\t", i0);
      	for(i1=N1; i1>=0; i1--){
      		fprintf(fp,"%3.15e\t", c[i0][i1] );
      	}
      	fprintf(fp,"\n");
      }
      fclose(fp);
      //end Data
      free_dmatrix(c, 0, N0,0, N1);      
    }
  }  

  return;
}
//------------------------------------------------------------------
void output_cheb_TimeEvolution(parameters par, derivs ID, derivs Sol){
  
  FILE *fp;
  char fn[1000];

  
  int iDom , iField,
      i0, i1, N0=par.N0, N1;
  
  for(iDom=0; iDom<nDoms; iDom++){
    N1=par.N1[iDom];
    for(iField=0; iField<nFields; iField++){
      sprintf(fn, "data/%s/Cheb_x1_TimeEvolution_field_%01d_dom_%01d.dat", par.SimName, iField, iDom);
      fp = par.ncount==0? fopen(fn, "w"): fopen(fn, "a"); 



      double *p=dvector(0, N1), *cp=dvector(0, N1), 
             *p_bound=dvector(0, N0), *cp_bound=dvector(0, N0),
             *p_part=dvector(0, N0), *cp_part=dvector(0, N0);

      if(par.ncount==0){
        double chi_0 = -1.;
        double tau = get_t_from_x0(par, chi_0);
        fprintf(fp, "\"Time = %3.20e\n", tau);
        
        for(i1=0; i1<=N1; i1++){    
    			int indx = Index_slice(par, iDom, iField, i1);
    			p[i1] = ID.X[indx];          
    		}
        Chebyshev_Coefficients(p, cp, N1, par.grid_x1[iDom]);

        for(i1=0; i1<=N1; i1++) fprintf(fp, "%d\t%3.20e\n",i1, fabs(cp[i1]));
        fprintf(fp, "\n");        
      }     
      
      
      for(i0=N0; i0>=0; i0--){
        double chi_0 = par.grid_points_x0[i0];
        double  tau = get_t_from_x0(par, chi_0);
        fprintf(fp, "\"Time = %3.20e\n", tau);

      	for(i1=0; i1<=N1; i1++){
      		int indx = Index(par, iDom,iField,i1,i0);
      		p[i1]=Sol.X[indx];
      	}
      	Chebyshev_Coefficients(p, cp, N1, par.grid_x1[iDom]);

        for(i1=0; i1<=N1; i1++) fprintf(fp, "%d\t%3.20e\n",i1, fabs(cp[i1]));
        fprintf(fp, "\n"); 

        int i1_bound = iDom==0? N1:0, i1_part = iDom==0? 0:N1;
        int indx_bound = Index(par, iDom,iField,i1_bound,i0),
            indx_part = Index(par, iDom,iField,i1_part,i0);
        
        p_bound[i0] = Sol.X[indx_bound];
        p_part[i0]  = Sol.X[indx_part];      	
      }
      free_dvector(p,0,N1); free_dvector(cp, 0, N1);
      fclose(fp);

      
      sprintf(fn, "data/%s/Cheb_x0_TimeEvolution_field_%01d_dom_%01d_bound.dat", par.SimName, iField, iDom);
      fp = par.ncount==0? fopen(fn, "w"): fopen(fn, "a"); 
      
      fprintf(fp, "\"Time = %3.20e\n", par.t1);
      Chebyshev_Coefficients(p_bound, cp_bound, N0, par.grid_x0);
      for(i0=0; i0<=N0; i0++) fprintf(fp, "%d\t%3.20e\n",i0, fabs(cp_bound[i0]));
      fprintf(fp, "\n");
      fclose(fp);

      sprintf(fn, "data/%s/Cheb_x0_TimeEvolution_field_%01d_dom_%01d_part.dat", par.SimName, iField, iDom);
      fp = par.ncount==0? fopen(fn, "w"): fopen(fn, "a"); 
      
      fprintf(fp, "\"Time = %3.20e\n", par.t1);
      Chebyshev_Coefficients(p_part, cp_part, N0, par.grid_x0);
      for(i0=0; i0<=N0; i0++) fprintf(fp, "%d\t%3.20e\n",i0, fabs(cp_part[i0]));
      fprintf(fp, "\n");
      fclose(fp);

          
    }
  }  

  return;
}
//----------------------------------------------------------------------------------
void output_Sol_particle(parameters par, derivs Sol){
 
  // int i0, N0 = par.N0;

  // char fn_dom0[500], fn_dom1[500],
  //      fn_dtau_dom0[500], fn_dtau_dom1[500];
  // sprintf(fn_dom0, "data/%s/psi_Part_dom0.dat", par.SimName);
  // sprintf(fn_dom1, "data/%s/psi_Part_dom1.dat", par.SimName);

  // sprintf(fn_dtau_dom0, "data/%s/Pi_Part_dom0.dat", par.SimName);
  // sprintf(fn_dtau_dom1, "data/%s/Pi_Part_dom1.dat", par.SimName);

  // FILE *fp_dom0, *fp_dom1, *fp_dom0_dtau, *fp_dom1_dtau;

  // if(par.ncount == 0){
  //   fp_dom0 = fopen(fn_dom0, "w");
  //   fprintf(fp_dom0,"#1:tau\t 2:Re(psi)\t 3:Im(phi)\t 4:Abs(psi)\n");
  //   fp_dom1 = fopen(fn_dom1, "w");
  //   fprintf(fp_dom1,"#1:tau\t 2:Re(psi)\t 3:Im(phi)\t 4:Abs(psi)\n");

  //   fp_dom0_dtau = fopen(fn_dtau_dom0, "w");
  //   fprintf(fp_dom0_dtau,"#1:tau\t 2:Re(Pi)\t 3:Im(Pi)\t 4:Abs(Pi)\n");
  //   fp_dom1_dtau = fopen(fn_dtau_dom1, "w");
  //   fprintf(fp_dom1_dtau,"#1:tau\t 2:Re(Pi)\t 3:Im(Pi)\t 4:Abs(Pi)\n");
  // }
  // else{
  //   fp_dom0 = fopen(fn_dom0, "a");
  //   fp_dom1 = fopen(fn_dom1, "a");
  //   fp_dom0_dtau = fopen(fn_dtau_dom0, "a");
  //   fp_dom1_dtau = fopen(fn_dtau_dom1, "a");
  // }

  // double *U_part_Re_dom0  = dvector(0,N0), *cU_part_Re_dom0  = dvector(0,N0),
  //        *U_part_Im_dom0  = dvector(0,N0), *cU_part_Im_dom0  = dvector(0,N0),  
  //        *U_part_Re_dom1  = dvector(0,N0), *cU_part_Re_dom1  = dvector(0,N0), 
  //        *U_part_Im_dom1  = dvector(0,N0), *cU_part_Im_dom1  = dvector(0,N0),
  //        *dU_part_Re_dom0 = dvector(0,N0), *cdU_part_Re_dom0 = dvector(0,N0),
  //        *dU_part_Im_dom0 = dvector(0,N0), *cdU_part_Im_dom0 = dvector(0,N0),  
  //        *dU_part_Re_dom1 = dvector(0,N0), *cdU_part_Re_dom1 = dvector(0,N0), 
  //        *dU_part_Im_dom1 = dvector(0,N0), *cdU_part_Im_dom1 = dvector(0,N0);

  // for(i0=0; i0<=N0; i0++){
  //   int Indx_part_Re_dom0 = Index(par, 0, 0, 0, i0),
  //       Indx_part_Im_dom0 = Index(par, 0, 1, 0, i0),
  //       Indx_part_Re_dom1 = Index(par, 1, 0, par.N1[1], i0),
  //       Indx_part_Im_dom1 = Index(par, 1, 1, par.N1[1], i0),
  //       Indx_part_Re_dom0_dtau = Index(par, 0, 2, 0, i0),
  //       Indx_part_Im_dom0_dtau = Index(par, 0, 3, 0, i0),
  //       Indx_part_Re_dom1_dtau = Index(par, 1, 2, par.N1[1], i0),
  //       Indx_part_Im_dom1_dtau = Index(par, 1, 3, par.N1[1], i0);

  //       U_part_Re_dom0[i0] = Sol.X[Indx_part_Re_dom0];
  //       U_part_Im_dom0[i0] = Sol.X[Indx_part_Im_dom0];
  //       U_part_Re_dom1[i0] = Sol.X[Indx_part_Re_dom1];
  //       U_part_Im_dom1[i0] = Sol.X[Indx_part_Im_dom1];

  //       dU_part_Re_dom0[i0] = Sol.X[Indx_part_Re_dom0_dtau];
  //       dU_part_Im_dom0[i0] = Sol.X[Indx_part_Im_dom0_dtau];
  //       dU_part_Re_dom1[i0] = Sol.X[Indx_part_Re_dom1_dtau];
  //       dU_part_Im_dom1[i0] = Sol.X[Indx_part_Im_dom1_dtau];
  // }
  // Chebyshev_Coefficients(U_part_Re_dom0, cU_part_Re_dom0, N0, par.grid_x0);
  // Chebyshev_Coefficients(U_part_Im_dom0, cU_part_Im_dom0, N0, par.grid_x0);
  // Chebyshev_Coefficients(U_part_Re_dom1, cU_part_Re_dom1, N0, par.grid_x0);
  // Chebyshev_Coefficients(U_part_Im_dom1, cU_part_Im_dom1, N0, par.grid_x0);
  // Chebyshev_Coefficients(dU_part_Re_dom0, cdU_part_Re_dom0, N0, par.grid_x0);
  // Chebyshev_Coefficients(dU_part_Im_dom0, cdU_part_Im_dom0, N0, par.grid_x0);
  // Chebyshev_Coefficients(dU_part_Re_dom1, cdU_part_Re_dom1, N0, par.grid_x0);
  // Chebyshev_Coefficients(dU_part_Im_dom1, cdU_part_Im_dom1, N0, par.grid_x0);

  // int j, j0 = par.ncount==0? 0:1;
  // for(j=j0; j<=jmax; j++){
  //   double x0 = -1.+2.*j/jmax, tau=get_t_from_x0(par, x0);
    
    
  //   double U_part_Re_dom0_output = Clenshaw_Chebyshev(cU_part_Re_dom0, N0, x0),
  //          U_part_Im_dom0_output = Clenshaw_Chebyshev(cU_part_Im_dom0, N0, x0),
  //          U_dom0_abs = sqrt( sqr(U_part_Re_dom0_output) + sqr(U_part_Im_dom0_output)),

  //          U_part_Re_dom1_output = Clenshaw_Chebyshev(cU_part_Re_dom1, N0, x0),
  //          U_part_Im_dom1_output = Clenshaw_Chebyshev(cU_part_Im_dom1, N0, x0),
  //          U_dom1_abs = sqrt( sqr(U_part_Re_dom1_output) + sqr(U_part_Im_dom1_output)),

  //          dU_part_Re_dom0_output = Clenshaw_Chebyshev(cdU_part_Re_dom0, N0, x0),
  //          dU_part_Im_dom0_output = Clenshaw_Chebyshev(cdU_part_Im_dom0, N0, x0),
  //          dU_dom0_abs = sqrt( sqr(dU_part_Re_dom0_output) + sqr(dU_part_Im_dom0_output)),

  //          dU_part_Re_dom1_output = Clenshaw_Chebyshev(cdU_part_Re_dom1, N0, x0),
  //          dU_part_Im_dom1_output = Clenshaw_Chebyshev(cdU_part_Im_dom1, N0, x0),
  //          dU_dom1_abs = sqrt( sqr(dU_part_Re_dom1_output) + sqr(dU_part_Im_dom1_output));

  //          fprintf(fp_dom0, "%3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, U_part_Re_dom0_output, U_part_Im_dom0_output, U_dom0_abs);
  //          fprintf(fp_dom1, "%3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, U_part_Re_dom1_output, U_part_Im_dom1_output, U_dom1_abs);

  //          fprintf(fp_dom0_dtau, "%3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, dU_part_Re_dom0_output, dU_part_Im_dom0_output, dU_dom0_abs);
  //          fprintf(fp_dom1_dtau, "%3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, dU_part_Re_dom1_output, dU_part_Im_dom1_output, dU_dom1_abs);
  // }
  // fclose(fp_dom0);
  // fclose(fp_dom1);
  // fclose(fp_dom0_dtau);
  // fclose(fp_dom1_dtau);

  // free_dvector(U_part_Re_dom0, 0, N0); free_dvector(cU_part_Re_dom0, 0, N0);
  // free_dvector(U_part_Im_dom0, 0, N0); free_dvector(cU_part_Im_dom0, 0, N0);
  // free_dvector(U_part_Re_dom1, 0, N0); free_dvector(cU_part_Re_dom1, 0, N0);
  // free_dvector(U_part_Im_dom1, 0, N0); free_dvector(cU_part_Im_dom1, 0, N0);
  return;
}
//----------------------------------------------------------------------------------
void output_Sol(parameters par, derivs Sol){
 
  int N1, N0=par.N0,  i0, i1, iField, idom;
  double **U, **cU, *UBound=dvector( 0,N0), *cUBound=dvector(0,N0);
  
  FILE *fpU, *fp_Ugnuplot, *fp2_Ugnuplot;
  char fnU[500], fn_Ugnuplot[500], fn2_Ugnuplot[500];
  
  for(iField=0; iField < nFields1D; iField++){
    sprintf(fnU, "data/%s/gnuplot_FreeBound_%d.dat", par.SimName,iField);
    for(i0=0; i0<=N0; i0++){
      int indx = Index_field1D(par, iField, i0);
      UBound[i0]    = Sol.X[indx];
    }
    Chebyshev_Coefficients(UBound, cUBound, N0, par.grid_x0);
    
    fpU = par.ncount==0? fopen(fnU, "w"): fopen(fnU, "a");
    int j,  j0 = par.ncount==0? 0:1;
    for(j=j0; j<=jmax; j++){
    	double x0 = -1.+2.*j/jmax, t=get_t_from_x0(par, x0), U_output;
      
      U_output = Clenshaw_Chebyshev(cUBound, N0, x0);
      fprintf(fpU, "%3.20e\t%3.20e\n", t, U_output);
    }
    fclose(fpU);
  }
  
  for(idom=0; idom<nDoms; idom++){
    N1=par.N1[idom];
    U=dmatrix(0,N1, 0,N0); 
    cU=dmatrix(0,N1, 0,N0);
    
    for(iField=0; iField<nFields; iField++){
    	for(i1=0; i1<=N1; i1++){
    		for(i0=0; i0<=N0; i0++){
    			int indx = Index(par, idom, iField, i1, i0);
    			U[i1][i0]=Sol.X[indx];
    		}
    	}
    	Chebyshev_Coefficients_2D(U, cU, N1, par.grid_x1[idom], N0, par.grid_x0);

    	sprintf(fn_Ugnuplot, "data/%s/gnuplot3D_field_%01d_dom_%01d.dat", par.SimName, iField, idom);
    	sprintf(fn2_Ugnuplot, "data/%s/gnuplot2D_field_%01d_dom_%01d.dat", par.SimName, iField, idom);

    	fp_Ugnuplot = par.ncount==0? fopen(fn_Ugnuplot, "w"): fopen(fn_Ugnuplot, "a");
    	fp2_Ugnuplot = par.ncount==0? fopen(fn2_Ugnuplot, "w"): fopen(fn2_Ugnuplot, "a");

    	int i, j, j0= par.ncount==0? 0:1;

    	//----gnuplot2D header
    	if(par.ncount==0){
    		fprintf(fp2_Ugnuplot, "#1:T\t");
    		for(i=0; i<=imax; i++){
    			double chi=-1.+2.*i/imax,
    						x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi), 
    						sigma = get_sigma_from_x1(par, idom, x1);
    			
    			fprintf(fp2_Ugnuplot, "#%d:%3.3e\t", i+2, sigma);
    		}
    		fprintf(fp2_Ugnuplot, "\n");
    	}
    	//----gnuplot2D header

    	//----gnuplot3D header
    	fprintf(fp_Ugnuplot, "#1:t\t2:s\t3:Field\n");
    	//----gnuplot3D header
    	for(j=j0; j<=jmax; j++){
    		double x0 = -1.+2.*j/jmax, t=get_t_from_x0(par, x0), U_output;

		  	fprintf(fp2_Ugnuplot, "%3.20e\t", t);

		  	for(i=0; i<=imax; i++){
		  		double chi=-1.+2.*i/imax,
    						x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi), 
    						sigma = get_sigma_from_x1(par, idom, x1);

		  		U_output = Clenshaw_Chebyshev_2D(cU, N1, N0, chi, x0);

		  		fprintf(fp_Ugnuplot, "%3.20e\t%3.20e\t%3.20e\n", t,sigma, U_output);
		  		fprintf(fp2_Ugnuplot, "%3.20e\t", U_output );
		  	}
		  	fprintf(fp_Ugnuplot, "\n");
		  	fprintf(fp2_Ugnuplot, "\n");
		  }
		  fclose(fp_Ugnuplot);
		  fclose(fp2_Ugnuplot);
		}
		free_dmatrix(U,0,N1, 0,N0 ); free_dmatrix(cU,0,N1, 0,N0 );
	}
	free_dvector(UBound, 0, N0);  free_dvector(cUBound, 0, N0);

	return;
}
//----------------------------------------------------------------------------------
void output_TimeEvolution(parameters par, derivs ID, derivs Sol){
 
  int N1, N0=par.N0,  i0, i1, iField, idom;
    
  FILE *fpU;
  char fnU[1000];
  
  for(idom=0; idom<nDoms; idom++){
    N1=par.N1[idom];
   
    for(iField=0; iField<nFields; iField++){
      sprintf(fnU, "data/%s/TimeEvolution_field_%01d_dom_%01d.dat", par.SimName, iField, idom);
      fpU = par.ncount==0? fopen(fnU, "w"): fopen(fnU, "a"); 
      
      
      if(par.ncount==0){
        double chi_0 = -1.;
        double tau = get_t_from_x0(par, chi_0);
        fprintf(fpU, "\"Time = %3.20e\n", tau);
        for(i1=0; i1<=N1; i1++){   
          double chi_1= par.grid_points_x1[idom][i1], U_output,
                 x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi_1), 
    						 sigma = get_sigma_from_x1(par, idom, x1);
    			int indx = Index_slice(par, idom, iField, i1);
    			U_output=ID.X[indx];
          fprintf(fpU, "%3.20e\t%3.20e\n",sigma, U_output);
    		}
        fprintf(fpU, "\n");        
      }
      
      
      
      for(i0=N0; i0>=0; i0--){
        double chi_0 = par.grid_points_x0[i0];
        double tau = get_t_from_x0(par, chi_0);
        fprintf(fpU, "\"Time = %3.20e\n", tau);

    	  for(i1=0; i1<=N1; i1++){   
          double chi_1= par.grid_points_x1[idom][i1], U_output,
                 x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi_1), 
    						 sigma = get_sigma_from_x1(par, idom, x1);
    			int indx = Index(par, idom, iField, i1, i0);
    			U_output=Sol.X[indx];
          fprintf(fpU, "%3.20e\t%3.20e\n",sigma, U_output);
    		}
        fprintf(fpU, "\n");
    	}
      fclose(fpU);    
		}
	}
	
	return;
}
//----------------------------------------------------------------------------------
void output_Sol_derv_sigma(parameters par, derivs Sol){
 
  // int N1, N0=par.N0,  i0, i1, iField, idom;
  // double **U, **cU, *UBound=dvector( 0,N0), *cUBound=dvector(0,N0);
  
  // FILE *fp_Ugnuplot, *fp2_Ugnuplot;
  // char fn_Ugnuplot[500], fn2_Ugnuplot[500];

  // derivs Sol_x0x1, Sol_sigtau;
  // allocate_derivs(&Sol_x0x1, 0, 1);
  // allocate_derivs(&Sol_sigtau, 0, 1);

  // for(idom=0; idom<nDoms; idom++){
  //   N1=par.N1[idom];
  //   U=dmatrix(0,N1, 0,N0); 
  //   cU=dmatrix(0,N1, 0,N0);
    
  //   for(iField=0; iField<nFields; iField++){
  //     for(i1=0; i1<=N1; i1++){
  //       for(i0=0; i0<=N0; i0++){
  //         int indx = Index(par, idom, iField, i1, i0);
  //         Sol_x0x1.X[0] = Sol.X[indx];     Sol_x0x1.Y[0]   = Sol.Y[indx];
  //         Sol_x0x1.X1[0] = Sol.X1[indx];   Sol_x0x1.Y1[0]  = Sol.Y1[indx];
  //         Sol_x0x1.X11[0] = Sol.X11[indx]; Sol_x0x1.Y11[0] = Sol.Y11[indx]; 
  //         scale_derivs_x0x1_to_tau_sigma(par, idom, i1, Sol_x0x1, Sol_sigtau, 1);
  //         U[i1][i0]=Sol_sigtau.X1[0];
  //       }
  //     }
  //     Chebyshev_Coefficients_2D(U, cU, N1, par.grid_x1[idom], N0, par.grid_x0);

  //     sprintf(fn_Ugnuplot, "data/%s/gnuplot3D_field_sigma_%01d_dom_%01d.dat", par.SimName, iField, idom);
  //     sprintf(fn2_Ugnuplot, "data/%s/gnuplot2D_field_sigma_%01d_dom_%01d.dat", par.SimName, iField, idom);

  //     fp_Ugnuplot = par.ncount==0? fopen(fn_Ugnuplot, "w"): fopen(fn_Ugnuplot, "a");
  //     fp2_Ugnuplot = par.ncount==0? fopen(fn2_Ugnuplot, "w"): fopen(fn2_Ugnuplot, "a");

  //     int i, j,  j0= par.ncount==0? 0:1;

  //     //----gnuplot2D header
  //     if(par.ncount==0){
  //       fprintf(fp2_Ugnuplot, "#1:T\t");
  //       for(i=0; i<=imax; i++){
  //         double chi=-1.+2.*i/imax,
  //               x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi), 
  //               sigma = get_sigma_from_x1(par, idom, x1);
          
  //         fprintf(fp2_Ugnuplot, "#%d:%3.3e\t", i+2, sigma);
  //       }
  //       fprintf(fp2_Ugnuplot, "\n");
  //     }
  //     //----gnuplot2D header

  //     //----gnuplot3D header
  //     fprintf(fp_Ugnuplot, "#1:t\t2:s\t3:Field\n");
  //     //----gnuplot3D header
  //     for(j=j0; j<=jmax; j++){
  //       double x0 = -1.+2.*j/jmax, t=get_t_from_x0(par, x0), U_output;

  //       fprintf(fp2_Ugnuplot, "%3.20e\t", t);

  //       for(i=0; i<=imax; i++){
  //         double chi=-1.+2.*i/imax,
  //               x1 = get_x_from_chi(par.AnMR_x_boundary[idom], par.AnMR_kappa[idom], chi), 
  //               sigma = get_sigma_from_x1(par, idom, x1);

  //         U_output = Clenshaw_Chebyshev_2D(cU, N1, N0, x1, x0);

  //         fprintf(fp_Ugnuplot, "%3.20e\t%3.20e\t%3.20e\n", t,sigma, U_output);
  //         fprintf(fp2_Ugnuplot, "%3.20e\t", U_output );
  //       }
  //       fprintf(fp_Ugnuplot, "\n");
  //       fprintf(fp2_Ugnuplot, "\n");
  //     }
  //     fclose(fp_Ugnuplot);
  //     fclose(fp2_Ugnuplot);
  //   }
  //   free_dmatrix(U,0,N1, 0,N0 ); free_dmatrix(cU,0,N1, 0,N0 );
  // }
  // free_dvector(UBound, 0, N0);  free_dvector(cUBound, 0, N0);
  // free_derivs(&Sol_x0x1, 0, 1);
  // free_derivs(&Sol_sigtau, 0, 1);

  return;
}
//---------------------------------------------------------------------
void output_Source(parameters par){
  int idom, i1, i0, j0= par.ncount==0? 0:1;

  FILE *fp_Re, *fp_Im;
  char fn_Re[500], fn_Im[500];
  

  for(idom=0; idom<nDoms; idom++){
    int N1=par.N1[idom];

    sprintf(fn_Re, "data/%s/Re_Source_dom_%d.dat", par.SimName, idom);
    sprintf(fn_Im, "data/%s/Im_Source_dom_%d.dat", par.SimName, idom);
    if(par.ncount==0){
      fp_Re   =fopen(fn_Re, "w");  
      fp_Im   =fopen(fn_Im, "w");  

      fprintf(fp_Re,   "#1:tau\t ");
      fprintf(fp_Im,   "#1:tau\t ");
      for(i1=0; i1<=N1; i1++) {
        fprintf(fp_Re,   "#%d:S^2(%3.2f)\t ", i1+2, par.grid_points_x1[idom][i1]);
        fprintf(fp_Im,   "#%d:S^2(%3.2f)\t ", i1+2, par.grid_points_x1[idom][i1]);
      }
      fprintf(fp_Re,   "\n");
      fprintf(fp_Im,   "\n");
    }
    else{
      fp_Re   = fopen(fn_Re, "a");
      fp_Im   = fopen(fn_Im, "a");
    }

    for(i0=j0; i0<=jmax; i0++){
        double x0 = -1.+2.*i0/jmax, t=get_t_from_x0(par, x0);
        fprintf(fp_Re,   "%3.15e\t ", t);
        fprintf(fp_Im,   "%3.15e\t ", t);

      for(i1=0; i1<=N1; i1++){
        double complex Source = Get_Source(par, idom,  x0, i1);
        fprintf(fp_Re,   "%3.15e\t ", creal(Source));
        fprintf(fp_Im,   "%3.15e\t ", cimag(Source));
      }
      fprintf(fp_Re,   "\n");
      fprintf(fp_Im,   "\n");
    }
    fclose(fp_Re);
    fclose(fp_Im);
  }

  return;
}
//--------------------------------------------------------------------
void output_Fr(parameters par, derivs Sol){
  // int i0, N0 = par.N0;
  // double complex Fr_SSF_lm_dom0, Fr_SSF_lm_dom1;

  // FILE *fp_Fr;
  // char fn_Fr[500];
  // sprintf(fn_Fr, "data/%s/Fr.dat", par.SimName);

  // if(par.ncount==0){
  //   fp_Fr   =fopen(fn_Fr, "w");  
  //   fprintf(fp_Fr,   "#1:tau\t 2:Re(Fr_SSF_dom0)\t 3:Re(Fr_SSF_dom1)\n");
  // }
  // else{
  //  fp_Fr   = fopen(fn_Fr, "a");
  // }

  // for(i0=N0; i0>=0; i0--){
  //   double tau = get_t_from_x0(par, par.grid_points_x0[i0]);
  //   get_Fr_SSF_lm(par, Sol, i0, &Fr_SSF_lm_dom0, &Fr_SSF_lm_dom1);    
  //   fprintf(fp_Fr, "%3.15e\t %3.15e\t %3.15e\n", tau, creal(Fr_SSF_lm_dom0), creal(Fr_SSF_lm_dom1) );
    

  // }
  // fclose(fp_Fr);

  // return;
}
//--------------------------------------------------------------------
void output_Fluxes_and_Ft(parameters par, derivs Sol){
  // int i0, N0 = par.N0;
  // double complex Ft_SSF_lm_dom0, Ft_SSF_lm_dom1;
  // double Flux_lm_hrzn, Flux_lm_scri, Flux_lm, Ft_lm_flux, rel_error_dom0, rel_error_dom1;

  // FILE *fp_flux, *fp_Ft;
  // char fn_flux[500], fn_Ft[500];

  // sprintf(fn_flux, "data/%s/EnergyFluxes.dat", par.SimName);
  // sprintf(fn_Ft, "data/%s/Ft.dat", par.SimName);

  // if(par.ncount==0){
  //   fp_flux = fopen(fn_flux, "w"); 
  //   fp_Ft   =fopen(fn_Ft, "w");
  //   fprintf(fp_flux, "#1:tau\t 2:Flux Scri\t 3:Flux Horizon\t 4:Total Flux\n");
  //   fprintf(fp_Ft,   "#1:tau\t 2:Re(Ft_SSF_dom0)\t 3:Error_balance_dom0\t 4:Re(Ft_SSF_dom1)\t 5:Error_balance_dom0\n");
  // }
  // else{
  //  fp_flux = fopen(fn_flux, "a"); 
  //  fp_Ft   = fopen(fn_Ft, "a");
  // }

  // for(i0=N0; i0>=0; i0--){
  //   double tau = get_t_from_x0(par, par.grid_points_x0[i0]);
  //   Ft_lm_flux = get_Ft_FLux_lm(par, Sol, i0);
    
  //   get_Flux_lm(par, Sol, i0, &Flux_lm_hrzn, &Flux_lm_scri, &Flux_lm);
  //   get_Ft_SSF_lm(par, Sol, i0, &Ft_SSF_lm_dom0, &Ft_SSF_lm_dom1);

  //   rel_error_dom0 = fabs(1. - creal(Ft_SSF_lm_dom0)/Ft_lm_flux );
  //   rel_error_dom1 = fabs(1. - creal(Ft_SSF_lm_dom1)/Ft_lm_flux );

  //   fprintf(fp_flux, "%3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, Flux_lm_scri, Flux_lm_hrzn, Flux_lm);
  //   fprintf(fp_Ft, "%3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\n", tau, creal(Ft_SSF_lm_dom0), rel_error_dom0, creal(Ft_SSF_lm_dom1), rel_error_dom1);
    

  // }

  // fclose(fp_flux);
  // fclose(fp_Ft);

  return;
}

