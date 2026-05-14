#include "1+1_Free_Boundaries.h"

//-----------------------------------------------
double get_grid_angle_incremente(char *grid, int N){
  
  double arg=0;
  
  if(strcmp( grid,"Radau_RHS" ) ==0)
    arg = 2*Pi/(2*N+1);   
  else if(strcmp(grid,"Radau_LHS")==0)
    arg = - 2*Pi/(2*N+1);
  else if(strcmp(grid,"Gauss")==0)
    arg = Pi/(N+1);
  else if(strcmp( grid,"Lobatto")==0)
    arg = Pi/N;
  else if(strcmp( grid,"Fourier")==0){
    arg = 2.*Pi/(2*N+1);
  }
  else{ 
    fprintf(stderr,
    "Error in get_grid_angle_incremente: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto / Fourier\n input grid was: %s\n", grid);
    exit(1);
  }
  
  return arg;
}
//------------------------------------------

double get_grid_point(char *grid, int i, int N){
  
  double arg=0, x=0;
  
  if(strcmp( grid,"Radau_RHS" ) ==0){
    arg = 2*Pi*i/(2*N+1);
    x = cos(arg); 
  }
  else if(strcmp(grid,"Radau_LHS")==0){
    arg = Pi - 2*Pi*i/(2*N+1);
    x = cos(arg);
  }
  else if(strcmp(grid,"Gauss")==0){
    arg = Pi*(i+0.5)/(N+1);
    x = cos(arg);
  }
  else if(strcmp( grid,"Lobatto")==0){
    arg = Pi*i/N;
    x = cos(arg);
  }
  else if(strcmp( grid,"Fourier")==0){
    x = 2.*Pi*i/(2*N+1);
  }
  else{ 
    fprintf(stderr,
    "Error in get_grid_point: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto / Fourier\n input grid was: %s\n", grid);
    exit(1);
  }  
  
  return x;
}
//----------------------------------------------------------------
double get_x_from_chi(double x_boundary, double kappa, double chi){
  //Input: x_boundary = -1 (Analytical Mesh Refinement at domain's Left Boundary)
  //       x_boundary =  1 (Analytical Mesh Refinement at domain's Right Boundary)

  double x;

  if(kappa==0.){
    x = chi;
  }
  else{
    if(1.-sqr(x_boundary)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundary was: %lf\n", x_boundary);
      exit(1);
    }
    x = x_boundary*( 1. - 2.* sinh(kappa* (1-chi*x_boundary))/sinh(2*kappa) );
  }
  return x;
}
//---------------------------------------------------------------------
double get_t(parameters par, int i0){
  double x0 = par.grid_points_x0[i0], t = par.t0 + 0.5*(par.t1-par.t0)*(1+x0);

  return t;
}
//----------------------------------------------------------------
void get_sigma(parameters par, int iDom, int i1, double *sigma, double *dx_dsigma){
  
  int N=par.N1[iDom];
  double chi, x, x_B = par.AnMR_x_boundary[iDom], kappa = par.AnMR_kappa[iDom];


  chi = get_grid_point(par.grid_x1[iDom], i1, N);
  x = get_x_from_chi(x_B, kappa, chi);

  *sigma =0.5*(par.sigma[iDom+1]*(1.+ x) + par.sigma[iDom]*(1.- x));
  *dx_dsigma=2./(par.sigma[iDom+1]-par.sigma[iDom]);
  
  return;
}

//---------------------------------------------------------------------
double get_sigma_from_x1(parameters par, int iDom, double x){
  
  double sigma;

  sigma =0.5*(par.sigma[iDom+1]*(1.+ x) + par.sigma[iDom]*(1.- x));
  return sigma;
}
//---------------------------------------------------------------------
double get_x1_from_sigma(parameters par, int iDom, double sigma){
  
  double x, 
        sigma_sum = par.sigma[iDom] + par.sigma[iDom-1] , 
        sigma_dif = par.sigma[iDom] - par.sigma[iDom-1] ; 

  x =  (2*sigma - sigma_sum)/sigma_dif;
  return x;
}
//---------------------------------------------------------------------
double get_t_from_x0(parameters par, double x0){
 
 return par.t0 + 0.5*(par.t1-par.t0)*(1+x0);
}
//-----------------------------------------------
void scale_derivs_x0x1_to_tau_sigma(parameters par, int iDom, int i1, derivs W_x0x1, derivs W_tau_sigma, int nsize){
  int N=par.N1[iDom];
  double x_B = par.AnMR_x_boundary[iDom], kappa = par.AnMR_kappa[iDom], chi, sigma, dx_dsigma, dx_dsigma_2,
         fact = 2./(par.t1-par.t0);

  chi = get_grid_point(par.grid_x1[iDom], i1, N);
  get_sigma(par, iDom, i1, &sigma, &dx_dsigma);
  dx_dsigma_2 = sqr(dx_dsigma);
  
 int i;

 for(i=0; i<nsize; i++ ){
    // //Map derivatives w.r.t to chi [spectral] into derivatives w.r.t to x [AnMr]
    if(kappa==0.){
      W_tau_sigma.X1[i]  = W_x0x1.X1[i];
      W_tau_sigma.X11[i] = W_x0x1.X11[i];

      W_tau_sigma.Y1[i]  = W_x0x1.Y1[i];
      W_tau_sigma.Y11[i] = W_x0x1.Y11[i];
    }
    else{
      if(1.-sqr(x_B)!= 0. ){
        fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundary was: %lf\n", x_B);
        exit(1);
      }
      double J, J2;
      J= 2.*kappa* cosh(kappa* (1-chi*x_B))/sinh(2*kappa);
      J2 = -2.*sqr(kappa)*x_B* sinh(kappa* (1-chi*x_B))/sinh(2*kappa);

      W_tau_sigma.X1[i] = W_x0x1.X1[i]/J;
      W_tau_sigma.X11[i] = (W_x0x1.X11[i] - J2* W_tau_sigma.X1[i] )/sqr(J);

      W_tau_sigma.Y1[i] = W_x0x1.Y1[i]/J;
      W_tau_sigma.Y11[i] = (W_x0x1.Y11[i] - J2* W_tau_sigma.Y1[i] )/sqr(J);
    }
    //Map derivatives w.r.t to x [AnMr] into derivatives w.r.t to sigma [Hyperboloidal]
    
    W_tau_sigma.X[i] = W_x0x1.X[i];
    W_tau_sigma.X1[i] *= dx_dsigma;
    W_tau_sigma.X11[i] *= dx_dsigma_2;

    W_tau_sigma.Y[i] = fact*W_x0x1.Y[i];
    W_tau_sigma.Y1[i] *= fact*dx_dsigma;
    W_tau_sigma.Y11[i] *= fact*dx_dsigma_2;
  }

  return; 
}
//-----------------------------------------------
void scale_derivs_x0_to_tau(parameters par, derivs W1D_x0, derivs W1D_tau, int nsize_1D){
  int i;
  double t0=par.t0, t1=par.t1, fact = 2./(t1-t0);

  for(i=0; i<nsize_1D; i++ ){
    W1D_tau.X[i] = W1D_x0.X[i];
    W1D_tau.Y[i] = fact*W1D_x0.Y[i];
  }
  return; 
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
double get_z_from_chi(parameters par, double chi){
  int flag_NewCoord = par.flag_NewCoord;
  double z;
  
  if(flag_NewCoord==0)
    z= chi;
  else
    z = chi==0.? 0.: exp(1. - 1./chi);
    
  return z;
  
}
//-----------------------------------------------
double get_chi_from_z(parameters par, double z){
  int flag_NewCoord = par.flag_NewCoord;
  double chi;
  
  if(flag_NewCoord==0)
    chi = z;
  else
    chi = z==0.? 0.: 1./(1.-log(z));
    
    
    
  return chi;
  
}