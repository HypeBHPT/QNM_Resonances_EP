#include "1+1_Free_Boundaries.h"

//--------------------------------------------------------------------------------------
double get_r_over_M_from_rStar(parameters par, double rStar){
  double M = par.Mass, rh=par.rh_over_M, W_arg, W, r, r_over_M;


    W_arg = exp(-1. + rStar/rh);
    W = gsl_sf_lambert_W0(W_arg);

    r = (1. + W)*rh;
    r_over_M = r/M;

    //CONSISTENT CHECK
   double rStar_over_lambda, sigma0, dsigma0_dr, x, dx_dsigma, lambda = par.lambda_over_M;
   rStar_over_lambda = rStar/lambda;
	 func_sigma_of_r(par, r_over_M, &sigma0, &dsigma0_dr);
	 func_tortoise_x(par, sigma0, &x, &dx_dsigma);
   if(fabs(1.- x/rStar_over_lambda)>1.e-15){
    fprintf(stderr, "Error in function get_r_over_M_from_rStar: rStar=%3.15e not being recoverded (x = %3.15e, error = %3.15e)\n",rStar_over_lambda,x, fabs(1.-x/rStar_over_lambda) );
    exit(-1);
   }
	
  return r_over_M;
}
//--------------------------------------------------------------------------------------
void func_rho(parameters par, double sigma, double *rho, double *drho_dsigma, double *d2rho_dsigma2){
  double rho_0 = par.rho_0, rho_1 = par.rho_1;


  *rho = rho_0 + rho_1*sigma;
  *drho_dsigma = rho_1;
  *d2rho_dsigma2 = 0.;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_beta(parameters par, double sigma, double *beta, double *dbeta_dsigma){
  double rho, drho, d2rho;
  func_rho(par, sigma, &rho, &drho, &d2rho);
  
  *beta = rho - sigma*drho;
  *dbeta_dsigma = - sigma*d2rho;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_r_of_sigma(parameters par, double sigma, double *r_over_M, double *dr_dsigma){
  double rho, drho, d2rho, beta, dbeta;
  
  func_rho(par, sigma, &rho, &drho, &d2rho);
  func_beta(par, sigma, &beta, &dbeta);
  
  
  *r_over_M = par.lambda_over_M*rho/sigma;
  *dr_dsigma = -par.lambda_over_M*beta/sqr(sigma);
  
  return;
}
//--------------------------------------------------------------------------------------
void func_sigma_of_r(parameters par, double r_over_M, double *sigma, double *dsigma_dr){
  
  *sigma = par.lambda_over_M*par.rho_0/r_over_M;
  *dsigma_dr = - par.lambda_over_M*par.rho_0/sqr(r_over_M);
  
  return;
}
//--------------------------------------------------------------------------------------
void func_f(double r_over_M, double *f, double *df_dr){     
  
  *f=(1.-2./r_over_M);
  *df_dr=2/sqr(r_over_M);

  return;
}
//--------------------------------------------------------------------------------------
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma){
  

  *x= 2*( 1./sigma + log(1.-sigma) - log(sigma) )/par.lambda_over_M;

  double r,dr_dsigma;  
  func_r_of_sigma(par, sigma, &r, &dr_dsigma);    
    
  double f, df_dr;
  func_f(r, &f, &df_dr);

  double beta, dbeta;
  func_beta(par, sigma, &beta, &dbeta);
  
  double drstar_dr = 1./f, dx_drstar = 1./par.lambda_over_M;
  
  *dx_dsigma = dx_drstar*drstar_dr*dr_dsigma;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_height(parameters par, double sigma, double *h, double *dh_dsigma){
 //HEIGHT FUNCTION h(sigma) WITH RESPECT TO INGOING EDDINGTON FINKELSTEIN TIME v

  *h = -2*par.rho_0*(1./sigma - 2*log(sigma)/(par.lambda_over_M*par.rho_0)  );
  *dh_dsigma = 2*par.rho_0*(1./sqr(sigma) + 2./(sigma*par.lambda_over_M*par.rho_0)  );
  
  return ;
}
//--------------------------------------------------------------------------------------
void func_Height(parameters par,double sigma, double *H, double *dH_dsigma){
//HEIGHT FUNCTION H(sigma) WITH RESPECT TO SCHWARZSCHILD TIME t
  double rh = par.rh_over_M, lambda = par.lambda_over_M, rh_over_lambda = rh/lambda;

  double p, dp;
  func_p(par, sigma, &p, &dp);
  double gamma, dgamma;
  func_gamma(par, sigma, &gamma, &dgamma);

  *H = rh_over_lambda * ( -1./sigma + log(sigma) + log(1-sigma)  );
  *dH_dsigma = gamma/p;
  
  
  return;
}
//--------------------------------------------------------------------------------------
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma){

  *Omega = sigma/par.lambda_over_M;
  *dOmega_dsigma = 1./par.lambda_over_M;
  return;
}
//--------------------------------------------------------------------------------------
void func_Z(parameters par, double sigma, double complex *Z, double complex *dlnZ_dsigma){
  double complex s=par.s;
   double H, dH;
  
  
  func_height(par, sigma, &H, &dH);
  
  *Z = cexp( s*(H) );
  
  *dlnZ_dsigma = s*dH; 
  
  return;
}
//------------------------------------------------------------
double func_alpha2(parameters par, double sigma){
  double a2, sigma2=sqr(sigma), rh = par.rh_over_M, lambda = par.lambda_over_M;

  if(strcmp( par.Equation,"BardeenPress") ==0){
    a2 = sigma2*(1.-sigma); 
  }
  else{
    a2 = lambda*sigma2*(1.-sigma)/rh;//*rh/lambda
  }
  
  return a2;
}//------------------------------------------------------------
double func_alpha1(parameters par, double sigma){
  double a1, spin = 1.*par.spin, rh = par.rh_over_M, lambda = par.lambda_over_M;

  if(strcmp( par.Equation,"BardeenPress") ==0){
    a1 = sigma*( 2. -3.*sigma +   spin*(2 - sigma)  ); 
  }
  else{
    a1 = lambda * sigma * (2 - 3 *sigma)/rh; //*rh/lambda
  }
  
  return a1;
}
//------------------------------------------------------------
double func_alpha0(parameters par, double sigma){
  double a0, rh = par.rh_over_M, lambda = par.lambda_over_M, spin = 1.*par.spin, l = 1.*par.ell, V_bump;

  func_V_bump(par, sigma, &V_bump);
  if(strcmp( par.Equation, "BardeenPress") ==0){
    a0 = - (   l*(l+1.) -  spin*(1.+spin) + sigma*(1+spin)  ) - V_bump; 
  }
  else if(strcmp( par.Equation,"ReggeWheeler") ==0){
    a0 = - lambda/rh * ( l*(l+1) +  (1- sqr(spin))* sigma    ) - V_bump ; //*rh/lambda

  }
  else{
    printf("Error if function func_alpha0: Potential for %s not implemented\n", par.Equation);
    exit(-1);
  }
  
  
  return a0;
}
//------------------------------------------------------------
double func_beta1(parameters par, double sigma){
  double b1, rh_over_lambda = par.rh_over_M/par.lambda_over_M, sigma2=sqr(sigma);

  if(strcmp( par.Equation,"BardeenPress") ==0){
    b1 = 2. * rh_over_lambda*( 1.-2*sigma2); 
  }
  else{
    b1 = 2*(1 - 2 * sigma2); //*rh/lambda
  }
  
  
  return b1;
}
//------------------------------------------------------------
double func_beta0(parameters par, double sigma){
  double b0, spin = 1.*par.spin, rh_over_lambda = par.rh_over_M/par.lambda_over_M;

  if(strcmp( par.Equation,"BardeenPress") ==0){
    b0 = - 2. * rh_over_lambda*( 2.*sigma -   spin *(1. - sigma) ); 
  }
  else{
    b0 = -4*sigma;//*rh/lambda
  }
  return b0;
}
//------------------------------------------------------------
double func_w(parameters par, double sigma){
  double w, rh_over_lambda = par.rh_over_M/par.lambda_over_M;

  if(strcmp( par.Equation,"BardeenPress") ==0){
      w = sqr( 2. * rh_over_lambda) * (1+sigma); 
  }
  else{
      w = 4. * rh_over_lambda * (1+sigma); //*rh/lamdda
  }

  return w;
}
//------------------------------------------------------------
void func_p(parameters par, double sigma, double *p, double *dp_dsigma){
  double rh = par.rh_over_M ,lambda=par.lambda_over_M, lambda_over_rh =lambda/rh, sigma2=sqr(sigma);

  *p = lambda_over_rh*sigma2*(1.-sigma);
  *dp_dsigma = lambda_over_rh*sigma*(2.-3.*sigma);

  return;
}
//------------------------------------------------------------
void func_gamma(parameters par, double sigma, double *gamma, double *dgamma_dsigma){
  double sigma2=sqr(sigma);


  *gamma = 1. - 2.*sigma2;
  *dgamma_dsigma = - 4.*sigma;

  return;
}
//------------------------------------------------------------
void func_w_RWZ(parameters par, double sigma, double *w){
  double rh = par.rh_over_M ,lambda=par.lambda_over_M, rh_over_lambda = rh/lambda;


  *w = 4. * rh_over_lambda * (1+sigma);

  return;
}
//------------------------------------------------------------
void func_V(parameters par, double sigma, double *V){
  double l = 1.*par.ell, rh = par.rh_over_M, lambda=par.lambda_over_M, lambda_over_rh =lambda/rh, s=par.spin, s2=sqr(s);


  *V = lambda_over_rh*(l*(l+1) + (1-s2)*sigma );

  return;
}
//------------------------------------------------------------
void func_V_bump(parameters par, double sigma, double *V){
  double eps = par.epsilon, a0 = par.a0_over_lambda, x, dx_dsigma, p, dp_dsigma;

  func_tortoise_x(par, sigma, &x, &dx_dsigma);
  func_p(par, sigma, &p, &dp_dsigma);

  if(sigma ==0. || sigma ==1.){
    *V = 0.;
  }
  else{
    *V = eps/(sqr(cosh( x - a0)) * p);
  }

  

  return;
}
//------------------------------------------------------------
double complex Diff_Operator(parameters par, complex_derivs W, int iDom, int i, int FLAG_dr0){
  
  // double sigma, dx_dsigma, sigma_2, dx_dsigma_2 ;
  // get_sigma(par, iDom, i, &sigma, &dx_dsigma);
  // sigma_2=sqr(sigma);
  // dx_dsigma_2=sqr(dx_dsigma);
  
  // double complex s, a2, a1, a0, da2, da2_dsig, da1, da0, A_f=0.;
  // complex_sigma_derivs f;

  // get_SigmaDerv_from_SpecDerv_complex(par, iDom, i, W, &f);
    
  // s = par.s;
  // func_alpha0(par, sigma, &a0, &da0);
  // func_alpha1(par, sigma, &a1, &da1);
  // func_alpha2(par, sigma, &a2, &da2_dsig, &da2);

  // if(FLAG_dr0==0) A_f = a2*f.d2sigma + a1*f.dsigma + a0*f.d0;  
  // else if(FLAG_dr0==1) A_f = da2*f.d2sigma + da1*f.dsigma + da0*f.d0;  
  // else {fprintf(stderr, "Error in function Diff_Operator\nFLAG_dr0 = %d not implementend\n", FLAG_dr0); exit(-1);}
  
  return 0.;
}
// -------------------------------------------------------------------------------
double complex Get_ModeSum_Source(parameters par, int iDom, int i){
  
  double complex HypSource = 0.;
  return HypSource;
}
// -------------------------------------------------------------------------------
double complex Get_Extended_Source(parameters par, complex_derivs W, int iDom, int i){
  
  double complex HypSource;
    
  HypSource = -Diff_Operator(par, W, iDom, i, 1);


  return HypSource;
}
//------------------------------------------------------------
void Alocate_Diff_Operator_Functions(parameters *par){
  int N = (*par).Ntotal;
  (*par).W = dvector(0,N);
  (*par).alpha_2 = dvector(0,N);
  (*par).alpha_1 = dvector(0,N);
  (*par).alpha_0 = dvector(0,N);
  (*par).beta_1 = dvector(0,N);
  (*par).beta_0 = dvector(0,N);
  

  return;
}
//------------------------------------------------------------
void Get_Diff_Operator_Functions(parameters *par){
  int i0, i1, idom, iF, N0, N1, N=(*par).Ntotal, II; 
  // func_hyperboloidal_coords w	, a2, a1, a0, b1, b0;
  double  sigma, dx, w, a2, a1, a0, b1, b0;

    
  fill0_dvector((*par).W, 0, N);    
  fill0_dvector((*par).alpha_2, 0, N);    
  fill0_dvector((*par).alpha_1, 0, N);    
  fill0_dvector((*par).alpha_0, 0, N);    
  fill0_dvector((*par).beta_1, 0, N);    
  fill0_dvector((*par).beta_0, 0, N);
  
  N0 = (*par).N0;
  for(idom=0; idom<nDoms; idom++){
    N1 = (*par).N1[idom];
    

    
      
    for(i1=0;i1<=N1; i1++){             
      get_sigma((*par), idom, i1, &sigma, &dx);
      a2 = func_alpha2((*par), sigma);
      a1 = func_alpha1((*par), sigma);
      a0 = func_alpha0((*par), sigma);
      b1 = func_beta1((*par), sigma);
      b0 = func_beta0((*par), sigma);
      w  = func_w((*par), sigma);

      for(i0=0;i0<=N0; i0++){
        for(iF = 0; iF<nFields; iF++){
          II = Index((*par), idom, iF, i1, i0);
          (*par).W[II] = w;
          (*par).alpha_2[II] = a2;
          (*par).alpha_1[II] = a1;
          (*par).alpha_0[II] = a0;
          (*par).beta_1[II]  = b1;
          (*par).beta_0[II]  = b0;

        }          
      }
    }
  }


  return;
}
//-----------------------------------------------
void free_Diff_Operator_Functions(parameters *par){
	int N = (*par).Ntotal;
 
  free_dvector((*par).W, 0, N);
  free_dvector((*par).alpha_2, 0, N);
  free_dvector((*par).alpha_1, 0, N);
  free_dvector((*par).alpha_0, 0, N);
  free_dvector((*par).beta_1, 0, N);
  free_dvector((*par).beta_0, 0, N);
  
  
	return;
}
//----------------------------
double particle_trajectory(parameters par, double tau){
//TO-DO: IMPLEMENTE CORRECT EXPRESSION FOR PARTICLES TRAJECTORY (24/09/24)
  double rp_over_M, r_initial = 1000, v=1.e-3;


  rp_over_M = r_initial + v * tau; //JUST A GENERIC FUNCTION


  return rp_over_M;
}

//----------------------------
double delta_func_Gaussian_approximation(parameters par, double tau, double sigma){
  double delta_gaussian, lambda_rStar = par.lambda_rStar, delta_rp, x, dx_dsigma, rStar_over_M, rStar_P, x_P, dxP_dsigma, rp, sigma_p, dsigmap_dr;


  func_tortoise_x(par,  sigma, &x, &dx_dsigma);
  rStar_over_M = par.lambda_over_M * x;


  rp = particle_trajectory(par, tau);
  func_sigma_of_r(par, rp, &sigma_p, &dsigmap_dr);

  func_tortoise_x(par,  sigma_p, &x_P, &dxP_dsigma);
  rStar_P = par.lambda_over_M * x_P;

  delta_rp = rStar_over_M - rStar_P;
  delta_gaussian = 1./ ( sqrt(2 * Pi ) * lambda_rStar)  * exp(  -0.5 * sqr( delta_rp/lambda_rStar )   )  ;

  //ATTENTION: ATTECIAPTING POSSIBLE PROBLEMS (24/09/24)
  // GRID HAS sigma = 0 (future null infintyt) and sigma = 1 (HORIZON)
  // TORTOISE COORDINATE x EXPLODES AT sigma = 0 and sigma = 1
  // BUT GAUSSIAN SHOULD GO TO ZERO AT sigma = 0 and sigma = 1

  return delta_gaussian;
}
//----------------------------
double complex Source_in_Schwarzschild_Coord(parameters par, double tau, double sigma){
  double complex Source, G=0;
//   // double delta, fp, fp_2, dfp_dr, rp ;
  
//   // rp = particle_trajectory(par, tau);
//   // func_f(rp, &fp, &dfp_dr);
//   // fp_2 = sqr(fp);
  
  
//   // G = - 4.*Pi*fp_2*cY/(E0*rp);
//   // //TO-DO LIST: 24/09/24: IMPLEMENT ENERGY E0, SPHERICAL HARMONIC;


//   // delta = delta_func_Gaussian_approximation(par, tau, sigma);
//   // Source = G * delta;
Source = G;
  return Source;
}