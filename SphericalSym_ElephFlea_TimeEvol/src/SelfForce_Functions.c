#include "1+1_Free_Boundaries.h"

//----------------------------
double get_clm(int l, int m){
  int x=0;
  double c_lm = 1.;

  while(x<= 2*m-1){

    c_lm/=sqrt(l + m - x);
    x++;
  }
  c_lm *=sqrt( (2.*l +1.)/(4*Pi));

  return c_lm;
}
//----------------------------

double Azimutal_phase_particle(parameters par, double t){
  double phi_p,  Omega=par.M_Omega0;

  phi_p = Omega*t;

  return phi_p;
}
//----------------------------
void get_SphericalHarmonicPhase(parameters par, double t, double complex *cY, double complex *dcY_dtheta){
	double complex Y, phi_p;
	double l = 1.*par.ell, m = 1.*par.m, 
          c_lm = get_clm(l, m), 
          P_lm = plgndr(l, m, 0.),
          P_lMINUSONE_m = plgndr(l-1, m, 0.),
          dP_lm_dx = (l+m)*P_lMINUSONE_m,
          dP_lm_dtheta = -dP_lm_dx;	

	phi_p = Azimutal_phase_particle(par, t);
	

  *cY           = c_lm * P_lm * cexp(- I * m * phi_p);
  *dcYlm_dtheta = c_lm * dP_lm_dtheta * cexp(- I * m * phi_p);

	
	return;
}
//----------------------------
void get_gotC_lm(parameters par, double complex *gotC_lm, double complex *d_gotC_lm_dr){
  double f0 = par.f0, L0 = par.L0_over_M, r0 = par.r0_over_M, r0_2 = sqr(r0), r0_3 = r0*r0_2, r0_4 = r0_3*r0, l =  1.*par.ell,
         cl     =   -16.*Pi*f0*L0/(  l*(l+1.) * r0_3    ),
         dcl_dr =   -16.*Pi*f0*L0*(-3.)/(  l*(l+1.) * r0_4    );

  get_SphericalHarmonicPhase(par, t, &cY, &dcY_dtheta);
  
  *gotC_lm = cl * dcY_dtheta;
  *d_gotC_lm_dr = dcl_dr * dcY_dtheta;
  return;
}
//----------------------------
double complex Source_G(parameters par, double t){
  double complex G=0., cY, dcY_dtheta, gotC_lm, d_gotC_lm_dr;
  double f0, df0_dr, f0_2=sqr(f0), E0 = par.E0, r0=par.r0_over_M, r0_2=sqr(r0);

  func_f(r0, &f0, &df0_dr);  
  get_SphericalHarmonicPhase(par, t, &cY, &dcY_dtheta);
  get_gotC_lm(par, &gotC_lm, &d_gotC_lm_dr);

  if(par.spin == 0){
		G = - 4.*Pi*f0_2*cY/(E0*r0);
	}
	else if (par)spin == 2){
    if ( par.parity==1 ){
      sprintf((*par).Equation, "ReggeWheeler");
      G = f0*r0_2*d_gotC_lm_dr + f0*r0*gotC_lm - (df0_dr*r0_2*gotC_lm + f0*2*r0*gotC_lm + f0*r0_2*d_gotC_lm_dr);
    }
	//else{
	// 		sprintf((*par).Equation, "Zerilli"); 
	// 	}
	}
	else{
		fprintf(stderr, "Error in set_parameters: source for spin = %d not implemented\n", par.spin );
		exit(-1);
	}


  return G;
}
//----------------------------
double complex Source_F(parameters par, double t){
  double complex F=0.;
  double f0 = par.f0, r0 = par.r0_over_M, r0_2 = sqr(r0), gotC_lm, d_gotC_lm_dr; 

  get_gotC_lm(par, &gotC_lm, &d_gotC_lm_dr);
  
  if(par.spin == 0){
		F = 0.;
	}
	else if (par)spin == 2){
    if ( par.parity==1 ){
      sprintf((*par).Equation, "ReggeWheeler");
      F = f0 * r0_2 * gotC_lm;
    }
  //else{
  //    sprintf((*par).Equation, "Zerilli"); 
  //  }
	}
	else{
		fprintf(stderr, "Error in set_parameters: source for spin = %d not implemented\n", par.spin );
		exit(-1);
	}


  return F;
}


//----------------------------------------------------------------------------------------
double complex Get_phi_From_hyp_phi(parameters par, double sigma, double complex hyp_phi){
 //MAPPING IS ONLY VALID FOR FREQUECY DOMAIN
  double complex phi, Z, dlnZ;
  func_Z(par, sigma, &Z, &dlnZ);

  phi = Z*hyp_phi;


  return phi;
}
//----------------------------------------------------------------------------------------
double complex Get_dphi_dr_From_dhyp_phi_dsigma(parameters par, double sigma, double complex hyp_psi, double complex dhyp_psi_dsigma, double complex dhyp_psi_dtau){

  double r, dr_dsigma, dsigma_dr, H, dH;
  double complex dphi_dr, Z, dlnZ_dsigma;
  func_Z(par, sigma, &Z, &dlnZ_dsigma);
  func_Height(par,sigma, &H, &dH);

  func_r_of_sigma(par, sigma, &r, &dr_dsigma);
  dsigma_dr = 1./dr_dsigma;
  dphi_dr = dsigma_dr * (dhyp_psi_dsigma + dH*dhyp_psi_dtau) ;
           

  return dphi_dr;
}
//---------------------------------
double get_Schwarzschild_t(parameters par, double tau, double sigma){
  double lambda = par.lambda_over_M, H, dH, t;

  func_Height(par,sigma, &H, &dH);

  t = lambda*(tau - H);

  return t;
}
//---------------------------------
void get_Fr_SSF_lm(parameters par, derivs Sol, int i0, double complex *Fr_SSF_lm_dom0, double complex *Fr_SSF_lm_dom1){
  double q = par.q,  Omega0 = par.M_Omega0, sigma0=par.sigma0,
         l=1.*par.ell, m=1.*par.m, tau, t, r0 = par.r0_over_M,      
         c_lm = sqrt( (2.* l+1.)*factrl(l-m)/(4*Pi*factrl(l+m))),
         P_lm = plgndr(l, m, 0.);

 
  tau = get_t(par, i0); //BAD NAME FOR FUNCTION: HERE IT CALCULATES THE HYPERBOLOIDAL TIME TAU
  t = get_Schwarzschild_t(par, tau, sigma0);

  double complex dphi_dr_dom0, dphi_dr_dom1,
                 hyp_phi_dom0, hyp_phi_dom1, 
                 dhyp_phi_dsigma_dom0, dhyp_phi_dsigma_dom1,
                 dhyp_phi_dtau_dom0, dhyp_phi_dtau_dom1, 
                 Y_lm = c_lm*P_lm*cexp(I*m*Omega0*t);

  derivs Sol_x0x1, Sol_sigtau;
  allocate_derivs(&Sol_x0x1, 0, 1);
  allocate_derivs(&Sol_sigtau, 0, 1);
  int idom, i1_part, iF, I_phi_particle;


  idom = 0; i1_part = 0; 

    iF = 0; //Re(phi)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    hyp_phi_dom0 = Sol_sigtau.X[0]; 
    dhyp_phi_dsigma_dom0 = Sol_sigtau.X1[0];

    iF = 1; //Re(phi)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    hyp_phi_dom0 += I*Sol_sigtau.X[0]; 
    dhyp_phi_dsigma_dom0 += I*Sol_sigtau.X1[0];

    iF = 2; //Re(phi,tau)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    dhyp_phi_dtau_dom0 = Sol_sigtau.X[0];

    iF = 3; //Re(phi,tau)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    dhyp_phi_dtau_dom0 += I*Sol_sigtau.X[0];

    dphi_dr_dom0 = Get_dphi_dr_From_dhyp_phi_dsigma(par, sigma0, hyp_phi_dom0, dhyp_phi_dsigma_dom0, dhyp_phi_dtau_dom0);

  idom = 1; i1_part = par.N1[1]; 

    iF = 0; //Re(phi)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    hyp_phi_dom1 = Sol_sigtau.X[0]; 
    dhyp_phi_dsigma_dom1 = Sol_sigtau.X1[0];

    iF = 1; //Re(phi)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    hyp_phi_dom1 += I*Sol_sigtau.X[0]; 
    dhyp_phi_dsigma_dom1 += I*Sol_sigtau.X1[0];

    iF = 2; //Re(phi,tau)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    dhyp_phi_dtau_dom1 = Sol_sigtau.X[0];

    iF = 3; //Re(phi,tau)
    I_phi_particle  = Index(par, idom, iF, i1_part, i0);
    Sol_x0x1.X[0] = Sol.X[I_phi_particle];     Sol_x0x1.Y[0]   = Sol.Y[I_phi_particle];
    Sol_x0x1.X1[0] = Sol.X1[I_phi_particle];   Sol_x0x1.Y1[0]  = Sol.Y1[I_phi_particle];
    Sol_x0x1.X11[0] = Sol.X11[I_phi_particle]; Sol_x0x1.Y11[0] = Sol.Y11[I_phi_particle]; 
    scale_derivs_x0x1_to_tau_sigma(par, idom, i1_part, Sol_x0x1, Sol_sigtau, 1);
    dhyp_phi_dtau_dom1 += I*Sol_sigtau.X[0];  
    
    dphi_dr_dom1 = Get_dphi_dr_From_dhyp_phi_dsigma(par, sigma0, hyp_phi_dom1, dhyp_phi_dsigma_dom1, dhyp_phi_dtau_dom1);


  *Fr_SSF_lm_dom0 = -q * (dphi_dr_dom0 - hyp_phi_dom0/r0)*Y_lm/r0;
  *Fr_SSF_lm_dom1 = -q * (dphi_dr_dom1 - hyp_phi_dom1/r0)*Y_lm/r0;

  free_derivs(&Sol_x0x1, 0, 1);
  free_derivs(&Sol_sigtau, 0, 1);
  return;
}
//----------------------------------------------------------------------------------------
void get_Ft_SSF_lm(parameters par, derivs Sol, int i0, double complex *Ft_SSF_lm_dom0, double complex *Ft_SSF_lm_dom1){
  double q = par.q,  lambda=par.lambda_over_M, Omega0 = par.M_Omega0, sigma0=par.sigma0,
         l=1.*par.ell, m=1.*par.m, tau, t, r0 = par.r0_over_M,      
         c_lm = sqrt( (2.* l+1.)*factrl(l-m)/(4*Pi*factrl(l+m))),
         P_lm = plgndr(l, m, 0.);

 
  tau = get_t(par, i0); //BAD NAME FOR FUNCTION: HERE IT CALCULATES THE HYPERBOLOIDAL TIME TAU
  t = get_Schwarzschild_t(par, tau, sigma0);

  double complex dv_dt_dom0, dv_dt_dom1, Y_lm = c_lm*P_lm*cexp(I*m*Omega0*t);

  int I_real_particle_dom0 = Index(par, 0, 2, 0, i0), //Index(par, 0, 3, 0, i0), 
      I_imag_particle_dom0 = Index(par, 0, 3, 0, i0), //Index(par, 0, 4, 0, i0), 
      I_real_particle_dom1 = Index(par, 1, 2, par.N1[1], i0), //Index(par, 1, 3, par.N1[1], i0), 
      I_imag_particle_dom1 = Index(par, 1, 3, par.N1[1], i0);//Index(par, 1, 4, par.N1[1], i0);

      dv_dt_dom0 = Sol.X[I_real_particle_dom0] + I*Sol.X[I_imag_particle_dom0];
      dv_dt_dom1 = Sol.X[I_real_particle_dom1] + I*Sol.X[I_imag_particle_dom1];
      

  *Ft_SSF_lm_dom0 = q*dv_dt_dom0 * Y_lm/(lambda*r0);
  *Ft_SSF_lm_dom1 = q*dv_dt_dom1 * Y_lm/(lambda*r0);

  return;
}
//----------------------------------------------------------------------------------------
void get_Flux_lm(parameters par, derivs Sol, int i0, double *Flux_lm_hrzn, double *Flux_lm_scri, double *Flux_lm){
  // int I_phi_real_scri = Index(par, 0, 0, par.N1[0], i0), 
  //     I_phi_imag_scri = Index(par, 0, 1, par.N1[0], i0),
  //     I_phi_real_hrz  = Index(par, 1, 0, 0, i0), 
  //     I_phi_imag_hrz  = Index(par, 1, 1, 0, i0),

     int I_psi_real_scri = Index(par, 0, 2, par.N1[0], i0), 
      I_psi_imag_scri = Index(par, 0, 3, par.N1[0], i0),
      I_psi_real_hrz  = Index(par, 1, 2, 0, i0), 
      I_psi_imag_hrz  = Index(par, 1, 3, 0, i0);

      
  double lambda = par.lambda_over_M;
  double complex dv_dt_scri, dv_dt_hrz;


  dv_dt_scri = Sol.X[I_psi_real_scri] + I*Sol.X[I_psi_imag_scri];
  dv_dt_hrz = Sol.X[I_psi_real_hrz] + I*Sol.X[I_psi_imag_hrz];

  *Flux_lm_hrzn = sqr(cabs(dv_dt_hrz))/(4*Pi*sqr(lambda));   ///(4*Pi*sqr(2*lambda));
  *Flux_lm_scri = sqr(cabs(dv_dt_scri))/(4*Pi*sqr(lambda));   ///(4*Pi*sqr(2*lambda));  

  *Flux_lm = *Flux_lm_hrzn + *Flux_lm_scri;
  return;
}
//----------------------------------------------------------------------------------------
double get_Ft_FLux_lm(parameters par, derivs Sol, int i0){
  double Flux_lm_hrzn, Flux_lm_scri, Flux_lm, ut , Ft_lm_flux;

  get_Flux_lm(par, Sol, i0, &Flux_lm_hrzn, &Flux_lm_scri, &Flux_lm);
  
  ut=par.E0/par.f0; 
  
  Ft_lm_flux = Flux_lm * ut;

  return Ft_lm_flux;
}