// source /home_local/nb-oma/Downloads/src_Intel/bin/compilervars.sh intel64
// source /home_local/ansorg/src/intel/bin/compilervars.sh intel64
// clear; icc -O2 *.c -o 1+1_Free_Boundaries

#include "utilities.h"
#include "spectral_utilities.h"

#define nDoms     4
#define nFields   2
#define nFields1D 0

#define Newton_itmin      2     
#define Newton_itmax      10     
#define Newton_tol        5.e-10
#define Newton_verb       1

#define bicgstab_decr     1.e-3    
#define bicgstab_verb     1 
#define bicgstab_itmax    10

#define Newton_SDIRK_itmin  2      
#define Newton_SDIRK_itmax  20      
#define Newton_SDIRK_tol    5.e-10 
#define Newton_SDIRK_verb   0

#define Newton_LinSolve_SDIRK_itmin  1      
#define Newton_LinSolve_SDIRK_itmax  5 
#define Newton_LinSolve_SDIRK_tol 5.e-10
#define Newton_LinSolve_SDIRK_verb 0

#define Newton_verb_ODE 0
#define Newton_tol_ODE  1.e-11




//Newton Constant of Gravitation [ in m^3/(kg s^2) ]
#define G_Grav  6.67545e-11 

//Polytropic Constant [ in m^5/(kg s^2)]
#define K_poly 4.5e-3

//Central density [ in kg/m^3]
#define rho_c 1.9891e18

typedef struct PARAMETERS{ 
	int N0, N1[nDoms], mA[nDoms+1],  Nmax,
		SDIRK_s, nmax, IG_FLAG, Solver_FLAG, EquationType_FLAG,
		ntotal_2Dslice,//nVslice, //TOTAL NUMBER OF SPACE GRID POINTS AT T=CONSTANT SLICE FOR THE 1+1D FIELDS
		ntotal_slice, //TOTAL NUMBER OF GRID POINTS AT A SLICE FOR THE 1+1D FIELDS PLUS 1D FIELD
		ntotal_2D, //nV2D, //TOTAL NUMBER OF SPACETIME GRID POINTS FOR THE 1+1D FIELDS
		ntotal_1D, //TOTAL NUMBER OF TIME GRID POINTS FOR THE 1D FIELDS
		// Nmax, nX2D, nXslice, 
		ntotal,// nVtotal, //TOTAL NUMBER OF SPACETIME GRID POINTS FOR THE 1+1D FIELDS PLUS 1D FIELDS
		Ntotal,
		// *iX_of_iV, *iV_of_iX, 
		N_ODE,
		*indx_LU;

	double *grid_points_x0, *grid_points_x1[nDoms],
			angle_incr_x0, angle_incr_x1[nDoms],
			Tmin, Tmax, DT, t0, t1, Del_t,
			 **Dx0, **Dx1[nDoms], bicgstab_tol,
			 **D0_2D, **D1_2D, **D11_2D, **D01_2D,
			 **Id_2D, **Dtau_2D, **Dsig_2D, **Dsigsig_2D, **Dtausig_2D,
			 *W, *alpha_2, *alpha_1, *alpha_0, *beta_1, *beta_0, 
			**SDIRK_A, *SDIRK_B, *SDIRK_C, gamma, gah,
			**J_LU;
		
	char	SimName[1000], grid_x0[50], grid_x1[nDoms][50],  Equation[50];

	FILE *fout;
	
	double AnMR_kappa[nDoms], AnMR_x_boundary[nDoms];

	double complex s, ss, Amp_RegularMode, Amp_MirrorMode;

	int ell, spin, m,  n_QNM, i_omp, parity, HomoEq;

	double Mass, rh_over_M, rho_0, rho_1,
		   lambda_over_M, sigma[nDoms+1],
		   *cheb_real_phi_ID[nDoms], *cheb_imag_phi_ID[nDoms],
		   //*cheb_real_dphi_dtau_ID[nDoms], //, *cheb_imag_dphi_dtau_ID[nDoms],
		   *cheb_real_Source[nDoms], *cheb_imag_Source[nDoms];
	
	int N_input_QNM[nDoms], N_input_Source[nDoms];
	double prec_input_QNM;
	
	int flag_NewCoord, flag_debug;
	
	int ncount;

	double lambda_rStar, ro_over_M;

	double epsilon, a0_over_lambda;


} parameters;

typedef struct DERIVS{ 
	double *X, *X1, *X11, *Y, *Y1, *Y11;
} derivs;

typedef struct CHEB_DERIVS{ 
	double **X, **X1, **X11, **Y, **Y1, **Y11;
} cheb_derivs;

typedef struct DERIVS_ODE{ 
	double *d0, *d1, *d11;
} derivs_ODE;

typedef struct COMPLEX_DERIVS_Space{
	double complex d0, d1, d11;

} complex_derivs;

typedef struct COMPLEX_SIGMA_DERIVS_Space{
	double complex d0, dsigma, d2sigma;

} complex_sigma_derivs;

// Subroutines in "derivatives.c"
void Get_DiffMatrices_2D(parameters *par);
void Allocate_PhysDiffMatrices_2D(parameters *par);
void Get_PhysDiffMatrices_2D(parameters *par);
void free_DiffMatrices_2D(parameters *par);
void Derivatives_x1(parameters par, derivs Wslice);
void Derivatives_x0x1(parameters par, derivs V);

// Subroutines in "parameters.c"
void set_parameters(parameters *par, int spin, int l, int m, double eps, double a0, double ID_r0_over_M, double ID_Gauss_width_lambda, int N_gaussian, int N_bump);
void free_parameters(parameters *par);

// Subroutines in "solve_equations.c"
void solve_equations(parameters par, derivs ID, double *X, derivs Sol);
void get_solution(parameters par, derivs ID, double *X, derivs W);

// Subroutines in "Get_Differentiation_Matrices"
void Get_Differentiation_Matrices(parameters *par);
void Free_Differentiation_Matrices(parameters *par);

// Subroutines in "allocate_and_free_derivs.c"
void allocate_derivs(derivs *V, int nl, int nh);
void allocate_derivs_NULL(derivs *V);
void fill0_derivs(derivs V, int nl, int nh);
void copy_derivs(derivs Vout, derivs Vin, int nl, int nh);
void free_derivs(derivs *V, int nl, int nh);
void allocate_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch);
void fill0_cheb_derivs(cheb_derivs V, int nrl, int nrh, int ncl, int nch);
void copy_cheb_derivs(cheb_derivs Vout, cheb_derivs Vin, int nrl, int nrh, int ncl, int nch);
void free_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch);

// Subroutines in "Coordinates.c"
double get_grid_angle_incremente(char *grid, int N);
double get_grid_point(char *grid, int i, int N);
double get_x_from_chi(double x_boundary, double kappa, double chi);
void get_sigma(parameters par, int iDom, int i1, double *sigma, double *dx_dsigma);
double get_t(parameters par, int i0);
double get_sigma_from_x1(parameters par, int iDom, double x);
double get_x1_from_sigma(parameters par, int iDom, double sigma);
double get_t_from_x0(parameters par, double x0);
void scale_derivs_x0_to_tau(parameters par, derivs W1D_x0, derivs W1D_tau, int nsize_1D);
void scale_derivs_x0x1_to_tau_sigma(parameters par, int iDom, int i1, derivs W_x0x1, derivs W_tau_sigma, int nsize);
double get_z_from_chi(parameters par, double chi);
double get_chi_from_z(parameters par, double z);

// Subroutines in "IndexRoutines.c"
int Index_i(parameters par, int iDom, int i);
int Index_j(parameters par, int iDom, int j);
int Index_slice(parameters par, int iDom, int iField, int i1);
int Index_field1D_slice(parameters par, int iField);
int Index_field1D(parameters par, int iField, int i0);
int Index(parameters par, int iDom, int iField, int i1, int i0);
void Get_Index_Arrays(parameters *par);
void Get_Indices_From_iV(parameters par, int iV, int *iDom, int *iField, int *i, int *j);
void Get_Indices_From_iVslice(parameters par, int iVslice, int *iDom, int *iField, int *i);
void Free_Index_Arrays(parameters *par);

// Subroutines in "Get_Initial_Guess.c"
void get_initial_guess(parameters par, derivs ID, double *X);
void get_Initial_Guess_SDIRK(parameters par, derivs ID, double *X);
void SDIRK_Step_IG(parameters par, double B, double h, double *Xn, double *Xnp1);
void get_Ws_from_M_and_k(parameters par, double *Mj, double *kj, derivs Wslice, derivs W1D);
void Get_Fields2D_Slice_atGrid(parameters par, int iDom, int i1, derivs Wslice, derivs W);
void F_of_kj_IG(parameters par, double x0, double *Mj, double *kj, double *F);
// void J_times_Dkj_IG(parameters par, double B, double *Dkj, derivs Wslice, derivs W1D, double *JDkj);
void J_times_Dkj_IG(parameters par, double x0, double *Mj, double *kj, double *Dkj, double *JDqj);

// Subroutines in "newton_SDIRK_IG.c"
int Newton_SDIRK_IG(parameters par, double x0, double *Mj, double *kj);
// void Jacobian_SDIRK_IG(parameters par, double B, derivs Wslice, derivs W1D, double **J);
void Jacobian_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J);
// void DF_of_X_SDIRK_IG(parameters par, double B, double *Mj, double *kj, derivs Wslice, derivs W1D, double **J);
void DF_of_X_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J);


// Subroutines in "FieldAndBoundEqns.c"
void LinearOperator_Matrix(parameters par, int idom,  int i0, int i1, double **M);
void LinearBoundaryConditions_Matrix(parameters par, int idom,  int i0, int i1, double **M);

double window_function(parameters par, double t);
double jump_field(parameters par, int iF, double x0);
double jump_field_dtau(parameters par, int iF, double x0);
double jump_derivative(parameters par, int iF, double x0);
double jump_derivative_dtau(parameters par, int iF, double x0);
double complex Get_Source(parameters par, int iDom,  double x0, int i1);
void F_of_XY(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, double *F);
void JDX_of_XY(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *J);
void Boundary_Data_Left(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, double *F);
void LinearBoundary_Data_Left(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *F);

void Boundary_Data_Right(parameters par, int iDom, double x0, int x1, derivs W, derivs W1D, double *F);
void LinearBoundary_Data_Right(parameters par, int iDom, double x0, int i1, derivs W, derivs W1D, derivs DW, derivs DW1D, double *F);


void FBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, double *F);
void JDXBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, derivs DWslice, derivs DW1D, double *J);

// Subroutines in "FuncAndJacobian.c"
void Get_V_From_X(parameters par, double *X, derivs V);
void Get_Fields1D(parameters par, int i0, derivs ID, derivs V, derivs W1D);
void Get_Fields2D_atGrid(parameters par, int iDom, int i0, int i1, derivs ID, derivs V, derivs W_grid);
void Get_Fields2D_atSlice(parameters par, int i0, int i1, derivs ID, derivs V, derivs Wslice);
void F_of_X(parameters par, double *X, derivs V, derivs ID, double *F);
void J_times_DX(parameters par, double *DX, derivs V, derivs ID, double *JDX);
void Get_LinearOperator(parameters par, double **M);


// Subroutines in "newton_direct.c"
void Jacobian(parameters par, derivs V, derivs ID, double **J);
void DF_of_X(parameters par, double *X, derivs V, derivs ID, double **J);
int newton_direct(parameters par, derivs ID, double *X);
void Allocate_LinearStatic_Jacobian(parameters *par);
void Free_LinearStatic_Jacobian(parameters *par);

// Subroutines in "get_Initial_Data.c"
double delta(double a, double b);
double ID_Function(parameters par, double sigma);
void get_Initial_Data(parameters par, derivs ID);
void get_Initial_Data_newTimeStep(parameters par, derivs W, derivs ID);
void Output_Initial_Data(parameters par, derivs ID);
void Output_Cheb_Initial_Data(parameters par,derivs ID);

// Subroutines in "Butcher_Tableau.c"
void Butcher_Tableau(parameters *par);
void free_Butcher_Tableau(parameters *par);

// Subroutines in "newton.c"
// void get_cheb(parameters par, int indx_slice, double *p, double *V, double **chebV);
int newton(parameters par, derivs ID, double *X);
int bicgstab(parameters par, derivs V, derivs ID, double *F, double *DX, double *normres);
void get_Chebyshev_Coefficients_derivs(derivs psi, derivs c, int N, char *grid);
void get_chebV(parameters par, derivs V, cheb_derivs chebV);

// Subroutines in "SDIRK_Preconditioner.c"
void get_chebF(parameters par, double *F, double **chebF);
void Get_Wslice_W1D_Fslice(parameters par, double x0, derivs ID, double **chebF, cheb_derivs chebV, derivs Wslice, derivs W1D, double *Fslice);
void F_of_kj(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F);
void LinSolve_SDIRK(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj);
void SDIRK_Step(parameters par, derivs ID, cheb_derivs chebV, double **chebF, double x0, double h, double *X_previous, double *X_next);
void PreCond(parameters par, cheb_derivs chebV, double *F, derivs ID, double *X);

//Routines in "hyperboloidal_functions.c"
double get_r_over_M_from_rStar(parameters par, double rStar);
void func_rho(parameters par, double sigma, double *rho, double *drho_dsigma, double *d2rho_dsigma2);
void func_beta(parameters par, double sigma, double *beta, double *dbeta_dsigma);
void func_r_of_sigma(parameters par, double sigma, double *r, double *dr_dsigma);
void func_sigma_of_r(parameters par, double r, double *sigma, double *dsigma_dr);
void func_f(double r, double *f, double *df_dr);
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma);
void func_height(parameters par, double sigma, double *h, double *dh_dsigma);
void func_Height(parameters par,double sigma, double *H, double *dH_dsigma);
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma);
void func_Z(parameters par, double sigma, double complex *Z, double complex *dlnZ_dsigma);
double func_alpha2(parameters par, double sigma);
double func_alpha1(parameters par, double sigma);
double func_alpha0(parameters par, double sigma);
double func_beta1(parameters par, double sigma);
double func_beta0(parameters par, double sigma);
double func_w(parameters par, double sigma);
void func_p(parameters par, double sigma, double *p, double *dp_dsigma);
void func_gamma(parameters par, double sigma, double *gamma, double *dgamma_dsigma);
void func_w_RWZ(parameters par, double sigma, double *w);
void func_V(parameters par, double sigma, double *V);
void func_V_bump(parameters par, double sigma, double *V);
double complex Diff_Operator(parameters par, complex_derivs W, int iDom, int i, int FLAG_dr0);
double complex Get_ModeSum_Source(parameters par, int iDom, int i);
double complex Get_Extended_Source(parameters par, complex_derivs W, int iDom, int i);
void Alocate_Diff_Operator_Functions(parameters *par);
void Get_Diff_Operator_Functions(parameters *par);
void free_Diff_Operator_Functions(parameters *par);

//Routines in "SelfForce_functions.c"
double get_clm(int l, int m);
void get_SphericalHarmonicPhase(parameters par, double t, double complex *cY, double complex *dcY_dtheta);
void get_gotC_lm(parameters par, double *gotC_lm, double *d_gotC_lm_dr);
double complex SphericalHarmonicPhase(parameters par, double t);
double complex Source_G(parameters par, double t);
double complex Source_F(parameters par, double t);
double complex Get_phi_From_hyp_phi(parameters par, double sigma, double complex hyp_phi);
double complex Get_dphi_dr_From_dhyp_phi_dsigma(parameters par, double sigma, double complex hyp_phi, double complex dhyp_phi_dsigma, double complex dhyp_phi_dtau);
double get_Schwarzschild_t(parameters par, double tau, double sigma);
void get_Ft_SSF_lm(parameters par, derivs Sol, int i0, double complex *Ft_SSF_lm_dom0, double complex *Ft_SSF_lm_dom1);
void get_Flux_lm(parameters par, derivs Sol, int i0, double *Flux_lm_hrzn, double *Flux_lm_scri, double *Flux_lm);
double get_Ft_FLux_lm(parameters par, derivs Sol, int i0);
void get_Fr_SSF_lm(parameters par, derivs Sol, int i0, double complex *Fr_SSF_lm_dom0, double complex *Fr_SSF_lm_dom1);

// Subroutines in "output.c"
void output_cheb(parameters par, derivs Sol);
void output_cheb_TimeEvolution(parameters par, derivs ID, derivs Sol);
void output_Sol_particle(parameters par, derivs Sol);
void output_Sol(parameters par, derivs Sol);
void output_TimeEvolution(parameters par, derivs ID, derivs Sol);
void output_Sol_derv_sigma(parameters par, derivs Sol);
void output_Fluxes_and_Ft(parameters par, derivs Sol);
void output_Fr(parameters par, derivs Sol);
void output_Source(parameters par);


// Subroutines in "load_data.c"
void read_header(parameters *par, FILE *fr, double complex *s, int *N);
void load_QNM_ID(parameters *par);
void load_Source(parameters *par);
void free_QNM_ID(parameters *par);
void free_Source(parameters *par);

// Subroutines in "problem_specific_functions.c"
void load_QNMSet(parameters par, char *fn, double complex *qnm_s, double complex *qnm_phi, double complex *qnm_psi);
double ID_Functions(double sigma);
double dID_Functions(double sigma);


// Subroutines in "debug.c"
void pause();
void PrintVector(double *J, int na, int nb);
void PrintMatrix(char *fn, double **J, int la, int lb, int ca, int cb);
void PrintMatrix_to_File(parameters par, char *fn, double **J, int la, int lb, int ca, int cb);
void Check_OS();

// Subroutines in "self_force_functions.c"
long double acoth(double x);
long double exact_solution_circular_orbit_static_mode(parameters par, int iDom, long double sigma);



// Subroutines in "1+1_2_Domains_SDIRK.c"
int main();

