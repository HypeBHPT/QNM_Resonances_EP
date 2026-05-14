#include "1+1_Free_Boundaries.h"

void set_parameters(parameters *par, int spin, int l, int m, double eps, double a0, double ID_r0_over_M, double ID_Gauss_width_lambda, int N_gaussian, int N_bump){

	//INPUT PHYSICAL PARAMETERS------------------------------------
	
	(*par).HomoEq = 0; //HomoEq = 0 (Vanishing source, non-vanishing ID), //HomoEq = 1 (non-vanishin source, vanishing ID)
	
	(*par).rh_over_M = 2.; //BH Horizon in units of M
	(*par).Mass = (*par).rh_over_M/2.;
	(*par).lambda_over_M = (*par).rh_over_M; //Hyperboloidal Length Scale in units of M

	(*par).rho_0 = (*par).rh_over_M/(*par).lambda_over_M; //Hyperboloidal Radial Gauge
	(*par).rho_1 = 0. ; //Hyperboloidal Radial Gauge
	
	(*par).spin = spin; //Field's spin weight (spin = 0 : Scalar field)
	(*par).ell = l; //Angular Mode
	(*par).m = m; //Angular Mode
	(*par).n_QNM = 0; //Individual QNM (CURRENTLY NOT BEING USED 03/09/24)

	//INPUT NUMERICAL FOR BUMP POTENTIAL------------------------------------
	(*par).epsilon = eps;
	(*par).a0_over_lambda = a0; 
	
	double rStar_0 = (*par).a0_over_lambda * (*par).lambda_over_M, r0_over_M, sigma_bump, dsigma0_dr;
	r0_over_M = get_r_over_M_from_rStar(*par, rStar_0);
	printf("r*_bump = %f M, r_bump = %f M\n", rStar_0* (*par).Mass, r0_over_M*(*par).Mass);
	
	func_sigma_of_r(*par, r0_over_M, &sigma_bump, &dsigma0_dr);
	

	//INPUT NUMERICAL PARAMETERS FOR GAUSSIAN ID------------------------------------
	(*par).lambda_rStar = ID_Gauss_width_lambda;
	(*par).ro_over_M = ID_r0_over_M;
	double sigma_0, dsigma_dr;
	func_sigma_of_r(*par, (*par).ro_over_M, &sigma_0, &dsigma_dr);
	printf("r0 = %f M , sigma_0 = %f \n", (*par).ro_over_M*(*par).Mass, sigma_0);
	printf("r0 = %f rh , sigma_0 = %f \n", (*par).ro_over_M / (*par).rh_over_M, sigma_0);
	exit(-1);

	int idom_bump, idom_gauss;
	if(sigma_0<sigma_bump){
		idom_gauss = 0;
		idom_bump  = 2;
	}
	else if(sigma_0>sigma_bump){
		idom_gauss = 2;
		idom_bump  = 0;
	}
	else{
		fprintf(stderr, "Error in funcion set_parameters: sigma_bump = %lf, sigma_gauss=%f conflicting\n",sigma_bump, sigma_0);
		exit(-1);
	}

	//DERIVED PHYSICAL PARAMETERS------------------------------------
	//Space

	(*par).sigma[0]    = 0.; //Location of Scri
	
	(*par).sigma[idom_gauss+1]  = sigma_0; //Location of Gauss peak
	
	(*par).sigma[2]    = 0.5*(sigma_bump+sigma_0); //Location of Bump
	
	(*par).sigma[idom_bump+1]    = sigma_bump; //Location of Bump
	
	(*par).sigma[nDoms]=1.; //Location of Horizon

	printf("sigma_bump = %10f, sigma_gauss=%10f\n",sigma_bump, sigma_0);
	exit(-1);
	
	
	(*par).parity = ( (*par).ell+(*par).m ) %2 ;

	sprintf((*par).Equation, "ReggeWheeler");
	// if((*par).spin == 0){
	// 		sprintf((*par).Equation, "ReggeWheeler"); //Zerilli; ;
	// 	}
	// 	else if((*par).spin == 2){
	// 		if( (*par).parity==0 ){
	// 			sprintf((*par).Equation, "Zerilli"); //Zerilli; ;
	// 		}
	// 		else{
	// 			sprintf((*par).Equation, "ReggeWheeler"); //Zerilli; ;
	// 		}
	// 	}
	// 	else{
	// 		fprintf(stderr, "Error in set_parameters: source for spin = %d not implemented\n", (*par).spin );
	// 		exit(-1);
	// }



	// (*par).N_input_QNM[0] = 150;
	// (*par).N_input_Source[0] = 200;
	// (*par).prec_input_QNM=159.546;

	

	//INPUT TIME STRIPS PARAMETERS------------------------------------
	(*par).Tmin = 0.; //Start Time
	(*par).Tmax = 500.; //Final Time
	(*par).nmax = 1000; //Division of Total Time interval
	(*par).DT = ((*par).Tmax-(*par).Tmin)/(*par).nmax;

	//INPUT NUMERICAL RESOLUTION------------------------------------
	int N0 = 5;


			
	//SET NUMERICAL RESOLUTION------------------------------------
	//TIME DIRECTION
	(*par).N0    = N0; //Time Resolution Domain
	sprintf((*par).grid_x0,"Radau_RHS");


	//SPACE DIRECTION
	int iDom;
	//DOM 0
	iDom = idom_gauss ;
	(*par).N1[iDom]  = N_gaussian; //Spatial Resolution Domain 0
	sprintf((*par).grid_x1[iDom],"Lobatto");
	(*par).AnMR_x_boundary[iDom]= 1.;
	(*par).AnMR_kappa[iDom] = 2.5;

	//DOM 1
	iDom = idom_gauss+1 ;
	(*par).N1[iDom]  = N_gaussian; //Spatial Resolution Domain 0
	sprintf((*par).grid_x1[iDom],"Lobatto");
	(*par).AnMR_x_boundary[iDom]= -1.;
	(*par).AnMR_kappa[iDom] = 2.5;//-log10(1/sqr((*par).lambda_rStar));

	//DOM 2
	iDom = idom_bump ;
	(*par).N1[iDom]  = N_bump; //Spatial Resolution Domain 0
	sprintf((*par).grid_x1[iDom],"Lobatto");
	(*par).AnMR_x_boundary[iDom]= -1.;
	(*par).AnMR_kappa[iDom] = 0;//-log10(1/sqr((*par).lambda_rStar));

	//DOM 3
	iDom = idom_bump+1;
	(*par).N1[iDom]  = N_bump; //Spatial Resolution Domain 0
	sprintf((*par).grid_x1[iDom],"Lobatto");
	(*par).AnMR_x_boundary[iDom]= -1.;
	(*par).AnMR_kappa[iDom] = 0;//-log10(1/sqr((*par).lambda_rStar));


	

	if((*par).HomoEq == 0 ){
		// sprintf( (*par).SimName, "CteField/spin_%d/l_%d/epsilon_%f_a0_%f/N0_%d_Ngauss_%d_Ngauss_%d_ndom_%d_nfields_%d",  (*par).spin,(*par).ell, (*par).epsilon, (*par).a0_over_lambda, (*par).N0, N_gaussian, N_bump ,nDoms, nFields);
		sprintf( (*par).SimName, "Gauss_ID_r0_over_M_%f_delta_%f/spin_%d/l_%d/epsilon_%f_a0_%f/N0_%d_Ngauss_%d_Ngauss_%d_ndom_%d_nfields_%d", (*par).ro_over_M, (*par).lambda_rStar, (*par).spin,(*par).ell, (*par).epsilon, (*par).a0_over_lambda, (*par).N0,  N_gaussian, N_bump ,nDoms, nFields);
			
		
	}
	else if((*par).HomoEq == 1 ){
		sprintf((*par).SimName, "Source_Evol" );
		// sprintf((*par).SimName, "Test");
	}


	// sprintf((*par).SimName, "Test_Function");


	//INPUT METHODS PARAMETERS---------------------------------------
	//0 : Initial Guess with null Vector 
	//1 : Initial Guess with SDIRK
	(*par).IG_FLAG = 0;

	//0 : Newton-Raphson LU Solver 
	//1 : Newton-Raphson BiCGStab + SDIRK Preconditioner
	(*par).Solver_FLAG = 0;

	//0 : Linear Equations with time-independent coefficients
	//1 : Equations are linear with time-dependent coefficients
	//2 : Non Linear Equations
	(*par).EquationType_FLAG = 0;  





	//----Fix Parameters for SDIRK Pre-conditioner
	(*par).SDIRK_s= 3;
	Butcher_Tableau(par);

	//DERIVED NUMERICAL PARAMETERS------------------------------
	Get_Index_Arrays(par);
	Get_Differentiation_Matrices(par);

	Get_DiffMatrices_2D(par);	
	Allocate_PhysDiffMatrices_2D(par);
	Alocate_Diff_Operator_Functions(par);
	
	if((*par).EquationType_FLAG == 0){
		Allocate_LinearStatic_Jacobian(par);
	}


	char mkdir[200];
	sprintf(mkdir, "mkdir -p data/%s",(*par).SimName) ;
	int func_out =  system(mkdir);
	if(func_out){ fprintf(stderr, "Error in set_parameters: Unable to create data directory\n"); exit(1);}

	(*par).flag_debug = 0;
	
	return;
}
//-----------------------------------------------
void free_parameters(parameters *par){

	Free_Index_Arrays(par);

	Free_Differentiation_Matrices(par);
	free_DiffMatrices_2D(par);
	free_Diff_Operator_Functions(par);



	free_Butcher_Tableau(par);
	if((*par).EquationType_FLAG == 0){
		Free_LinearStatic_Jacobian(par);
	}

	// if((*par).PertOrder==1)
	// 	{free_QNM_ID(par);}
	// if((*par).PertOrder==2)
	// 	{free_Source(par);}

	return;
}