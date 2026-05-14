#include "1+1_Free_Boundaries.h"

void load_QNMSet(parameters par, char *fn, double complex *qnm_s, double complex *qnm_phi, double complex *qnm_psi){

	// char fn = "InputData/QNMSet0_N_200_Freq_eps_1._ksig_10_Prec_159.dat";
	FILE *fr=fopen(fn, "r"); 

	long double Real_read, Imag_read;

	fscanf(fr,"%Lf\t%Lf", &Real_read, &Imag_read);

	*qnm_s = Real_read + I*Imag_read;
	

	int count = 0;
	while (fscanf(fr,"%Lf\t%Lf", &Real_read, &Imag_read)!=EOF){
		qnm_phi[count] = Real_read + I*Imag_read;
		qnm_psi[count] = (*qnm_s)*qnm_phi[count];
		
		
		count++; 
		if(count>par.NA[0]+1){
			printf("Error in load_QNMSet: Data is bigger than resolution. Interpolantion not implemented\n");
			exit(1);
		}
		 
	}


	if(count<par.NA[0]+1){
			printf("Error in load_QNMSet: Data is smaller than resolution. Interpolantion not implemented\n");
			exit(1);
		} 

	fclose(fr);
	return;
}
//------------------------------
double ID_Functions(double sigma){
	double f;

	f = sigma*(1-sigma);

	return f;
}
//------------------------------
double dID_Functions(double sigma){
	double f;

	f = 0.;

	return f;
}