#include "1+1_Free_Boundaries.h"

void read_header(parameters *par, FILE *fr, double complex *s, int *N){

	double Re_s, Im_s;
	int N_read;//, QNM_L, QNM_M;
	// double Re_Amp_RegMode, Im_Amp_RegMode,  Re_Amp_MirrorMode, Im_Amp_MirrorMode;

	// if((*par).PertOrder ==2){
	// 	fscanf(fr, "QQNM$l\t%d\n", &QNM_L);	
	// 	fscanf(fr, "QQNM$m\t%d\n", &QNM_M);
	// 	if(QNM_L != (*par).ell || QNM_M != (*par).m){
	// 		fprintf(stderr,"Error in function read_header\n");
	// 		fprintf(stderr,"read (L,M) = (%d,%d) doesn't agree with parameter (L,M) = (%d,%d)\n", QNM_L,QNM_M, (*par).ell, (*par).m);
	// 		exit(-1);
	// 	}
	// }
	

	fscanf(fr, "s\t%lf\t%lf\n", &Re_s, &Im_s);
	*s = Re_s + I * Im_s;

	// if((*par).PertOrder ==2){
	// 	fscanf(fr, "AmplitudeReg\t%lf\t%lf\n", &Re_Amp_RegMode, &Im_Amp_RegMode);	
	// 	fscanf(fr, "AmplitudeMirror\t%lf\t%lf\n", &Re_Amp_MirrorMode, &Im_Amp_MirrorMode);
	// 	(*par).Amp_RegularMode = Re_Amp_RegMode + I*Im_Amp_RegMode;
	// 	(*par).Amp_MirrorMode = Re_Amp_MirrorMode + I*Im_Amp_MirrorMode;
	// }
	
	fscanf(fr, "N\t%d\n\n", &N_read);
	*N = N_read;


	if(nDoms>1 ){
		fprintf(stderr,"Error in function read_header\n");
		fprintf(stderr,"Multi-domain not implemented: nDoms = %d\n", nDoms);
		exit(-1);
	}
	
	return;
}
//----------------------------------------------
void load_QNM_ID(parameters *par){
	FILE *fr;

	char fn[500];
	int iDom=0;
	if(nDoms>1 ){
		fprintf(stderr,"Error in function load_QNM_ID\n");
		fprintf(stderr,"Multi-domain not implemented: nDoms = %d\n", nDoms);
		exit(-1);
	}

	// for(iDom=0; iDom<nDoms; iDom++){

		sprintf(fn,"input_data/Cheb_QNM_N_%d_n%d_l%d_spin%d_Prec_%3.3f.dat", (*par).N_input_QNM[iDom], (*par).n_QNM, (*par).ell, (*par).spin, (*par).prec_input_QNM);

		if( (fr = fopen(fn,"r")) == NULL){
			fprintf(stderr, "Error in load_data: file %s not found\n", fn);
			exit(-1);
		}

		int i1, N_read;
		double complex s;
		read_header(par, fr, &s, &N_read);
		(*par).s = s;
		(*par).N_input_QNM[iDom] = N_read;


		(*par).cheb_real_phi_ID[iDom] = dvector(0, (*par).N_input_QNM[iDom]);
		(*par).cheb_imag_phi_ID[iDom] = dvector(0, (*par).N_input_QNM[iDom]);

		fscanf(fr,"#1:i	 #2:Real(Cheb_phi_QNM)	 #3:Imag(Cheb_phi_QNM)\n");

		for(i1=0; i1<=(*par).N_input_QNM[iDom]; i1++){
			int i1_read;
			double c_Re_phi, c_Im_phi;

			fscanf(fr,"%d\t%lf\t%lf\n", &i1_read, &c_Re_phi, &c_Im_phi);
			(*par).cheb_real_phi_ID[iDom][i1] = c_Re_phi; 
			(*par).cheb_imag_phi_ID[iDom][i1] = c_Im_phi;				
		}
		fclose(fr);

	// }
	
	return;
}
//----------------------------------------------
void load_Source(parameters *par){
	FILE *fr;

	char fn[500];
	int iDom=0;
	if(nDoms>1 ){
		fprintf(stderr,"Error in function load_QNM_ID\n");
		fprintf(stderr,"Multi-domain not implemented: nDoms = %d\n", nDoms);
		exit(-1);
	}

	// for(iDom=0; iDom<nDoms; iDom++){

		sprintf(fn,"input_data/ReRegMode_%3.5f_ImRegMode_%3.5f_ReMirrorMode_%3.5f_ImMirrorMode_%3.5f/l_%d/m_%d/2ndOrderSource_N%d.dat", creal((*par).Amp_RegularMode), cimag((*par).Amp_RegularMode), creal((*par).Amp_MirrorMode) ,cimag((*par).Amp_RegularMode), (*par).ell, (*par).m,(*par).N_input_Source[iDom]);
		if( (fr = fopen(fn,"r")) == NULL){
			fprintf(stderr, "Error in load_Source: file %s not found\n", fn);
			exit(-1);
		}

		int i1, N_read;
		double complex s;
		
		read_header(par, fr, &s, &N_read);
		(*par).ss = s;
		(*par).N_input_Source[iDom] = N_read;
		// printf("%d %3.15e %3.15e\n", N_read, creal(s), cimag(s));		
		// pause();

		(*par).cheb_real_Source[iDom] = dvector(0, (*par).N_input_Source[iDom]);
		(*par).cheb_imag_Source[iDom] = dvector(0, (*par).N_input_Source[iDom]);
		fscanf(fr, "#1:i	 #2:Real(Cheb_S2nd)	 #3:Imag(Cheb_S2nd)\n");
		for(i1=0; i1<=(*par).N_input_Source[iDom]; i1++){
			int i1_read;
			double c_Re_source, c_Im_source;
			
			fscanf(fr,"%d\t%lf\t%lf\n", &i1_read, &c_Re_source, &c_Im_source);			
			(*par).cheb_real_Source[iDom][i1] = c_Re_source; 
			(*par).cheb_imag_Source[iDom][i1] = c_Im_source;			
		}
		fclose(fr);
	
		
	// }
	
	//------------
	// for(i=0; i<=(*par).N_input_Source[iDom]; i++){
	// 	printf("%d %3.15e %3.15e\n", i, (*par).cheb_real_Source_pm[iDom][i], (*par).cheb_imag_Source_pm[iDom][i]);
	// }
	// exit(-1);
	//------------

	return;
}
// //----------------------------------------------
// void load_Source(parameters *par){
// 	FILE  *fr, *fr_pp, *fr_pm;

// 	char fn[500], fn_pp[500], fn_pm[500];
// 	int iDom=0, N_read;
	
// 	// //WARNING: HARDCODED FOR nDoms = 1
// 	sprintf(fn,"input_data/Cheb_QNM_N_%d_n%d_l%d_spin%d_Prec_%3.3f.dat", (*par).N_input_QNM[0], (*par).n_QNM, (*par).ell, (*par).spin, (*par).prec_input_QNM);
// 	if( (fr = fopen(fn,"r")) == NULL){
// 			fprintf(stderr, "Error in load_data: file %s not found\n", fn);
// 			exit(-1);
// 	}	
// 	read_header(par, fr);
// 	fclose(fr);

// 	(*par).s_pp = (*par).s + (*par).s;
// 	(*par).s_pm = (*par).s + conj((*par).s);
// 	double *a_pp, *b_pp, *a_pm, *b_pm;
// 	// for(iDom=0; iDom<nDoms; iDom++){

// 		sprintf(fn_pp,"input_data/2ndOrderSource_++.dat");
// 		sprintf(fn_pm,"input_data/2ndOrderSource_+-.dat");

// 		if( (fr_pp = fopen(fn_pp,"r")) == NULL){
// 			fprintf(stderr, "Error in load_Source: file %s not found\n", fn_pp);
// 			exit(-1);
// 		}
// 		if( (fr_pm = fopen(fn_pm,"r")) == NULL){
// 			fprintf(stderr, "Error in load_Source: file %s not found\n", fn_pm);
// 			exit(-1);
// 		}

		
// 		fscanf(fr_pp, "N\t%d\n\n", &N_read);
// 		printf("%d \n", N_read);
// 		fscanf(fr_pm, "N\t%d\n\n", &N_read);
// 		printf("%d \n", N_read);
// 		(*par).N_input_Source[iDom] = N_read;
// 		pause();

// 		a_pp=dvector(0,N_read); b_pp=dvector(0,N_read);
// 		(*par).cheb_real_Source_pp[iDom] = dvector(0, (*par).N_input_QNM[iDom]);
// 		(*par).cheb_imag_Source_pp[iDom] = dvector(0, (*par).N_input_QNM[iDom]);

// 		a_pm=dvector(0,N_read); b_pm=dvector(0,N_read);
// 		(*par).cheb_real_Source_pm[iDom] = dvector(0, (*par).N_input_QNM[iDom]);
// 		(*par).cheb_imag_Source_pm[iDom] = dvector(0, (*par).N_input_QNM[iDom]);

// 		fscanf(fr_pp,"#1:i	 #2:Real(Cheb_S2nd_++)	 #3:Imag(Cheb_S2nd_++)\n");
// 		fscanf(fr_pm,"#1:i	 #2:Real(Cheb_S2nd_+-)	 #3:Imag(Cheb_S2nd_+-)\n");
		
// 		int i1=0;
// 		for(i1=0; i1<=(*par).N_input_Source[iDom]; i1++){
// 			int i1_read;
// 			double c_Re_source_pp, c_Im_source_pp, c_Re_source_pm, c_Im_source_pm;

// 			fscanf(fr_pp,"%d\t%lf\t%lf\n", &i1_read, &c_Re_source_pp, &c_Im_source_pp);
// 			// printf("%d %d %3.15e %3.15e\n", i1, i1_read, c_Re_source_pp, c_Im_source_pp );
// 			a_pp[i1]=c_Re_source_pp;
// 			b_pp[i1]=c_Im_source_pp;
// 			// printf("%d %d %3.15e %3.15e\n", i1, -1, a_pp[i1], b_pp[i1] );
			
// 			fscanf(fr_pm,"%d\t%lf\t%lf\n", &i1_read, &c_Re_source_pm, &c_Im_source_pm);
// 			printf("%d %d %3.15e %3.15e\n", i1, i1_read, c_Re_source_pm, c_Im_source_pm );
// 			a_pm[i1]=c_Re_source_pm;
// 			b_pm[i1]=c_Im_source_pm;
// 			printf("%d %d %3.15e %3.15e\n", i1, -1, a_pm[i1], b_pm[i1] );			

				
// 		}
// 		fclose(fr_pp);
// 		fclose(fr_pm);
// 		pause();

// 	// }

// 		//---------
// 	printf("In load Source\n");
// 	int i; 
// 	for(i=(*par).N_input_Source[0]; i>=0; i--){
// 		(*par).cheb_real_Source_pp[0][i] = a_pp[i];
// 		(*par).cheb_imag_Source_pp[0][i] = b_pp[i];

// 		printf("%d %3.15e %3.15e\n", i, (*par).cheb_real_Source_pp[0][i],(*par).cheb_imag_Source_pp[0][i] );
// 		printf("%d %3.15e %3.15e\n", i, a_pp[i], b_pp[i] );

// 	}
// 	pause();
// 	for(i=(*par).N_input_Source[0]; i>=0; i--){
// 		(*par).cheb_real_Source_pm[0][i] = a_pm[i];
// 		(*par).cheb_imag_Source_pm[0][i] = b_pm[i];

// 		printf("%d %3.15e %3.15e\n", i, (*par).cheb_real_Source_pm[0][i],(*par).cheb_imag_Source_pm[0][i] );
// 		printf("%d %3.15e %3.15e\n", i, a_pm[i], b_pm[i] );

// 	}
// 	exit(-1);
// 	//----------

	
// 	return;
// }
//--------------------------------
void free_QNM_ID(parameters *par){
	int iDom;

	for(iDom=0; iDom<nDoms; iDom++){
		free_dvector( (*par).cheb_real_phi_ID[iDom],       0, (*par).N_input_QNM[iDom]);
		free_dvector( (*par).cheb_imag_phi_ID[iDom], 0, (*par).N_input_QNM[iDom]);
	}

	return;
}
//--------------------------------
void free_Source(parameters *par){
	int iDom;

	for(iDom=0; iDom<nDoms; iDom++){
		free_dvector( (*par).cheb_real_Source[iDom], 0, (*par).N_input_Source[iDom]);
		free_dvector( (*par).cheb_imag_Source[iDom], 0, (*par).N_input_Source[iDom]);
	}

	return;
}