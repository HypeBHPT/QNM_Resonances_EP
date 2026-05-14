#include "1+1_Free_Boundaries.h"

void pause() {

    printf("\nPress any key to continue\n(or 'q' to quit)\n");

    int pause = getchar(); // Read a character from input

    if (pause == EOF) {
        printf("Error in 'Pause'\n");
        return;
    }

    if (pause == 'q' || pause == 'Q') {
        printf("Exiting...\n");
        exit(0); // Exit the program
    }

    // int c;
    // while ((c = getchar()) != '\n' && c != EOF); // Consume characters until newline or EOF

    printf("\n");
    return;
}
//------------------------------------------------------------------------------
// void pause(){
//   struct utsname unameData;
//   uname(&unameData);

//   char pause;
//   printf("\nPress any key to continue\n(or 'q' to quit)\n");
  
//   // __fpurge(stdin);
//   fpurge(stdin);

//   int ret;
//   ret = scanf("%c", &pause);
//   if (ret == EOF)
//   {
//     printf("...Error in 'Pause'\n");
//     exit(1);
//   }
//   if(pause=='q')
//     exit(1);
//   // __fpurge(stdin);
//   fpurge(stdin);
//   printf("\n");

//   return;
// }
//------------------------------------------------------------------------------
void PrintVector(double *J, int na, int nb){
 int i;
 
 for(i=na; i<=nb; i++){
   printf("vec[%d]=%3.15e\n",i, J[i]);
//    pause();
 }
}
//------------------------------------------------------------------------------
void PrintMatrix(char *fn, double **J, int la, int lb, int ca, int cb){
 int i, j;

 FILE *fp = fopen(fn,"w");
 
 for(i=la; i<=lb; i++){
   for(j=ca; j<=cb; j++){
      double out = J[i][j];

      // if(isnan(out)){
        // printf("%d, %d, %3.2e\n", i, j, J[i][j]);
        // pause();
      // }

     fprintf(fp, "%3.2e ", out );
   }
   fprintf(fp, "\n");
 }
 fclose(fp);

 return;
}
//------------------------------------------------------------------------------
void PrintMatrix_to_File(parameters par, char *fn, double **J, int la, int lb, int ca, int cb){
 int i, j; 
 
 FILE *fp = fopen(fn, "w");

 
 for(i=la; i<=lb; i++){
    // int j1_grid, j2_grid, iF, iDom;
        // Get_Indices_From_Index(par, i, &iDom, &iF, &j1_grid, &j2_grid);
        // printf("I=%d, idom=%d, iF=%d, j1=%d, j2=%d\n", i, iDom, iF, j1_grid, j2_grid);
   for(j=ca; j<=cb; j++){
     double out = (fabs(J[i][j])<1.0e-12 || J[i][j] == -0. ) ? 0.:J[i][j];
//      printf(fp,"%3.5e ", out);
     fprintf(fp,"%3.2lf ", out);
       
     // if(fabs( J[i][j])>5.0e-8 /*isnan(J[i])!=0*/ ){

        
     //    fprintf(fp, "\t\tJ[%d][%d]=%3.15e : Equation (iDom = %d, iF = %d, j1 = %d, j2 = %d)\n", i, j, out, iDom, iF, j1_grid, j2_grid);
       
     //   // pause();
     // }    
     
   }
//    printf("\n");
   fprintf(fp, "\n");
 }
 fclose(fp);
}
//------------------------------------------------------------------------------
void Check_OS(){
  struct utsname unameData;
  uname(&unameData);

  // Check if the operating system is Mac (Darwin)
  if (strcmp(unameData.sysname, "Darwin") == 0) {
      printf("Running on Mac\n\n");
  }
  // Check if the operating system is Linux
  else if (strcmp(unameData.sysname, "Linux") == 0) {
      printf("Running on Linux\n\n");
  }
  // If it's neither Mac nor Linux
  else {
      printf("Unknown OS\n\n");
  }

  return;
}
