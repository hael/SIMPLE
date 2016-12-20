/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of May 2015
 *      Monash University
 *      May 2015
 *
 *      Routine which screens the vectorised index of 2 and 3D arrays
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/
#include "simple.h"

int main() {

  int n = 3,m=3,l=3; 
  int nBytes_2D, nBytes_3D;
  int lda,ldb;
  int i,j,k;
  float *h_a_2D; //host array
  float *h_b_2D; //host array
  float *h_a_3D; //host array
  float *h_b_3D; //host array

  lda =n;
  ldb =n;

  nBytes_2D = n*n*sizeof(float);
  h_a_2D = (float *)malloc(nBytes_2D);
  h_b_2D = (float *)malloc(nBytes_2D);

  //printf(ANSI_COLOR_BRIGHT_RED"2D representation [i][j]\n"ANSI_COLOR_RESET);
  printf(ANSI_COLOR_BRIGHT_YELLOW
         "[i][j] (h_a_2D[i][j]), (h_b_2D[i][j]), at Line %i in file %s %s\n"
         ANSI_COLOR_RESET,__LINE__,__FILE__,__FUNCTION__);
  for (i=0; i<n ; i++) {
    for (j=0; j<m ; j++) {
      h_a_2D[(i)*n + (j)] = (i)*n + (j) + 1; // +1 to match Fortran
      h_b_2D[(j)*n + (i)] = (j)*n + (i) + 1; // +1 to match Fortran  
      printf(ANSI_COLOR_BRIGHT_GREEN
             "[%i][%i] (%f) , (%f) , at Line %i in file %s %s\n" ANSI_COLOR_RESET,
             i,j,h_a_2D[(i)*n + (j)], h_b_2D[(j)*n + (i)],
             __LINE__,__FILE__,__FUNCTION__);
    }
  }

  free(h_a_2D);
  free(h_b_2D);

  nBytes_3D = n*n*n*sizeof(float);
  h_a_3D = (float *)malloc(nBytes_3D);
  h_b_3D = (float *)malloc(nBytes_3D);
  printf(ANSI_COLOR_BRIGHT_RED);
  printf("3D representation [i][j][k]\n"ANSI_COLOR_RESET);

  for (i=0; i<n ; i++) {
    for (j=0; j<m ; j++) {
      for (k=0; k<l ; k++) {
        h_a_3D[(i+m*k)*n + (j)] = (i+m*k)*n + (j) + 1; // +1 to match Fortran
        h_b_3D[(j+m*k)*n + (i)] = (j+m*k)*n + (i) + 1; // +1 to match Fortran  
        printf(ANSI_COLOR_BRIGHT_BLUE
               "[%i][%i][%i] (%f) , (%f) , at Line %i in file %s %s\n" ANSI_COLOR_RESET,
               i,j,k,h_a_3D[(i+m*k)*n + (j)], h_b_3D[(j+m*k)*n + (i)],
               __LINE__,__FILE__,__FUNCTION__);
      }
    } 
  }

  free(h_a_3D);
  free(h_b_3D);

} /*end main*/
