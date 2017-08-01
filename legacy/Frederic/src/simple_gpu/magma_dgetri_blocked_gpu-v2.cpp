/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 1st Apr 2015
 *      Monash University
 *      Apr 2015
 *
 *      Routine which calculates the inverse of a matrix A.
 * @precisions normal z -> s d c
 */

/*preprossing for the CUDA and MAGMA environment */
#if defined (CUDA) && defined (MAGMA)

//MAGMA include file from {path}/magma-1.6.1/control
#include "common_magma.h"

//#include "magma_dgetri_blocked_gpu-v2.h"
#include "simple_math_gpu.h"

// includes, project
#include "flops.h"

// === Define what BLAS to use ============================================
#define PRECISION_d
#define BLOCK_SIZE 128

#define getri_min(a,b) ((a) <= (b) ? (a) : (b))
#define getri_max(a,b) ((a) >= (b) ? (a) : (b))

#if (defined(PRECISION_s) || defined(PRECISION_z))
  #define cublasDgemm magmablas_dgemm
  #define cublasDtrsm magmablas_dtrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_dgetri_block_gpu_v2( magma_int_t n,
			   double *dA,
			   magma_int_t ldda, magma_int_t *ipiv,
			   double *work, magma_int_t *lwork,
			   magma_int_t *info) {

  //printf("At entry of magma_dgetri_block_gpu_v2\n");
  //printf("---------------------------------------------------\n");

  double d_one     = MAGMA_D_ONE;
  double d_zero    = MAGMA_D_ZERO;

  //variables used by magma
  magma_int_t lquery;             //query variables
  magma_int_t nb;                 //block size
  magma_int_t s;                  //number of blocks   
  magma_int_t nn;                 //number of blocks   

  double *dAT;           //pointer for A matrix 
  double *work_dgetri;   //pointer for the intermediate variables unblocked
  double *work_dtrti2;   //pointer for the intermediate variables block
  double *work_work;     //pointer for the temporary matrix

  magma_int_t nbmin;              //minimum value for the block
  magma_int_t ldwork;
  magma_int_t iws;
  
  magma_int_t mat_size = n * n * sizeof(double); //setting the size of the matrices

  //counters and variables used for loops and the like
  int i__1;
  int i, j, jj;
  int jb, jp;

  //timers variables

  struct timeval start_mainTime,end_mainTime;
  struct timeval start_dtrti2Time,end_dtrti2Time;
  struct timezone tz;
  double elapsed_timeMain;
  double elapsed_dtrti2Time;

  FILE * timeFile;

  //struct timeval start_magmaDZgemm,end_magmaDZgemm;
  double elapsed_magmaDZgemm;
  double flops;

#ifdef BENCH
  timeFile = fopen("Times_magma_dgetri_block_gpu_v2.log","a");
#endif

  // Function Body

  gettimeofday(&start_mainTime, &tz);

  *info = 0;
  work[1] = n * nb * d_one;
  lquery = *lwork == -1;  
  if( n < 0 )
    {
      *info = -1;
      printf("*info: %i\n",*info);
    }
  else if( ldda < getri_max( 1, n ) )
    {
      *info = -3;
      printf("*info: %i, *lwork: %i, n: %i, ldda: %i\n",*info,*lwork,n,ldda);
    }
  else if( *lwork < getri_max( 1, n ) && ! lquery )
    {
      *info = -6;
      printf("*info: %i, *lwork: %i, n: %i, ldda: %i\n",*info,*lwork,n,ldda);
    }

  if( *info != 0 )
    {
      i__1 = -(*info);
      cublasXerbla ("Magma_DGETRI", i__1);
      return 0;
    }
  else if( lquery )
    {
      return 0;
    }

  // Quick return if possible
  if ( n == 0 ) return MAGMA_SUCCESS;

  // Function Body
  nb = magma_get_getri_nb_gpu(n);   //Definition of blocking sizes 
  s = n / nb;                    //calculating the number of blocks

  nbmin = 2;
  if ( nb <= nbmin || nb >= n )
    {
      /*
	the CPU code for small size n <= 128
      */
      printf("*        Use unblocked code on CPU.\n");
      printf("*        Using the standard lapack CPU code\n");

      work_dgetri = (double*)malloc((n) * (n) * sizeof(double));
      cublasGetMatrix(ldda, n, sizeof(double), dA, ldda, work_dgetri, n);
      printf("---------------------------------------------------------------------\n");
      lapackf77_dgetri( &n, work_dgetri, &ldda, ipiv, work, lwork, info );
      cublasSetMatrix(ldda, n, sizeof(double), work_dgetri, n, dA, ldda);
      free(work_dgetri);
    }
  else
    {
      /*
	the GPU code for size n > 128
      */
      /*
      printf("*\n");
      printf("*        Use block matrix to solve AX = I on GPU where A is obbtained from\n");
      printf("*        magma_dgetrf_gpu.cpp and A = P * L * U.\n");
      printf("*        Using the standard CuBlas code for GPU code\n");
      printf("*        First get inv(U) then\n");
      printf("*        Solve the linear system inv(A)*L = inv(U)\n");
      printf("*\n");
      */
      dAT = dA;

      if( n == 0 ) return   MAGMA_SUCCESS; //check if the matrix size n>0 
      //
      //compute the inverse of the upper part, inv(U)
      //
      work_dtrti2 = (double*)malloc((nb) * (nb) * sizeof(double)); //allc temp mem. var. nb blocks

      for ( j = 1 ; j <= n ; j += nb ) //loop over the blocks 
	{
	  jb = getri_min( nb , n - j + 1 );

	  // compute rows 0:j of current block column

	  cublasDtrmm( SimpleLeft, SimpleUpper, SimpleNoTrans, SimpleNonUnit, 
		       j-1, jb,
		       d_one,
		       dAT, ldda,
		       &dAT[(j-1) * ldda], ldda );
	  cublasDtrsm( SimpleRight, SimpleUpper, SimpleNoTrans, SimpleNonUnit,
		       j-1, jb,
		       -d_one,
		       &dAT[(j-1) + (j-1) * ldda], ldda,
		       &dAT[(j-1) * ldda], ldda );
	  
	  //making sure that the GPU is synchronized with CPU before proceeding
	  cuCtxSynchronize();

	  // Compute inverse of current diagonal block using
	  // the standard lapack77 for small blocks on the CPU

	  gettimeofday(&start_dtrti2Time, &tz);

	  cublasGetMatrix(jb, jb, sizeof(double), &dAT[(j-1) + (j-1) * ldda], ldda, &work_dtrti2[0], jb);
	  lapackf77_dtrti2( "Upper", "Non-unit", &jb, &work_dtrti2[0], &jb, info );
	  cublasSetMatrix(jb, jb, sizeof(double), &work_dtrti2[0], jb, &dAT[(j-1)+(j-1)*ldda], ldda);

	  gettimeofday(&end_dtrti2Time, &tz);

	} //end of loop over the blocks 
      
      nbmin = 2;
      ldwork = n;

      //
      // Solve the equation inv(A)*L = inv(U) for inv(A).
      //

      nn = ( ( n - 1 ) / nb )*nb + 1;
      cudaMalloc((void **) &work_work, (nb)*(ldwork)*sizeof(double));

      int count = 0;
      for ( j = nn ; j >= 1 ; j -= nb ) // loop over the blocks 
	{
	  jb = getri_min( nb , n - j + 1);

	  // Copy current block column of L to work_work and replace A(i,jj) with zeros.

	  int counter = 0;
	  for ( jj = j ; jj <= j + jb - 1 ; ++jj )
	    {
	      counter = counter + 1;
	      
	      cudaMemcpy(&work_work[ jj + (jj - j)*ldwork ],

			 &dAT[(  jj + (jj-1) * ldda )],

			 (n-(jj+1)+1)*sizeof(double), cudaMemcpyDeviceToDevice);

	      cudaMemset(&dAT[ jj + (jj-1) * ldda ], 0, (n-(jj+1)+1)*sizeof(double));

	    }

	  count = count + counter;
	  //printf("When (j): %i , the after loop (counter): %i, tot (count): %i\n",j,counter,count);

	  // Compute current block column of inv(A).

	  if ( j + jb <= n )
	    {
	      cublasDgemm( SimpleNoTrans, SimpleNoTrans,
			   n, jb, n-j-jb+1,
			   -d_one,
			   &dAT[ ( (j-1) + jb) * ldda ], ldda,
			   &work_work[ (j-1) + jb ], ldwork,
			   d_one,
			   &dAT[(j-1) * ldda ], ldda );

	    }
          cublasDtrsm( SimpleRight, SimpleLower, SimpleNoTrans, SimpleUnit,
		       n, jb,
		       d_one,
		       &work_work[(j-1)], ldwork,
		       &dAT[(j-1) * ldda ], ldda );

	  //making sure that the GPU is synchronized with CPU before proceeding
	  cuCtxSynchronize();

	} //end of loop over block (j=nn, 1, -nb)
      

      //
      // Apply column interchanges.
      //

      for (j = n - 1; j >= 1; --j)
	{
	  jp = ipiv[j-1];
	  if (jp != j)
	    {
	      cublasDswap(n, &dAT[ (j-1) * ldda ], 1, &dAT[ (jp-1) * ldda ], 1 );
	    }
	}

      //
      // Freeing memory at the end.
      //

      //on the CPU 
      free(work_dtrti2);
      
    } //end of the else if for blocked code

#ifdef BENCH

  gettimeofday(&end_mainTime, &tz);

  elapsed_timeMain = (end_mainTime.tv_sec  - start_mainTime.tv_sec) + 
                     (end_mainTime.tv_usec - start_mainTime.tv_usec)/1000000.0 ;

  elapsed_dtrti2Time =  (end_dtrti2Time.tv_sec  - start_dtrti2Time.tv_sec) + 
                     (end_dtrti2Time.tv_usec - start_dtrti2Time.tv_usec)/1000000.0 ;

  elapsed_magmaDZgemm = (end_mainTime.tv_sec  - start_mainTime.tv_sec)*1000.0     +
                        (end_mainTime.tv_usec - start_mainTime.tv_usec)*0.001     ;

  //flops = FLOPS((double)(*m), (double)(*n), (double)(*k)) / 1000000;
  flops = ( 6. * FMULS_GETRI((double)(n)) + 
	    2. * FADDS_GETRI((double)(n)) ) / 1000000;


  if ( n <= 128 )
    {
      printf("     [nxn]      Elements    Size (MB)      time (sec)   GFLop/s\n");
      printf("-----------------------------------------------------------------------------------------------\n");
      printf("   %i x %i    %i        %.2f   %15.8f          %6.2f\n",
	     n,n,mat_size, mat_size/1.e6,elapsed_timeMain, ( flops / elapsed_magmaDZgemm ));

      fprintf(timeFile,"%i, %.2f, %15.8f, %15.8f, %8.4f, %6.2f\n",
	      n,mat_size/1.e6,
	      elapsed_timeMain, 
	      ( flops / elapsed_magmaDZgemm ));
    }
  else
    {
      printf("     [nxn]        Elements        Size (MB)      time (sec)     time(Dtrti2) (sec)     (%%) of Dtrti2        GFLop/s\n");
      printf("-------------------------------------------------------------------------------------------------------------------\n");
      printf("   %i x %i    %i        %.2f   %15.8f   %15.8f           %8.4f             %6.2f\n",
	     n,n,mat_size, mat_size/1.e6,elapsed_timeMain,elapsed_dtrti2Time,
	     (elapsed_dtrti2Time/elapsed_timeMain)*100.0, ( flops / elapsed_magmaDZgemm ) );

      fprintf(timeFile,"%i, %.2f, %15.8f, %15.8f, %8.4f, %6.2f\n",
	      n,mat_size/1.e6,
	      elapsed_timeMain,
	      elapsed_dtrti2Time, 
	      (elapsed_dtrti2Time/elapsed_timeMain)*100.0, 
	      ( flops / elapsed_magmaDZgemm ) );
    }

  fclose(timeFile);

#endif

  //printf("Bye magma_dgetri_block_gpu_v2\n");    

  return MAGMA_SUCCESS;

  /* End of magma_dgetri_gpu */
}
#undef inAT

#endif /* (CUDA) && (MAGMA) */
