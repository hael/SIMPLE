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

//#include "magma_zgetri_blocked_gpu-v2.h"
#include "simple_math_gpu.h"

// includes, project
#include "flops.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#define BLOCK_SIZE 128

#define getri_min(a,b) ((a) <= (b) ? (a) : (b))
#define getri_max(a,b) ((a) >= (b) ? (a) : (b))

#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_zgetri_block_gpu_v2( magma_int_t n,
			   cuDoubleComplex *dA,
			   magma_int_t ldda, magma_int_t *ipiv,
			   cuDoubleComplex *work, magma_int_t *lwork,
			   magma_int_t *info) {

  //printf("At entry of magma_zgetri_block_gpu_v2\n");
  //printf("---------------------------------------------------\n");

  cuDoubleComplex c_one     = MAGMA_Z_ONE;
  cuDoubleComplex c_zero    = MAGMA_Z_ZERO;

  //variables used by magma
  magma_int_t lquery;             //query variables
  magma_int_t nb;                 //block size
  magma_int_t s;                  //number of blocks   
  magma_int_t nn;                 //number of blocks   

  cuDoubleComplex *dAT;           //pointer for A matrix 
  cuDoubleComplex *work_zgetri;   //pointer for the intermediate variables unblocked
  cuDoubleComplex *work_ztrti2;   //pointer for the intermediate variables block
  cuDoubleComplex *work_work;     //pointer for the temporary matrix

  magma_int_t nbmin;              //minimum value for the block
  magma_int_t ldwork;
  magma_int_t iws;
  
  magma_int_t mat_size = n * n * sizeof(cuDoubleComplex); //setting the size of the matrices

  //counters and variables used for loops and the like
  int i__1;
  int i, j, jj;
  int jb, jp;

  //timers variables

  struct timeval start_mainTime,end_mainTime;
  struct timeval start_ztrti2Time,end_ztrti2Time;
  struct timezone tz;
  double elapsed_timeMain;
  double elapsed_ztrti2Time;

  FILE * timeFile;

  //struct timeval start_magmaDZgemm,end_magmaDZgemm;
  double elapsed_magmaDZgemm;
  double flops;

#ifdef BENCH
  timeFile = fopen("Times_magma_zgetri_block_gpu_v2.log","a");
#endif

  // Function Body

  gettimeofday(&start_mainTime, &tz);

  *info = 0;
  work[1] = n * nb * c_one;
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
      cublasXerbla ("Magma_ZGETRI", i__1);
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

      work_zgetri = (cuDoubleComplex*)malloc((n) * (n) * sizeof(cuDoubleComplex));
      cublasGetMatrix(ldda, n, sizeof(cuDoubleComplex), dA, ldda, work_zgetri, n);
      printf("---------------------------------------------------------------------\n");
      lapackf77_zgetri( &n, work_zgetri, &ldda, ipiv, work, lwork, info );
      cublasSetMatrix(ldda, n, sizeof(cuDoubleComplex), work_zgetri, n, dA, ldda);
      free(work_zgetri);
    }
  else
    {
      /*
	the GPU code for size n > 128
      */
      /*
      printf("*\n");
      printf("*        Use block matrix to solve AX = I on GPU where A is obbtained from\n");
      printf("*        magma_zgetrf_gpu.cpp and A = P * L * U.\n");
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
      work_ztrti2 = (cuDoubleComplex*)malloc((nb) * (nb) * sizeof(cuDoubleComplex)); //allc temp mem. var. nb blocks

      for ( j = 1 ; j <= n ; j += nb ) //loop over the blocks 
	{
	  jb = getri_min( nb , n - j + 1 );

	  // compute rows 0:j of current block column

	  cublasZtrmm( SimpleLeft, SimpleUpper, SimpleNoTrans, SimpleNonUnit, 
		       j-1, jb,
		       c_one,
		       dAT, ldda,
		       &dAT[(j-1) * ldda], ldda );
	  cublasZtrsm( SimpleRight, SimpleUpper, SimpleNoTrans, SimpleNonUnit,
		       j-1, jb,
		       -c_one,
		       &dAT[(j-1) + (j-1) * ldda], ldda,
		       &dAT[(j-1) * ldda], ldda );

	  //making sure that the GPU is synchronized with CPU before proceeding
	  cuCtxSynchronize();

	  // Compute inverse of current diagonal block using
	  // the standard lapack77 for small blocks on the CPU

	  gettimeofday(&start_ztrti2Time, &tz);

	  cublasGetMatrix(jb, jb, sizeof(cuDoubleComplex), &dAT[(j-1) + (j-1) * ldda], ldda, &work_ztrti2[0], jb);
	  lapackf77_ztrti2( "Upper", "Non-unit", &jb, &work_ztrti2[0], &jb, info );
	  cublasSetMatrix(jb, jb, sizeof(cuDoubleComplex), &work_ztrti2[0], jb, &dAT[(j-1)+(j-1)*ldda], ldda);

	  gettimeofday(&end_ztrti2Time, &tz);

	} //end of loop over the blocks 
      
      nbmin = 2;
      ldwork = n;

      //
      // Solve the equation inv(A)*L = inv(U) for inv(A).
      //

      nn = ( ( n - 1 ) / nb )*nb + 1;
      cudaMalloc((void **) &work_work, (nb)*(ldwork)*sizeof(cuDoubleComplex));

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

			 (n-(jj+1)+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);

	      cudaMemset(&dAT[ jj + (jj-1) * ldda ], 0, (n-(jj+1)+1)*sizeof(cuDoubleComplex));

	    }

	  count = count + counter;
	  //printf("When (j): %i , the after loop (counter): %i, tot (count): %i\n",j,counter,count);

	  // Compute current block column of inv(A).

	  if ( j + jb <= n )
	    {
	      cublasZgemm( SimpleNoTrans, SimpleNoTrans,
			   n, jb, n-j-jb+1,
			   -c_one,
			   &dAT[ ( (j-1) + jb) * ldda ], ldda,
			   &work_work[ (j-1) + jb ], ldwork,
			   c_one,
			   &dAT[(j-1) * ldda ], ldda );

	    }
	  cublasZtrsm( SimpleRight, SimpleLower, SimpleNoTrans, SimpleUnit,
		       n, jb,
		       c_one,
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
	      cublasZswap(n, &dAT[ (j-1) * ldda ], 1, &dAT[ (jp-1) * ldda ], 1 );
	    }
	}

      //
      // Freeing memory at the end.
      //

      //on the CPU 
      free(work_ztrti2);
      
    } //end of the else if for blocked code

#ifdef BENCH

  gettimeofday(&end_mainTime, &tz);

  elapsed_timeMain = (end_mainTime.tv_sec  - start_mainTime.tv_sec) + 
                     (end_mainTime.tv_usec - start_mainTime.tv_usec)/1000000.0 ;

  elapsed_ztrti2Time =  (end_ztrti2Time.tv_sec  - start_ztrti2Time.tv_sec) + 
                     (end_ztrti2Time.tv_usec - start_ztrti2Time.tv_usec)/1000000.0 ;

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
      printf("     [nxn]        Elements        Size (MB)      time (sec)     time(Ztrti2) (sec)     (%%) of Ztrti2        GFLop/s\n");
      printf("-------------------------------------------------------------------------------------------------------------------\n");
      printf("   %i x %i    %i        %.2f   %15.8f   %15.8f           %8.4f             %6.2f\n",
	     n,n,mat_size, mat_size/1.e6,elapsed_timeMain,elapsed_ztrti2Time,
	     (elapsed_ztrti2Time/elapsed_timeMain)*100.0, ( flops / elapsed_magmaDZgemm ) );

      fprintf(timeFile,"%i, %.2f, %15.8f, %15.8f, %8.4f, %6.2f\n",
	      n,mat_size/1.e6,
	      elapsed_timeMain,
	      elapsed_ztrti2Time, 
	      (elapsed_ztrti2Time/elapsed_timeMain)*100.0, 
	      ( flops / elapsed_magmaDZgemm ) );
    }

  fclose(timeFile);

#endif

  //printf("Bye magma_zgetri_block_gpu_v2\n");    

  return MAGMA_SUCCESS;

  /* End of magma_zgetri_gpu */
}
#undef inAT

#endif /* (CUDA) && (MAGMA) */
