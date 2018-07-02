/*
!................................................................................
! Copyright (C) 2009 david.car, david.car7@gmail.com

! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.

! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
................................................................................
*/

/*
 * Modified by Michael Eager (michael.eager@monash.edu) 2018
 */

#include "cuda.h"
#include <stdio.h>

__global__ void vecAddInt( const int *A, const int *B, int *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B[i];
  }
}
__global__ void vecAddFloat( const float *A, const float *B, float *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B[i];
  }
}
__global__ void vecAddConstInt( const int *A, const int B, int *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B;
  }
}
__global__ void vecAddConstFloat( const float *A, const float B, float *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B;
  }
}


extern "C"
{
  void vecadd_float_(float *A, float *B, float *C, dim3 *dimGrid, dim3 *dimBlk,
                   int N, cudaStream_t *stream)
  {
    vecAddFloat<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, N);
  }
  void vecadd_int_(int *A, int *B, int *C, dim3 *dimGrid, dim3 *dimBlk,
                 int N, cudaStream_t *stream)
  {
    vecAddInt<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, N);
  }
  void vecaddconst_int_(int *A, int *B, int *C, dim3 *dimGrid, dim3 *dimBlk,
                 int N, cudaStream_t *stream)
  {
    vecAddConstInt<<<*dimGrid, *dimBlk, 0, *stream>>>(A, *B, C, N);
  }
  void vecaddconst_float_(float *A, float *B, float *C, dim3 *dimGrid, dim3 *dimBlk,
                   int N, cudaStream_t *stream)
  {
    vecAddConstFloat<<<*dimGrid, *dimBlk, 0, *stream>>>(A, *B, C, N);
  }
}
