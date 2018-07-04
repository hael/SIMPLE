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
 * 2018, modified by Michael Eager (michael.eager@monash.edu)
 */

#include "cuda.h"
#include <stdio.h>

__global__ void vec_Add_Int( const int *A, const int *B, int *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B[i];
  }
}
__global__ void vec_Add_Float( const float *A, const float *B, float *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B[i];
  }
}
__global__ void vec_Add_ConstInt( const int *A, const int B, int *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B;
  }
}
__global__ void vec_Add_ConstFloat( const float *A, const float B, float *C, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ) {
    C[i] = A[i] + B;
  }
}


#define VECADDFLOAT vecaddfloat_
#define VECADDINT vecaddint_
#define VECADDCONSTFLOAT vecaddconstfloat_
#define VECADDCONSTINT vecaddconstint_

extern "C"
{
  void vecaddfloat(float *A, float *B, float *C, dim3 *dimGrid, dim3 *dimBlk,
                   int N, cudaStream_t *stream)
  {
    vec_Add_Float<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, N);
  }
  void vecaddint(int *A, int *B, int *C, dim3 *dimGrid, dim3 *dimBlk,
                 int N, cudaStream_t *stream)
  {
    vec_Add_Int<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, N);
  }
  void vecaddconstint(int *A, int *B, int *C, dim3 *dimGrid, dim3 *dimBlk,
                 int N, cudaStream_t *stream)
  {
    vec_Add_ConstInt<<<*dimGrid, *dimBlk, 0, *stream>>>(A, *B, C, N);
  }
  void vecaddconstfloat(float *A, float *B, float *C, dim3 *dimGrid, dim3 *dimBlk,
                   int N, cudaStream_t *stream)
  {
    vec_Add_ConstFloat<<<*dimGrid, *dimBlk, 0, *stream>>>(A, *B, C, N);
  }
}
