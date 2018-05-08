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

#include "cuda.h"
#include <stdio.h>

__global__ void vecAdd( const float *A, const float *B, float *C, int N )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < N ) {
	C[i] = A[i] + B[i];
    }
}

extern "C"
{
    void vecaddf95(float *A, float *B, float *C, dim3 *dimGrid, dim3 *dimBlk,
	           int N, cudaStream_t *stream)
    {
	vecAdd<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, N);
    }

}
