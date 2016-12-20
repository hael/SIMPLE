
/* summing the 3D matrix treating it a 1D vec */
extern "C" __global__ void
psum_1D( float *A, float *partial, int N) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int sharedIndex = threadIdx.x;
  //allocating the shared memory array
  __shared__ float shared_A[256];

  //summing over the 
  float temp = 0.0;
  while ( i < N ){
    temp += A[i];
    i += blockDim.x * gridDim.x;
  } 
  shared_A[sharedIndex] = temp;
  //syncronizing the threads
  __syncthreads();

  //practising the 1/2 reducion sum at each step
  int j = blockDim.x / 2.0;
  while ( j != 0 ) {
    if ( sharedIndex < j ) shared_A[sharedIndex] += shared_A[sharedIndex + j];
    __syncthreads();
    j /= 2;
  }
  if ( sharedIndex == 0 ) partial[blockIdx.x] = shared_A[0];
}

/* summing the 3D matrix in the nrot and nk dimension */
extern "C" __global__ void
psum_23D(float *partial, float *A, int npart, int nrot, int nk) {

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;
  int    ik = blockIdx.z * blockDim.z + threadIdx.z;

  int sharedIndex_x = threadIdx.x;
  int sharedIndex_y = threadIdx.y;
  int sharedIndex_z = threadIdx.z;
  //allocating the shared memory array
  __shared__ float shared_A[16][16][4];
  __shared__ float temp_A[16];

  temp_A[sharedIndex_x] = 0.0;
  while ( ipart < npart ){
    while ( irot < nrot ){
      while ( ik < nk ){
        temp_A[sharedIndex_x] += A[(irot+nrot*ik)*npart+ipart];
        ik += blockDim.z * gridDim.z;
      }
      irot += blockDim.y * gridDim.y;
    }
    ipart += blockDim.x * gridDim.x;
  }
  shared_A[sharedIndex_x][sharedIndex_y][sharedIndex_z] = temp_A[sharedIndex_x];
  //syncronizing the threads
  __syncthreads();

  //practising the 1/2 reducion sum at each step
  int j = blockDim.y / 2.0;
  while ( j != 0 ) {
    if ( sharedIndex_y < j ) {
      int k = blockDim.z / 2.0;
      while ( k != 0 ) {
        if ( sharedIndex_z < k ) shared_A[sharedIndex_x][sharedIndex_y][sharedIndex_z] += 
                                 shared_A[sharedIndex_x][sharedIndex_y + j][sharedIndex_z + k];
        __syncthreads();
        k /= 2;
      }
    }
    __syncthreads();
    j /= 2;
  }

  if ( sharedIndex_y == 0 ) {
    if ( sharedIndex_z == 0 ) {
      if ( ipart < npart ) partial[gridDim.x*blockIdx.x + ipart] = shared_A[sharedIndex_x][0][0];
    }
  }

}

/* summing the 3D matrix in the nrot and nk dimension in a congruent way*/
extern "C" __global__ void
psum_2D(float *partial, float *A, int npart, int nrot, int nk) {

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot_ik = blockIdx.y * blockDim.y + threadIdx.y;
  
  int sharedIndex_x = threadIdx.x;
  int sharedIndex_y = threadIdx.y;
  //allocating the shared memory array
  __shared__ float shared_A[16][64];
  float temp_A;

  temp_A = 0.0;
  while ( ipart < npart ){
    while ( irot_ik < nrot*nk ){
      temp_A += A[irot_ik * npart + ipart];
      irot_ik += blockDim.y * gridDim.y;
    }
    ipart += blockDim.x * gridDim.x;
  }
  shared_A[sharedIndex_x][sharedIndex_y] = temp_A;
  //syncronizing the threads
  __syncthreads();

  //practising the 1/2 reducion sum at each step
  int j = blockDim.y / 2.0;
  while ( j != 0 ) {
    if ( sharedIndex_y < j ) shared_A[sharedIndex_x][sharedIndex_y] += shared_A[sharedIndex_x][sharedIndex_y + j];
    __syncthreads();
    j /= 2;
  }
  
  if ( sharedIndex_y == 0 ) {
    if ( ipart < npart ) 
      partial[gridDim.x*blockIdx.x + ipart] = shared_A[sharedIndex_x][0];
  }

}
/* summing the 3D matrix in the nrot and nk dimension in a congruent way*/
/* used for debugging purpose*/
extern "C" __global__ void
psum_2SD(float *partial, float *A, int npart, int nrot, int nk) {

  int   ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int irot_ik = blockIdx.y * blockDim.y + threadIdx.y;

  int sharedIndex_x = threadIdx.x;
  int sharedIndex_y = threadIdx.y;
  //allocating the shared memory array
  __shared__ float shared_A[16][64];
  float temp_A;

  temp_A = 0.0;

  while ( ipart < npart ){
    while ( irot_ik < nrot*nk ){
      temp_A += A[irot_ik * npart + ipart];
      partial[ipart] = temp_A;//A[irot_ik * npart + ipart];
      irot_ik += blockDim.y * gridDim.y;
    }
    ipart += blockDim.x * gridDim.x;
  }
  shared_A[sharedIndex_x][sharedIndex_y] = 1.0;
  //syncronizing the threads
  __syncthreads();
  partial[ipart] = shared_A[sharedIndex_x][0];

}
/* summing the 3D matrix in the nrot and nk dimension in a congruent way*/
/* used for debugging purpose*/
extern "C" __global__ void
mysum_23D(float *sumA, float *A, int n, int npart){
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int sharedIdx = threadIdx.x;
  __shared__ float shared_A[4];
  float temp = 0.0;
  while ( i < n) {
    float *endptr = A;// + i*npart; 
    temp += *endptr;
    i += blockDim.x * gridDim.x;
  }
  shared_A[sharedIdx] = temp;
  __syncthreads();
  int j = blockDim.x / 2;
  while ( j != 0 ) {
    if ( sharedIdx < j ) shared_A[sharedIdx] += shared_A[sharedIdx + j];
    __syncthreads();
    j /= 2;
  }
  if ( sharedIdx == 0 ) sumA[blockIdx.x] = shared_A[0];
}


/* doing the sum(a,b)sq product Re(A*conjg(A)) */
extern "C" __global__ void
pabsq_3D_mat( float *C,
             const cuFloatComplex *A,
             int npart, int nrot, int nk,
             float alpha, float beta)
{

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;
  int    ik = blockIdx.z * blockDim.z + threadIdx.z;

  if ( ipart < npart) {
    if (irot < nrot ){
      if (ik < nk ){
        C[(irot+nrot*ik)*npart+ipart ] = cuReCCstarf( A[(irot+nrot*ik)*npart+ipart]);
      }
    }
  }
  __syncthreads();
}


