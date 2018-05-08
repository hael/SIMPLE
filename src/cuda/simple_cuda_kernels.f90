module FortCUDA_kernels
implicit none
  interface

     subroutine vecAddF95(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddf95")
       use, intrinsic :: ISO_C_BINDING
       use CUDA, only : dim3, cudaStream_t
       type (c_ptr), value :: a, b, c
       type (dim3) :: dimGrid
       type (dim3) :: dimBlk
       integer(c_int), value :: N
       type (cudaStream_t) :: stream
     end subroutine vecAddF95

  end interface

end module FortCUDA_kernels
