module simple_cuda_kernels
    include 'simple_lib.f08'
    use , intrinsic :: ISO_C_BINDING
    use CUDA
implicit none
interface
    !! c_kernels.cu
     subroutine vecAddF(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecadd_float")
       use, intrinsic :: ISO_C_BINDING
       use CUDA, only : dim3, cudaStream_t
       type (c_ptr), value :: a, b, c
       type (dim3) :: dimGrid
       type (dim3) :: dimBlk
       integer(c_int), value :: N
       type (cudaStream_t) :: stream
   end subroutine vecAddF
  subroutine vecAddI(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecadd_int")
       use, intrinsic :: ISO_C_BINDING
       use CUDA, only : dim3, cudaStream_t
       type (c_ptr), value :: a, b, c
       type (dim3) :: dimGrid
       type (dim3) :: dimBlk
       integer(c_int), value :: N
       type (cudaStream_t) :: stream
   end subroutine vecAddI
     subroutine vecAddConstF(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddconst_float")
       use, intrinsic :: ISO_C_BINDING
       use CUDA, only : dim3, cudaStream_t
       type (c_ptr), value :: a, b, c ! b is a const
       type (dim3) :: dimGrid
       type (dim3) :: dimBlk
       integer(c_int), value :: N
       type (cudaStream_t) :: stream
   end subroutine vecAddConstF
  subroutine vecAddConstI(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddconst_int")
       use, intrinsic :: ISO_C_BINDING
       use CUDA, only : dim3, cudaStream_t
       type (c_ptr), value :: a, b, c ! b is a const
       type (dim3) :: dimGrid
       type (dim3) :: dimBlk
       integer(c_int), value :: N
       type (cudaStream_t) :: stream
   end subroutine vecAddConstI


  end interface

  !! Matmul.cu
  interface
      subroutine multiply_by_element( dimGrid, dimBlk, a, b, c, N, stream) bind(C, name="multiply_by_element")
          use, intrinsic :: ISO_C_BINDING
          use CUDA, only : dim3, cudaStream_t
          type (c_ptr), value :: a,b, c
          type (dim3) :: dimGrid
          type (dim3) :: dimBlk
          integer(c_int), value :: N
          type (cudaStream_t) :: stream
      end subroutine multiply_by_element
  end interface

  interface
      subroutine multiply_by_block(dimGrid, threads, a, b, c, N, stream) bind(C, name="multiply_by_block")
          use, intrinsic :: ISO_C_BINDING
          use CUDA, only : dim3, cudaStream_t
          type (c_ptr), value :: a, b,c
          type (dim3) :: dimGrid
          type (dim3) :: threads
          integer(c_int), value :: N
          type (cudaStream_t) :: stream
      end subroutine multiply_by_block
  end interface
   interface
      subroutine  transpose_by_block(dimGrid, dimBlk, a,  c, N, stream) bind(C, name="transpose_by_block")
          use, intrinsic :: ISO_C_BINDING
          use CUDA, only : dim3, cudaStream_t
          type (c_ptr), value :: a, c
          type (dim3) :: dimGrid
          type (dim3) :: dimBlk
          integer(c_int), value :: N
          type (cudaStream_t) :: stream
      end subroutine transpose_by_block


  end interface



  !! filter_kernels.cu
   interface
      subroutine filter_gaussKernel(U,V,W,Z, s, dimGrid, dimBlk, N, stream) bind(C, name="gaussKernel")
          use, intrinsic :: ISO_C_BINDING
          use CUDA, only : dim3, cudaStream_t
          type (c_ptr), value ::U,V,W !grid
          type(c_ptr) :: Z ! return value
          real(c_float) :: s ! sigma
          type (dim3) :: dimGrid
          type (dim3) :: dimBlk
          integer(c_int), value :: N
          type (cudaStream_t) :: stream
      end subroutine
  end interface

contains


    subroutine test_FortCUDA_kernels (arg)
        implicit none
        real, intent(in) :: arg
#include "simple_cuda_handle.inc"
        type :: HostPtr
            real, pointer, dimension(:) :: host_ptr => null()
        end type HostPtr

        real :: diff
        integer :: i, j
        integer, parameter :: num_streams = 5
        integer(c_int), parameter :: NTHREADS = 2000000

        type (cudaEvent_t) :: e_start, e_stop
        type (dim3) :: nblk, nthrd
        type (c_ptr) :: err_ptr
        type (cudaStream_t) :: stream(num_streams)
        type (cudaStream_t) :: all_streams
        type (c_ptr) :: d_A(num_streams), d_B(num_streams), d_C(num_streams)
        type (c_ptr) :: h_c_A(num_streams), h_c_B(num_streams), h_c_C(num_streams)
        type (HostPtr) :: h_A(num_streams), h_B(num_streams), h_C(num_streams)
        integer(kind(cudaSuccess)) :: err

        type (cudaDeviceProp) :: prop

        integer (c_int) :: zero = 0
        real (c_float) :: c_time
        integer(c_int) :: sz, flags, device_count
        integer(c_int) :: NTHREADS_PER_BLOCK = 256
        integer(c_int) :: NUM_BLOCKS

        real(KIND=KIND(1.d0)), parameter :: REAL_DP_ELEMENT = 0.d0
        real, parameter :: REAL_ELEMENT = 0.0
        integer, parameter :: INTEGER_ELEMENT = 0
        character(len=1), parameter :: BYTE(1) = 'A'

        integer, parameter :: SIZE_OF_REAL_DP = SIZE( TRANSFER( REAL_DP_ELEMENT, ['a'] ) )
        integer, parameter :: SIZE_OF_REAL = SIZE( TRANSFER( REAL_ELEMENT, ['a'] ) )
        integer, parameter :: SIZE_OF_INTEGER = SIZE( TRANSFER( INTEGER_ELEMENT, ['a'] ) )

        logical :: b, passing

        print*, "Cuda Kernels test:  test_fortcuda_kernels "
        sz = NTHREADS
        ! allocate(h_A(sz), h_B(sz), h_C(sz))
        !..........................................
        ! Create the streams and events for timing
        !..........................................
        stream = cudaStream_t(C_NULL_PTR)
        all_streams = cudaStream_t(C_NULL_PTR)
        e_start = cudaEvent_t(C_NULL_PTR)
        e_stop = cudaEvent_t(C_NULL_PTR)
#define SIZE_OF( element ) SIZE( TRANSFER( element, ['a'] ) )

        print *, 'SIZE_OF_REAL_DP = ', SIZE_OF_REAL_DP, SIZE_OF( REAL_DP_ELEMENT )
        print *, 'SIZE_OF_REAL = ', SIZE_OF_REAL, SIZE_OF( REAL_ELEMENT )
        print *, 'SIZE_OF_INTEGER = ' , SIZE_OF_INTEGER, SIZE_OF( INTEGER_ELEMENT )

        HANDLE_CUDAERROR( cudaGetDeviceCount( device_count ) , passing, 'cudaGetDeviceCount failed ' )
        do i = 0, device_count - 1
            HANDLE_CUDAERROR( cudaGetDeviceProperties( prop, i ), passing, 'cudaGetProperties failed ' )
            print *, '-------------------------------------------'
            print *, 'Device #', i
            print *, 'Total Global Memory =', prop % totalGlobalMem
            print *, 'Max Threads Per Block =', prop % maxThreadsPerBlock
            print *, 'Max Grid Size =', prop % maxGridSize( 1:3 )
        enddo

        NTHREADS_PER_BLOCK = MAX( NTHREADS_PER_BLOCK, prop % maxThreadsPerBlock )
        NUM_BLOCKS = MIN( prop % maxGridSize( 1 ), &
            & ( NTHREADS + NTHREADS_PER_BLOCK - 1 ) / NTHREADS_PER_BLOCK )

        print *, ''
        print *, '-----------------------------------------'
        print *, '  Dimensions for this Example'
        print *, '-----------------------------------------'
        print *, 'NTHREADS =', NTHREADS
        print *, 'NTHREADS_PER_BLOCK =', NTHREADS_PER_BLOCK
        print *, 'NUM_BLOCKS = ', NUM_BLOCKS
        print *, '-----------------------------------------'

        HANDLE_CUDAERROR( cudaSetDevice( 0 ) , passing, 'cudaSetDevice failed ')
        HANDLE_CUDAERROR( cudaEventCreate( e_start ), passing, 'cudaEventCreate  failed' )
        HANDLE_CUDAERROR( cudaEventCreate( e_stop ) , passing, 'cudaEventCreate failed ')

        do i = 1, num_streams
            HANDLE_CUDAERROR( cudaStreamCreate( stream(i) ), passing , 'cudaStreamCreate failed')
            HANDLE_CUDAERROR( cudaMalloc( d_A(i), sz*4 ), passing, 'cudaMalloc failed' )
            HANDLE_CUDAERROR( cudaMalloc( d_B(i), sz*4 ), passing, 'cudaMalloc failed'  )
            HANDLE_CUDAERROR( cudaMalloc( d_C(i), sz*4 ), passing, 'cudaMalloc failed'  )
        end do


        nblk = dim3( NUM_BLOCKS ,1, 1 )
        nthrd = dim3( NTHREADS_PER_BLOCK, 1, 1 )

        do i = 1, num_streams
            !........................................
            ! Create pinned host memory
            !........................................
            HANDLE_CUDAERROR( cudaMallocHost( h_c_A(i), sz*4 ) , passing, 'cudaMallocHost failed' )
            HANDLE_CUDAERROR( cudaMallocHost( h_c_B(i), sz*4 ) , passing, 'cudaMallocHost failed' )
            HANDLE_CUDAERROR( cudaMallocHost( h_c_C(i), sz*4 ) , passing, 'cudaMallocHost failed' )

            call C_F_POINTER( h_c_A(i), h_A(i) % host_ptr, [sz] )
            call C_F_POINTER( h_c_B(i), h_B(i) % host_ptr, [sz] )
            call C_F_POINTER( h_c_C(i), h_C(i) % host_ptr, [sz] )

            h_A(i) % host_ptr = 10.2
            h_B(i) % host_ptr = 20.1
            h_C(i) % host_ptr  = 0.

        enddo

        HANDLE_CUDAERROR( cudaEventRecord( e_start, all_streams ), passing , 'cudaEventRecord  failed')

        do i = 1, num_streams

            err = cudaMemcpyAsync( d_A(i), &
                & c_loc( h_A(i) % host_ptr(1) ), &
                & sz*4, &
                & cudaMemCpyHostToDevice, &
                & stream(i) )
            HANDLE_CUDAERROR( err , passing, 'cudaMemcpyAsync failed')
            err = cudaMemcpyAsync( d_B(i), &
                & c_loc( h_B(i) % host_ptr(1) ), &
                & sz*4, &
                & cudaMemCpyHostToDevice, &
                & stream(i) )
            HANDLE_CUDAERROR( err , passing, 'cudaMemcpyAsync failed')
            err = cudaMemcpyAsync( d_C(i), &
                & c_loc( h_C(i) % host_ptr(1) ), &
                & sz*4, &
                & cudaMemCpyHostToDevice, &
                & stream(i) )
            HANDLE_CUDAERROR( err, passing , 'cudaMemcpyAsync failed')

        enddo

        do i = 1, num_streams

            call vecAddF(d_A(i), d_B(i), d_C(i), nblk, nthrd, NTHREADS, stream(i))

        enddo

        do i = 1, num_streams

            err = cudaMemcpyAsync( c_loc( h_C(i) % host_ptr(1) ), &
                & d_C(i), &
                & sz*4, &
                & cudaMemCpyDeviceToHost, &
                & stream(i) )
            HANDLE_CUDAERROR( err, passing, 'cudaMemcpyAsync failed' )
            HANDLE_CUDAERROR( cudaStreamQuery( stream(i) ) , passing, 'cudaStreamQuery failed')

        end do

        HANDLE_CUDAERROR( cudaEventRecord( e_stop, all_streams ) , passing, 'cudaEventRecord  failed')
        HANDLE_CUDAERROR( cudaEventSynchronize( e_stop ) , passing ,  'cudaEventSynchronize failed')

        do i = 1, num_streams
            HANDLE_CUDAERROR( cudaStreamQuery( stream(i) ), passing , 'cudaStreamQuery failed')
        enddo

        HANDLE_CUDAERROR( cudaEventElapsedTime( c_time, e_start, e_stop ) , passing, 'cudaEventElapsedTime  failed')

        print *, 'Elapsed time = ', c_time
        print *,'*--------------------------------------------------------------------------------*'
        if ( sz < 51 ) then
            do j = 1, num_streams
                do i=1,sz
                    diff = h_A(j) % host_ptr(i) + h_B(j) % host_ptr(i) - h_C(j) % host_ptr(i)
                    print *, h_A(j) % host_ptr(i), h_B(j) % host_ptr(i), h_C(j) % host_ptr(i), diff
                enddo
            end do
        else
            do j = 1, num_streams
                passing = .true.
                diff = 0.
                do i = 1, sz
                    diff = h_A(j) % host_ptr(i) + h_B(j) % host_ptr(i) - h_C(j) % host_ptr(i)
                    if ( abs( diff ) > 1.e-6 ) then
                        print *, 'diff = ', diff
                        passing = .false.
                        exit
                    endif
                end do

                if ( .not. passing ) then
                    print *, 'Stream (', j, ') kernel resulted in incorrect values.'
                else
                    print *, 'Stream (', j, ') passinged'
                end if
            end do
        end if

        print *,'*--------------------------------------------------------------------------------*'

        HANDLE_CUDAERROR( cudaEventDestroy( e_start ) , passing, 'cudaEventDestroy  failed')
        HANDLE_CUDAERROR( cudaEventDestroy( e_stop ), passing, 'cudaEventDestroy  failed' )
        do i = 1, num_streams
            HANDLE_CUDAERROR( cudaFree( d_A(i) ), passing , 'cudaFree  failed')
            HANDLE_CUDAERROR( cudaFree( d_B(i) ), passing, 'cudaFree  failed' )
            HANDLE_CUDAERROR( cudaFree( d_C(i) ), passing, 'cudaFree  failed' )
            HANDLE_CUDAERROR( cudaStreamDestroy( stream(i)), passing , 'cudaStreamDestroy failed')

            HANDLE_CUDAERROR( cudaFreeHost( h_c_A(i) ), passing, 'cudaFreeHost  failed' )
            HANDLE_CUDAERROR( cudaFreeHost( h_c_B(i) ), passing , 'cudaFreeHost  failed')
            HANDLE_CUDAERROR( cudaFreeHost( h_c_C(i) ), passing , 'cudaFreeHost  failed')
        enddo

    end subroutine test_FortCUDA_kernels



    subroutine test_fortran_mul2dComplex_kernels
        !$ use omp_lib
        implicit none

        !define the floating point kind to be single precision
        integer, parameter :: fp_kind = kind(0.0)

        !define length of the array
        integer, parameter :: Nmin=4
        integer, parameter :: Nmax=14
        integer, parameter :: iterations= 10


        complex(fp_kind), dimension(:,:), allocatable :: c, c2
        integer(8) :: i,j, it,N, M,B
        integer(8) :: t1, t2, t3, t4, t5,t6, t7,t8, t9,t10,crate, cmax
        real(8) :: ftime, cutime, mmtime, vectime,omptime, cuda_su, cuda_per, vec_su,vec_per, omp_su, omp_per
        logical :: printc = .false.
        ftime=0.
        cutime=0.
        mmtime=0.
        vectime=0.
        print *," CUDA 2D matrix multiplication testing"
        print *," Iteration  Matrix-size     { Timing  %Reduced    SpeedupX}"
        print *,"                 Fortran(serial)                 Fortran(vectorised)         Fortran(omp)              CUDA   "
        call system_clock(count_max=cmax, count_rate=crate)
        M=1024
        ! Initialize array c, compute c2=c*c
        do it=1,iterations
            N=2**(Nmin +it)
            if(N > 2**Nmax) exit
            M=N
            allocate(c(N,M), c2(N,M))
            do i = 1, N
                do j=1, M
                    c(i,j) = cmplx(i,2*j)
                end do
            end do


            c2=cmplx(0.,0.)
            call system_clock(t1)
            do i = 1, N
                do j=1,M
                    c2(i,j)= c(i,j)*c(i,j)
                end do
            end do
            call system_clock(t2)
            ftime = ftime + REAL(t2-t1)/REAL(crate)
            if(printc)then
                print *, "Results from Fortran"
                do i = N-3, N
                    do j=M-3,M
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if
            c2=cmplx(0.,0.)
            call system_clock(t9)
            !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close)
            do i = 1, N
                do j=1,M
                    c2(i,j)= c(i,j)*c(i,j)
                end do
            end do
            !$omp end parallel do
            call system_clock(t10)
            omptime = omptime + REAL(t10-t9)/REAL(crate)
            omp_per=(REAL(t10-t9) -  REAL(t2-t1))*100/(REAL(t2-t1))
            omp_su =(REAL(t2-t1) -  REAL(t10-t9))/(REAL(t10-t9))


            if(printc)then
                print *, "Results from Fortran"
                do i = N-3, N
                    do j=M-3,M
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if




            c2=cmplx(0.,0.)
            call system_clock(t5)
            if(N==M)  c2= c*c
            call system_clock(t6)
            vectime = vectime +REAL(t6-t5)/REAL(crate)
            vec_per=(REAL(t6-t5) -  REAL(t2-t1))*100/(REAL(t2-t1))
            vec_su =(REAL(t2-t1) -  REAL(t6-t5))/(REAL(t6-t5))

            !     c2=cmplx(0.,0.)
            !     call system_clock(t7)
            ! !    c2= matmul(c,c)
            !     call system_clock(t8)
            !     mmtime = mmtime +REAL(t8-t7)/REAL(crate)

            ! Print results from Fortran
            if(printc)then
                print *, "Results from Fortran"
                do i = N-3, N
                    do j=M-3,M
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if



            ! Put
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmul2dcomplex(c,c,c2,N,M,64)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)

            !Results from CUDA
            if(printc)then
                print *, "Results from CUDA"
                do i = N-3, N
                    do j=M-3,M
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if
            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            write(*,'(I4,I10,F15.5,F15.5,F15.5,F15.5,F15.3,F15.3,F15.5,F15.5,F15.3,F15.3)')it,N,&
                &REAL(t2-t1)/REAL(crate),REAL(t6-t5)/REAL(crate),vec_per,vec_su,&
                REAL(t10-t9)/REAL(crate), omp_per, omp_su,&
                & REAL(t4-t3)/REAL(crate), cuda_per, cuda_su


            deallocate(c,c2)

        enddo
        print *, " Fortran time (serial) :", ftime
        print *, " Fortran time (vectorised):", vectime
        print *, " Fortran time (openmp):", omptime
        print *, " CUDA time    :", cutime


        print *, ""
        print *," Block_size testing"
        ftime=0.
        cutime=0.
        call system_clock(count_max=cmax, count_rate=crate)
        ! Initialize array c, compute c2=c*c
        N=(2_8)**12
        allocate(c(N,N), c2(N,N))
        do i = 1, N
            do j=1, N
                c(i,j) = cmplx(i,2*j)
            end do
        end do


        c2=cmplx(0.,0.)
        call system_clock(t1)
        do i = 1, N
            do j=1,N
                c2(i,j)= c(i,j)*c(i,j)
            end do
        end do
        call system_clock(t2)
        ftime = ftime + REAL(t2-t1)/REAL(crate)



        ! Print results from Fortran
        if(printc)then
            print *, "Results from Fortran"
            do i = N-3, N
                do j=N-3,N
                    print *,i, c(i,j),c2(i,j)
                end do
            end do
        end if

        deallocate(c,c2)
        do it=1,7
            B= 8 + 2**(it-1)
            allocate(c(N,N), c2(N,N))
            do i = 1, N
                do j=1, N
                    c(i,j) = cmplx(i,2*j)
                end do
            end do

            ! Put
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call Kernelmul2DComplex(c,c,c2,N,N,B)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)

            !Results from CUDA
            if(printc)then
                print *, "Results from CUDA"
                do i = N-3, N
                    do j=N-3,N
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if


            print *, it, B,  REAL(t2-t1)/REAL(crate), REAL(t4-t3)/REAL(crate), &
                (REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1)),   (REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))


            deallocate(c,c2)

        enddo
        print *, " Fortran time :", ftime
        print *, " CUDA time    :", cutime

    end subroutine test_fortran_mul2dComplex_kernels


end module simple_cuda_kernels
