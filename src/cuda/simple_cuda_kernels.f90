module simple_cuda_kernels
include 'simple_lib.f08'
use , intrinsic :: ISO_C_BINDING
use CUDA
implicit none
#include "simple_cuda_handle.inc"


interface
    !! c_kernels.cu
    subroutine vecAddF(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddfloat")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value :: a, b, c
        type (dim3) :: dimGrid
        type (dim3) :: dimBlk
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine vecAddF
    subroutine vecAddI(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddint")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value :: a, b, c
        type (dim3) :: dimGrid
        type (dim3) :: dimBlk
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine vecAddI
    subroutine vecAddConstF(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddconstfloat")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value :: a, b, c ! b is a const
        type (dim3) :: dimGrid
        type (dim3) :: dimBlk
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine vecAddConstF
    subroutine vecAddConstI(a, b, c, dimGrid, dimBlk, N, stream) bind(C, name="vecaddconstint")
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
    subroutine multiply_by_block(dimGrid, threads, a, b, c, N, stream) bind(C, name="multiply_by_block")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value :: a, b,c
        type (dim3) :: dimGrid
        type (dim3) :: threads
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine multiply_by_block
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

!! Mul2D.cu  -- independent 2D multiplication kernels
! interface
!     subroutine mul2dfloat( a,b, c, Np,Ns, Bs) bind(C, name="kernelmul2dfloat")
!         use, intrinsic :: ISO_C_BINDING
!         real(c_float),dimension(:,:) :: a, b,c
!         integer(c_int), value :: Np,Ns, Bs
!     end subroutine mul2dfloat
!     subroutine mul2dcomplex( a,b, c, Np,Ns, Bs) bind(C, name="kernelmul2dcomplex")
!         use, intrinsic :: ISO_C_BINDING
!         complex(C_FLOAT_COMPLEX),dimension(:,:) :: a,b,c
!         integer(c_int), value :: Np, Ns, Bs
!     end subroutine mul2dcomplex
!     subroutine muladd2dcomplex( a,b, c,d, Np,Ns, Bs) bind(C, name="kernelmul2dcomplex")
!         use, intrinsic :: ISO_C_BINDING
!         complex(C_FLOAT_COMPLEX),dimension(:,:) :: a,b,c,d
!         integer(c_int), value :: Np, Ns, Bs
!     end subroutine muladd2dcomplex! end interface

!! filter_kernels.cu
interface
    !! Obtain
    subroutine filter_gaussKernel(U,V,W,Z, sigma, dimGrid, dimBlk, N, stream) bind(C, name="gausskernelelementwise")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value ::U,V,W !grid
        type(c_ptr) :: Z ! return value
        real(c_float) ::  sigma
        type (dim3) :: dimGrid
        type (dim3) :: dimBlk
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine filter_gaussKernel
    subroutine mulgaussKernel(A,B , sigma, dimGrid, dimBlk, N, stream) bind(C, name="gaussconvolution3d")
        use, intrinsic :: ISO_C_BINDING
        use CUDA, only : dim3, cudaStream_t
        type (c_ptr), value ::A  !! B
        type(c_ptr) :: B ! return value B= GaussKernel x A
        real(c_float) ::  sigma
        type (dim3) :: dimGrid
        type (dim3) :: dimBlk
        integer(c_int), value :: N
        type (cudaStream_t) :: stream
    end subroutine mulgaussKernel
end interface

contains



    !>  \brief corr is for correlating two images
    function corr_cuda( self1, self2, lp_dyn, hp_dyn ) result( r )
        use simple_image, only: image
        type(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        complex, allocatable :: cmat1(:,:,:), cmat2(:,:,:), vec1(:), vec2(:)
        integer, allocatable :: bpmask(:,:,:)
        real    :: r, sumasq, sumbsq
        integer :: h, k, l, phys(3), lims(3,2),ldim(3), sqarg, sqlp, sqhp, npoints, npointsNP2
        logical :: didft1, didft2
        integer(kind=timer_int_kind) :: t1,t2
#include "simple_local_flags.inc"
        r = 0.
#ifdef USING_CUDA
        t1=tic()
        sumasq = 0.
        sumbsq = 0.
        if( self1.eqdims.self2 )then
            didft1 = .false.
            if( .not. self1%is_ft() )then
                call self1%fft()
                didft1 = .true.
            endif
            didft2 = .false.
            if( .not. self2%is_ft() )then
                call self2%fft()
                didft2 = .true.
            endif

            if( present(lp_dyn) )then
                lims = self1%loop_lims(1,lp_dyn)
            else
                lims = self1%loop_lims(2) ! Nyqvist default low-pass limit
            endif
            sqlp = (maxval(lims(:,2)))**2
            if( present(hp_dyn) )then
                sqhp = max(2,self1%get_find(hp_dyn))**2
            else
                sqhp = 2 ! index 2 default high-pass limit
            endif
            ldim=self1%get_ldim()
            allocate(cmat1(1:ldim(1),1:ldim(2),1:ldim(3)),cmat2(1:ldim(1),1:ldim(2),1:ldim(3)))
            allocate(bpmask(1:ldim(1),1:ldim(2),1:ldim(3)))
            bpmask=0

            !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,phys)&
            !$omp  schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sqarg = h*h + k*k + l*l
                        if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                            phys = self1%comp_addr_phys(h,k,l)
                            bpmask(h-lims(1,1)+1, k-lims(2,1)+1, l-lims(3,1)+1) = 1
                            cmat1(h-lims(1,1)+1, k-lims(2,1)+1, l-lims(3,1)+1) = self1%get_cmat_at(phys(1),phys(2),phys(3))
                            cmat2(h-lims(1,1)+1, k-lims(2,1)+1, l-lims(3,1)+1) = self2%get_cmat_at(phys(1),phys(2),phys(3))
                        else
                            cmat1(h-lims(1,1)+1, k-lims(2,1)+1, l-lims(3,1)+1) = cmplx(0.,0.)
                            cmat2(h-lims(1,1)+1, k-lims(2,1)+1, l-lims(3,1)+1) = cmplx(0.,0.)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do

            npoints =sum(sum(sum(bpmask,3),2),1)
            print *," Number of valid points ",npoints, ", in ", product(ldim)
            npointsNP2 = nextPow2(npoints)
            print *," Number of valid points (nextpow of 2)",npointsNP2
            allocate(vec1(npointsNP2), source=cmplx(0.,0.))
            allocate(vec2(npointsNP2), source=cmplx(0.,0.))
            vec1=pack(cmat1,MASK=bpmask.ne.0, VECTOR=vec1)
            vec2=pack(cmat2,MASK=bpmask.ne.0, VECTOR=vec2)
            print *, " corr_cuda initialisation ", toc(t1)
            t2 = tic()
            !! Single pass
            !           call kernelcrosscorr(cmat1, cmat2, ldim(1),ldim(2), sqlp, sqhp, )

            !! Stages
            call kernelrcorr(cmat1,cmat2, r, ldim, sqlp, sqhp, 64,8)
            print *, " corr_cuda r correlation ",r, toc(t2)
            call kernelsumcsq(cmat1, sumasq, ldim, 64, 8 )
            print *, " corr_cuda A sum squared ", sumasq, toc()
            call kernelsumcsq(cmat2, sumbsq, ldim, 64, 8 )
            print *, " corr_cuda A sum squared ",sumbsq,  toc()


            if( sumasq < TINY .or. sumbsq < TINY )then
                r = 0.
            else
                r = r / sqrt(sumasq * sumbsq)
            endif

            print *, " Image corr ", r, " took ", toc(t2)

        else
            write(*,*) 'self1%ldim:', self1%get_ldim()
            write(*,*) 'self2%ldim:', self2%get_ldim()
            THROW_HARD('images to be correlated need to have same dimensions')
        endif
        if(allocated(cmat1)) deallocate(cmat1)
        if(allocated(cmat2)) deallocate(cmat2)
        if(allocated(bpmask)) deallocate(bpmask)
        if(allocated(vec1)) deallocate(vec1)
        if(allocated(vec2)) deallocate(vec2)

#else
        print *, " CUDA corr in image class not available for this build"
#endif
    end function corr_cuda

    subroutine test_FortCUDA_kernels (arg)
        implicit none
        real, intent(in) :: arg

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


    subroutine test_fortran_mul1dComplex_kernels
        !$ use omp_lib
        implicit none

        !define the floating point kind to be single precision
        integer, parameter :: fp_kind = kind(0.0)

        !define length of the array
        integer, parameter :: Nmin=14
        integer, parameter :: Nmax=24
        integer, parameter :: iterations= 30


        complex(fp_kind), dimension(:), allocatable :: c, c2
        complex(fp_kind) :: cserial
        integer(8) :: i, it,iloop
        integer :: N,  B
        integer(8) :: t1, t2, t3, t4, t5,t6, t7,t8, t9,t10,crate, cmax
        real(8) :: ftime, cutime, mmtime, vectime,omptime, cuda_su, cuda_per, vec_su,vec_per, omp_su, omp_per
        logical :: printc = .false.
        ftime=0.
        cutime=0.
        mmtime=0.
        vectime=0.

        print *," CUDA 1D matrix multiplication testing "
        print *," Iteration  N                          { Timing      %Reduced     SpeedupX }"

        print *,"                      Fortran(serial)   Fortran(omp)                                   CUDA   "
        call system_clock(count_max=cmax, count_rate=crate)

        ! Initialize array c, compute c2=c*c
        do iloop=1,30
            N=2**(Nmin +iloop-1)
            if(N > 2**Nmax) exit
            allocate(c(N), c2(N))
            do i = 1, N
                c(i) = cmplx(i,2*real(i)/real(N))
            end do

            !! SERIAL
            do  it=1,iterations
                c2=cmplx(0.,0.)
                call system_clock(t1)
                do i = 1, N
                    c2(i)= c(i)*c(i)
                end do
                call system_clock(t2)

                if(it/=1)ftime = ftime + REAL(t2-t1)/REAL(crate)
            end do
            cserial = c2(N-1)
            if(printc)then
                print *, "Results from Fortran"
                do i = N-3, N

                    print *,i, c(i),c2(i)

                end do
            end if

            !! OpenMP
            do  it=1,iterations
                c2=cmplx(0.,0.)
                call system_clock(t9)
                !$omp parallel do  private(i) default(shared) proc_bind(close)
                do i = 1, N

                    c2(i)= c(i)*c(i)

                end do
                !$omp end parallel do
                call system_clock(t10)
                if(it/=1)omptime = omptime + REAL(t10-t9)/REAL(crate)
            end do
            omptime = omptime + REAL(t10-t9)/REAL(crate)


            omp_per=(REAL(t10-t9) -  REAL(t2-t1))*100/(REAL(t2-t1))
            omp_su =(REAL(t2-t1) -  REAL(t10-t9))/(REAL(t10-t9))
            if(cserial .ne. c2(N-1))then
                print *,' OpenMP [',N,N, '] failed'
                do i =1, N

                    if(abs(c2(i)) == 0.)then
                        print *,i, c(i),c2(i),  c(i)*c(i)
                        stop
                    endif


                end do
                stop
            endif
            if(printc)then
                print *, "Results from Fortran"
                do i = N-3, N

                    print *,i, c(i),c2(i)

                end do
            end if

            ! CUDA
            do  it=1,iterations
                c2=cmplx(0.,0.)
                call system_clock(t3)
                ! Do the same computation with CUDA.
                ! Fortran -> C -> CUDA ->C ->Fortran
                call kernelmul1dcomplex(c,c,c2,N,8)
                call system_clock(t4)
                if(it/=1)cutime = cutime + REAL(t4-t3)/REAL(crate)

            end do
            if(cserial .ne. c2(N-1))then
                print *,' CUDA [',N,'] failed', cserial ,c2(N-1)
                do i =1, N
                    if(abs(c2(i)) == 0.)then
                        print *,i, c(i),c2(i),  c(i)*c(i)
                        stop
                    endif
                end do
                stop
                print *, c2
                stop
            endif

            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            write(*,'(4x,I4,I10,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5)')iloop,N,&
                &REAL(t2-t1)/REAL(crate),&
                &REAL(t10-t9)/REAL(crate), omp_per, omp_su,&
                & REAL(t4-t3)/REAL(crate), cuda_per, cuda_su


            deallocate(c,c2)

        enddo
        print *, " Fortran time (serial) :", ftime
        !   print *, " Fortran time (vectorised):", vectime
        print *, " Fortran time (openmp):", omptime
        print *, " CUDA time    :", cutime



        ! print *, ""
        ! print *," Block_size testing  1D array"
        ! ftime=0.
        ! cutime=0.
        ! call system_clock(count_max=cmax, count_rate=crate)
        ! ! Initialize array c, compute c2=c*c
        ! N=(2_8)**18
        ! allocate(c(N), c2(N))
        ! do i = 1, N
        !     c(i) = cmplx(i,2*real(i)/real(N))
        ! end do

        ! !! SERIAL
        ! c2=cmplx(0.,0.)
        ! call system_clock(t1)
        ! do i = 1, N
        !     c2(i)= c(i)*c(i)
        ! end do
        ! call system_clock(t2)
        ! ftime = ftime + REAL(t2-t1)/REAL(crate)
        ! cserial = c2(N-1)




        ! deallocate(c,c2)
        ! do it=1,7
        !     B= 2**(it+2)
        !     allocate(c(N), c2(N))
        !     do i = 1, N
        !         c(i) = cmplx(i,2*real(i)/real(N))
        !     end do

        !     ! Put
        !     c2=cmplx(0.,0.)
        !     call system_clock(t3)
        !     ! Do the same computation with CUDA.
        !     ! Fortran -> C -> CUDA ->C ->Fortran
        !     call kernelmul2dcomplex(c,c,c2,N,B)
        !     call system_clock(t4)
        !     cutime = cutime + REAL(t4-t3)/REAL(crate)
        !     if(cserial .ne. c2(N-1)) then
        !         print *,' CUDA N=',N,' failed', cserial, c2(N-1), N-1
        !         stop
        !     endif


        !     print *, it+2, B,  REAL(t2-t1)/REAL(crate), REAL(t4-t3)/REAL(crate), &
        !         (REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1)),   (REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))


        !     deallocate(c,c2)

        ! enddo
        ! print *, " Fortran time :", ftime
        ! print *, " CUDA time    :", cutime

    end subroutine test_fortran_mul1dComplex_kernels



    subroutine test_fortran_squaremul2dComplex_kernels
        !$ use omp_lib
        implicit none

        !define the floating point kind to be single precision
        integer, parameter :: fp_kind = kind(0.0)

        !define length of the array
        integer, parameter :: NNmin=8
        integer, parameter :: NNmax=14
        integer, parameter :: iterations= 30


        complex(fp_kind), dimension(:,:), allocatable :: c, c2
        complex(fp_kind) :: cserial
        integer(8) :: i,j, it
        integer :: N,  B
        integer(8) :: t1, t2, t3, t4, t5,t6, t7,t8, t9,t10,crate, cmax
        real(8) :: ftime, cutime, mmtime, vectime,omptime, cuda_su, cuda_per, vec_su,vec_per, omp_su, omp_per
        logical :: printc = .false.
        ftime=0.
        cutime=0.
        mmtime=0.
        vectime=0.


        print *," CUDA 2D matrix multiplication testing -- Square 2D matrix"
        print *," Iteration  NxN                       { Timing      %Reduced      SpeedupX }"

        print *,"                      Fortran(serial)  Fortran(omp)                                   CUDA   "
        call system_clock(count_max=cmax, count_rate=crate)

        ! Initialize array c, compute c2=c*c
        do it=1,iterations
            N=2**(NNmin +it-1)
            if(N > 2**NNmax) exit

            allocate(c(N,N), c2(N,N))
            do i = 1, N
                do j=1, N
                    c(i,j) = cmplx(i,2*real(j)/real(N))
                end do
            end do

            !! SERIAL
            c2=cmplx(0.,0.)
            call system_clock(t1)
            do i = 1, N
                do j=1,N
                    c2(i,j)= c(i,j)*c(i,j)
                end do
            end do
            call system_clock(t2)
            cserial = c2(N-1,N-1)
            ftime = ftime + REAL(t2-t1)/REAL(crate)
            ! if(printc)then
            !     print *, "Results from Fortran Serial"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if

            !! OpenMP
            c2=cmplx(0.,0.)
            call system_clock(t9)
            !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close)
            do i = 1, N
                do j=1,N
                    c2(i,j)= c(i,j)*c(i,j)
                end do
            end do
            !$omp end parallel do
            call system_clock(t10)
            omptime = omptime + REAL(t10-t9)/REAL(crate)
            omp_per=(REAL(t10-t9) -  REAL(t2-t1))*100/(REAL(t2-t1))
            omp_su =(REAL(t2-t1) -  REAL(t10-t9))/(REAL(t10-t9))
            if(cserial .ne. c2(N-1,N-1)) print *,' OpenMP N=',N, ' failed ', cserial, c2(N-1,N-1)


            ! if(printc)then
            !     print *, "Results from Fortran"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if




            ! CUDA
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmul2dcomplex(c,c,c2,N,N,8)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,N-1)) then
                print *,' CUDA N=',N,' failed', cserial, c2(N-1,N-1), N-1,N-1
                !  print * ,cmplx(N-1,2*(N-1)), cmplx(N-1,2*(N-1))*cmplx(N-1,2*(N-1))
                print *, c2
                stop
            endif


            ! !Results from CUDA
            ! if(printc)then
            !     print *, "Results from CUDA"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if
            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            write(*,'(4x,I4,I10,F15.5,F15.5,F15.5,F15.5,F15.3,F15.5,F15.5,F15.5)')it,N,&
                &REAL(t2-t1)/REAL(crate),&
                REAL(t10-t9)/REAL(crate), omp_per, omp_su,&
                & REAL(t4-t3)/REAL(crate), cuda_per, cuda_su


            deallocate(c,c2)

        enddo
        print *, " Fortran time (serial) :", ftime
        print *, " Fortran time (openmp):", omptime
        print *, " CUDA time    :", cutime


        print *, ""
        print *," Block_size testing  square 2D array"
        ftime=0.
        cutime=0.
        call system_clock(count_max=cmax, count_rate=crate)
        ! Initialize array c, compute c2=c*c
        N=(2_8)**12
        allocate(c(N,N), c2(N,N))
        do i = 1, N
            do j=1, N
                c(i,j) = cmplx(i,2*real(j)/real(N))
            end do
        end do

        !! SERIAL
        c2=cmplx(0.,0.)
        call system_clock(t1)
        do i = 1, N
            do j=1,N
                c2(i,j)= c(i,j)*c(i,j)
            end do
        end do
        call system_clock(t2)
        ftime = ftime + REAL(t2-t1)/REAL(crate)
        cserial = c2(N-1,N-1)


        ! ! Print results from Fortran
        ! if(printc)then
        !     print *, "Results from Fortran"
        !     do i = N-3, N
        !         do j=N-3,N
        !             print *,i, c(i,j),c2(i,j)
        !         end do
        !     end do
        ! end if

        deallocate(c,c2)
        do it=1,4
            B= 2**(it+2)
            allocate(c(N,N), c2(N,N))
            do i = 1, N
                do j=1, N
                    c(i,j) = cmplx(i,2*real(j)/real(N))
                end do
            end do

            ! Put
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmul2dcomplex(c,c,c2,N,N,B)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,N-1)) then
                print *,' CUDA N=',N,' failed', cserial, c2(N-1,N-1), N-1,N-1
                stop
            endif

            ! !Results from CUDA
            ! if(printc)then
            !     print *, "Results from CUDA"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if
            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            print *, it+2, B,  REAL(t2-t1)/REAL(crate), REAL(t4-t3)/REAL(crate), &
                cuda_per,   cuda_su


            deallocate(c,c2)

        enddo

        ! print *, " Fortran time :", ftime
        ! print *, " CUDA time    :", cutime

    end subroutine test_fortran_squaremul2dComplex_kernels

    subroutine test_fortran_squaremuladd2dComplex_kernels
        !$ use omp_lib
        implicit none

        !define the floating point kind to be single precision
        integer, parameter :: fp_kind = kind(0.0)

        !define length of the array
        integer, parameter :: NNmin=8
        integer, parameter :: NNmax=14
        integer, parameter :: iterations= 30


        complex(fp_kind), dimension(:,:), allocatable :: c, c2
        complex(fp_kind) :: cserial
        integer(8) :: i,j, it
        integer :: N,  B
        integer(8) :: t1, t2, t3, t4, t5,t6, t7,t8, t9,t10,crate, cmax
        real(8) :: ftime, cutime, mmtime, vectime,omptime, cuda_su, cuda_per, vec_su,vec_per, omp_su, omp_per
        logical :: printc = .false.
        ftime=0.
        cutime=0.
        mmtime=0.
        vectime=0.


        print *," CUDA 2D matrix muladd testing -- Square 2D matrix"
        print *," Iteration  NxN                       { Timing      %Reduced      SpeedupX }"

        print *,"                      Fortran(serial)  Fortran(omp)                                   CUDA   "
        call system_clock(count_max=cmax, count_rate=crate)

        ! Initialize array c, compute c2=c*c
        do it=1,iterations
            N=2**(NNmin +it-1)
            if(N > 2**NNmax) exit

            allocate(c(N,N), c2(N,N))
            do i = 1, N
                do j=1, N
                    c(i,j) = cmplx(i,2*real(j)/real(N))
                end do
            end do

            !! SERIAL
            c2=cmplx(0.,0.)
            call system_clock(t1)
            do i = 1, N
                do j=1,N
                    c2(i,j)= c(i,j)*c(i,j) + c(i,j)
                end do
            end do
            call system_clock(t2)
            cserial = c2(N-1,N-1)
            ftime = ftime + REAL(t2-t1)/REAL(crate)
            ! if(printc)then
            !     print *, "Results from Fortran Serial"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if

            !! OpenMP
            c2=cmplx(0.,0.)
            call system_clock(t9)
            !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close)
            do i = 1, N
                do j=1,N
                    c2(i,j)= c(i,j)*c(i,j) + c(i,j)
                end do
            end do
            !$omp end parallel do
            call system_clock(t10)
            omptime = omptime + REAL(t10-t9)/REAL(crate)
            omp_per=(REAL(t10-t9) -  REAL(t2-t1))*100/(REAL(t2-t1))
            omp_su =(REAL(t2-t1) -  REAL(t10-t9))/(REAL(t10-t9))
            if(cserial .ne. c2(N-1,N-1)) print *,' OpenMP N=',N, ' failed ', cserial, c2(N-1,N-1)


            ! if(printc)then
            !     print *, "Results from Fortran"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if




            ! CUDA
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmuladd2dcomplex(c,c,c,c2,N,N,8)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,N-1)) then
                print *,' CUDA N=',N,' failed', cserial, c2(N-1,N-1), N-1,N-1
                !  print * ,cmplx(N-1,2*(N-1)), cmplx(N-1,2*(N-1))*cmplx(N-1,2*(N-1))
                print *, c2
                stop
            endif


            ! !Results from CUDA
            ! if(printc)then
            !     print *, "Results from CUDA"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if
            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            write(*,'(4x,I4,I10,F15.5,F15.5,F15.5,F15.5,F15.3,F15.5,F15.5,F15.5)')it,N,&
                &REAL(t2-t1)/REAL(crate),&
                REAL(t10-t9)/REAL(crate), omp_per, omp_su,&
                & REAL(t4-t3)/REAL(crate), cuda_per, cuda_su


            deallocate(c,c2)

        enddo
        print *, " Fortran time (serial) :", ftime
        print *, " Fortran time (openmp):", omptime
        print *, " CUDA time    :", cutime


        print *, ""
        print *," Block_size testing  square 2D array"
        ftime=0.
        cutime=0.
        call system_clock(count_max=cmax, count_rate=crate)
        ! Initialize array c, compute c2=c*c
        N=(2_8)**12
        allocate(c(N,N), c2(N,N))
        do i = 1, N
            do j=1, N
                c(i,j) = cmplx(i,2*real(j)/real(N))
            end do
        end do

        !! SERIAL
        c2=cmplx(0.,0.)
        call system_clock(t1)
        do i = 1, N
            do j=1,N
                c2(i,j)= c(i,j)*c(i,j)
            end do
        end do
        call system_clock(t2)
        ftime = ftime + REAL(t2-t1)/REAL(crate)
        cserial = c2(N-1,N-1)


        ! ! Print results from Fortran
        ! if(printc)then
        !     print *, "Results from Fortran"
        !     do i = N-3, N
        !         do j=N-3,N
        !             print *,i, c(i,j),c2(i,j)
        !         end do
        !     end do
        ! end if

        deallocate(c,c2)
        do it=1,4
            B= 2**(it+2)
            allocate(c(N,N), c2(N,N))
            do i = 1, N
                do j=1, N
                    c(i,j) = cmplx(i,2*real(j)/real(N))
                end do
            end do

            ! Put
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmuladd2dcomplex(c,c,c,c2,N,N,B)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,N-1)) then
                print *,' CUDA N=',N,' failed', cserial, c2(N-1,N-1), N-1,N-1
                stop
            endif

            ! !Results from CUDA
            ! if(printc)then
            !     print *, "Results from CUDA"
            !     do i = N-3, N
            !         do j=N-3,N
            !             print *,i, c(i,j),c2(i,j)
            !         end do
            !     end do
            ! end if
            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            print *, it+2, B,  REAL(t2-t1)/REAL(crate), REAL(t4-t3)/REAL(crate), &
                cuda_per,   cuda_su


            deallocate(c,c2)

        enddo

        ! print *, " Fortran time :", ftime
        ! print *, " CUDA time    :", cutime

    end subroutine test_fortran_squaremuladd2dComplex_kernels



    subroutine test_fortran_mul2dComplex_kernels
        !$ use omp_lib
        implicit none

        !define the floating point kind to be single precision
        integer, parameter :: fp_kind = kind(0.0)

        !define length of the array
        integer, parameter :: Nmin=2
        integer, parameter :: NNmin=2
        integer, parameter :: Nmax=16
        integer, parameter :: NNmax=14
        integer, parameter :: iterations= 10


        complex(fp_kind), dimension(:,:), allocatable :: c, c2
        complex(fp_kind) :: cserial
        integer(8) :: i,j, it
        integer :: N, M, B
        integer(8) :: t1, t2, t3, t4, t5,t6, t7,t8, t9,t10,crate, cmax
        real(8) :: ftime, cutime, mmtime, vectime,omptime, cuda_su, cuda_per, vec_su,vec_per, omp_su, omp_per
        logical :: printc = .false.
        ftime=0.
        cutime=0.
        mmtime=0.
        vectime=0.

        print *," CUDA 2D matrix multiplication testing -- non-square array "
        print *," Iteration  Nx1024                     { Timing        %Reduced       SpeedupX }"

        print *,"                      Fortran(serial)  Fortran(OpenMP)                                CUDA   "
        call system_clock(count_max=cmax, count_rate=crate)
        M=1024
        ! Initialize array c, compute c2=c*c
        do it=1,iterations
            N=2**(Nmin +it-1) + 1
            if(N > 2**Nmax) exit
            allocate(c(N,M), c2(N,M))
            do i = 1, N
                do j=1, M
                    c(i,j) = cmplx(i,2*real(j)/real(M))
                end do
            end do

            !! SERIAL
            c2=cmplx(0.,0.)
            call system_clock(t1)
            do i = 1, N
                do j=1,M
                    c2(i,j)= c(i,j)*c(i,j)
                end do
            end do
            call system_clock(t2)
            cserial = c2(N-1,M-1)
            ftime = ftime + REAL(t2-t1)/REAL(crate)
            !! OpenMP
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
            if(cserial .ne. c2(N-1,M-1))then
                print *,' OpenMP [',N,M, '] failed'
                do i =1, N
                    do j=1,M
                        if(abs(c2(i,j)) == 0.)then
                            print *,i,j, c(i,j),c2(i,j),  c(i,j)*c(i,j)
                            stop
                        endif
                    end do
                end do
                stop
            endif

            ! CUDA
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmul2dcomplex(c,c,c2,N,M,8)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,M-1))then
                print *,' CUDA [',N,M,'] failed', cserial ,c2(N-1,M-1)
                do i =1, N
                    do j=1,M
                        if(abs(c2(i,j)) == 0.)then
                            print *,i,j, c(i,j),c2(i,j),  c(i,j)*c(i,j)
                            stop
                        endif
                    end do
                end do
                stop

            endif

            cuda_per=(REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1))
            cuda_su =(REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))

            write(*,'(4x,I4,I10,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5)')it,N,&
                &REAL(t2-t1)/REAL(crate),&
                REAL(t10-t9)/REAL(crate), omp_per, omp_su,&
                & REAL(t4-t3)/REAL(crate), cuda_per, cuda_su


            deallocate(c,c2)

        enddo
        print *, " Fortran time (serial) :", ftime
        !   print *, " Fortran time (vectorised):", vectime
        print *, " Fortran time (openmp):", omptime
        print *, " CUDA time    :", cutime


        print *, ""
        print *," Block_size testing  non-square 2D array"
        ftime=0.
        cutime=0.
        call system_clock(count_max=cmax, count_rate=crate)
        ! Initialize array c, compute c2=c*c
        N=(2_8)**12
        M=(2_8)**8
        allocate(c(N,M), c2(N,M))
        do i = 1, N
            do j=1, M
                c(i,j) = cmplx(i,2*real(j)/real(M))
            end do
        end do

        !! SERIAL
        c2=cmplx(0.,0.)
        call system_clock(t1)
        do i = 1, N
            do j=1,M
                c2(i,j)= c(i,j)*c(i,j)
            end do
        end do
        call system_clock(t2)
        ftime = ftime + REAL(t2-t1)/REAL(crate)
        cserial = c2(N-1,M-1)


        ! Print results from Fortran
        if(printc)then
            print *, "Results from Fortran"
            do i = N-3, N
                do j=M-3,M
                    print *,i, c(i,j),c2(i,j)
                end do
            end do
        end if

        deallocate(c,c2)
        do it=1,7
            B= 2**(it+2)

            allocate(c(N,M), c2(N,M))
            do i = 1, N
                do j=1, M
                    c(i,j) = cmplx(i,2*real(j)/real(M))
                end do
            end do

            ! Put
            c2=cmplx(0.,0.)
            call system_clock(t3)
            ! Do the same computation with CUDA.
            ! Fortran -> C -> CUDA ->C ->Fortran
            call kernelmul2dcomplex(c,c,c2,N,M,B)
            call system_clock(t4)
            cutime = cutime + REAL(t4-t3)/REAL(crate)
            if(cserial .ne. c2(N-1,M-1)) then
                print *,' CUDA N=',N,' failed', cserial, c2(N-1,M-1), N-1,M-1
                stop
            endif

            !Results from CUDA
            if(printc)then
                print *, "Results from CUDA"
                do i = N-3, N
                    do j=N-3,N
                        print *,i, c(i,j),c2(i,j)
                    end do
                end do
            end if

            print *, it+2, B,  REAL(t2-t1)/REAL(crate), REAL(t4-t3)/REAL(crate), &
                (REAL(t4-t3) -  REAL(t2-t1))*100/(REAL(t2-t1)),   (REAL(t2-t1) -  REAL(t4-t3))/(REAL(t4-t3))


            deallocate(c,c2)

        enddo
        print *, " Fortran time :", ftime
        print *, " CUDA time    :", cutime

    end subroutine test_fortran_mul2dComplex_kernels


end module simple_cuda_kernels
