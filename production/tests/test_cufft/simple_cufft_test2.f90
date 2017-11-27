module simple_cufft_test2
    include 'simple_lib.f08'
    use simple_image,       only: image
    use simple_fftshifter
    use gnufor2
    use simple_cufft
    implicit none

#ifdef PGI
    public ::  test_pgi_ifft, test_pgi, test_pgi1,&
        test_pgi2, test_pgi3, test_pgi4
#endif
    private

    ! module global constants
    integer, parameter :: SQRAD=60, NTST=50, NNOISY=20
    real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2

    real(kind=sp), parameter :: tolerance = 1.e-05_sp
#include "simple_local_flags.inc"
contains

    !> \brief fft  forward Fourier transform
    !!
    subroutine fft_pgi_cuda_test( self )
        ! use simple_fftshifter
        use openacc
        use cudafor
        !! IN-PLACE 3D FFT using CUDA's cufft
        use simple_cufft
        implicit none
        type(image), intent(inout)      :: self
        real(sp),    allocatable         :: rinput(:,:,:)
        !real(sp),    allocatable, device :: rinput_d(:,:,:)
        !complex(sp), allocatable, device :: coutput_d(:,:,:)
        complex(sp), allocatable, device :: inplace_d(:,:,:) !< inplace real and complex FFT data
        complex(sp), allocatable         :: coutput(:,:,:)
        integer      :: plan, planType, ierr
        integer      :: i,j,k,h,l,n,istat
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        integer(timer_int_kind) :: t1


        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::fft start ",__FILENAME__,__LINE__)
        endif


        verbose=.false.
        if(verbose) t1=tic()
        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        nscale =  real(product(ldim))
        !   lims = self%loop_lims(3)
        VerbosePrint "In simple_cufft::fft_pgi_cuda_test  PGI test ldim ",ldim
        VerbosePrint "In simple_cufft::fft_pgi_cuda_test  PGI test cdim ", cdim
        VerbosePrint "In simple_cufft::fft_pgi_cuda_test  PGI test nc ", nscale
        !    VerbosePrint "In simple_cufft::fft_pgi  PGI test loop lims ", lims
        ierr=0
        ! allocate arrays on the host
        !  if(allocated(coutput))deallocate(coutput)
        allocate(coutput(ldim(1),ldim(2),ldim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat /= 0)call allocchk("In fft_cuda_test coutput",istat)
        !    if(allocated(rinput))deallocate(rinput)
        !  allocate( rinput(ldim(1), ldim(2), ldim(3)),source=0.,stat=istat)
        !    if(istat.ne.0)call allocchk("In simple_image::pgi_fft rinput",istat)

        allocate( rinput(ldim(1), ldim(2), ldim(3)),stat=istat)
        rinput = self%get_rmat()
        if(istat /= 0)call allocchk("In fft_cuda_test rinput ",istat)

        ! call gnufor_image(rinput(:ldim(1),:ldim(2),1),palette='gray')
        !coutput(:ldim(1), :ldim(2), :ldim(3))=cmplx(rinput(:ldim(1), :ldim(2), :ldim(3)),0.)

        !call fftshift(rinput)
        !VerbosePrint "In simple_cufft::fft_pgi  PGI test ", shape(coutput), shape(rinput)
        ! istat=cudaDeviceSynchronize()
        ! istat=istat+cudaThreadSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::fft sync ",__FILENAME__,__LINE__)
        ! endif
        !        call fftshift(rinput)
        VerbosePrint "In simple_cufft::fft_pgi  PGI start "
        ! if(.not.allocated(coutput)) call simple_stop(' coutput failed ')
        ! if(.not.allocated(rinput)) call simple_stop(' rinput failed ')
        ! allocate arrays on the device
        allocate(inplace_d(ldim(1),ldim(2),ldim(3)),stat=istat)
        if(istat /= 0) call allocchk("fft_cuda_test inplace_d",stat)
        !if(allocated(rinput_d))deallocate( rinput_d)
        !allocate( rinput_d(ldim(1),ldim(2),ldim(3)),source=0.,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::fft rinput_d",istat)
        ! if(allocated(coutput_d))deallocate(coutput_d)
        !allocate(coutput_d(cdim(1),cdim(2),cdim(3)),stat=istat)
        !if(istat /= 0) call allocchk("In simple_image::fft rinput_d",stat)
        ! set planType to either single or double precision
        planType = CUFFT_C2C    !! (fp_kind == singlePrecision)
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI device allocated"

        !copy arrays to device
        inplace_d = cmplx(rinput) !(1:ldim(1)-1:2,:ldim(2),:ldim(3)),rinput(2:ldim(1):2,:ldim(2),:ldim(3)))
        !rinput_d = rinput
        VerbosePrint "In simple_cufft::fft_pgi  PGI  arrays copied to device "


        ! Initialize the plan for real to complex transform
        call cufftPlan3D(plan,ldim(1),ldim(2),ldim(3), planType)
        VerbosePrint "In fft_pgi_cuda_test  PGI cufftPlan completed"

        ! Execute  Forward transform in place
        !call cufftExec(plan,planType,rinput_d,coutput_d,CUFFT_FORWARD )
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_FORWARD )
        VerbosePrint "In fft_pgi_cuda_test  PGI cufftexec completed "

        ! now scale the values so that a ifft of the output yields the
        ! original image back, rather than a scaled version
        !$acc kernels
        inplace_d = inplace_d / nscale
        !$acc end kernels
        VerbosePrint "In simple_cufft::fft_pgi  PGI acc kernel completed"

        ! Copy results back to host
        coutput = inplace_d
        VerbosePrint "In fft_pgi_cuda_test PGI copied results back to host"
        ! istat=cudaThreadSynchronize()
        ! istat=istat+cudaDeviceSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
        ! endif
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI synchronized "
        ! call ifftshift(coutput,lims)
        ! VerbosePrint " COUTPUT ", real(coutput(1:3,1,1)), cabs(coutput(1:3,1,1))
        VerbosePrint "In fft_pgi_cuda_test  set cmat ", shape(coutput), LBOUND(coutput), UBOUND(coutput)
        ! call self%set_cmat(coutput(:ldim(1),:ldim(2),:ldim(3)))
        ! VerbosePrint "In simple_cufft::fft_pgi  PGI destroyed "

        !      call self%set_cmat( coutput )

        call self%set_ft(.true.)
        call self%set_cmat(coutput(:cdim(1),:cdim(2),:cdim(3)))
        ! release memory on the host and on the device
        !  if(allocated(rinput_d))then
        !      deallocate(rinput_d)!,stat=istat)
        !if(istat.ne.0)call allocchk("In simple_image::fft dealloc device coutput_d ",istat)
        !  endif
        !  if(allocated(coutput_d))then
        istat=0
        ! deallocate(coutput_d, rinput_d,stat=istat)
        deallocate(inplace_d,stat=istat)
        VerbosePrint "In fft_pgi_cuda_test Released memory on the device PGI  ",istat
        if(istat /= 0)call allocchk("In fft_cuda_test dealloc device coutput_d ",istat)
        !   endif



        ! call self%set_ft(.true.)
        ! call self%set_cmat(coutput(:ldim(1), :ldim(2), :ldim(3)))
        ! ! omp parallel do collapse(3) schedule(static) default(shared) private(h,k,l,phys) proc_bind(close)
        ! do h=lims(1,1),lims(1,2)
        !     do l=lims(2,1),lims(2,2)
        !         do k=lims(3,1),lims(3,2)
        !             phys = self%comp_addr_phys([h,l,k])
        !             ! if (h > 0) then
        !             !     phys(1) = h + 1
        !             !     phys(2) = l + 1 + MERGE(ldim(2),0,l  < 0)
        !             !     phys(3) = k + 1 + MERGE(ldim(3),0,k  < 0)
        !             ! else
        !             !     phys(1) = -h + 1
        !             !     phys(2) = -l + 1 + MERGE(ldim(2),0,-l  < 0)
        !             !     phys(3) = -k + 1 + MERGE(ldim(3),0,-k  < 0)
        !             ! end if
        !             call self%set_cmat_at_ind(phys(1),phys(2),phys(3),coutput(h-lims(1,1)+1,l-lims(2,1)+1,k-lims(3,1)+1))
        !         end do
        !     end do
        ! end do
        ! ! omp end parallel do
        !   call gnufor_image(real(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(aimag(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(cabs(coutput(:,:,1)),text='gpo')
        !   call gnufor_image(atan2(real(coutput(:,:,1)),aimag(coutput(:,:,1))),text='gpo')

        VerbosePrint "In simple_cufft::fft_pgi  deallocating host arrays"
        ! Remove host storage
        ! if(allocated(rinput))then
        !     VerbosePrint "In simple_cufft::fft_pgi  deallocating rinput"
        !     deallocate(rinput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::fft dealloc host rinput",istat)
        ! end if
        ! if(allocated(rinput).and.allocated(coutput))then
        !     VerbosePrint "In simple_cufft::fft_pgi  deallocating coutput and rinput"
        !     deallocate(rinput, coutput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::fft dealloc host arrays",istat)
        ! end if
        ! if(allocated(coutput))then

        if(allocated(coutput))deallocate(coutput,stat=istat)
        VerbosePrint "In fft_pgi_cuda_test  deallocated coutput", istat
        if(istat /= 0) call allocchk("In fft_cuda_test dealloc host coutput",istat)
        ! end if
        ! if(allocated(rinput))then
        if(allocated(rinput))deallocate(rinput,stat=istat)
        VerbosePrint "In fft_pgi_cuda_test  deallocated rinput", istat
        if(istat /= 0)call allocchk("In fft_cuda_test dealloc host rinput",istat)
        ! end if
        VerbosePrint "In fft_pgi_cuda_test deallocating done "

        ! Destroy the plans
        call cufftDestroy(plan)
        !if (ierr /= 0) print *, " fft_pgi_cuda_test error: ", ierr
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In fft_cuda_test end ",__FILENAME__,__LINE__)
        endif
        VerbosePrint "In fft_pgi_cuda_test finished "
        !return
    end subroutine fft_pgi_cuda_test

     subroutine fft_pgi_stream_test (self)

        use simple_cufft
        use cudafor
        implicit none
        type(image), intent(inout)  :: self
        complex(fp_kind), allocatable, dimension(:,:,:),pinned :: A,B,C
        complex(fp_kind), allocatable, dimension(:,:,:),device :: A_d,B_d
        integer, parameter :: num_streams=4
        integer:: nx, ny, nomega, ifr, i,j, istat,stream_index
        integer(8):: clock_start,clock_end,clock_rate
        integer(kind=cuda_stream_kind) :: stream(num_streams)
        integer :: plan
        real:: elapsed_time
        real(fp_kind):: scale

        nx=512; ny=512;  nomega=196
        scale=1./real(nx*ny,fp_kind)

        ! Initialize FFT plan
        call cufftPlan2D(plan,nx,ny,CUFFT_C2C)

        ! Create streams
        do i = 1,num_streams
            istat= cudaStreamCreate(stream(i))
        end do

        call SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate

        ! Allocate arrays on CPU and GPU
        allocate(A(nx,ny,nomega),B(nx,ny,nomega),C(nx,ny,nomega))
        allocate(A_d(nx,ny,num_streams),B_d(nx,ny,num_streams))

        ! Initialize arrays on CPU
        A=cmplx(1.,1.,fp_kind); B=cmplx(1.,1.,fp_kind); C=cmplx(0.,0.,fp_kind)

        ! Measure only the transfer time
        istat=cudaThreadSynchronize()

        print *,"I/O only"
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

        do ifr=1,nomega
            istat= cudaMemcpy(A_d(1,1,1),A(1,1,ifr),nx*ny)
            istat= cudaMemcpy(B_d(1,1,1),B(1,1,ifr),nx*ny)
            istat= cudaMemcpy(C(1,1,ifr),A_d(1,1,1),nx*ny)
        end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! End timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        ! Measure the transfer time H2D, FFT , IFFT and transfer time D2H

        print '(/a)',"Single stream  loop"
        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
        stream_index = 1
    !     call cufftSetStream(plan,stream(stream_index))
    !     do ifr=1,nomega
    !         istat= cudaMemcpy(A_d(1,1,stream_index),A(1,1,ifr),nx*ny)
    !         istat= cudaMemcpy(B_d(1,1,stream_index),B(1,1,ifr),nx*ny)
    !         call cufftExec(plan,CUFFT_C2C ,A_d(1,1,stream_index),&
    !             A_d(1,1,stream_index),CUFFT_FORWARD)
    !         call cufftExec(plan, CUFFT_C2C,B_d(1,1,stream_index),&
    !             B_d(1,1,stream_index),CUFFT_FORWARD)

    !         ! Convolution and scaling of the  arrays
    !         !$cuf kernel do(2) <<<*,(16,16),stream=stream(stream_index)>>>
    !         do j=1,ny
    !             do i=1,nx
    !                 B_d(i,j,stream_index)= A_d(i,j,stream_index)*&
    !                     B_d(i,j,stream_index)*scale
    !             end do
    !         end do

    !         call cufftExec(plan, CUFFT_C2C,B_d(1,1,stream_index),&
    !             B_d(1,1,stream_index),CUFFT_INVERSE)
    !         istat=cudaMemcpy( C(1,1,ifr),B_d(1,1,stream_index),nx*ny)
    !     end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! End timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        ! Overlap I/O and compute using multiple streams and async copies
        print '(/a)',"Do loop with multiple streams"
        call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

        do ifr=1,nomega

            ! assign a stream for the current plan
            stream_index = mod(ifr,num_streams)+1

    !         ! Set the stream used by CUFFT
    !         call cufftSetStream(plan,stream(stream_index))

    !         ! Send A to GPU
    !         istat= cudaMemcpyAsync(A_d(1,1,stream_index),A(1,1,ifr),&
    !             nx*ny, stream(stream_index))

    !         ! Execute forward FFTs on GPU
    !         call cufftExec(plan ,CUFFT_C2C,A_d(1,1,stream_index),&
    !             A_d(1,1,stream_index),CUFFT_FORWARD)

    !         ! Send B to GPU
    !         istat= cudaMemcpyAsync(B_d(1,1,stream_index), &
    !             B(1,1,ifr),nx*ny, stream(stream_index))

    !         ! Execute forward FFTs on GPU
    !         call cufftExec(plan, CUFFT_C2C,B_d(1,1,stream_index),&
    !             B_d(1,1,stream_index),CUFFT_FORWARD)

    !         ! Convolution and scaling of the  arrays
    !         !$cuf kernel do(2) <<<*,(16,16),stream=stream(stream_index)>>>
    !         do j=1,ny
    !             do i=1,nx
    !                 B_d(i,j,stream_index)= A_d(i,j,stream_index)* &
    !                     B_d(i,j,stream_index)*scale
    !             end do
    !         end do

    !         ! Execute inverse FFTs on GPU
    !         call cufftExec(plan, CUFFT_C2C ,B_d(1,1,stream_index), &
    !             B_d(1,1,stream_index),CUFFT_INVERSE)

    !         ! Copy results back
    !         istat=cudaMemcpyAsync( C(1,1,ifr),B_d(1,1,stream_index), &
    !             nx*ny, stream=stream(stream_index))

        end do

        istat=cudaThreadSynchronize()
        call SYSTEM_CLOCK(COUNT=clock_end) ! Start timing
        elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)
        print *,"Elapsed time :",elapsed_time

        if (istat .eq. 0) then
            print *,"Test Passed"
        else
            print *,"Test Failed"
        endif

         deallocate(A,B,C); deallocate(A_d,B_d)
         call cufftDestroy(plan)

     end subroutine fft_pgi_stream_test

    subroutine ifft_pgi_cuda_test( self )
        !#ifndef PGI
        !        class(image), intent(inout) :: self
        !        print *,'CUDA ifft cannot be performed with current compiler'
        !#else
        !! IN-PLACE 3D FFT using CUDA's cufft
        use simple_cufft
        use simple_fftshifter
        use openacc
        use cudafor

        type(image), intent(inout)       :: self
        real(dp),    allocatable         :: routput(:,:,:),rtmp(:,:,:)
        real(dp),    allocatable, device :: routput_d(:,:,:)
        complex(dp), allocatable, device :: cinput_d(:,:,:)
        complex(dp), allocatable         :: cinput(:,:,:)
        complex(sp), allocatable         :: ctmp(:,:,:)
        ! type(image)  :: c1,c2,c3
        integer      :: plan, planType
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(dp)      :: comp, cswap
        integer(timer_int_kind) :: t1

        verbose=.false.
        !   if(verbose) t1=tic()
        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        fsize = self%get_filtsz()
        lims = self%loop_lims(2)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI cmat lims", lims
        nscale =  real(product(ldim))
        VerbosePrint "In simple_cufft::ifft_pgi  PGI LDIM/CDIM", ldim, cdim

        ! allocate arrays on the host
        if(allocated(routput))deallocate(routput)
        allocate(routput(ldim(1), ldim(2), ldim(3)),source=0._dp,stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft routput",istat)
        if(allocated(cinput))deallocate( cinput )
        allocate(cinput((ldim(1)/2)+1, ldim(2), ldim(3)),source=dcmplx(0._dp,0._dp),stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft cinput",istat)
        allocate(ctmp(cdim(1),cdim(2),cdim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat.ne.0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !rtmp = self%get_rmat()
        ctmp = self%get_cmatfull()
        !! VerbosePrint "In simple_cufft::ifft_pgi  PGI routput ", shape(routput), LBOUND(routput), UBOUND(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI cinput ", shape(cinput), LBOUND(cinput), UBOUND(cinput)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI ctmp ", shape(ctmp), LBOUND(ctmp), UBOUND(ctmp)
        !call ifftshift(ctmp)
        !call gnufor_image(rtmp(:,:,1), palette='gray')
        !call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        !call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')

        cinput=dcmplx(ctmp(:(ldim(1)/2)+1,:ldim(2),:ldim(3)))
        print *,"shape cinput", shape(cinput)
        !cinput(:cdim(1),:ldim(2),:ldim(3)) = cmplx(rtmp(1:ldim(1):2,:ldim(2),:ldim(3)),rtmp(2:ldim(1):2,:ldim(2),:ldim(3)))
        ! cinput(cdim(1)+2:ldim(1),1:ldim(2)-1,:ldim(3)) =  cmplx(rtmp(:cdim(1)-2,2:ldim(2),:ldim(3))
        !  cinput(cdim(1)+1:ldim(1),:cdim(2),:cdim(3)) = ctmp(cdim(1):-1:2, :cdim(2), :cdim(3))
        !conjg(ctmp(2:cdim(1), :cdim(2),:cdim(3)))
        ! !omp parallel do collapse(3) default(shared) private(cswap,i,j,k)&
        ! !omp schedule(static) proc_bind(close)
        !         do i=1,ldim(1)/2
        !             do j=1,ldim(2)/2
        !                 do k=1,ldim(3)/2
        !                     !(1)
        !                     cswap = cinput(i,j,k)
        !                     cinput(i,j,k) = cinput(ldim(1)/2+i,ldim(2)/2+j,ldim(3)/2+k)
        !                     cinput(ldim(1)/2+i,ldim(2)/2+j,ldim(3)/2+k) = cswap
        !                     !(2)
        !                     cswap = cinput(i,ldim(2)/2+j,ldim(3)/2+k)
        !                     cinput(i,ldim(2)/2+j,ldim(3)/2+k) = cinput(ldim(1)/2+i,j,k)
        !                     cinput(ldim(1)/2+i,j,k) = cswap
        !                     !(3)
        !                     cswap = cinput(ldim(1)/2+i,j,ldim(3)/2+k)
        !                     cinput(ldim(1)/2+i,j,ldim(3)/2+k) = cinput(i,ldim(2)/2+j,k)
        !                     cinput(i,ldim(2)/2+j,k) = cswap
        !                     !(4)
        !                     cswap = cinput(i,j,ldim(3)/2+k)
        !                     cinput(i,j,ldim(3)/2+k) = cinput(ldim(1)/2+i,ldim(2)/2+j,k)
        !                     cinput(ldim(1)/2+i,ldim(2)/2+j,k) = cswap
        !                 end do
        !             end do
        !         end do
        !         !omp end parallel do

        ! call ifftshift(cinput)
        !     call gnufor_image((real(cinput(:,:,1),sp)),palette='gray')
        !      call gnufor_image(real(aimag(cinput(:,:,1)),sp),palette='gray')

        ! !omp parallel do collapse(3) schedule(static) default(shared) private(h,k,l,phys,comp) proc_bind(close)
        ! do h=lims(1,1),lims(1,2)
        !     do k=lims(2,1),lims(2,2)
        !         do l=lims(3,1),lims(3,2)
        !             !                        phys = self%comp_addr_phys([h,k,l])
        !             if (h > 0) then
        !                 phys(1) = h + 1
        !                 phys(2) = k + 1 + MERGE(ldim(2),0,k  < 0)
        !                 phys(3) = l + 1 + MERGE(ldim(3),0,l  < 0)
        !             else
        !                 phys(1) = -h + 1
        !                 phys(2) = -k + 1 + MERGE(ldim(2),0,-k  < 0)
        !                 phys(3) = -l + 1 + MERGE(ldim(3),0,-l  < 0)
        !             end if
        !             comp = self%get_cmat_at_ind(phys(1),phys(2),phys(3))
        !             ! if(h<0)comp=conjg(comp)
        !             cinput(phys(1),phys(2),phys(3))=comp
        !         end do
        !     end do
        ! end do
        ! !omp end parallel do
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        endif
        allocate(routput_d(ldim(1), ldim(2), ldim(3)),stat=istat)
        if(istat.ne.0) call allocchk("In simple_image::ifft routput_d",istat)

        allocate(cinput_d((ldim(1)/2)+1, ldim(2), ldim(3)),stat=istat)
        if(istat.ne.0) call allocchk("In simple_image::ifft cinput_d",istat)
        ! set planType to either single or double precision
        planType = CUFFT_Z2D    !! (fp_kind == singlePrecision)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI device allocated ", size(cinput_d), size(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput), LBOUND(cinput), UBOUND(cinput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput_d), LBOUND(cinput_d), UBOUND(cinput_d)
        ! copy input to device
        cinput_d = cinput
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI "
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI host -> device ", shape(cinput_d), lbound(cinput_d), ubound(cinput_d)


        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        VerbosePrint "In simple_cufft::ifft_pgi  PGI "

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,cinput_d,routput_d,CUFFT_INVERSE )
        VerbosePrint "In simple_cufft::ifft_pgi  PGI exec completed "

        ! ! Copy results back to host
        ! istat=cudaThreadSynchronize()
        ! istat=istat+cudaDeviceSynchronize()
        ! if(istat.ne.0)then
        !     call simple_cuda_stop("In simple_image::ifft post fft sync ",__FILENAME__,__LINE__)
        ! endif
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI copy results back to host"
        ! Copy results back to host
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI device -> host ", size(cinput_d), size(routput)
        !VerbosePrint "In simple_cufft::ifft_pgi  PGI host shape ", shape(routput), LBOUND(routput), UBOUND(routput)
        routput = routput_d
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI device -> host ", size(cinput_d), size(routput)
        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI deallocate arrays on device"

        ! release memory on the host and on the device
        ! if(allocated(routput_d))then
        deallocate(routput_d)!,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc routput_d",istat)
        !  end if
        !  if(allocated(cinput_d))then
        deallocate(cinput_d)!,stat=istat)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc cinput_d",istat)
        !  end if

        ! VerbosePrint "In simple_cufft::ifft_pgi  PGI destroy cufft plans"
        ! Destroy the plans
        call cufftDestroy(plan)
        istat=cudaThreadSynchronize()
        istat=istat+cudaDeviceSynchronize()
        if(istat /= 0) then
            call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        endif
        VerbosePrint "In simple_cufft::ifft_pgi  PGIC UFFT vs FFTW3"
        !call fftshift(routput)
        call gnufor_image((real(routput(:,:,1))),text='go',palette='hot')
        call self%set_rmat( real(routput(:ldim(1),:ldim(2),:ldim(3))) )
        call self%set_ft(.false.)
        call self%shift_phorig
        ! Cleanup
        if(allocated(rtmp))then
            deallocate( rtmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating rtmp",istat)
        end if
        if(allocated(ctmp))then
            deallocate( ctmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating ctmp",istat)
        end if
        if(allocated(routput))then
            ! deallocate(routput,stat=istat)
            ! if(istat /= 0)call allocchk("In simple_image::ifft deallocatiing routput",istat)
        end if
        if(allocated(cinput))then
            deallocate(cinput,stat=istat)
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinput",istat)
        end if

        !  #endif
        return
    end subroutine ifft_pgi_cuda_test


    subroutine ifft_pgi_cuda_test2( self )

        use simple_cufft
        use simple_fftshifter
        use openacc
        use cudafor

        type(image), intent(inout)       :: self
        ! real(sp),    allocatable         :: routput(:,:,:)
        !complex(sp),    allocatable, device :: routput_d(:,:,:)
        !  complex(sp), allocatable, device :: cinput_d(:,:,:)
        complex(sp), allocatable, device :: inplace_d(:,:,:)
        ! complex(sp), allocatable         :: cinput(:,:,:), coutput(:,:,:)
        complex(sp), allocatable         :: cinout(:,:,:),ctmp(:,:,:)
        ! type(image)  :: c1,c2,c3
        integer      :: plan, planType, xdim,ydim
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(sp)      :: comp, cswap
        integer(timer_int_kind) :: t1
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        verbose=.true.
        if(verbose)then
            t1=tic()
            print *,"In ifft_pgi_cuda_test2 "
        end if
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ldim = self%get_ldim()
        cdim = self%get_array_shape()
        cdim(1)=cdim(1)
        fsize = self%get_filtsz()
        lims = self%loop_lims(3)
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cmat lims, fsize", lims, fsize
        nscale =  real(product(ldim))
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI LDIM/nscale", ldim, nscale
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI CDIM       ", cdim
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ! allocate arrays on the host
        ! if(allocated(coutput))deallocate(coutput)
        ! allocate(coutput(ldim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        ! if(istat.ne.0)call allocchk("In simple_image::ifft routput",istat)
        ! if(allocated(cinout))deallocate( cinout )
        ! print *,"Half plus 1", ceiling(real(ldim(1))/2.)+1
        !allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)
        ! if(istat.ne.0)call allocchk("In simple_image::ifft cinout",istat)
        ! allocate(ctmp(lims(1,1):lims(1,2),lims(2,1):lims(2,2),lims(3,1):lims(3,2)),source=cmplx(0.,0.),stat=istat)
        allocate(ctmp(cdim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        if(istat /= 0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !  rtmp = self%get_rmat()
        ctmp = self%get_cmatfull()
        !call ifftshift(ctmp)
        call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp),lbound(ctmp), ubound(ctmp)
        k=1

        if(is_even(ldim(1)))then
            xdim = ldim(1)/2
            ydim = ldim(2)/2
        else
            xdim = (ldim(1)-1)/2
            ydim = (ldim(2)-1)/2
        endif
        !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)
        allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)
        cinout = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,phys)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                if(h<0)then
                    cinout(h+xdim+1,k+ydim+1,1) = ctmp(phys(1),phys(2),1)
                else
                    cinout(h+xdim+1,k+ydim+1,1) = conjg(ctmp(phys(1),phys(2),1))
                endif

            end do
        end do
        !$omp end parallel do
        ! call fftshift(cinout)
        call gnufor_image(real(cinout(:,:,1)), palette='gray')
        call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !    read(*,*)
        !         deallocate( cinout )
        !         allocate(cinout(ldim(1), ldim(2), ldim(3)),stat=istat)

        !         !$omp parallel do collapse(2) default(shared) private(i,j)
        !   do i=1,cdim(1)
        !          do j=1,floor(real(ldim(2))/2.)
        !              k=1
        ! !             do k=1,floor(real(ldim(3))/2.)+1
        !              cinout(i,j,k) =( ctmp(i,j,k))

        !                !if(i/=cdim(1))cinout(ldim(1)-i+1,j,k) =  conjg(ctmp(i,ldim(2)-j+1,k))

        !                !  cinout(i,ldim(2)-j+1,k) = ctmp(i,ldim(2)-j+1,k)
        !             !     if(i/=1)cinout(ldim(1)-i+2,ldim(2)-j+1,k) =  ( ctmp(i,j,k)) !ldim(2)-j+1,k) )
        !                  ! if(i/=cdim(1))cinout(ldim(1)-i+1,ldim(2)-j+1,k) =  ( ctmp(i,ldim(2)-j+1,k) )
        !         !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !         !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !  !                end do
        !              end do
        !              do j=floor(real(ldim(2))/2.)+1,ldim(2)
        !              k=1
        ! !             do k=1,floor(real(ldim(3))/2.)+1
        !             ! cinout(i,j,k) = ( ctmp(i,j,k))

        !                !if(i/=cdim(1))cinout(ldim(1)-i+1,j,k) =  conjg(ctmp(i,ldim(2)-j+1,k))

        !                !  cinout(i,ldim(2)-j+1,k) = ctmp(i,ldim(2)-j+1,k)
        !              !    if(i/=1)cinout(ldim(1)-i+2,ldim(2)-j+1,k) =  ( ctmp(i,j,k)) !ldim(2)-j+1,k) )

        !         !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !         !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !  !                end do
        !                     end do
        !          end do
        !          !$omp end parallel do
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), 1:ceiling(cdim(2)/2), 1:ceiling(cdim(3)/2)
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), ceiling(cdim(2)/2)+1:cdim(2), ceiling(cdim(3)/2)+1:cdim(2)
        !         ! cinout(cdim(1)-1:ldim(1), cdim(2)-1:ldim(2), cdim(3)-1:ldim(3))=ctmp(1:cdim(1), 1:ceiling(cdim(2)/2), 1:ceiling(cdim(3)/2)


        !         VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp), lbound(ctmp), ubound(ctmp)
        !         !call ifftshift(cinout)

        !         !call gnufor_image(cabs(ctmp(:cdim(1),:ldim(2),1)), text='go',palette='hot')
        !         call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !         call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !         VerbosePrint "In simple_cufft::ifft_pgi2  ctmp shape", shape(ctmp)
        !        ! cinout(:cdim(1), :ldim(2), :ldim(3)) = cmplx(ctmp(:cdim(1), :ldim(2), :ldim(3)))
        !        ! VerbosePrint "In simple_cufft::ifft_pgi2  cinout shape", shape(cinout)

        !         call fftshift(cinout)
        !         call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !         call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        ! Synchronize
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"

        !allocate on device
        ! if(allocated(routput_d)) deallocate( routput_d )
        ! VerbosePrint "In simple_cufft::ifft_pgi2 device allocation"
        ! call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        allocate(inplace_d(ldim(1), ldim(2), ldim(3)))!,stat=istat)
        call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)
        !allocate(routput_d(cdim(1), ldim(2), ldim(3)))!,stat=istat)
        !call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)
        !if(istat.ne.0) call allocchk("In simple_image::ifft routput_d",istat)
        VerbosePrint "In simple_cufft::ifft_pgi2 routput_d allocation completed"
        ! if(allocated(cinout_d)) deallocate( cinout_d )
        ! allocate(cinout_d(cdim(1), cdim(2), cdim(3)))!,stat=istat)
        ! call simple_cuda_stop("In ifft TEST2 cinout_d",__FILENAME__,__LINE__)
        !if(istat.ne.0) call allocchk("In simple_image::ifft cinout_d",istat)
        VerbosePrint "In simple_cufft::ifft_pgi2 device allocation completed"

        ! set planType to either single or double precision
        planType = CUFFT_C2C
        VerbosePrint "In simple_cufft::ifft_pgi2 set planType completed"

        ! copy input to device
        ! routput_d = cinout
        inplace_d = cinout
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 copy input to device completed"

        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Initialize the plan completed"

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_INVERSE )
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Execute  Backward transform completed"

        ! Copy results back to host
        cinout = inplace_d

        VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp), ubound(ctmp),lbound(ctmp)
        ! !$omp parallel do default(shared) private(i,j,k)
        ! !$acc kernels loop
        ! do j=1,cdim(2)
        !     do k=1,cdim(3)
        !         !$acc loop gang(16), vector(32)
        !         do i=1,cdim(1)-1
        !             routput((i-1)*2+1,j,k) = real(cinout(i,j,k))
        !             routput((i-1)*2+2,j,k) = aimag(cinout(i,j,k))
        !         end do
        !         !$acc end loop
        !     end do
        ! end do
        ! !$acc end kernels
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(i)
        !  do i=1,cdim(1)-1
        !      routput((i-1)*2+1, :ldim(2), :ldim(3)) = real(coutput(:cdim(1),:cdim(2),:cdim(3)))
        !      routput((i-1)*2+2, :ldim(2), :ldim(3)) = aimag(coutput(:cdim(1),:cdim(2),:cdim(3)))
        ! end do
        ! !$omp end parallel do
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Copy results back to host completed"

        ! release memory on the host and on the device
        !  if(allocated(routput_d))then
        !      call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        deallocate(inplace_d)!,stat=istat)
        call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !if(istat /= 0)call allocchk("In simple_image::ifft dealloc routput_d",istat)
        !  end if
        ! if(allocated(cinout_d))then
        !     call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !     deallocate(cinout_d)!,stat=istat)
        !     call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        !     !if(istat /= 0)call allocchk("In simple_image::ifft dealloc cinout_d",istat)
        ! end if
        !        call fftshift(ctmp)
        call gnufor_image(real(cinout(:,:,1)), palette='gray')
        call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        call gnufor_image(cabs(cinout(:,:,1)), palette='gray')
        !call gnufor_image(aimag(cinout(:,:,1)), palette='gray')
        !call fftshift(routput)
        !        call self%set_cmat( cinout )
        call self%set_ft(.false.)
        call self%set_rmat( real(cinout) )
        call self%shift_phorig
        ! Cleanup
        ! if(allocated(rtmp))then
        !     deallocate( rtmp, stat=istat )
        !     if(istat /= 0)call allocchk("In simple_image::ifft deallocating rtmp",istat)
        ! end if
        if(allocated(ctmp))then
            deallocate( ctmp, stat=istat )
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating ctmp",istat)
        end if
        ! if(allocated(coutput))then
        !     deallocate(coutput,stat=istat)
        !     if(istat /= 0)call allocchk("In simple_image::ifft deallocatiing coutput",istat)
        ! end if
        if(allocated(cinout))then
            deallocate(cinout,stat=istat)
            if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinout",istat)
        end if
        ! Destroy the plans
        call cufftDestroy(plan)
        !if (ierr /= 0) print *, " fft_pgi_cuda_test error: ", ierr
        istat=cudaDeviceSynchronize()
        istat=istat+cudaThreadSynchronize()
        if(istat.ne.0)then
            call simple_cuda_stop("In simple_image::ifft end ",__FILENAME__,__LINE__)
        endif
        call self%set_ft(.false.)
        VerbosePrint "In simple_cufft::ifft_pgi finished "
    end subroutine ifft_pgi_cuda_test2

    subroutine ifft_pgi_cuda_test3( self )

        use simple_cufft
        !  use simple_fftshifter
        !  use openacc
        !  use cudafor

        type(image), intent(inout)       :: self

        complex(sp), allocatable, device :: inplace_d(:,:,:)

        complex(sp), allocatable         :: cinout(:,:,:) !,ctmp(:,:,:)
        !   complex(sp), allocatable         :: coutput(:,:,:)
        integer      :: plan, planType, xdim,ydim
        integer      :: i,j,k,h,l,n,istat,mk,mh,ml,fsize
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        complex(sp)      :: comp, cswap
        integer(timer_int_kind) :: t1
        !   call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        verbose=.true.
        ! if(verbose)then
        !     t1=tic()
        !     print *,"In ifft_pgi_cuda_test2 "
        ! end if
        ! call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)
        call self%set_ft(.false.) ;print *,"In simple_image::ifft set ft"
        ldim = self%get_ldim()
        allocate(cinout(ldim(1),ldim(2),ldim(3)),inplace_d(ldim(1),ldim(2),ldim(3)))
        ! cdim = self%get_array_shape()
        ! cdim(1)=cdim(1)
        ! fsize = self%get_filtsz()
        ! lims = self%loop_lims(3)
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI cmat lims, fsize", lims, fsize
        ! nscale =  real(product(ldim))
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI LDIM/nscale", ldim, nscale
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI CDIM       ", cdim
        !  call simple_cuda_stop("In ifft TEST2 ",__FILENAME__,__LINE__)

        ! allocate(ctmp(cdim(1), ldim(2), ldim(3)),source=cmplx(0.,0.),stat=istat)
        ! if(istat /= 0)call allocchk("In simple_image::ifft ctmp",istat)

        ! allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        !  rtmp = self%get_rmat()
        ! ctmp = self%get_cmatfull()
        !call ifftshift(ctmp)
        !   call gnufor_image(real(ctmp(:,:,1)), palette='gray')
        !  call gnufor_image(aimag(ctmp(:,:,1)), palette='gray')
        ! VerbosePrint "In simple_cufft::ifft_pgi2  PGI ctmp ",size(ctmp),lbound(ctmp), ubound(ctmp)
        k=1

        if(is_even(ldim(1)))then
            xdim = ldim(1)/2
            ydim = ldim(2)/2
        else
            xdim = (ldim(1)-1)/2
            ydim = (ldim(2)-1)/2
        endif
        !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)

        cinout = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,phys,comp)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                comp= self%get_cmat_at_ind(phys(1),phys(2),1)
                if(h>=0)then
                    cinout(h+xdim+1,k+ydim+1,1) = comp
                else
                    cinout(h+xdim+1,k+ydim+1,1) = conjg(comp)
                endif

            end do
        end do
        !$omp end parallel do
        ! call fftshift(cinout)
        !call gnufor_image(real(cinout(:,:,1)), palette='gray')
        !call gnufor_image(aimag(cinout(:,:,1)), palette='gray')

        ! Synchronize
        !  istat=cudaDeviceSynchronize()
        !  istat=istat+cudaThreadSynchronize()
        !  call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"


        ! allocate(inplace_d(ldim(1), ldim(2), ldim(3)),stat=istat)
        ! call simple_cuda_stop("In ifft TEST2 routput_d",__FILENAME__,__LINE__)



        ! set planType to either single or double precision
        planType = CUFFT_C2C
        !VerbosePrint "In simple_cufft::ifft_pgi2 set planType completed"

        ! copy input to device
        ! routput_d = cinout
        inplace_d = cinout
        VerbosePrint "In simple_cufft::ifft_pgi2 copy input to device completed"

        ! Initialize the plan for complex to real transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        VerbosePrint "In simple_cufft::ifft_pgi2 Initialize the plan completed"

        ! Execute  Backward transform in place
        call cufftExec(plan,planType,inplace_d,inplace_d,CUFFT_INVERSE )
        VerbosePrint "In simple_cufft::ifft_pgi2 Execute  Backward transform completed"

        ! Copy results back to host
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI inplace_d ",size(inplace_d), ubound(inplace_d),lbound(inplace_d)
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cinout ",size(cinout), ubound(cinout),lbound(cinout)
        cinout = inplace_d
        VerbosePrint "In simple_cufft::ifft_pgi2 Copied results back to host completed"
        VerbosePrint "In simple_cufft::ifft_pgi2  PGI cinout ",size(cinout), ubound(cinout),lbound(cinout)
        ! Save to image obj

        !call self%set_rmat( cabs(cinout) ) ;print *,"In simple_image::ifft cabs inout"
        !$omp parallel do collapse(2) default(shared) private(i,j)
        do i=1,ldim(1)
            do j=1,ldim(2)
                call self%set_rmat_at_ind(i, j, 1, cabs( cinout(i,j,1)) )
            end do
        end do
        !$omp end parallel do
        !  call self%shift_phorig ;print *,"In simple_image::ifft shift"
        ! Destroy the plans
        call cufftDestroy(plan)

        ! release memory on the host and on the device
        print *,"In simple_image::ifft deallocating cinout"

        ! deallocate(cinout ) ; print *," deallocated cinout"
        !   if(istat /= 0)call allocchk("In simple_image::ifft deallocating cinout",istat)
        !   deallocate(ctmp); print *," deallocated ctmp"
        deallocate(cinout, inplace_d ) ; print *," deallocated cinout"
        !deallocate( inplace_d) ; print *," deallocated inplace_d"

        ! Synchronize
        !  istat=cudaDeviceSynchronize()
        !  istat=istat+cudaThreadSynchronize()
        !   call simple_cuda_stop("In simple_image::ifft sync ",__FILENAME__,__LINE__)
        VerbosePrint "In simple_cufft::ifft_pgi2 Synchronize completed"

        VerbosePrint "In simple_cufft::ifft_pgi finished "
    end subroutine ifft_pgi_cuda_test3

    subroutine test_pgi( ld1, ld2, ld3, doplot,iter )
 !       use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        integer, intent(in),optional:: iter
        type(image)          :: img, img_2
        integer              :: i, j, k, iter_here,cnt, lfny, ldim(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        real                 :: smpd=2.
        logical              :: passed, test(6)

        iter_here=1
        if(present(iter))iter_here=iter
        verbose=.false.

        write(*,'(a)') '**info(unit_test, PGI CUFFT: test_pgi  '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)

        call img%gauimg(10)
        call img%gauimg2(5, floor(real(ld1)/5.),floor(real(ld2)/3.))
        call img%gauimg2(8, -floor(real(ld1)/5.),floor(real(ld2)/6.))
        call img%add_gauran(0.5)

        img_2 = img

        do i=1,iter_here
            call fft_pgi_cuda_test(img_2)
            call ifft_pgi_cuda_test4(img_2)
        end do

        write(*,'(a,ES20.10)') "**info(simple_image PGI test): L2 norm sum ", img.lonesum.img_2




        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        if(doplot)call img_2%vis(geomorsphr=.false.)
        !if(doplot)call img%vis(geomorsphr=.false.)
        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)

        !     write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '
        !     call img_3%fwd_ft()
        !     call img_3%ft2img('amp',img_4)
        !     if(doplot)call img_4%vis(geomorsphr=.false.)
        !     call img_3%bwd_ft()
        ! !    write(*,'(a)') '**info(simple_image PGI test): testing subtr '
        !   !  read *
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !     print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2

        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '

        !   call img%add(img_2)
        ! read(*,*)
        !  write(*,'(a)') '**info(simple_image PGI test): comparing cufft and fftw3 '



        ! do i=1,ldim(1)
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs(cabs(img%get_cmat_at_ind(i,j,k)) - cabs(img_2%get_cmat_at_ind(i,j,k))) >= TINY )then
        !                 print *,'wrong FFT i/j/k fftw cufft', i,j,k,real(img%get_cmat_at_ind(i,j,k)) ,&
        !                     &aimag(img%get_cmat_at_ind(i,j,k)) , real( img_2%get_cmat_at_ind(i,j,k)), &
        !                     & aimag( img_2%%get_cmat_at_ind(i,j,k))
        !                 call simple_stop("image FFT inverse PGI test  stopped")
        !             end if
        !         end do
        !     end do
        ! end do
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L1 norm sum", img.lonesum.img_2

        ! do i=1,ldim(1)
        !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
        !     do j=1,ldim(2),2
        !         do k=1,ldim(3)
        !             if (abs( img%%get_rmat_at_ind(i,j,k) - img_2%get_rmat_at_ind(i,j,k)) >= 1.e-06 )then
        !                 print *,'wrong FFT i/j/k fftw /cufft', i,j,k,img%%get_rmat_at_ind(i,j,k) , img_2%%get_rmat_at_ind(i,j,k)
        !                 call simple_stop("image FFT PGI test stopped")
        !             end if
        !         end do

        !     end do
        ! end do
        call img_2%kill
        call img%kill


    end subroutine test_pgi

    !> test_pgi1 is a simple fwd/bwd FFT
    subroutine test_pgi1( ld1, ld2, ld3, doplot ,iter)
!        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        integer, intent(in),optional:: iter
        type(image)          :: img, img_2
        integer              :: i, j, k, iter_here,cnt, lfny, ldim(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        real                 :: smpd=2.
        logical              :: passed, test(6)
        integer(timer_int_kind) :: t1
        iter_here=1
        if(present(iter))iter_here=iter

        verbose=.false.
        write(*,'(a)') '**info(unit_test, PGI CUFFT: test_pgi  '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)

        call img%gauimg(10)
        call img%gauimg2(5, floor(real(ld1)/5.),floor(real(ld2)/3.))
        call img%gauimg2(8, -floor(real(ld1)/5.),floor(real(ld2)/6.))
        call img%add_gauran(0.5)
        if(doplot)call img%vis(geomorsphr=.false.)
        img_2 = img

        do i=1,iter_here
            call img_2%fft_pgi_cuda()
            call img_2%ifft_pgi_cuda()
        end do
        write(*,'(a,ES20.10)') "**info(simple_image PGI test): L1 norm sum ", img.lonesum.img_2
        if(doplot)call img_2%vis()


        !  write(*,'(a)') '**info(simple_image PGI test): comparing cufft and fftw3 '



        ! do i=1,ldim(1)
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs(cabs(img%get_cmat_at_ind(i,j,k)) - cabs(img_2%get_cmat_at_ind(i,j,k))) >= TINY )then
        !                 print *,'wrong FFT i/j/k fftw cufft', i,j,k,real(img%get_cmat_at_ind(i,j,k)) ,&
        !                     &aimag(img%get_cmat_at_ind(i,j,k)) , real( img_2%get_cmat_at_ind(i,j,k)), &
        !                     & aimag( img_2%%get_cmat_at_ind(i,j,k))
        !                 call simple_stop("image FFT inverse PGI test  stopped")
        !             end if
        !         end do
        !     end do
        ! end do
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L2 norm sum", img.lonesum.img_2

        ! do i=1,ldim(1)
        !     do j=1,ldim(2),2
        !         do k=1,ldim(3)
        !             if (abs( img%get_rmat_at_ind(i,j,k) - img_2%get_rmat_at_ind(i,j,k)) >= 1.e-05_sp )then
        !                 print *,'test_pgi1: mismatch FFT i/j/k fftw /cufft', i, j, k,&
        !                     img%get_rmat_at_ind(i,j,k) , img_2%get_rmat_at_ind(i,j,k)
        !               !  call simple_stop("image FFT PGI test stopped")
        !             end if
        !         end do

        !     end do
        ! end do
        call img_2%kill
        call img%kill

    end subroutine test_pgi1


    subroutine test_pgi2( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4!, img_5, img3d
        !    type(image)          :: imgs(20)
        complex, allocatable :: ctmp(:,:,:), ctmp2(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), lims(3,2), phys(3)
        real                 :: input, msk, ave, sdev, var, med, xyz(3), pow
        real                 :: imcorr, recorr, corr, corr_lp, maxv, minv
        !   real, allocatable    :: pcavec1(:), pcavec2(:), spec(:), res(:),rsh(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)


        verbose=.false.
        !debug=.true.
        !write(*,'(a)') '**info(simple_image PGI unit_test, part 1): testing basal constructors'

        ! call img%new([ld1,ld2,1], 1.)

        ! call img_3%new([ld1,ld2,1], 1.)

        ! call img3d%new([ld1,ld2,1], 1.)
        ! call img_4%new([ld1,ld2,1], 1.)

        ! call img_4%set_rmat_at_ind(5,5,1, 1.)
        ! call img_4%set_rmat_at_ind(10,3,1, 1.)
        ! call img_4%set_rmat_at_ind(27,26,1, 1.)
        ! call img_4%set_rmat_at_ind(29,25,1, 1.)


        ! if( .not. img%exists() ) call simple_stop('ERROR, in constructor or in exists function, 1')
        ! if( .not. img3d%exists() ) call simple_stop('ERROR, in constructor or in exists function, 2')


        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        !  call img_4%new(ldim, smpd)
        !     call img_5%new(ldim, smpd)
        call img%gauimg(10)
        call img%gauimg2(10, -floor(ld1/4.),-floor(ld2/3.))
        call img%gauimg2(10, floor(ld1/5.),-floor(ld2/8.))
        call img%gauimg2(10, floor(ld1/3.),floor(ld2/6.))
        call img%add_gauran(1.)
        !   call img%set_rmat_at_ind(floor(ld1/4.),floor(ld2/5.),1, 1.)
        !   call img%set_rmat_at_ind(floor(ld1/6.),floor(ld2/3.),1, 1.)
        if(doplot)call img%vis()
        img_2 = img
        img_3 = img
        !     if(doplot)call img_2%vis
        !  ave = (img.lonesum.img_2)
        !  print *," L1 norm sum", ave
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftshifter bwd'
        !     write(*,'(a)') 'Press enter to test ifftshift '
        !     read(*,*)
        !  call img_2%bwd_fftshift()
        !  img_4= img_2
        !  if(doplot)call img_4%vis(geomorsphr=.false.)


        !write(*,'(a)') 'Press enter to test cufft '
        !read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): testing cufft '

        call fft_pgi_cuda_test(img_2)
        call simple_cuda_stop("In PGI TEST2 ",__FILENAME__,__LINE__)
        write(*,'(a)') '**info(simple_image PGI test): testing fftw '
        call img_3%fwd_ft()
        print *,"FFTW"
        call img_3%print_cmat()
        !call img_2%ft2img('amp',img_4)
        print *,"CUFFT"
        call img_2%print_cmat()


        !     print *,"Diff"
        !     call img_4%print_cmat()

        if(doplot)call img_3%vis(geomorsphr=.false.)
        if(doplot)call img_2%vis(geomorsphr=.false.)

        !   img_4 = img_2
        !   img_5 = img_3
        print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2
        !  lims=img_3%loop_lims(2)
        !   write(*,'(a)') '**info(simple_image PGI test): cufft returned '
        !   write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        !     call ifft_pgi_cuda_test(img_2)
        ! read(*,*)
        ! if(doplot)call img_2%vis
        ! print *," L1 norm sum", img.lonesum.img_2

        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        !    if(doplot)call img_2%vis(geomorsphr=.false.)


        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)


        !        call img_3%ft2img('amp',img_4)
        ! ! if(doplot)call img_3%vis(geomorsphr=.false.)
        !   call img_3%bwd_ft()
        !   if(doplot)call img_3%vis()
        !   if(doplot)call img_2%vis()
        ! ! !    write(*,'(a)') '**info(simple_image PGI test): testing subtr '
        ! ! !  read *
        !   print *," L1 norm sum Orig v CUFFT", img.lonesum.img_2
        !  print *," L1 norm sum Orig v FFTW", img.lonesum.img_3
        !   print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2

        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '

        !   call img%add(img_2)
        !read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): comparing FWD cufft and fftw3 '
        passed=.true.
        if(allocated(ctmp)) deallocate(ctmp)
        if( allocated(ctmp2) ) deallocate(ctmp2)
        !$ allocate(ctmp(ldim(1)/2 +1,ldim(2),ldim(3)),ctmp2(ldim(1)/2 +1,ldim(2),ldim(3)) )
        ctmp2= img_2%get_cmatfull(); ctmp= img_3%get_cmatfull()
        VerbosePrint " Ldim/shapes ", ldim, shape(ctmp), shape(ctmp2)

        do i=1,ldim(1)/2 +1 !lims(1,1),lims(1,2)
            do j=1,ldim(2)!lims(2,1),lims(2,2)
                do k=1,ldim(3)!lims(3,1),lims(3,2)
                    ! phys = img_3%comp_addr_phys([h,k,l])
                    !if (abs(cabs(img_2%get_cmat_at(phys)) - cabs(img_3%get_cmat_at(phys))) >= 1.e-06 )then
                    !   print *,'wrong FFT i/j/k fftw cufft', h,k,l,img_3%get_cmat_at(phys) ,&
                    !      & img_2%get_cmat_at(phys)
                    ! passed = .false.
                    !                 call simple_stop("image FFT inverse PGI test  stopped")
                    if ((abs(real(ctmp2(i,j,k)) - real(ctmp(i,j,k))) >= tolerance ) .or. &
                        &(abs(aimag(ctmp2(i,j,k)) - aimag(ctmp(i,j,k))) >= tolerance)) then
                        print *,'wrong IFFT i/j/k fftw cufft', i,j,k,ctmp2(i,j,k),ctmp(i,j,k)
                    end if
                end do
            end do
        end do
        write(*,'(a)') '**info(simple_image PGI test): compare fftw and cufft (images) '
        ctmp = ctmp - ctmp2
        call gnufor_image(real(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(aimag(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(cabs(ctmp(:ldim(1),:ldim(2),1)), palette='gray')
        call gnufor_image(atan2(real(ctmp(:ldim(1),:ldim(2),1)),aimag(ctmp(:ldim(1),:ldim(2),1))), palette='gray')
        ! read(*,*)
        img_4 = img_2 - img_3
        call img_4%set_ft(.true.)
        if(doplot)call img_4%vis(geomorsphr=.false.)
        if(doplot)call img_4%vis(geomorsphr=.true.)
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftw3 inverse fft'

        !   call img%ifft
        !  write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '
        ! call ifft_pgi_cuda_test(img_2)
        !  print *," L1 norm sum", img.lonesum.img_2

        !  do i=1,ldim(1)
        ! !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
        !      do j=1,ldim(2)
        !          do k=1,ldim(3)
        !              if (abs( img_2%get([i,j,k]) - img_3%get([i,j,k])) >= 1.e-06 )then
        !                 print *,'wrong IFFT i/j/k fftw /cufft', i,j,img_3%get([i,j,k]) , img_2%get([i,j,k])
        ! !                 call simple_stop("image FFT PGI test stopped")
        !              end if
        !          end do
        !      end do
        !  end do
        !     VerbosePrint "In simple_cufft::ifft_pgi testing done"
        !read(*,*)
        call img_4%kill
        !      call img_5%kill
        call img_2%kill
        call img_3%kill
        call img%kill        ! write(*,'(a)') '**info(simple_image_unit_test, part 7): testing fftshift'

        if( allocated(ctmp) ) deallocate(ctmp)
        if( allocated(ctmp2) ) deallocate(ctmp2)
        ! passed = .false.
        ! msk=50
        ! call img%gauimg(10)
        ! write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
        ! if( doplot ) call img%vis
        ! write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
        ! call img%serialize(pcavec1, msk)
        ! call img%shift([-9.345,-5.786,0.])
        ! if( doplot ) call img%vis
        ! call img%shift([9.345,5.786,0.])
        ! call img%serialize(pcavec2, msk)
        ! if( doplot ) call img%vis
        ! if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
        if( .not. passed )  call simple_stop('pgi fft test failed')

    end subroutine test_pgi2


    subroutine test_pgi3( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        use gnufor2
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4!, img_5
        complex, allocatable ::  ctmp(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), lims(3,2), phys(3)
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)

        verbose=doplot

        VerbosePrint "In simple_cufft::ifft_pgi testing start"
        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,ld3]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        call img_4%new(ldim, smpd)
        ! call img_5%new(ldim, smpd)
        img=2.; img_2=0.; img_3=0.; img_4=0.
        call img%gauimg(10)
        call img%add_gauran(.5)
        ! call img%set_rmat_at_ind(8,5,1, 1.5)
        ! call img%set_rmat_at_ind(6,3,1, 1.)
        !$ allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        rtmp = img%get_rmat()
        call fftshift(rtmp)
        VerbosePrint "In simple_cufft::ifft_pgi testing fftshift"
        !   rtmp(8,5,1) = 1.5
        !  rtmp(6,3,1)= 1.
        !  rtmp(6:8,:,1)=0.
        call img_2%set_rmat(rtmp); call img_2%set_ft(.false.)
        call gnufor_image(rtmp(:,:,1),palette='gray')
        !call gnufor_image(rtmp(ldim(1):-1:1,ldim(2):-1:1,1),palette='gray')
        allocate(ctmp(ldim(1),ldim(2),ldim(3)))
        ctmp = cmplx(rtmp,abs(1.0-rtmp))

        call img_2%set_cmat(ctmp); call img_2%set_ft(.true.)
        call img_3%set_cmat(ctmp); call img_3%set_ft(.true.)
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        !     if(doplot)call img_2%vis
        !  ave = (img.lonesum.img_2)
        !  print *," L1 norm sum", ave
        !   write(*,'(a)') '**info(simple_image PGI test): testing fftshifter bwd'
        !     write(*,'(a)') 'Press enter to test ifftshift '
        !     read(*,*)
        !  call img_2%bwd_fftshift()
        !  img_4= img_2
        !  if(doplot)call img_4%vis(geomorsphr=.false.)
        !
        !  write(*,'(a)') 'Press enter to test cufft '
        ! read(*,*)
        write(*,'(a)') '**info(simple_image PGI test): testing cufft '

        call ifft_pgi_cuda_test(img_2)
        call simple_cuda_stop("In PGI TEST2 ",__FILENAME__,__LINE__)
        call img_3%bwd_ft()
        if(doplot)call img_3%vis()
        print *,"FFTW"
        call img_3%print_rmat
        !call img_2%ft2img('amp',img_4)
        if(doplot)call img_2%vis()
        print *,"CUFFT"
        call img_2%print_rmat()
        print *," test_pgi3 L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2



        !   write(*,'(a)') '**info(simple_image PGI test): cufft returned '
        !   write(*,'(a)') '**info(simple_image PGI test): testing cufft inverse fft '


        ! if(doplot)call img_2%vis
        ! print *," L1 norm sum", img.lonesum.img_2

        ! call img_3%kill
        ! call img_4%kill
        ! call img3d%kill

        !img_4 = img_2
        !    call img_2%ft2img('amp',img_4)
        !    if(doplot)call img_4%vis(geomorsphr=.false.)
        !  write(*,'(a)') '**info(simple_image PGI test): testing fftshifter fwd'
        !  write(*,'(a)') 'Press enter to test fftshift '
        !    read(*,*)
        !     call img_3%fwd_fftshift()
        !    img_4 = img_2
        !    if(doplot)call img_2%vis(geomorsphr=.false.)


        !     write(*,'(a)') 'Press enter to test fftw3 '
        !    read(*,*)


        !        call img_3%ft2img('amp',img_4)
        ! if(doplot)call img_3%vis(geomorsphr=.false.)

        lims=img_3%loop_lims(2)
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        write(*,'(a)') 'Press enter to test fftw '
        VerbosePrint "In simple_cufft::ifft_pgi testing set_cmat"
        !read(*,*)


        !  call img%subtr(img_2)
        !    call img%ft2img('amp',img_3)
        !    if(doplot)call img_3%vis(geomorsphr=.false.)
        !    write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '
        ldim = img_2%get_ldim()
        if(allocated(rtmp)) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        !$ allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)) )
        rtmp2= img_2%get_rmat(); rtmp= img_3%get_rmat()
        print *," test_pgi3 Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
        do i=1,ldim(1)
            !     VerbosePrint "In simple_cufft::ifft_pgi testing ",i
            do j=1,ldim(2)
                do k=1,ldim(3)
                    if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= tolerance )then
                        write(*,'(A,3I3,1x,2F20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
                        !                 call simple_stop("image FFT PGI test stopped")
                    end if
                end do
            end do
        end do

        write(*,'(a)') 'Press enter to test difference '
        !read(*,*)
        img_4= img_2 - img_3
        call img_4%print_rmat()
        if(doplot) call img_4%vis()




        if( allocated(rtmp) ) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        if( allocated(ctmp) ) deallocate(ctmp)

        ! call img_5%kill
        call img_4%kill
        call img_3%kill
        call img_2%kill
        call img%kill
        VerbosePrint "In simple_cufft::ifft_pgi testing done"
    end subroutine test_pgi3


    subroutine test_pgi4( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        use gnufor2
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3
        complex, allocatable :: ctmp(:,:,:)
        integer              :: i, j, k,h,l, cnt, lfny, ldim(3), phys(3), cdim(3)
        !   complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null()
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)
        integer(timer_int_kind):: t4
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)

        verbose=doplot
        if(verbose)t4=tic()
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        VerbosePrint "In test_pgi4 start"
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
        passed = .false.
        ldim = [ld1,ld2,1]
        !cdim = [9,16,1]
        call img%new(ldim, smpd)
        call img_2%new(ldim, smpd)
        call img_3%new(ldim, smpd)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        img=2.; img_2=0.; img_3=0.
        call img%gauimg(10)
        call img%gauimg2(5, floor(real(ld1)/5.),floor(real(ld2)/3.))
        call img%gauimg2(8, -floor(real(ld1)/5.),floor(real(ld2)/6.))

        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        allocate(rtmp(ldim(1),ldim(2),ldim(3)))
        rtmp = img%get_rmat()
        call ifftshift(rtmp)
        call img%set_rmat(rtmp)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        if(verbose)print *, "testing fftshift",toc()
        !call fftshift(rtmp)
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        VerbosePrint "In test_pgi4 testing fftshift"
        !call img%set_rmat_at_ind(17,16,1,10.)
        !rtmp(6,3,1) = 1.
        if(doplot)call img%vis()

        call img%set_ft(.true.)
        ! call img%fwd_ft()
        img_2 = img
        img_3 = img
        cdim = img%get_array_shape()
        !call img_2%set_rmat(rtmp); call img_2%set_ft(.false.)
        !call gnufor_image(rtmp(:,:,1),palette='gray')
        ! if(verbose)read(*,*)
        !call gnufor_image(rtmp(ldim(1):-1:1,ldim(2):-1:1,1),palette='gray')
        !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)

        if(doplot)call img%vis(geomorsphr=.true.)
        !  if(doplot)call img_2%vis(geomorsphr=.true.)
        !  if(doplot)call img_3%vis(geomorsphr=.true.)
        print *,' pgi4  ctmp array shape ', cdim
        allocate(ctmp(cdim(1),cdim(2),cdim(3)))
        ctmp = img%get_cmatfull()
        !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        !        call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !         call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !         call ifftshift(ctmp)
        !    call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        !    call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
        ! read(*,*)
        VerbosePrint "In test_pgi4 testing cufft ifft_pgi_cuda"
        call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
        ! if(verbose)read(*,*)
        if(verbose)t4=tic()
        call ifft_pgi_cuda_test3(img_2)
        if(img_2%is_ft())print *, " cufft ifft test failed",toc(t4)
        if(doplot)call img_2%vis()
        if(verbose)print *, "test completed backward cufft ifft",toc(t4)
        !read(*,*)
        if(doplot)call img_2%vis()
        if(verbose)print *,"CUFFT "
        !if(verbose)call img_2%print_rmat()
        !if(verbose)read(*,*)
        if(verbose)t4=tic()
        call img_3%bwd_ft()
        if(verbose)print *, "testing backward cufft ifft",toc(t4)
        if(doplot)call img_3%vis()
        if(verbose)print *,"FFTW"
        !f(verbose)call img_3%print_rmat

        print *," L1 norm sum FFTW v CUFFT", img_3.lonesum.img_2
        call img%set_ft(.false.)
        img = img_3-img_2
        if(doplot)call img%vis()
        ! if(verbose)read(*,*)

        if(allocated(rtmp)) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)) )
        rtmp2= img_2%get_rmat(); rtmp= img_3%get_rmat()
        ! print *," Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
        ! cnt=0
        ! do i=1,ldim(1)
        !     !     VerbosePrint "In test_pgi4 testing diff",i
        !     do j=1,ldim(2)
        !         do k=1,ldim(3)
        !             if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= 1.e-06 )then
        !                 write(*,'(A,3I0,1x,2ES20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
        !                 !                 call simple_stop("image FFT PGI test stopped")
        !                 cnt=cnt+1
        !             end if
        !         end do
        !     end do
        ! end do
        ! write(*,'(A,3I0,1x,2ES20.10)'), ' fftw /cufft', i,j,k,rtmp(i,j,k) , rtmp2(i,j,k)
        ! VerbosePrint "In test_pgi4 testing done"

        if( allocated(rtmp) ) deallocate(rtmp)
        if( allocated(rtmp2) ) deallocate(rtmp2)
        if( allocated(ctmp) ) deallocate(ctmp)

        call img_3%kill
        call img_2%kill
        call img%kill
        return
    end subroutine test_pgi4

    function ifft_pgi_cuda_test3a( cin ,ldim) result (routput)
        !    use simple_cufft
        use gnufor2
        complex,intent(in) ::cin(:,:,:)
        integer, intent(in) :: ldim(3)
        !       type(image), intent(inout)       :: self
        complex, device, allocatable :: inplace_d(:,:,:)
        complex, allocatable         :: cinput(:,:,:),coutput(:,:,:)
        complex, allocatable         :: routput(:,:,:)

        integer      :: i,j,k,h,l,f, xdim(2),ydim(2),zdim(2),n,m,p, lowerbounds(3), upperbounds(3)
        integer      :: nerrors, cdim(3),lims(3,2),phys(3)
        real         :: nscale
        complex      :: comp, cswap
        integer      :: planType, plan


        ! ldim = self%get_ldim()
        !do p=1,2
        n=ldim(1)
        m=ldim(2)
        p=ldim(3)
        ! allocate arrays on the host
        if (allocated(routput))deallocate(routput)
        if (allocated(cinput))deallocate(cinput)
        if (allocated(coutput))deallocate(coutput)
        allocate(routput(ldim(1),ldim(2),ldim(3)))
        allocate(cinput(n,m,p), coutput(n,m,p))
        ! allocate arrays on the device
        allocate(inplace_d(n,m,p))
        ! k=1
        ! write(*,*)"*"
        ! if(is_even(n))then
        lowerbounds(1) = lbound(cin,1);lowerbounds(2) = lbound(cin,2);lowerbounds(3) = lbound(cin,3);
        upperbounds(1) = ubound(cin,1); upperbounds(2) = ubound(cin,2); upperbounds(3) = ubound(cin,3)
        !     xbounds = (\ lbound(cin,1), ubound(cin,1) \) ! n/2
        !     ybounds = (\ lbound(cin,2), ubound(cin,2) \) !  m/2
        ! else
        xdim = [ -(n-1)/2, (n-1)/2 -1 ]
        ydim = [ -(m)/2, (m)/2 -1 ]
        zdim =  [ -(p)/2, (p)/2 -1 ]
        ! endif
        ! !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)
        !    print *,"test3 CIN    ", xdim, ydim, lowerbounds, upperbounds
        !    print *,"test3 CINPUT ", shape(cinput)
        cinput = cmplx(0.,0.)
        ! !  omp parallel do collapse(2) default(shared) private(h,k,phys,comp)
        do h=xdim(1),xdim(2)
            do k=ydim(1),ydim(2)
                do l=zdim(1),zdim(2)
                    if (h .ge. 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(m,0, k < 0)
                        phys(3) = l + 1 + MERGE(p,0, l < 0)
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(m,0, -k < 0)
                        phys(2) = -l + 1 + MERGE(p,0, -l < 0)

                    endif
                    if (phys(1)< lowerbounds(1) .or. phys(1) > upperbounds(1)) print *,"test3 CIN ind 1 err ", phys, h, k
                    if (phys(2)< lowerbounds(2) .or. phys(2) > upperbounds(2)) print *,"test3 CIN ind 2 err ", phys, h, k
                    if (phys(3)< lowerbounds(3) .or. phys(3) > upperbounds(3)) print *,"test3 CIN ind 3 err ", phys, h, k, l
                    comp=cin(phys(1),phys(2),phys(3))
                    if (h-xdim(1)+1 < 1 .or. h-xdim(1)+1 > n) print *,"test3 CINOUT ind 1 err ", h-xdim(1)+1, h, n
                    if (k-ydim(1)+1 < 1 .or. k-ydim(1)+1 > m) print *,"test3 CINOUT ind 2 err ", k-ydim(1)+1, k, m
                    if (l-zdim(1)+1 < 1 .or. l-zdim(1)+1 > p) print *,"test3 CINOUT ind 2 err ", l-zdim(1)+1, l, p

                    if(h>=0)then
                        cinput(h-xdim(1)+1,k-ydim(1)+1, l-zdim(1)+1) = comp
                    else
                        cinput(h-xdim(1)+1,k-ydim(1)+1, l-zdim(1)+1) = conjg(comp)
                    endif
                end do
            end do
        end do
        ! !  omp end parallel do

        ! ! copy input to device
        inplace_d = cinput
        ! ! write(*,*)"*Initialize the plan for complex to complex transform"
        call cufftPlan2D(plan,n,n,CUFFT_C2C)
        ! ! write(*,*)" Execute  Backward transform in place"
        call cufftExec(plan,CUFFT_C2C,inplace_d,inplace_d,CUFFT_INVERSE)


        ! ! Copy results back to host
        coutput = inplace_d

        ! ! Save to image obj

        ! ! omp parallel do collapse(2) default(shared) private(i,j)
        do i=1,ldim(1)
            do j=1,ldim(2)
                do k=1,ldim(3)
                    !       call self%set_rmat_at_ind(i, j, 1, cabs( cinout(i,j)) )
                    routput(i, j, k)= cabs( coutput(i,j,k))
                    !  if (j==m)print *," ",  cabs(cinout(i,j))
                end do
            end do
        end do

        ! ! omp end parallel do
        ! !call self%set_ft(.false.)

        ! !  call self%shift_phorig

        ! ! release memory on the host and on the device


        deallocate( cinput, coutput, inplace_d ,stat=alloc_stat)

        ! ! Destroy the plans
        call cufftDestroy(plan)
        print *,"In simple_image::ifft_pgi done "

    end function ifft_pgi_cuda_test3a



    subroutine ifft_pgi_cuda_test4( self) result (routput)
        !    use simple_cufft
        use gnufor2
        !        complex,intent(in) ::cin(:,:,:)
        !        integer, intent(in) :: ldim(3)
        type(image), intent(inout)       :: self
        complex, device, allocatable :: inplace_d(:,:,:)
        complex, allocatable         :: cinout(:,:,:), cin(:,:,:)
        real, allocatable         :: routput(:,:,:)

        integer      :: i,j,k,g,f,h, xdim(2),ydim(2),zdim(2),ii,jj,kk
        integer      :: nerrors,  ldim(3), cdim(3),lims(3,2),phys(3), cmat_shape(3)
        real         :: nscale
        complex      :: comp, cswap
        integer      :: planType, plan


        ldim = self%get_ldim()
        !do p=1,2
        cmat_shape = self%get_array_shape()

        ! allocate arrays on the host
        if (allocated(routput))deallocate(routput)
        if (allocated(cinout))deallocate(cinout)
        allocate(cin(cmat_shape(1),cmat_shape(2),cmat_shape(3)))
        cin = self%get_cmatfull()
        allocate(routput(ldim(1),ldim(2),ldim(3)))
        allocate(cinout(ldim(1),ldim(2),ldim(3)))
        ! allocate arrays on the device
        allocate(inplace_d(ldim(1),ldim(2),ldim(3)))
        ! k=1
        ! write(*,*)"*"
        ! if(is_even(n))then
        ! lowerbounds(1) = lbound(cin,1);lowerbounds(2) = lbound(cin,2);
        ! upperbounds(1) = ubound(cin,1); upperbounds(2) = ubound(cin,2)
        !     xbounds = (\ lbound(cin,1), ubound(cin,1) \) ! n/2
        !     ybounds = (\ lbound(cin,2), ubound(cin,2) \) !  m/2
        ! else
        xdim = [ -floor(real(ldim(1)-1)/2.) , ceiling(real(ldim(1)-1)/2.) - 1  ]
        ydim = [ -floor(real(ldim(2))/2.)   , ceiling(real(ldim(2))/2.) - 1      ]
        zdim = [ -floor(real(ldim(3))/2.)   , ceiling(real(ldim(3))/2.) - 1      ]

        ! endif
        ! !allocate(cinout(-xdim:xdim,-ydim:ydim,1:1),stat=istat)
        print *,"test4 CIN    ", xdim, ydim, cmat_shape
        print *,"test4 CINPUT ", shape(cinout)
        cinout = cmplx(0.,0.)
        phys=1
        ! omp parallel do collapse(3) default(shared)  private(i,j,k,f,g,h,phys,comp)
        do i=1,ldim(1)
            do j=1,ldim(2)
                do k=1,ldim(3)

                    h=xdim(1)+(i-1)
                    g=ydim(1)+(j-1)
                    f=zdim(1)+(k-1)
                    phys=1
                    if (h .ge. 0) then
                        phys(1) = h + 1
                        phys(2) = g + 1 + MERGE(ldim(2),0, g < 0)
                        phys(3) = f + 1 + MERGE(ldim(3),0, f < 0)
                    else
                        phys(1) = -h + 1
                        phys(2) = -g + 1 + MERGE(ldim(2),0, -g < 0)
                        phys(3) = -f + 1 + MERGE(ldim(3),0, -f < 0)
                    endif
                    ! #ifdef _DEBUG
                    if (phys(1) < 1 .or. phys(1) > cmat_shape(1) ) print *," CIN ind 1 err ", i, h, phys(1)
                    if (phys(2) < 1 .or. phys(2) > cmat_shape(2)) print *," CIN ind 2 err ", j, g, phys(2)
                    if (phys(3) < 1 .or. phys(3) > cmat_shape(3)) print *," CIN ind 3 err ", k,f, phys(3)
                    ! #endif
                    comp=cin(phys(1),phys(2),phys(3))
                    ! #ifdef _DEBUG
                    if (i < 1 .or. i > ldim(1)) print *," CINOUT ind 1 err ", i, h
                    if (j < 1 .or. j > ldim(2)) print *," CINOUT ind 2 err ", j, g
                    if (k < 1 .or. k > ldim(3)) print *," CINOUT ind 3 err ", k, f
                    ! #endif
                    if(h>=0)then
                        cinout(i,j,k) = comp
                    else
                        cinout(i,j,k) = conjg(comp)
                    endif
                end do
            end do
        end do
        ! omp end parallel do

        ! ! copy input to device
        inplace_d = cinout
        ! ! write(*,*)"*Initialize the plan for complex to complex transform"
        call cufftPlan3D(plan,ldim(1),ldim(2),ldim(3),CUFFT_C2C)
        ! ! write(*,*)" Execute  Backward transform in place"
        call cufftExec(plan,CUFFT_C2C,inplace_d,inplace_d,CUFFT_INVERSE)


        ! ! Copy results back to host
        cinout = inplace_d

        ! ! Save to image obj

        !$omp parallel do collapse(3) default(shared) private(i,j,k)
        do i=1,ldim(1)
            do j=1,ldim(2)
                do k=1,ldim(3)
                    !       call self%set_rmat_at_ind(i, j, 1, cabs( cinout(i,j)) )
                    routput(i, j, k)= cabs( cinout(i,j,k))
                    !  if (j==m)print *," ",  cabs(cinout(i,j))
                end do
            end do
        end do
        !$omp end parallel do
        call self%set_ft(.false.)
        call self%set_rmat(routput)
        call self%shift_phorig

        ! ! release memory on the host and on the device
        deallocate( cin,cinout, routput, inplace_d ,stat=alloc_stat)

        ! ! Destroy the plans
        call cufftDestroy(plan)
        print *,"In simple_image::ifft_pgi done "

    end subroutine ifft_pgi_cuda_test4


    subroutine test_pgi_ifft(  doplot )
        use simple_fftshifter
        use simple_image, only: image
        use gnufor2

        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4
        complex, allocatable :: ctmp(:,:,:)
        integer              :: i, j, k,h,l,p, cnt, lfny, ldim(3), phys(3), cdim(3)
        !   complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null()
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:),rtmp3(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)
        !   type(timer_cuda) :: ctimer
        !    type(cudaEvent) :: t4

        integer(8) :: t1

        do p=5,8
            t1=tic()
            ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
            passed = .false.
            ldim = [2**p ,2**(p) ,1]
            !cdim = [9,16,1]
            call img%new(ldim, smpd)
            call img_2%new(ldim, smpd)
            call img_3%new(ldim, smpd)
            call img_4%new(ldim, smpd)

            ! call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            img=2.; img_2=0.; img_3=0.
            call img%gauimg(10)
            call img%gauimg2(5, floor(real(ldim(1))/5.),floor(real(ldim(2))/3.))
            call img%gauimg2(8, -floor(real(ldim(1))/5.),floor(real(ldim(2))/6.))
            call img%gauimg2(2, -floor(real(ldim(1))/8.),-floor(real(ldim(2))/9.))
            call img%gauimg2(3, floor(real(ldim(1))/5.),-floor(real(ldim(2))/4.))
            call img%add_gauran(100.)
            ! call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            ! allocate(rtmp(n,m,ldim(3)))
            ! rtmp = img%get_rmat()
            !  call ifftshift(rtmp)
            !  call img%set_rmat(rtmp)
            !  call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            !     print *, "testing fft ",ctimer%tocU()
            !call fftshift(rtmp)
            !  call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            ! print *,"In test_pgi4 testing fftshift"
            !call img%set_rmat_at_ind(17,16,1,10.)
            !rtmp(6,3,1) = 1.
            !     if(doplot)call img%vis()

            ! call img%set_ft(.true.)
            call img%fwd_ft()
            img_2 = img
            img_3 = img
            img_4 = img
            cdim = img%get_array_shape()
            !    if(doplot)call img%vis(geomorsphr=.false.)
            ! !call img_2%set_rmat(rtmp); call img_2%set_ft(.false.)
            ! !call gnufor_image(rtmp(:,:,1),palette='gray')
            ! if(verbose)read(*,*)
            ! !call gnufor_image(rtmp(n:-1:1,m:-1:1,1),palette='gray')
            ! !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            print *, "testing backward cufft and fftw  dim size:", ldim(1)

            ! !  if(doplot)call img%vis(geomorsphr=.true.)
            ! !  if(doplot)call img_2%vis(geomorsphr=.true.)
            ! !  if(doplot)call img_3%vis(geomorsphr=.true.)
            ! print *,' pgi4  ctmp array shape ', cdim
            t1=tic()
            call img%bwd_ft()
            ! print *, "testing backward fftw ifft time(s):",toc(t4)
            print *, "test completed fftw ifft as image subroutine", toc(t1)
            if(doplot)call img%vis()

            t1=tic()
            allocate(ctmp(cdim(1),cdim(2),cdim(3)))
            ctmp = img%get_cmatfull()
            ! !call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            ! !        call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
            ! !         call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
            ! !         call ifftshift(ctmp)
            ! !    call gnufor_image(real(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
            ! !    call gnufor_image(aimag(ctmp(:cdim(1),:cdim(2),1)),palette='gray')
            ! ! read(*,*)
            ! print *,"In test_pgi4 testing cufft ifft_pgi_cuda"
            ! !  call simple_cuda_stop("In TEST_PGI4 ",__FILENAME__,__LINE__)
            ! ! if(verbose)read(*,*)
            ! ! t4=tic()
            allocate(rtmp(ldim(1),ldim(2),ldim(3)))
            rtmp = ifft_pgi_cuda_test3a(ctmp,ldim)
            ! ! ! call gnufor_image(rtmp(:ldim(1),:ldim(2),1),palette='gray')
            call img_2%set_rmat(rtmp)
            call img_2%shift_phorig
            if(doplot)call img_2%vis()

            t1=tic()
            call ifft_pgi_cuda_test4(img_3)
            ! if(img_2%is_ft())then
            !     print *, " cufft ifft test failed"
            ! ! ! if(doplot)call img_2%vis()
            ! else
            print *, "test completed cufft ifft as function", toc(t1)
            ! end if
            if(doplot)call img_3%vis()


            t1=tic()
            call img_4%ifft_pgi_cuda()
            print *, "test completed cufft ifft as image subroutine", toc(t1)

            ! ! call gnufor_image(rtmp(:ldim(1),:ldim(2),1),palette='gray')
            if(doplot)call img_4%vis()

            print *," L2 norm sum FFTW v CUFFT func", img.lonesum.img_2
            print *," L2 norm sum FFTW v CUFFT subr", img.lonesum.img_4
            print *," L2 norm sum FFTW v CUFFT image class", img.lonesum.img_4

            if(allocated(rtmp)) deallocate(rtmp)
            if( allocated(rtmp2) ) deallocate(rtmp2)
            if( allocated(rtmp3) ) deallocate(rtmp3)
            allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)),rtmp3(ldim(1),ldim(2),ldim(3)) )
            rtmp = img%get_rmat();  rtmp2= img_4%get_rmat(); rtmp3 = img_4%get_rmat()

             if(doplot)call gnufor_image(rtmp(:ldim(1),:ldim(2),1),rtmp2(:ldim(1),:ldim(2),1), palette='hot')
             if(doplot)call gnufor_image(rtmp(:ldim(1),:ldim(2),1),rtmp3(:ldim(1),:ldim(2),1), palette='ocean')

            !read(*,*)
            !  print *," Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
            cnt=0
            do i=1,ldim(1)
                ! !     !     VerbosePrint "In test_pgi4 testing diff",i
                do j=1,ldim(2)
                    do k=1,ldim(3)
                        if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= 5.e-05 )then
                            write(*,'(A,I5,I5,I5,1x,2ES20.8)'), ' fftw /cufft diff ', i, j, k,&
                                rtmp(i,j,k) , rtmp2(i,j,k)
                            !            call simple_stop("image FFT PGI test stopped")
                            cnt=cnt+1
                            !    rtmp(i,j,k)= rtmp(i,j,k) - rtmp2(i,j,k)
                        end if

                    end do
                end do
            end do
            if( cnt > 0 ) write(*,'(A,I0,1x,F10.2,"%")'), ' fftw differs from cufft ', cnt, &
                real(cnt*100)/real(ldim(1)*ldim(2)*ldim(3))

            cnt=0
            do i=1,ldim(1)
                ! !     !     VerbosePrint "In test_pgi4 testing diff",i
                do j=1,ldim(2)
                    do k=1,ldim(3)
                        if (abs( rtmp(i,j,k) - rtmp3(i,j,k)) >= 5.e-05 )then
                            write(*,'(A,I5,I5,I5,1x,2ES20.8)'), ' fftw /cufft diff ', &
                                i, j, k,rtmp(i,j,k) , rtmp3(i,j,k)
                            !            call simple_stop("image FFT PGI test stopped")
                            cnt=cnt+1
                            !    rtmp(i,j,k)= rtmp(i,j,k) - rtmp2(i,j,k)
                        end if

                    end do
                end do
            end do
            if( cnt > 0 ) write(*,'(A,I0,1x,F10.2,"%")'), ' fftw differs from cufft ', cnt, &
                real(cnt*100)/real(ldim(1)*ldim(2)*ldim(3))

            !  call gnufor_image(rtmp(:ldim(1),:ldim(2),1),palette='gray')
            ! ! VerbosePrint "In test_pgi4 testing done"

            call img%kill
            call img_3%kill
            call img_2%kill
            call img_4%kill
            if( allocated(rtmp) ) deallocate(rtmp)
            if( allocated(rtmp2) ) deallocate(rtmp2)
            if( allocated(rtmp3) ) deallocate(rtmp3)
            if( allocated(ctmp) ) deallocate(ctmp)

            write(*,*) "In test_pgi4 testing done"
            ! read(*,*)
            call exec_cmdline("pkill gnuplot_x11")
        end do
    end subroutine test_pgi_ifft

end module simple_cufft_test2
