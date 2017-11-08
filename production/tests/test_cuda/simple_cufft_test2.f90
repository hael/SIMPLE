module simple_cufft_test2
    include 'simple_lib.f08'
    use simple_image,       only: image
    use simple_fftshifter
    use gnufor2
    use simple_cufft
    implicit none

#ifdef PGI
    public :: cufftTest, test_pgi4
#endif
#include "simple_local_flags.inc"
contains
    !     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
    !
    ! NVIDIA CORPORATION and its licensors retain all intellectual property
    ! and proprietary rights in and to this software, related documentation
    ! and any modifications thereto.  Any use, reproduction, disclosure or
    ! distribution of this software and related documentation without an express
    ! license agreement from NVIDIA CORPORATION is strictly prohibited.
    !
    subroutine cufftTest
        !  use simple_cufft

        complex, allocatable :: a(:,:),b(:,:)
        complex, device, allocatable :: a_d(:,:), b_d(:,:)
        real, allocatable :: rtimes(:)
        integer i,j,it
        integer :: nerrors
        integer(8)::t1
        real :: minr, mini
        integer,parameter :: n=600, iters=30
        integer :: plan, planType
        allocate(rtimes(iters-1))
        do it=1,iters
            t1=tic()
            ! allocate arrays on the host
            allocate(a(n,n), b(n,n))

            ! allocate arrays on the device
            allocate(a_d(n,n), b_d(n,n))

            !initialize arrays on host
            do i = 1, n
                do j=1,n
                    a(i,j) = cmplx(cos((i-1) * atan2(0.0,-1.0) / n),cos((j-1) * atan2(0.0,-1.0) / n))
                end do
            end do

            !copy arrays to device
            a_d = a

            ! Print initial array
            !print *, "Array A:"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") a

            ! set planType to either single or double precision
            !if (sp == singlePrecision) then
            planType = CUFFT_C2C
            !else
            !   planType = CUFFT_Z2Z
            !endif

            ! initialize the plan and execute the FFTs.

            call cufftPlan2D(plan,n,n,planType)
            call cufftExec(plan,planType,a_d,b_d,CUFFT_INVERSE)

            ! Copy results back to host
            b = b_d
            !print *, "Forward B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            call cufftExec(plan,planType,b_d,b_d,CUFFT_FORWARD)

            ! Copy results back to host
            b = b_d
            !print *, "Inverse B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b

            ! Scale
            b = b / (n*n)
            !print *, "Scaled B"
            !write (*,"(8('(',f6.3,',',f6.3,')',1x))") b
            if (it .ne. 1 ) rtimes(it-1) =  toc(t1)
            nerrors = 0;minr=1e20;mini=1e20
            do i = 1, n
                do j=1,n
                    if ( abs(real(a(i,j)) - real(b(i,j))) .gt. 2.0e-6 ) then
                        if ( abs(real(a(i,j)) - real(b(i,j))) < minr ) then
                            minr = abs(real(a(i,j)) - real(b(i,j)))
                        end if
                        nerrors = nerrors + 1
                    endif
                    if ( abs(imag(a(i,j)) - imag(b(i,j))) .gt. 2.0e-6 ) then
                        if ( abs(imag(a(i,j)) - imag(b(i,j))) < mini ) then
                            mini = abs(imag(a(i,j)) - imag(b(i,j)))
                        end if
                        nerrors = nerrors + 1
                    endif
                end do
            end do
            if (nerrors .ne. 0) then
                print *, "Test FAILED", minr, mini
            else
                if (it == 1 ) then
                    print *, "Test PASSED (first)",toc(t1)
                else
                    print *, "Test PASSED", rtimes(it-1)
                endif

            endif

            !release memory on the host and device
            deallocate(a, b, a_d, b_d)

            ! Destroy the plan
            call cufftDestroy(plan)
        end do
        print *," cufftTest  600x600 Cmplx FWD/BWD In-place, average runtime (excl first)", sum(rtimes)/real(iters-1)
        deallocate(rtimes)
    end subroutine cufftTest

    function ifft_pgi_cuda_test3( cin ,ldim) result (routput)
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

    end function ifft_pgi_cuda_test3



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




    subroutine test_pgi4( ld1, ld2, ld3, doplot )
        use simple_fftshifter
        use simple_image, only: image
        !    use simple_timer_cuda
        use gnufor2
        integer, intent(in)  :: ld1, ld2, ld3
        logical, intent(in)  :: doplot
        type(image)          :: img, img_2, img_3, img_4
        complex, allocatable :: ctmp(:,:,:)
        integer              :: i, j, k,h,l,p, cnt, lfny, ldim(3), phys(3), cdim(3)
        !   complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null()
        real, allocatable    :: rtmp(:,:,:),rtmp2(:,:,:),rtmp3(:,:,:)
        real                 :: smpd=2.
        logical              :: passed, test(6)
        !    type(timer_cuda)::ctimer
        !    type(cudaEvent):: t4

        integer(8) :: t1


        do p=5,8
            t1=tic()
            ! write(*,'(a)') '**info(simple_image_unit_test, PGI CUFFT 1): testing '
            passed = .false.
            ldim = [2**p ,2**(p) ,8]
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
              rtmp = ifft_pgi_cuda_test3(ctmp,ldim)
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

            print *," L2 norm sum FFTW v CUFFT", img.lonesum.img_2
            print *," L2 norm sum FFTW v CUFFT", img.lonesum.img_4
            print *," L2 norm sum FFTW v CUFFT", img.lonesum.img_4

            if(allocated(rtmp)) deallocate(rtmp)
            if( allocated(rtmp2) ) deallocate(rtmp2)
            if( allocated(rtmp3) ) deallocate(rtmp3)
            allocate(rtmp(ldim(1),ldim(2),ldim(3)),rtmp2(ldim(1),ldim(2),ldim(3)),rtmp3(ldim(1),ldim(2),ldim(3)) )
            rtmp = img%get_rmat();  rtmp2= img_4%get_rmat(); rtmp3 = img_4%get_rmat()

            call gnufor_image(rtmp(:ldim(1),:ldim(2),1),rtmp2(:ldim(1),:ldim(2),1), palette='hot')
            call gnufor_image(rtmp(:ldim(1),:ldim(2),1),rtmp3(:ldim(1),:ldim(2),1), palette='ocean')

            read(*,*)
            !  print *," Ldim/shapes ", ldim, shape(rtmp), shape(rtmp2)
            cnt=0
            do i=1,ldim(1)
                ! !     !     VerbosePrint "In test_pgi4 testing diff",i
                do j=1,ldim(2)
                    do k=1,ldim(3)
                        if (abs( rtmp(i,j,k) - rtmp2(i,j,k)) >= 5.e-05 )then
                            write(*,'(A,I5,I5,I5,1x,2ES20.8)'), ' fftw /cufft diff ', i, j, k,rtmp(i,j,k) , rtmp2(i,j,k)
                            !            call simple_stop("image FFT PGI test stopped")
                            cnt=cnt+1
                            !    rtmp(i,j,k)= rtmp(i,j,k) - rtmp2(i,j,k)
                        end if

                    end do
                end do
            end do
            if( cnt > 0 ) write(*,'(A,I0,1x,F10.2,"%")'), ' fftw differs from cufft ', cnt, real(cnt*100)/real(ldim(1)*ldim(2)*ldim(3))

            cnt=0
            do i=1,ldim(1)
                ! !     !     VerbosePrint "In test_pgi4 testing diff",i
                do j=1,ldim(2)
                    do k=1,ldim(3)
                        if (abs( rtmp(i,j,k) - rtmp3(i,j,k)) >= 5.e-05 )then
                            write(*,'(A,I5,I5,I5,1x,2ES20.8)'), ' fftw /cufft diff ', i, j, k,rtmp(i,j,k) , rtmp3(i,j,k)
                            !            call simple_stop("image FFT PGI test stopped")
                            cnt=cnt+1
                            !    rtmp(i,j,k)= rtmp(i,j,k) - rtmp2(i,j,k)
                        end if

                    end do
                end do
            end do
            if( cnt > 0 ) write(*,'(A,I0,1x,F10.2,"%")'), ' fftw differs from cufft ', cnt, real(cnt*100)/real(ldim(1)*ldim(2)*ldim(3))

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
    end subroutine test_pgi4

end module simple_cufft_test2
