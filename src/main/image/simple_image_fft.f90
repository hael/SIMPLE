submodule (simple_image) simple_image_fft
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
use, intrinsic :: iso_c_binding
implicit none
contains

    !===========================
    ! Core FFT routines
    !===========================

    module subroutine fwd_ft(self)
        class(image), intent(inout) :: self
        if( self%ft ) return
        if( shift_to_phase_origin ) call self%shift_phorig
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        ! now scale the values so that a ifft() of the output yields the
        ! original image back following FFTW
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fwd_ft

    module subroutine bwd_ft( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
            self%ft = .false.
            if( shift_to_phase_origin ) call self%shift_phorig
        endif
    end subroutine bwd_ft

    module subroutine fft_noshift( self )
        class(image), intent(inout) :: self
        if( self%ft ) return
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fft_noshift

    !===========================
    ! FT / image conversion
    !===========================

    module subroutine img2ft( self, img )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: img
        integer :: h,k,l,lims(3,2),logi(3),phys(3)
        integer :: xcnt,ycnt,zcnt
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; img2ft')
        endif
        call img%zero_and_flag_ft
        xcnt = 0
        ycnt = 0
        zcnt = 0
        lims = self%loop_lims(3)
        do h=lims(1,1),lims(1,2)
            xcnt = xcnt + 1
            if( xcnt > self%ldim(1) ) cycle
            ycnt = 0
            do k=lims(2,1),lims(2,2)
                ycnt = ycnt + 1
                if( ycnt > self%ldim(2) ) cycle
                zcnt = 0
                do l=lims(3,1),lims(3,2)
                    zcnt = zcnt + 1
                    if( zcnt > self%ldim(3) ) cycle
                    logi = [h,k,l]
                    phys = self%comp_addr_phys(logi)
                    call img%set_fcomp(logi, phys, cmplx(self%rmat(xcnt,ycnt,zcnt),0.))
                end do
            end do
        end do
    end subroutine img2ft

    module subroutine ft2img( self, which, img )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: which
        class(image),     intent(inout) :: img
        integer :: h,mh,k,mk,l,ml,lims(3,2),inds(3),phys(3)
        integer :: which_flag
        logical :: didft
        complex :: comp
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; ft2img')
        endif
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        select case(which)
            case ('real')
                which_flag = 0
            case('power')
                which_flag = 1
            case('sqrt')
                which_flag = 2
            case ('log')
                which_flag = 3
            case('phase')
                which_flag = 4
            case DEFAULT
                THROW_HARD('unsupported mode: '//trim(which)//'; ft2img')
        end select
        call img%zero_and_unflag_ft
        lims = self%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        ml   = abs(lims(3,1))
        if( .not.self%wthreads .and. self%is_2d() )then
            do k=lims(2,1),lims(2,2)
                inds(2) = min(max(1,k+mk+1),self%ldim(2))
                do h=lims(1,1),lims(1,2)
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    comp    = self%get_fcomp2D(h,k)
                    select case(which_flag)
                        case(0)
                            img%rmat(inds(1),inds(2),1) = real(comp)
                        case(1)
                            img%rmat(inds(1),inds(2),1) = csq(comp)
                        case(2)
                            img%rmat(inds(1),inds(2),1) = sqrt(csq(comp))
                        case(3)
                            img%rmat(inds(1),inds(2),1) = log(csq(comp))
                        case(4)
                            img%rmat(inds(1),inds(2),1) = phase_angle(comp)
                    end select
                end do
            end do
        else
            !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,comp,inds)&
            !$omp schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self%get_fcomp([h,k,l],phys)
                        inds(1) = min(max(1,h+mh+1),self%ldim(1))
                        inds(2) = min(max(1,k+mk+1),self%ldim(2))
                        inds(3) = min(max(1,l+ml+1),self%ldim(3))
                        select case(which_flag)
                        case (0)
                            call img%set(inds,real(comp))
                        case(1)
                            call img%set(inds,csq(comp))
                        case(2)
                            call img%set(inds,sqrt(csq(comp)))
                        case (3)
                            call img%set(inds,log(csq(comp)))
                        case(4)
                            call img%set(inds,phase_angle(comp))
                        end select
                    end do
                end do
            end do
            !$omp end parallel do
        endif
        if( didft ) call self%ifft()
    end subroutine ft2img

    !===========================
    ! Padding / normalization
    !===========================

    module subroutine pad_fft( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine pad_fft

    subroutine norm_noise_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        real    :: sdev_noise
        call self%norm_noise(lmsk, sdev_noise)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine norm_noise_pad_fft

    !>  \brief  expand_ft is for getting a Fourier plane using the old SIMPLE logics
    module function expand_ft( self ) result( fplane )
        class(image), intent(in) :: self
        complex, allocatable :: fplane(:,:)
        integer :: xdim, ydim, h, k, phys(3)
        if(is_even(self%ldim(1)))then
            xdim = self%ldim(1)/2
            ydim = self%ldim(2)/2
        else
            xdim = (self%ldim(1)-1)/2
            ydim = (self%ldim(2)-1)/2
        endif
        allocate(fplane(-xdim:xdim,-ydim:ydim))
        fplane = cmplx(0.,0.)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                fplane(h,k) = self%get_fcomp([h,k,0],phys)
            end do
        end do
    end function expand_ft

end submodule simple_image_fft
