! projection of 3D volumes in the Fourier domain by convolution interpolation
! to generate band-pass limited Cartesian and polar 2D Fourier transforms,
! high-level routines
module simple_projector_hlev
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_projector,  only: projector
use simple_parameters, only: params_glob
implicit none

contains

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function reproject( vol, o, top ) result( imgs )
        use simple_oris, only: oris
        use simple_ori,  only: ori
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        integer, optional, intent(in)    :: top     !< stop index
        type(image),       allocatable :: imgs(:)   !< resulting images
        type(image),       allocatable :: imgs_pad(:)
        type(ori)        :: o2
        type(projector)  :: vol_pad
        integer          :: n, i, ithr, ldim(3), boxpd, ldim_pd(3), box
        real             :: smpd
        ldim    = vol%get_ldim()
        box     = ldim(1)
        smpd    = vol%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(box / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        ! padding & fft
        call vol_pad%new(ldim_pd, smpd)
        call vol%pad(vol_pad)
        if( params_glob%gridding.eq.'yes' )then
            ! corrects for interpolation function
            call vol_pad%div_w_instrfun(params_glob%interpfun, alpha=KBALPHA, padded_dim=boxpd)
        endif
        call vol_pad%fft
        call vol_pad%mul(real(boxpd)) ! correct for FFTW convention
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), imgs_pad(nthr_glob) )
        ! construct thread safe images
        do i=1,n
            call imgs(i)%new([box,box,1], smpd, wthreads=.false.)
        end do
        do ithr=1,nthr_glob
            call imgs_pad(ithr)%new([boxpd,boxpd,1], smpd, wthreads=.false.)
        end do
        ! prepare for projection
        call vol_pad%expand_cmat(KBALPHA)
        write(logfhandle,'(A)') '>>> GENERATING PROJECTIONS'
        !$omp parallel do schedule(static) default(shared)&
        !$omp private(i,ithr,o2) proc_bind(close)
        do i=1,n
            if(o%isthere('state'))then
                if(o%get_state(i)==0)cycle
            endif
            ! get thread index
            ithr = omp_get_thread_num() + 1
            ! extract central secion
            call o%get_ori(i, o2)
            call vol_pad%fproject_serial(o2, imgs_pad(ithr))
            call o2%kill
            ! back FT
            call imgs_pad(ithr)%ifft()
            ! clip
            call imgs_pad(ithr)%clip(imgs(i))
        end do
        !$omp end parallel do
        ! destruct
        do ithr=1,nthr_glob
            call imgs_pad(ithr)%kill
        end do
        deallocate(imgs_pad)
        call vol_pad%kill_expanded
        call vol_pad%kill
    end function reproject

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( vol, o, shvec ) result( rovol )
        use simple_ori, only: ori
        class(image),   intent(inout) :: vol      !< volume to project
        class(ori),     intent(inout) :: o        !< orientation
        real, optional, intent(in)    :: shvec(3) !< 3D shift vector
        type(projector)  :: vol_pad
        type(image)      :: rovol_pad, rovol
        integer          :: ldim(3), ldim_pd(3), boxpd
        real             :: smpd
        ldim    = vol%get_ldim()
        smpd    = vol%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(ldim(1) / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        call rovol%new(ldim, smpd)
        call rovol_pad%new(ldim_pd, smpd)
        call vol_pad%new(ldim_pd, smpd)
        call vol%pad(vol_pad)
        if( params_glob%gridding.eq.'yes' )then
            ! corrects for interpolation function
            call vol_pad%div_w_instrfun(params_glob%interpfun, alpha=KBALPHA, padded_dim=boxpd)
        endif
        call vol_pad%fft
        call vol_pad%expand_cmat(KBALPHA)
        call rotvol_slim( vol_pad, rovol_pad, rovol, o, shvec )
        call vol_pad%kill_expanded
        call vol_pad%kill
        call rovol_pad%kill
    end function rotvol

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    subroutine rotvol_slim( vol_pad, rovol_pad, rovol, o, shvec )
        use simple_ori, only: ori
        class(projector), intent(inout) :: vol_pad
        class(image),     intent(inout) :: rovol_pad, rovol
        class(ori),       intent(inout) :: o
        real, optional,   intent(in)    :: shvec(3)
        integer :: h,k,l,nyq,lims(3,2),logi(3),phys(3),lim_sq
        real    :: loc(3), rmat(3,3)
        logical :: l_shvec_present
        l_shvec_present = present(shvec)
        rovol_pad = cmplx(0.,0.)
        rovol     = 0.
        lims      = vol_pad%loop_lims(2)
        nyq       = vol_pad%get_lfny(1)
        lim_sq    = (nyq + 1)**2
        rmat      = o%get_mat()
        !$omp parallel do collapse(2) default(shared) private(h,k,l,loc,logi,phys)&
        !$omp schedule(dynamic) proc_bind(close)
        do l=lims(3,1),lims(3,2)
            do k=lims(2,1),lims(2,2)
                if(l*l+k*k > lim_sq) cycle
                do h=lims(1,1),lims(1,2)
                    if(h*h+k*k+l*l > lim_sq) cycle
                    logi = [h,k,l]
                    phys = rovol_pad%comp_addr_phys(logi)
                    loc  = matmul(real(logi), rmat)
                    if( l_shvec_present )then
                        call rovol_pad%set_fcomp(logi, phys, vol_pad%interp_fcomp(loc) * rovol_pad%oshift(loc, shvec))
                    else
                        call rovol_pad%set_fcomp(logi, phys, vol_pad%interp_fcomp(loc))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        call rovol_pad%ifft()
        call rovol_pad%clip(rovol)
    end subroutine rotvol_slim

end module simple_projector_hlev
