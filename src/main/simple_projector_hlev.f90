! projection of 3D volumes in the Fourier domain by convolution interpolation
! to generate band-pass limited Cartesian and polar 2D Fourier transforms,
! high-level routines
module simple_projector_hlev
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_image,      only: image
use simple_oris,       only: oris
use simple_params,     only: params
use simple_gridding,   only: prep4cgrid
use simple_ori,        only: ori
use simple_projector,  only: projector
use simple_kbinterpol, only: kbinterpol
implicit none

contains

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function project( vol, o, p, top ) result( imgs )
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(inout) :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        type(image),       allocatable :: imgs(:)   !< resulting images
        type(image),       allocatable :: imgs_pad(:)
        type(projector)  :: vol_pad
        type(kbinterpol) :: kbwin
        integer          :: n, i, ithr
        kbwin = kbinterpol(KBWINSZ, p%alpha)
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call prep4cgrid(vol, vol_pad, kbwin)
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), imgs_pad(p%nthr), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('project; simple_projector')
        ! construct thread safe images
        do i=1,n
            call imgs(i)%new([p%box,p%box,1], p%smpd, wthreads=.false.)
        end do
        do ithr=1,p%nthr
            call imgs_pad(ithr)%new([p%boxpd,p%boxpd,1], p%smpd, wthreads=.false.)
        end do
        ! prepare for projection
        call vol_pad%expand_cmat(p%alpha)
        write(*,'(A)') '>>> GENERATES PROJECTIONS'
        !$omp parallel do schedule(static) default(shared)&
        !$omp private(i,ithr) proc_bind(close)
        do i=1,n
            ! get thread index
            ithr = omp_get_thread_num() + 1
            ! extract central secion
            call vol_pad%fproject_serial(o%get_ori(i), imgs_pad(ithr))
            ! back FT
            call imgs_pad(ithr)%ifft()
            ! clip
            call imgs_pad(ithr)%clip(imgs(i))
            ! normalise
            call imgs(i)%norm()
        end do
        !$omp end parallel do
        ! destruct
        do ithr=1,p%nthr
            call imgs_pad(ithr)%kill
        end do
        deallocate(imgs_pad)
        call vol_pad%kill_expanded
        call vol_pad%kill
    end function project

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( vol, o, p, shvec ) result( rovol )
        class(image),   intent(inout) :: vol      !< volume to project
        class(ori),     intent(inout) :: o        !< orientation
        class(params),  intent(in)    :: p        !< parameters
        real, optional, intent(in)    :: shvec(3) !< 3D shift vector
        type(projector)  :: vol_pad
        type(image)      :: rovol_pad, rovol
        type(kbinterpol) :: kbwin
        integer          :: sh,h,k,l,nyq,lims(3,2),logi(3),phys(3),ldim(3),ldim_pd(3)
        real             :: loc(3)
        logical          :: l_shvec_present
        kbwin           = kbinterpol(KBWINSZ, p%alpha)
        ldim            = vol%get_ldim()
        ldim_pd         = [p%boxpd,p%boxpd,p%boxpd]
        l_shvec_present = present(shvec)
        call vol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%set_ft(.true.)
        call rovol%new(ldim, p%smpd)
        call prep4cgrid(vol, vol_pad, kbwin)
        call vol_pad%expand_cmat(p%alpha)
        lims = vol_pad%loop_lims(2)
        nyq  = vol_pad%get_lfny(1)
        write(*,'(A)') '>>> ROTATING VOLUME'
        !$omp parallel do collapse(3) default(shared) private(sh,h,k,l,loc,logi,phys)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    logi = [h,k,l]
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh > nyq + 1 )cycle
                    phys = rovol_pad%comp_addr_phys(logi)
                    loc  = matmul(real(logi), o%get_mat())
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
        call rovol%norm()
        call vol_pad%kill_expanded
        call vol_pad%kill
        call rovol_pad%kill
    end function rotvol

end module simple_projector_hlev
