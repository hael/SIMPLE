! projection of 3D volumes in the Fourier domain by convolution interpolation
! to generate band-pass limited Cartesian and polar 2D Fourier transforms,
! high-level routines
module simple_projector_hlev
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_syslib
use simple_image,      only: image
use simple_oris,       only: oris
use simple_params,     only: params
use simple_gridding,   only: prep4cgrid
use simple_ori,        only: ori
use simple_math,       only: rotmat2d
use simple_projector,  only: projector
use simple_kbinterpol, only: kbinterpol
use simple_jiffys,     only: progress ! use all in there
implicit none

contains

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function projvol( vol, o, p, top, lp ) result( imgs )
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(inout) :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        real,    optional, intent(inout) :: lp      !< low-pass
        type(image),       allocatable :: imgs(:)   !< resulting images
        type(projector)  :: vol_pad, img_pad
        type(kbinterpol) :: kbwin
        integer          :: n, i, alloc_stat
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call prep4cgrid(vol, vol_pad, p%msk, kbwin)
        call img_pad%new([p%boxpd,p%boxpd,1], p%smpd)
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), stat=alloc_stat )
        call alloc_errchk('projvol; simple_projector', alloc_stat)
        call vol_pad%expand_cmat
        write(*,'(A)') '>>> GENERATES PROJECTIONS' 
        do i=1,n
            call progress(i, n)
            call imgs(i)%new([p%box,p%box,1], p%smpd)
            if( present(lp) )then
                call vol_pad%fproject( o%get_ori(i), img_pad, lp=lp )
            else
                call vol_pad%fproject( o%get_ori(i), img_pad )
            endif
            call img_pad%bwd_ft
            call img_pad%clip(imgs(i))
        end do
        call vol_pad%kill_expanded
        call vol_pad%kill
        call img_pad%kill
    end function projvol

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( vol, o, p, shvec ) result( rovol )
        class(image),   intent(inout) :: vol      !< volume to project
        class(ori),     intent(inout) :: o        !< orientation
        class(params),  intent(in)    :: p        !< parameters
        real, optional, intent(in)    :: shvec(3) !< 3D shift vector
        type(projector)  :: vol_pad
        type(image)      :: rovol_pad, rovol 
        type(kbinterpol) :: kbwin
        integer          :: h,k,l,lims(3,2),logi(3),phys(3),ldim(3),ldim_pd(3)
        real             :: loc(3)
        logical          :: l_shvec_present
        kbwin           = kbinterpol(KBWINSZ, KBALPHA)
        ldim            = vol%get_ldim()
        ldim_pd         = nint(KBALPHA)*ldim
        l_shvec_present = present(shvec)
        call vol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%set_ft(.true.)
        call rovol%new(ldim, p%smpd)
        call prep4cgrid(vol, vol_pad, p%msk, kbwin)
        lims = vol_pad%loop_lims(2)
        write(*,'(A)') '>>> ROTATING VOLUME'
        !$omp parallel do collapse(3) default(shared) private(h,k,l,loc,logi,phys)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    logi = [h,k,l]
                    phys = rovol_pad%comp_addr_phys(logi)
                    loc  = matmul(real(logi), o%get_mat())
                    if( l_shvec_present )then
                        call rovol_pad%set_fcomp(logi, phys, vol_pad%extr_gridfcomp(loc) * rovol_pad%oshift(loc, shvec))
                    else
                        call rovol_pad%set_fcomp(logi, phys, vol_pad%extr_gridfcomp(loc))
                    endif
                end do 
            end do
        end do
        !$omp end parallel do
        call rovol_pad%bwd_ft
        call rovol_pad%clip(rovol)
        call rovol%norm()
        call vol_pad%kill
        call rovol_pad%kill
    end function rotvol

    !>  \brief  rotates an image by angle ang using Fourier gridding
    subroutine rotimg( img, ang, msk, roimg )
        use simple_math, only: hyp
        class(image),     intent(inout) :: img   !< image to rotate
        real,             intent(in)    :: ang   !< angle of rotation
        real,             intent(in)    :: msk   !< mask radius (in pixels)
        class(image),     intent(out)   :: roimg !< rotated image
        type(projector)  :: img_pad 
        type(image)      :: roimg_pad
        type(kbinterpol) :: kbwin
        integer          :: h,k,lims(3,2),ldim(3),ldim_pd(3),logi(3),phys(3),sh,nyq
        real             :: loc(3),mat(2,2),smpd
        kbwin      = kbinterpol(KBWINSZ, KBALPHA)
        ldim       = img%get_ldim()
        ldim_pd    = 2*ldim
        ldim_pd(3) = 1
        smpd       = img%get_smpd()
        call roimg%new(ldim, smpd)
        call img_pad%new(ldim_pd, smpd)
        call roimg_pad%new(ldim_pd, smpd)
        nyq        = img_pad%get_nyq()
        roimg_pad  = cmplx(0.,0.)
        call prep4cgrid(img, img_pad, msk, kbwin)
        lims       = img_pad%loop_lims(2)
        mat        = rotmat2d(ang) 
        !$omp parallel do collapse(2) default(shared) private(h,k,loc,logi,phys,sh)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)                
                loc(:2) = matmul(real([h,k]),mat)
                loc(3)  = 0.
                logi    = [h,k,0]
                phys    = img_pad%comp_addr_phys(logi)
                call roimg_pad%set_fcomp(logi, phys, img_pad%extr_gridfcomp(loc))
            end do
        end do
        !$omp end parallel do
        call roimg_pad%bwd_ft
        call roimg_pad%clip(roimg)
        call img_pad%kill
        call roimg_pad%kill
    end subroutine rotimg

end module simple_projector_hlev
