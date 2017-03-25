module simple_projector_hlev
use simple_image,     only: image
use simple_oris,      only: oris
use simple_params,    only: params
use simple_gridding,  only: prep4cgrid
use simple_ori,       only: ori
use simple_math,      only: rotmat2d
use simple_projector, only: projector
use simple_jiffys     ! use all in there
implicit none

contains

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function projvol( vol, o, p, top ) result( imgs )   
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(in)    :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        type(image),      allocatable :: imgs(:)    !< resulting images
        character(len=:), allocatable :: imgk
        type(projector) :: vol_pad
        type(image)     :: img_pad
        integer         :: n, i, alloc_stat
        imgk = VOL%get_imgkind()
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, imgk)
        if( imgk .eq. 'xfel' )then
            call vol%pad(vol_pad)
        else
            call prep4cgrid(vol, vol_pad, p%msk)
        endif
        call img_pad%new([p%boxpd,p%boxpd,1], p%smpd, imgk)
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), stat=alloc_stat )
        call alloc_err('projvol; simple_projector', alloc_stat)
        write(*,'(A)') '>>> GENERATES PROJECTIONS' 
        do i=1,n
            call progress(i, n)
            call imgs(i)%new([p%box,p%box,1], p%smpd, imgk)
            call vol_pad%fproject( o%get_ori(i), img_pad )
            if( imgk .eq. 'xfel' )then
                call img_pad%clip(imgs(i))
            else
                call img_pad%bwd_ft
                call img_pad%clip(imgs(i))
                ! HAD TO TAKE OUT BECAUSE PGI COMPILER BAILS
                ! call imgs(i)%norm
            endif
        end do
        call vol_pad%kill
        call img_pad%kill
    end function projvol

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function projvol_expanded( vol, o, p, top, lp ) result( imgs )
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(inout) :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        real,    optional, intent(inout) :: lp      !< low-pass
        type(image),      allocatable :: imgs(:) !< resulting images
        character(len=:), allocatable :: imgk
        type(projector) :: vol_pad, img_pad
        integer         :: n, i, alloc_stat
        imgk = vol%get_imgkind()
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, imgk)
        if( imgk .eq. 'xfel' )then
            call vol%pad(vol_pad)
        else
            call prep4cgrid(vol, vol_pad, p%msk)
        endif
        call img_pad%new([p%boxpd,p%boxpd,1], p%smpd, imgk)
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), stat=alloc_stat )
        call vol_pad%expand_cmat
        call alloc_err('projvol; simple_projector', alloc_stat)
        write(*,'(A)') '>>> GENERATES PROJECTIONS' 
        do i=1,n
            call progress(i, n)
            call imgs(i)%new([p%box,p%box,1], p%smpd, imgk)
            if( present(lp) )then
                call vol_pad%fproject_expanded( o%get_ori(i), img_pad, lp=lp )
            else
                call vol_pad%fproject_expanded( o%get_ori(i), img_pad )
            endif
            if( imgk .eq. 'xfel' )then
                call img_pad%clip(imgs(i))
            else
                call img_pad%bwd_ft
                call img_pad%clip(imgs(i))
                ! HAD TO TAKE OUT BECAUSE PGI COMPILER BAILS
                ! call imgs(i)%norm
            endif
        end do
        call vol_pad%kill
        call img_pad%kill
    end function projvol_expanded

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( vol, o, p ) result( rovol )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(image),     intent(inout) :: vol  !< volume to project
        class(ori),       intent(inout) :: o    !< orientation
        class(params),    intent(in)    :: p    !< parameters
        type(projector) :: vol_pad
        type(image)     :: rovol_pad, rovol  
        integer         :: h,k,l,lims(3,2),logi(3),phys(3)
        real            :: loc(3)
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        rovol_pad = vol_pad
        call rovol_pad%set_ft(.true.)
        call rovol%new([p%box,p%box,p%box], p%smpd)
        call prep4cgrid(vol, vol_pad, p%msk)
        lims = vol_pad%loop_lims(2)
        write(*,'(A)') '>>> ROTATING VOLUME'
        !$omp parallel do default(shared) private(h,k,l,loc,logi,phys) schedule(auto)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    logi = [h,k,l]
                    phys = rovol_pad%comp_addr_phys([h,k,l])
                    loc  = matmul(real(logi),o%get_mat())
                    call rovol_pad%set_fcomp(logi, phys, vol_pad%extr_gridfcomp(loc))
                end do 
            end do
        end do
        !$omp end parallel do
        call rovol_pad%bwd_ft
        call rovol_pad%clip(rovol)
        ! HAD TO TAKE OUT BECAUSE PGI COMPILER BAILS
        ! call rovol%norm
        call vol_pad%kill
        call rovol_pad%kill
    end function rotvol

    !>  \brief  rotates an image by angle ang using Fourier gridding
    subroutine rotimg( img, ang, msk, roimg, shellw )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: hyp
        class(image),     intent(inout) :: img       !< image to rotate
        real,             intent(in)    :: ang       !< angle of rotation
        real,             intent(in)    :: msk       !< mask radius (in pixels)
        class(image),     intent(out)   :: roimg     !< rotated image
        real, optional,   intent(in)    :: shellw(:) !< shell-weights
        type(projector) :: img_pad 
        type(image)     :: roimg_pad, img_copy
        integer         :: h,k,lims(3,2),ldim(3),ldim_pd(3),logi(3),phys(3),sh,nyq
        real            :: loc(3),mat(2,2),smpd,fwght,wzero
        logical         :: doshellw
        doshellw   = present(shellw)
        ldim       = img%get_ldim()
        ldim_pd    = 2*ldim
        ldim_pd(3) = 1
        smpd       = img%get_smpd()
        call roimg%new(ldim, smpd)
        call img_pad%new(ldim_pd, smpd)
        call roimg_pad%new(ldim_pd, smpd)
        nyq        = img_pad%get_nyq()
        roimg_pad  = cmplx(0.,0.)
        img_copy   = img
        call img_copy%mask(msk, 'soft')
        call prep4cgrid(img, img_pad, msk)
        lims       = img_pad%loop_lims(2)
        mat        = rotmat2d(ang)
        wzero      = 1.0
        if( doshellw ) wzero = maxval(shellw)
        !$omp parallel do collapse(2) default(shared) private(h,k,loc,logi,phys,fwght,sh) schedule(auto)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)                
                loc(:2) = matmul(real([h,k]),mat)
                loc(3)  = 0.
                logi    = [h,k,0]
                phys    = img_pad%comp_addr_phys(logi)
                fwght   = 1.0
                if( doshellw )then
                    sh = nint(hyp(real(h),real(k)))
                    if( sh > nyq )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = shellw(sh)
                    endif
                endif
                call roimg_pad%set_fcomp(logi, phys, fwght * img_pad%extr_gridfcomp(loc))
            end do
        end do
        !$omp end parallel do
        call roimg_pad%bwd_ft
        call roimg_pad%clip(roimg)
        call img_pad%kill
        call roimg_pad%kill
        call img_copy%kill
    end subroutine rotimg

end module simple_projector_hlev
