module simple_projector_hlev
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,      only: image
use simple_oris,       only: oris
use simple_params,     only: params
use simple_gridding,   only: prep4cgrid
use simple_ori,        only: ori
use simple_math,       only: rotmat2d
use simple_projector,  only: projector
use simple_kbinterpol, only: kbinterpol
use simple_jiffys      ! use all in there
implicit none

contains

    !>  \brief  generates an array of projection images of volume vol in orientations o
    function projvol( vol, o, p, top, lp ) result( imgs )
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(inout) :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        real,    optional, intent(inout) :: lp      !< low-pass
        type(image),      allocatable :: imgs(:) !< resulting images
        character(len=:), allocatable :: imgk
        type(projector)  :: vol_pad, img_pad
        type(kbinterpol) :: kbwin
        integer          :: n, i, alloc_stat
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        imgk  = vol%get_imgkind()
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, imgk)
        if( imgk .eq. 'xfel' )then
            call vol%pad(vol_pad)
        else
            call prep4cgrid(vol, vol_pad, p%msk, kbwin)
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
                call vol_pad%fproject( o%get_ori(i), img_pad, lp=lp )
            else
                call vol_pad%fproject( o%get_ori(i), img_pad )
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
    end function projvol

    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( vol, o, p ) result( rovol )
        class(image),     intent(inout) :: vol  !< volume to project
        class(ori),       intent(inout) :: o    !< orientation
        class(params),    intent(in)    :: p    !< parameters
        type(projector)  :: vol_pad
        type(image)      :: rovol_pad, rovol 
        type(kbinterpol) :: kbwin
        integer          :: h,k,l,lims(3,2),logi(3),phys(3),ldim(3),ldim_pd(3)
        real             :: loc(3)
        kbwin   = kbinterpol(KBWINSZ, KBALPHA)
        ldim    = vol%get_ldim()
        ldim_pd = nint(KBALPHA)*ldim
        call vol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%new(ldim_pd, p%smpd)
        call rovol_pad%set_ft(.true.)
        call rovol%new(ldim, p%smpd)
        call prep4cgrid(vol, vol_pad, p%msk, kbwin)
        ! call vol_pad%expand_cmat REMOVED BECAUSE SEE BELOW
        lims = vol_pad%loop_lims(2)
        write(*,'(A)') '>>> ROTATING VOLUME'
        !$omp parallel do collapse(3) default(shared) private(h,k,l,loc,logi,phys)&
        !$omp schedule(static) proc_bind(close)
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
        use simple_math, only: hyp
        class(image),     intent(inout) :: img       !< image to rotate
        real,             intent(in)    :: ang       !< angle of rotation
        real,             intent(in)    :: msk       !< mask radius (in pixels)
        class(image),     intent(out)   :: roimg     !< rotated image
        real, optional,   intent(in)    :: shellw(:) !< shell-weights
        type(projector)  :: img_pad 
        type(image)      :: roimg_pad
        type(kbinterpol) :: kbwin
        integer          :: h,k,lims(3,2),ldim(3),ldim_pd(3),logi(3),phys(3),sh,nyq
        real             :: loc(3),mat(2,2),smpd,fwght,wzero
        logical          :: doshellw
        kbwin      = kbinterpol(KBWINSZ, KBALPHA)
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
        call prep4cgrid(img, img_pad, msk, kbwin)
        lims       = img_pad%loop_lims(2)
        mat        = rotmat2d(ang)
        wzero      = 1.0
        if( doshellw ) wzero = maxval(shellw)
        !$omp parallel do collapse(2) default(shared) private(h,k,loc,logi,phys,fwght,sh)&
        !$omp schedule(static) proc_bind(close)
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
    end subroutine rotimg

    !>  \brief  rotates an image batch by angle ang using Fourier gridding
    !>          the weighted images sum is returned
    subroutine rot_imgbatch( imgs, os, imgsum, msk )
        use simple_math, only: cyci_1d, sqwin_2d, rotmat2d
        class(image), intent(inout) :: imgs(:)   !< images to rotate
        type(oris),   intent(inout) :: os        !< orientations
        class(image), intent(inout) :: imgsum    !< reconstituted image
        real,         intent(in)    :: msk       !< mask radius (in pixels)
        type(projector), allocatable :: padded_imgs(:)
        complex,         allocatable :: cmat(:,:), comps(:,:)
        real,            allocatable :: w(:,:)
        integer,         allocatable :: cyc1(:), cyc2(:)
        type(kbinterpol) :: kbwin
        complex          :: comp, zero
        integer          :: lims(3,2), ldim(3), ldim_pd(3), logi(3), phys(3), win(2,2), cyc_lims(3,2)
        integer          :: alloc_stat, wdim, incr, nimgs, i, h, k,l,m
        real             :: loc(2), mat(2,2), smpd, winsz, pw
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        ! init
        zero       = cmplx(0.,0.)
        winsz      = KBWINSZ
        ldim       = imgs(1)%get_ldim()
        ldim_pd    = nint(KBALPHA)*ldim
        ldim_pd(3) = 1
        smpd       = imgs(1)%get_smpd()
        nimgs      = size(imgs)
        wdim       = ceiling(KBALPHA*winsz) + 1
        ! gridding
        allocate(padded_imgs(nimgs))
        do i = 1, nimgs
            call padded_imgs(i)%new(ldim_pd, smpd)
            call prep4cgrid(imgs(i), padded_imgs(i), msk, kbwin)
        enddo
        ! rotation
        lims     = padded_imgs(1)%loop_lims(2)
        cyc_lims = padded_imgs(1)%loop_lims(3)
        allocate(cmat(lims(1,1):lims(1,2), lims(2,1):lims(2,2)), cyc1(wdim), cyc2(wdim),&
        &w(wdim, wdim), comps(wdim, wdim), stat=alloc_stat)
        cmat = zero
        !$omp parallel do default(shared) private(i,h,k,l,m,loc,mat,logi,phys,cyc1,cyc2,w,comps,win,incr,pw)&
        !$omp schedule(static) reduction(+:cmat) proc_bind(close)
        do i = 1, nimgs
            pw  = os%get(i,'w')
            if( pw > TINY )then
                mat = rotmat2d( -os%e3get(i) )
                do h = lims(1,1), lims(1,2)
                    do k = lims(2,1), lims(2,2)
                        loc   = matmul(real([h,k]),mat)
                        win   = sqwin_2d(loc(1),loc(2), winsz)
                        comps = zero
                        w     = 1.
                        do l = 1, wdim
                            incr = l-1
                            ! circular addresses
                            cyc1(l)  = cyci_1d(cyc_lims(1,:), win(1,1)+incr)
                            cyc2(l)  = cyci_1d(cyc_lims(2,:), win(2,1)+incr)
                            ! interpolation kernel matrix
                            w(l,:) = w(l,:) * kbwin%apod( real(win(1,1)+incr)-loc(1) )
                            w(:,l) = w(:,l) * kbwin%apod( real(win(2,1)+incr)-loc(2) )
                        enddo
                        ! fetch fourier components
                        do l=1, wdim
                            do m=1, wdim
                                if(w(l,m) == 0.)cycle
                                logi       = [cyc1(l), cyc2(m), 0]
                                phys       = padded_imgs(i)%comp_addr_phys(logi)
                                comps(l,m) = padded_imgs(i)%get_fcomp(logi, phys)
                            end do
                        end do
                        ! SUM( kernel x components )
                        cmat(h,k) = cmat(h,k) + pw * sum(w * comps)
                        ! above is an optimized version of:
                        ! cmat(h,k) = cmat(h,k) + pw * padded_imgs(i)%extr_gridfcomp( [loc(1),loc(2),0.] )
                    end do
                end do
            endif
            ! cleanup
            call padded_imgs(i)%kill
        enddo
        !$omp end parallel do
        ! transfer to image
        call imgsum%new(ldim_pd, smpd)
        call imgsum%set_ft(.true.)
        imgsum = zero
        !$omp parallel do collapse(2) default(shared) private(h,k,logi,phys,comp)&
        !$omp schedule(static) proc_bind(close)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                comp = cmat(h, k)
                if( comp.eq.zero )cycle
                logi = [h, k, 0]
                phys = imgsum%comp_addr_phys(logi)
                call imgsum%set_fcomp(logi, phys, comp)
            end do 
        end do
        !$omp end parallel do
        ! real space & clipping
        call imgsum%bwd_ft
        call imgsum%clip_inplace(ldim)
        ! cleanup
        deallocate(padded_imgs,comps,w,cyc1,cyc2,cmat)
    end subroutine rot_imgbatch

end module simple_projector_hlev
