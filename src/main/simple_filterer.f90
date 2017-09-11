! Fourier filtering routines
module simple_filterer
use simple_defs   ! use all in there
use simple_image, only: image
!$ use omp_lib
!$ use omp_lib_kinds
implicit none

contains

    !>  \brief generates an optimal low-pass filter (3D) from FRCs (2D) by 
    !!         (1) taking the average filter coefficients for frequencies 1 & 2
    !!         (2) minimizing the angle between the vector defined by the 3D Fourier index
    !!         and the planes used for calculating FRCs to find a matching filter coeff
    subroutine gen_anisotropic_optlp( vol_filter, projfrcs, e_space, state, pgrp )
        use simple_estimate_ssnr,   only: fsc2optlp
        use simple_projection_frcs, only: projection_frcs
        use simple_oris,            only: oris
        use simple_math,            only: hyp
        use simple_sym,             only: sym
        use simple_ori,             only: ori
        class(image),           intent(inout) :: vol_filter
        class(projection_frcs), intent(in)    :: projfrcs
        class(oris),            intent(inout) :: e_space
        integer,                intent(in)    :: state
        character(len=*),       intent(in)    :: pgrp
        type(sym)         :: se
        type(ori)         :: orientation, o_sym
        integer           :: noris, iori, nprojs, ldim(3), lims(3,2), isym, nsym
        integer           :: sh, logi(3), phys(3), imatch, filtsz, h, k, l
        real, allocatable :: plane_normals(:,:,:), plane_normals_L2(:,:), filters2D(:,:)
        real, allocatable :: frc(:)
        real              :: fwght, fwght_find0, fwght_find1, fwght_find2
        write(*,'(a)') '>>> GENERATING ANISOTROPIC OPTIMAL 3D LOW-PASS FILTER'
        ! sanity checking
        noris  = e_space%get_noris()
        nprojs = projfrcs%get_nprojs()
        if( noris /= nprojs )&
        &stop 'e_space & projfrcs objects non-conforming; filterer :: gen_anisotropic_optlp'
        ldim   = vol_filter%get_ldim()
        if( ldim(3) == 1 ) stop 'only for 3D filter generation; filterer :: gen_anisotropic_optlp'
        frc    = projfrcs%get_frc(1, ldim(1), state)
        filtsz = size(frc)
        lims   = vol_filter%loop_lims(2)
        if( pgrp .eq. 'c1' )then
            nsym = 1
            ! extract plane normals and L2 norms from e_space & extract 2D filters from projfrcs
            allocate( plane_normals(1,noris,3), plane_normals_L2(1,noris), filters2D(noris,filtsz) )
            plane_normals    = 0.0
            plane_normals_L2 = 0.0 
            filters2D        = 0.0 
            do iori=1,noris
                ! plane normals & L2 norms
                plane_normals(1,iori,:)  = e_space%get_normal(iori)
                plane_normals_L2(1,iori) = sum(plane_normals(1,iori,:) * plane_normals(1,iori,:))
                if( plane_normals_L2(1,iori) > TINY )then
                    plane_normals_L2(1,iori) = sqrt(plane_normals_L2(1,iori))
                endif
                ! 2D filters
                frc = projfrcs%get_frc(iori, ldim(1), state)
                filters2D(iori,:) = fsc2optlp(frc)
            end do
        else
            ! we need to expand over the symmetry group
            call se%new(pgrp)
            nsym = se%get_nsym()
            ! extract plane normals and L2 norms from e_space & extract 2D filters from projfrcs
            allocate( plane_normals(nsym,noris,3), plane_normals_L2(nsym,noris), filters2D(noris,filtsz) )
            plane_normals    = 0.0
            plane_normals_L2 = 0.0 
            filters2D        = 0.0 
            do iori=1,noris
                orientation = e_space%get_ori(iori)
                do isym=1,nsym
                    o_sym = se%apply(orientation, isym)
                    ! plane normals & L2 norms
                    plane_normals(isym,iori,:) = o_sym%get_normal()
                    plane_normals_L2(isym,iori) = sum(plane_normals(isym,iori,:) * plane_normals(isym,iori,:))
                    if( plane_normals_L2(isym,iori) > TINY )then
                        plane_normals_L2(isym,iori) = sqrt(plane_normals_L2(isym,iori))
                    endif
                end do
                ! 2D filters
                frc = projfrcs%get_frc(iori, ldim(1), state)
                filters2D(iori,:) = fsc2optlp(frc)
            end do
            call se%kill
        endif
        ! generate the 3D filter
        fwght_find0 = maxval(filters2D)
        fwght_find1 = sum(filters2D(:,1)) / real(noris)
        fwght_find2 = sum(filters2D(:,2)) / real(noris)
        vol_filter  = cmplx(0.0,0.0)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,logi,phys,imatch,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! shell
                    sh   = nint(hyp(real(h),real(k),real(l)))
                    ! logical index
                    logi = [h, k, l]
                    ! physical index
                    phys = vol_filter%comp_addr_phys(logi)
                    ! set 3D filter coeff
                    select case(sh)
                        case(0)
                            call vol_filter%set_fcomp(logi, phys, cmplx(fwght_find0,0.0))
                        case(1)
                            call vol_filter%set_fcomp(logi, phys, cmplx(fwght_find1,0.0))
                        case(2)
                            call vol_filter%set_fcomp(logi, phys, cmplx(fwght_find2,0.0))
                        case DEFAULT
                            if( sh <= filtsz )then
                                imatch = find2Dmatch(real([h,k,l]))
                                fwght  = filters2D(imatch,sh) ! filter coeff
                                call vol_filter%set_fcomp(logi, phys, cmplx(fwght,0.0))
                            endif
                    end select
                end do
            end do
        end do
        !$omp end parallel do

        contains

            ! index of matching 2D filter
            integer function find2Dmatch( vec )
                real, intent(in) :: vec(3)
                integer :: iori, loc(1)
                real    :: angles(noris)
                do iori=1,noris
                    angles(iori) = angle_btw_vec_and_normal( vec, iori )
                end do
                loc = minloc(angles)
                find2Dmatch = loc(1)
            end function find2Dmatch

            ! angle minimised over symmetry group
            real function angle_btw_vec_and_normal( vec, iori )
                real,    intent(in) :: vec(3)
                integer, intent(in) :: iori
                real    :: vec_L2, x, angle
                integer :: isym
                vec_L2 = sum(vec * vec)
                angle_btw_vec_and_normal = huge(x)
                if( vec_L2 > TINY )then
                    vec_L2 = sqrt(vec_L2)
                    do isym=1,nsym
                        if( plane_normals_L2(isym,iori) > TINY )then
                            angle = asin( abs(sum(plane_normals(isym,iori,:) * vec(:))) / (vec_L2 * plane_normals_L2(isym,iori)) )
                            if( angle < angle_btw_vec_and_normal ) angle_btw_vec_and_normal = angle
                        endif
                    end do
                endif
            end function angle_btw_vec_and_normal

    end subroutine gen_anisotropic_optlp

    !> \brief fits B-factor, untested
    subroutine fit_bfac( img_ref, img_ptcl, o, tfplan, bfac_range, lp, msk, bfac_best )
        use simple_ctf,   only: ctf
        use simple_ori,   only: ori
        class(image),  intent(in)    :: img_ref, img_ptcl
        class(ori),    intent(inout) :: o
        type(ctfplan), intent(in)    :: tfplan
        real,          intent(in)    :: bfac_range(2), lp, msk
        real,          intent(out)   :: bfac_best
        character(len=STDLEN) :: mode
        type(ctf)             :: tfun
        type(image)           :: ref, ptcl, maskimg, ref_tmp
        real                  :: dfx, dfy, angast, bfac_range_refine(2), smpd, corr_best, bfac, cc
        integer               :: npix, ldim(3), imode
        real, parameter       :: BFAC_STEPSZ = 5.0, BFAC_STEPSZ_REFINE = 1.0
        if( .not. (img_ref.eqdims.img_ptcl)) stop 'ref & ptcl imgs not of same dims; filterer :: fit_bfac'
        ! extract image info
        ldim = img_ref%get_ldim()
        smpd = img_ref%get_smpd()
        ! extract CTF info
        select case(tfplan%flag)
            case('yes')  ! multiply with CTF
                mode  ='ctf'
                imode = 1
            case('flip') ! multiply with abs(CTF)
                mode  = 'abs'
                imode = 2
            case('mul','no')
                mode  = ''
                imode = 3
        end select
        dfx    = o%get('dfx')
        dfy    = o%get('dfy')
        angast = o%get('angast')
        ! make hard mask
        call maskimg%disc(ldim, smpd, msk, npix)
        ! make CTF object
        tfun = ctf(smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
        ! prepare particle image
        ptcl = img_ptcl
        call ptcl%bp(0.,lp)
        call ptcl%bwd_ft
        ! prepare reference image
        ref = img_ref
        call ref%fwd_ft
        ! init
        corr_best = -1.
        bfac_best =  0.
        ! coarse search
        bfac = bfac_range(1)
        do while( bfac < bfac_range(2) )
            call iteration
            bfac = bfac + BFAC_STEPSZ
        end do
        ! refinement
        bfac_range_refine(1) = bfac_best - BFAC_STEPSZ
        bfac_range_refine(2) = bfac_best + BFAC_STEPSZ
        bfac = bfac_range_refine(1)
        do while( bfac < bfac_range_refine(2) )
            call iteration
             bfac = bfac + BFAC_STEPSZ_REFINE
        end do

        contains

            subroutine iteration
                ref_tmp = ref
                if( imode < 3 )then
                    call tfun%apply(ref_tmp, dfx, trim(mode), dfy, angast, bfac)
                else
                    call ref_tmp%apply_bfac(bfac)
                endif
                call ref_tmp%bp(0.,lp)
                call ref_tmp%bwd_ft
                cc = ref_tmp%real_corr(ptcl, maskimg)
                if( cc > corr_best )then
                    corr_best = cc
                    bfac_best = bfac
                endif
            end subroutine iteration

    end subroutine fit_bfac

end module simple_filterer
