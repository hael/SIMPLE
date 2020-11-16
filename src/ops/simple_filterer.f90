! Fourier filtering routines
module simple_filterer
!$ use omp_lib
!$ use omp_lib_kinds
use simple_math
use simple_image, only: image
use simple_oris,  only: oris
use simple_ori,   only: ori
use simple_error, only: simple_exception

implicit none

public :: gen_anisotropic_optlp
private
#include "simple_local_flags.inc"

contains

    !>  \brief generates an optimal low-pass filter (3D) from FRCs (2D) by
    !!         (1) taking the average filter coefficients for frequencies 1 & 2
    !!         (2) minimizing the angle between the vector defined by the 3D Fourier index
    !!         and the planes used for calculating FRCs to find a matching filter coeff
    subroutine gen_anisotropic_optlp( vol_filter, projfrcs, e_space, state, pgrp, hpind_fsc, phaseplate )
        use simple_sym,       only: sym
        use simple_estimate_ssnr,   only: fsc2optlp
        use simple_projection_frcs, only: projection_frcs
        class(image),           intent(inout) :: vol_filter
        class(projection_frcs), intent(in)    :: projfrcs
        class(oris),            intent(inout) :: e_space
        integer,                intent(in)    :: state
        character(len=*),       intent(in)    :: pgrp
        integer,                intent(in)    :: hpind_fsc
        logical,                intent(in)    :: phaseplate
        complex(kind=c_float_complex), pointer :: pcmat(:,:,:) => null()
        type(sym)         :: se
        type(ori)         :: orientation, o_sym
        integer           :: noris, iori, nprojs, ldim(3), lims(3,2), isym, nsym
        integer           :: sh, logi(3), phys(3), imatch, filtsz, h, k, l
        real, allocatable :: plane_normals(:,:,:), plane_normals_L2(:,:), filters2D(:,:)
        real, allocatable :: frc(:)
        real              :: loginorm(3),fwght, fwght_find0, fwght_find1, fwght_find2
        write(logfhandle,'(a)') '>>> GENERATING ANISOTROPIC OPTIMAL 3D LOW-PASS FILTER'
        ! sanity checking
        noris  = e_space%get_noris()
        nprojs = projfrcs%get_nprojs()
        if( noris /= nprojs ) THROW_HARD('e_space & projfrcs objects non-conforming; gen_anisotropic_optlp')
        ldim   = vol_filter%get_ldim()
        if( ldim(3) == 1 ) THROW_HARD('only for 3D filter generation; gen_anisotropic_optlp')
        filtsz = vol_filter%get_filtsz()
        allocate(frc(filtsz))
        lims   = vol_filter%loop_lims(2)
        if( pgrp .eq. 'c1' )then
            nsym = 1
            ! extract plane normals and L2 norms from e_space & extract 2D filters from projfrcs
            allocate( plane_normals(3,noris,1), plane_normals_L2(noris,1), filters2D(noris,filtsz) )
            plane_normals    = 0.0
            plane_normals_L2 = 0.0
            filters2D        = 0.0
            do iori=1,noris
                ! normalized plane normals & L2 norms
                plane_normals(:,iori,1)  = e_space%get_normal(iori)
                plane_normals_L2(iori,1) = sum(plane_normals(:,iori,1) * plane_normals(:,iori,1))
                if( plane_normals_L2(iori,1) > TINY )then
                    plane_normals_L2(iori,1) = sqrt(plane_normals_L2(iori,1))
                    plane_normals(:,iori,1)  = plane_normals(:,iori,1) / plane_normals_L2(iori,1)
                endif
                ! 2D filters
                call projfrcs%frc_getter(iori, hpind_fsc, phaseplate, frc, state)
                filters2D(iori,:) = fsc2optlp(frc)
            end do
        else
            ! we need to expand over the symmetry group
            call se%new(pgrp)
            nsym = se%get_nsym()
            ! extract plane normals and L2 norms from e_space & extract 2D filters from projfrcs
            allocate( plane_normals(3,noris,nsym), plane_normals_L2(noris,nsym), filters2D(noris,filtsz) )
            plane_normals    = 0.0
            plane_normals_L2 = 0.0
            filters2D        = 0.0
            do iori=1,noris
                call e_space%get_ori(iori, orientation)
                do isym=1,nsym
                    call se%apply(orientation, isym, o_sym)
                    ! plane normals & L2 norms
                    plane_normals(:,iori,isym)  = o_sym%get_normal()
                    plane_normals_L2(iori,isym) = sum(plane_normals(:,iori,isym) * plane_normals(:,iori,isym))
                    if( plane_normals_L2(iori,isym) > TINY )then
                        plane_normals_L2(iori,isym) = sqrt(plane_normals_L2(iori,isym))
                        plane_normals(:,iori,isym)  = plane_normals(:,iori,isym) / plane_normals_L2(iori,isym)
                    endif
                end do
                ! 2D filters
                call projfrcs%frc_getter(iori, hpind_fsc, phaseplate, frc, state)
                filters2D(iori,:) = fsc2optlp(frc)
            end do
            call se%kill
        endif
        ! generate the 3D filter
        fwght_find0 = maxval(filters2D)
        fwght_find1 = sum(filters2D(:,1)) / real(noris)
        fwght_find2 = sum(filters2D(:,2)) / real(noris)
        call vol_filter%zero_and_flag_ft
        call vol_filter%get_cmat_ptr(pcmat)
        ! dynamic scheduling works better here because of the imbalanced amount of compute
        ! it cuts execution time by half as compared to static scheduling
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,logi,phys,imatch,fwght,loginorm)&
        !$omp schedule(dynamic) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! shell
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh <= filtsz )then
                        ! logical index
                        logi = [h, k, l]
                        ! physical index
                        phys = vol_filter%comp_addr_phys(logi)
                        ! set 3D filter coeff
                        select case(sh)
                            case(0)
                                pcmat(phys(1),phys(2),phys(3)) = cmplx(fwght_find0)
                            case(1)
                                pcmat(phys(1),phys(2),phys(3)) = cmplx(fwght_find1)
                            case(2)
                                pcmat(phys(1),phys(2),phys(3)) = cmplx(fwght_find2)
                            case DEFAULT
                                loginorm = real(logi) / sum(real(logi*logi)) ! normalize
                                imatch   = find2Dmatch(loginorm)
                                fwght    = filters2D(imatch,sh) ! filter coeff
                                pcmat(phys(1),phys(2),phys(3)) = cmplx(fwght)
                        end select
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        nullify(pcmat)
        call orientation%kill
        call o_sym%kill
        contains

            ! index of matching 2D filter
            integer function find2Dmatch( vec )
                real, intent(in) :: vec(3)
                real    :: sines(noris)
                integer :: iori
                do iori=1,noris
                    sines(iori) = sine_btw_vec_and_normal(vec, iori)
                end do
                find2Dmatch = minloc(sines,dim=1)
            end function find2Dmatch

            ! sine minimised over symmetry group (= angle minimized)
            real function sine_btw_vec_and_normal( vec, iori )
                real,    intent(in) :: vec(3)
                integer, intent(in) :: iori
                real    :: sine
                integer :: isym
                sine_btw_vec_and_normal = huge(sine)
                do isym=1,nsym
                    if( plane_normals_L2(iori,isym) > TINY )then
                        sine = abs( sum(plane_normals(:,iori,isym) * vec) )
                        sine_btw_vec_and_normal = min(sine, sine_btw_vec_and_normal)
                    endif
                end do
            end function sine_btw_vec_and_normal

    end subroutine gen_anisotropic_optlp

end module simple_filterer
