module simple_polarops
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_builder,           only: builder, build_glob
use simple_parameters,        only: params_glob
use simple_sp_project,        only: sp_project
use simple_image,             only: image
use simple_polarft_corrcalc,  only: polarft_corrcalc, pftcc_glob
use simple_strategy2D_utils
implicit none

! Restoration
public :: polar_cavger_new, polar_cavger_dims_from_header, polar_cavger_update_sums
public :: polar_cavger_calc_and_write_frcs_and_eoavg, polar_cavger_merge_eos_and_norm
public :: polar_cavger_assemble_sums_from_parts, polar_cavger_set_ref_pftcc
public :: polar_cavger_zero_pft_refs, polar_cavger_kill
! I/O
public :: polar_cavger_refs2cartesian, polar_cavger_write_cartrefs, polar_cavger_writeall_pftccrefs
public :: polar_cavger_write, polar_cavger_writeall, polar_cavger_writeall_cartrefs
public :: polar_cavger_readwrite_partial_sums, polar_cavger_read_all
! Book-keeping
public :: polar_cavger_gen2Dclassdoc, polar_cavger_calc_pops
! Alignment
public :: polar_prep2Dref, polar_prep3Dref
private
#include "simple_local_flags.inc"


complex(dp), allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:) ! PFTs arrays
real(dp),    allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:)                   ! PFT-size CTF2 arrays
integer,     allocatable :: prev_eo_pops(:,:), eo_pops(:,:)                     ! Class populations
real                     :: smpd       = 0.                                     ! Pixel size
integer                  :: ncls       = 0                                      ! # classes
integer                  :: kfromto(2) = 0                                      ! Resolution range
integer                  :: pftsz      = 0                                      ! Size of PFT in pftcc along rotation dimension
integer                  :: nrots      = 0                                      ! # of in-plane rotations; =2*pftsz
logical                  :: l_comlin   = .false.                                ! Whether we operate in 2D/3D(=added common lines)
contains

    !> Module initialization
    subroutine polar_cavger_new( pftcc, comlin, nrefs )
        class(polarft_corrcalc), intent(in) :: pftcc
        logical,                 intent(in) :: comlin
        integer,       optional, intent(in) :: nrefs
        call polar_cavger_kill
        nrots   = pftcc%get_nrots()
        pftsz   = pftcc%get_pftsz()
        kfromto = pftcc%get_kfromto()
        if( present(nrefs) )then
            ncls = nrefs
        else
            ncls = pftcc%get_nrefs()
        endif
        l_comlin = comlin
        ! dimensions
        smpd    = params_glob%smpd
        allocate(prev_eo_pops(2,ncls), eo_pops(2,ncls), source=0)
        ! Arrays        
        allocate(pfts_even(pftsz,kfromto(1):kfromto(2),ncls),pfts_odd(pftsz,kfromto(1):kfromto(2),ncls),&
                &ctf2_even(pftsz,kfromto(1):kfromto(2),ncls),ctf2_odd(pftsz,kfromto(1):kfromto(2),ncls),&
                &pfts_merg(pftsz,kfromto(1):kfromto(2),ncls))
        call polar_cavger_zero_pft_refs
        pfts_merg = DCMPLX_ZERO
    end subroutine polar_cavger_new

    subroutine polar_cavger_zero_pft_refs
        pfts_even = DCMPLX_ZERO
        pfts_odd  = DCMPLX_ZERO
        ctf2_even = 0.d0
        ctf2_odd  = 0.d0
    end subroutine polar_cavger_zero_pft_refs

    subroutine polar_cavger_set_ref_pftcc( icls, which, pftcc )
        integer,                 intent(in)    :: icls
        character(len=*),        intent(in)    :: which
        class(polarft_corrcalc), intent(inout) :: pftcc
        select case(trim(which))
        case('merged')
            call pftcc%set_ref_pft(icls, cmplx(pfts_merg(:,:,icls),kind=sp), .true.)
        case('even')
            call pftcc%set_ref_pft(icls, cmplx(pfts_even(:,:,icls),kind=sp), .true.)
        case('odd')
            call pftcc%set_ref_pft(icls, cmplx(pfts_odd(:,:,icls),kind=sp), .false.)
        end select
    end subroutine polar_cavger_set_ref_pftcc

    subroutine polar_cavger_calc_pops( spproj )
        class(sp_project), target, intent(in) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer :: i, icls, iptcl, eo
        logical :: l_3D
        l_3D = .false.
        select case(trim(params_glob%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(params_glob%oritype))
        end select
        eo_pops = 0
        !$omp parallel do schedule(guided) proc_bind(close) default(shared)&
        !$omp private(iptcl,eo,icls)&
        !$omp reduction(+:eo_pops)
        do iptcl = 1,ptcl_field%get_noris()
            if( ptcl_field%get_state(iptcl) == 0  ) cycle
            if( ptcl_field%get(iptcl,'w') < SMALL ) cycle
            eo = ptcl_field%get_eo(iptcl)+1
            if( l_3D )then
                icls = ptcl_field%get_proj(iptcl)
            else
                icls = ptcl_field%get_class(iptcl)
            endif
            eo_pops(eo,icls) = eo_pops(eo,icls) + 1
        enddo
        !$omp end parallel do
        prev_eo_pops = 0
        if( cls_field%get_noris() == ncls )then
            do i = 1,ncls
                if( l_3D )then
                    icls = ptcl_field%get_proj(i)
                else
                    icls = ptcl_field%get_class(i)
                endif
                if( .not.cls_field%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(1,icls) = cls_field%get_int(i,'prev_pop_even')
                prev_eo_pops(2,icls) = cls_field%get_int(i,'prev_pop_odd')
            enddo
        endif
        eo_pops = eo_pops + prev_eo_pops
    end subroutine polar_cavger_calc_pops

    !>  \brief  Updates Fourier components and normalization matrices with new particles
    subroutine polar_cavger_update_sums( nptcls, pinds, spproj, pftcc, incr_shifts, is3D )
        integer,                         intent(in)    :: nptcls
        integer,                         intent(in)    :: pinds(nptcls)
        class(sp_project),               intent(inout) :: spproj
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        real,                  optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,               optional, intent(in)    :: is3d
        class(oris), pointer :: spproj_field
        complex(sp), pointer :: pptcls(:,:,:), rptcl(:,:)
        real(sp),    pointer :: pctfmats(:,:,:), rctf(:,:)
        real(dp) :: w
        real     :: incr_shift(2)
        integer  :: eopops(2,ncls), i, icls, iptcl, irot
        logical  :: l_ctf, l_even, l_3D, l_shift
        l_3D = .false.
        if( present(is3D) ) l_3D = is3D
        l_shift = .false.
        if( present(incr_shifts) ) l_shift = .true.
        ! retrieve particle info & pointers
        call spproj%ptr2oritype(params_glob%oritype, spproj_field)
        call pftcc%get_ptcls_ptr(pptcls)
        l_ctf = pftcc%is_with_ctf()
        if( l_ctf ) call pftcc%get_ctfmats_ptr(pctfmats)
        ! update classes
        eopops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static) reduction(+:eopops)&
        !$omp private(i,iptcl,w,l_even,icls,irot,incr_shift,rptcl,rctf)
        do i = 1,nptcls
            ! particles parameters
            iptcl = pinds(i)
            if( spproj_field%get_state(iptcl) == 0  ) cycle
            w = real(spproj_field%get(iptcl,'w'),dp)
            if( w < DSMALL ) cycle
            l_even = spproj_field%get_eo(iptcl)==0
            if( l_3D )then
                icls = spproj_field%get_proj(iptcl)
            else
                icls = spproj_field%get_class(iptcl)
            endif
            irot = pftcc%get_roind_fast(spproj_field%e3get(iptcl))
            if( l_shift )then
                incr_shift = incr_shifts(:,i)
                ! weighted restoration
                if( any(abs(incr_shift) > 1.e-6) ) call pftcc%shift_ptcl(iptcl, -incr_shift)
            endif
            call pftcc%get_work_pft_ptr(rptcl)
            call pftcc%rotate_pft(pptcls(:,:,i), irot, rptcl)
            if( l_ctf )then
                call pftcc%get_work_rpft_ptr(rctf)
                call pftcc%rotate_pft(pctfmats(:,:,i), irot, rctf)
                if( l_even )then
                    !$omp critical
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + w * cmplx(rptcl,kind=dp) * real(rctf,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w * real(rctf,kind=dp)**2
                    !$omp end critical
                else
                    !$omp critical
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + w * cmplx(rptcl,kind=dp) * real(rctf,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w * real(rctf,kind=dp)**2
                    !$omp end critical
                endif
            else
                if( l_even )then
                    !$omp critical
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + w * cmplx(rptcl,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w
                    !$omp end critical
                else
                    !$omp critical
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + w * cmplx(rptcl,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w
                    !$omp end critical
                endif
            endif
            ! total population
            if( l_even )then
                eopops(1,icls) = eopops(1,icls) + 1
            else
                eopops(2,icls) = eopops(2,icls) + 1
            endif
        enddo
        !$omp end parallel do
        eo_pops = eo_pops + eopops
        ! cleanup
        nullify(spproj_field,rptcl,rctf,pptcls,pctfmats)
    end subroutine polar_cavger_update_sums

    !>  \brief  Restores class-averages
    subroutine polar_cavger_merge_eos_and_norm( reforis )
        type(oris), optional, intent(in) :: reforis
        real,          parameter :: EPSILON = 0.1
        logical,       parameter :: l_kb = .true.
        complex(dp), allocatable :: pfts_cavg(:,:,:), pfts_clin(:,:,:), clin_odd(:,:,:), clin_even(:,:,:)
        ! real,        allocatable :: res(:)
        complex(dp) :: pfts_clin_even(pftsz,kfromto(1):kfromto(2),ncls),pfts_clin_odd(pftsz,kfromto(1):kfromto(2),ncls),&
                      &pft(pftsz,kfromto(1):kfromto(2)),pfte(pftsz,kfromto(1):kfromto(2)),pfto(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)),ctf2e(pftsz,kfromto(1):kfromto(2)),ctf2o(pftsz,kfromto(1):kfromto(2)),&
                      &ctf2_clin_even(pftsz,kfromto(1):kfromto(2),ncls), ctf2_clin_odd(pftsz,kfromto(1):kfromto(2),ncls)
        integer     :: icls, eo_pop(2), pop, k, pops(ncls), npops, m
        ! real        :: res_fsc05, res_fsc0143, min_res_fsc0143, max_res_fsc0143, avg_res_fsc0143, avg_res_fsc05
        real        :: cavg_clin_frcs(pftsz,kfromto(1):kfromto(2),ncls), psi!, dfrcs(kfromto(1):kfromto(2),ncls)
        logical     :: l_rotm
        if( l_comlin )then
            ! 3D-related tasks
            if( .not. present(reforis) )THROW_HARD('Reference orientations need be inputted in polar_cavger_merge_eos_and_norm')
            ! Mirroring etc.
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! Common-lines conribution
            call calc_comlin_contrib(reforis, build_glob%pgrpsyms)
            ! need to compute the cavg/clin frcs prior to their summing
            if( trim(params_glob%polar_frcs) .eq. 'yes' )then
                allocate(pfts_cavg(pftsz,kfromto(1):kfromto(2),ncls),&
                        &pfts_clin(pftsz,kfromto(1):kfromto(2),ncls),&
                        &clin_even(pftsz,kfromto(1):kfromto(2),ncls),&
                        &clin_odd( pftsz,kfromto(1):kfromto(2),ncls),source=DCMPLX_ZERO)
                cavg_clin_frcs = 1.
                call get_cavg_clin
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,k)
                do icls = 1, ncls
                    if( pops(icls) < 2 )cycle
                    do k = kfromto(1), kfromto(2)
                        call nonuni_frc(pfts_clin(:,k,icls), pfts_cavg(:,k,icls), cavg_clin_frcs(:,k,icls))
                    enddo
                enddo
                !$omp end parallel do
                deallocate(pfts_cavg,pfts_clin,clin_odd,clin_even)
            endif
        endif
        pfts_merg = DCMPLX_ZERO
        select case(trim(params_glob%ref_type))
            case('cavg')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,eo_pop,pop,pft,ctf2)
                do icls = 1,ncls
                    eo_pop = prev_eo_pops(:,icls) + eo_pops(:,icls) ! eo_pops has to be calculated differently
                    pop    = sum(eo_pop)
                    if(pop == 0)then
                        pfts_even(:,:,icls) = DCMPLX_ZERO
                        pfts_odd(:,:,icls)  = DCMPLX_ZERO
                        ctf2_even(:,:,icls) = 0.d0
                        ctf2_odd(:,:,icls)  = 0.d0
                    else
                        if(pop > 1)then
                            pft  = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                            ctf2 = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                            if( pop <= 5 ) ctf2 = ctf2 + real(EPSILON/real(pop),dp)
                            where( ctf2 > DSMALL ) pfts_merg(:,:,icls) = pft / ctf2
                        endif
                        if(eo_pop(1) > 1)then
                            where( ctf2_even(:,:,icls) > DSMALL ) pfts_even(:,:,icls) = pfts_even(:,:,icls) / ctf2_even(:,:,icls)
                        endif
                        if(eo_pop(2) > 1)then
                            where( ctf2_odd(:,:,icls) > DSMALL )  pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  / ctf2_odd(:,:,icls)
                        endif
                    endif
                end do
                !$omp end parallel do
            case('clin')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
                do icls = 1,ncls/2
                    ! already mirrored common-line contribution
                    pfte  = pfts_clin_even(:,:,icls)
                    pfto  = pfts_clin_odd(:,:,icls)
                    ctf2e = ctf2_clin_even(:,:,icls)
                    ctf2o = ctf2_clin_odd(:,:,icls)
                    ! merged then e/o
                    pft   = pfte  + pfto
                    ctf2  = ctf2e + ctf2o
                    call safe_norm(pft,  ctf2,  pfts_merg(:,:,icls))
                    call safe_norm(pfte, ctf2e, pfts_even(:,:,icls))
                    call safe_norm(pfto, ctf2o, pfts_odd(:,:,icls))
                    ! mirroring the restored images
                    m = reforis%get_int(icls,'mirr')
                    call mirror_pft(pfts_merg(:,:,icls), pfts_merg(:,:,m))
                    call mirror_pft(pfts_even(:,:,icls), pfts_even(:,:,m))
                    call mirror_pft(pfts_odd(:,:,icls),  pfts_odd(:,:,m))
                    psi    = abs(reforis%get(m, 'psi'))
                    l_rotm = (psi > 0.1) .and. (psi < 359.9)
                    if( l_rotm )then
                        pfts_merg(:,:,m) = conjg(pfts_merg(:,:,m))
                        pfts_even(:,:,m) = conjg(pfts_even(:,:,m))
                        pfts_odd(:,:,m)  = conjg(pfts_odd(:,:,m))
                    endif
                enddo
                !$omp end parallel do
            case('vol')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
                do icls = 1,ncls/2
                    ! calculate sum of already mirrored class + common-line contribution
                    pfte  = pfts_even(:,:,icls) + pfts_clin_even(:,:,icls)
                    pfto  = pfts_odd(:,:,icls)  + pfts_clin_odd(:,:,icls)
                    ctf2e = ctf2_even(:,:,icls) + ctf2_clin_even(:,:,icls)
                    ctf2o = ctf2_odd(:,:,icls)  + ctf2_clin_odd(:,:,icls)
                    ! merged then e/o
                    pft   = pfte  + pfto
                    ctf2  = ctf2e + ctf2o
                    call safe_norm(pft,  ctf2,  pfts_merg(:,:,icls))
                    call safe_norm(pfte, ctf2e, pfts_even(:,:,icls))
                    call safe_norm(pfto, ctf2o, pfts_odd(:,:,icls))
                    ! mirroring the restored images
                    m = reforis%get_int(icls,'mirr')
                    call mirror_pft(pfts_merg(:,:,icls), pfts_merg(:,:,m))
                    call mirror_pft(pfts_even(:,:,icls), pfts_even(:,:,m))
                    call mirror_pft(pfts_odd(:,:,icls),  pfts_odd(:,:,m))
                    psi    = abs(reforis%get(m, 'psi'))
                    l_rotm = (psi > 0.1) .and. (psi < 359.9)
                    if( l_rotm )then
                        pfts_merg(:,:,m) = conjg(pfts_merg(:,:,m))
                        pfts_even(:,:,m) = conjg(pfts_even(:,:,m))
                        pfts_odd(:,:,m)  = conjg(pfts_odd(:,:,m))
                    endif
                enddo
                !$omp end parallel do
            case DEFAULT
                THROW_HARD('Unsupported ref_type mode. It should be cavg, clin, or vol')
        end select
        ! Directional FRCs
        if( trim(params_glob%polar_frcs) .eq. 'yes' )then
            pfts_merg = pfts_merg * cavg_clin_frcs
            ! res = get_resarr(params_glob%box_crop, params_glob%smpd_crop)
            ! ! min/max frc between cavg and clin
            ! min_res_fsc0143 = HUGE(min_res_fsc0143)
            ! max_res_fsc0143 = 0.
            ! avg_res_fsc0143 = 0.
            ! avg_res_fsc05   = 0.
            ! npops           = 0
            ! do icls = 1, ncls
            !     if( pops(icls) < 2 )cycle
            !     npops = npops + 1
            !     call get_resolution_at_fsc(cavg_clin_frcs(:,icls), res, 0.5,   res_fsc05)
            !     call get_resolution_at_fsc(cavg_clin_frcs(:,icls), res, 0.143, res_fsc0143)
            !     avg_res_fsc0143 = avg_res_fsc0143 + res_fsc0143
            !     avg_res_fsc05   = avg_res_fsc05   + res_fsc05
            !     if( res_fsc0143 < min_res_fsc0143 ) min_res_fsc0143 = res_fsc0143
            !     if( res_fsc0143 > max_res_fsc0143 ) max_res_fsc0143 = res_fsc0143
            ! enddo
            ! avg_res_fsc05   = avg_res_fsc05   / real(npops)
            ! avg_res_fsc0143 = avg_res_fsc0143 / real(npops)
            ! write(logfhandle,'(A,2F8.2)')'>>> AVG CAVG/CLIN DIRECTIONAL RESOLUTION @ FSC=0.5/0.143: ',avg_res_fsc05, avg_res_fsc0143
            ! write(logfhandle,'(A,2F8.2)')'>>> MIN/MAX CAVG/CLIN DIRECTIONAL RESOLUTION @ FSC=0.143: ',min_res_fsc0143, max_res_fsc0143
            ! ! 3D FSC = AVERAGE RESOLUTION
            ! write(logfhandle,'(A,2F8.2)')'>>> 3D CAVG/CLIN RESOLUTION @ FSC=0.5/0.143             : ',avg_res_fsc05,avg_res_fsc0143
            ! ! computing directional frcs
            ! dfrcs = 0.
            ! !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            ! !$omp private(icls,k,numer,denom1,denom2)
            ! do icls = 1, ncls
            !     do k = kfromto(1), kfromto(2)
            !         numer  = real(sum(    pfts_even(:,k,icls) * conjg(pfts_odd(:,k,icls))), dp)
            !         denom1 =      sum(csq_fast(pfts_even(:,k,icls)))
            !         denom2 =      sum(csq_fast(pfts_odd( :,k,icls)))
            !         if( denom1*denom2 > DTINY ) dfrcs(k,icls) = real(numer / dsqrt(denom1*denom2))
            !     enddo
            ! enddo
            ! !$omp end parallel do
            ! ! min/max directional frcs
            ! min_res_fsc0143 = HUGE(min_res_fsc0143)
            ! max_res_fsc0143 = 0.
            ! avg_res_fsc0143 = 0.
            ! avg_res_fsc05   = 0.
            ! npops           = 0
            ! do icls = 1, ncls
            !     if( pops(icls) < 2 )cycle
            !     npops = npops + 1
            !     call get_resolution_at_fsc(dfrcs(:,icls), res, 0.5,   res_fsc05)
            !     call get_resolution_at_fsc(dfrcs(:,icls), res, 0.143, res_fsc0143)
            !     avg_res_fsc0143 = avg_res_fsc0143 + res_fsc0143
            !     avg_res_fsc05   = avg_res_fsc05   + res_fsc05
            !     if( res_fsc0143 < min_res_fsc0143 ) min_res_fsc0143 = res_fsc0143
            !     if( res_fsc0143 > max_res_fsc0143 ) max_res_fsc0143 = res_fsc0143
            ! enddo
            ! avg_res_fsc05   = avg_res_fsc05   / real(npops)
            ! avg_res_fsc0143 = avg_res_fsc0143 / real(npops)
            ! write(logfhandle,'(A,2F8.2)')'>>> AVG EVEN/ODD DIRECTIONAL RESOLUTION @ FSC=0.5/0.143 : ',avg_res_fsc05,avg_res_fsc0143
            ! write(logfhandle,'(A,2F8.2)')'>>> MIN/MAX EVEN/ODD DIRECTIONAL RESOLUTION @ FSC=0.143 : ', min_res_fsc0143,max_res_fsc0143
            ! ! 3D FSC = AVERAGE RESOLUTION
            ! write(logfhandle,'(A,2F8.2)')'>>> 3D EVEN/ODD RESOLUTION @ FSC=0.5/0.143              : ',avg_res_fsc05,avg_res_fsc0143
        endif

      contains

        ! nonuniform frc filter
        subroutine nonuni_frc(pft1, pft2, pfrc)
            complex(dp), intent(in)    :: pft1(pftsz)
            complex(dp), intent(in)    :: pft2(pftsz)
            real,        intent(inout) :: pfrc(pftsz)
            integer,     parameter     :: NKER = 3, HALF = NKER/2, MID = HALF + 1
            integer  :: i
            real(dp) :: numer, denom1, denom2
            do i = MID, pftsz - HALF
                pfrc(i) = 1._dp
                numer   = real(sum(         pft1(i-HALF:i+HALF) * conjg(pft2(i-HALF:i+HALF))),dp)
                denom1  =      sum(csq_fast(pft1(i-HALF:i+HALF)))
                denom2  =      sum(csq_fast(pft2(i-HALF:i+HALF)))
                if( denom1*denom2 > DTINY ) pfrc(i) = real(numer / dsqrt(denom1*denom2))
            enddo
            do i = 1, MID-1
                pfrc(i) = 1._dp
                numer   = real(sum(         pft1(           1:i+HALF) * conjg(pft2(           1:i+HALF))),dp)
                numer   = real(sum(         pft1(pftsz-HALF+i:pftsz ) * conjg(pft2(pftsz-HALF+i:pftsz ))),dp) + numer
                denom1  =      sum(csq_fast(pft1(1:i+HALF))) + sum(csq_fast(pft1(pftsz-HALF+i:pftsz)))
                denom2  =      sum(csq_fast(pft2(1:i+HALF))) + sum(csq_fast(pft2(pftsz-HALF+i:pftsz)))
                if( denom1*denom2 > DTINY ) pfrc(i) = real(numer / dsqrt(denom1*denom2))
            enddo
            do i = pftsz-HALF+1, pftsz
                pfrc(i) = 1._dp
                numer   = real(sum(         pft1(     1:HALF-(pftsz-i)) * conjg(pft2(     1:HALF-(pftsz-i)))),dp)
                numer   = real(sum(         pft1(i-HALF:pftsz         ) * conjg(pft2(i-HALF:pftsz         ))),dp) + numer
                denom1  =      sum(csq_fast(pft1(1:HALF-(pftsz-i)))) + sum(csq_fast(pft1(i-HALF:pftsz)))
                denom2  =      sum(csq_fast(pft2(1:HALF-(pftsz-i)))) + sum(csq_fast(pft2(i-HALF:pftsz)))
                if( denom1*denom2 > DTINY ) pfrc(i) = real(numer / dsqrt(denom1*denom2))
            enddo
        end subroutine nonuni_frc

        ! Deals with summing slices and their mirror
        subroutine mirror_slices( ref_space, symop )
            type(oris), intent(in) :: ref_space
            type(sym),  intent(in) :: symop
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: iref, m
            logical     :: l_rotm
            if( .not.ref_space%isthere('mirr') )then
                THROW_HARD('Mirror index missing in reference search space')
            endif
            !$omp parallel do default(shared) proc_bind(close) private(iref,m,pft,ctf2,psi,l_rotm)
            do iref = 1,ncls/2
                m      = ref_space%get_int(iref,'mirr')
                psi    = abs(ref_space%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                ! Fourier components
                if( l_rotm )then
                    call mirror_pft(pfts_even(:,:,m), pft)
                    pfts_even(:,:,iref) = pfts_even(:,:,iref) + conjg(pft)
                    call mirror_pft(pfts_odd(:,:,m), pft)
                    pfts_odd(:,:,iref)  = pfts_odd(:,:,iref)  + conjg(pft)
                    call mirror_pft(conjg(pfts_even(:,:,iref)), pfts_even(:,:,m))
                    call mirror_pft(conjg(pfts_odd(:,:,iref)),  pfts_odd(:,:,m))
                else
                    call mirror_pft(pfts_even(:,:,m), pft)
                    pfts_even(:,:,iref) = pfts_even(:,:,iref) + pft
                    call mirror_pft(pfts_odd(:,:,m), pft)
                    pfts_odd(:,:,iref)  = pfts_odd(:,:,iref)  + pft
                    call mirror_pft(pfts_even(:,:,iref), pfts_even(:,:,m))
                    call mirror_pft(pfts_odd(:,:,iref),  pfts_odd(:,:,m))
                endif
                ! CTF
                call mirror_ctf2(ctf2_even(:,:,m), ctf2)
                ctf2_even(:,:,iref) = ctf2_even(:,:,iref) + ctf2
                call mirror_ctf2(ctf2_odd(:,:,m),  ctf2)
                ctf2_odd(:,:,iref)  = ctf2_odd(:,:,iref)  + ctf2
                call mirror_ctf2(ctf2_even(:,:,iref), ctf2_even(:,:,m))
                call mirror_ctf2(ctf2_odd(:,:,iref),  ctf2_odd(:,:,m))
            enddo
            !$omp end parallel do
        end subroutine mirror_slices

        subroutine get_cavg_clin
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,eo_pop,pop,pft,ctf2)
            do icls = 1,ncls
                eo_pop     = prev_eo_pops(:,icls) + eo_pops(:,icls) ! eo_pops has to be calculated differently
                pop        = sum(eo_pop)
                pops(icls) = pop
                if(pop == 0)then
                    pfts_cavg(:,:,icls) = DCMPLX_ZERO
                else
                    if(pop > 1)then
                        pft  = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                        ctf2 = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                        if( pop <= 5 ) ctf2 = ctf2 + real(EPSILON/real(pop),dp)
                        call safe_norm(pft, ctf2, pfts_cavg(:,:,icls))
                    endif
                endif
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,pft,ctf2)
            do icls = 1,ncls
                call safe_norm(pfts_clin_even(:,:,icls), ctf2_clin_even(:,:,icls), clin_even(:,:,icls))
                call safe_norm(pfts_clin_odd( :,:,icls), ctf2_clin_odd( :,:,icls), clin_odd( :,:,icls))
                pft  = pfts_clin_even(:,:,icls) + pfts_clin_odd(:,:,icls)
                ctf2 = ctf2_clin_even(:,:,icls) + ctf2_clin_odd(:,:,icls)
                call safe_norm(pft, ctf2, pfts_clin(:,:,icls))
            enddo
            !$omp end parallel do
        end subroutine get_cavg_clin

        ! Extrapolate cline, rline to pfts_clin and ctf2_clin
        subroutine extrapolate_line(ref, rot, w, cline_e, cline_o, rline_e, rline_o)
            integer,     intent(in) :: ref, rot
            real(dp),    intent(in) :: w
            complex(dp), intent(in) :: cline_e(kfromto(1):kfromto(2))
            complex(dp), intent(in) :: cline_o(kfromto(1):kfromto(2))
            real(dp),    intent(in) :: rline_e(kfromto(1):kfromto(2))
            real(dp),    intent(in) :: rline_o(kfromto(1):kfromto(2))
            integer :: irot
            irot = rot
            if( irot < 1 )then
                irot = irot + pftsz
                pfts_clin_even(irot,:,ref) = pfts_clin_even(irot,:,ref) + w * conjg(cline_e)
                pfts_clin_odd( irot,:,ref) = pfts_clin_odd( irot,:,ref) + w * conjg(cline_o)
            elseif( irot > pftsz )then
                irot = irot - pftsz
                pfts_clin_even(irot,:,ref) = pfts_clin_even(irot,:,ref) + w * conjg(cline_e)
                pfts_clin_odd( irot,:,ref) = pfts_clin_odd( irot,:,ref) + w * conjg(cline_o)
            else
                pfts_clin_even(irot,:,ref) = pfts_clin_even(irot,:,ref) + w * cline_e
                pfts_clin_odd( irot,:,ref) = pfts_clin_odd( irot,:,ref) + w * cline_o
            endif
            ctf2_clin_even(irot,:,ref) = ctf2_clin_even(irot,:,ref) + w * rline_e
            ctf2_clin_odd( irot,:,ref) = ctf2_clin_odd( irot,:,ref) + w * rline_o
        end subroutine extrapolate_line

        subroutine calc_comlin_contrib( ref_space, symop )
            type(oris), intent(in) :: ref_space
            type(sym),  intent(in) :: symop
            type(kbinterpol)  :: kbwin
            real, allocatable :: Rsym(:,:,:)
            complex(dp) :: cline_l(kfromto(1):kfromto(2)), cline_r(kfromto(1):kfromto(2))
            complex(dp) :: cline_e(kfromto(1):kfromto(2)), cline_o(kfromto(1):kfromto(2))
            real(dp)    :: rline_l(kfromto(1):kfromto(2)), rline_r(kfromto(1):kfromto(2))
            real(dp)    :: rline_e(kfromto(1):kfromto(2)), rline_o(kfromto(1):kfromto(2))
            real(dp)    :: dd, w, wl, wr, sumw
            real        :: eulers(3),R(3,3,ncls),Rj(3,3),tRi(3,3),psi,drot,d,targ_w,self_w
            integer     :: rotl,rotr, iref, jref, m, self_irot, targ_irot, isym, nsym
            logical     :: l_rotm
            if( .not.ref_space%isthere('mirr') )then
                THROW_HARD('Mirror index missing in reference search space')
            endif
            drot = pftcc_glob%get_dang()
            if( l_kb ) kbwin = kbinterpol(1.5, KBALPHA)
            pfts_clin_even = DCMPLX_ZERO; pfts_clin_odd = DCMPLX_ZERO
            ctf2_clin_even = 0.d0; ctf2_clin_odd  = 0.d0
            ! Symmetry rotation matrices
            nsym = symop%get_nsym()
            allocate(Rsym(3,3,nsym),source=0.)
            !$omp parallel default(shared) proc_bind(close)&
            !$omp private(iref,jref,m,tRi,Rj,eulers,targ_irot,self_irot,targ_w,self_w)&
            !$omp& private(d,dd,w,wl,wr,sumw,psi,l_rotm,cline_l,cline_r,cline_e,cline_o)&
            !$omp& private(rline_l,rline_r,rline_e,rline_o,pft,ctf2,rotl,rotr,isym)
            ! Caching rotation matrices
            !$omp do schedule(static)
            do iref = 1,ncls
                R(:,:,iref) = ref_space%get_mat(iref)
            enddo
            !$omp end do nowait
            !$omp do schedule(static)
            do isym = 1,nsym
                call symop%get_sym_rmat(isym, Rsym(:,:,isym))
            end do
            !$omp end do
            ! Common lines contribution
            !$omp do schedule(static)
            do iref = 1,ncls/2
                tRi = transpose(R(:,:,iref))
                m   = ref_space%get_int(iref,'mirr')
                do jref = 1,ncls/2
                    do isym = 1,nsym
                        if( isym == 1 )then
                            if( jref == iref ) cycle    ! self   exclusion
                            if( jref == m )    cycle    ! mirror exclusion
                            ! Rotation of both planes by transpose of Ri (tRixRi=I & Rsym=I)
                            Rj = matmul(R(:,:,jref), tRi)
                        else
                            ! Symmetry operator
                            Rj = matmul(R(:,:,jref), Rsym(:,:,isym))
                            ! Rotation of both planes by transpose of Ri (tRixRi -> I)
                            Rj = matmul(Rj, tRi)
                        endif
                        ! Euler angles identification
                        eulers = m2euler_fast(Rj)
                        ! Interpolation
                        if( l_kb )then
                            ! KB
                            ! in plane rotation index of jref slice intersecting iref
                            psi       = 360.0 - eulers(3)
                            targ_irot = pftcc_glob%get_roind_fast(psi)
                            d         = psi - pftcc_glob%get_rot(targ_irot)
                            if( d > drot ) d = d - 360.0
                            targ_w    = d / drot
                            ! in plane rotation index of iref slice
                            psi       = eulers(1)
                            self_irot = pftcc_glob%get_roind_fast(psi)
                            d         = psi - pftcc_glob%get_rot(self_irot)
                            if( d > drot ) d = d - 360.0
                            self_w    = d / drot
                            ! intepolate common line in jref-th slice
                            rotl = targ_irot - 1; rotr = targ_irot + 1
                            dd   = real(targ_w,dp)
                            w    = kbwin%apod_dp(dd); wl = kbwin%apod_dp(dd-1.d0); wr = kbwin%apod_dp(dd+1.d0)
                            sumw = wl + w + wr
                            w    = w / sumw; wl = wl / sumw; wr = wr / sumw
                            call get_line(jref, targ_irot, .true., cline_e, rline_e)
                            call get_line(jref, rotl,      .true., cline_l, rline_l)
                            call get_line(jref, rotr,      .true., cline_r, rline_r)
                            cline_e = wl*cline_l + w*cline_e + wr*cline_r
                            rline_e = wl*rline_l + w*rline_e + wr*rline_r
                            call get_line(jref, targ_irot, .false., cline_o, rline_o)
                            call get_line(jref, rotl,      .false., cline_l, rline_l)
                            call get_line(jref, rotr,      .false., cline_r, rline_r)
                            cline_o = wl*cline_l + w*cline_o + wr*cline_r
                            rline_o = wl*rline_l + w*rline_o + wr*rline_r
                            ! extrapolate the common line to iref-th slice
                            rotl = self_irot - 1; rotr = self_irot + 1
                            if( rotr > nrots ) rotr = rotr - nrots
                            dd   = real(self_w,dp)
                            w    = kbwin%apod_dp(dd); wl = kbwin%apod_dp(dd-1.d0); wr = kbwin%apod_dp(dd+1.d0)
                            sumw = wl + w + wr
                            w    = w / sumw; wl = wl / sumw; wr = wr / sumw
                            ! leftmost line
                            call extrapolate_line(iref, rotl,      wl, cline_e, cline_o, rline_e, rline_o)
                            ! nearest line
                            call extrapolate_line(iref, self_irot, w,  cline_e, cline_o, rline_e, rline_o)
                            ! rightmost line
                            call extrapolate_line(iref, rotr,      wr, cline_e, cline_o, rline_e, rline_o)
                        else
                            ! Linear interpolation
                            ! in plane rotation index of jref slice intersecting iref
                            psi       = 360.0 - eulers(3)
                            targ_irot = pftcc_glob%get_roind_fast(psi)
                            d         = psi - pftcc_glob%get_rot(targ_irot)
                            if( d > drot ) d = d - 360.0
                            if( d < 0. )then
                                targ_irot = targ_irot - 1
                                if( targ_irot < 1 ) targ_irot = targ_irot + nrots
                                d = d + drot
                            endif
                            targ_w = d / drot
                            ! in plane rotation index of iref slice
                            psi       = eulers(1)
                            self_irot = pftcc_glob%get_roind_fast(psi)
                            d         = psi - pftcc_glob%get_rot(self_irot)
                            if( d > drot ) d = d - 360.0
                            if( d < 0. )then
                                self_irot = self_irot - 1
                                if( self_irot < 1 ) self_irot = self_irot + nrots
                                d = d + drot
                            endif
                            self_w = d / drot
                            ! intepolate common line in jref-th slice
                            rotl = targ_irot; rotr = rotl+1
                            wl   = real(targ_w,dp); wr = 1.d0 - wl
                            call get_line(jref, rotl, .true., cline_e, rline_e)
                            call get_line(jref, rotr, .true., cline_r, rline_r)
                            cline_e = wr*cline_e + wl*cline_r
                            rline_e = wr*rline_e + wl*rline_r
                            call get_line(jref, rotl, .false., cline_o, rline_o)
                            call get_line(jref, rotr, .false., cline_r, rline_r)
                            cline_o = wr*cline_o + wl*cline_r
                            rline_o = wr*rline_o + wl*rline_r
                            ! extrapolate the common line to iref-th slice
                            rotl = self_irot; rotr = rotl + 1
                            if( rotr > nrots ) rotr = rotr - nrots
                            w  = real(self_w,dp)
                            call extrapolate_line(iref, rotl, 1.d0-w, cline_e, cline_o, rline_e, rline_o)
                            call extrapolate_line(iref, rotr,      w, cline_e, cline_o, rline_e, rline_o)
                        endif
                    enddo
                enddo
            enddo
            !$omp end do
            ! Mirroring contributions
            !$omp do schedule(static)
            do iref = 1,ncls/2
                m      = ref_space%get_int(iref,'mirr')
                psi    = abs(ref_space%get(m, 'psi'))
                l_rotm = (psi>0.1) .and. (psi<359.9)
                if( l_rotm )then
                    call mirror_pft(conjg(pfts_clin_even(:,:,iref)), pfts_clin_even(:,:,m))
                    call mirror_pft(conjg(pfts_clin_odd(:,:,iref)),  pfts_clin_odd(:,:,m))
                else
                    call mirror_pft(pfts_clin_even(:,:,iref), pfts_clin_even(:,:,m))
                    call mirror_pft(pfts_clin_odd(:,:,iref),  pfts_clin_odd(:,:,m))
                endif
                call mirror_ctf2(ctf2_clin_even(:,:,iref), ctf2_clin_even(:,:,m))
                call mirror_ctf2(ctf2_clin_odd(:,:,iref),  ctf2_clin_odd(:,:,m))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine calc_comlin_contrib

    end subroutine polar_cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    subroutine polar_cavger_calc_and_write_frcs_and_eoavg( fname, cline )
        use simple_cmdline, only: cmdline
        character(len=*), intent(in) :: fname
        type(cmdline),    intent(in) :: cline
        complex(dp), allocatable :: prev_prefs(:,:,:)
        real,        allocatable :: frc(:)
        real(dp) :: ufrac_trec
        integer  :: icls, find, pop, filtsz
        filtsz = fdim(params_glob%box_crop) - 1
        allocate(frc(filtsz),source=0.)
        !$omp parallel do default(shared) private(icls,frc,find,pop) schedule(static) proc_bind(close)
        do icls = 1,ncls
            if( l_comlin )then
                ! calculate FRC (pseudo-cavgs are never empty)
                call calc_frc(pfts_even(:,:,icls), pfts_odd(:,:,icls), filtsz, frc)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = min(kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                if( find >= kfromto(1) )then
                    pfts_even(:,kfromto(1):find,icls) = pfts_merg(:,kfromto(1):find,icls)
                    pfts_odd(:,kfromto(1):find,icls)  = pfts_merg(:,kfromto(1):find,icls)
                endif
            else
                pop = sum(prev_eo_pops(:,icls) + eo_pops(:,icls))
                if( pop == 0 )then
                    frc = 0.
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                else
                    ! calculate FRC
                    call calc_frc(pfts_even(:,:,icls), pfts_odd(:,:,icls), filtsz, frc)
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                    ! average low-resolution info between eo pairs to keep things in register
                    find = min(kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                    if( find >= kfromto(1) )then
                        pfts_even(:,kfromto(1):find,icls) = pfts_merg(:,kfromto(1):find,icls)
                        pfts_odd(:,kfromto(1):find,icls)  = pfts_merg(:,kfromto(1):find,icls)
                    endif
                endif
            endif
        end do
        !$omp end parallel do
        ! write FRCs
        call build_glob%clsfrcs%write(fname)
        ! e/o Trailing reconstruction
        if( params_glob%l_trail_rec )then
            if( cline%defined('ufrac_trec') )then
                ufrac_trec = real(params_glob%ufrac_trec,dp)
            else
                ufrac_trec = real(build_glob%spproj_field%get_update_frac(),dp)
            endif
            call read_pft_array(POLAR_REFS_FBODY//'_even'//trim(BIN_EXT), prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_even = ufrac_trec * pfts_even + (1.d0-ufrac_trec) * prev_prefs
            !$omp end parallel workshare
            call read_pft_array(POLAR_REFS_FBODY//'_odd'//trim(BIN_EXT), prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_odd  = ufrac_trec * pfts_odd   + (1.d0-ufrac_trec) * prev_prefs
            pfts_merg = 0.5d0 * (pfts_even + pfts_odd)
            !$omp end parallel workshare
            deallocate(prev_prefs)
        endif
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    !>  \brief  Converts the polar references to a cartesian grid
    subroutine polar_cavger_refs2cartesian( pftcc, cavgs, which, pfts_in )
        use simple_image
        class(polarft_corrcalc), intent(in)    :: pftcc
        type(image),             intent(inout) :: cavgs(ncls)
        character(len=*),        intent(in)    :: which
        complex(dp),   optional, intent(in)    :: pfts_in(1:pftsz,kfromto(1):kfromto(2),1:ncls)
        complex, allocatable :: cmat(:,:)
        real,    allocatable :: norm(:,:)
        complex :: pft(1:pftsz,kfromto(1):kfromto(2)), fc
        real    :: phys(2), dh,dk,mdk,mdh
        integer :: k,c,irot,physh,physk,box,icls
        box = params_glob%box_crop
        c   = box/2+1
        allocate(cmat(c,box),norm(c,box))
        !$omp parallel do schedule(guided) proc_bind(close) default(shared)&
        !$omp private(icls,pft,cmat,norm,irot,k,phys,fc,physh,physk,dh,dk,mdh,mdk)
        do icls = 1, ncls
            if( present(pfts_in) )then
                pft = cmplx(pfts_in(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
            else
                select case(trim(which))
                case('even')
                    pft = cmplx(pfts_even(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                case('odd')
                    pft = cmplx(pfts_odd(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                case('merged')
                    pft = cmplx(pfts_merg(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                end select
            endif
            ! Bi-linear interpolation
            cmat = CMPLX_ZERO
            norm = 0.0
            do irot = 1,pftsz
                do k = kfromto(1),kfromto(2)
                    phys  = pftcc%get_coord(irot,k) + [1.,real(c)]
                    fc    = pft(irot,k)
                    physh = floor(phys(1))
                    physk = floor(phys(2))
                    dh    = phys(1) - real(physh)
                    dk    = phys(2) - real(physk)
                    mdh   = 1.0 - dh
                    mdk   = 1.0 - dk
                    if( physh > 0 .and. physh <= c )then
                        if( physk <= box .and. physk >= 1 )then
                            cmat(physh,physk) = cmat(physh,physk) + mdh*mdk*fc
                            norm(physh,physk) = norm(physh,physk) + mdh*mdk
                            if( physk+1 <= box .and. physk+1 >= 1 )then
                                cmat(physh,physk+1) = cmat(physh,physk+1) + mdh*dk*fc
                                norm(physh,physk+1) = norm(physh,physk+1) + mdh*dk
                            endif
                        endif
                    endif
                    physh = physh + 1
                    if( physh > 0 .and. physh <= c )then
                        if( physk <= box .and. physk >= 1 )then
                            cmat(physh,physk) = cmat(physh,physk) + dh*mdk*fc
                            norm(physh,physk) = norm(physh,physk) + dh*mdk
                            if( physk+1 <= box .and. physk+1 >= 1 )then
                                cmat(physh,physk+1) = cmat(physh,physk+1) + dh*dk*fc
                                norm(physh,physk+1) = norm(physh,physk+1) + dh*dk
                            endif
                        endif
                    endif
                end do
            end do
            where( norm > TINY )
                cmat = cmat / norm
            elsewhere
                cmat = 0.0
            end where
            ! irot = self%pftsz+1, eg. angle=180.
            do k = 1,box/2-1
                cmat(1,k+c) = conjg(cmat(1,c-k))
            enddo
            ! arbitrary magnitude
            cmat(1,c) = CMPLX_ZERO
            ! set image
            call cavgs(icls)%new([box,box,1], smpd, wthreads=.false.)
            call cavgs(icls)%set_cmat(cmat)
            call cavgs(icls)%shift_phorig()
            call cavgs(icls)%ifft
        enddo
        !$omp end parallel do
    end subroutine polar_cavger_refs2cartesian

    !>  \brief  Reads in and reduces partial matrices prior to restoration
    subroutine polar_cavger_assemble_sums_from_parts( reforis )
        type(oris), optional, intent(in) :: reforis
        character(len=:), allocatable :: cae, cao, cte, cto
        complex(dp),      allocatable :: pfte(:,:,:), pfto(:,:,:)
        real(dp),         allocatable :: ctf2e(:,:,:), ctf2o(:,:,:)
        integer :: ipart
        allocate(pfte(pftsz,kfromto(1):kfromto(2),ncls),  pfto(pftsz,kfromto(1):kfromto(2),ncls),&
               &ctf2e(pftsz,kfromto(1):kfromto(2),ncls), ctf2o(pftsz,kfromto(1):kfromto(2),ncls))
        call polar_cavger_zero_pft_refs
        do ipart = 1,params_glob%nparts
            cae = 'cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            call read_pft_array(cae, pfte)
            call read_pft_array(cao, pfto)
            call read_ctf2_array(cte, ctf2e)
            call read_ctf2_array(cto, ctf2o)
            !$omp parallel workshare proc_bind(close)
            pfts_even = pfts_even + pfte
            pfts_odd  = pfts_odd  + pfto
            ctf2_even = ctf2_even + ctf2e
            ctf2_odd  = ctf2_odd  + ctf2o
            !$omp end parallel workshare
        enddo
        ! merge eo-pairs and normalize
        call polar_cavger_merge_eos_and_norm(reforis=reforis)
    end subroutine polar_cavger_assemble_sums_from_parts

    ! I/O

    ! Writes cavgs PFT array
    subroutine polar_cavger_write( fname, which )
        character(len=*),  intent(in) :: fname, which
        character(len=:), allocatable :: fname_here
        fname_here  = trim(fname)
        select case(which)
            case('even')
                call write_pft_array(pfts_even, fname_here)
            case('odd')
                call write_pft_array(pfts_odd,  fname_here)
            case('merged')
                call write_pft_array(pfts_merg, fname_here)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_write

    ! Writes all cavgs PFT arrays
    subroutine polar_cavger_writeall( tmpl_fname )
        character(len=*),  intent(in) :: tmpl_fname
        call polar_cavger_write(trim(tmpl_fname)//'_even'//BIN_EXT,'even')
        call polar_cavger_write(trim(tmpl_fname)//'_odd'//BIN_EXT, 'odd')
        call polar_cavger_write(trim(tmpl_fname)//BIN_EXT,         'merged')
    end subroutine polar_cavger_writeall

    ! Write references contained in PFTCC
    subroutine polar_cavger_writeall_pftccrefs( tmpl_fname )
        character(len=*),  intent(in) :: tmpl_fname
        complex(sp),  pointer :: ptre(:,:,:), ptro(:,:,:)
        call pftcc_glob%get_refs_ptr( ptre, ptro )
        pfts_even = cmplx(ptre,kind=dp)
        pfts_odd  = cmplx(ptro,kind=dp)
        call polar_cavger_write(trim(tmpl_fname)//'_even'//BIN_EXT,'even')
        call polar_cavger_write(trim(tmpl_fname)//'_odd'//BIN_EXT, 'odd')
        call polar_cavger_zero_pft_refs !! removed after writing
        nullify(ptre, ptro)
    end subroutine polar_cavger_writeall_pftccrefs

    ! Converts cavgs PFTS to cartesian grids and writes them
    subroutine polar_cavger_write_cartrefs( pftcc, tmpl_fname, which )
        class(polarft_corrcalc), intent(in) :: pftcc
        character(len=*),  intent(in)       :: tmpl_fname, which
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(ncls, [params_glob%box_crop, params_glob%box_crop,1], smpd, imgs)
        select case(trim(which))
        case('even','odd')
            call polar_cavger_refs2cartesian( pftcc, imgs, trim(which) )
            call write_cavgs(imgs, trim(tmpl_fname)//'_'//trim(which)//params_glob%ext)
        case('merged')
            call polar_cavger_refs2cartesian( pftcc, imgs, 'merged' )
            call write_cavgs(imgs, trim(tmpl_fname)//params_glob%ext)
        end select
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_write_cartrefs

    ! Converts all cavgs PFTS to cartesian grids and writes them
    subroutine polar_cavger_writeall_cartrefs( pftcc, tmpl_fname )
        class(polarft_corrcalc), intent(in) :: pftcc
        character(len=*),  intent(in)       :: tmpl_fname
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(ncls, [params_glob%box_crop, params_glob%box_crop,1], smpd, imgs)
        call polar_cavger_refs2cartesian( pftcc, imgs, 'even' )
        call write_cavgs(imgs, trim(tmpl_fname)//'_even'//params_glob%ext)
        call polar_cavger_refs2cartesian( pftcc, imgs, 'odd' )
        call write_cavgs(imgs, trim(tmpl_fname)//'_odd'//params_glob%ext)
        call polar_cavger_refs2cartesian( pftcc, imgs, 'merged' )
        call write_cavgs(imgs, trim(tmpl_fname)//params_glob%ext)
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_writeall_cartrefs

    ! Read cavgs PFT array
    subroutine polar_cavger_read( fname, which )
        character(len=*),  intent(in) :: fname, which
        character(len=:), allocatable :: fname_here
        fname_here  = trim(fname)
        select case(which)
            case('even')
                call read_pft_array(fname_here, pfts_even)
            case('odd')
                call read_pft_array(fname_here, pfts_odd)
            case('merged')
                call read_pft_array(fname_here, pfts_merg)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_read

    ! Reads all cavgs PFT arrays
    subroutine polar_cavger_read_all( fname )
        character(len=*),  intent(in) :: fname
        character(len=:), allocatable :: refs, refs_even, refs_odd, ext
        ext = '.'//trim(fname2ext(fname))
        if( ext == trim(params_glob%ext) )then
            refs = trim(get_fbody(trim(fname), params_glob%ext, separator=.false.))//BIN_EXT
        elseif( ext == BIN_EXT )then
            refs = trim(fname)
        else
            THROW_HARD('Unsupported file format: '//ext)
        endif
        refs_even = trim(get_fbody(refs,BIN_EXT,separator=.false.))//'_even'//BIN_EXT
        refs_odd  = trim(get_fbody(refs,BIN_EXT,separator=.false.))//'_odd'//BIN_EXT
        if( .not. file_exists(refs) )then
            THROW_HARD('Polar references do not exist in cwd: '//trim(refs))
        endif
        call polar_cavger_read(refs, 'merged')
        if( file_exists(refs_even) )then
            call polar_cavger_read(refs_even, 'even')
        else
            call polar_cavger_read(refs, 'even')
        endif
        if( file_exists(refs_odd) )then
            call polar_cavger_read(refs_odd, 'odd')
        else
            call polar_cavger_read(refs, 'odd')
        endif
    end subroutine polar_cavger_read_all

    !>  \brief  writes partial class averages (PFTS + CTF2) to disk (distributed execution)
    subroutine polar_cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        character(len=:), allocatable :: cae, cao, cte, cto
        cae = 'cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cao = 'cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cte = 'ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cto = 'ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        select case(trim(which))
            case('read')
                call read_pft_array(cae, pfts_even)
                call read_pft_array(cao, pfts_odd)
                call read_ctf2_array(cte, ctf2_even)
                call read_ctf2_array(cto, ctf2_odd)
            case('write')
                call write_pft_array(pfts_even, cae)
                call write_pft_array(pfts_odd,  cao)
                call write_ctf2_array(ctf2_even, cte)
                call write_ctf2_array(ctf2_odd,  cto)
            case DEFAULT
                THROW_HARD('unknown which flag; only read & write supported; cavger_readwrite_partial_sums')
        end select
        deallocate(cae, cao, cte, cto)
    end subroutine polar_cavger_readwrite_partial_sums

    ! BOOK-KEEPING

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    subroutine polar_cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), target, intent(inout) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(ncls)
        real(dp) :: corrs(ncls), ws(ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: iptcl, icls, pop, nptcls
        logical  :: l_3D
        l_3D = .false.
        select case(trim(params_glob%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(params_glob%oritype))
        end select
        nptcls = ptcl_field%get_noris()
        pops   = 0
        corrs  = 0.d0
        ws     = 0.d0
        !$omp parallel do default(shared) private(iptcl,rstate,icls,w) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs,ws)
        do iptcl=1,nptcls
            rstate = ptcl_field%get(iptcl,'state')
            if( rstate < 0.5 ) cycle
            w = ptcl_field%get(iptcl,'w')
            if( w < SMALL ) cycle
            if( l_3D )then
                icls = ptcl_field%get_proj(iptcl)
            else
                icls = ptcl_field%get_class(iptcl)
            endif
            if( icls<1 .or. icls>params_glob%ncls )cycle
            pops(icls)  = pops(icls)  + 1
            corrs(icls) = corrs(icls) + real(ptcl_field%get(iptcl,'corr'),dp)
            ws(icls)    = ws(icls)    + real(w,dp)
        enddo
        !$omp end parallel do
        where(pops>1)
            corrs = corrs / real(pops)
            ws    = ws / real(pops)
        elsewhere
            corrs = -1.
            ws    = 0.
        end where
        call cls_field%new(ncls, is_ptcl=.false.)
        do icls=1,ncls
            pop = pops(icls)
            call build_glob%clsfrcs%estimate_res(icls, frc05, frc0143)
            if( l_3D )then
                call cls_field%set(icls, 'proj',  icls)
            else
                call cls_field%set(icls, 'class', icls)
            endif
            call cls_field%set(icls, 'pop',       pop)
            call cls_field%set(icls, 'res',       frc0143)
            call cls_field%set(icls, 'corr',      corrs(icls))
            call cls_field%set(icls, 'w',         ws(icls))
            call cls_field%set_state(icls, 1) ! needs to be default val if no selection has been done
            if( pop == 0 )call cls_field%set_state(icls, 0)
        end do
    end subroutine polar_cavger_gen2Dclassdoc

    ! ALIGNEMENT

    !>  \brief  prepares one polar cluster centre image for alignment
    subroutine polar_prep2Dref( icls, cavg, center, xyz )
        integer,                intent(in)    :: icls
        class(image), optional, intent(inout) :: cavg
        logical,      optional, intent(in)    :: center
        real,         optional, intent(out)   :: xyz(3)
        real, allocatable :: frc(:), filter(:), gaufilter(:)
        real    :: xy_cavg(2), crop_factor
        integer :: filtsz, k
        logical :: do_center
        ! Centering
        if( present(cavg).and.present(center).and.present(xyz) )then
            ! The three optional arguments are required
            do_center = .false.
            if( present(center) )then
                xyz       = 0.
                do_center = center.and.(params_glob%center .eq. 'yes')
            endif
            if( do_center )then
                crop_factor = real(params_glob%box_crop) / real(params_glob%box)
                select case(trim(params_glob%center_type))
                case('params')
                    ! offset from document in original pixel unit
                    call build_glob%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                    if( arg(xy_cavg) < CENTHRESH )then
                        xyz = 0.
                    else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                        xyz = [xy_cavg(1), xy_cavg(2), 0.]
                    else
                        xyz = cavg%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
                        xyz = xyz / crop_factor         ! scaled pixel unit
                        if( arg(xyz(1:2) - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                    endif
                case('seg')
                    call calc_cavg_offset(cavg, params_glob%cenlp, params_glob%msk_crop, xy_cavg)
                    xyz(1:2) = xy_cavg / crop_factor    ! scaled pixel unit
                case('mass')
                    xyz = cavg%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
                    xyz = xyz / crop_factor             ! scaled pixel unit
                end select
            endif
        endif
        ! Filtering
        if( params_glob%l_ml_reg )then
            ! no filtering, not supported yet
        else
            ! FRC-based optimal filter
            filtsz = build_glob%clsfrcs%get_filtsz()
            allocate(frc(filtsz),filter(filtsz),source=0.)
            call build_glob%clsfrcs%frc_getter(icls, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
            else
                filter = 1.0
            endif
            ! gaussian filter
            if(trim(params_glob%gauref).eq.'yes')then
                allocate(gaufilter(filtsz),source=0.)
                call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
                ! take the minimum of FRC-based & gaussian filters
                forall(k = 1:filtsz) filter(k) = min(filter(k), gaufilter(k))
                deallocate(gaufilter)
            endif
            call filterrefs(icls, filter)
            deallocate(frc,filter)
        endif
    end subroutine polar_prep2Dref

    !>  \brief  prepares one polar cluster centre image for alignment
    subroutine polar_prep3Dref( icls )
        integer, intent(in) :: icls
        real, allocatable   :: frc(:), filter(:), gaufilter(:)
        integer :: filtsz, k
        ! Filtering
        if( params_glob%l_ml_reg )then
            ! no filtering, not supported yet
        else
            filtsz = build_glob%clsfrcs%get_filtsz()
            allocate(filter(filtsz),source=1.)
            ! FRC-based optimal filter
            if(trim(params_glob%frcref).eq.'yes')then
                allocate(frc(filtsz),source=0.)
                call build_glob%clsfrcs%frc_getter(icls, frc)
                if( any(frc > 0.143) )then
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
                endif
                deallocate(frc)
            endif
            ! optional gaussian filter
            if(trim(params_glob%gauref).eq.'yes')then
                allocate(gaufilter(filtsz),source=0.)
                call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
                ! take the minimum of FRC-based & gaussian filters: relevant in 2.5D?
                forall(k = 1:filtsz) filter(k) = min(filter(k), gaufilter(k))
                deallocate(gaufilter)
            endif
            call filterrefs(icls, filter)
            deallocate(filter)
        endif
    end subroutine polar_prep3Dref

    ! DESTRUCTOR

    subroutine polar_cavger_kill
        if( allocated(pfts_even) )then
            deallocate(pfts_even,pfts_odd,ctf2_even,ctf2_odd,pfts_merg,eo_pops,prev_eo_pops)
        endif
        smpd       = 0.
        ncls       = 0
        nrots      = 0
        kfromto(2) = 0
        pftsz      = 0
        l_comlin   = .false.
    end subroutine polar_cavger_kill

    ! PRIVATE UTILITIES

    ! produces y-mirror of real (reciprocal) matrix 
    pure subroutine mirror_ctf2( ctf2in, ctf2out )
        real(dp), intent(in)    :: ctf2in(pftsz,kfromto(1):kfromto(2))
        real(dp), intent(inout) :: ctf2out(pftsz,kfromto(1):kfromto(2))
        integer :: i,j
        ctf2out(1,:) = ctf2in(1,:)
        do i = 2,pftsz/2
            j = pftsz-i+2
            ctf2out(i,:) = ctf2in(j,:)
            ctf2out(j,:) = ctf2in(i,:)
        enddo
        i = pftsz/2 + 1
        ctf2out(i,:) = ctf2in(i,:)
    end subroutine mirror_ctf2

    ! produces y-mirror of complex (reciprocal) matrix 
    pure subroutine mirror_pft( pftin, pftout )
        complex(dp), intent(in)    :: pftin(pftsz,kfromto(1):kfromto(2))
        complex(dp), intent(inout) :: pftout(pftsz,kfromto(1):kfromto(2))
        integer :: i,j
        pftout(1,:) = conjg(pftin(1,:))
        do i = 2,pftsz/2
            j = pftsz-i+2
            pftout(i,:) = pftin(j,:)
            pftout(j,:) = pftin(i,:)
        enddo
        i = pftsz/2 + 1
        pftout(i,:) = pftin(i,:)
    end subroutine mirror_pft

    ! Private utility
    pure subroutine safe_norm( Mnum, Mdenom, Mout )
        complex(dp), intent(in)    :: Mnum(pftsz,kfromto(1):kfromto(2))
        real(dp),    intent(inout) :: Mdenom(pftsz,kfromto(1):kfromto(2))
        complex(dp), intent(inout) :: Mout(pftsz,kfromto(1):kfromto(2))
        logical  :: msk(pftsz)
        real(dp) :: avg, t
        integer  :: k
        do k = kfromto(1),kfromto(2)
            msk = Mdenom(:,k) > DSMALL
            avg = sum(Mdenom(:,k),mask=msk) / real(count(msk),dp)
            t   = avg/50.d0
            where((Mdenom(:,k) < t).and.msk) Mdenom(:,k) = Mdenom(:,k) + t
        enddo
        where( Mdenom > DSMALL)
            Mout = Mnum / Mdenom
        elsewhere
            Mout = DCMPLX_ZERO
        end where
    end subroutine safe_norm

    ! Returns complex and ctf2 polar lines given ref and rotational indices
    pure subroutine get_line( ref, rot, even, pftline, ctf2line )
        integer,     intent(in)    :: ref, rot
        logical,     intent(in)    :: even
        complex(dp), intent(out)   :: pftline(kfromto(1):kfromto(2))
        real(dp),    intent(out)   :: ctf2line(kfromto(1):kfromto(2))
        integer :: irot
        if( rot >  nrots )then
            irot = rot - nrots
        else
            irot = rot
        endif
        if( even )then
            if( irot < 1 )then
                irot    = irot + pftsz
                pftline = conjg(pfts_even(irot,:,ref))
            elseif( irot > pftsz )then
                irot    = irot - pftsz
                pftline = conjg(pfts_even(irot,:,ref))
            else
                pftline = pfts_even(irot,:,ref)
            endif
            ctf2line = ctf2_even(irot,:,ref)
        else
            if( irot < 1 )then
                irot    = irot + pftsz
                pftline = conjg(pfts_odd(irot,:,ref))
            elseif( irot > pftsz )then
                irot    = irot - pftsz
                pftline = conjg(pfts_odd(irot,:,ref))
            else
                pftline = pfts_odd(irot,:,ref)
            endif
            ctf2line = ctf2_odd(irot,:,ref)
        endif
    end subroutine get_line

    !>  \brief  Filter references
    subroutine filterrefs( icls, filter )
        integer, intent(in) :: icls
        real,    intent(in) :: filter(:)
        integer :: n, k
        n = size(filter)
        if( n < kfromto(2) )then
            THROW_HARD('Incompatible filter size!; polar_cavger_filterref')
        endif
        do k = kfromto(1),kfromto(2)
            pfts_merg(:,k,icls) = filter(k) * pfts_merg(:,k,icls)
            pfts_even(:,k,icls) = filter(k) * pfts_even(:,k,icls)
            pfts_odd(:,k,icls)  = filter(k) * pfts_odd(:,k,icls)
        enddo
    end subroutine filterrefs

    subroutine calc_frc( pft1, pft2, n, frc )
        complex(dp), intent(in)    :: pft1(pftsz,kfromto(1):kfromto(2)), pft2(pftsz,kfromto(1):kfromto(2))
        integer,     intent(in)    :: n
        real(sp),    intent(inout) :: frc(1:n)
        real(dp) :: var1, var2, denom
        integer  :: k
        frc(1:kfromto(1)-1) = 0.999
        do k = kfromto(1), kfromto(2)
            var1  = sum(csq_fast(pft1(:,k)))
            var2  = sum(csq_fast(pft2(:,k)))
            if( (var1>DTINY) .and. (var2>DTINY) )then
                denom  = sqrt(var1) * sqrt(var2)
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / denom, sp)
            else
                frc(k) = 0.0
            endif
        enddo
        if( kfromto(2) < n ) frc(kfromto(2)+1:) = 0.0
    end subroutine calc_frc

    ! Format for PFT I/O
    ! First  integer: PFTSZ
    ! Second integer: KFROMTO(1)
    ! Third  integer: KFROMTO(2)
    ! Fourth integer: NCLS
    ! input/ouput in kind=dp but read/written in kind=sp
    subroutine write_pft_array( array, fname )
        complex(dp),      intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        character(len=*), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) cmplx(array,kind=sp)
        call fclose(funit)
    end subroutine write_pft_array

    subroutine write_ctf2_array( array, fname )
        real(dp),         intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        character(len=*), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) real(array,kind=sp)
        call fclose(funit)
    end subroutine write_ctf2_array

    subroutine read_pft_array( fname, array )
        character(len=*),         intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        complex(sp), allocatable :: tmp(:,:,:)
        integer :: dims(4), funit,io_stat, k
        logical :: samedims
        if( .not.file_exists(trim(fname)) ) THROW_HARD(trim(fname)//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//trim(fname), io_stat)
        read(unit=funit,pos=1) dims
        if( .not.allocated(array) )then
            allocate(array(pftsz,kfromto(1):kfromto(2),ncls))
        endif
        samedims = all(dims == [pftsz, kfromto(1), kfromto(2), ncls])
        if( samedims )then
            allocate(tmp(pftsz,kfromto(1):kfromto(2),ncls))
            read(unit=funit, pos=(sizeof(dims)+1)) tmp
            array(:,:,:) = cmplx(tmp(:,:,:),kind=dp)
        else
            if( pftsz /= dims(1) )then
                THROW_HARD('Incompatible PFT size in '//trim(fname)//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//trim(fname)//': '//int2str(ncls)//' vs '//int2str(dims(4)))
            endif
            allocate(tmp(dims(1),dims(2):dims(3),dims(4)))
            read(unit=funit, pos=(sizeof(dims)+1)) tmp
            do k = kfromto(1),kfromto(2)
                if( (k >= dims(2)) .and. (k <= dims(3)) )then
                    array(:,k,:) = cmplx(tmp(:,k,:),kind=dp)    ! from stored array
                else
                    array(:,k,:) = 0.d0                         ! pad with zeros
                endif
            enddo
        endif
        deallocate(tmp)
        call fclose(funit)
    end subroutine read_pft_array

    subroutine read_ctf2_array( fname, array )
        character(len=*),      intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        real(sp), allocatable :: tmp(:,:,:)
        integer :: dims(4), funit,io_stat, k
        logical :: samedims
        if( .not.file_exists(trim(fname)) ) THROW_HARD(trim(fname)//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//trim(fname), io_stat)
        read(unit=funit,pos=1) dims
        if( .not.allocated(array) )then
            allocate(array(pftsz,kfromto(1):kfromto(2),ncls))
        endif
        samedims = all(dims == [pftsz, kfromto(1), kfromto(2), ncls])
        if( samedims )then
            allocate(tmp(pftsz,kfromto(1):kfromto(2),ncls))
            read(unit=funit, pos=(sizeof(dims)+1)) tmp
            array(:,:,:) = real(tmp(:,:,:),dp)
        else
            if( pftsz /= dims(1) )then
                THROW_HARD('Incompatible real array size in '//trim(fname)//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//trim(fname)//': '//int2str(ncls)//' vs '//int2str(dims(4)))
            endif
            allocate(tmp(dims(1),dims(2):dims(3),dims(4)))
            read(unit=funit, pos=(sizeof(dims)+1)) tmp
            do k = kfromto(1),kfromto(2)
                if( (k >= dims(2)) .or. (k <= dims(3)) )then
                    array(:,k,:) = real(tmp(:,k,:),dp)  ! from stored array
                else
                    array(:,k,:) = 0.d0                 ! pad with zeros
                endif
            enddo
        endif
        deallocate(tmp)
        call fclose(funit)
    end subroutine read_ctf2_array

    subroutine polar_cavger_dims_from_header( fname, pftsz_here, kfromto_here, ncls_here )
        character(len=*), intent(in)    :: fname
        integer,          intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
        integer :: dims(4), funit, io_stat
        if( .not.file_exists(trim(fname)) ) THROW_HARD(trim(fname)//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('dims_from_header; fopen failed: '//trim(fname), io_stat)
        read(unit=funit,pos=1) dims
        call fclose(funit)
        pftsz_here   = dims(1)
        kfromto_here = dims(2:3)
        ncls_here    = dims(4)
    end subroutine polar_cavger_dims_from_header

    ! TEST UNIT

    subroutine test_polarops
        use simple_cmdline,    only: cmdline
        use simple_parameters, only: parameters
        integer,     parameter :: N=128
        integer,     parameter :: NIMGS=200
        integer,     parameter :: NCLS=5
        type(image)            :: tmpl_img, img, cavgs(NCLS)
        type(cmdline)          :: cline
        type(polarft_corrcalc) :: pftcc
        type(parameters)       :: p
        type(builder)          :: b
        real    :: ang, shift(2), shifts(2,NIMGS)
        integer :: pinds(NIMGS), i, eo, icls
        ! dummy structure
        call tmpl_img%soft_ring([N,N,1], 1., 8.)
        call tmpl_img%fft
        call tmpl_img%shift2Dserial([ 8.,-16.])
        call img%soft_ring([N,N,1], 1., 12.)
        call img%fft
        call img%shift2Dserial([ 32., 0.])
        call tmpl_img%add(img)
        call img%soft_ring([N,N,1], 1., 16.)
        call img%fft
        call img%shift2Dserial([ -16., 8.])
        call tmpl_img%add(img)
        call img%soft_ring([N,N,1], 1., 32.)
        call img%fft
        call tmpl_img%add(img)
        call tmpl_img%ifft
        call tmpl_img%write('template.mrc')
        ! init of options & parameters
        call cline%set('prg',    'xxx')
        call cline%set('objfun', 'cc')
        call cline%set('smpd',   1.0)
        call cline%set('box',    N)
        call cline%set('ctf',    'no')
        call cline%set('oritype','ptcl2D')
        call cline%set('ncls',    NCLS)
        call cline%set('nptcls',  NIMGs)
        call cline%set('lp',      3.)
        call cline%set('nthr',    8)
        call cline%set('mskdiam', real(N)/2-10.)
        ! Calculators
        call b%init_params_and_build_strategy2D_tbox(cline, p)
        call pftcc%new(NCLS, [1,NIMGS], p%kfromto)
        pinds = (/(i,i=1,NIMGS)/)
        call b%img_crop_polarizer%init_polarizer(pftcc, p%alpha)
        do i = 1,NIMGS
            shift = 10.*[ran3(), ran3()] - 5.
            ! ang   = 360. * ran3()
            ang   = 0.
            eo    = 0
            if( .not.is_even(i) ) eo = 1
            icls  = ceiling(ran3()*4.)
            call img%copy_fast(tmpl_img)
            call img%fft
            call img%shift2Dserial(-shift)
            call img%ifft
            call img%rtsq(ang, 0.,0.)
            call img%add_gauran(2.)
            call img%write('rotimgs.mrc', i)
            call img%fft
            call b%spproj_field%set_euler(i, [0.,0.,ang])
            call b%spproj_field%set_shift(i, shift)
            call b%spproj_field%set(i,'w',1.0)
            call b%spproj_field%set(i,'state',1)
            call b%spproj_field%set(i,'class', icls)
            call b%spproj_field%set(i,'eo',eo)
            shifts(:,i) = -shift
            call b%img_crop_polarizer%polarize(pftcc, img, i, isptcl=.true., iseven=eo==0, mask=b%l_resmsk)
        enddo
        call polar_cavger_new(pftcc,.false.)
        call polar_cavger_update_sums(NIMGS, pinds, b%spproj, pftcc, shifts)
        call polar_cavger_merge_eos_and_norm
        call polar_cavger_calc_and_write_frcs_and_eoavg(FRCS_FILE, cline)
        ! write
        call polar_cavger_write('cavgs_even.bin', 'even')
        call polar_cavger_write('cavgs_odd.bin',  'odd')
        call polar_cavger_write('cavgs.bin',      'merged')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'even')
        call write_cavgs(cavgs, 'cavgs_even.mrc')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'odd')
        call write_cavgs(cavgs, 'cavgs_odd.mrc')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'merged')
        call write_cavgs(cavgs, 'cavgs_merged.mrc')
        call polar_cavger_kill
        ! read & write again
        call polar_cavger_new(pftcc,.false.)
        call polar_cavger_read('cavgs_even.bin', 'even')
        call polar_cavger_read('cavgs_odd.bin',  'odd')
        call polar_cavger_read('cavgs.bin',      'merged')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'even')
        call write_cavgs(cavgs, 'cavgs2_even.mrc')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'odd')
        call write_cavgs(cavgs, 'cavgs2_odd.mrc')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'merged')
        call write_cavgs(cavgs, 'cavgs2_merged.mrc')
        call polar_cavger_kill
    end subroutine test_polarops

end module simple_polarops
