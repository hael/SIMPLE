module simple_polarops
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_builder,          only: builder, build_glob
use simple_image,            only: image
use simple_parameters,       only: params_glob
use simple_polarft_calc,     only: polarft_calc, pftc_glob
use simple_sp_project,       only: sp_project
use simple_strategy2D_utils
implicit none

! Restoration
public :: polar_cavger_new, polar_cavger_dims_from_header, polar_cavger_update_sums
public :: polar_cavger_calc_and_write_frcs_and_eoavg
public :: polar_cavger_merge_eos_and_norm2D, polar_cavger_merge_eos_and_norm
public :: polar_cavger_assemble_sums_from_parts, polar_cavger_set_ref_pftc
public :: polar_cavger_zero_pft_refs, polar_cavger_kill
! I/O
public :: polar_cavger_refs2cartesian, polar_cavger_write_cartrefs, polar_cavger_writeall_pftcrefs
public :: polar_cavger_write, polar_cavger_writeall, polar_cavger_writeall_cartrefs
public :: polar_cavger_readwrite_partial_sums, polar_cavger_read_all
! Book-keeping
public :: polar_cavger_gen2Dclassdoc, polar_cavger_calc_pops
! Alignment
public :: polar_prep2Dref, polar_prep3Dref
! Public Utils
public :: center_3Dpolar_refs
private
#include "simple_local_flags.inc"

complex(dp), allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:) ! PFTs arrays
real(dp),    allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:)                   ! PFT-size CTF2 arrays
integer,     allocatable :: prev_eo_pops(:,:), eo_pops(:,:)                     ! Class populations
real                     :: smpd       = 0.                                     ! Pixel size
integer                  :: ncls       = 0                                      ! # classes
integer                  :: kfromto(2) = 0                                      ! Resolution range
integer                  :: pftsz      = 0                                      ! Size of PFT in pftc along rotation dimension
integer                  :: nrots      = 0                                      ! # of in-plane rotations; =2*pftsz
logical                  :: l_comlin   = .false.                                ! Whether we operate in 2D/3D(=added common lines)

contains

    !> Module initialization
    subroutine polar_cavger_new( pftc, comlin, nrefs )
        class(polarft_calc), intent(in) :: pftc
        logical,                 intent(in) :: comlin
        integer,       optional, intent(in) :: nrefs
        call polar_cavger_kill
        nrots    = pftc%get_nrots()
        pftsz    = pftc%get_pftsz()
        kfromto  = pftc%get_kfromto()
        l_comlin = comlin
        smpd     = params_glob%smpd
        if( present(nrefs) )then
            ncls = nrefs
        else
            ncls = pftc%get_nrefs()
        endif
        allocate(prev_eo_pops(2,ncls), eo_pops(2,ncls), source=0)
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

    subroutine polar_cavger_set_ref_pftc( icls, which, pftc )
        integer,                 intent(in)    :: icls
        character(len=*),        intent(in)    :: which
        class(polarft_calc), intent(inout) :: pftc
        select case(trim(which))
        case('merged')
            call pftc%set_ref_pft(icls, cmplx(pfts_merg(:,:,icls),kind=sp), .true.)
        case('even')
            call pftc%set_ref_pft(icls, cmplx(pfts_even(:,:,icls),kind=sp), .true.)
        case('odd')
            call pftc%set_ref_pft(icls, cmplx(pfts_odd(:,:,icls),kind=sp), .false.)
        end select
    end subroutine polar_cavger_set_ref_pftc

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
    subroutine polar_cavger_update_sums( nptcls, pinds, spproj, pftc, incr_shifts, is3D )
        use simple_euclid_sigma2, only: eucl_sigma2_glob
        integer,                         intent(in)    :: nptcls
        integer,                         intent(in)    :: pinds(nptcls)
        class(sp_project),               intent(inout) :: spproj
        class(polarft_calc), target, intent(inout) :: pftc
        real,                  optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,               optional, intent(in)    :: is3d
        class(oris), pointer :: spproj_field
        complex(sp), pointer :: pptcls(:,:,:), rptcl(:,:)
        real(sp),    pointer :: pctfmats(:,:,:), rctf(:,:)
        real(dp),    pointer :: rctf2(:,:)
        real(dp) :: w
        real     :: sigma2(kfromto(1):kfromto(2)), incr_shift(2)
        integer  :: eopops(2,ncls), i, icls, iptcl, irot, k
        logical  :: l_ctf, l_even, l_3D, l_shift
        l_3D = .false.
        if( present(is3D) ) l_3D = is3D
        l_shift = .false.
        if( present(incr_shifts) ) l_shift = .true.
        ! retrieve particle info & pointers
        call spproj%ptr2oritype(params_glob%oritype, spproj_field)
        call pftc%get_ptcls_ptr(pptcls)
        l_ctf = pftc%is_with_ctf()
        if( l_ctf ) call pftc%get_ctfmats_ptr(pctfmats)
        ! update classes
        eopops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static) reduction(+:eopops)&
        !$omp private(i,iptcl,w,l_even,icls,irot,incr_shift,rptcl,rctf,rctf2,k,sigma2)
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
            irot = pftc%get_roind_fast(spproj_field%e3get(iptcl))
            if( l_shift )then
                incr_shift = incr_shifts(:,i)
                ! weighted restoration
                if( any(abs(incr_shift) > 1.e-6) ) call pftc%shift_ptcl(iptcl, -incr_shift)
            endif
            call pftc%get_work_pft_ptr(rptcl)
            ! Particle rotation
            call pftc%rotate_pft(pptcls(:,:,i), irot, rptcl)
            ! Particle weight
            rptcl = real(w) * rptcl
            ! Particle ML regularization
            if( params_glob%l_ml_reg )then
                sigma2 = eucl_sigma2_glob%sigma2_noise(kfromto(1):kfromto(2),iptcl)
                do k = kfromto(1),kfromto(2)
                    rptcl(:,k) = rptcl(:,k) / sigma2(k)
                enddo
            endif
            ! Array updates
            if( l_ctf )then
                call pftc%get_work_rpft_ptr(rctf)
                call pftc%get_work_rpft8_ptr(rctf2)
                ! weighted CTF2
                call pftc%rotate_pft(pctfmats(:,:,i), irot, rctf)
                rctf2 = w * real(rctf,kind=dp)**2
                rptcl = rptcl * rctf    ! PhFlip(X).|CTF|
                ! CTF2 ML regularization
                if( params_glob%l_ml_reg )then
                    do k = kfromto(1),kfromto(2)
                        rctf2(:,k) = rctf2(:,k) / real(sigma2(k),dp)
                    enddo
                endif
                if( l_even )then
                    !$omp critical
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + rctf2
                    !$omp end critical
                else
                    !$omp critical
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + rctf2
                    !$omp end critical
                endif
            else
                if( params_glob%l_ml_reg )then
                    ! CTF2=1 & ML regularization
                    call pftc%get_work_rpft8_ptr(rctf2)
                    do k = kfromto(1),kfromto(2)
                        rctf2(:,k) = w / real(sigma2(k),dp)
                    enddo
                    if( l_even )then
                        !$omp critical
                        pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + rctf2
                        !$omp end critical
                    else
                        !$omp critical
                        pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + rctf2
                        !$omp end critical
                    endif
                else
                    if( l_even )then
                        !$omp critical
                        pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w
                        !$omp end critical
                    else
                        !$omp critical
                        pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w
                        !$omp end critical
                    endif
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

    !>  \brief  Restores 2D class-averages only
    subroutine polar_cavger_merge_eos_and_norm2D
        complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2))
        integer     :: icls, eo_pop(2), pop
        select case(trim(params_glob%ref_type))
        case('polar_cavg')
            ! all good
        case DEFAULT
            THROW_HARD('polar_cavger_merge_eos_and_norm2D only for 2D cavgs restoration')
        end select
        pfts_merg = DCMPLX_ZERO
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,eo_pop,pop,pft,ctf2)
        do icls = 1,ncls
            eo_pop = prev_eo_pops(:,icls) + eo_pops(:,icls)
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
                    call safe_norm(pft, ctf2, pfts_merg(:,:,icls))
                endif
                if(eo_pop(1) > 1)then
                    pft = pfts_even(:,:,icls)
                    call safe_norm(pft, ctf2_even(:,:,icls), pfts_even(:,:,icls))
                endif
                if(eo_pop(2) > 1)then
                    pft = pfts_odd(:,:,icls)
                    call safe_norm(pft, ctf2_odd(:,:,icls), pfts_odd(:,:,icls))
                endif
            endif
        end do
        !$omp end parallel do
    end subroutine polar_cavger_merge_eos_and_norm2D

    !>  \brief  Restores 3D slices
    subroutine polar_cavger_merge_eos_and_norm( reforis, cl_weight )
        use simple_class_frcs
        type(oris),           intent(in) :: reforis
        real,       optional, intent(in) :: cl_weight
        type(class_frcs)   :: cavg2clin_frcs
        real,  allocatable :: cavg_clin_frcs(:,:,:)
        complex(dp) :: pfts_clin_even(pftsz,kfromto(1):kfromto(2),ncls)
        complex(dp) :: pfts_clin_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: ctf2_clin_even(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: ctf2_clin_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: fsc(kfromto(1):kfromto(2)), clw
        clw = 1.d0
        select case(trim(params_glob%ref_type))
        case('comlin_hybrid')
            ! 2.5D: cavgs + variable CLs contribution
            ! Common-line contribution weight
            if( present(cl_weight) ) clw = min(max(0.d0,real(cl_weight,dp)),1.d0)
            if( clw > 1.d-6 )then
                ! Mirroring slices
                call mirror_slices(reforis, build_glob%pgrpsyms)
                ! Common-lines conribution
                call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
                &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            endif
        case('comlin_noself', 'comlin')
            ! Mirroring slices
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! Common-lines conribution
            call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
            &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
        case DEFAULT
            THROW_HARD('Invalid REF_TYPE in polar_cavger_merge_eos_and_norm')
        end select
        ! ML-regularization
        if( params_glob%l_ml_reg ) call add_invtausq2rho
        ! Restoration of references
        pfts_merg = DCMPLX_ZERO
        select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                if( clw > 1.d-6 )then
                    call calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call restore_cavgs_comlins(clw)
                else
                    call mirr_and_calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call polar_cavger_merge_eos_and_norm2D
                endif
                call cavg2clin_frcs%write(string('frcs_cavg2clin'//BIN_EXT))
            case('comlin_noself')
                call restore_comlins
            case('comlin')
                call calc_cavg_comlin_frcs( cavg2clin_frcs )
                call cavg2clin_frcs%write(string('frcs_cavg2clin'//BIN_EXT))
                call restore_cavgs_comlins( 1.d0 )
        end select
        ! cleanup
        call cavg2clin_frcs%kill
        if( allocated(cavg_clin_frcs) ) deallocate(cavg_clin_frcs)
      contains

        subroutine add_invtausq2rho
            logical, parameter :: DEBUG = .false.
            complex(dp) :: even(pftsz,kfromto(1):kfromto(2)), odd(pftsz,kfromto(1):kfromto(2))
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: vare(kfromto(1):kfromto(2)), varo(kfromto(1):kfromto(2))
            real(dp)    :: sig2e(kfromto(1):kfromto(2)), sig2o(kfromto(1):kfromto(2))
            real(dp)    :: ssnr(kfromto(1):kfromto(2)), tau2(kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)), cc, fudge, invtau2, s, a,b
            integer     :: icls, k, kstart, p
            ! FSC per slice contribution
            fsc  = 0.d0
            vare = 0.d0; varo = 0.d0
            sig2e = 0.d0; sig2o = 0.d0
            select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                THROW_HARD('Not supported yet')
            case('comlin_noself')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo,sig2e,sig2o)
                do icls = 1,ncls/2
                    ! e/o restoration
                    pft  = pfts_clin_even(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls)
                    call safe_norm(pft, ctf2, even)
                    pft  = pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = kfromto(1),kfromto(2)
                        fsc(k)   = fsc(k)   + sum(real(even(:,k) * conjg(odd(:,k)), dp))
                        vare(k)  = vare(k)  + sum(real(even(:,k) * conjg(even(:,k)),dp))
                        varo(k)  = varo(k)  + sum(real(odd(:,k)  * conjg(odd(:,k)), dp))
                        sig2e(k) = sig2e(k) + sum(ctf2_clin_even(:,k,icls))
                        sig2o(k) = sig2o(k) + sum(ctf2_clin_odd(:,k,icls))
                    enddo
                enddo
                !$omp end parallel do
            case('comlin')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo,sig2e,sig2o)
                do icls = 1,ncls/2
                    ! e/o restoration
                    pft  = pfts_even(:,:,icls) + pfts_clin_even(:,:,icls)
                    ctf2 = ctf2_even(:,:,icls) + ctf2_clin_even(:,:,icls)
                    call safe_norm(pft, ctf2, even)
                    pft  = pfts_odd(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_odd(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = kfromto(1),kfromto(2)
                        fsc(k)   = fsc(k)   + real(sum(even(:,k) * conjg(odd(:,k))) ,dp)
                        vare(k)  = vare(k)  + real(sum(even(:,k) * conjg(even(:,k))),dp)
                        varo(k)  = varo(k)  + real(sum(odd(:,k)  * conjg(odd(:,k))) ,dp)
                        sig2e(k) = sig2e(k) + sum(ctf2_clin_even(:,k,icls))
                        sig2o(k) = sig2o(k) + sum(ctf2_clin_odd(:,k,icls))
                    enddo
                enddo
                !$omp end parallel do
            end select
            ! Variances
            vare = vare * varo
            ! FSC
            where( vare > DTINY )
                fsc = fsc / sqrt(vare)
            elsewhere
                fsc = 0.d0
            end where
            ! SSNR
            fudge = real(params_glob%tau,dp)
            do k = kfromto(1),kfromto(2)
                cc = max(0.001d0,min(0.999d0,fsc(k)))
                ssnr(k) = fudge * cc / (1.d0 - cc)
            enddo
            ! Add Tau2 inverse to denominators
            ! because signal assumed infinite at very low resolution there is no addition
            kstart = max(6, calc_fourier_index(params_glob%hp, params_glob%box_crop, params_glob%smpd_crop))
            ! Even
            where( sig2e > DTINY )
                sig2e = real(ncls*pftsz,dp) / sig2e
            elsewhere
                sig2e = 0.d0
            end where
            tau2 = ssnr * sig2e
            if( DEBUG )then
                do k = kfromto(1),kfromto(2)
                    cc = max(0.001d0,min(0.999d0,fsc(k)))
                    s = cc / (1.d0 - cc)
                    print *, k, real(s), real(fsc(k)),&
                    &real(fudge**2*s/(fudge**2*s+1.d0)), tau2(k), 1.d0 / (fudge * tau2(k))
                enddo
                do k = kfromto(1),kfromto(2)
                    if( tau2(k) > DTINY )then
                        invtau2 = 1.d0 / (fudge * tau2(k))
                        a = sum(sqrt(real(pfts_clin_even(:,k,:)*conjg(pfts_clin_even(:,k,:))))) / real(pftsz*ncls,dp)
                        b =      sum(ctf2_clin_even(:,k,:))                                     / real(pftsz*ncls,dp)
                        print *,k,(a/(b+invtau2)) / (a/b)
                    endif
                enddo
            endif
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,kfromto(2)
                if( tau2(k) > DTINY )then
                    ! CTF2 <- CTF2 + avgCTF2/(tau*SSNR)
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_even(:,k,:) = ctf2_clin_even(:,k,:) + invtau2
                else
                    do icls = 1,ncls
                        do p = 1,pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_clin_even(p,k,icls))
                            ctf2_clin_even(p,k,icls) = ctf2_clin_even(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
            ! Odd
            where( sig2o > DTINY )
                sig2o = real(ncls*pftsz,dp) / sig2o
            elsewhere
                sig2o = 0.d0
            end where
            tau2 = ssnr * sig2o
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,kfromto(2)
                if( tau2(k) > DTINY )then
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_odd(:,k,:) = ctf2_clin_odd(:,k,:) + invtau2
                else
                    do icls = 1,ncls
                        do p = 1,pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_clin_odd(p,k,icls))
                            ctf2_clin_odd(p,k,icls) = ctf2_clin_odd(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
        end subroutine add_invtausq2rho

        ! Mirror 2D classes, calculate CL contribution & FRCS
        ! Module arrays are untouched on exit
        subroutine mirr_and_calc_cavg_comlin_frcs( frcs )
            class(class_frcs), intent(inout) :: frcs
            complex(dp) :: pfte_backup(pftsz,kfromto(1):kfromto(2),ncls)
            complex(dp) :: pfto_backup(pftsz,kfromto(1):kfromto(2),ncls)
            real(dp)    :: ctf2e_backup(pftsz,kfromto(1):kfromto(2),ncls)
            real(dp)    :: ctf2o_backup(pftsz,kfromto(1):kfromto(2),ncls)
            ! backup classes
            pfte_backup(:,:,:)  = pfts_even(:,:,:); pfto_backup(:,:,:)  = pfts_odd(:,:,:)
            ctf2e_backup(:,:,:) = ctf2_even(:,:,:); ctf2o_backup(:,:,:) = ctf2_odd(:,:,:)
            ! mirror classes belonging to the same slice
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! calculate the per slice CLs
            call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
            &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            ! CLs vs. CLS FRCs
            call calc_cavg_comlin_frcs(frcs)
            ! restores classes
            pfts_even(:,:,:) = pfte_backup(:,:,:);  pfts_odd(:,:,:) = pfto_backup(:,:,:)
            ctf2_even(:,:,:) = ctf2e_backup(:,:,:); ctf2_odd(:,:,:) = ctf2o_backup(:,:,:)
        end subroutine mirr_and_calc_cavg_comlin_frcs

        ! Calculate CL contribution & FRCS from mirrored 2D classes
        ! Module arrays are untouched on exit
        subroutine calc_cavg_comlin_frcs( frcs )
            class(class_frcs), intent(inout) :: frcs
            real, allocatable :: frc(:)
            complex(dp)       :: cavg(pftsz,kfromto(1):kfromto(2)), clin(pftsz,kfromto(1):kfromto(2))
            complex(dp)       :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)          :: ctf2(pftsz,kfromto(1):kfromto(2))
            integer           :: icls, pop, nk
            call frcs%new(ncls, params_glob%box, params_glob%smpd, 1)
            nk = frcs%get_filtsz()
            allocate(frc(1:nk),source=0.)
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,pop,pft,ctf2,cavg,clin,frc)
            do icls = 1,ncls
                pop = sum(prev_eo_pops(:,icls) + eo_pops(:,icls))
                if( pop > 1 )then
                    ! cavg
                    pft  = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                    ctf2 = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                    call safe_norm(pft, ctf2, cavg)
                    ! comlin
                    pft  = pfts_clin_even(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, clin)
                    ! FRC
                    call calc_frc(cavg, clin, nk, frc)
                else
                    frc = 0.0
                endif
                call frcs%set_frc(icls, frc)
            end do
            !$omp end parallel do
            deallocate(frc)
        end subroutine calc_cavg_comlin_frcs

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

        ! Restore common lines contributions only
        subroutine restore_comlins
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2)),pfte(pftsz,kfromto(1):kfromto(2)),pfto(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)),ctf2e(pftsz,kfromto(1):kfromto(2)),ctf2o(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
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
        end subroutine restore_comlins

        ! Restores slices (=cavgs + comlin)
        subroutine restore_cavgs_comlins( clw )
            real(dp), intent(in) :: clw
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2)),pfte(pftsz,kfromto(1):kfromto(2))
            complex(dp) :: pfto(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)),ctf2e(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2o(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            ! Restoration
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,ncls/2
                ! calculate sum of already mirrored class + common-line contribution
                pfte  = pfts_even(:,:,icls) + clw * pfts_clin_even(:,:,icls)
                pfto  = pfts_odd(:,:,icls)  + clw * pfts_clin_odd(:,:,icls)
                ctf2e = ctf2_even(:,:,icls) + clw * ctf2_clin_even(:,:,icls)
                ctf2o = ctf2_odd(:,:,icls)  + clw * ctf2_clin_odd(:,:,icls)
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
        end subroutine restore_cavgs_comlins

    end subroutine polar_cavger_merge_eos_and_norm

    ! Calculate common-lines contributions from all the slices
    subroutine calc_comlin_contrib( ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        logical, parameter :: L_KB = .true.
        type(oris),       intent(in)    :: ref_space
        type(sym),        intent(in)    :: symop
        complex(kind=dp), intent(inout) :: pfts_cl_even(pftsz,kfromto(1):kfromto(2),ncls)
        complex(kind=dp), intent(inout) :: pfts_cl_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(kind=dp),    intent(inout) :: ctf2_cl_even(pftsz,kfromto(1):kfromto(2),ncls)
        real(kind=dp),    intent(inout) :: ctf2_cl_odd(pftsz,kfromto(1):kfromto(2),ncls)
        type(kbinterpol)  :: kbwin
        real, allocatable :: Rsym(:,:,:)
        complex(dp) :: cl_l(kfromto(1):kfromto(2)), cl_r(kfromto(1):kfromto(2)), cl_e(kfromto(1):kfromto(2))
        complex(dp) :: cl_o(kfromto(1):kfromto(2)), pft(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: rl_l(kfromto(1):kfromto(2)), rl_r(kfromto(1):kfromto(2)), rl_e(kfromto(1):kfromto(2))
        real(dp)    :: rl_o(kfromto(1):kfromto(2)), ctf2(pftsz,kfromto(1):kfromto(2)), dd, w, wl, wr, sumw
        real        :: eulers(3),R(3,3,ncls),Rj(3,3),tRi(3,3),psi,drot,d,targ_w,self_w
        integer     :: rotl,rotr, iref, jref, m, self_irot, targ_irot, isym, nsym
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        drot = pftc_glob%get_dang()
        if( L_KB ) kbwin = kbinterpol(1.5, KBALPHA)
        pfts_cl_even = DCMPLX_ZERO; pfts_cl_odd = DCMPLX_ZERO
        ctf2_cl_even = 0.d0; ctf2_cl_odd  = 0.d0
        ! Symmetry rotation matrices
        nsym = symop%get_nsym()
        allocate(Rsym(3,3,nsym),source=0.)
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(iref,jref,m,tRi,Rj,eulers,targ_irot,self_irot,targ_w,self_w)&
        !$omp& private(d,dd,w,wl,wr,sumw,psi,l_rotm,cl_l,cl_r,cl_e,cl_o)&
        !$omp& private(rl_l,rl_r,rl_e,rl_o,pft,ctf2,rotl,rotr,isym)
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
                    if( L_KB )then
                        ! KB
                        ! in plane rotation index of jref slice intersecting iref
                        psi       = 360.0 - eulers(3)
                        targ_irot = pftc_glob%get_roind_fast(psi)
                        d         = psi - pftc_glob%get_rot(targ_irot)
                        if( d > drot ) d = d - 360.0
                        targ_w    = d / drot
                        ! in plane rotation index of iref slice
                        psi       = eulers(1)
                        self_irot = pftc_glob%get_roind_fast(psi)
                        d         = psi - pftc_glob%get_rot(self_irot)
                        if( d > drot ) d = d - 360.0
                        self_w    = d / drot
                        ! intepolate common line in jref-th slice
                        rotl = targ_irot - 1; rotr = targ_irot + 1
                        dd   = real(targ_w,dp)
                        w    = kbwin%apod_dp(dd); wl = kbwin%apod_dp(dd-1.d0); wr = kbwin%apod_dp(dd+1.d0)
                        sumw = wl + w + wr
                        w    = w / sumw; wl = wl / sumw; wr = wr / sumw
                        call get_line(jref, targ_irot, .true., cl_e, rl_e)
                        call get_line(jref, rotl,      .true., cl_l, rl_l)
                        call get_line(jref, rotr,      .true., cl_r, rl_r)
                        cl_e = wl*cl_l + w*cl_e + wr*cl_r
                        rl_e = wl*rl_l + w*rl_e + wr*rl_r
                        call get_line(jref, targ_irot, .false., cl_o, rl_o)
                        call get_line(jref, rotl,      .false., cl_l, rl_l)
                        call get_line(jref, rotr,      .false., cl_r, rl_r)
                        cl_o = wl*cl_l + w*cl_o + wr*cl_r
                        rl_o = wl*rl_l + w*rl_o + wr*rl_r
                        ! extrapolate the common line to iref-th slice
                        rotl = self_irot - 1; rotr = self_irot + 1
                        if( rotr > nrots ) rotr = rotr - nrots
                        dd   = real(self_w,dp)
                        w    = kbwin%apod_dp(dd); wl = kbwin%apod_dp(dd-1.d0); wr = kbwin%apod_dp(dd+1.d0)
                        sumw = wl + w + wr
                        w    = w / sumw; wl = wl / sumw; wr = wr / sumw
                        ! leftmost line
                        call extrapolate_line(iref, rotl,      wl, cl_e, cl_o, rl_e, rl_o)
                        ! nearest line
                        call extrapolate_line(iref, self_irot, w,  cl_e, cl_o, rl_e, rl_o)
                        ! rightmost line
                        call extrapolate_line(iref, rotr,      wr, cl_e, cl_o, rl_e, rl_o)
                    else
                        ! Linear interpolation
                        ! in plane rotation index of jref slice intersecting iref
                        psi       = 360.0 - eulers(3)
                        targ_irot = pftc_glob%get_roind_fast(psi)
                        d         = psi - pftc_glob%get_rot(targ_irot)
                        if( d > drot ) d = d - 360.0
                        if( d < 0. )then
                            targ_irot = targ_irot - 1
                            if( targ_irot < 1 ) targ_irot = targ_irot + nrots
                            d = d + drot
                        endif
                        targ_w = d / drot
                        ! in plane rotation index of iref slice
                        psi       = eulers(1)
                        self_irot = pftc_glob%get_roind_fast(psi)
                        d         = psi - pftc_glob%get_rot(self_irot)
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
                        call get_line(jref, rotl, .true., cl_e, rl_e)
                        call get_line(jref, rotr, .true., cl_r, rl_r)
                        cl_e = wr*cl_e + wl*cl_r
                        rl_e = wr*rl_e + wl*rl_r
                        call get_line(jref, rotl, .false., cl_o, rl_o)
                        call get_line(jref, rotr, .false., cl_r, rl_r)
                        cl_o = wr*cl_o + wl*cl_r
                        rl_o = wr*rl_o + wl*rl_r
                        ! extrapolate the common line to iref-th slice
                        rotl = self_irot; rotr = rotl + 1
                        if( rotr > nrots ) rotr = rotr - nrots
                        w  = real(self_w,dp)
                        call extrapolate_line(iref, rotl, 1.d0-w, cl_e, cl_o, rl_e, rl_o)
                        call extrapolate_line(iref, rotr,      w, cl_e, cl_o, rl_e, rl_o)
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
                call mirror_pft(conjg(pfts_cl_even(:,:,iref)), pfts_cl_even(:,:,m))
                call mirror_pft(conjg(pfts_cl_odd(:,:,iref)),  pfts_cl_odd(:,:,m))
            else
                call mirror_pft(pfts_cl_even(:,:,iref), pfts_cl_even(:,:,m))
                call mirror_pft(pfts_cl_odd(:,:,iref),  pfts_cl_odd(:,:,m))
            endif
            call mirror_ctf2(ctf2_cl_even(:,:,iref), ctf2_cl_even(:,:,m))
            call mirror_ctf2(ctf2_cl_odd(:,:,iref),  ctf2_cl_odd(:,:,m))
        enddo
        !$omp end do
        !$omp end parallel

    contains

        ! Extrapolate cline, rline to pfts_clin and ctf2_clin
        subroutine extrapolate_line(ref, rot, weight, cle, clo, rle, rlo)
            integer,     intent(in) :: ref, rot
            real(dp),    intent(in) :: weight
            complex(dp), intent(in) :: cle(kfromto(1):kfromto(2)), clo(kfromto(1):kfromto(2))
            real(dp),    intent(in) :: rle(kfromto(1):kfromto(2)), rlo(kfromto(1):kfromto(2))
            integer :: irot
            irot = rot
            if( irot < 1 )then
                irot = irot + pftsz
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * conjg(cle)
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * conjg(clo)
            elseif( irot > pftsz )then
                irot = irot - pftsz
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * conjg(cle)
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * conjg(clo)
            else
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * cle
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * clo
            endif
            ctf2_cl_even(irot,:,ref) = ctf2_cl_even(irot,:,ref) + weight * rle
            ctf2_cl_odd( irot,:,ref) = ctf2_cl_odd( irot,:,ref) + weight * rlo
        end subroutine extrapolate_line

    end subroutine calc_comlin_contrib

    !>  \brief  calculates Fourier ring correlations
    subroutine polar_cavger_calc_and_write_frcs_and_eoavg( fname, cline )
        use simple_cmdline, only: cmdline
        class(string), intent(in) :: fname
        type(cmdline), intent(in) :: cline
        complex(dp), allocatable :: prev_prefs(:,:,:)
        real,        allocatable :: frc(:)
        real(dp) :: ufrac_trec
        integer  :: icls, find, pop, filtsz
        filtsz = fdim(params_glob%box_crop) - 1
        allocate(frc(filtsz),source=0.)
        ! In case nspace/ncls has changed
        if( build_glob%clsfrcs%get_ncls() /= ncls )then
            call build_glob%clsfrcs%new(ncls, params_glob%box_crop, params_glob%smpd_crop, params_glob%nstates)
        endif
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
            call read_pft_array(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_even = ufrac_trec * pfts_even + (1.d0-ufrac_trec) * prev_prefs
            !$omp end parallel workshare
            call read_pft_array(string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_odd  = ufrac_trec * pfts_odd   + (1.d0-ufrac_trec) * prev_prefs
            pfts_merg = 0.5d0 * (pfts_even + pfts_odd)
            !$omp end parallel workshare
            deallocate(prev_prefs)
        endif
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    !>  \brief  Converts the polar references to a cartesian grid
    subroutine polar_cavger_refs2cartesian( pftc, cavgs, which, pfts_in )
        use simple_image
        class(polarft_calc), intent(in)    :: pftc
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
                    phys  = pftc%get_coord(irot,k) + [1.,real(c)]
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
    subroutine polar_cavger_assemble_sums_from_parts( reforis, clin_anneal )
        type(oris), optional, intent(in) :: reforis
        real,       optional, intent(in) :: clin_anneal
        complex(dp), allocatable :: pfte(:,:,:), pfto(:,:,:)
        real(dp),    allocatable :: ctf2e(:,:,:), ctf2o(:,:,:)
        type(string) :: cae, cao, cte, cto
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
        call polar_cavger_merge_eos_and_norm(reforis=reforis, cl_weight=clin_anneal)
    end subroutine polar_cavger_assemble_sums_from_parts

    ! I/O

    ! Writes cavgs PFT array
    subroutine polar_cavger_write( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        select case(which)
            case('even')
                call write_pft_array(pfts_even, fname)
            case('odd')
                call write_pft_array(pfts_odd,  fname)
            case('merged')
                call write_pft_array(pfts_merg, fname)
                ! call pft2img
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_write

    ! writes out |PFTs_MERG| as mrcs
    subroutine pft2img()
        type(image) :: img
        integer :: nk,i,k,icls
        nk = kfromto(2)
        if( .not.is_even(nk) ) nk = nk+1
        call img%new([pftsz,nk,1],1.0)
        do icls = 1,ncls
            img = 0.0
            do i = 1,pftsz
                do k = kfromto(1),kfromto(2)
                    call img%set([i,k,1], real(abs(pfts_merg(i,k,icls))))
                enddo
            enddo
            call img%write(string('pfts_it'//int2str(params_glob%which_iter)//'.mrc'),icls)
        enddo
        call img%kill
    end subroutine pft2img


    ! Writes all cavgs PFT arrays
    subroutine polar_cavger_writeall( tmpl_fname )
        class(string), intent(in) :: tmpl_fname
        call polar_cavger_write(tmpl_fname//'_even'//BIN_EXT,'even')
        call polar_cavger_write(tmpl_fname//'_odd'//BIN_EXT, 'odd')
        call polar_cavger_write(tmpl_fname//BIN_EXT,         'merged')
    end subroutine polar_cavger_writeall

    ! Write references contained in pftc
    subroutine polar_cavger_writeall_pftcrefs( tmpl_fname )
        class(string),  intent(in) :: tmpl_fname
        complex(sp),  pointer :: ptre(:,:,:), ptro(:,:,:)
        call pftc_glob%get_refs_ptr( ptre, ptro )
        pfts_even = cmplx(ptre,kind=dp)
        pfts_odd  = cmplx(ptro,kind=dp)
        call polar_cavger_write(tmpl_fname//'_even'//BIN_EXT,'even')
        call polar_cavger_write(tmpl_fname//'_odd'//BIN_EXT, 'odd')
        call polar_cavger_zero_pft_refs !! removed after writing
        nullify(ptre, ptro)
    end subroutine polar_cavger_writeall_pftcrefs

    ! Converts cavgs PFTS to cartesian grids and writes them
    subroutine polar_cavger_write_cartrefs( pftc, tmpl_fname, which )
        class(polarft_calc), intent(in) :: pftc
        class(string),           intent(in) :: tmpl_fname
        character(len=*),        intent(in) :: which
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(ncls, [params_glob%box_crop, params_glob%box_crop,1], smpd, imgs)
        select case(trim(which))
            case('even','odd')
                call polar_cavger_refs2cartesian( pftc, imgs, trim(which) )
                call write_imgarr(imgs, tmpl_fname//'_'//trim(which)//params_glob%ext%to_char())
            case('merged')
                call polar_cavger_refs2cartesian( pftc, imgs, 'merged' )
                call write_imgarr(imgs, tmpl_fname//params_glob%ext)
        end select
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_write_cartrefs

    ! Converts all cavgs PFTS to cartesian grids and writes them
    subroutine polar_cavger_writeall_cartrefs( pftc, tmpl_fname )
        class(polarft_calc), intent(in) :: pftc
        class(string),           intent(in) :: tmpl_fname
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(ncls, [params_glob%box_crop, params_glob%box_crop,1], smpd, imgs)
        call polar_cavger_refs2cartesian( pftc, imgs, 'even' )
        call write_imgarr(imgs, tmpl_fname//'_even'//params_glob%ext%to_char())
        call polar_cavger_refs2cartesian( pftc, imgs, 'odd' )
        call write_imgarr(imgs, tmpl_fname//'_odd'//params_glob%ext%to_char())
        call polar_cavger_refs2cartesian( pftc, imgs, 'merged' )
        call write_imgarr(imgs, tmpl_fname//params_glob%ext)
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_writeall_cartrefs

    ! Read cavgs PFT array
    subroutine polar_cavger_read( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        select case(which)
            case('even')
                call read_pft_array(fname, pfts_even)
            case('odd')
                call read_pft_array(fname, pfts_odd)
            case('merged')
                call read_pft_array(fname, pfts_merg)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_read

    ! Reads all cavgs PFT arrays
    subroutine polar_cavger_read_all( fname )
        class(string),  intent(in) :: fname
        type(string) :: refs, refs_even, refs_odd, ext
        ext = string('.')//fname2ext(fname)
        if( ext == params_glob%ext )then
            refs = get_fbody(fname, params_glob%ext, separator=.false.)//BIN_EXT
        elseif( ext == BIN_EXT )then
            refs = fname
        else
            THROW_HARD('Unsupported file format: '//ext%to_char())
        endif
        refs_even = get_fbody(refs,BIN_EXT,separator=.false.)//'_even'//BIN_EXT
        refs_odd  = get_fbody(refs,BIN_EXT,separator=.false.)//'_odd'//BIN_EXT
        if( .not. file_exists(refs) )then
            THROW_HARD('Polar references do not exist in cwd: '//refs%to_char())
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
        type(string) :: cae, cao, cte, cto
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
        call cae%kill
        call cao%kill
        call cte%kill
        call cto%kill
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
                if( params_glob%beta > 0.001 )then
                    call fsc2boostfilter(params_glob%beta, filtsz, frc, filter, merged=params_glob%l_lpset)
                else
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
                endif
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
        real, allocatable   :: gaufilter(:)
        integer :: filtsz
        ! Gaussian filter if needed
        if(trim(params_glob%gauref).eq.'yes')then
            filtsz = build_glob%clsfrcs%get_filtsz()
            allocate(gaufilter(filtsz),source=0.)
            call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
            call filterrefs(icls, gaufilter)
            deallocate(gaufilter)
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

    ! PUBLIC UTILITIES

    ! Estimates the center of the volume based on the distribution of
    ! the individual particles in-plane offsets and map the shifts to both
    ! the particles and the references stored in the pftc
    subroutine center_3Dpolar_refs( pftc, algndoc, algnrefs )
        class(polarft_calc), intent(inout) :: pftc
        class(oris),             intent(inout) :: algndoc
        class(oris),             intent(in)    :: algnrefs
        real    :: R(3,3),offset3D(3), offset2D(3)
        integer :: iref
        ! estimate 3D offset from particle alignement parameters
        call algndoc%calc_avg_offset3D(offset3D, state=1)
        ! report 3D offset to particles 2D offsets
        call algndoc%map3dshift22d(-offset3D, state=1)
        ! report 3D offset to alignment references
        !$omp parallel do proc_bind(close) default(shared) private(iref,R,offset2D)
        do iref = 1,ncls
            ! Projection direction rotation matrix
            R = euler2m([algnrefs%e1get(iref), algnrefs%e2get(iref), 0.0])
            ! 3D Shift rotated with respect to projection direction
            offset2D = matmul(R, offset3D)
            ! Apply offset to e/o references
            call pftc%shift_ref(iref, offset2D(1:2))
        enddo
        !$omp end parallel do
    end subroutine center_3Dpolar_refs

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
                denom  = sqrt(var1*var2)
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
        complex(dp),   intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//fname%to_char(),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) cmplx(array,kind=sp)
        call fclose(funit)
    end subroutine write_pft_array

    subroutine write_ctf2_array( array, fname )
        real(dp),      intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//fname%to_char(),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) real(array,kind=sp)
        call fclose(funit)
    end subroutine write_ctf2_array

    subroutine read_pft_array( fname, array )
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        complex(sp), allocatable :: tmp(:,:,:)
        integer :: dims(4), funit,io_stat, k
        logical :: samedims
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//fname%to_char(), io_stat)
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
                THROW_HARD('Incompatible PFT size in '//fname%to_char()//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//fname%to_char()//': '//int2str(ncls)//' vs '//int2str(dims(4)))
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
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        real(sp), allocatable :: tmp(:,:,:)
        integer :: dims(4), funit,io_stat, k
        logical :: samedims
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//fname%to_char(), io_stat)
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
                THROW_HARD('Incompatible real array size in '//fname%to_char()//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//fname%to_char()//': '//int2str(ncls)//' vs '//int2str(dims(4)))
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
        class(string), intent(in)    :: fname
        integer,       intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
        integer :: dims(4), funit, io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('dims_from_header; fopen failed: '//fname%to_char(), io_stat)
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
        type(polarft_calc)     :: pftc
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
        call tmpl_img%write(string('template.mrc'))
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
        call pftc%new(NCLS, [1,NIMGS], p%kfromto)
        pinds = (/(i,i=1,NIMGS)/)
        call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
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
            call img%write(string('rotimgs.mrc'), i)
            call img%fft
            call b%spproj_field%set_euler(i, [0.,0.,ang])
            call b%spproj_field%set_shift(i, shift)
            call b%spproj_field%set(i,'w',1.0)
            call b%spproj_field%set(i,'state',1)
            call b%spproj_field%set(i,'class', icls)
            call b%spproj_field%set(i,'eo',eo)
            shifts(:,i) = -shift
            call b%img_crop_polarizer%polarize(pftc, img, i, isptcl=.true., iseven=eo==0, mask=b%l_resmsk)
        enddo
        call polar_cavger_new(pftc,.false.)
        call polar_cavger_update_sums(NIMGS, pinds, b%spproj, pftc, shifts)
        call polar_cavger_merge_eos_and_norm2D
        call polar_cavger_calc_and_write_frcs_and_eoavg(string(FRCS_FILE), cline)
        ! write
        call polar_cavger_write(string('cavgs_even.bin'), 'even')
        call polar_cavger_write(string('cavgs_odd.bin'),  'odd')
        call polar_cavger_write(string('cavgs.bin'),      'merged')
        call polar_cavger_refs2cartesian(pftc, cavgs, 'even')
        call write_imgarr(cavgs, string('cavgs_even.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'odd')
        call write_imgarr(cavgs, string('cavgs_odd.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'merged')
        call write_imgarr(cavgs, string('cavgs_merged.mrc'))
        call polar_cavger_kill
        ! read & write again
        call polar_cavger_new(pftc,.false.)
        call polar_cavger_read(string('cavgs_even.bin'), 'even')
        call polar_cavger_read(string('cavgs_odd.bin'),  'odd')
        call polar_cavger_read(string('cavgs.bin'),      'merged')
        call polar_cavger_refs2cartesian(pftc, cavgs, 'even')
        call write_imgarr(cavgs, string('cavgs2_even.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'odd')
        call write_imgarr(cavgs, string('cavgs2_odd.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'merged')
        call write_imgarr(cavgs, string('cavgs2_merged.mrc'))
        call polar_cavger_kill
    end subroutine test_polarops

end module simple_polarops
