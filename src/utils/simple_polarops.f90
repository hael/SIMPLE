module simple_polarops
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_builder,           only: builder, build_glob
use simple_parameters,        only: params_glob
use simple_sp_project,        only: sp_project
use simple_image,             only: image
use simple_stack_io,          only: stack_io
! use simple_ctf,               only: ctf
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_euclid_sigma2
use simple_progress
implicit none

public :: polar_cavger_new, polar_cavger_update_sums, polar_cavger_merge_eos_and_norm
public :: polar_cavger_write, polar_cavger_kill
public :: polar_cavger_restore_classes
public :: test_polarops
private
#include "simple_local_flags.inc"


complex(dp), allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_refs(:,:,:)
real(dp),    allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:) 
integer,     allocatable :: prev_eo_pops(:,:), eo_pops(:,:)
real                     :: smpd       = 0.          !< sampling distance
integer                  :: ncls       = 0           !< # classes
integer                  :: kfromto(2) = 0
integer                  :: pftsz = 0

contains

    subroutine polar_cavger_new( pftcc )
        class(polarft_corrcalc), target, intent(in) :: pftcc
        ncls    = params_glob%ncls
        pftsz   = pftcc%get_pftsz()
        kfromto = pftcc%get_kfromto()
        ! dimensions
        smpd          = params_glob%smpd
        allocate(prev_eo_pops(ncls,2), eo_pops(ncls,2), source=0)
        ! Arrays        
        allocate(pfts_even(pftsz,kfromto(1):kfromto(2),ncls),pfts_odd(pftsz,kfromto(1):kfromto(2),ncls),&
            &ctf2_even(pftsz,kfromto(1):kfromto(2),ncls),ctf2_odd(pftsz,kfromto(1):kfromto(2),ncls),&
            &pfts_refs(pftsz,kfromto(1):kfromto(2),ncls))
        call polar_cavger_zero_pft_refs
        pfts_refs = DCMPLX_ZERO
    end subroutine polar_cavger_new

    subroutine polar_cavger_zero_pft_refs
        pfts_even = DCMPLX_ZERO
        pfts_odd  = DCMPLX_ZERO
        ctf2_even = 0.d0
        ctf2_odd  = 0.d0
    end subroutine polar_cavger_zero_pft_refs

    subroutine polar_cavger_update_sums( nptcls, pinds, spproj, pftcc, incr_shifts )
        integer,                         intent(in)    :: nptcls
        integer,                         intent(in)    :: pinds(nptcls)
        class(sp_project),               intent(inout) :: spproj
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        real,                            intent(in)    :: incr_shifts(2,nptcls)
        class(oris), pointer :: spproj_field
        complex(sp), pointer :: pptcls(:,:,:), rptcl(:,:)
        real(sp),    pointer :: pctfmats(:,:,:), rctf(:,:)
        real(dp) :: w
        real     :: incr_shift(2)
        integer  :: i, icls, iptcl, irot
        logical  :: l_ctf, l_even
        ! retrieve particle info & pointers
        call spproj%ptr2oritype(params_glob%oritype, spproj_field)
        l_ctf = pftcc%is_with_ctf()
        call pftcc%get_ptcls_ptr(pptcls)
        if( l_ctf )call pftcc%get_ctfmats_ptr(pctfmats)
        ! update classes
        do i = 1,nptcls
            ! particles parameters
            iptcl = pinds(i)
            if( spproj_field%get_state(iptcl) == 0  ) cycle
            w = real(spproj_field%get(iptcl,'w'),dp)
            if( w < DSMALL ) cycle
            l_even = spproj_field%get_eo(iptcl)==0
            icls   = spproj_field%get_class(iptcl)
            irot   = pftcc%get_roind(spproj_field%e3get(iptcl))
            incr_shift = incr_shifts(:,i)
            ! weighted restoration
            ! if(abs(incr_shift(1)) > 1.e-6 .or. abs(incr_shift(2)) > 1.e-6 )then
            !     call pftcc%shift_ptcl(iptcl, incr_shift)
            ! endif
            call pftcc%get_work_pft_ptr(rptcl)
            call pftcc%rotate_pft(pptcls(:,:,i), irot, rptcl)
            if( l_ctf )then
                call pftcc%get_work_rpft_ptr(rctf)
                call pftcc%rotate_pft(pctfmats(:,:,i), irot, rctf)
                if( l_even )then
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + w * cmplx(rptcl,kind=dp) * real(rctf,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w * real(rctf,kind=dp)**2
                else
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + w * cmplx(rptcl,kind=dp) * real(rctf,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w * real(rctf,kind=dp)**2
                endif
            else
                if( l_even )then
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + w * cmplx(rptcl,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w
                else
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + w * cmplx(rptcl,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w
                endif
            endif
            ! total population
            if( l_even )then
                eo_pops(icls,1) = eo_pops(icls,1) + 1
            else
                eo_pops(icls,2) = eo_pops(icls,2) + 1
            endif
        enddo
        ! cleanup
        nullify(rptcl,rctf,pptcls,pctfmats)
    end subroutine polar_cavger_update_sums

    subroutine polar_cavger_merge_eos_and_norm
        real, parameter :: EPSILON = 0.1
        complex(dp) :: numerator(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: denominator(pftsz,kfromto(1):kfromto(2))
        integer     :: icls, eo_pop(2), pop
        pfts_refs = CMPLX_ZERO
        !$omp parallel do default(shared), schedule(static) proc_bind(close)&
        !$omp private(icls,eo_pop,pop,numerator,denominator)
        do icls=1,ncls
            eo_pop = prev_eo_pops(icls,:) + eo_pops(icls,:) ! eo_pops has to be calculated differently
            pop    = sum(eo_pop)
            if(pop == 0)then
                pfts_even(:,:,icls) = DCMPLX_ZERO
                pfts_odd(:,:,icls)  = DCMPLX_ZERO
                ctf2_even(:,:,icls) = 0.d0
                ctf2_odd(:,:,icls)  = 0.d0
            else
                ! w*CTF**2 density correction
                if(pop > 1)then
                    numerator   = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                    denominator = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                    if( pop <= 5 ) denominator = denominator + real(EPSILON/real(pop),dp)
                    where( abs(denominator) > DSMALL ) pfts_refs(:,:,icls) = numerator / denominator
                endif
                if(eo_pop(1) > 1)then
                    where( abs(ctf2_even(:,:,icls)) > DSMALL ) pfts_even(:,:,icls) = pfts_even(:,:,icls) / ctf2_even(:,:,icls)
                endif
                if(eo_pop(2) > 1)then
                    where( abs(ctf2_odd(:,:,icls)) > DSMALL )  pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  / ctf2_odd(:,:,icls)
                endif
            endif
        end do
        !$omp end parallel do
    end subroutine polar_cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    subroutine polar_cavger_calc_and_write_frcs_and_eoavg( fname )
        use simple_builder, only: build_glob
        character(len=*), intent(in) :: fname
        real, allocatable :: frc(:)
        integer           :: eo_pop(2), icls, find, pop, filtsz
        filtsz = fdim(params_glob%box) - 1
        allocate(frc(filtsz),source=0.)
        !$omp parallel do default(shared) private(icls,frc,find,pop,eo_pop) schedule(static) proc_bind(close)
        do icls = 1,ncls
            eo_pop = prev_eo_pops(icls,:) + eo_pops(icls,:) ! eo_pops has to be calculated differently
            pop    = sum(eo_pop)
            if( pop == 0 )then
                frc = 0.
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
            else
                ! calculate FRC
                call calc_frc(pfts_even(:,:,icls), pfts_odd(:,:,icls), filtsz, frc)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                if( find >= kfromto(1) )then
                    pfts_even(:,kfromto(1):find,icls) = pfts_refs(:,kfromto(1):find,icls) / 2.d0
                    pfts_odd(:,kfromto(1):find,icls)  = pfts_even(:,kfromto(1):find,icls)
                endif
            endif
        end do
        !$omp end parallel do
        ! write FRCs
        call build_glob%clsfrcs%write(fname)
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    subroutine polar_cavger_refs2cartesian( pftcc, cavgs, which )
        use simple_image
        class(polarft_corrcalc), intent(in)    :: pftcc
        type(image),             intent(inout) :: cavgs(ncls)
        character(len=*),        intent(in)    :: which
        complex(dp), allocatable :: cmat(:,:)
        real(dp),    allocatable :: norm(:,:)
        complex :: fc
        real    :: phys(2), dh,dk
        integer :: k,c,irot,physh,physk,box,icls,sz
        box = params_glob%box
        c   = box/2+1
        allocate(cmat(c,box),norm(c,box))
        do icls = 1, ncls
            cmat = DCMPLX_ZERO
            norm = 0.d0
            select case(trim(which))
                case('even')
                    do irot = 1,pftsz
                        do k = kfromto(1),kfromto(2)
                            phys  = pftcc%get_coord(irot,k)
                            fc    = cmplx(pfts_even(irot,k,icls),kind=dp)
                            phys  = phys + [1.,real(c)]
                            physh = floor(phys(1))
                            physk = floor(phys(2))
                            dh = phys(1) - real(physh) 
                            dk = phys(2) - real(physk)
                            if( physh > 0 .and. physh <= c )then
                                if( physk <= box )then
                                    cmat(physh,physk) = cmat(physh,physk) + (1.-dh)*(1-dk)*fc
                                    norm(physh,physk) = norm(physh,physk) + (1.-dh)*(1-dk)
                                    if( physk+1 <= box )then
                                        cmat(physh,physk+1) = cmat(physh,physk+1) + (1.-dh)*dk*fc
                                        norm(physh,physk+1) = norm(physh,physk+1) + (1.-dh)*dk
                                    endif
                                endif
                            endif
                            physh = physh + 1
                            if( physh > 0 .and. physh <= c )then
                                if( physk <= box )then
                                    cmat(physh,physk) = cmat(physh,physk) + dh*(1-dk)*fc
                                    norm(physh,physk) = norm(physh,physk) + dh*(1-dk)
                                    if( physk+1 <= box )then
                                        cmat(physh,physk+1) = cmat(physh,physk+1) + dh*dk*fc
                                        norm(physh,physk+1) = norm(physh,physk+1) + dh*dk
                                    endif
                                endif
                            endif
                        end do
                    end do
                case('odd')
                    do irot = 1,pftsz
                        do k = kfromto(1),kfromto(2)
                            phys = pftcc%get_coord(irot,k)
                            physh = nint(phys(1)) + 1
                            physk = nint(phys(2)) + c
                            if( physk > box ) cycle
                            cmat(physh,physk) = cmat(physh,physk) + pfts_odd(irot,k,icls)
                            norm(physh,physk) = norm(physh,physk) + 1.d0
                        end do
                    end do
                case('merged')
                    do irot = 1,pftsz
                        do k = kfromto(1),kfromto(2)
                            phys  = pftcc%get_coord(irot,k)
                            physh = nint(phys(1)) + 1
                            physk = nint(phys(2)) + c
                            if( physk > box ) cycle
                            cmat(physh,physk) = cmat(physh,physk) + pfts_refs(irot,k,icls)
                            norm(physh,physk) = norm(physh,physk) + 1
                        end do
                    end do
            end select
            where( norm > DTINY )
                cmat = cmat / norm
            elsewhere
                cmat = 0.d0
            end where
            ! irot = self%pftsz+1, eg. angle=180.
            do k = 1,box/2-1
                cmat(1,k+c) = conjg(cmat(1,c-k))
            enddo
            ! arbitrary magnitude
            cmat(1,c) = DCMPLX_ZERO
            call cavgs(icls)%new([box,box,1],smpd)
            call cavgs(icls)%set_cmat(cmplx(cmat,kind=sp))
            call cavgs(icls)%shift_phorig()
            call cavgs(icls)%ifft
            ! call cavgs(icls)%div_w_instrfun('nn')
        enddo
    end subroutine polar_cavger_refs2cartesian

    ! I/O

    subroutine polar_cavger_write( fname, which )
        character(len=*),  intent(in) :: fname, which
        character(len=:), allocatable :: fname_here
        fname_here  = trim(fname)
        select case(which)
            case('even')
                call write_pft_array(pfts_even, fname_here)
            case('odd')
                call write_pft_array(pfts_odd, fname_here)
            case('merged')
                call write_pft_array(pfts_refs, fname_here)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_write

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
                call read_pft_array(fname_here, pfts_refs)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine polar_cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        character(len=:), allocatable :: cae, cao, cte, cto
        allocate(cae, source='cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT)
        allocate(cao, source='cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT)
        allocate(cte, source='ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT)
        allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT)
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

    subroutine polar_cavger_restore_classes( pinds )
        use simple_ctf,                 only: ctf
        use simple_builder,             only: build_glob
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, killimgbatch, prepimgbatch
        integer,                 intent(in)    :: pinds(:)
        type(ctfparams)          :: ctfparms
        type(polarft_corrcalc)   :: pftcc
        type(ctf)                :: tfun
        type(image), allocatable :: cavgs(:)
        real, allocatable :: incr_shifts(:,:)
        real    :: sdevnoise
        integer :: iptcl, ithr, icls, i, nptcls
        logical :: eo, l_ctf
        nptcls = size(pinds)
        params_glob%kfromto = [2, nint(real(params_glob%box)/2.2)]
        call pftcc%new(params_glob%ncls, [1,nptcls], params_glob%kfromto)
        l_ctf = build_glob%spproj%get_ctfflag('ptcl2D',iptcl=params_glob%fromp).ne.'no'
        call prepimgbatch(nptcls)
        call discrete_read_imgbatch(nptcls, pinds, [1,nptcls])
        call pftcc%reallocate_ptcls(nptcls, pinds)
        if( l_ctf ) call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl2D')
        call build_glob%img_crop_polarizer%init_polarizer(pftcc, params_glob%alpha)
        allocate(incr_shifts(2,nptcls),source=0.)
        !$omp parallel do default(shared) private(iptcl,i,ithr,eo,sdevnoise,ctfparms,tfun)&
        !$omp schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(i)
            ! normalization
            call build_glob%imgbatch(i)%norm_noise(build_glob%lmsk, sdevnoise)
            if( trim(params_glob%gridding).eq.'yes' )then
                call build_glob%img_crop_polarizer%div_by_instrfun(build_glob%imgbatch(i))
            endif
            call build_glob%imgbatch(i)%fft
            ! shift
            incr_shifts(:,i) = build_glob%spproj_field%get_2Dshift(iptcl)
            call build_glob%imgbatch(i)%shift2Dserial(-incr_shifts(:,i))
            ! phase-flipping
            ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
            select case(ctfparms%ctfflag)
                case(CTFFLAG_NO, CTFFLAG_FLIP)
                case(CTFFLAG_YES)
                    tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                    call tfun%apply_serial(build_glob%imgbatch(i), 'flip', ctfparms)
                case DEFAULT
                    THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' polar_cavger_restore_classes')
            end select
            ! even/odd
            eo = build_glob%spproj_field%get_eo(iptcl) < 0.5
            ! polar transform
            call build_glob%img_crop_polarizer%polarize(pftcc, build_glob%imgbatch(i), iptcl, .true., eo, mask=build_glob%l_resmsk)
        end do
        !$omp end parallel do
        call killimgbatch
        call polar_cavger_new(pftcc)
        call polar_cavger_update_sums(nptcls, pinds, build_glob%spproj, pftcc, incr_shifts)
        call polar_cavger_merge_eos_and_norm
        call polar_cavger_calc_and_write_frcs_and_eoavg(FRCS_FILE)
        call polar_cavger_write('cavgs_even.bin', 'even')
        call polar_cavger_write('cavgs_odd.bin',  'odd')
        call polar_cavger_write('cavgs.bin',      'merged')
        allocate(cavgs(params_glob%ncls))
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'even')
        do icls = 1,params_glob%ncls
            call cavgs(icls)%write('cavgs_even.mrc', icls)
        enddo
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'odd')
        do icls = 1,params_glob%ncls
            call cavgs(icls)%write('cavgs_odd.mrc', icls)
        enddo
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'odd')
        do icls = 1,params_glob%ncls
            call cavgs(icls)%write('cavgs_odd.mrc', icls)
        enddo
        call pftcc%kill
        call polar_cavger_kill
    end subroutine polar_cavger_restore_classes

    subroutine polar_cavger_kill
        if( allocated(pfts_even) ) deallocate(pfts_even,pfts_odd,ctf2_even,ctf2_odd,pfts_refs)
        smpd       = 0.
        ncls       = 0
        kfromto(2) = 0
        pftsz = 0
    end subroutine polar_cavger_kill

    ! PRIVATE UTILITIES

    subroutine calc_frc( pft1, pft2, n, frc )
        complex(dp), intent(in)    :: pft1(pftsz,kfromto(1):kfromto(2)), pft2(pftsz,kfromto(1):kfromto(2))
        integer,     intent(in)    :: n
        real(sp),    intent(inout) :: frc(1:n)
        real(dp) :: denom
        integer  :: k
        frc(1:kfromto(1)-1) = 0.999
        do k = kfromto(1), kfromto(2)
            denom = sum(csq_fast(pft1(:,k))) * sum(csq_fast(pft2(:,k)))
            if( denom > DTINY )then
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / sqrt(denom), sp)
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
    subroutine write_pft_array( array, fname )
        complex(dp),      intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        character(len=*), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) array
        call fclose(funit)
    end subroutine write_pft_array

    subroutine write_ctf2_array( array, fname )
        real(dp),         intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        character(len=*), intent(in) :: fname
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) array
        call fclose(funit)
    end subroutine write_ctf2_array

    subroutine read_pft_array( fname, array )
        character(len=*),         intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        complex(dp), allocatable :: tmp(:,:,:)
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
            read(unit=funit, pos=(sizeof(dims)+1)) array
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
                if( (k >= dims(2)) .or. (k <= dims(3)) )then
                    array(:,k,:) = tmp(:,k,:)   ! from stored array
                else
                    array(:,k,:) = 0.d0         ! pad with zeros
                endif
            enddo
            deallocate(tmp)
        endif
        call fclose(funit)
    end subroutine read_pft_array

    subroutine read_ctf2_array( fname, array )
        character(len=*),      intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        real(dp), allocatable :: tmp(:,:,:)
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
            read(unit=funit, pos=(sizeof(dims)+1)) array
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
                if( (k >= dims(2)) .or. (k <= dims(3)) )then
                    array(:,k,:) = tmp(:,k,:)   ! from stored array
                else
                    array(:,k,:) = 0.d0         ! pad with zeros
                endif
            enddo
            deallocate(tmp)
        endif
        call fclose(funit)
    end subroutine read_ctf2_array

    ! TEST UNIT

    subroutine test_polarops
        use simple_cmdline,    only: cmdline
        use simple_parameters, only: parameters
        use simple_builder,    only: builder
        integer,     parameter :: N=128
        integer,     parameter :: NIMGS=100
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
            ! shift = 10.*[ran3(), ran3()] - 5.
            shift = 0.
            ang   = 360. * ran3()
            ! ang   = 0.
            eo    = 0
            if( .not.is_even(i) ) eo = 1
            icls  = ceiling(ran3()*4.)

            call img%copy_fast(tmpl_img)
            call img%fft
            call img%shift2Dserial(-shift)
            call img%ifft
            call img%rtsq(ang, 0.,0.)
            call img%write('rotimgs.mrc', i)
            call img%fft

            call b%spproj_field%set_euler(i, [0.,0.,ang])
            call b%spproj_field%set_shift(i, [0.,0.])
            call b%spproj_field%set(i,'w',1.0)
            call b%spproj_field%set(i,'state',1)
            call b%spproj_field%set(i,'class', icls)
            call b%spproj_field%set(i,'eo',eo)
            shifts(:,i) = shift
            call b%img_crop_polarizer%polarize(pftcc, img, i, isptcl=.true., iseven=eo==0, mask=b%l_resmsk)
        enddo
        !!
        shifts = 0.
        !!
        call polar_cavger_new(pftcc)
        call polar_cavger_update_sums(NIMGS, pinds, b%spproj, pftcc, shifts)
        call polar_cavger_merge_eos_and_norm
        call polar_cavger_calc_and_write_frcs_and_eoavg(FRCS_FILE)
        call polar_cavger_write('cavgs_even.bin', 'even')
        call polar_cavger_write('cavgs_odd.bin', 'odd')
        call polar_cavger_write('cavgs.bin', 'merged')
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'even')
        do icls = 1,NCLS
            call cavgs(icls)%write('cavgs_even.mrc', icls)
        enddo
        call polar_cavger_refs2cartesian(pftcc, cavgs, 'odd')
        do icls = 1,NCLS
            call cavgs(icls)%write('cavgs_odd.mrc', icls)
        enddo
    end subroutine test_polarops

end module simple_polarops