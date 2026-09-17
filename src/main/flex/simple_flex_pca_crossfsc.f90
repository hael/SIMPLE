!@descr: versioned cross-fit-FSC artifact (flex_pca_crossfsc.bin): writer, reader, SSNR conversion, series restart
!! Spec: doc/for_developers/ideas/flex_pca_crossfsc_shrinkage_marching_spec.md (drop-in spec, 2026-08-31).
!! The paired flex_pca engine's lasting output is the per-component, per-shell cross-fit FSC between
!! two independent half-fits. It is persisted every iteration as a first-class versioned artifact.
!! Of the spec's three consumers only (1) SSNR/tau^2 shrinkage (spec par.3-4) remains; frequency
!! marching (par.5), the stopping series (par.6) and the single-fit SCAFFOLDING writer were removed
!! 2026-09-07 (never adopted by a recipe). The record layout is unchanged: `paired`, `khi_shared`,
!! `march_on` and `s_stop` stay as file-format fields (march_on is always 0 now), and
!! crossfsc_assert_paired still THROWs when a paired=0 record from an old artifact is asked to drive
!! rank selection.
module simple_flex_pca_crossfsc
use simple_core_module_api
use simple_flex_reconstructor_latent_ops, only: pair_index
implicit none
private
#include "simple_local_flags.inc"

public :: crossfsc_record, crossfsc_file
public :: crossfsc_load, crossfsc_write, crossfsc_append, crossfsc_latest_upto
public :: crossfsc_assert_paired, crossfsc_kill, crossfsc_kill_record
public :: crossfsc_to_invtau2, crossfsc_harvest_h, crossfsc_stop_stat
public :: crossfsc_inband_mean, crossfsc_khi_deepest
public :: COV_XFSC_FNAME

character(len=8), parameter :: COV_XFSC_MAGIC   = 'SIMPLFXF'
integer,          parameter :: COV_XFSC_VERSION = 1
character(len=*), parameter :: COV_XFSC_FNAME        = 'flex_pca_crossfsc.bin'
character(len=*), parameter :: COV_XFSC_SERIES_FNAME = 'flex_pca_crossfsc_series.txt'
!> per-shell sampling below this is a dead shell, matching add_invtausq2rho's rsum floor
!! (simple_reconstructor.f90:1128)
real(dp),         parameter :: COV_XFSC_H_FLOOR = 1.0d-10

!> One per-iteration record (spec par.2.4, version 1). Unmatched components appear in the per-fit
!! blocks (fsc_int/h/eigvals) but not in the matched blocks (match_*/fsc_cross).
type crossfsc_record
    integer :: it_eff     = 0   !< global iteration stamp (never a worker-local counter)
    integer :: ncomp_a    = 0, ncomp_b = 0   !< per-fit delivered ranks
    integer :: kmatch     = 0   !< number of matched pairs
    integer :: khi_a      = 0, khi_b = 0     !< each fit's internal-FSC band (par.5.2 criterion)
    integer :: khi_shared = 0   !< the shared band in force (par.5)
    integer :: reg_mode   = 0   !< the par.4 arm ACTIVE this iteration (0 when degraded)
    integer :: march_on   = 0   !< 1 when the par.5 shared band was APPLIED this iteration
    integer,  allocatable :: match_a(:), match_b(:)  !< signed-permutation pairing (component indices)
    integer,  allocatable :: match_sign(:)           !< +1/-1
    real,     allocatable :: match_cos(:)            !< masked matched cosine (par.5 criterion scalar)
    real,     allocatable :: fsc_cross(:,:)          !< (filtsz,kmatch) cross-fit FSC curves
    real,     allocatable :: fsc_int_a(:,:)          !< (filtsz,ncomp_a) fit A internal e/o FSC
    real,     allocatable :: fsc_int_b(:,:)          !< (filtsz,ncomp_b)
    real,     allocatable :: h_a(:,:)                !< (filtsz,ncomp_a) fit A sampling profiles (par.2.3)
    real,     allocatable :: h_b(:,:)                !< (filtsz,ncomp_b)
    integer,  allocatable :: cnt(:)                  !< (filtsz) shared per-shell voxel counts
    real(dp), allocatable :: eigvals_a(:), eigvals_b(:)  !< latent variances (amplitude diagnostic)
    real(dp) :: s_stop = 0.d0   !< stopping statistic S(t) at khi_cmp (par.6), precomputed
end type crossfsc_record

!> The artifact: file header + all records so far. Full rewrite per iteration (status='replace',
!! stream unformatted -- the write_embedding_cache idiom), so the artifact is restart-complete:
!! on master restart the series reloads exactly.
type crossfsc_file
    integer :: paired     = 0   !< 1 = honest paired fits; 0 = scaffolding (internal e/o)
    integer :: pairing_id = 0   !< balanced mod-4 pairing id (1..3); 0 = no pairing (scaffolding)
    integer :: box_crop   = 0
    integer :: filtsz     = 0   !< fdim(box_crop)-1
    real    :: smpd_crop  = 0.
    integer :: khi_full   = 0   !< full band cap (covariance_kfromto at the production lp)
    integer :: khi_cmp    = 0   !< FIXED comparison band (par.6), frozen at first record
    integer :: nrec       = 0
    type(crossfsc_record), allocatable :: recs(:)
end type crossfsc_file

contains

    ! ============ artifact I/O ============

    !> Load flex_pca_crossfsc.bin from the working directory. found=.false. when absent (fresh run).
    !! Hard-fails on magic and version mismatch: any layout change bumps COV_XFSC_VERSION.
    subroutine crossfsc_load( self, found )
        type(crossfsc_file), intent(inout) :: self
        logical,             intent(out)   :: found
        character(len=len(COV_XFSC_MAGIC)) :: magic
        integer :: funit, io_stat, ver, irec
        call crossfsc_kill(self)
        found = .false.
        if( .not. file_exists(string(COV_XFSC_FNAME)) ) return
        call fopen(funit, file=string(COV_XFSC_FNAME), access='STREAM', action='READ', &
            &status='OLD', iostat=io_stat)
        call fileiochk('crossfsc_load; open '//COV_XFSC_FNAME, io_stat)
        read(funit, iostat=io_stat) magic, ver
        call fileiochk('crossfsc_load; magic', io_stat)
        if( magic /= COV_XFSC_MAGIC ) THROW_HARD('not a flex_pca crossfsc artifact: '//COV_XFSC_FNAME)
        if( ver /= COV_XFSC_VERSION )then
            write(logfhandle,'(A,I0,A,I0)') 'crossfsc artifact version found: ',ver,'  expected: ', &
                &COV_XFSC_VERSION
            THROW_HARD('flex_pca crossfsc artifact version mismatch; any layout change bumps the &
                &version -- delete '//COV_XFSC_FNAME//' to restart the series')
        endif
        read(funit, iostat=io_stat) self%paired, self%pairing_id, self%box_crop, self%filtsz, &
            &self%smpd_crop, self%khi_full, self%khi_cmp, self%nrec
        call fileiochk('crossfsc_load; header', io_stat)
        if( self%filtsz < 1 .or. self%nrec < 0 ) THROW_HARD('corrupt crossfsc header: '//COV_XFSC_FNAME)
        allocate(self%recs(max(1,self%nrec)))
        do irec = 1, self%nrec
            call read_record(funit, self%filtsz, self%recs(irec))
        end do
        call fclose(funit)
        found = .true.
    end subroutine crossfsc_load

    subroutine read_record( funit, filtsz, rec )
        integer,               intent(in)    :: funit, filtsz
        type(crossfsc_record), intent(inout) :: rec
        integer :: io_stat
        call crossfsc_kill_record(rec)
        read(funit, iostat=io_stat) rec%it_eff, rec%ncomp_a, rec%ncomp_b, rec%kmatch, &
            &rec%khi_a, rec%khi_b, rec%khi_shared, rec%reg_mode, rec%march_on
        call fileiochk('crossfsc read_record; scalars', io_stat)
        if( rec%kmatch < 0 .or. rec%ncomp_a < 0 .or. rec%ncomp_b < 0 ) &
            &THROW_HARD('corrupt crossfsc record')
        allocate(rec%match_a(rec%kmatch), rec%match_b(rec%kmatch), rec%match_sign(rec%kmatch))
        allocate(rec%match_cos(rec%kmatch), rec%fsc_cross(filtsz,rec%kmatch))
        allocate(rec%fsc_int_a(filtsz,rec%ncomp_a), rec%fsc_int_b(filtsz,rec%ncomp_b))
        allocate(rec%h_a(filtsz,rec%ncomp_a), rec%h_b(filtsz,rec%ncomp_b))
        allocate(rec%cnt(filtsz), rec%eigvals_a(rec%ncomp_a), rec%eigvals_b(rec%ncomp_b))
        if( rec%kmatch > 0 )then
            read(funit, iostat=io_stat) rec%match_a, rec%match_b, rec%match_sign, rec%match_cos, &
                &rec%fsc_cross
            call fileiochk('crossfsc read_record; matched blocks', io_stat)
        endif
        if( rec%ncomp_a > 0 )then
            read(funit, iostat=io_stat) rec%fsc_int_a, rec%h_a
            call fileiochk('crossfsc read_record; fit A blocks', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            read(funit, iostat=io_stat) rec%fsc_int_b, rec%h_b
            call fileiochk('crossfsc read_record; fit B blocks', io_stat)
        endif
        read(funit, iostat=io_stat) rec%cnt
        call fileiochk('crossfsc read_record; cnt', io_stat)
        if( rec%ncomp_a > 0 )then
            read(funit, iostat=io_stat) rec%eigvals_a
            call fileiochk('crossfsc read_record; eigvals_a', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            read(funit, iostat=io_stat) rec%eigvals_b
            call fileiochk('crossfsc read_record; eigvals_b', io_stat)
        endif
        read(funit, iostat=io_stat) rec%s_stop
        call fileiochk('crossfsc read_record; s_stop', io_stat)
    end subroutine read_record

    !> Full rewrite of the artifact (all records so far) + the human-greppable series mirror.
    !! Write-to-tmp-and-rename, so a reader never sees a torn file (the part-file idiom).
    subroutine crossfsc_write( self )
        type(crossfsc_file), intent(in) :: self
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat, irec
        fname     = string(COV_XFSC_FNAME)
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', &
            &iostat=io_stat)
        call fileiochk('crossfsc_write; open', io_stat)
        write(funit, iostat=io_stat) COV_XFSC_MAGIC, COV_XFSC_VERSION
        call fileiochk('crossfsc_write; magic', io_stat)
        write(funit, iostat=io_stat) self%paired, self%pairing_id, self%box_crop, self%filtsz, &
            &self%smpd_crop, self%khi_full, self%khi_cmp, self%nrec
        call fileiochk('crossfsc_write; header', io_stat)
        do irec = 1, self%nrec
            call write_record(funit, self%recs(irec))
        end do
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
        call write_series_mirror(self)
    end subroutine crossfsc_write

    subroutine write_record( funit, rec )
        integer,               intent(in) :: funit
        type(crossfsc_record), intent(in) :: rec
        integer :: io_stat
        write(funit, iostat=io_stat) rec%it_eff, rec%ncomp_a, rec%ncomp_b, rec%kmatch, &
            &rec%khi_a, rec%khi_b, rec%khi_shared, rec%reg_mode, rec%march_on
        call fileiochk('crossfsc write_record; scalars', io_stat)
        if( rec%kmatch > 0 )then
            write(funit, iostat=io_stat) rec%match_a, rec%match_b, rec%match_sign, rec%match_cos, &
                &rec%fsc_cross
            call fileiochk('crossfsc write_record; matched blocks', io_stat)
        endif
        if( rec%ncomp_a > 0 )then
            write(funit, iostat=io_stat) rec%fsc_int_a, rec%h_a
            call fileiochk('crossfsc write_record; fit A blocks', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            write(funit, iostat=io_stat) rec%fsc_int_b, rec%h_b
            call fileiochk('crossfsc write_record; fit B blocks', io_stat)
        endif
        write(funit, iostat=io_stat) rec%cnt
        call fileiochk('crossfsc write_record; cnt', io_stat)
        if( rec%ncomp_a > 0 )then
            write(funit, iostat=io_stat) rec%eigvals_a
            call fileiochk('crossfsc write_record; eigvals_a', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            write(funit, iostat=io_stat) rec%eigvals_b
            call fileiochk('crossfsc write_record; eigvals_b', io_stat)
        endif
        write(funit, iostat=io_stat) rec%s_stop
        call fileiochk('crossfsc write_record; s_stop', io_stat)
    end subroutine write_record

    !> Human-greppable mirror of the stopping series only; the .bin is authoritative.
    subroutine write_series_mirror( self )
        type(crossfsc_file), intent(in) :: self
        integer :: u, irec, k
        call del_file(COV_XFSC_SERIES_FNAME)
        open(newunit=u, file=COV_XFSC_SERIES_FNAME, status='replace', action='write')
        write(u,'(A,I0,A,I0,A,I0,A,I0)') '# flex_pca crossfsc series  paired=',self%paired, &
            &'  pairing_id=',self%pairing_id,'  khi_cmp=',self%khi_cmp,'  khi_full=',self%khi_full
        write(u,'(A)') '# it  S(t)  khi_a  khi_b  khi_shared  reg_mode  march_on  inband_mean_per_matched_component'
        do irec = 1, self%nrec
            write(u,'(I6,1X,ES14.6,5(1X,I5))', advance='no') self%recs(irec)%it_eff, &
                &self%recs(irec)%s_stop, self%recs(irec)%khi_a, self%recs(irec)%khi_b, &
                &self%recs(irec)%khi_shared, self%recs(irec)%reg_mode, self%recs(irec)%march_on
            do k = 1, self%recs(irec)%kmatch
                write(u,'(1X,F8.4)', advance='no') &
                    &crossfsc_inband_mean(self%recs(irec)%fsc_cross(:,k), self%khi_cmp)
            end do
            write(u,*)
        end do
        close(u)
    end subroutine write_series_mirror

    !> Append one record, keeping the series stamp-monotonic: any trailing records with
    !! it_eff >= the new stamp are dropped first (a re-run over an existing artifact replaces
    !! from its own first iteration on, which is what makes the full-rewrite restart-complete).
    subroutine crossfsc_append( self, rec )
        type(crossfsc_file),   intent(inout) :: self
        type(crossfsc_record), intent(in)    :: rec
        type(crossfsc_record), allocatable   :: tmp(:)
        integer :: nkeep, irec
        nkeep = 0
        do irec = 1, self%nrec
            if( self%recs(irec)%it_eff < rec%it_eff )then
                nkeep = irec
            else
                exit
            endif
        end do
        allocate(tmp(nkeep+1))
        do irec = 1, nkeep
            tmp(irec) = self%recs(irec)
        end do
        tmp(nkeep+1) = rec
        if( allocated(self%recs) )then
            do irec = 1, size(self%recs)
                call crossfsc_kill_record(self%recs(irec))
            end do
            deallocate(self%recs)
        endif
        call move_alloc(tmp, self%recs)
        self%nrec = nkeep + 1
    end subroutine crossfsc_append

    !> Index of the latest record with it_eff <= it (the par.4.2 timing rule: iteration t may
    !! consume only records stamped <= t-1, so callers pass it = t-1). 0 when none qualifies.
    integer function crossfsc_latest_upto( self, it ) result( irec )
        type(crossfsc_file), intent(in) :: self
        integer,             intent(in) :: it
        integer :: i
        irec = 0
        do i = 1, self%nrec
            if( self%recs(i)%it_eff <= it ) irec = i
        end do
    end function crossfsc_latest_upto

    !> The scaffolding guard: a paired=0 record is not a cross-fit statistic. Rank selection and
    !! stopping THROW here; the in-loop shrinkage arms and marching degrade silently instead
    !! (their own call sites log the degradation).
    subroutine crossfsc_assert_paired( self, purpose )
        type(crossfsc_file), intent(in) :: self
        character(len=*),    intent(in) :: purpose
        if( self%paired /= 1 )then
            write(logfhandle,'(A,A)') 'crossfsc consumer refused for purpose: ', trim(purpose)
            THROW_HARD('flex_pca crossfsc artifact holds SCAFFOLDING (paired=0) records, which are &
                &internal even/odd curves, not cross-fit statistics; refusing to use them for the &
                &purpose named above')
        endif
    end subroutine crossfsc_assert_paired

    subroutine crossfsc_kill_record( rec )
        type(crossfsc_record), intent(inout) :: rec
        if( allocated(rec%match_a)    ) deallocate(rec%match_a, rec%match_b, rec%match_sign)
        if( allocated(rec%match_cos)  ) deallocate(rec%match_cos)
        if( allocated(rec%fsc_cross)  ) deallocate(rec%fsc_cross)
        if( allocated(rec%fsc_int_a)  ) deallocate(rec%fsc_int_a)
        if( allocated(rec%fsc_int_b)  ) deallocate(rec%fsc_int_b)
        if( allocated(rec%h_a)        ) deallocate(rec%h_a)
        if( allocated(rec%h_b)        ) deallocate(rec%h_b)
        if( allocated(rec%cnt)        ) deallocate(rec%cnt)
        if( allocated(rec%eigvals_a)  ) deallocate(rec%eigvals_a)
        if( allocated(rec%eigvals_b)  ) deallocate(rec%eigvals_b)
        rec%it_eff = 0; rec%ncomp_a = 0; rec%ncomp_b = 0; rec%kmatch = 0
        rec%khi_a = 0; rec%khi_b = 0; rec%khi_shared = 0
        rec%reg_mode = 0; rec%march_on = 0; rec%s_stop = 0.d0
    end subroutine crossfsc_kill_record

    subroutine crossfsc_kill( self )
        type(crossfsc_file), intent(inout) :: self
        integer :: irec
        if( allocated(self%recs) )then
            do irec = 1, size(self%recs)
                call crossfsc_kill_record(self%recs(irec))
            end do
            deallocate(self%recs)
        endif
        self%paired = 0; self%pairing_id = 0; self%box_crop = 0; self%filtsz = 0
        self%smpd_crop = 0.; self%khi_full = 0; self%khi_cmp = 0; self%nrec = 0
    end subroutine crossfsc_kill

    ! ============ the S.11/S.13 conversion (spec par.3.1) ============

    !> Convert one component's cross-fit FSC curve F plus that fit's own per-shell sampling
    !! profile H into the per-shell inverse prior variance invtau2, mirroring
    !! simple_reconstructor::add_invtausq2rho step for step -- H IS the rsum/cnt pass
    !! precomputed, which is the entire reason it lives in the artifact.
    !! fudge = params%tau (the ml_reg knob); k_lo = max(6, covariance_kfromto(1)), the
    !! reslim_ind analog: no addition at very low resolution (signal assumed infinite).
    subroutine crossfsc_to_invtau2( fsc, h, fudge, k_lo, invtau2 )
        real,    intent(in)  :: fsc(:)      !< one component's cross-fit FSC (1:filtsz)
        real,    intent(in)  :: h(:)        !< that fit's own sampling profile (1:filtsz)
        real,    intent(in)  :: fudge       !< params%tau
        integer, intent(in)  :: k_lo        !< low-resolution exemption index
        real,    intent(out) :: invtau2(:)  !< (1:filtsz)
        real    :: fc, ssnr, sig2, tau2
        integer :: sh, n
        n = min(size(fsc), min(size(h), size(invtau2)))
        invtau2 = 0.0
        do sh = 1, n
            if( sh < k_lo ) cycle                          ! no addition below k_lo (:1137,:1144)
            fc   = min(0.999, max(0.001, fsc(sh)))
            ssnr = fc / (1.0 - fc)                         ! (:1108-1112)
            if( real(h(sh),dp) > COV_XFSC_H_FLOOR )then    ! (:1128-1132)
                sig2 = 1.0 / h(sh)
            else
                sig2 = 0.0                                 ! unpopulated shell
            endif
            tau2 = ssnr * sig2                             ! (:1134)
            if( tau2 > TINY .and. fsc(sh) > 0.0 )then
                invtau2(sh) = 1.0 / (max(fudge, TINY)*tau2)          ! (:1147)
            else
                ! kill branch, per-shell here rather than per-voxel: the shell-mean H stands in
                ! for the per-voxel rho of (:1149), within a factor wherever the shell is populated
                invtau2(sh) = min(1.0e3, 1.0e3 * h(sh))
            endif
        end do
    end subroutine crossfsc_to_invtau2

    ! ============ the sampling profile H (spec par.2.3) ============

    !> Per-component, per-shell mean of the DIAGONAL rows of one packed coupled density --
    !! the exact analog of the rsum/cnt pass in add_invtausq2rho, run over pair_index(q,q).
    !! `lb` are the frequency-space lower bounds of the exp lattice (lbound(cmat_exp)) and
    !! `nyq` its spherical Nyquist (get_lfny(1)): the shell convention and support rule of
    !! solve_coupled_basis_exp itself. `cnt` is shared across components (stored once per record).
    !! Harvest BEFORE any invtau2 is added: H is the fit's own ACCUMULATED sampling.
    subroutine crossfsc_harvest_h( rho, npairs, ncomp, lb, nyq, filtsz, h, cnt )
        real,    intent(in)  :: rho(:,:,:,:)     !< packed coupled density (npairs, nh, nk, nm)
        integer, intent(in)  :: npairs, ncomp, lb(3), nyq, filtsz
        real,    intent(out) :: h(filtsz,ncomp)
        integer, intent(out) :: cnt(filtsz)
        real(dp) :: rsum(filtsz,ncomp)
        integer  :: nh, nk, nm, ih, ik, im, hf, kf, mf, sh, q, shmax
        if( size(rho,1) /= npairs ) THROW_HARD('crossfsc_harvest_h: rho leading extent /= npairs')
        nh = size(rho,2); nk = size(rho,3); nm = size(rho,4)
        rsum  = 0.d0
        cnt   = 0
        shmax = min(filtsz, nyq)
        !$omp parallel do collapse(2) default(shared) schedule(static) &
        !$omp private(im,ik,ih,hf,kf,mf,sh,q) reduction(+:cnt,rsum) proc_bind(close)
        do im = 1, nm
            do ik = 1, nk
                do ih = 1, nh
                    hf = lb(1) + ih - 1
                    kf = lb(2) + ik - 1
                    mf = lb(3) + im - 1
                    sh = nint(sqrt(real(hf*hf + kf*kf + mf*mf)))
                    if( sh < 1 .or. sh > shmax ) cycle
                    cnt(sh) = cnt(sh) + 1
                    do q = 1, ncomp
                        rsum(sh,q) = rsum(sh,q) + real(rho(pair_index(q,q),ih,ik,im),dp)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        do q = 1, ncomp
            do sh = 1, filtsz
                if( cnt(sh) > 0 )then
                    h(sh,q) = real(rsum(sh,q) / real(cnt(sh),dp))
                else
                    h(sh,q) = 0.0
                endif
            end do
        end do
    end subroutine crossfsc_harvest_h

    ! ============ aggregations (spec par.5.2 / par.6.1) ============

    !> Mean FSC over shells 2..khi (the fmean_dg aggregation shape applied from shell 2,
    !! matching the stopping statistic's support).
    real function crossfsc_inband_mean( fsc, khi ) result( fmean )
        real,    intent(in) :: fsc(:)
        integer, intent(in) :: khi
        integer :: k_hi
        k_hi  = max(2, min(khi, size(fsc)))
        fmean = sum(fsc(2:k_hi)) / real(k_hi - 1)
    end function crossfsc_inband_mean

    !> The stopping statistic S(t) at the FIXED comparison band khi_cmp (spec par.6.1):
    !! mean over matched pairs of the in-band mean cross-fit FSC. Computed over shells up to
    !! khi_cmp even while the working band is below it (those shells score ~0, honestly).
    real(dp) function crossfsc_stop_stat( rec, khi_cmp ) result( s )
        type(crossfsc_record), intent(in) :: rec
        integer,               intent(in) :: khi_cmp
        integer :: k
        s = 0.d0
        if( rec%kmatch < 1 ) return
        do k = 1, rec%kmatch
            s = s + real(crossfsc_inband_mean(rec%fsc_cross(:,k), khi_cmp), dp)
        end do
        s = s / real(rec%kmatch, dp)
    end function crossfsc_stop_stat

    !> Criterion crossing of one FSC curve -- the get_find_at_crit analog
    !! (simple_math_ft.f90:188-210): first shell h >= 3 with fsc(h) < crit returns h-1;
    !! a curve that never crosses returns size-1.
    integer function crossfsc_find_at_crit( fsc, crit ) result( find )
        real, intent(in) :: fsc(:)
        real, intent(in) :: crit
        integer :: n, h
        n    = size(fsc)
        find = n - 1
        do h = 3, n - 1
            if( fsc(h) >= crit )then
                cycle
            else
                find = h - 1
                exit
            endif
        end do
        find = max(1, min(find, n - 1))
    end function crossfsc_find_at_crit

    !> The band driver over a set of curves: the DEEPEST-crossing component (spec par.5.2 item 2,
    !! the "best resolved state drives" analog). Not the mean: one honest component earning band
    !! is the point of marching.
    integer function crossfsc_khi_deepest( curves, ncurves, crit ) result( khi )
        real,    intent(in) :: curves(:,:)   !< (filtsz, ncurves)
        integer, intent(in) :: ncurves
        real,    intent(in) :: crit
        integer :: k
        khi = 1
        do k = 1, ncurves
            khi = max(khi, crossfsc_find_at_crit(curves(:,k), crit))
        end do
    end function crossfsc_khi_deepest

end module simple_flex_pca_crossfsc
