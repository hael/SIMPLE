!@descr: versioned cross-fit-FSC artifact (flex_pca_crossfsc.bin): writer, reader, SSNR conversion, series restart
!! The paired flex_pca engine persists the per-component, per-shell cross-fit FSC between its two
!! half-fits every iteration. Its only consumer is the SSNR/tau^2 ridge (crossfsc_to_invtau2), so a
!! record holds what that ridge reads: the pairing, the cross-fit FSC curves and each fit's sampling H.
module simple_flex_pca_crossfsc
use simple_core_module_api, only: dp, fclose, file_exists, fileiochk, fopen, logfhandle, simple_exception, &
    &simple_rename, string, tiny
use simple_flex_reconstructor_latent_ops, only: pair_index
implicit none
private
#include "simple_local_flags.inc"

public :: crossfsc_record, crossfsc_file
public :: crossfsc_to_invtau2, crossfsc_harvest_h, crossfsc_inband_mean
public :: COV_XFSC_FNAME

character(len=8), parameter :: COV_XFSC_MAGIC        = 'SIMPLFXF'
integer,          parameter :: COV_XFSC_VERSION      = 2   ! 2: only what the ridge reads
character(len=*), parameter :: COV_XFSC_FNAME        = 'flex_pca_crossfsc.bin'
!> per-shell sampling below this is a dead shell, matching add_invtausq2rho's rsum floor
real(dp),         parameter :: COV_XFSC_H_FLOOR      = 1.0d-10

!> One per-iteration record (layout COV_XFSC_VERSION). Unmatched components appear in the per-fit
!! sampling blocks (h_a/h_b) but not in the matched blocks (match_*/fsc_cross).
type crossfsc_record
    integer :: it_eff  = 0                   !< global iteration stamp (never a worker-local counter)
    integer :: ncomp_a = 0, ncomp_b = 0      !< per-fit delivered ranks
    integer :: kmatch  = 0                   !< number of matched pairs
    integer, allocatable :: match_a(:), match_b(:)  !< pairing (component indices)
    real,    allocatable :: fsc_cross(:,:)          !< (filtsz,kmatch) cross-fit FSC curves (sign resolved)
    real,    allocatable :: h_a(:,:)                !< (filtsz,ncomp_a) fit A sampling profiles
    real,    allocatable :: h_b(:,:)                !< (filtsz,ncomp_b)
  contains
    procedure :: kill => crossfsc_kill_record
end type crossfsc_record

!> The artifact: file header + all records so far. Full rewrite per iteration (status='replace',
!! stream unformatted -- the write_embedding_cache idiom), so the artifact is restart-complete:
!! on master restart the series reloads exactly.
type crossfsc_file
    integer :: box_crop = 0
    integer :: filtsz   = 0   !< fdim(box_crop)-1
    integer :: nrec     = 0
    type(crossfsc_record), allocatable :: recs(:)
  contains
    procedure :: load        => crossfsc_load
    procedure :: write       => crossfsc_write
    procedure :: append      => crossfsc_append
    procedure :: latest_upto => crossfsc_latest_upto
    procedure :: kill        => crossfsc_kill
end type crossfsc_file

contains

    ! ============ artifact I/O ============

    !> Load flex_pca_crossfsc.bin from the working directory. found=.false. when absent (fresh run).
    !! Hard-fails on magic and version mismatch: any layout change bumps COV_XFSC_VERSION.
    subroutine crossfsc_load( self, found )
        class(crossfsc_file), intent(inout) :: self
        logical,             intent(out)   :: found
        character(len=len(COV_XFSC_MAGIC)) :: magic
        integer :: funit, io_stat, ver, irec
        call self%kill
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
        read(funit, iostat=io_stat) self%box_crop, self%filtsz, self%nrec
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
        call rec%kill
        read(funit, iostat=io_stat) rec%it_eff, rec%ncomp_a, rec%ncomp_b, rec%kmatch
        call fileiochk('crossfsc read_record; scalars', io_stat)
        if( rec%kmatch < 0 .or. rec%ncomp_a < 0 .or. rec%ncomp_b < 0 ) &
            &THROW_HARD('corrupt crossfsc record')
        allocate(rec%match_a(rec%kmatch), rec%match_b(rec%kmatch), rec%fsc_cross(filtsz,rec%kmatch))
        allocate(rec%h_a(filtsz,rec%ncomp_a), rec%h_b(filtsz,rec%ncomp_b))
        if( rec%kmatch > 0 )then
            read(funit, iostat=io_stat) rec%match_a, rec%match_b, rec%fsc_cross
            call fileiochk('crossfsc read_record; matched blocks', io_stat)
        endif
        if( rec%ncomp_a > 0 )then
            read(funit, iostat=io_stat) rec%h_a
            call fileiochk('crossfsc read_record; fit A sampling', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            read(funit, iostat=io_stat) rec%h_b
            call fileiochk('crossfsc read_record; fit B sampling', io_stat)
        endif
    end subroutine read_record

    !> Full rewrite of the artifact (all records so far). Write-to-tmp-and-rename, so a reader never
    !! sees a torn file (the part-file idiom).
    subroutine crossfsc_write( self )
        class(crossfsc_file), intent(in) :: self
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat, irec
        fname     = string(COV_XFSC_FNAME)
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', &
            &iostat=io_stat)
        call fileiochk('crossfsc_write; open', io_stat)
        write(funit, iostat=io_stat) COV_XFSC_MAGIC, COV_XFSC_VERSION
        call fileiochk('crossfsc_write; magic', io_stat)
        write(funit, iostat=io_stat) self%box_crop, self%filtsz, self%nrec
        call fileiochk('crossfsc_write; header', io_stat)
        do irec = 1, self%nrec
            call write_record(funit, self%recs(irec))
        end do
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine crossfsc_write

    subroutine write_record( funit, rec )
        integer,               intent(in) :: funit
        type(crossfsc_record), intent(in) :: rec
        integer :: io_stat
        write(funit, iostat=io_stat) rec%it_eff, rec%ncomp_a, rec%ncomp_b, rec%kmatch
        call fileiochk('crossfsc write_record; scalars', io_stat)
        if( rec%kmatch > 0 )then
            write(funit, iostat=io_stat) rec%match_a, rec%match_b, rec%fsc_cross
            call fileiochk('crossfsc write_record; matched blocks', io_stat)
        endif
        if( rec%ncomp_a > 0 )then
            write(funit, iostat=io_stat) rec%h_a
            call fileiochk('crossfsc write_record; fit A sampling', io_stat)
        endif
        if( rec%ncomp_b > 0 )then
            write(funit, iostat=io_stat) rec%h_b
            call fileiochk('crossfsc write_record; fit B sampling', io_stat)
        endif
    end subroutine write_record

    !> Append one record, keeping the series stamp-monotonic: any trailing records with
    !! it_eff >= the new stamp are dropped first (a re-run over an existing artifact replaces
    !! from its own first iteration on, which is what makes the full-rewrite restart-complete).
    subroutine crossfsc_append( self, rec )
        class(crossfsc_file),  intent(inout) :: self
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
                call self%recs(irec)%kill
            end do
            deallocate(self%recs)
        endif
        call move_alloc(tmp, self%recs)
        self%nrec = nkeep + 1
    end subroutine crossfsc_append

    !> Index of the latest record with it_eff <= it (the timing rule: iteration t may
    !! consume only records stamped <= t-1, so callers pass it = t-1). 0 when none qualifies.
    integer function crossfsc_latest_upto( self, it ) result( irec )
        class(crossfsc_file), intent(in) :: self
        integer,             intent(in) :: it
        integer :: i
        irec = 0
        do i = 1, self%nrec
            if( self%recs(i)%it_eff <= it ) irec = i
        end do
    end function crossfsc_latest_upto

    subroutine crossfsc_kill_record( rec )
        class(crossfsc_record), intent(inout) :: rec
        if( allocated(rec%match_a)   ) deallocate(rec%match_a, rec%match_b)
        if( allocated(rec%fsc_cross) ) deallocate(rec%fsc_cross)
        if( allocated(rec%h_a)       ) deallocate(rec%h_a)
        if( allocated(rec%h_b)       ) deallocate(rec%h_b)
        rec%it_eff = 0; rec%ncomp_a = 0; rec%ncomp_b = 0; rec%kmatch = 0
    end subroutine crossfsc_kill_record

    subroutine crossfsc_kill( self )
        class(crossfsc_file), intent(inout) :: self
        integer :: irec
        if( allocated(self%recs) )then
            do irec = 1, size(self%recs)
                call self%recs(irec)%kill
            end do
            deallocate(self%recs)
        endif
        self%box_crop = 0; self%filtsz = 0; self%nrec = 0
    end subroutine crossfsc_kill

    ! ============ the S.11/S.13 conversion ============

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
            if( sh < k_lo ) cycle                          ! no addition below k_lo
            fc   = min(0.999, max(0.001, fsc(sh)))
            ssnr = fc / (1.0 - fc)
            if( real(h(sh),dp) > COV_XFSC_H_FLOOR )then
                sig2 = 1.0 / h(sh)
            else
                sig2 = 0.0                                 ! unpopulated shell
            endif
            tau2 = ssnr * sig2
            if( tau2 > TINY .and. fsc(sh) > 0.0 )then
                invtau2(sh) = 1.0 / (max(fudge, TINY)*tau2)
            else
                ! kill branch, per-shell here rather than per-voxel: the shell-mean H stands in
                ! for add_invtausq2rho's per-voxel rho, within a factor wherever the shell is populated
                invtau2(sh) = min(1.0e3, 1.0e3 * h(sh))
            endif
        end do
    end subroutine crossfsc_to_invtau2

    ! ============ the sampling profile H ============

    !> Per-component, per-shell mean of the DIAGONAL rows of one packed coupled density --
    !! the exact analog of the rsum/cnt pass in add_invtausq2rho, run over pair_index(q,q).
    !! `lb` are the frequency-space lower bounds of the exp lattice (lbound(cmat_exp)) and
    !! `nyq` its spherical Nyquist (get_lfny(1)): the shell convention and support rule of
    !! solve_coupled_basis_exp itself. Harvest BEFORE any invtau2 is added: H is the fit's own
    !! ACCUMULATED sampling.
    subroutine crossfsc_harvest_h( rho, npairs, ncomp, lb, nyq, filtsz, h )
        real,    intent(in)  :: rho(:,:,:,:)     !< packed coupled density (npairs, nh, nk, nm)
        integer, intent(in)  :: npairs, ncomp, lb(3), nyq, filtsz
        real,    intent(out) :: h(filtsz,ncomp)
        real(dp) :: rsum(filtsz,ncomp)
        integer  :: cnt(filtsz)
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

    ! ============ aggregations ============

    !> Mean FSC over shells 2..khi (the fmean_dg aggregation shape applied from shell 2).
    real function crossfsc_inband_mean( fsc, khi ) result( fmean )
        real,    intent(in) :: fsc(:)
        integer, intent(in) :: khi
        integer :: k_hi
        k_hi  = max(2, min(khi, size(fsc)))
        fmean = sum(fsc(2:k_hi)) / real(k_hi - 1)
    end function crossfsc_inband_mean

end module simple_flex_pca_crossfsc
