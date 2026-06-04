!@descr: statistics/reporting helpers for 2D nonuniform filtering
module simple_nu_filter2D_stats
use simple_core_module_api
use simple_image, only: image
implicit none

public :: NU2D_LABEL_KIND, nu_filter2D_stats
public :: init_nu_filter2D_stats, kill_nu_filter2D_stats, merge_nu_filter2D_stats
public :: count_nu_filter2D_labels, accumulate_nu_filter2D_stats
public :: print_nu_filter2D_stats, build_nu2D_local_resolution_image
private
#include "simple_local_flags.inc"

integer, parameter :: NU2D_LABEL_KIND    = selected_int_kind(4)
integer, parameter :: NU2D_WEAK_LP_STEP  = 1

type :: nu_filter2D_stats
    real,    allocatable :: lowpass_limits(:)
    integer, allocatable :: label_counts_raw(:), label_counts_potts(:), stepdiff_counts(:)
    integer :: nbase=0, nlabels=0
    integer :: nclasses=0, npix_total=0
    integer :: potts_iters_total=0, potts_label_changes=0
    integer :: neighbor_pairs=0, identical_pairs=0, weak_step_pairs=0
    integer :: strong_jump_pairs=0, transition_pairs=0, pixels_strong_jump=0
end type nu_filter2D_stats

contains

    subroutine init_nu_filter2D_stats( stats, lowpass_limits )
        type(nu_filter2D_stats), intent(inout) :: stats
        real,                    intent(in)    :: lowpass_limits(:)
        integer :: max_step
        call kill_nu_filter2D_stats(stats)
        stats%nbase     = size(lowpass_limits)
        stats%nlabels   = stats%nbase
        if( stats%nbase < 1 ) THROW_HARD('empty 2D NU low-pass bank; init_nu_filter2D_stats')
        allocate(stats%lowpass_limits(stats%nbase), source=lowpass_limits)
        allocate(stats%label_counts_raw(stats%nlabels),   source=0)
        allocate(stats%label_counts_potts(stats%nlabels), source=0)
        max_step = max(1, stats%nlabels - 1)
        allocate(stats%stepdiff_counts(max_step), source=0)
    end subroutine init_nu_filter2D_stats

    subroutine kill_nu_filter2D_stats( stats )
        type(nu_filter2D_stats), intent(inout) :: stats
        if( allocated(stats%lowpass_limits)    ) deallocate(stats%lowpass_limits)
        if( allocated(stats%label_counts_raw)  ) deallocate(stats%label_counts_raw)
        if( allocated(stats%label_counts_potts)) deallocate(stats%label_counts_potts)
        if( allocated(stats%stepdiff_counts)   ) deallocate(stats%stepdiff_counts)
        stats%nbase = 0
        stats%nlabels = 0
        stats%nclasses = 0
        stats%npix_total = 0
        stats%potts_iters_total = 0
        stats%potts_label_changes = 0
        stats%neighbor_pairs = 0
        stats%identical_pairs = 0
        stats%weak_step_pairs = 0
        stats%strong_jump_pairs = 0
        stats%transition_pairs = 0
        stats%pixels_strong_jump = 0
    end subroutine kill_nu_filter2D_stats

    subroutine merge_nu_filter2D_stats( dst, src )
        type(nu_filter2D_stats), intent(inout) :: dst
        type(nu_filter2D_stats), intent(in)    :: src
        if( src%nlabels == 0 ) return
        if( dst%nlabels == 0 ) call init_nu_filter2D_stats(dst, src%lowpass_limits)
        if( dst%nlabels /= src%nlabels ) THROW_HARD('2D NU stats label-count mismatch; merge_nu_filter2D_stats')
        if( any(abs(dst%lowpass_limits - src%lowpass_limits) > TINY) ) &
            &THROW_HARD('2D NU stats low-pass bank mismatch; merge_nu_filter2D_stats')
        dst%label_counts_raw   = dst%label_counts_raw   + src%label_counts_raw
        dst%label_counts_potts = dst%label_counts_potts + src%label_counts_potts
        dst%stepdiff_counts    = dst%stepdiff_counts    + src%stepdiff_counts
        dst%nclasses           = dst%nclasses           + src%nclasses
        dst%npix_total         = dst%npix_total         + src%npix_total
        dst%potts_iters_total  = dst%potts_iters_total  + src%potts_iters_total
        dst%potts_label_changes = dst%potts_label_changes + src%potts_label_changes
        dst%neighbor_pairs        = dst%neighbor_pairs        + src%neighbor_pairs
        dst%identical_pairs       = dst%identical_pairs       + src%identical_pairs
        dst%weak_step_pairs       = dst%weak_step_pairs       + src%weak_step_pairs
        dst%strong_jump_pairs     = dst%strong_jump_pairs     + src%strong_jump_pairs
        dst%transition_pairs      = dst%transition_pairs      + src%transition_pairs
        dst%pixels_strong_jump    = dst%pixels_strong_jump    + src%pixels_strong_jump
    end subroutine merge_nu_filter2D_stats

    subroutine count_nu_filter2D_labels( labelmap, pix, nlabels, counts )
        integer(kind=NU2D_LABEL_KIND), intent(in)  :: labelmap(:,:,:)
        integer,                       intent(in)  :: pix(:,:)
        integer,                       intent(in)  :: nlabels
        integer,                       intent(out) :: counts(:)
        integer :: ipix, ilabel, i, j
        if( size(counts) /= nlabels ) THROW_HARD('counts size mismatch; count_nu_filter2D_labels')
        counts = 0
        !$omp parallel do schedule(static) default(shared) private(ipix,ilabel,i,j) reduction(+:counts) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            ilabel = int(labelmap(i,j,1))
            if( ilabel >= 1 .and. ilabel <= nlabels ) counts(ilabel) = counts(ilabel) + 1
        end do
        !$omp end parallel do
    end subroutine count_nu_filter2D_labels

    subroutine accumulate_nu_filter2D_stats( stats, labelmap, pix, raw_counts, potts_iters, potts_changes )
        type(nu_filter2D_stats),        intent(inout) :: stats
        integer(kind=NU2D_LABEL_KIND),  intent(in)    :: labelmap(:,:,:)
        integer,                        intent(in)    :: pix(:,:)
        integer,                        intent(in)    :: raw_counts(:)
        integer,                        intent(in)    :: potts_iters, potts_changes
        integer :: final_counts(stats%nlabels)
        if( stats%nlabels == 0 ) THROW_HARD('2D NU stats not initialized; accumulate_nu_filter2D_stats')
        if( size(raw_counts) /= stats%nlabels ) THROW_HARD('raw_counts size mismatch; accumulate_nu_filter2D_stats')
        call count_nu_filter2D_labels(labelmap, pix, stats%nlabels, final_counts)
        stats%nclasses = stats%nclasses + 1
        stats%npix_total            = stats%npix_total + size(pix,2)
        stats%label_counts_raw      = stats%label_counts_raw   + raw_counts
        stats%label_counts_potts    = stats%label_counts_potts + final_counts
        stats%potts_iters_total     = stats%potts_iters_total + potts_iters
        stats%potts_label_changes   = stats%potts_label_changes + potts_changes
        call accumulate_neighbor_continuity(stats, labelmap, pix)
    end subroutine accumulate_nu_filter2D_stats

    subroutine print_nu_filter2D_stats( stats )
        type(nu_filter2D_stats), intent(in) :: stats
        integer :: ilabel, nraw, nfinal
        real    :: raw_pct, final_pct, pct
        if( stats%nlabels == 0 ) return
        if( stats%npix_total == 0 )then
            write(logfhandle,'(A)') '>>> 2D NU FILTER: no pixels analyzed'
            return
        endif
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> 2D NU FILTER LOW-PASS ASSIGNMENTS'
        write(logfhandle,'(A,I12)') '    Class averages analyzed: ', stats%nclasses
        write(logfhandle,'(A,I12)') '    Pixels analyzed:         ', stats%npix_total
        write(logfhandle,'(A)') &
            &'    Source       Label  LP limit (A)       Raw    Raw %      Final  Final %'
        do ilabel = 1, stats%nbase
            nraw = stats%label_counts_raw(ilabel)
            nfinal = stats%label_counts_potts(ilabel)
            raw_pct = 100. * real(nraw) / real(stats%npix_total)
            final_pct = 100. * real(nfinal) / real(stats%npix_total)
            write(logfhandle,'(4X,A10,2X,I5,2X,F12.3,2X,I10,2X,F7.2,2X,I10,2X,F7.2)') &
                &'Base', ilabel, stats%lowpass_limits(ilabel), nraw, raw_pct, nfinal, final_pct
        end do
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> 2D NU POTTS PRIOR BEHAVIOR'
        write(logfhandle,'(A,I12)') '    Total ICM sweeps:         ', stats%potts_iters_total
        write(logfhandle,'(A,I12)') '    Label changes from Potts: ', stats%potts_label_changes
        if( stats%npix_total > 0 )then
            pct = 100. * real(stats%potts_label_changes) / real(stats%npix_total)
            write(logfhandle,'(A,F8.3,A)') '    Label changes per analyzed pixel: ', pct, '%'
        endif
        call print_neighbor_continuity(stats)
    end subroutine print_nu_filter2D_stats

    subroutine build_nu2D_local_resolution_image( labelmap, pix, lowpass_limits, ldim, smpd, resmap )
        integer(kind=NU2D_LABEL_KIND), intent(in)    :: labelmap(:,:,:)
        integer,                       intent(in)    :: pix(:,:)
        real,                          intent(in)    :: lowpass_limits(:)
        integer,                       intent(in)    :: ldim(3)
        real,                          intent(in)    :: smpd
        class(image),                  intent(inout) :: resmap
        real(kind=c_float), pointer :: rmat(:,:,:)
        integer :: ipix, i, j, ilabel
        call resmap%kill
        call resmap%new(ldim, smpd, wthreads=.false.)
        call resmap%get_rmat_ptr(rmat)
        rmat(:ldim(1),:ldim(2),1) = 0.
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            ilabel = int(labelmap(i,j,1))
            if( ilabel < 1 .or. ilabel > size(lowpass_limits) ) cycle
            if( lowpass_limits(ilabel) <= TINY ) cycle
            rmat(i,j,1) = 1. / lowpass_limits(ilabel)
        end do
    end subroutine build_nu2D_local_resolution_image

    subroutine accumulate_neighbor_continuity( stats, labelmap, pix )
        type(nu_filter2D_stats),       intent(inout) :: stats
        integer(kind=NU2D_LABEL_KIND), intent(in)    :: labelmap(:,:,:)
        integer,                       intent(in)    :: pix(:,:)
        integer :: ipix, i, j, di, dj, ni, nj, ilabel, jlabel, step_diff
        integer :: n_strong_jump_neighbors, local_pairs, local_identical, local_weak, local_strong
        integer :: local_transitions, local_pixels_strong
        integer :: local_steps(size(stats%stepdiff_counts))
        local_pairs = 0
        local_identical = 0
        local_weak = 0
        local_strong = 0
        local_transitions = 0
        local_pixels_strong = 0
        local_steps = 0
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            ilabel = int(labelmap(i,j,1))
            if( ilabel < 1 .or. ilabel > stats%nlabels ) cycle
            n_strong_jump_neighbors = 0
            do dj = -1, 1
                do di = -1, 1
                    if( di == 0 .and. dj == 0 ) cycle
                    ni = i + di
                    nj = j + dj
                    if( ni < 1 .or. ni > size(labelmap,1) ) cycle
                    if( nj < 1 .or. nj > size(labelmap,2) ) cycle
                    jlabel = int(labelmap(ni,nj,1))
                    if( jlabel < 1 .or. jlabel > stats%nlabels ) cycle
                    step_diff = abs(ilabel - jlabel)
                    if( step_diff > NU2D_WEAK_LP_STEP ) n_strong_jump_neighbors = n_strong_jump_neighbors + 1
                    if( .not.(dj > 0 .or. (dj == 0 .and. di > 0)) ) cycle
                    local_pairs = local_pairs + 1
                    if( step_diff == 0 )then
                        local_identical = local_identical + 1
                    else
                        local_transitions = local_transitions + 1
                        if( step_diff <= NU2D_WEAK_LP_STEP )then
                            local_weak = local_weak + 1
                        else
                            local_strong = local_strong + 1
                        endif
                        if( step_diff >= 1 .and. step_diff <= size(local_steps) ) &
                            &local_steps(step_diff) = local_steps(step_diff) + 1
                    endif
                end do
            end do
            if( n_strong_jump_neighbors > 0 ) local_pixels_strong = local_pixels_strong + 1
        end do
        stats%neighbor_pairs        = stats%neighbor_pairs        + local_pairs
        stats%identical_pairs       = stats%identical_pairs       + local_identical
        stats%weak_step_pairs       = stats%weak_step_pairs       + local_weak
        stats%strong_jump_pairs     = stats%strong_jump_pairs     + local_strong
        stats%transition_pairs      = stats%transition_pairs      + local_transitions
        stats%pixels_strong_jump    = stats%pixels_strong_jump    + local_pixels_strong
        stats%stepdiff_counts       = stats%stepdiff_counts       + local_steps
    end subroutine accumulate_neighbor_continuity

    subroutine print_neighbor_continuity( stats )
        type(nu_filter2D_stats), intent(in) :: stats
        integer :: istep
        real :: pct, pct_penalized
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> 2D NU NEIGHBOR CONTINUITY'
        write(logfhandle,'(A,I0,A)') '    Weakly penalized label transition: <= ', &
            &NU2D_WEAK_LP_STEP, ' coordinate step(s)'
        if( stats%neighbor_pairs == 0 )then
            write(logfhandle,'(A)') '    Continuity assessment: no neighbor pairs found'
            return
        endif
        pct = 100. * real(stats%pixels_strong_jump) / real(stats%npix_total)
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Pixels touching a strong jump: ', &
            &stats%pixels_strong_jump, ' (', pct, '%)'
        write(logfhandle,'(A,I14)') '    Neighbor pairs examined: ', stats%neighbor_pairs
        pct = 100. * real(stats%identical_pairs) / real(stats%neighbor_pairs)
        write(logfhandle,'(A,I14,A,F8.2,A)') '      identical:        ', stats%identical_pairs, ' (', pct, '%)'
        pct = 100. * real(stats%weak_step_pairs) / real(stats%neighbor_pairs)
        write(logfhandle,'(A,I14,A,F8.2,A)') '      weak step:        ', stats%weak_step_pairs, ' (', pct, '%)'
        pct_penalized = 100. * real(stats%strong_jump_pairs) / real(stats%neighbor_pairs)
        write(logfhandle,'(A,I14,A,F8.2,A)') '      strong jump:      ', stats%strong_jump_pairs, ' (', pct_penalized, '%)'
        pct = 100. * real(stats%transition_pairs) / real(stats%neighbor_pairs)
        write(logfhandle,'(A,I14,A,F8.2,A)') '      transition total: ', stats%transition_pairs, ' (', pct, '%)'
        write(logfhandle,'(A)') '    Label coordinate-step distribution:'
        pct = 100. * real(stats%identical_pairs) / real(stats%neighbor_pairs)
        write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') 0, stats%identical_pairs, pct, 'identical'
        do istep = 1, size(stats%stepdiff_counts)
            if( stats%stepdiff_counts(istep) == 0 ) cycle
            pct = 100. * real(stats%stepdiff_counts(istep)) / real(stats%neighbor_pairs)
            if( istep <= NU2D_WEAK_LP_STEP )then
                write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') &
                    &istep, stats%stepdiff_counts(istep), pct, 'weak step'
            else
                write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') &
                    &istep, stats%stepdiff_counts(istep), pct, 'strong jump'
            endif
        end do
        if( pct_penalized <= 1. )then
            write(logfhandle,'(A)') '    Continuity assessment: low strong-jump rate'
        else if( pct_penalized <= 5. )then
            write(logfhandle,'(A)') '    Continuity assessment: moderate strong-jump rate'
        else
            write(logfhandle,'(A)') '    Continuity assessment: high strong-jump rate'
        endif
        write(logfhandle,'(A)') ''
    end subroutine print_neighbor_continuity

end module simple_nu_filter2D_stats
