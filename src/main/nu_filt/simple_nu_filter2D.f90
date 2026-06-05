!@descr: 2D nonuniform low-pass filtering of even/odd class-average pairs
module simple_nu_filter2D
use simple_core_module_api
use simple_image,       only: image
use simple_butterworth, only: butterworth_filter
use simple_nu_filter2D_stats, only: NU2D_LABEL_KIND, nu_filter2D_stats, init_nu_filter2D_stats, &
    &kill_nu_filter2D_stats, count_nu_filter2D_labels, accumulate_nu_filter2D_stats
implicit none

public :: nu_filter2D_state, print_nu_filter2D_policy
private
#include "simple_local_flags.inc"

real,    parameter :: NU2D_LOW_PASS_LIMITS(8) = [30.,20.,15.,12.,8.,6.,5.,4.]
real,    parameter :: NU2D_TEST_AUX_LOW_PASS_LIMITS(2) = [20.,8.]
real,    parameter :: NU2D_TEST_AUX_EQUIV_LOW_PASS_LIMIT = 8.
logical, parameter :: SIMPLE_NU2D_TEST_AUX_BANK = .true.
character(len=*), parameter :: NU2D_TEST_AUX_BANK_ENVVAR = 'SIMPLE_NU2D_TEST_AUX_BANK'
real,    parameter :: NU2D_OBJECTIVE_SMOOTH_AWF          = 3.0
real,    parameter :: NU2D_OBJECTIVE_SMOOTH_RADIUS_FRAC  = 1.0
real,    parameter :: NU2D_OBJECTIVE_SMOOTH_MAX_RADIUS_A = 30.0
integer, parameter :: NU2D_LABEL_SMOOTH_MAXITS           = 8
integer, parameter :: NU2D_LABEL_SMOOTH_ADJACENT_STEPS  = 1
integer, parameter :: NU2D_LABEL_SMOOTH_RADIUS           = 2
integer, parameter :: NU2D_LABEL_SMOOTH_NNEIGH           = (2 * NU2D_LABEL_SMOOTH_RADIUS + 1)**2 - 1
integer, parameter :: NU2D_LABEL_SMOOTH_NCOLORS          = (NU2D_LABEL_SMOOTH_RADIUS + 1)**2
real,    parameter :: NU2D_LABEL_SMOOTH_BETA_FRAC        = 2.5
real,    parameter :: NU2D_LABEL_HIST_BETA_FRAC          = 8.0
real,    parameter :: NU2D_FOREGROUND_POTTS_BETA_SCALE   = 1.25
real,    parameter :: NU2D_BACKGROUND_POTTS_BETA_SCALE   = 1.0
real,    parameter :: NU2D_REGION_LOGODDS_BETA_FRAC      = 0.2
real,    parameter :: NU2D_REGION_PRIOR_PSEUDOCOUNT      = 0.5
real,    parameter :: NU2D_LABEL_SMOOTH_QUAD_FRAC        = 1.0
real,    parameter :: NU2D_LABEL_SMOOTH_TIE_EPS          = 1.e-6
real,    parameter :: NU2D_LABEL_SMOOTH_RING1_WEIGHT     = 1.0
real,    parameter :: NU2D_LABEL_SMOOTH_RING2_WEIGHT     = 0.5
real,    parameter :: NU2D_LABEL_SMOOTH_ADJACENT_FRAC    = 0.25
real,    parameter :: NU2D_SYNTH_LABEL_SMOOTH_RADIUS_A   = 10.0
integer, parameter :: NU2D_HIST_REGION_FOREGROUND        = 1
integer, parameter :: NU2D_HIST_REGION_BACKGROUND        = 2

type :: nu_filter2D_state
    real, allocatable :: dmats(:,:), candidate_coords(:), dmat(:,:,:), tmp(:,:,:), mask_tmp(:,:,:)
    real, allocatable :: bwfilters(:,:), lowpass_limits(:)
    integer, allocatable :: cutoff_finds(:), pix(:,:)
    integer(kind=NU2D_LABEL_KIND), allocatable :: filtmap(:,:,:)
    type(nu_filter2D_stats) :: stats
    integer :: ldim(3) = [0,0,0], box = 0, npix = 0
    real    :: smpd = 0.
    logical :: l_setup = .false., l_test_aux_bank = .false.
contains
    procedure :: setup => setup_nu_filter2D
    procedure :: apply => nu_filter2D_classavg
    procedure :: kill  => cleanup_nu_filter2D
end type nu_filter2D_state

contains

    subroutine print_nu_filter2D_policy()
        if( nu2D_test_aux_bank_enabled() )then
            write(logfhandle,'(A,A,A)') &
                &'>>> 2D nonuniform filter: TEST gate ', NU2D_TEST_AUX_BANK_ENVVAR, &
                &' active; candidate bank 20 A, 8 A, aux'
            write(logfhandle,'(A)') &
                &'>>> 2D nonuniform filter: aux candidate uses the supplied auxiliary even/odd pair'
            write(logfhandle,'(A)') &
                &'>>> 2D nonuniform filter: all three candidates compete over the whole image'
            write(logfhandle,'(A)') &
                &'>>> 2D nonuniform filter: Potts and histogram priors are disabled under this test gate'
            write(logfhandle,'(A,F4.1,A,F4.1,A)') &
                &'>>> 2D nonuniform filter: objective smoothing radius = ', &
                &NU2D_OBJECTIVE_SMOOTH_AWF * NU2D_OBJECTIVE_SMOOTH_RADIUS_FRAC, &
                &' * max(LP,smpd), capped at ', NU2D_OBJECTIVE_SMOOTH_MAX_RADIUS_A, ' A'
            write(logfhandle,'(A,F4.1,A)') &
                &'>>> 2D nonuniform filter: hard labels are tent-smoothed; output blends candidates over a ', &
                &NU2D_SYNTH_LABEL_SMOOTH_RADIUS_A, ' A tent-smoothed label field'
            return
        endif
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: low-pass bank 30,20,15,12,8,6,5,4 A; binary automask support'
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: all bank members compete over the whole image'
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: automask only defines foreground/background Potts regions'
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: automask failure filters every pixel with the lowest-resolution bank member'
        write(logfhandle,'(A,F4.1,A,F4.1,A)') &
            &'>>> 2D nonuniform filter: objective smoothing radius = ', &
            &NU2D_OBJECTIVE_SMOOTH_AWF * NU2D_OBJECTIVE_SMOOTH_RADIUS_FRAC, &
            &' * max(LP,smpd), capped at ', NU2D_OBJECTIVE_SMOOTH_MAX_RADIUS_A, ' A'
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: foreground/background log-odds prior uses raw label histograms'
        write(logfhandle,'(A,F4.2,A,F4.2)') &
            &'>>> 2D nonuniform filter: log-odds prior scale = ', NU2D_REGION_LOGODDS_BETA_FRAC, &
            &' * Potts beta; pseudocount = ', NU2D_REGION_PRIOR_PSEUDOCOUNT
        write(logfhandle,'(A)') &
            &'>>> 2D nonuniform filter: background histogram target is the lowest-resolution label'
        write(logfhandle,'(A,F4.2,A,F4.2)') &
            &'>>> 2D nonuniform filter: Potts beta scale foreground/background = ', &
            &NU2D_FOREGROUND_POTTS_BETA_SCALE, '/', NU2D_BACKGROUND_POTTS_BETA_SCALE
        write(logfhandle,'(A,I0,A)') &
            &'>>> 2D nonuniform filter: Potts ICM weights one- and two-pixel neighborhoods; weak step <= ', &
            &NU2D_LABEL_SMOOTH_ADJACENT_STEPS, ' label coordinate'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: Potts penalty includes weak 1-label jumps'
        write(logfhandle,'(A,F4.1,A)') &
            &'>>> 2D nonuniform filter: output blends bank members over a ', &
            &NU2D_SYNTH_LABEL_SMOOTH_RADIUS_A, ' A tent-smoothed label field'
    end subroutine print_nu_filter2D_policy

    subroutine setup_nu_filter2D( self, ldim_in, smpd_in )
        class(nu_filter2D_state), intent(inout) :: self
        integer,                   intent(in)    :: ldim_in(3)
        real,                      intent(in)    :: smpd_in
        logical, allocatable :: is_aux_label(:)
        integer :: nfilters, nlabels, icand
        call self%kill
        self%ldim = ldim_in
        self%box  = self%ldim(1)
        self%smpd = smpd_in
        self%l_test_aux_bank = nu2D_test_aux_bank_enabled()
        if( self%ldim(3) /= 1 ) THROW_HARD('2D images only; setup_nu_filter2D')
        if( self%ldim(1) /= self%ldim(2) ) THROW_HARD('square 2D images expected; setup_nu_filter2D')
        if( self%box < 2 ) THROW_HARD('invalid box; setup_nu_filter2D')
        if( self%smpd <= TINY ) THROW_HARD('invalid smpd; setup_nu_filter2D')
        call setup_pixels(self%ldim, self%pix, self%npix)
        if( self%l_test_aux_bank )then
            call setup_cutoff_finds(self%box, self%smpd, NU2D_TEST_AUX_LOW_PASS_LIMITS, self%cutoff_finds)
        else
            call setup_cutoff_finds(self%box, self%smpd, NU2D_LOW_PASS_LIMITS, self%cutoff_finds)
        endif
        nfilters = size(self%cutoff_finds)
        nlabels  = nfilters
        if( self%l_test_aux_bank ) nlabels = nlabels + 1
        call setup_candidate_coords(nlabels, self%candidate_coords)
        allocate(self%dmats(self%npix,nlabels), source=huge(0.))
        allocate(self%dmat(self%ldim(1),self%ldim(2),1))
        allocate(self%tmp(self%ldim(1),self%ldim(2),1))
        allocate(self%mask_tmp(self%ldim(1),self%ldim(2),1))
        allocate(self%filtmap(self%ldim(1),self%ldim(2),1), source=1_NU2D_LABEL_KIND)
        allocate(self%bwfilters(self%box,nfilters), source=0.)
        allocate(self%lowpass_limits(nlabels), source=0.)
        allocate(is_aux_label(nlabels), source=.false.)
        do icand = 1, nfilters
            call butterworth_filter(self%cutoff_finds(icand), self%bwfilters(:,icand))
            self%lowpass_limits(icand) = calc_lowpass_lim(self%cutoff_finds(icand), self%box, self%smpd)
        end do
        if( self%l_test_aux_bank )then
            self%lowpass_limits(nlabels) = NU2D_TEST_AUX_EQUIV_LOW_PASS_LIMIT
            is_aux_label(nlabels) = .true.
        endif
        call init_nu_filter2D_stats(self%stats, self%lowpass_limits, is_aux_label)
        deallocate(is_aux_label)
        self%l_setup = .true.
    end subroutine setup_nu_filter2D

    subroutine nu_filter2D_classavg( self, even_raw, odd_raw, aux_even, aux_odd, aux_avg, even_out, odd_out, &
            &avg_out, align_lp, locres_out, support_pix )
        class(nu_filter2D_state), intent(inout) :: self
        class(image),             intent(in)    :: even_raw, odd_raw, aux_even, aux_odd, aux_avg
        class(image),             intent(out)   :: even_out, odd_out, avg_out
        real,           optional, intent(out)   :: align_lp
        class(image),   optional, intent(inout) :: locres_out
        integer,        optional, intent(in)    :: support_pix(:,:)
        call validate_inputs()
        if( present(support_pix) )then
            call validate_support_pixels(support_pix)
            if( size(support_pix,2) > 0 )then
                call run_with_pixels(self%pix, support_pix)
            else if( self%l_test_aux_bank )then
                call run_with_pixels(self%pix)
            else
                call run_lowest_resolution_only()
            endif
        else
            call run_with_pixels(self%pix)
        endif

    contains

        subroutine run_with_pixels( eval_pix, hist_support_pix )
            integer, intent(in) :: eval_pix(:,:)
            integer, optional, intent(in) :: hist_support_pix(:,:)
            type(image), allocatable :: even_bank(:), odd_bank(:)
            type(image) :: even_ft, odd_ft
            integer, allocatable :: label_counts_raw(:)
            integer :: nfilters, nlabels, nout, icand, winsz, potts_iters, potts_changes, npix_eval
            real    :: noise_sigma, edge_mean
            npix_eval = size(eval_pix,2)
            if( npix_eval < 1 )then
                call run_lowest_resolution_only()
                return
            endif
            nfilters = size(self%cutoff_finds)
            nlabels = size(self%lowpass_limits)
            nout = nlabels
            self%dmats(:npix_eval,:nlabels) = huge(0.)
            self%filtmap = 0_NU2D_LABEL_KIND
            allocate(even_bank(nout), odd_bank(nout))
            allocate(label_counts_raw(nlabels), source=0)
            if( present(hist_support_pix) )then
                call build_nu2D_support_mask(hist_support_pix, self%mask_tmp, self%ldim)
                noise_sigma = nu2D_objective_noise_scale(even_raw, odd_raw, hist_support_pix)
            else
                noise_sigma = nu2D_objective_noise_scale(even_raw, odd_raw, eval_pix)
            endif
            winsz = nint(COSMSKHALFWIDTH)
            call even_ft%copy(even_raw)
            call odd_ft%copy(odd_raw)
            call even_ft%taper_edges_particle(winsz, edge_mean)
            call odd_ft%taper_edges_particle(winsz, edge_mean)
            call even_ft%fft
            call odd_ft%fft
            do icand = 1, nfilters
                call even_bank(icand)%new(self%ldim, self%smpd, wthreads=.false.)
                call odd_bank(icand)%new(self%ldim, self%smpd, wthreads=.false.)
                call even_bank(icand)%copy_fast(even_ft)
                call odd_bank(icand)%copy_fast(odd_ft)
                call even_bank(icand)%apply_filter(odd_bank(icand), self%bwfilters(:,icand))
                call even_bank(icand)%ifft
                call odd_bank(icand)%ifft
                call calc_nu2D_huber_objective(even_raw, even_bank(icand), odd_raw, odd_bank(icand), &
                    &noise_sigma, self%dmat, eval_pix)
                call smooth_nu2D_objective(self%dmat, self%tmp, self%mask_tmp, self%ldim, self%smpd, &
                    &calc_lowpass_lim(self%cutoff_finds(icand), self%box, self%smpd), eval_pix)
                call pack_nu2D_dmat(self%dmat, eval_pix, self%dmats(:npix_eval,icand))
            end do
            if( self%l_test_aux_bank )then
                icand = nlabels
                call even_bank(icand)%new(self%ldim, self%smpd, wthreads=.false.)
                call odd_bank(icand)%new(self%ldim, self%smpd, wthreads=.false.)
                call even_bank(icand)%copy(aux_even)
                call odd_bank(icand)%copy(aux_odd)
                call calc_nu2D_huber_objective(even_raw, even_bank(icand), odd_raw, odd_bank(icand), &
                    &noise_sigma, self%dmat, eval_pix)
                call smooth_nu2D_objective(self%dmat, self%tmp, self%mask_tmp, self%ldim, self%smpd, &
                    &NU2D_TEST_AUX_EQUIV_LOW_PASS_LIMIT, eval_pix)
                call pack_nu2D_dmat(self%dmat, eval_pix, self%dmats(:npix_eval,icand))
            endif
            call select_nu2D_labels(self%dmats(:npix_eval,:nlabels), eval_pix, self%filtmap)
            call count_nu_filter2D_labels(self%filtmap, eval_pix, nlabels, label_counts_raw)
            if( self%l_test_aux_bank )then
                potts_iters = 0
                potts_changes = 0
            else
                if( present(hist_support_pix) )then
                    call refine_nu2D_labels(self%filtmap, self%dmats(:npix_eval,:nlabels), eval_pix, &
                        &self%candidate_coords(:nlabels), self%ldim, potts_iters, potts_changes, region_mask=self%mask_tmp)
                else
                    call refine_nu2D_labels(self%filtmap, self%dmats(:npix_eval,:nlabels), eval_pix, &
                        &self%candidate_coords(:nlabels), self%ldim, potts_iters, potts_changes)
                endif
            endif
            call accumulate_nu_filter2D_stats(self%stats, self%filtmap, eval_pix, label_counts_raw, &
                &potts_iters, potts_changes)
            if( present(align_lp) ) align_lp = nu2D_finest_selected_lp(self%filtmap, eval_pix, self%lowpass_limits)
            call smooth_nu2D_label_coords(self%filtmap, eval_pix, self%dmat, self%tmp, self%mask_tmp, &
                &self%ldim, self%smpd, nout)
            if( present(locres_out) ) call build_nu2D_soft_local_resolution_image(self%dmat, self%pix, &
                &self%lowpass_limits, self%ldim, self%smpd, locres_out)
            call synthesize_nu2D_soft(even_bank, odd_bank, self%dmat, self%pix, self%ldim, self%smpd, &
                &even_out, odd_out, avg_out)
            do icand = 1, nout
                call even_bank(icand)%kill
                call odd_bank(icand)%kill
            end do
            call even_ft%kill
            call odd_ft%kill
            deallocate(even_bank, odd_bank)
            deallocate(label_counts_raw)
        end subroutine run_with_pixels

        subroutine run_lowest_resolution_only()
            type(image) :: even_ft, odd_ft
            type(image), allocatable :: even_bank(:), odd_bank(:)
            integer :: winsz
            real    :: edge_mean
            if( present(align_lp) ) align_lp = 0.
            allocate(even_bank(1), odd_bank(1))
            winsz = nint(COSMSKHALFWIDTH)
            call even_ft%copy(even_raw)
            call odd_ft%copy(odd_raw)
            call even_ft%taper_edges_particle(winsz, edge_mean)
            call odd_ft%taper_edges_particle(winsz, edge_mean)
            call even_ft%fft
            call odd_ft%fft
            call even_bank(1)%new(self%ldim, self%smpd, wthreads=.false.)
            call odd_bank(1)%new(self%ldim, self%smpd, wthreads=.false.)
            call even_bank(1)%copy_fast(even_ft)
            call odd_bank(1)%copy_fast(odd_ft)
            call even_bank(1)%apply_filter(odd_bank(1), self%bwfilters(:,1))
            call even_bank(1)%ifft
            call odd_bank(1)%ifft
            self%dmat(:self%ldim(1),:self%ldim(2),1) = 1.
            if( present(locres_out) ) call build_nu2D_soft_local_resolution_image(self%dmat, self%pix, &
                &self%lowpass_limits(:1), self%ldim, self%smpd, locres_out)
            call synthesize_nu2D_soft(even_bank, odd_bank, self%dmat, self%pix, self%ldim, self%smpd, &
                &even_out, odd_out, avg_out)
            call even_bank(1)%kill
            call odd_bank(1)%kill
            call even_ft%kill
            call odd_ft%kill
            deallocate(even_bank, odd_bank)
        end subroutine run_lowest_resolution_only

        subroutine validate_inputs()
            integer :: img_ldim(3)
            real    :: img_smpd
            if( .not.self%l_setup ) THROW_HARD('setup must be called first; nu_filter2D_classavg')
            img_ldim = even_raw%get_ldim()
            img_smpd = even_raw%get_smpd()
            if( even_raw%is_3d() ) THROW_HARD('2D images only; nu_filter2D_classavg')
            if( any(img_ldim /= self%ldim) ) THROW_HARD('even_raw dimensions differ from setup; nu_filter2D_classavg')
            if( any(odd_raw%get_ldim() /= self%ldim) ) THROW_HARD('odd_raw dimensions differ; nu_filter2D_classavg')
            if( any(aux_even%get_ldim() /= self%ldim) ) THROW_HARD('aux_even dimensions differ; nu_filter2D_classavg')
            if( any(aux_odd%get_ldim()  /= self%ldim) ) THROW_HARD('aux_odd dimensions differ; nu_filter2D_classavg')
            if( any(aux_avg%get_ldim()  /= self%ldim) ) THROW_HARD('aux_avg dimensions differ; nu_filter2D_classavg')
            if( abs(img_smpd - self%smpd) > TINY ) THROW_HARD('image smpd differs from setup; nu_filter2D_classavg')
            if( abs(odd_raw%get_smpd() - img_smpd) > TINY ) THROW_HARD('odd_raw smpd differs; nu_filter2D_classavg')
            if( abs(aux_even%get_smpd() - img_smpd) > TINY ) THROW_HARD('aux_even smpd differs; nu_filter2D_classavg')
            if( abs(aux_odd%get_smpd()  - img_smpd) > TINY ) THROW_HARD('aux_odd smpd differs; nu_filter2D_classavg')
            if( abs(aux_avg%get_smpd()  - img_smpd) > TINY ) THROW_HARD('aux_avg smpd differs; nu_filter2D_classavg')
        end subroutine validate_inputs

        subroutine validate_support_pixels( pix )
            integer, intent(in) :: pix(:,:)
            integer :: ipix
            if( size(pix,1) /= 2 ) THROW_HARD('support pixel list must have two rows; nu_filter2D_classavg')
            do ipix = 1, size(pix,2)
                if( pix(1,ipix) < 1 .or. pix(1,ipix) > self%ldim(1) ) &
                    &THROW_HARD('support pixel row out of bounds; nu_filter2D_classavg')
                if( pix(2,ipix) < 1 .or. pix(2,ipix) > self%ldim(2) ) &
                    &THROW_HARD('support pixel column out of bounds; nu_filter2D_classavg')
            end do
        end subroutine validate_support_pixels

    end subroutine nu_filter2D_classavg

    subroutine cleanup_nu_filter2D( self )
        class(nu_filter2D_state), intent(inout) :: self
        if( allocated(self%dmats)            ) deallocate(self%dmats)
        if( allocated(self%candidate_coords) ) deallocate(self%candidate_coords)
        if( allocated(self%dmat)             ) deallocate(self%dmat)
        if( allocated(self%tmp)              ) deallocate(self%tmp)
        if( allocated(self%mask_tmp)         ) deallocate(self%mask_tmp)
        if( allocated(self%bwfilters)        ) deallocate(self%bwfilters)
        if( allocated(self%lowpass_limits)   ) deallocate(self%lowpass_limits)
        if( allocated(self%cutoff_finds)     ) deallocate(self%cutoff_finds)
        if( allocated(self%pix)              ) deallocate(self%pix)
        if( allocated(self%filtmap)          ) deallocate(self%filtmap)
        call kill_nu_filter2D_stats(self%stats)
        self%ldim = [0,0,0]
        self%box = 0
        self%npix = 0
        self%smpd = 0.
        self%l_setup = .false.
        self%l_test_aux_bank = .false.
    end subroutine cleanup_nu_filter2D

    subroutine setup_pixels( ldim, pix, npix )
        integer, intent(in) :: ldim(3)
        integer, allocatable, intent(out) :: pix(:,:)
        integer, intent(out) :: npix
        integer :: i, j, ipix
        npix = ldim(1) * ldim(2)
        if( npix < 1 ) THROW_HARD('empty image; setup_pixels')
        allocate(pix(2,npix), source=0)
        ipix = 0
        do j = 1, ldim(2)
            do i = 1, ldim(1)
                ipix = ipix + 1
                pix(:,ipix) = [i,j]
            end do
        end do
        if( ipix /= npix ) THROW_HARD('pixel count mismatch; setup_pixels')
    end subroutine setup_pixels

    subroutine setup_cutoff_finds( box, smpd, lp_limits, cutoff_finds )
        integer, intent(in) :: box
        real,    intent(in) :: smpd, lp_limits(:)
        integer, allocatable, intent(out) :: cutoff_finds(:)
        integer :: i, find, nvalid
        integer :: tmp(size(lp_limits))
        nvalid = 0
        do i = 1, size(lp_limits)
            find = calc_fourier_index(lp_limits(i), box, smpd)
            find = max(1, min(box/2, find))
            if( nvalid > 0 )then
                if( any(tmp(:nvalid) == find) ) cycle
            endif
            nvalid = nvalid + 1
            tmp(nvalid) = find
        end do
        if( nvalid < 1 ) THROW_HARD('empty cutoff bank; setup_cutoff_finds')
        allocate(cutoff_finds(nvalid))
        cutoff_finds = tmp(:nvalid)
    end subroutine setup_cutoff_finds

    subroutine setup_candidate_coords( nbase, coords )
        integer, intent(in) :: nbase
        real, allocatable, intent(out) :: coords(:)
        integer :: i
        allocate(coords(nbase), source=0.)
        do i = 1, nbase
            coords(i) = real(i)
        end do
    end subroutine setup_candidate_coords

    logical function nu2D_test_aux_bank_enabled()
        character(len=32) :: env_value
        integer :: status
        nu2D_test_aux_bank_enabled = SIMPLE_NU2D_TEST_AUX_BANK
        call get_environment_variable(NU2D_TEST_AUX_BANK_ENVVAR, env_value, status=status)
        if( status /= 0 ) return
        select case( trim(adjustl(env_value)) )
        case( '1', 'y', 'Y', 'yes', 'Yes', 'YES', 'true', 'True', 'TRUE', 'on', 'On', 'ON' )
            nu2D_test_aux_bank_enabled = .true.
        case( '0', 'n', 'N', 'no', 'No', 'NO', 'false', 'False', 'FALSE', 'off', 'Off', 'OFF' )
            nu2D_test_aux_bank_enabled = .false.
        end select
    end function nu2D_test_aux_bank_enabled

    real function nu2D_objective_noise_scale( even_raw, odd_raw, pix )
        class(image), intent(in) :: even_raw, odd_raw
        integer,      intent(in) :: pix(:,:)
        real(kind=c_float), pointer :: r_even(:,:,:), r_odd(:,:,:)
        real, allocatable :: vals(:)
        real :: med
        integer :: i, j, ipix, npix
        npix = size(pix,2)
        if( npix < 1 ) THROW_HARD('empty image; nu2D_objective_noise_scale')
        allocate(vals(npix))
        call even_raw%get_rmat_ptr(r_even)
        call odd_raw%get_rmat_ptr(r_odd)
        do ipix = 1, npix
            i = pix(1,ipix)
            j = pix(2,ipix)
            vals(ipix) = real(r_even(i,j,1) - r_odd(i,j,1))
        end do
        med = median_nocopy(vals)
        nu2D_objective_noise_scale = mad_gau(vals, med)
        if( nu2D_objective_noise_scale <= TINY )then
            nu2D_objective_noise_scale = sqrt(sum(vals * vals) / real(npix))
        endif
        if( nu2D_objective_noise_scale <= TINY ) nu2D_objective_noise_scale = 1.
        deallocate(vals)
    end function nu2D_objective_noise_scale

    subroutine calc_nu2D_huber_objective( even_raw, even_cand, odd_raw, odd_cand, noise_sigma, dmat, pix )
        class(image), intent(in)  :: even_raw, even_cand, odd_raw, odd_cand
        real,         intent(in)  :: noise_sigma
        integer,      intent(in)  :: pix(:,:)
        real,         intent(out) :: dmat(:,:,:)
        real(kind=c_float), pointer :: r_eraw(:,:,:), r_ecand(:,:,:), r_oraw(:,:,:), r_ocand(:,:,:)
        real :: sigma, r1, r2
        integer :: i, j, nx, ny, ipix
        nx = size(dmat,1)
        ny = size(dmat,2)
        sigma = max(noise_sigma, TINY)
        call even_raw%get_rmat_ptr(r_eraw)
        call even_cand%get_rmat_ptr(r_ecand)
        call odd_raw%get_rmat_ptr(r_oraw)
        call odd_cand%get_rmat_ptr(r_ocand)
        dmat(:nx,:ny,1) = 0.
        if( size(pix,2) == nx * ny )then
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,r1,r2) proc_bind(close)
            do j = 1, ny
                do i = 1, nx
                    r1 = real(r_eraw(i,j,1) - r_ocand(i,j,1)) / sigma
                    r2 = real(r_ecand(i,j,1) - r_oraw(i,j,1)) / sigma
                    dmat(i,j,1) = huber_loss(r1) + huber_loss(r2)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) private(ipix,i,j,r1,r2) proc_bind(close)
            do ipix = 1, size(pix,2)
                i = pix(1,ipix)
                j = pix(2,ipix)
                r1 = real(r_eraw(i,j,1) - r_ocand(i,j,1)) / sigma
                r2 = real(r_ecand(i,j,1) - r_oraw(i,j,1)) / sigma
                dmat(i,j,1) = huber_loss(r1) + huber_loss(r2)
            end do
            !$omp end parallel do
        endif
    contains
        elemental real function huber_loss( r ) result( loss )
            real, intent(in) :: r
            real, parameter :: HUBER_DELTA = 1.345
            real :: ar
            ar = abs(r)
            if( ar <= HUBER_DELTA )then
                loss = 0.5 * r * r
            else
                loss = HUBER_DELTA * (ar - 0.5 * HUBER_DELTA)
            endif
        end function huber_loss
    end subroutine calc_nu2D_huber_objective

    subroutine smooth_nu2D_objective( dmat, tmp, mask_tmp, ldim, smpd, lp_angstrom, pix )
        real,    intent(inout) :: dmat(:,:,:), tmp(:,:,:), mask_tmp(:,:,:)
        integer, intent(in)    :: ldim(3), pix(:,:)
        real,    intent(in)    :: smpd, lp_angstrom
        integer :: radius_px, ipix, i, j
        real :: radius_angstrom, denom
        radius_angstrom = NU2D_OBJECTIVE_SMOOTH_RADIUS_FRAC * NU2D_OBJECTIVE_SMOOTH_AWF * max(lp_angstrom, smpd)
        if( NU2D_OBJECTIVE_SMOOTH_MAX_RADIUS_A > TINY ) radius_angstrom = min(radius_angstrom, NU2D_OBJECTIVE_SMOOTH_MAX_RADIUS_A)
        radius_px = max(1, nint(max(smpd, radius_angstrom) / smpd))
        if( size(pix,2) == ldim(1) * ldim(2) )then
            call tent_smooth_2d(dmat(:,:,1), tmp(:,:,1), ldim(1), ldim(2), radius_px)
            return
        endif
        mask_tmp(:ldim(1),:ldim(2),1) = 0.
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            mask_tmp(i,j,1) = 1.
        end do
        call tent_smooth_2d(dmat(:,:,1), tmp(:,:,1), ldim(1), ldim(2), radius_px)
        call tent_smooth_2d(mask_tmp(:,:,1), tmp(:,:,1), ldim(1), ldim(2), radius_px)
        !$omp parallel do schedule(static) default(shared) private(ipix,i,j,denom) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            denom = max(mask_tmp(i,j,1), TINY)
            dmat(i,j,1) = dmat(i,j,1) / denom
        end do
        !$omp end parallel do
    end subroutine smooth_nu2D_objective

    subroutine tent_smooth_2d( img, tmp, nx, ny, radius )
        integer, intent(in)    :: nx, ny, radius
        real,    intent(inout) :: img(nx,ny)
        real,    intent(inout) :: tmp(nx,ny)
        integer :: box_width, lext1, rext1, lext2, rext2
        box_width = radius + 1
        lext1 = (box_width - 1) / 2
        rext1 = (box_width - 1) - lext1
        lext2 = rext1
        rext2 = lext1
        call box_filter_x_2d(img, tmp, nx, ny, lext1, rext1)
        call box_filter_x_2d(tmp, img, nx, ny, lext2, rext2)
        call box_filter_y_2d(img, tmp, nx, ny, lext1, rext1)
        call box_filter_y_2d(tmp, img, nx, ny, lext2, rext2)
    end subroutine tent_smooth_2d

    subroutine box_filter_x_2d( src, dst, nx, ny, lext, rext )
        integer, intent(in)  :: nx, ny, lext, rext
        real,    intent(in)  :: src(nx,ny)
        real,    intent(out) :: dst(nx,ny)
        integer :: j
        !$omp parallel do schedule(static) default(shared) private(j) proc_bind(close)
        do j = 1, ny
            call box_filter_1d(src(:,j), dst(:,j), nx, lext, rext)
        end do
        !$omp end parallel do
    end subroutine box_filter_x_2d

    subroutine box_filter_y_2d( src, dst, nx, ny, lext, rext )
        integer, intent(in)  :: nx, ny, lext, rext
        real,    intent(in)  :: src(nx,ny)
        real,    intent(out) :: dst(nx,ny)
        integer :: i
        real :: buf_in(ny), buf_out(ny)
        !$omp parallel do schedule(static) default(shared) private(i,buf_in,buf_out) proc_bind(close)
        do i = 1, nx
            buf_in = src(i,:)
            call box_filter_1d(buf_in, buf_out, ny, lext, rext)
            dst(i,:) = buf_out
        end do
        !$omp end parallel do
    end subroutine box_filter_y_2d

    subroutine box_filter_1d( src, dst, n, lext, rext )
        integer, intent(in)  :: n, lext, rext
        real,    intent(in)  :: src(n)
        real,    intent(out) :: dst(n)
        integer :: i, j, box_width
        real :: running_sum
        box_width = lext + rext + 1
        running_sum = 0.
        do j = -lext, rext
            running_sum = running_sum + src(clamp_index(1 + j, 1, n))
        end do
        dst(1) = running_sum / real(box_width)
        do i = 2, n
            running_sum = running_sum &
                &+ src(clamp_index(i + rext,     1, n)) &
                &- src(clamp_index(i - 1 - lext, 1, n))
            dst(i) = running_sum / real(box_width)
        end do
    end subroutine box_filter_1d

    pure integer function clamp_index( i, lo, hi )
        integer, intent(in) :: i, lo, hi
        clamp_index = max(lo, min(hi, i))
    end function clamp_index

    subroutine pack_nu2D_dmat( dmat, pix, dvec )
        real,    intent(in)  :: dmat(:,:,:)
        integer, intent(in)  :: pix(:,:)
        real,    intent(out) :: dvec(:)
        integer :: ipix, i, j
        !$omp parallel do schedule(static) default(shared) private(ipix,i,j) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            dvec(ipix) = dmat(i,j,1)
        end do
        !$omp end parallel do
    end subroutine pack_nu2D_dmat

    subroutine build_nu2D_support_mask( pix, mask, ldim )
        integer, intent(in)    :: pix(:,:), ldim(3)
        real,    intent(inout) :: mask(:,:,:)
        integer :: ipix, i, j
        if( any(shape(mask) /= ldim) ) THROW_HARD('support mask shape mismatch; build_nu2D_support_mask')
        mask(:ldim(1),:ldim(2),1) = 0.
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            mask(i,j,1) = 1.
        end do
    end subroutine build_nu2D_support_mask

    subroutine smooth_nu2D_label_coords( labelmap, pix, label_coord, work, mask_work, ldim, smpd, nlabels )
        integer(kind=NU2D_LABEL_KIND), intent(in)    :: labelmap(:,:,:)
        integer,                       intent(in)    :: pix(:,:), ldim(3), nlabels
        real,                          intent(in)    :: smpd
        real,                          intent(inout) :: label_coord(:,:,:), work(:,:,:), mask_work(:,:,:)
        integer :: ipix, i, j, ilabel, radius_px
        if( nlabels < 1 ) THROW_HARD('empty label bank; smooth_nu2D_label_coords')
        if( smpd <= TINY ) THROW_HARD('invalid smpd; smooth_nu2D_label_coords')
        radius_px = max(1, nint(NU2D_SYNTH_LABEL_SMOOTH_RADIUS_A / smpd))
        if( size(pix,2) == ldim(1) * ldim(2) )then
            label_coord(:ldim(1),:ldim(2),1) = 1.
            mask_work(:ldim(1),:ldim(2),1) = 1.
            !$omp parallel do schedule(static) default(shared) private(ipix,i,j,ilabel) proc_bind(close)
            do ipix = 1, size(pix,2)
                i = pix(1,ipix)
                j = pix(2,ipix)
                ilabel = max(1, min(nlabels, int(labelmap(i,j,1))))
                label_coord(i,j,1) = real(ilabel)
            end do
            !$omp end parallel do
            call tent_smooth_2d(label_coord(:,:,1), work(:,:,1), ldim(1), ldim(2), radius_px)
            !$omp parallel do schedule(static) default(shared) private(ipix,i,j) proc_bind(close)
            do ipix = 1, size(pix,2)
                i = pix(1,ipix)
                j = pix(2,ipix)
                label_coord(i,j,1) = max(1., min(real(nlabels), label_coord(i,j,1)))
            end do
            !$omp end parallel do
            return
        endif
        label_coord(:ldim(1),:ldim(2),1) = 0.
        mask_work(:ldim(1),:ldim(2),1) = 0.
        !$omp parallel do schedule(static) default(shared) private(ipix,i,j,ilabel) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            ilabel = max(1, min(nlabels, int(labelmap(i,j,1))))
            label_coord(i,j,1) = real(ilabel)
            mask_work(i,j,1) = 1.
        end do
        !$omp end parallel do
        call tent_smooth_2d(label_coord(:,:,1), work(:,:,1), ldim(1), ldim(2), radius_px)
        call tent_smooth_2d(mask_work(:,:,1), work(:,:,1), ldim(1), ldim(2), radius_px)
        work(:ldim(1),:ldim(2),1) = 0.
        !$omp parallel do schedule(static) default(shared) private(ipix,i,j) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            label_coord(i,j,1) = max(1., min(real(nlabels), label_coord(i,j,1) / max(mask_work(i,j,1), TINY)))
            work(i,j,1) = max(0., min(1., mask_work(i,j,1)))
        end do
        !$omp end parallel do
        label_coord(:ldim(1),:ldim(2),1) = merge(label_coord(:ldim(1),:ldim(2),1), 1., work(:ldim(1),:ldim(2),1) > 0.)
        mask_work(:ldim(1),:ldim(2),1) = work(:ldim(1),:ldim(2),1)
    end subroutine smooth_nu2D_label_coords

    subroutine synthesize_nu2D_soft( even_bank, odd_bank, label_coord, pix, ldim, smpd, even_out, odd_out, avg_out )
        type(image), intent(inout) :: even_bank(:), odd_bank(:)
        real,        intent(in)    :: label_coord(:,:,:), smpd
        integer,     intent(in)    :: pix(:,:), ldim(3)
        class(image), intent(out)  :: even_out, odd_out, avg_out
        real(kind=c_float), pointer :: r_even(:,:,:), r_odd(:,:,:), r_avg(:,:,:)
        real(kind=c_float), pointer :: r_be(:,:,:), r_bo(:,:,:)
        integer :: ncand, icand, ipix, i, j
        real    :: weight
        ncand = size(even_bank)
        if( ncand < 1 ) THROW_HARD('empty filter bank; synthesize_nu2D_soft')
        if( size(odd_bank) /= ncand ) THROW_HARD('bank size mismatch; synthesize_nu2D_soft')
        call even_out%new(ldim, smpd, wthreads=.false.)
        call odd_out%new(ldim, smpd, wthreads=.false.)
        call avg_out%new(ldim, smpd, wthreads=.false.)
        call even_out%get_rmat_ptr(r_even)
        call odd_out%get_rmat_ptr(r_odd)
        call avg_out%get_rmat_ptr(r_avg)
        r_even(:ldim(1),:ldim(2),1) = 0.
        r_odd (:ldim(1),:ldim(2),1) = 0.
        r_avg (:ldim(1),:ldim(2),1) = 0.
        do icand = 1, ncand
            call even_bank(icand)%get_rmat_ptr(r_be)
            call odd_bank(icand)%get_rmat_ptr(r_bo)
            !$omp parallel do schedule(static) default(shared) private(ipix,i,j,weight) proc_bind(close)
            do ipix = 1, size(pix,2)
                i = pix(1,ipix)
                j = pix(2,ipix)
                weight = nu2D_soft_label_weight(label_coord(i,j,1), icand)
                if( weight > TINY )then
                    r_even(i,j,1) = r_even(i,j,1) + weight * r_be(i,j,1)
                    r_odd (i,j,1) = r_odd (i,j,1) + weight * r_bo(i,j,1)
                    r_avg (i,j,1) = r_avg (i,j,1) + 0.5 * weight * (r_be(i,j,1) + r_bo(i,j,1))
                endif
            end do
            !$omp end parallel do
        end do
    end subroutine synthesize_nu2D_soft

    subroutine build_nu2D_soft_local_resolution_image( label_coord, pix, lowpass_limits, ldim, smpd, resmap )
        real,         intent(in)    :: label_coord(:,:,:), lowpass_limits(:), smpd
        integer,      intent(in)    :: pix(:,:), ldim(3)
        class(image), intent(inout) :: resmap
        real(kind=c_float), pointer :: rmat(:,:,:)
        real    :: coord, frac, lp
        integer :: ipix, i, j, ilow, ihigh, nlabels
        nlabels = size(lowpass_limits)
        if( nlabels < 1 ) THROW_HARD('empty filter bank; build_nu2D_soft_local_resolution_image')
        call resmap%kill
        call resmap%new(ldim, smpd, wthreads=.false.)
        call resmap%get_rmat_ptr(rmat)
        rmat(:ldim(1),:ldim(2),1) = 0.
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            coord = max(1., min(real(nlabels), label_coord(i,j,1)))
            ilow  = max(1, min(nlabels, int(coord)))
            ihigh = min(nlabels, ilow + 1)
            frac  = coord - real(ilow)
            if( ilow == ihigh ) frac = 0.
            lp = (1. - frac) * lowpass_limits(ilow) + frac * lowpass_limits(ihigh)
            if( lp > TINY ) rmat(i,j,1) = 1. / lp
        end do
    end subroutine build_nu2D_soft_local_resolution_image

    real function nu2D_soft_label_weight( label_coord, icand )
        real,    intent(in) :: label_coord
        integer, intent(in) :: icand
        nu2D_soft_label_weight = max(0., 1. - abs(label_coord - real(icand)))
    end function nu2D_soft_label_weight

    subroutine select_nu2D_labels( dmats, pix, labelmap )
        real,    intent(in) :: dmats(:,:)
        integer, intent(in) :: pix(:,:)
        integer(kind=NU2D_LABEL_KIND), intent(inout) :: labelmap(:,:,:)
        integer :: ipix, icand, cur_icand, i, j
        real    :: cur_dmat
        !$omp parallel do schedule(static) default(shared) private(ipix,icand,cur_icand,cur_dmat,i,j) proc_bind(close)
        do ipix = 1, size(pix,2)
            cur_icand = 1
            cur_dmat  = huge(cur_dmat)
            do icand = 1, size(dmats,2)
                if( dmats(ipix,icand) < cur_dmat )then
                    cur_dmat  = dmats(ipix,icand)
                    cur_icand = icand
                endif
            end do
            i = pix(1,ipix)
            j = pix(2,ipix)
            labelmap(i,j,1) = int(cur_icand, kind=NU2D_LABEL_KIND)
        end do
        !$omp end parallel do
    end subroutine select_nu2D_labels

    real function nu2D_finest_selected_lp( labelmap, pix, lowpass_limits )
        integer(kind=NU2D_LABEL_KIND), intent(in) :: labelmap(:,:,:)
        integer, intent(in) :: pix(:,:)
        real,    intent(in) :: lowpass_limits(:)
        integer :: ipix, i, j, ilabel
        real :: finest_lp, lp
        finest_lp = huge(finest_lp)
        !$omp parallel do schedule(static) default(shared) private(ipix,i,j,ilabel,lp) reduction(min:finest_lp) proc_bind(close)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            ilabel = int(labelmap(i,j,1))
            if( ilabel >= 1 .and. ilabel <= size(lowpass_limits) )then
                lp = lowpass_limits(ilabel)
                if( lp > TINY ) finest_lp = min(finest_lp, lp)
            endif
        end do
        !$omp end parallel do
        nu2D_finest_selected_lp = 0.
        if( finest_lp < huge(finest_lp) ) nu2D_finest_selected_lp = finest_lp
    end function nu2D_finest_selected_lp

    subroutine refine_nu2D_labels( labelmap, dmats, pix, coords, ldim, potts_iters, potts_changes, region_mask )
        integer(kind=NU2D_LABEL_KIND), intent(inout) :: labelmap(:,:,:)
        real,                          intent(in)    :: dmats(:,:), coords(:)
        integer,                       intent(in)    :: pix(:,:), ldim(3)
        integer,                       intent(out)   :: potts_iters, potts_changes
        real, optional,                intent(in)    :: region_mask(:,:,:)
        integer, allocatable :: counts(:,:), target_counts(:,:), region_sizes(:), pixel_regions(:)
        real,    allocatable :: hist_scales(:), potts_scales(:), region_priors(:,:)
        integer :: iter, color, ipix, icand, cur_icand, best_icand, i, j, neigh(3,NU2D_LABEL_SMOOTH_NNEIGH), nsz
        integer :: nchanged, ncand, iregion
        real    :: beta, region_beta, e, best_e, neigh_w(NU2D_LABEL_SMOOTH_NNEIGH)
        potts_iters   = 0
        potts_changes = 0
        ncand = size(dmats,2)
        if( ncand < 2 ) return
        if( present(region_mask) )then
            if( any(shape(region_mask) /= ldim) ) THROW_HARD('region mask shape mismatch; refine_nu2D_labels')
        endif
        beta = estimate_nu2D_beta(dmats)
        if( beta <= TINY ) return
        call init_nu2D_histogram_targets(labelmap, pix, ncand, beta, target_counts, counts, &
            &region_sizes, hist_scales, pixel_regions, region_priors, region_mask)
        call init_nu2D_potts_scales(potts_scales, region_mask)
        do iter = 1, NU2D_LABEL_SMOOTH_MAXITS
            potts_iters = potts_iters + 1
            nchanged = 0
            do color = 0, NU2D_LABEL_SMOOTH_NCOLORS - 1
                do ipix = 1, size(pix,2)
                    i = pix(1,ipix)
                    j = pix(2,ipix)
                    if( nu2D_label_color(i,j) /= color ) cycle
                    call collect_nu2D_neighbors(ldim, i, j, neigh, neigh_w, nsz)
                    cur_icand = max(1, min(ncand, int(labelmap(i,j,1))))
                    if( cur_icand /= int(labelmap(i,j,1)) ) labelmap(i,j,1) = int(cur_icand, kind=NU2D_LABEL_KIND)
                    iregion = pixel_regions(ipix)
                    region_beta = beta * potts_scales(iregion)
                    best_icand = cur_icand
                    best_e     = dmats(ipix,cur_icand) + region_priors(cur_icand,iregion) + region_beta * &
                        &nu2D_neighborhood_cost(cur_icand, labelmap, neigh, neigh_w, nsz, coords)
                    do icand = 1, ncand
                        if( icand == cur_icand ) cycle
                        e = dmats(ipix,icand) + region_priors(icand,iregion) + region_beta * &
                            &nu2D_neighborhood_cost(icand, labelmap, neigh, neigh_w, nsz, coords) + &
                            &nu2D_histogram_delta(icand, cur_icand, counts, target_counts, &
                            &iregion, hist_scales(iregion))
                        if( nu2D_is_better(e, best_e) )then
                            best_e     = e
                            best_icand = icand
                        endif
                    end do
                    if( best_icand /= cur_icand )then
                        nchanged = nchanged + 1
                        counts(cur_icand,iregion)  = counts(cur_icand,iregion) - 1
                        counts(best_icand,iregion) = counts(best_icand,iregion) + 1
                        labelmap(i,j,1) = int(best_icand, kind=NU2D_LABEL_KIND)
                    endif
                end do
            end do
            potts_changes = potts_changes + nchanged
            if( nchanged == 0 ) exit
        end do
        deallocate(counts, target_counts, region_sizes, hist_scales, pixel_regions, potts_scales, region_priors)
    end subroutine refine_nu2D_labels

    subroutine init_nu2D_histogram_targets( labelmap, pix, ncand, beta, target_counts, counts, &
            &region_sizes, hist_scales, pixel_regions, region_priors, region_mask )
        integer(kind=NU2D_LABEL_KIND), intent(in) :: labelmap(:,:,:)
        integer, intent(in) :: pix(:,:), ncand
        real,    intent(in) :: beta
        integer, allocatable, intent(out) :: target_counts(:,:), counts(:,:), region_sizes(:)
        integer, allocatable, intent(out) :: pixel_regions(:)
        real,    allocatable, intent(out) :: hist_scales(:), region_priors(:,:)
        real, optional, intent(in) :: region_mask(:,:,:)
        integer :: nregions, ipix, i, j, ilabel, iregion
        nregions = 1
        if( present(region_mask) ) nregions = 2
        allocate(target_counts(ncand,nregions), source=0)
        allocate(counts(ncand,nregions),        source=0)
        allocate(region_sizes(nregions),        source=0)
        allocate(hist_scales(nregions),         source=0.)
        allocate(pixel_regions(size(pix,2)),    source=1)
        allocate(region_priors(ncand,nregions), source=0.)
        do ipix = 1, size(pix,2)
            i = pix(1,ipix)
            j = pix(2,ipix)
            iregion = NU2D_HIST_REGION_FOREGROUND
            if( present(region_mask) )then
                if( region_mask(i,j,1) <= 0.5 ) iregion = NU2D_HIST_REGION_BACKGROUND
            endif
            pixel_regions(ipix) = iregion
            ilabel = int(labelmap(i,j,1))
            if( ilabel < 1 .or. ilabel > ncand ) cycle
            target_counts(ilabel,iregion) = target_counts(ilabel,iregion) + 1
            region_sizes(iregion) = region_sizes(iregion) + 1
        end do
        counts = target_counts
        if( present(region_mask) )then
            call init_nu2D_region_log_priors(target_counts, region_sizes, beta, region_priors)
            target_counts(:,NU2D_HIST_REGION_BACKGROUND) = 0
            target_counts(1,NU2D_HIST_REGION_BACKGROUND) = region_sizes(NU2D_HIST_REGION_BACKGROUND)
        endif
        do iregion = 1, nregions
            hist_scales(iregion) = NU2D_LABEL_HIST_BETA_FRAC * beta / real(max(1, region_sizes(iregion)))
        end do
    end subroutine init_nu2D_histogram_targets

    subroutine init_nu2D_potts_scales( potts_scales, region_mask )
        real, allocatable, intent(out) :: potts_scales(:)
        real, optional, intent(in) :: region_mask(:,:,:)
        integer :: nregions
        nregions = 1
        if( present(region_mask) ) nregions = 2
        allocate(potts_scales(nregions), source=1.)
        if( present(region_mask) )then
            potts_scales(NU2D_HIST_REGION_FOREGROUND) = NU2D_FOREGROUND_POTTS_BETA_SCALE
            potts_scales(NU2D_HIST_REGION_BACKGROUND) = NU2D_BACKGROUND_POTTS_BETA_SCALE
        endif
    end subroutine init_nu2D_potts_scales

    subroutine init_nu2D_region_log_priors( raw_counts, region_sizes, beta, region_priors )
        integer, intent(in)    :: raw_counts(:,:), region_sizes(:)
        real,    intent(in)    :: beta
        real,    intent(inout) :: region_priors(:,:)
        integer :: ilabel, ncand
        real    :: alpha, fg_denom, bg_denom, p_fg, p_bg, logodds, prior_scale
        ncand = size(raw_counts,1)
        if( size(raw_counts,2) < 2 ) return
        alpha = NU2D_REGION_PRIOR_PSEUDOCOUNT
        fg_denom = real(region_sizes(NU2D_HIST_REGION_FOREGROUND)) + alpha * real(ncand)
        bg_denom = real(region_sizes(NU2D_HIST_REGION_BACKGROUND)) + alpha * real(ncand)
        prior_scale = NU2D_REGION_LOGODDS_BETA_FRAC * beta
        do ilabel = 1, ncand
            p_fg = (real(raw_counts(ilabel,NU2D_HIST_REGION_FOREGROUND)) + alpha) / fg_denom
            p_bg = (real(raw_counts(ilabel,NU2D_HIST_REGION_BACKGROUND)) + alpha) / bg_denom
            logodds = log(max(p_fg,TINY) / max(p_bg,TINY))
            region_priors(ilabel,NU2D_HIST_REGION_FOREGROUND) = -prior_scale * logodds
            region_priors(ilabel,NU2D_HIST_REGION_BACKGROUND) =  prior_scale * logodds
        end do
    end subroutine init_nu2D_region_log_priors

    real function nu2D_histogram_delta( icand, cur_icand, counts, target_counts, iregion, hist_scale )
        integer, intent(in) :: icand, cur_icand, counts(:,:), target_counts(:,:), iregion
        real,    intent(in) :: hist_scale
        real :: cur_before, cand_before
        if( icand == cur_icand )then
            nu2D_histogram_delta = 0.
            return
        endif
        cur_before  = real(counts(cur_icand,iregion) - target_counts(cur_icand,iregion))
        cand_before = real(counts(icand,iregion)     - target_counts(icand,iregion))
        nu2D_histogram_delta = hist_scale * ((cur_before - 1.)**2 - cur_before**2 + &
            &(cand_before + 1.)**2 - cand_before**2)
    end function nu2D_histogram_delta

    real function estimate_nu2D_beta( dmats )
        real,    intent(in) :: dmats(:,:)
        integer :: ipix, icand, npix
        real :: best_e, second_e, cur_e, beta_sum
        estimate_nu2D_beta = 0.
        beta_sum = 0.
        npix = 0
        !$omp parallel do schedule(static) default(shared) private(ipix,icand,best_e,second_e,cur_e) &
        !$omp reduction(+:beta_sum,npix) proc_bind(close)
        do ipix = 1, size(dmats,1)
            best_e   = huge(best_e)
            second_e = huge(second_e)
            do icand = 1, size(dmats,2)
                cur_e = dmats(ipix,icand)
                if( cur_e < best_e )then
                    second_e = best_e
                    best_e   = cur_e
                else if( cur_e < second_e )then
                    second_e = cur_e
                endif
            end do
            if( second_e < huge(second_e) )then
                beta_sum = beta_sum + max(0., second_e - best_e)
                npix = npix + 1
            endif
        end do
        !$omp end parallel do
        if( npix > 0 ) estimate_nu2D_beta = NU2D_LABEL_SMOOTH_BETA_FRAC * beta_sum / real(npix)
    end function estimate_nu2D_beta

    real function nu2D_neighborhood_cost( icand, labelmap, neigh, neigh_w, nsz, coords )
        integer, intent(in) :: icand, neigh(3,NU2D_LABEL_SMOOTH_NNEIGH), nsz
        integer(kind=NU2D_LABEL_KIND), intent(in) :: labelmap(:,:,:)
        real,    intent(in) :: neigh_w(NU2D_LABEL_SMOOTH_NNEIGH), coords(:)
        integer :: ineigh, ni, nj, jlabel, degree
        real    :: wsum
        nu2D_neighborhood_cost = 0.
        degree = 0
        wsum = 0.
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            jlabel = int(labelmap(ni,nj,1))
            if( jlabel < 1 .or. jlabel > size(coords) ) cycle
            degree = degree + 1
            wsum = wsum + neigh_w(ineigh)
            nu2D_neighborhood_cost = nu2D_neighborhood_cost + &
                &neigh_w(ineigh) * nu2D_pair_cost(coords(icand), coords(jlabel))
        end do
        if( degree > 0 .and. wsum > TINY ) nu2D_neighborhood_cost = nu2D_neighborhood_cost / wsum
    end function nu2D_neighborhood_cost

    real function nu2D_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
        real :: adjacent_jump, excess_jump, label_jump
        if( abs(icoord - jcoord) <= TINY )then
            nu2D_pair_cost = 0.
        else
            label_jump = abs(icoord - jcoord)
            adjacent_jump = min(label_jump, real(NU2D_LABEL_SMOOTH_ADJACENT_STEPS))
            excess_jump = max(0., label_jump - real(NU2D_LABEL_SMOOTH_ADJACENT_STEPS))
            nu2D_pair_cost = NU2D_LABEL_SMOOTH_ADJACENT_FRAC * adjacent_jump + &
                &excess_jump + NU2D_LABEL_SMOOTH_QUAD_FRAC * excess_jump * excess_jump
        endif
    end function nu2D_pair_cost

    integer function nu2D_label_color( i, j )
        integer, intent(in) :: i, j
        nu2D_label_color = mod(i - 1, NU2D_LABEL_SMOOTH_RADIUS + 1) + &
            &(NU2D_LABEL_SMOOTH_RADIUS + 1) * mod(j - 1, NU2D_LABEL_SMOOTH_RADIUS + 1)
    end function nu2D_label_color

    subroutine collect_nu2D_neighbors( ldim, i, j, neigh, neigh_w, nsz )
        integer, intent(in)  :: ldim(3), i, j
        integer, intent(out) :: neigh(3,NU2D_LABEL_SMOOTH_NNEIGH), nsz
        real,    intent(out) :: neigh_w(NU2D_LABEL_SMOOTH_NNEIGH)
        integer :: di, dj, ni, nj, ring
        nsz = 0
        do dj = -NU2D_LABEL_SMOOTH_RADIUS, NU2D_LABEL_SMOOTH_RADIUS
            do di = -NU2D_LABEL_SMOOTH_RADIUS, NU2D_LABEL_SMOOTH_RADIUS
                if( di == 0 .and. dj == 0 ) cycle
                ring = max(abs(di), abs(dj))
                ni = i + di
                nj = j + dj
                if( ni < 1 .or. ni > ldim(1) ) cycle
                if( nj < 1 .or. nj > ldim(2) ) cycle
                nsz = nsz + 1
                neigh(:,nsz) = [ni,nj,1]
                if( ring == 1 )then
                    neigh_w(nsz) = NU2D_LABEL_SMOOTH_RING1_WEIGHT
                else
                    neigh_w(nsz) = NU2D_LABEL_SMOOTH_RING2_WEIGHT
                endif
            end do
        end do
    end subroutine collect_nu2D_neighbors

    logical function nu2D_is_better( e, best_e )
        real, intent(in) :: e, best_e
        real :: tol
        tol = NU2D_LABEL_SMOOTH_TIE_EPS * max(1., abs(best_e))
        nu2D_is_better = e < best_e - tol
    end function nu2D_is_better

end module simple_nu_filter2D
