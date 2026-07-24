module simple_motion_gain_analysis
use simple_core_module_api
use simple_image, only: image
implicit none
private
#include "simple_local_flags.inc"

public :: gain_flip_analyzer
public :: GAIN_FIRST_ANALYSIS_AT

real,    parameter :: GAIN_LP_RES_A              = 200.
integer, parameter :: GAIN_FIRST_ANALYSIS_AT     = 10
integer, parameter :: GAIN_ANALYSIS_STEP         = 5
integer, parameter :: GAIN_CONV_STREAK_REQ       = 2
real,    parameter :: GAIN_CONV_EPS_DELTA        = 1.e-2
real,    parameter :: GAIN_CONV_EPS_SMOOTH_CORR  = 3.e-3
real,    parameter :: GAIN_CONV_EPS_SMOOTH_DELTA = 3.e-3
real,    parameter :: GAIN_CONV_SMOOTH_ALPHA     = 0.4

type :: gain_flip_analyzer

    type(string) :: gainref_fname
    type(image)  :: gain_ref, gain_flip_x, gain_flip_y, gain_flip_xy
    type(image)  :: sum_img
    logical      :: initialized           = .false.
    logical      :: have_prev_analysis    = .false.
    logical      :: have_smooth_state     = .false.
    logical      :: converged             = .false.
    integer      :: best_idx              = 0
    integer      :: analysis_count        = 0
    integer      :: conv_streak           = 0
    integer      :: prev_best_idx         = 0
    integer      :: total_frames          = 0
    integer      :: total_movies          = 0
    real         :: smpd                  = 0.
    real         :: best_corr             = 0.
    real         :: second_corr           = 0.
    real         :: delta_r               = 0.
    real         :: corrs(4)              = 0.
    real         :: smooth_best_corr      = 0.
    real         :: smooth_delta_r        = 0.
    real         :: prev_smooth_best_corr = 0.
    real         :: prev_smooth_delta_r   = 0.

contains
    procedure :: new => analyzer_new
    procedure :: analyze_if_due
    procedure :: get_converged
    procedure :: kill => analyzer_kill
end type gain_flip_analyzer

contains

    subroutine analyzer_new(self, gainref_fname, smpd)
        class(gain_flip_analyzer), intent(inout) :: self
        class(string),             intent(in)    :: gainref_fname
        real,                      intent(in)    :: smpd
        integer :: ldim_gain(3), ifoo

        call self%kill()

        if( .not. file_exists(gainref_fname) )then
            THROW_HARD('Gain reference file does not exist: '//gainref_fname%to_char())
        endif

        self%gainref_fname = gainref_fname
        self%smpd          = smpd

        call find_ldim_nptcls(self%gainref_fname, ldim_gain, ifoo)
        ldim_gain(3) = 1

        call self%gain_ref%new(ldim_gain, self%smpd, wthreads=.false.)
        call self%gain_ref%read(self%gainref_fname)
        call self%gain_ref%bp(0., GAIN_LP_RES_A)

        call self%gain_flip_x%new(ldim_gain, self%smpd, wthreads=.false.)
        call self%gain_flip_x%copy(self%gain_ref)
        call self%gain_flip_x%flip('X')

        call self%gain_flip_y%new(ldim_gain, self%smpd, wthreads=.false.)
        call self%gain_flip_y%copy(self%gain_ref)
        call self%gain_flip_y%flip('Y')

        call self%gain_flip_xy%new(ldim_gain, self%smpd, wthreads=.false.)
        call self%gain_flip_xy%copy(self%gain_ref)
        call self%gain_flip_xy%flip('XY')

        call self%sum_img%new(ldim_gain, self%smpd, wthreads=.false.)
        call self%sum_img%zero()

        self%initialized = .true.

    end subroutine analyzer_new


    subroutine analyze_if_due(self, sum_part_img, part_frames, part_movies, ran_analysis)
        class(gain_flip_analyzer), intent(inout) :: self
        class(image),              intent(in)    :: sum_part_img
        integer,                   intent(in)    :: part_frames
        integer,                   intent(in)    :: part_movies
        logical,                   intent(out)   :: ran_analysis
        type(image) :: tmp_avg, tmp_gain
        integer :: iori
        integer :: ldim_sum(3), ldim_sum_part(3)
        logical :: idx_stable

        ran_analysis = .false.
        if( .not. self%initialized )then
            THROW_HARD('gain_flip_analyzer is not initialized; call new first')
        endif

        if( self%converged  ) return
        if( part_frames < 1 ) return

        ldim_sum      = self%sum_img%get_ldim()
        ldim_sum_part = sum_part_img%get_ldim()
        if( ldim_sum(1) /= ldim_sum_part(1) .or. ldim_sum(2) /= ldim_sum_part(2) )then
            THROW_HARD('Sum image dimensions differ from gain reference dimensions')
        endif

        call self%sum_img%add_workshare(sum_part_img)
        self%total_frames = self%total_frames + part_frames
        self%total_movies = self%total_movies + part_movies

        if( self%total_movies < GAIN_FIRST_ANALYSIS_AT ) return
        if( mod(self%total_movies - GAIN_FIRST_ANALYSIS_AT, GAIN_ANALYSIS_STEP) /= 0 ) return

        self%analysis_count = self%analysis_count + 1
        write(logfhandle,'(a)') '>>> ------------------------------------------------------------'
        write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> ANALYSIS', self%analysis_count, 'at movie count', self%total_movies

        call tmp_avg%new(ldim_sum, self%smpd, wthreads=.false.)
        call tmp_gain%new(ldim_sum, self%smpd, wthreads=.false.)

        call tmp_avg%copy(self%sum_img)
        call tmp_avg%div(real(self%total_frames))
        call tmp_avg%bp(0., GAIN_LP_RES_A)
        write(logfhandle,'(a,1x,i0)') '>>> TOTAL FRAMES AVERAGED:', self%total_frames
        write(logfhandle,'(a,1x,f8.2,a)') '>>> Applied low-pass filter to average:', GAIN_LP_RES_A, 'A'
        call tmp_avg%write(string('average_all_movies.mrc'))

        call tmp_gain%copy(self%gain_ref)
        self%corrs(1) = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(self%gain_flip_x)
        self%corrs(2) = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(self%gain_flip_y)
        self%corrs(3) = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(self%gain_flip_xy)
        self%corrs(4) = tmp_avg%corr(tmp_gain)

        write(logfhandle,'(a)') '>>> CORRELATION TABLE: average_all_movies vs gain variants'
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_unchanged :', self%corrs(1)
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_x    :', self%corrs(2)
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_y    :', self%corrs(3)
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_xy   :', self%corrs(4)

        self%best_idx = 1
        self%best_corr = self%corrs(1)
        if( self%corrs(2) < self%best_corr )then
            self%best_idx = 2
            self%best_corr = self%corrs(2)
        endif
        if( self%corrs(3) < self%best_corr )then
            self%best_idx = 3
            self%best_corr = self%corrs(3)
        endif
        if( self%corrs(4) < self%best_corr )then
            self%best_idx = 4
            self%best_corr = self%corrs(4)
        endif

        self%second_corr = huge(self%second_corr)
        do iori=1,4
            if( iori == self%best_idx ) cycle
            if( self%corrs(iori) < self%second_corr ) self%second_corr = self%corrs(iori)
        enddo
        self%delta_r = self%second_corr - self%best_corr
        write(logfhandle,'(a,1x,f12.6)') '>>> CONVERGENCE METRIC delta_r (2nd-best - best):', self%delta_r

        select case(self%best_idx)
        case(1)
            call self%gain_ref%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): unchanged, corr=', self%best_corr
        case(2)
            call self%gain_flip_x%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_x, corr=', self%best_corr
        case(3)
            call self%gain_flip_y%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_y, corr=', self%best_corr
        case(4)
            call self%gain_flip_xy%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_xy, corr=', self%best_corr
        end select

        if( .not. self%have_smooth_state )then
            self%smooth_best_corr = self%best_corr
            self%smooth_delta_r   = self%delta_r
            self%have_smooth_state = .true.
        else
            self%smooth_best_corr = GAIN_CONV_SMOOTH_ALPHA * self%best_corr + (1.0 - GAIN_CONV_SMOOTH_ALPHA) * self%smooth_best_corr
            self%smooth_delta_r   = GAIN_CONV_SMOOTH_ALPHA * self%delta_r   + (1.0 - GAIN_CONV_SMOOTH_ALPHA) * self%smooth_delta_r
        endif

        if( self%have_prev_analysis )then
            idx_stable = (self%best_idx == self%prev_best_idx) .or. (abs(self%delta_r) <= GAIN_CONV_EPS_DELTA)
            if( idx_stable .and. abs(self%smooth_best_corr - self%prev_smooth_best_corr) <= GAIN_CONV_EPS_SMOOTH_CORR .and. &
                abs(self%smooth_delta_r - self%prev_smooth_delta_r) <= GAIN_CONV_EPS_SMOOTH_DELTA )then
                self%conv_streak = self%conv_streak + 1
            else
                self%conv_streak = 0
            endif
        endif

        self%prev_best_idx = self%best_idx
        self%prev_smooth_best_corr = self%smooth_best_corr
        self%prev_smooth_delta_r = self%smooth_delta_r
        self%have_prev_analysis = .true.

        write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> CONVERGENCE streak:', self%conv_streak, 'required=', GAIN_CONV_STREAK_REQ
        write(logfhandle,'(a,1x,f12.6,1x,f12.6)') '>>> CONVERGENCE smooth best/delta:', self%smooth_best_corr, self%smooth_delta_r

        if( self%conv_streak >= GAIN_CONV_STREAK_REQ )then
            self%converged = .true.
            write(logfhandle,'(a,1x,i0)') '>>> CONVERGED at movie count:', part_movies
        endif

        call tmp_gain%kill()
        call tmp_avg%kill()

        ran_analysis = .true.
    end subroutine analyze_if_due


    logical function get_converged(self) result(is_converged)
        class(gain_flip_analyzer), intent(in) :: self
        is_converged = self%converged
    end function get_converged


    subroutine analyzer_kill(self)
        class(gain_flip_analyzer), intent(inout) :: self
        if( self%initialized )then
            call self%gain_ref%kill()
            call self%gain_flip_x%kill()
            call self%gain_flip_y%kill()
            call self%gain_flip_xy%kill()
            call self%sum_img%kill()
            self%initialized = .false.
        endif
        call self%gainref_fname%kill()
        self%best_idx              = 0
        self%best_corr             = 0.
        self%second_corr           = 0.
        self%delta_r               = 0.
        self%corrs                 = 0.
        self%analysis_count        = 0
        self%conv_streak           = 0
        self%prev_best_idx         = 0
        self%total_frames          = 0
        self%total_movies          = 0
        self%smooth_best_corr      = 0.
        self%smooth_delta_r        = 0.
        self%prev_smooth_best_corr = 0.
        self%prev_smooth_delta_r   = 0.
        self%have_prev_analysis    = .false.
        self%have_smooth_state     = .false.
        self%converged             = .false.
        self%initialized           = .false.
    end subroutine analyzer_kill

end module simple_motion_gain_analysis
