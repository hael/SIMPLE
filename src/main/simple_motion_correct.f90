! motion_correct does motion correction, dose-weighting and frame-weighting of DDD movies
module simple_motion_correct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ft_expanded,   only: ft_expanded, ft_exp_reset_tmp_pointers
use simple_image,         only: image
use simple_parameters,    only: params_glob
use simple_estimate_ssnr, only: acc_dose2filter
use simple_oris,          only: oris
implicit none

public :: motion_correct_movie, motion_correct_calc_sums, motion_correct_calc_sums_tomo
private
#include "simple_local_flags.inc"

interface motion_correct_calc_sums
    module procedure motion_correct_calc_sums_1
    module procedure motion_correct_calc_sums_2
end interface

type(ft_expanded), allocatable :: movie_frames_ftexp(:)             !< movie frames
type(ft_expanded), allocatable :: movie_frames_ftexp_sh(:)          !< shifted movie frames
type(ft_expanded), allocatable :: movie_sum_global_ftexp_threads(:) !< array of global movie sums for parallel refinement
type(image),       allocatable :: movie_frames_scaled(:)            !< scaled movie frames
type(image),       allocatable :: movie_frames_shifted(:)           !< shifted movie frames
type(image)                    :: movie_sum_global                  !< global movie sum for output
real, allocatable              :: corrs(:)                          !< per-frame correlations
real, allocatable              :: frameweights(:)                   !< array of frameweights
real, allocatable              :: frameweights_saved(:)             !< array of frameweights
real, allocatable              :: opt_shifts(:,:)                   !< optimal shifts identified
real, allocatable              :: opt_shifts_saved(:,:)             !< optimal shifts for local opt saved
real, allocatable              :: acc_doses(:)                      !< accumulated doses
complex,          allocatable  :: cmat(:,:,:), cmat_sum(:,:,:)      !< complex matrices for OpenMP reduction
integer                        :: nframes        = 0                !< number of frames
integer                        :: fixed_frame    = 0                !< fixed frame of reference (0,0)
integer                        :: ldim(3)        = [0,0,0]          !< logical dimension of frame
integer                        :: ldim_scaled(3) = [0,0,0]          !< shrunken logical dimension of frame
real                           :: hp             = 0.               !< high-pass limit
real                           :: lp             = 0.               !< low-pass limit
real                           :: resstep        = 0.               !< resolution step size (in angstrom)
real                           :: smpd           = 0.               !< sampling distance
real                           :: smpd_scaled    = 0.               !< sampling distance
real                           :: corr_saved     = 0.               !< opt corr for local opt saved
real                           :: kV             = 300.             !< acceleration voltage
real                           :: dose_rate      = 0.               !< dose rate
real                           :: nsig_here      = 6.0              !< nr of sigmas (for outlier removal)
logical                        :: do_dose_weight = .false.          !< dose weight or not
logical                        :: doscale        = .false.          !< scale or not
logical                        :: doprint        = .true.           !< print out correlations
logical                        :: existence      = .false.          !< to indicate existence

integer, parameter :: MITSREF    = 30 !< max nr iterations of refinement optimisation
real,    parameter :: SMALLSHIFT = 2. !< small initial shift to blur out fixed pattern noise
logical, parameter :: DEBUG_HERE = .true.

contains

    !> motion_correct DDD movie
    subroutine motion_correct_movie( movie_stack_fname, ctfvars, corr, shifts, err, gainref_fname, nsig )
        use simple_ftexp_shsrch, only: ftexp_shsrch
        character(len=*),           intent(in)    :: movie_stack_fname !< filename
        type(ctfparams),            intent(inout) :: ctfvars           !< CTF params
        real,                       intent(out)   :: corr              !< ave correlation per frame
        real,          allocatable, intent(out)   :: shifts(:,:)       !< the nframes shifts identified
        logical,                    intent(out)   :: err               !< error flag
        character(len=*), optional, intent(in)    :: gainref_fname     !< gain reference filename
        real,             optional, intent(in)    :: nsig              !< # sigmas (for outlier removal)
        type(ftexp_shsrch), allocatable :: ftexp_srch(:)
        real    :: ave, sdev, var, minw, maxw
        real    :: cxy(3), corr_prev, frac_improved, corrfrac
        integer :: iframe, iter, nimproved, updateres, i
        logical :: didsave, didupdateres, err_stat
        ! initialise
        nsig_here = 6.0
        if( present(nsig) ) nsig_here = nsig
        call motion_correct_init(movie_stack_fname, ctfvars, gainref_fname)
        err = .false.
        if( nframes < 2 )then
            err = .true.
            write(*,*) 'movie: ', trim(movie_stack_fname)
            THROW_WARN('nframes of movie < 2, aborting motion_correct')
            return
        endif
        ! make search objects for parallel execution
        allocate(ftexp_srch(nframes))
        do iframe=1,nframes
            call ftexp_srch(iframe)%new(movie_sum_global_ftexp_threads(iframe), movie_frames_ftexp(iframe),&
            &params_glob%scale * params_glob%trs, motion_correct_ftol = params_glob%motion_correctftol,&
            &motion_correct_gtol = params_glob%motion_correctgtol)
            ! initialise with small random shifts (to average out dead/hot pixels)
            opt_shifts(iframe,1) = ran3()*2.*SMALLSHIFT-SMALLSHIFT
            opt_shifts(iframe,2) = ran3()*2.*SMALLSHIFT-SMALLSHIFT
        end do
        ! generate movie sum for refinement
        call shift_wsum_calc_corrs(opt_shifts, 1)
        ! calc avg corr to weighted avg
        corr = sum(corrs)/real(nframes)
        if( doprint ) write(*,'(a)') '>>> WEIGHTED AVERAGE-BASED REFINEMENT'
        iter       = 0
        corr_saved = -1.
        didsave    = .false.
        updateres  = 0
        do i=1,MITSREF
            iter = iter+1
            nimproved = 0
            !$omp parallel do schedule(static) default(shared) private(iframe,cxy) proc_bind(close) reduction(+:nimproved)
            do iframe=1,nframes
                ! subtract the movie frame being aligned to reduce bias
                call movie_sum_global_ftexp_threads(iframe)%subtr(movie_frames_ftexp_sh(iframe), w=frameweights(iframe))
                ! set pointers
                call  ftexp_srch(iframe)%set_ptrs(movie_sum_global_ftexp_threads(iframe), movie_frames_ftexp(iframe))
                ! optimise shifts
                cxy = ftexp_srch(iframe)%minimize(corrs(iframe), opt_shifts(iframe,:))
                ! reset pointers used for calculation bookkeeping in ft_expanded
                call ft_exp_reset_tmp_pointers ! done for this thread only
                ! count # improved
                if( cxy(1) > corrs(iframe) ) nimproved = nimproved + 1
                ! update parameter arrays
                opt_shifts(iframe,:) = cxy(2:3)
                corrs(iframe)        = cxy(1)
                ! no need to add the frame back to the weighted sum since the sum will be updated after the loop
                ! (see below)
            end do
            !$omp end parallel do
            frac_improved = real(nimproved) / real(nframes) * 100.
            if( doprint ) write(*,'(a,1x,f4.0)') 'This % of frames improved their alignment: ', frac_improved
            call shift_wsum_calc_corrs( opt_shifts, 2 )
            corr_prev = corr
            corr      = sum(corrs) / real(nframes)
            if( corr >= corr_saved )then ! save the local optimum
                corr_saved         = corr
                opt_shifts_saved   = opt_shifts
                frameweights_saved = frameweights
                didsave = .true.
            endif
            corrfrac = corr_prev / corr
            didupdateres = .false.
            select case(updateres)
                case(0)
                    call update_res( 0.96, 70., updateres )
                case(1)
                    call update_res( 0.97, 60., updateres )
                case(2)
                    call update_res( 0.98, 50., updateres )
                case DEFAULT
                    ! nothing to do
            end select
            if( updateres > 2 .and. .not. didupdateres )then ! at least one iteration with new lim
                if( nimproved == 0 .and. i > 2 )  exit
                if( i > 10 .and. corrfrac > 0.9999 ) exit
            endif
        end do
        ! put the best local optimum back
        corr         = corr_saved
        opt_shifts   = opt_shifts_saved
        frameweights = frameweights_saved
        call shift_frames(opt_shifts)
        ! output shifts
        if( allocated(shifts) ) deallocate(shifts)
        allocate(shifts(nframes,2), source=opt_shifts)
        ! print
        if( corr < 0. )then
            if( doprint ) write(*,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
            if( doprint ) THROW_WARN('OPTIMAL CORRELATION < 0.0')
        endif
        call moment(frameweights, ave, sdev, var, err_stat)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT     :', ave
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS    :', sdev
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> MIN WEIGHT         :', minw
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> MAX WEIGHT         :', maxw
        ! report the sampling distance of the possibly scaled movies
        ctfvars%smpd = smpd_scaled
        ! destruct
        do iframe=1,nframes
            call ftexp_srch(iframe)%kill
        end do
        deallocate(ftexp_srch)

    contains

            subroutine update_res( thres_corrfrac, thres_frac_improved, which_update )
                real,    intent(in) :: thres_corrfrac, thres_frac_improved
                integer, intent(in) :: which_update
                if( corrfrac > thres_corrfrac .and. frac_improved <= thres_frac_improved&
                .and. updateres == which_update )then
                    lp = lp - resstep
                    if( doprint )  write(*,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', lp
                    ! need to re-make the ftexps
                    do iframe=1,nframes
                        call movie_frames_ftexp(iframe)%new(movie_frames_scaled(iframe), hp, lp, .true.)
                        call movie_frames_ftexp_sh(iframe)%new(movie_frames_scaled(iframe), hp, lp, .false.)
                        call movie_sum_global_ftexp_threads(iframe)%new(movie_frames_scaled(iframe), hp, lp, .false.)
                    end do
                    call shift_wsum_calc_corrs(opt_shifts,2)
                    ! need to indicate that we updated resolution limit
                    updateres  = updateres + 1
                    ! need to destroy all previous knowledge about correlations
                    corr       = sum(corrs)/real(nframes)
                    corr_prev  = corr
                    corr_saved = corr
                    ! indicate that reslim was updated
                    didupdateres = .true.
                endif
            end subroutine update_res

    end subroutine motion_correct_movie

    subroutine motion_correct_init( movie_stack_fname, ctfvars, gainref_fname )
        character(len=*),           intent(in) :: movie_stack_fname  !< input filename of stack
        type(ctfparams),            intent(in) :: ctfvars            !< CTF parameters
        character(len=*), optional, intent(in) :: gainref_fname      !< gain reference filename
        type(image), allocatable :: movie_frames(:)
        real,        allocatable :: rmat(:,:,:), rmat_pad(:,:), win(:,:), rmat_sum(:,:,:)
        logical,     allocatable :: outliers(:,:)
        type(image)        :: tmpmovsum, gainref
        real               :: moldiam, dimo4, time_per_frame, current_time
        integer            :: iframe, ncured, deadhot(2), i, j, winsz, ldim_scaled_tmp(3), shp(3), sz
        integer, parameter :: HWINSZ = 6
        logical            :: l_gain, do_alloc
        ! get number of frames & dim from stack
        call find_ldim_nptcls(movie_stack_fname, ldim, nframes)
        if( nframes < 2 ) return
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        if( params_glob%scale < 0.99 )then
            ldim_scaled(1) = round2even(real(ldim(1))*params_glob%scale)
            ldim_scaled(2) = round2even(real(ldim(2))*params_glob%scale)
            ldim_scaled(3) = 1
            doscale        = .true.
        else
            ldim_scaled = ldim
            doscale     = .false.
        endif
        ! set sampling distance
        smpd        = ctfvars%smpd
        smpd_scaled = ctfvars%smpd/params_glob%scale
        ! set fixed frame (all others are shifted by reference to this at 0,0)
        fixed_frame = nint(real(nframes)/2.)
        ! set reslims
        dimo4   = (real(minval(ldim_scaled(1:2)))*smpd_scaled)/4.
        moldiam = 0.7*real(params_glob%box)*smpd_scaled
        hp      = min(dimo4,2000.)
        lp      = params_glob%lpstart
        resstep = (params_glob%lpstart-params_glob%lpstop)/3.
        ! allocate
        do_alloc = .true.
        ! call motion_correct_kill
        if( allocated(movie_frames_scaled) )then
            sz = size(movie_frames_scaled)
            ldim_scaled_tmp = movie_frames_scaled(1)%get_ldim()
            if( sz == nframes .and. all(ldim_scaled == ldim_scaled_tmp) )then
                do_alloc = .false.
            else
                call motion_correct_kill
            endif
        endif
        if( do_alloc )then
            allocate( movie_frames_ftexp(nframes), movie_frames_scaled(nframes),&
                movie_frames_ftexp_sh(nframes), corrs(nframes), opt_shifts(nframes,2),&
                opt_shifts_saved(nframes,2), frameweights(nframes),&
                frameweights_saved(nframes), movie_sum_global_ftexp_threads(nframes),&
                movie_frames_shifted(nframes), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init; simple_motion_correct')
            do iframe=1,nframes
                call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
                call movie_frames_shifted(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
            end do
            ! retrieve array shape
            shp = movie_frames_scaled(1)%get_array_shape()
            ! create image objects for the shifted frames + cmats
            allocate(cmat(shp(1),shp(2),shp(3)), cmat_sum(shp(1),shp(2),shp(3)))
        endif
        ! init
        corrs              = 0.
        opt_shifts         = 0.
        opt_shifts_saved   = 0.
        frameweights       = 1./real(nframes)
        frameweights_saved = frameweights
        ! local allocations
        winsz = 2 * HWINSZ + 1
        allocate(rmat(ldim(1),ldim(2),1), rmat_sum(ldim(1),ldim(2),1), rmat_pad(1-HWINSZ:ldim(1)+HWINSZ,&
        1-HWINSZ:ldim(2)+HWINSZ), win(winsz,winsz), movie_frames(nframes), stat=alloc_stat)
        if( alloc_stat.ne.0 )call allocchk('motion_correct_init; simple_motion_correct, rmat etc.')
        rmat     = 0.
        rmat_sum = 0.
        rmat_pad = 0.
        win      = 0.
        ! gain reference
        l_gain = .false.
        if( present(gainref_fname) )then
            if( file_exists(gainref_fname) )then
                l_gain = .true.
                call gainref%new(ldim, smpd)
                call gainref%read(gainref_fname)
            else
                THROW_HARD('gain reference: '//trim(gainref_fname)//' not found; motion_correct_init')
            endif
        endif
        ! allocate & read frames
        do iframe=1,nframes
            call movie_frames(iframe)%new(ldim, smpd, wthreads=.false.)
            call movie_frames(iframe)%read(movie_stack_fname, iframe)
        end do
        ! calculate image sum and identify outliers
        if( doprint ) write(*,'(a)') '>>> REMOVING DEAD/HOT PIXELS & FOURIER TRANSFORMING FRAMES'
        ! gain correction, generation of temporary movie sum and outlier detection
        rmat_sum = 0.
        !$omp parallel do schedule(static) default(shared) private(iframe,rmat) proc_bind(close) reduction(+:rmat_sum)
        do iframe=1,nframes
            if( l_gain ) call movie_frames(iframe)%mul(gainref) ! gain correction
            call movie_frames(iframe)%get_rmat_sub(rmat)
            rmat_sum = rmat_sum + rmat
        end do
        !$omp end parallel do
        if( l_gain ) call gainref%kill
        call tmpmovsum%new(ldim, smpd)
        call tmpmovsum%set_rmat(rmat_sum)
        call tmpmovsum%cure_outliers(ncured, nsig_here, deadhot, outliers)
        call tmpmovsum%kill
        write(*,'(a,1x,i7)') '>>> # DEAD PIXELS:', deadhot(1)
        write(*,'(a,1x,i7)') '>>> # HOT  PIXELS:', deadhot(2)
        if( any(outliers) )then ! remove the outliers & do the rest
            !$omp parallel do schedule(static) default(shared) private(iframe,rmat,rmat_pad,i,j,win) proc_bind(close)
            do iframe=1,nframes
                call movie_frames(iframe)%get_rmat_sub(rmat)
                rmat_pad = median( reshape(rmat(:,:,1),(/(ldim(1)*ldim(2))/)) )
                rmat_pad(1:ldim(1),1:ldim(2)) = rmat(1:ldim(1),1:ldim(2),1)
                do i=1,ldim(1)
                    do j=1,ldim(2)
                        if( outliers(i,j) )then
                            win = rmat_pad( i-HWINSZ:i+HWINSZ, j-HWINSZ:j+HWINSZ )
                            rmat(i,j,1) = median( reshape(win,(/winsz**2/)) )
                        endif
                    enddo
                enddo
                call movie_frames(iframe)%set_rmat(rmat)
                call movie_frames(iframe)%fft()
                call movie_frames(iframe)%clip(movie_frames_scaled(iframe))
            enddo
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
            do iframe=1,nframes
                call movie_frames(iframe)%fft()
                call movie_frames(iframe)%clip(movie_frames_scaled(iframe))
            end do
            !$omp end parallel do
        endif
        ! create the expanded data structures
        do iframe=1,nframes
            call movie_frames_ftexp(iframe)%new(movie_frames_scaled(iframe), hp, lp, .true.)
            call movie_frames_ftexp_sh(iframe)%new(movie_frames_scaled(iframe), hp, lp, .false.)
            call movie_sum_global_ftexp_threads(iframe)%new(movie_frames_scaled(iframe), hp, lp, .false.)
        end do
        ! check if we are doing dose weighting
        if( params_glob%l_dose_weight )then
            do_dose_weight = .true.
            if( allocated(acc_doses) ) deallocate(acc_doses)
            allocate( acc_doses(nframes), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init; simple_motion_correct, acc_doses')
            kV = ctfvars%kv
            time_per_frame = params_glob%exp_time/real(nframes)  ! unit: s
            dose_rate      = params_glob%dose_rate
            do iframe=1,nframes
                current_time      = real(iframe)*time_per_frame ! unit: s
                acc_doses(iframe) = dose_rate*current_time      ! unit: e/A2/s * s = e/A2
            end do
        endif
        ! local deallocation
        do iframe=1,nframes
            call movie_frames(iframe)%kill
        end do
        deallocate(rmat, rmat_sum, rmat_pad, win, outliers, movie_frames)
        existence = .true.
    end subroutine motion_correct_init

    subroutine shift_frames( shifts )
        real, intent(in) :: shifts(nframes,2)
        integer :: iframe
        real    :: shvec(3)
        !$omp parallel do schedule(static) default(shared) private(iframe,shvec) proc_bind(close)
        do iframe=1,nframes
            shvec(1) = -shifts(iframe,1)
            shvec(2) = -shifts(iframe,2)
            shvec(3) = 0.0
            call movie_frames_ftexp(iframe)%shift(shvec, movie_frames_ftexp_sh(iframe))
        end do
        !$omp end parallel do
    end subroutine shift_frames

    subroutine corrmat2weights
        integer :: iframe, jframe
        real    :: corrmat(nframes,nframes)
        corrmat = 1. ! diagonal elements are 1
        corrs   = 0.
        !$omp parallel default(shared) private(iframe,jframe) proc_bind(close)
        !$omp do schedule(guided)
        do iframe=1,nframes-1
            do jframe=iframe+1,nframes
                corrmat(iframe,jframe) = movie_frames_ftexp_sh(iframe)%corr(movie_frames_ftexp_sh(jframe))
                corrmat(jframe,iframe) = corrmat(iframe,jframe)
            end do
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do iframe=1,nframes
            do jframe=1,nframes
                if( jframe == iframe ) cycle
                corrs(iframe) = corrs(iframe)+corrmat(iframe,jframe)
            end do
            corrs(iframe) = corrs(iframe)/real(nframes-1)
        end do
        !$omp end do nowait
        !$omp end parallel
        frameweights = corrs2weights(corrs)
    end subroutine corrmat2weights

    subroutine shift_wsum_calc_corrs( shifts, imode )
        real,     intent(inout) :: shifts(nframes,2)
        integer,  intent(in)    :: imode
        complex, allocatable    :: cmat_sum(:,:,:), cmat(:,:,:)
        integer :: iframe, flims(3,2)
        real    :: shvec(3), w, xsh, ysh
        ! CENTER SHIFTS
        xsh = -shifts(fixed_frame,1)
        ysh = -shifts(fixed_frame,2)
        do iframe=1,nframes
            shifts(iframe,1) = shifts(iframe,1) + xsh
            shifts(iframe,2) = shifts(iframe,2) + ysh
            if( abs(shifts(iframe,1)) < 1e-6 ) shifts(iframe,1) = 0.
            if( abs(shifts(iframe,2)) < 1e-6 ) shifts(iframe,2) = 0.
        end do
        if( imode == 1 )then
            call corrmat2weights ! initialisation
        else if( imode == 2 )then
            frameweights = corrs2weights(corrs) ! update
        endif
        ! allocate matrices for reduction
        flims = movie_sum_global_ftexp_threads(1)%get_flims()
        allocate(cmat(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)),&
             cmat_sum(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)), source=cmplx(0.,0.))
        ! FIRST LOOP TO OBTAIN WEIGHTED SUM
        !$omp parallel default(shared) private(iframe,shvec,cmat) proc_bind(close)
        !$omp do schedule(static) reduction(+:cmat_sum)
        do iframe=1,nframes
            shvec(1) = -shifts(iframe,1)
            shvec(2) = -shifts(iframe,2)
            shvec(3) = 0.0
            call movie_frames_ftexp(iframe)%shift(shvec, movie_frames_ftexp_sh(iframe))
            call movie_frames_ftexp_sh(iframe)%get_cmat(cmat)
            cmat_sum = cmat_sum + cmat * frameweights(iframe)
        end do
        !$omp end do nowait
        ! SECOND LOOP TO UPDATE movie_sum_global_ftexp_threads AND CALCULATE CORRS
        !$omp do schedule(static)
        do iframe=1,nframes
            ! update array of sums (for future parallel exec)
            call movie_sum_global_ftexp_threads(iframe)%set_cmat(cmat_sum)
            ! subtract the movie frame being correlated to reduce bias
            call movie_sum_global_ftexp_threads(iframe)%subtr(movie_frames_ftexp_sh(iframe), w=frameweights(iframe))
            ! calc corr
            corrs(iframe) = movie_sum_global_ftexp_threads(iframe)%corr(movie_frames_ftexp_sh(iframe))
            ! add the subtracted movie frame back to the weighted sum
            call movie_sum_global_ftexp_threads(iframe)%add(movie_frames_ftexp_sh(iframe), w=frameweights(iframe))
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine shift_wsum_calc_corrs

    subroutine motion_correct_calc_sums_1( movie_sum, movie_sum_corrected, movie_sum_ctf )
        type(image), intent(out) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        ! calculate the sum for CTF estimation
        call sum_movie_frames(opt_shifts)
        movie_sum_ctf = movie_sum_global
        call movie_sum_ctf%ifft()
        ! re-calculate the weighted sum
        call wsum_movie_frames(opt_shifts)
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%ifft()
        ! generate straight integrated movie frame for comparison
        call sum_movie_frames
        movie_sum = movie_sum_global
        call movie_sum%ifft()
    end subroutine motion_correct_calc_sums_1

    subroutine motion_correct_calc_sums_2( movie_sum_corrected, fromto )
        type(image), intent(out) :: movie_sum_corrected
        integer,     intent(in)  :: fromto(2)
        logical :: l_tmp
        ! re-calculate the weighted sum with dose_weighting turned off
        l_tmp = do_dose_weight
        do_dose_weight = .false.
        call wsum_movie_frames(opt_shifts, fromto)
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%ifft()
        do_dose_weight = l_tmp
    end subroutine motion_correct_calc_sums_2

    subroutine motion_correct_calc_sums_tomo( frame_counter, time_per_frame, movie_sum, movie_sum_corrected, movie_sum_ctf )
        integer,     intent(inout) :: frame_counter  !< frame counter
        real,        intent(in)    :: time_per_frame !< time resolution
        type(image), intent(out)   :: movie_sum, movie_sum_corrected, movie_sum_ctf
        ! calculate the sum for CTF estimation
        call sum_movie_frames(opt_shifts)
        movie_sum_ctf = movie_sum_global
        call movie_sum_ctf%ifft()
        ! re-calculate the weighted sum
        call wsum_movie_frames_tomo(opt_shifts, frame_counter, time_per_frame)
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%ifft()
        ! generate straight integrated movie frame for comparison
        call sum_movie_frames
        movie_sum = movie_sum_global
        call movie_sum%ifft()
    end subroutine motion_correct_calc_sums_tomo

    subroutine sum_movie_frames( shifts )
        real, intent(in), optional :: shifts(nframes,2)
        integer     :: iframe
        real        :: w
        logical     :: doshift
        type(image) :: frame_tmp
        doshift = present(shifts)
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
        call frame_tmp%new(ldim_scaled, smpd_scaled)
        w = 1./real(nframes)
        cmat_sum = cmplx(0.,0.)
        !$omp parallel do default(shared) private(iframe,cmat) proc_bind(close) schedule(static) reduction(+:cmat_sum)
        do iframe=1,nframes
            if( doshift )then
                call movie_frames_scaled(iframe)%shift2Dserial([-shifts(iframe,1),-shifts(iframe,2)], movie_frames_shifted(iframe))
                call movie_frames_shifted(iframe)%get_cmat_sub(cmat)
            else
                call movie_frames_scaled(iframe)%get_cmat_sub(cmat)
            endif
            cmat_sum = cmat_sum + cmat * w
        end do
        !$omp end parallel do
        call movie_sum_global%set_cmat(cmat_sum)
    end subroutine sum_movie_frames

    subroutine wsum_movie_frames( shifts, fromto )
        real,              intent(in) :: shifts(nframes,2)
        integer, optional, intent(in) :: fromto(2)
        real, allocatable :: filtarr(:,:)
        integer :: iframe, ffromto(2), sz
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
        if( do_dose_weight )then
            sz = movie_frames_scaled(1)%get_filtsz()
            allocate(filtarr(nframes,sz))
            do iframe=1,nframes
                filtarr(iframe,:) = acc_dose2filter(movie_frames_scaled(1), acc_doses(iframe), kV)
            end do
        endif
        ffromto(1) = 1
        ffromto(2) = nframes
        if( present(fromto) ) ffromto = fromto
        cmat_sum = cmplx(0.,0.)
        !$omp parallel do default(shared) private(iframe,cmat) proc_bind(close) schedule(static) reduction(+:cmat_sum)
        do iframe=ffromto(1),ffromto(2)
            if( frameweights(iframe) > 0. )then
                call movie_frames_scaled(iframe)%shift2Dserial([-shifts(iframe,1),-shifts(iframe,2)], movie_frames_shifted(iframe))
                if( do_dose_weight ) call movie_frames_shifted(iframe)%apply_filter_serial(filtarr(iframe,:))
                call movie_frames_shifted(iframe)%get_cmat_sub(cmat)
                cmat_sum = cmat_sum + cmat * frameweights(iframe)
            endif
        end do
        !$omp end parallel do
        call movie_sum_global%set_cmat(cmat_sum)
        if( allocated(filtarr) ) deallocate(filtarr)
    end subroutine wsum_movie_frames

    subroutine wsum_movie_frames_tomo( shifts, frame_counter, time_per_frame )
        real,    intent(in)    :: shifts(nframes,2)
        integer, intent(inout) :: frame_counter
        real,    intent(in)    :: time_per_frame
        real, allocatable :: filtarr(:,:)
        type(image)       :: frame_tmp
        integer           :: iframe, sz
        real              :: current_time, acc_dose
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
        sz = movie_frames_scaled(1)%get_filtsz()
        allocate(filtarr(nframes,sz))
        do iframe=1,nframes
            frame_counter     = frame_counter + 1
            current_time      = real(frame_counter)*time_per_frame ! unit: seconds
            acc_dose          = dose_rate*current_time             ! unit e/A2
            filtarr(iframe,:) = acc_dose2filter(movie_frames_scaled(iframe), acc_dose, kV)
        end do
        !$omp parallel do default(shared) private(iframe,cmat) proc_bind(close) schedule(static) reduction(+:cmat_sum)
        do iframe=1,nframes
            if( frameweights(iframe) > 0. )then
                call movie_frames_scaled(iframe)%shift2Dserial([-shifts(iframe,1),-shifts(iframe,2)], movie_frames_shifted(iframe))
                call movie_frames_shifted(iframe)%apply_filter_serial(filtarr(iframe,:))
                call movie_frames_shifted(iframe)%get_cmat_sub(cmat)
                cmat_sum = cmat_sum + cmat * frameweights(iframe)
            endif
        end do
        !$omp end parallel do
        call movie_sum_global%set_cmat(cmat_sum)
        if( allocated(filtarr) ) deallocate(filtarr)
    end subroutine wsum_movie_frames_tomo

    subroutine motion_correct_kill
        integer :: iframe
        if( existence )then
            do iframe=1,nframes
                call movie_frames_ftexp(iframe)%kill
                call movie_frames_ftexp_sh(iframe)%kill
                call movie_frames_scaled(iframe)%kill
                call movie_frames_shifted(iframe)%kill
                call movie_sum_global_ftexp_threads(iframe)%kill
            end do
            call movie_sum_global%kill
            deallocate( movie_frames_ftexp, movie_frames_ftexp_sh, movie_frames_scaled,&
            frameweights, frameweights_saved, corrs, opt_shifts, opt_shifts_saved,&
            movie_sum_global_ftexp_threads, movie_frames_shifted, cmat, cmat_sum)
            if( allocated(acc_doses) ) deallocate(acc_doses)
            existence = .false.
        endif
    end subroutine motion_correct_kill

end module simple_motion_correct
