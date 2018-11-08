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

type(ft_expanded), allocatable  :: movie_frames_ftexp(:)             !< movie frames
type(ft_expanded), allocatable  :: movie_frames_ftexp_sh(:)          !< shifted movie frames
type(ft_expanded), allocatable  :: movie_sum_global_ftexp_threads(:) !< array of global movie sums for parallel refinement
type(image),       allocatable  :: movie_frames(:)                   !< temporary frame
type(image),       allocatable  :: movie_frames_scaled(:)            !< scaled movie frames
type(image)                     :: movie_sum_global                  !< global movie sum for output
real, allocatable               :: corrs(:)                          !< per-frame correlations
real, allocatable               :: frameweights(:)                   !< array of frameweights
real, allocatable               :: frameweights_saved(:)             !< array of frameweights
real, allocatable               :: opt_shifts(:,:)                   !< optimal shifts identified
real, allocatable               :: opt_shifts_saved(:,:)             !< optimal shifts for local opt saved
real, allocatable               :: acc_doses(:)                      !< accumulated doses
integer                         :: nframes        = 0                !< number of frames
integer                         :: fixed_frame    = 0                !< fixed frame of reference (0,0)
integer                         :: ldim(3)        = [0,0,0]          !< logical dimension of frame
integer                         :: ldim_scaled(3) = [0,0,0]          !< shrunken logical dimension of frame
real                            :: hp             = 0.               !< high-pass limit
real                            :: lp             = 0.               !< low-pass limit
real                            :: resstep        = 0.               !< resolution step size (in angstrom)
real                            :: smpd           = 0.               !< sampling distance
real                            :: smpd_scaled    = 0.               !< sampling distance
real                            :: corr_saved     = 0.               !< opt corr for local opt saved
real                            :: kV             = 300.             !< acceleration voltage
real                            :: dose_rate      = 0.               !< dose rate
real                            :: nsig_here      = 6.0              !< nr of sigmas (for outlier removal)
logical                         :: do_dose_weight = .false.          !< dose weight or not
logical                         :: doscale        = .false.          !< scale or not
logical                         :: doprint        = .true.           !< print out correlations
logical                         :: existence      = .false.          !< to indicate existence

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
        real    :: ave, sdev, var, minw, maxw, cxy(3), corr_prev, frac_improved, corrfrac
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
        type(ctfparams),            intent(in) :: ctfvars
        character(len=*), optional, intent(in) :: gainref_fname     !< gain reference filename
        type(image)          :: tmpmovsum, gainref
        real                 :: moldiam, dimo4, time_per_frame, current_time
        integer              :: iframe, ncured, deadhot(2), i, j, winsz, ldim_tmp(3), ldim_scaled_tmp(3)
        integer, parameter   :: HWINSZ = 6
        real,    allocatable :: rmat(:,:,:), rmat_pad(:,:), win(:,:), rmat_sum(:,:,:)
        logical, allocatable :: outliers(:,:)
        logical              :: l_gain, do_alloc
        if( DEBUG_HERE ) print *, '(DEBUG) ALLOCATING AND READING'
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
        if( allocated(movie_frames_ftexp) )then
            if( size(movie_frames_ftexp) /= nframes )then
                call motion_correct_kill
            else
                ldim_tmp = movie_frames(1)%get_ldim()
                ldim_scaled_tmp = movie_frames_scaled(1)%get_ldim()
                if( all(ldim == ldim_tmp) .and. all(ldim_scaled == ldim_scaled_tmp) )then
                    do_alloc = .false.
                else
                    call motion_correct_kill
                endif
            endif
        endif
        if( do_alloc )then
            allocate( movie_frames_ftexp(nframes), movie_frames_scaled(nframes),&
                movie_frames_ftexp_sh(nframes), corrs(nframes), opt_shifts(nframes,2),&
                opt_shifts_saved(nframes,2), frameweights(nframes), movie_frames(nframes),&
                frameweights_saved(nframes), movie_sum_global_ftexp_threads(nframes), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init; simple_motion_correct')
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
        1-HWINSZ:ldim(2)+HWINSZ), win(winsz,winsz), source=0., stat=alloc_stat)
        if( alloc_stat.ne.0 )call allocchk('motion_correct_init; simple_motion_correct, rmat etc.')
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
        ! create objects (if needed)
        if( do_alloc )then
            do iframe=1,nframes
                call movie_frames(iframe)%new(ldim, smpd, wthreads=.false.)
                call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
            end do
        endif
        ! read frames
        do iframe=1,nframes
            call movie_frames(iframe)%read(movie_stack_fname, iframe)
        end do
        if( DEBUG_HERE ) print *, '(DEBUG) DONE, ALLOCATING AND READING'
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
        deallocate(rmat, rmat_sum, rmat_pad, win, outliers)
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

    subroutine motion_correct_calc_sums( movsum, movsum_ctf, movsum_corrected, movsum_fromto, fromto, filtarr_in )
        class(image),           intent(out) :: movsum, movsum_ctf, movsum_corrected
        class(image), optional, intent(out) :: movsum_fromto
        integer,      optional, intent(in)  :: fromto(2)
        real,         optional, intent(in)  :: filtarr_in(:,:)
        type(image), allocatable :: frames_shifted(:)
        complex,     allocatable :: cmat_sum_noshift(:,:,:), cmat_sum(:,:,:), cmat(:,:,:)
        complex,     allocatable :: cmat_shifted(:,:,:), cmat_wsum(:,:,:), cmat_wsum_fromto(:,:,:)
        real,        allocatable :: filtarr(:,:)
        integer :: iframe, sz, shp(3)
        real    :: w_scalar
        logical :: do_fromto
        do_fromto = .false.
        if( present(movsum_fromto) )then
            if( .not. present(fromto) )then
                THROW_HARD('need fromto optional input in conjunction with movie_wsum_fromto optional input; sum_movie_frames')
            endif
            do_fromto = .true.
        endif
        ! create image objects for the shifted frames
        allocate(frames_shifted(nframes))
        do iframe=1,nframes
            call frames_shifted(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
            call frames_shifted(iframe)%set_ft(.true.)
        end do
        ! create output images
        call movsum%new(ldim_scaled,           smpd_scaled, wthreads=.false.)
        call movsum_ctf%new(ldim_scaled,       smpd_scaled, wthreads=.false.)
        call movsum_corrected%new(ldim_scaled, smpd_scaled, wthreads=.false.)
        ! retrieve array shape
        shp = frames_shifted(1)%get_array_shape()
        ! create filter array for dose-weighting
        if( do_dose_weight )then
            sz = frames_shifted(1)%get_filtsz()
            allocate(filtarr(nframes,sz))
            if( present(filtarr_in) )then
                filtarr = filtarr_in
            else
                do iframe=1,nframes
                    filtarr(iframe,:) = acc_dose2filter(frames_shifted(1), acc_doses(iframe), kV)
                end do
            endif
        endif
        ! set scalar weight
        w_scalar = 1./real(nframes)
        ! allocate matrices
        allocate(cmat(shp(1),shp(2),shp(3)),             cmat_shifted(shp(1),shp(2),shp(3)), &
                &cmat_sum_noshift(shp(1),shp(2),shp(3)), cmat_sum(shp(1),shp(2),shp(3)),     &
                &cmat_wsum(shp(1),shp(2),shp(3)),        cmat_wsum_fromto(shp(1),shp(2),shp(3)), source=cmplx(0.,0.))
        !$omp parallel default(shared) private(iframe,cmat,cmat_shifted) proc_bind(close)
        !$omp do schedule(static) reduction(+:cmat_sum_noshift,cmat_sum,cmat_wsum,cmat_wsum_fromto)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%get_cmat_sub(cmat)
            call movie_frames_scaled(iframe)%shift2Dserial(-opt_shifts(iframe,:), frames_shifted(iframe))
            call frames_shifted(iframe)%get_cmat_sub(cmat_shifted)
            cmat_sum_noshift = cmat_sum_noshift + w_scalar * cmat         ! movie_sum
            cmat_sum         = cmat_sum         + w_scalar * cmat_shifted ! movie_sum_ctf
            if( do_dose_weight )then
                call frames_shifted(iframe)%apply_filter_serial(filtarr(iframe,:))
                call frames_shifted(iframe)%get_cmat_sub(cmat_shifted)
            endif
            cmat_wsum = cmat_wsum + frameweights(iframe) * cmat_shifted   ! movie_sum_corrected
            if( do_fromto )then
                if( iframe >= fromto(1) .and. iframe <= fromto(2) )then
                    cmat_wsum_fromto = cmat_wsum_fromto + frameweights(iframe) * cmat_shifted
                endif
            endif
        end do
        !$omp end do nowait
        !$omp sections
        !$omp section
        call movsum%set_cmat(cmat_sum_noshift)
        call movsum%ifft
        !$omp section
        call movsum_ctf%set_cmat(cmat_sum)
        call movsum_ctf%ifft
        !$omp section
        call movsum_corrected%set_cmat(cmat_wsum)
        call movsum_corrected%ifft
        !$omp section
        if( do_fromto )then
            call movsum_fromto%new(ldim_scaled, smpd_scaled)
            call movsum_fromto%set_cmat(cmat_wsum_fromto)
        endif
        !$omp end sections nowait
        !$omp end parallel
        ! destroy shifted frames
        do iframe=1,nframes
            call frames_shifted(iframe)%kill
        end do
        ! deallocate
        if( allocated(filtarr) ) deallocate(filtarr)
        deallocate(frames_shifted, cmat, cmat_shifted, cmat_sum_noshift, cmat_sum, cmat_wsum, cmat_wsum_fromto)
    end subroutine motion_correct_calc_sums

    subroutine motion_correct_calc_sums_tomo( movsum, movsum_ctf, movsum_corrected, frame_counter, time_per_frame )
        class(image),           intent(out)   :: movsum, movsum_ctf, movsum_corrected
        integer,                intent(inout) :: frame_counter
        real,                   intent(in)    :: time_per_frame
        real, allocatable :: filtarr(:,:)
        integer :: sz, iframe
        real    :: current_time, acc_dose
        sz = movie_frames_scaled(1)%get_filtsz()
        allocate(filtarr(nframes,sz))
        do iframe=1,nframes
            frame_counter     = frame_counter + 1
            current_time      = real(frame_counter) * time_per_frame ! unit: seconds
            acc_dose          = dose_rate * current_time             ! unit e/A2
            filtarr(iframe,:) = acc_dose2filter(movie_frames_scaled(1), acc_doses(iframe), kV)
        end do
        call  motion_correct_calc_sums( movsum, movsum_ctf, movsum_corrected )
    end subroutine motion_correct_calc_sums_tomo

    ! subroutine sum_movie_frames( shifts )
    !     real, intent(in), optional :: shifts(nframes,2)
    !     integer     :: iframe
    !     real        :: w
    !     logical     :: doshift
    !     type(image) :: frame_tmp
    !     doshift = present(shifts)
    !     call movie_sum_global%new(ldim_scaled, smpd_scaled)
    !     call movie_sum_global%set_ft(.true.)
    !     call frame_tmp%new(ldim_scaled, smpd_scaled)
    !     w = 1./real(nframes)
    !     do iframe=1,nframes
    !         if( doshift )then
    !             frame_tmp = movie_frames_scaled(iframe)
    !             call frame_tmp%shift([-shifts(iframe,1),-shifts(iframe,2),0.])
    !             call movie_sum_global%add(frame_tmp, w=w)
    !         else
    !             call movie_sum_global%add(movie_frames_scaled(iframe), w=w)
    !         endif
    !     end do
    !     call frame_tmp%kill
    ! end subroutine sum_movie_frames
    !
    ! subroutine wsum_movie_frames( shifts, fromto )
    !     real,              intent(in) :: shifts(nframes,2)
    !     integer, optional, intent(in) :: fromto(2)
    !     real, allocatable :: filter(:)
    !     type(image), allocatable :: frames_tmp(:)
    !     complex,     allocatable :: cmat(:,:,:), cmat_sum(:,:,:)
    !     integer :: iframe, ffromto(2), shp(3)
    !     call movie_sum_global%new(ldim_scaled, smpd_scaled)
    !     call movie_sum_global%set_ft(.true.)
    !     ffromto(1) = 1
    !     ffromto(2) = nframes
    !     if( present(fromto) ) ffromto = fromto
    !     shp = movie_sum_global%get_array_shape()
    !     allocate(frames_tmp(ffromto(1):ffromto(2)), cmat(shp(1),shp(2),shp(3)), cmat_sum(shp(1),shp(2),shp(3)))
    !     do iframe=ffromto(1),ffromto(2)
    !         call frames_tmp(iframe)%new(ldim_scaled, smpd_scaled)
    !     end do
    !     cmat_sum = cmplx(0.,0.)
    !     !$omp parallel do schedule(static) default(shared) private(iframe,filter,cmat) proc_bind(close) reduction(+:cmat_sum)
    !     do iframe=ffromto(1),ffromto(2)
    !         if( frameweights(iframe) > 0. )then
    !             call movie_frames_scaled(iframe)%shift2Dserial([-shifts(iframe,1),-shifts(iframe,2)], frames_tmp(iframe))
    !             if( do_dose_weight )then
    !                 filter = acc_dose2filter(frames_tmp(iframe), acc_doses(iframe), kV)
    !                 call frames_tmp(iframe)%apply_filter_serial(filter)
    !             endif
    !             call frames_tmp(iframe)%get_cmat_sub(cmat)
    !             cmat_sum = cmat_sum + cmat * frameweights(iframe)
    !         endif
    !     end do
    !     !$omp end parallel do
    !     call movie_sum_global%set_cmat(cmat_sum)
    !     do iframe=ffromto(1),ffromto(2)
    !         call frames_tmp(iframe)%kill
    !     end do
    !     deallocate(frames_tmp, cmat, cmat_sum)
    ! end subroutine wsum_movie_frames
    !
    ! subroutine wsum_movie_frames_tomo( shifts, frame_counter, time_per_frame )
    !     real,    intent(in)    :: shifts(nframes,2)
    !     integer, intent(inout) :: frame_counter
    !     real,    intent(in)    :: time_per_frame
    !     real, allocatable :: filter(:)
    !     type(image) :: frame_tmp
    !     integer :: iframe
    !     real    :: current_time, acc_dose
    !     call movie_sum_global%new(ldim_scaled, smpd_scaled)
    !     call movie_sum_global%set_ft(.true.)
    !     call frame_tmp%new(ldim_scaled, smpd_scaled)
    !     do iframe=1,nframes
    !         frame_counter = frame_counter + 1
    !         current_time  = real(frame_counter)*time_per_frame ! unit: seconds
    !         acc_dose      = dose_rate*current_time             ! unit e/A2
    !         if( frameweights(iframe) > 0. )then
    !             frame_tmp = movie_frames_scaled(iframe)
    !             call frame_tmp%shift([-shifts(iframe,1),-shifts(iframe,2),0.])
    !             filter = acc_dose2filter(movie_frames_scaled(iframe), acc_dose, kV)
    !             call frame_tmp%apply_filter(filter)
    !             call movie_sum_global%add(frame_tmp, w=frameweights(iframe))
    !             deallocate(filter)
    !         endif
    !     end do
    !     call frame_tmp%kill
    ! end subroutine wsum_movie_frames_tomo

    subroutine motion_correct_kill
        integer :: iframe
        if( existence )then
            do iframe=1,nframes
                call movie_frames_ftexp(iframe)%kill
                call movie_frames_ftexp_sh(iframe)%kill
                call movie_frames(iframe)%kill
                call movie_frames_scaled(iframe)%kill
                call movie_sum_global_ftexp_threads(iframe)%kill
            end do
            call movie_sum_global%kill
            deallocate( movie_frames_ftexp, movie_frames_ftexp_sh, movie_frames, movie_frames_scaled,&
            frameweights, frameweights_saved, corrs, opt_shifts, opt_shifts_saved,&
            movie_sum_global_ftexp_threads)
            if( allocated(acc_doses) ) deallocate(acc_doses)
            existence = .false.
        endif
    end subroutine motion_correct_kill

end module simple_motion_correct
