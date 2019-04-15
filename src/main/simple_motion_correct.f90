! motion_correct does motion correction, dose-weighting and frame-weighting of DDD movies
module simple_motion_correct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ft_expanded,         only: ft_expanded
use simple_motion_anisocor,     only: motion_anisocor, POLY_DIM
use simple_motion_patched,      only: motion_patched
use simple_motion_anisocor_dbl, only: motion_anisocor_dbl
use simple_motion_align_iso,    only: motion_align_iso
use simple_image,               only: image
use simple_parameters,          only: params_glob
use simple_estimate_ssnr,       only: acc_dose2filter
use simple_oris,                only: oris
use simple_opt_lbfgsb,          only: PRINT_NEVALS
implicit none

public :: motion_correct_iso, motion_correct_iso_calc_sums, motion_correct_iso_calc_sums_tomo, motion_correct_iso_kill
public :: motion_correct_aniso, motion_correct_aniso_calc_sums, motion_correct_aniso_kill, motion_correct_kill_common
public :: motion_correct_with_aniso, motion_correct_with_patched
public :: motion_correct_patched, motion_correct_patched_calc_sums, motion_correct_patched_kill
public :: patched_shift_fname
private
#include "simple_local_flags.inc"

interface motion_correct_iso_calc_sums
    module procedure motion_correct_iso_calc_sums_1
    module procedure motion_correct_iso_calc_sums_2
end interface

interface motion_correct_aniso_calc_sums
    module procedure motion_correct_aniso_calc_sums_1
    module procedure motion_correct_aniso_calc_sums_2
end interface

interface motion_correct_patched_calc_sums
    module procedure motion_correct_patched_calc_sums_1
    module procedure motion_correct_patched_calc_sums_2
end interface

! data structures for isotropic correction
type(image), target, allocatable :: movie_frames_scaled(:)            !< scaled movie frames
real,                allocatable :: opt_shifts(:,:)                   !< optimal shifts identified
real,                allocatable :: shifts_toplot(:,:)                !< shifts for plotting (in patch-based mc)
integer                          :: updateres
logical                          :: didupdateres

! data structures for patch-based motion correction
type(motion_patched) :: motion_patch
type(image),            allocatable :: movie_frames_shifted_patched(:) !< shifted movie frames

! data structures for anisotropic correction
class(motion_anisocor), allocatable :: motion_aniso(:)
type(image),            allocatable :: movie_frames_shifted_aniso(:) !< shifted movie frames
type(image),            allocatable :: movie_frames_shifted_saved(:) !< shifted movie frames
type(image),            allocatable :: movie_sum_global_threads(:)   !< array of global movie sums for parallel refinement (image)
real(dp),               allocatable :: opt_shifts_aniso(:,:)         !< optimal anisotropic shifts identified
real,                   allocatable :: opt_shifts_aniso_sp(:,:)      !< optimal anisotropic shifts identified, single precision
real(dp),               allocatable :: opt_shifts_aniso_saved(:,:)   !< optimal anisotropic shifts for local opt saved

! data structures used by both isotropic, anisotropic & patch-based correction
type(image), allocatable :: movie_frames_shifted(:)      !< shifted movie frames
type(image)              :: movie_sum_global             !< global movie sum for output
real,        allocatable :: corrs(:)                     !< per-frame correlations
real,        allocatable :: frameweights(:)              !< array of frameweights
real,        allocatable :: frameweights_saved(:)        !< array of frameweights
real,        allocatable :: acc_doses(:)                 !< accumulated doses
complex,     allocatable :: cmat(:,:,:), cmat_sum(:,:,:) !< complex matrices for OpenMP reduction

! module global variables
integer :: nframes        = 0        !< number of frames
integer :: fixed_frame    = 0        !< fixed frame of reference (0,0)
integer :: ldim(3)        = [0,0,0]  !< logical dimension of frame
integer :: ldim_scaled(3) = [0,0,0]  !< shrunken logical dimension of frame
real    :: hp             = 0.       !< high-pass limit
real    :: lp             = 0.       !< low-pass limit
real    :: resstep        = 0.       !< resolution step size (in angstrom)
real    :: smpd           = 0.       !< sampling distance
real    :: smpd_scaled    = 0.       !< sampling distance
real    :: corr_saved     = 0.       !< opt corr for local opt saved
real    :: kV             = 300.     !< acceleration voltage
real    :: dose_rate      = 0.       !< dose rate
real    :: nsig_here      = 6.0      !< nr of sigmas (for outlier removal)
logical :: do_dose_weight = .false.  !< dose weight or not
logical :: do_scale       = .false.  !< scale or not
logical :: motion_correct_with_aniso     = .false.  !< run aniso or not
logical :: motion_correct_with_patched   = .false.  !< run patch-based aniso or not
character(len=:), allocatable :: patched_shift_fname    !< file name for shift plot for patched-based alignment

! module global constants
integer, parameter :: MITSREF       = 30   !< max # iterations of refinement optimisation
integer, parameter :: MITSREF_ANISO = 9      !< max # iterations of anisotropic refinement optimisation
real,    parameter :: SMALLSHIFT    = 2.     !< small initial shift to blur out fixed pattern noise
logical, parameter :: ANISO_DBL     = .true. !< use double-precision for anisotropic motion correct?

contains

    ! PUBLIC METHODS, ISOTROPIC MOTION CORRECTION

    subroutine motion_correct_init( movie_stack_fname, ctfvars, err, gainref_fname )
        character(len=*),           intent(in)  :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(in)  :: ctfvars           !< CTF parameters
        logical,                    intent(out) :: err               !< error flag
        character(len=*), optional, intent(in)  :: gainref_fname     !< gain reference filename
        type(image), allocatable :: movie_frames(:)
        real,        allocatable :: rmat(:,:,:), rmat_pad(:,:), win(:,:), rmat_sum(:,:,:)
        logical,     allocatable :: outliers(:,:)
        type(image)        :: tmpmovsum, gainref
        real               :: moldiam, dimo4, time_per_frame, current_time
        integer            :: iframe, ncured, deadhot(2), i, j, winsz, ldim_scaled_tmp(3), shp(3)
        integer, parameter :: HWINSZ = 6
        logical            :: l_gain
        ! get number of frames & dim from stack
        call find_ldim_nptcls(movie_stack_fname, ldim, nframes)
        err = .false.
        if( nframes < 2 )then
            err = .true.
            write(logfhandle,*) 'movie: ', trim(movie_stack_fname)
            THROW_WARN('nframes of movie < 2, aborting motion_correct')
            return
        endif
        if( nframes < 2 ) return
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        if( params_glob%scale < 0.99 )then
            ldim_scaled(1) = round2even(real(ldim(1))*params_glob%scale)
            ldim_scaled(2) = round2even(real(ldim(2))*params_glob%scale)
            ldim_scaled(3) = 1
            do_scale        = .true.
        else
            ldim_scaled = ldim
            do_scale     = .false.
        endif
        ! set sampling distance
        smpd        = ctfvars%smpd
        smpd_scaled = ctfvars%smpd/params_glob%scale
        ! set fixed frame (all others are shifted by reference to this at 0,0)
        ! fixed_frame = nint(real(nframes)/2.) align the frames wrt the central one
        fixed_frame = 1                       !align the frames wrt the first one
        ! set reslims
        dimo4   = (real(minval(ldim_scaled(1:2))) * smpd_scaled) / 4.
        moldiam = 0.7 * real(params_glob%box) * smpd_scaled
        hp      = min(dimo4,2000.)
        lp      = params_glob%lpstart
        resstep = (params_glob%lpstart - params_glob%lpstop) / 3.
        ! allocate abstract data structures
        allocate( movie_frames(nframes), movie_frames_scaled(nframes), movie_frames_shifted(nframes),&
                  stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 1; simple_motion_correct')
        if ( motion_correct_with_aniso ) then
            allocate( movie_frames_shifted_saved(nframes), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init 2; simple_motion_correct')
        end if
        if ( motion_correct_with_patched ) then
            allocate( shifts_toplot(nframes, 2), source=0., stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init 3; simple_motion_correct')
        end if
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
            call movie_frames_shifted(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
        end do
        ! retrieve array shapes
        shp   = movie_frames_scaled(1)%get_array_shape()
        winsz = 2 * HWINSZ + 1
        ! additional allocations
        allocate(cmat(shp(1),shp(2),shp(3)), cmat_sum(shp(1),shp(2),shp(3)), rmat(ldim(1),ldim(2),1),&
        &rmat_sum(ldim(1),ldim(2),1), rmat_pad(1-HWINSZ:ldim(1)+HWINSZ,1-HWINSZ:ldim(2)+HWINSZ),&
        &win(winsz,winsz), corrs(nframes), opt_shifts(nframes,2), &
        &frameweights(nframes), frameweights_saved(nframes), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 4; simple_motion_correct')
        ! init
        cmat               = cmplx(0.,0.)
        cmat_sum           = cmplx(0.,0.)
        rmat               = 0.
        rmat_sum           = 0.
        rmat_pad           = 0.
        win                = 0.
        corrs              = 0.
        opt_shifts         = 0.
        frameweights       = 1./real(nframes)
        frameweights_saved = frameweights
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
        write(logfhandle,'(a)') '>>> REMOVING DEAD/HOT PIXELS & FOURIER TRANSFORMING FRAMES'
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
        !!!!!!!!!!!!!!!!!!!!!!!!!CHIARA!!!!!!!!!!!!!!!!!!!!!!!!
        !call tmpmovsum%write('TMPmovSUM.mrc')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call tmpmovsum%kill
        write(logfhandle,'(a,1x,i7)') '>>> # DEAD PIXELS:', deadhot(1)
        write(logfhandle,'(a,1x,i7)') '>>> # HOT  PIXELS:', deadhot(2)
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
    end subroutine motion_correct_init

    subroutine motion_correct_iso_callback( aPtr, align_iso, converged )
        class(*), pointer,       intent(inout) :: aPtr
        class(motion_align_iso), intent(inout) :: align_iso
        logical,                 intent(out)   :: converged
        integer :: nimproved, iter
        real    :: corrfrac, frac_improved
        corrfrac      = align_iso%get_corrfrac()
        frac_improved = align_iso%get_frac_improved()
        nimproved     = align_iso%get_nimproved()
        iter          = align_iso%get_iter()
        didupdateres  = .false.
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
            if( nimproved == 0 .and. iter > 2 )     converged = .true.
            if( iter > 10 .and. corrfrac > 0.9999 ) converged = .true.
        else
            converged = .false.
        end if

    contains

        subroutine update_res( thres_corrfrac, thres_frac_improved, which_update )
            real,    intent(in) :: thres_corrfrac, thres_frac_improved
            integer, intent(in) :: which_update
            if( corrfrac > thres_corrfrac .and. frac_improved <= thres_frac_improved&
                .and. updateres == which_update )then
                lp = lp - resstep
                call align_iso%set_hp_lp(hp, lp)
                write(logfhandle,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', lp
                ! need to indicate that we updated resolution limit
                updateres  = updateres + 1
                ! indicate that reslim was updated
                didupdateres = .true.
            endif
        end subroutine update_res

    end subroutine motion_correct_iso_callback

    !> isotropic motion_correction of DDD movie
    subroutine motion_correct_iso( movie_stack_fname, ctfvars, shifts, gainref_fname, nsig )
        character(len=*),           intent(in)    :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(inout) :: ctfvars           !< CTF params
        real,          allocatable, intent(out)   :: shifts(:,:)       !< the nframes shifts identified
        character(len=*), optional, intent(in)    :: gainref_fname     !< gain reference filename
        real,             optional, intent(in)    :: nsig              !< # sigmas (for outlier removal)
        real    :: ave, sdev, var, minw, maxw, corr
        logical :: err, err_stat
        type(motion_align_iso) :: align_iso
        class(*), pointer :: callback_ptr
        callback_ptr => null()
        ! initialise
        call align_iso%new
        nsig_here = 6.0
        if( present(nsig) ) nsig_here = nsig
        call motion_correct_init(movie_stack_fname, ctfvars, err, gainref_fname)
        if( err ) return
        call align_iso%set_frames(movie_frames_scaled, nframes)
        call align_iso%set_hp_lp(hp,lp)
        updateres = 0
        call align_iso%set_trs(params_glob%scale*params_glob%trs)
        call align_iso%set_smallshift(2.)
        call align_iso%set_rand_init_shifts(.true.)
        call align_iso%set_ftol_gtol(1e-7, 1e-7)
        call align_iso%set_shsrch_tol(1e-5)
        call align_iso%set_callback( motion_correct_iso_callback )
        call align_iso%align(callback_ptr)
        corr = align_iso%get_corr()
        if( corr < 0. )then
            write(logfhandle,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
            THROW_WARN('OPTIMAL CORRELATION < 0.0')
        endif
        call align_iso%get_opt_shifts(opt_shifts)
        if( allocated(shifts) ) deallocate(shifts)
        allocate(shifts(nframes,2), source=opt_shifts)
        call align_iso%get_weights(frameweights)
        call align_iso%get_shifts_toplot(shifts_toplot)
        call moment(frameweights, ave, sdev, var, err_stat)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        write(logfhandle,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT :', ave
        write(logfhandle,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS:', sdev
        write(logfhandle,'(a,7x,f7.4)') '>>> MIN WEIGHT     :', minw
        write(logfhandle,'(a,7x,f7.4)') '>>> MAX WEIGHT     :', maxw
        ! report the sampling distance of the possibly scaled movies
        ctfvars%smpd = smpd_scaled
        call align_iso%kill
    end subroutine motion_correct_iso

    subroutine motion_correct_iso_calc_sums_1( movie_sum, movie_sum_corrected, movie_sum_ctf )
        type(image), intent(out) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        integer :: iframe
        ! generate straight integrated movie frame for comparison
        call sum_movie_frames
        movie_sum = movie_sum_global
        call movie_sum%ifft()
        ! calculate the sum for CTF estimation
        call sum_movie_frames(opt_shifts)
        movie_sum_ctf = movie_sum_global
        ! save the isotropically corrected movie stack to disk for anisotropic movie alignment
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%ifft
            if (motion_correct_with_aniso) then
                movie_frames_shifted_saved(iframe) = movie_frames_shifted(iframe)
            end if
        end do
        call movie_sum_ctf%ifft()
        ! re-calculate the weighted sum
        call wsum_movie_frames(opt_shifts)
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%ifft()
    end subroutine motion_correct_iso_calc_sums_1

    subroutine motion_correct_iso_calc_sums_2( movie_sum_corrected, fromto )
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
    end subroutine motion_correct_iso_calc_sums_2

    subroutine motion_correct_iso_calc_sums_tomo( frame_counter, time_per_frame, movie_sum, movie_sum_corrected, movie_sum_ctf )
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
    end subroutine motion_correct_iso_calc_sums_tomo

    subroutine motion_correct_iso_kill
        integer :: iframe
        if( allocated(movie_frames_scaled) )then
            do iframe=1,size(movie_frames_scaled)
                call movie_frames_scaled(iframe)%kill
            end do
            deallocate(movie_frames_scaled)
        endif
        if( allocated(opt_shifts) )       deallocate(opt_shifts)
    end subroutine motion_correct_iso_kill

    ! PUBLIC METHODS, ANISOTROPIC MOTION CORRECTION

    !> anisotropic motion_correction of DDD movie
    subroutine motion_correct_aniso( shifts )
        real, allocatable, intent(out) :: shifts(:,:) !< the nframes polynomial models identifie
        real    :: corr, cxy(POLY_DIM + 1), ave, sdev, var, minw, maxw, corr_prev, frac_improved, corrfrac
        real    :: scale, smpd4scale
        integer :: iframe, nimproved, i, ldim4scale(3)
        logical :: didsave, err_stat, doscale
        real(dp):: acorr, aregu
        smpd4scale = params_glob%lpstop / 3.
        doscale = .false.
        if( smpd4scale > params_glob%smpd )then
            doscale = .true.
            scale   = smpd_scaled / smpd4scale
            ldim4scale(1) = round2even(scale * real(ldim_scaled(1)))
            ldim4scale(2) = round2even(scale * real(ldim_scaled(2)))
            ldim4scale(3) = 1
        else
            ldim4scale = ldim_scaled
            smpd4scale = smpd_scaled
        endif
        if( allocated(shifts) ) deallocate(shifts)
        ! construct
        if (ANISO_DBL) then
            allocate(motion_anisocor_dbl :: motion_aniso(nframes))
        else
            !allocate(motion_anisocor_sgl :: motion_aniso(nframes))
        end if
        allocate(movie_frames_shifted_aniso(nframes), movie_sum_global_threads(nframes),&
            &opt_shifts_aniso(nframes,POLY_DIM), opt_shifts_aniso_sp(nframes,POLY_DIM), &
            &opt_shifts_aniso_saved(nframes,POLY_DIM), shifts(nframes,POLY_DIM))
        opt_shifts_aniso       = 0._dp
        opt_shifts_aniso_sp    = 0.
        opt_shifts_aniso_saved = 0._dp
        shifts                 = 0.
        ! prepare the images for anisotropic correction by scaling & band-pass filtering & create motion_aniso objs
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%new(ldim_scaled, smpd_scaled)
            call movie_frames_shifted(iframe)%copy(movie_frames_shifted_saved(iframe))
            call movie_frames_shifted(iframe)%fft
            if( doscale ) call movie_frames_shifted(iframe)%clip_inplace(ldim4scale)
            call movie_frames_shifted(iframe)%bp(hp, lp)
            call movie_frames_shifted(iframe)%ifft
            call movie_sum_global_threads(iframe)%new(ldim4scale, smpd4scale, wthreads=.false.)
            call movie_frames_shifted_aniso(iframe)%copy(movie_frames_shifted(iframe))
            call motion_aniso(iframe)%new
        end do
        ! generate first average
        call gen_aniso_wsum_calc_corrs(0)
        ! calc avg corr to weighted avg
        corr = sum(corrs)/real(nframes)
        write(logfhandle,'(a)') '>>> ANISOTROPIC REFINEMENT'
        corr_saved = -1.
        didsave    = .false.
        PRINT_NEVALS = .true.
        do i=1,MITSREF_ANISO
            nimproved = 0
            !$omp parallel do schedule(static) default(shared) private(iframe,cxy) proc_bind(close) reduction(+:nimproved)
            do iframe=1,nframes
                ! subtract the movie frame being correlated to reduce bias
                call movie_sum_global_threads(iframe)%subtr(movie_frames_shifted_aniso(iframe), w=frameweights(iframe))
                ! optimise deformation
                cxy = real(motion_aniso(iframe)%minimize(ref=movie_sum_global_threads(iframe), &
                    frame=movie_frames_shifted(iframe), corr=acorr, regu=aregu))
                ! update parameter arrays
                opt_shifts_aniso(iframe,:) = cxy(2:POLY_DIM + 1)
                ! apply deformation
                call motion_aniso(iframe)%calc_T_out(opt_shifts_aniso(iframe,:), &
                    frame_in=movie_frames_shifted(iframe), frame_out=movie_frames_shifted_aniso(iframe)) !revisit this
                ! no need to add the frame back to the weighted sum since the sum will be updated after the loop
                ! (see below)
            end do
            !$omp end parallel do
            corr_prev = corr
            call gen_aniso_wsum_calc_corrs(i)
            frac_improved = real(nimproved) / real(nframes) * 100.
            write(logfhandle,'(a,1x,f4.0)') 'This % of frames improved their alignment: ', frac_improved
            corr = sum(corrs) / real(nframes)
            if( corr >= corr_saved )then ! save the local optimum
                corr_saved             = corr
                opt_shifts_aniso_saved = opt_shifts_aniso
                frameweights_saved     = frameweights
                didsave                = .true.
            endif
            corrfrac = corr_prev / corr
            if( nimproved == 0 .and. i > 2 )  exit
            if( i > 5 .and. corrfrac > 0.9999 ) exit
        end do
        PRINT_NEVALS = .false.
        ! put the best local optimum back
        corr             = corr_saved
        opt_shifts_aniso = opt_shifts_aniso_saved
        frameweights     = frameweights_saved
        ! output shifts
        if( allocated(shifts) ) deallocate(shifts)
        opt_shifts_aniso_sp = real(opt_shifts_aniso)
        allocate(shifts(nframes,POLY_DIM), source=opt_shifts_aniso_sp)
        ! print
        if( corr < 0. )then
            write(logfhandle,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
            THROW_WARN('OPTIMAL CORRELATION < 0.0')
        endif
        call moment(frameweights, ave, sdev, var, err_stat)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        write(logfhandle,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT :', ave
        write(logfhandle,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS:', sdev
        write(logfhandle,'(a,7x,f7.4)') '>>> MIN WEIGHT     :', minw
        write(logfhandle,'(a,7x,f7.4)') '>>> MAX WEIGHT     :', maxw
        ! put back the unfiltered isotropically corrected frames for final anisotropic sum generation
        do iframe=1,nframes
            call movie_frames_shifted_aniso(iframe)%copy(movie_frames_shifted_saved(iframe))
        end do
        ! apply nonlinear transformations
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call motion_aniso(iframe)%calc_T_out(opt_shifts_aniso(iframe,:), &
                frame_in=movie_frames_shifted_saved(iframe), frame_out=movie_frames_shifted_aniso(iframe))
        end do
        !$omp end parallel do
        ! we do not destruct search objects here (needed for sum production)

    contains

        subroutine gen_aniso_wsum_calc_corrs( imode )
            integer, intent(in) :: imode
            real, allocatable   :: rmat(:,:,:), rmat_sum(:,:,:)
            real    :: corr
            allocate( rmat(ldim4scale(1),ldim4scale(2),1), rmat_sum(ldim4scale(1),ldim4scale(2),1), source=0. )
            nimproved = 0
            ! FIRST LOOP TO OBTAIN WEIGHTED SUM
            !$omp parallel default(shared) private(iframe,rmat,corr) proc_bind(close)
            !$omp do schedule(static) reduction(+:rmat_sum)
            do iframe=1,nframes
                call movie_frames_shifted_aniso(iframe)%get_rmat_sub(rmat)
                rmat_sum = rmat_sum + rmat * frameweights(iframe)
            end do
            !$omp end do
            ! SECOND LOOP TO UPDATE movie_sum_global_threads AND CALCULATE CORRS
            !$omp do schedule(static) reduction(+:nimproved)
            do iframe=1,nframes
                ! update array of sums (for future parallel exec)
                call movie_sum_global_threads(iframe)%set_rmat(rmat_sum)
                ! subtract the movie frame being correlated to reduce bias
                call movie_sum_global_threads(iframe)%subtr(movie_frames_shifted_aniso(iframe), w=frameweights(iframe))
                ! calc corr
                if( imode == 0 )then
                    corrs(iframe) = movie_frames_shifted_aniso(iframe)%real_corr(movie_sum_global_threads(iframe))
                else
                    corr = movie_frames_shifted_aniso(iframe)%real_corr(movie_sum_global_threads(iframe))
                    if( corr > corrs(iframe) ) nimproved = nimproved + 1
                    corrs(iframe) = corr
                endif
                ! add the subtracted movie frame back to the weighted sum
                call movie_sum_global_threads(iframe)%add(movie_frames_shifted_aniso(iframe), w=frameweights(iframe))
            end do
            !$omp end do
            !$omp end parallel
            ! update frame weights
            frameweights = corrs2weights(corrs)
        end subroutine gen_aniso_wsum_calc_corrs

    end subroutine motion_correct_aniso

    subroutine motion_correct_aniso_calc_sums_1( movie_sum_corrected, movie_sum_ctf )
        type(image), intent(out) :: movie_sum_corrected, movie_sum_ctf
        call gen_aniso_sum(scalar_weight=1./real(nframes))
        movie_sum_ctf = movie_sum_global
        ! re-calculate the weighted sum
        call gen_aniso_sum
        movie_sum_corrected = movie_sum_global
    end subroutine motion_correct_aniso_calc_sums_1

    subroutine motion_correct_aniso_calc_sums_2( movie_sum_corrected, fromto )
        type(image), intent(out) :: movie_sum_corrected
        integer,     intent(in)  :: fromto(2)
        logical :: l_tmp
        ! re-calculate the weighted sum with dose_weighting turned off
        l_tmp = do_dose_weight
        do_dose_weight = .false.
        call gen_aniso_sum(fromto=fromto)
        movie_sum_corrected = movie_sum_global
        do_dose_weight = l_tmp
    end subroutine motion_correct_aniso_calc_sums_2

    subroutine motion_correct_aniso_kill
        integer :: iframe
        if( allocated(movie_frames_shifted_aniso) )then
            do iframe=1,size(movie_frames_shifted_aniso)
                call motion_aniso(iframe)%kill
                call movie_frames_shifted_aniso(iframe)%kill
                call movie_sum_global_threads(iframe)%kill
                call movie_frames_shifted_saved(iframe)%kill
            end do
            deallocate(motion_aniso, movie_frames_shifted_aniso, movie_sum_global_threads, movie_frames_shifted_saved)
        endif
        if( allocated(opt_shifts_aniso)       ) deallocate(opt_shifts_aniso)
        if( allocated(opt_shifts_aniso_sp)    ) deallocate(opt_shifts_aniso_sp)
        if( allocated(opt_shifts_aniso_saved) ) deallocate(opt_shifts_aniso_saved)
    end subroutine motion_correct_aniso_kill

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! PUBLIC METHODS, PATCH-BASED MOTION CORRECTION

   !> patch-based motion_correction of DDD movie
    subroutine motion_correct_patched()
        real    :: corr, corr_prev, corrfrac
        real    :: scale, smpd4scale
        integer :: iframe, ldim4scale(3)
        logical :: didsave, doscale
        call motion_patch%new(motion_correct_ftol = params_glob%motion_correctftol, &
            motion_correct_gtol = params_glob%motion_correctgtol, trs = params_glob%scale * params_glob%trs)
        smpd4scale = params_glob%smpd
        doscale = .false.
        if( smpd4scale > params_glob%smpd )then
            doscale = .true.
            scale   = smpd_scaled / smpd4scale
            ldim4scale(1) = round2even(scale * real(ldim_scaled(1)))
            ldim4scale(2) = round2even(scale * real(ldim_scaled(2)))
            ldim4scale(3) = 1
        else
            ldim4scale = ldim_scaled
            scale = 1.
            smpd4scale = smpd_scaled
        endif
        allocate(movie_frames_shifted_patched(nframes), movie_sum_global_threads(nframes))
        ! prepare the images for patch-based correction by scaling & band-pass filtering & create motion_aniso objs
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%new(ldim_scaled, smpd_scaled)
            call movie_frames_shifted(iframe)%copy(movie_frames_shifted_saved(iframe))
            call movie_sum_global_threads(iframe)%new(ldim4scale, smpd4scale, wthreads=.false.)
            call movie_frames_shifted_patched(iframe)%copy(movie_frames_shifted(iframe))
        end do
        ! generate first average
        call gen_patched_wsum_calc_corrs(0)
        ! calc avg corr to weighted avg
        corr = sum(corrs)/real(nframes)
        write(logfhandle,'(a)') '>>> PATCH-BASED REFINEMENT'
        corr_saved = -1.
        didsave    = .false.
        PRINT_NEVALS = .true.
        ! apply deformation
        call motion_patch%set_frameweights( frameweights )
        call motion_patch%correct( hp, lp, resstep, movie_frames_shifted, movie_frames_shifted_patched, patched_shift_fname, shifts_toplot )
        ! no need to add the frame back to the weighted sum since the sum will be updated after the loop
        ! (see below)
        corr_prev = corr
        call gen_patched_wsum_calc_corrs(0)
        corr = sum(corrs) / real(nframes)
        if( corr >= corr_saved )then ! save the local optimum
            corr_saved             = corr
            frameweights_saved     = frameweights
            didsave                = .true.
        endif
        corrfrac = corr_prev / corr

    contains

        subroutine gen_patched_wsum_calc_corrs( imode )
            integer, intent(in) :: imode
            real, allocatable   :: rmat(:,:,:), rmat_sum(:,:,:)
            real    :: corr
            allocate( rmat(ldim4scale(1),ldim4scale(2),1), rmat_sum(ldim4scale(1),ldim4scale(2),1), source=0. )
            ! FIRST LOOP TO OBTAIN WEIGHTED SUM
            !$omp parallel default(shared) private(iframe,rmat,corr) proc_bind(close)
            !$omp do schedule(static) reduction(+:rmat_sum)
            do iframe=1,nframes
                call movie_frames_shifted_patched(iframe)%get_rmat_sub(rmat)
                rmat_sum = rmat_sum + rmat * frameweights(iframe)
            end do
            !$omp end do
            ! SECOND LOOP TO UPDATE movie_sum_global_threads AND CALCULATE CORRS
            !$omp do schedule(static)
            do iframe=1,nframes
                ! update array of sums (for future parallel exec)
                call movie_sum_global_threads(iframe)%set_rmat(rmat_sum)
                ! subtract the movie frame being correlated to reduce bias
                call movie_sum_global_threads(iframe)%subtr(movie_frames_shifted_patched(iframe), w=frameweights(iframe))
                ! calc corr
                if( imode == 0 )then
                    corrs(iframe) = movie_frames_shifted_patched(iframe)%real_corr(movie_sum_global_threads(iframe))
                else
                    corr = movie_frames_shifted_patched(iframe)%real_corr(movie_sum_global_threads(iframe))
                endif
            end do
            !$omp end do
            !$omp end parallel
        end subroutine gen_patched_wsum_calc_corrs
    end subroutine motion_correct_patched

    subroutine motion_correct_patched_calc_sums_1( movie_sum_corrected, movie_sum_ctf )
        type(image), intent(out) :: movie_sum_corrected, movie_sum_ctf
        call gen_patched_sum(scalar_weight=1./real(nframes))
        movie_sum_ctf = movie_sum_global
        ! re-calculate the weighted sum
        call gen_patched_sum
        movie_sum_corrected = movie_sum_global
    end subroutine motion_correct_patched_calc_sums_1

    subroutine motion_correct_patched_calc_sums_2( movie_sum_corrected, fromto )
        type(image), intent(out) :: movie_sum_corrected
        integer,     intent(in)  :: fromto(2)
        logical :: l_tmp
        ! re-calculate the weighted sum with dose_weighting turned off
        l_tmp = do_dose_weight
        do_dose_weight = .false.
        call gen_patched_sum(fromto=fromto)
        movie_sum_corrected = movie_sum_global
        do_dose_weight = l_tmp
    end subroutine motion_correct_patched_calc_sums_2

    subroutine motion_correct_patched_kill
        integer :: iframe
        if( allocated(movie_frames_shifted_patched) )then
            do iframe=1,size(movie_frames_shifted_patched)
                call movie_frames_shifted_patched(iframe)%kill
                call movie_sum_global_threads(iframe)%kill
                call movie_frames_shifted_saved(iframe)%kill
            end do
            deallocate(movie_frames_shifted_patched, movie_sum_global_threads, movie_frames_shifted_saved)
        endif
        if ( allocated(shifts_toplot) ) deallocate(shifts_toplot)
    end subroutine motion_correct_patched_kill

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine motion_correct_kill_common
        integer :: iframe
        if( allocated(movie_frames_shifted) )then
            do iframe=1,size(movie_frames_shifted)
                call movie_frames_shifted(iframe)%kill
            end do
            deallocate(movie_frames_shifted)
        endif
        call movie_sum_global%kill
        if( allocated(corrs)              ) deallocate(corrs)
        if( allocated(frameweights)       ) deallocate(frameweights)
        if( allocated(frameweights_saved) ) deallocate(frameweights_saved)
        if( allocated(acc_doses)          ) deallocate(acc_doses)
        if( allocated(cmat)               ) deallocate(cmat)
        if( allocated(cmat_sum)           ) deallocate(cmat_sum)
    end subroutine motion_correct_kill_common

    ! PRIVATE UTILITY METHODS, ISO FIRST THEN ANISO

    subroutine sum_movie_frames( shifts )
        real, intent(in), optional :: shifts(nframes,2)
        integer :: iframe
        real    :: w
        logical :: doshift
        doshift = present(shifts)
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
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
        integer :: iframe, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = nframes
        if( present(fromto) ) ffromto = fromto
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
        if( do_dose_weight )then
            call gen_dose_weight_filter(filtarr)
        endif
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
    end subroutine wsum_movie_frames

    subroutine wsum_movie_frames_tomo( shifts, frame_counter, time_per_frame )
        real,    intent(in)    :: shifts(nframes,2)
        integer, intent(inout) :: frame_counter
        real,    intent(in)    :: time_per_frame
        real, allocatable :: filtarr(:,:)
        integer           :: iframe
        real              :: current_time, acc_dose
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        call movie_sum_global%set_ft(.true.)
        do iframe=1,nframes
            frame_counter     = frame_counter + 1
            current_time      = real(frame_counter)*time_per_frame ! unit: seconds
            acc_dose          = dose_rate*current_time             ! unit e/A2
        end do
        call gen_dose_weight_filter( filtarr )
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
    end subroutine wsum_movie_frames_tomo

    subroutine gen_dose_weight_filter( filtarr, movie_frames )
        real, allocatable,                          intent(inout) :: filtarr(:,:)
        type(image), allocatable, target, optional, intent(in)    :: movie_frames(:)
        type(image), pointer :: movie_frames_here(:)
        integer :: sz, filtsz, iframe
        if( allocated(filtarr) )deallocate(filtarr)
        if( present(movie_frames) ) then
            movie_frames_here => movie_frames
        else
            movie_frames_here => movie_frames_scaled
        end if
        sz = movie_frames_here(1)%get_filtsz()
        filtsz = 2*sz ! the filter goes beyond nyquist
        allocate(filtarr(nframes,filtsz))
        do iframe=1,nframes
            filtarr(iframe,:) = acc_dose2filter(movie_frames_here(1), acc_doses(iframe), kV, filtsz)
        end do
    end subroutine gen_dose_weight_filter

    ! PRIVATE UTILITY METHODS, ANISO

    subroutine gen_aniso_sum( fromto, scalar_weight )
        integer, optional, intent(in) :: fromto(2)
        real,    optional, intent(in) :: scalar_weight
        real, allocatable :: rmat(:,:,:), rmat_sum(:,:,:), filtarr(:,:)
        integer :: iframe, ffromto(2)
        logical :: l_w_scalar
        ffromto(1) = 1
        ffromto(2) = nframes
        if( present(fromto) ) ffromto = fromto
        l_w_scalar = present(scalar_weight)
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        if( do_dose_weight )then
            call gen_dose_weight_filter(filtarr)
        endif
        allocate( rmat(ldim_scaled(1),ldim_scaled(2),1), rmat_sum(ldim_scaled(1),ldim_scaled(2),1), source=0. )
        !$omp parallel do default(shared) private(iframe,rmat) proc_bind(close) schedule(static) reduction(+:rmat_sum)
        do iframe=ffromto(1),ffromto(2)
            if( do_dose_weight )then
                call movie_frames_shifted_aniso(iframe)%fft
                call movie_frames_shifted_aniso(iframe)%apply_filter_serial(filtarr(iframe,:))
                call movie_frames_shifted_aniso(iframe)%ifft
            endif
            call movie_frames_shifted_aniso(iframe)%get_rmat_sub(rmat)
            if( l_w_scalar )then
                rmat_sum = rmat_sum + rmat * scalar_weight
            else
                rmat_sum = rmat_sum + rmat * frameweights(iframe)
            endif
        end do
        !$omp end parallel do
        call movie_sum_global%set_rmat(rmat_sum)
    end subroutine gen_aniso_sum

    ! PRIVATE UTILITY METHODS, PATCH-BASED

    subroutine gen_patched_sum( fromto, scalar_weight )
        integer, optional, intent(in) :: fromto(2)
        real,    optional, intent(in) :: scalar_weight
        real, allocatable :: rmat(:,:,:), rmat_sum(:,:,:), filtarr(:,:)
        integer :: iframe, ffromto(2)
        logical :: l_w_scalar
        logical :: ft_flag
        ffromto(1) = 1
        ffromto(2) = nframes
        if( present(fromto) ) ffromto = fromto
        l_w_scalar = present(scalar_weight)
        call movie_sum_global%new(ldim_scaled, smpd_scaled)
        if (.false.) then
            call movie_sum_global%set_ft(.true.)
            cmat_sum = cmplx(0.,0.)
            !$omp parallel do default(shared) private(iframe,cmat) proc_bind(close) schedule(static) reduction(+:cmat_sum)
            do iframe=ffromto(1),ffromto(2)
                if (.not. movie_frames_shifted_patched(iframe)%is_ft()) call movie_frames_shifted_patched(iframe)%fft()
                call movie_frames_shifted_patched(iframe)%get_cmat_sub(cmat)
                if( l_w_scalar )then
                    cmat_sum = cmat_sum + cmat * scalar_weight
                else
                    cmat_sum = cmat_sum + cmat * frameweights(iframe)
                endif
                call movie_frames_shifted_patched(iframe)%ifft()
            end do
            !$omp end parallel do
            call movie_sum_global%set_cmat(cmat_sum)
            call movie_sum_global%ifft()
        end if
        allocate( rmat(ldim_scaled(1),ldim_scaled(2),1), rmat_sum(ldim_scaled(1),ldim_scaled(2),1), source=0. )
        if ( do_dose_weight ) then
            ft_flag = movie_sum_global%is_ft()
            call movie_sum_global%set_ft(.true.)
            call gen_dose_weight_filter(filtarr, movie_frames_shifted_patched)
            call movie_sum_global%set_ft(ft_flag)
        end if
        !$omp parallel do default(shared) private(iframe,rmat) proc_bind(close) schedule(static) reduction(+:rmat_sum)
        do iframe=ffromto(1),ffromto(2)
            if( do_dose_weight )then
                call movie_frames_shifted_patched(iframe)%fft
                call movie_frames_shifted_patched(iframe)%apply_filter_serial(filtarr(iframe,:))
                call movie_frames_shifted_patched(iframe)%ifft
            endif
            call movie_frames_shifted_patched(iframe)%get_rmat_sub(rmat)
            if( l_w_scalar )then
                rmat_sum = rmat_sum + rmat * scalar_weight
            else
                rmat_sum = rmat_sum + rmat * frameweights(iframe)
            endif
        end do
        !$omp end parallel do
        call movie_sum_global%set_rmat(rmat_sum)
    end subroutine gen_patched_sum

end module simple_motion_correct
