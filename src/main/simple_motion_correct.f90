! motion_correct does motion correction, dose-weighting and frame-weighting of DDD movies
module simple_motion_correct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ft_expanded,                   only: ft_expanded, ftexp_transfmat_init, ftexp_transfmat_kill
use simple_motion_patched,                only: motion_patched
use simple_motion_align_hybrid,           only: motion_align_hybrid
use simple_motion_align_iso_polyn_direct, only: motion_align_iso_polyn_direct
!use simple_motion_align_iso_direct,       only: motion_align_iso_direct
use simple_image,                         only: image
use simple_parameters,                    only: params_glob
use simple_estimate_ssnr,                 only: acc_dose2filter
use simple_opt_lbfgsb,                    only: PRINT_NEVALS
use simple_starfile_wrappers
implicit none

! Stage drift
public :: motion_correct_iso, motion_correct_iso_calc_sums, motion_correct_iso_calc_sums_tomo
public :: write_iso2star, motion_correct_iso_kill
! Beam-induced motion correction
public :: motion_correct_patched, motion_correct_patched_calc_sums, motion_correct_patched_kill
public :: write_aniso2star, motion_correct_with_patched
! Common & convenience
public :: motion_correct_kill_common, close_starfile, patched_shift_fname
private
#include "simple_local_flags.inc"

! data structures for isotropic correction
type(image), target, allocatable :: movie_frames_scaled(:)            !< scaled movie frames
real,                allocatable :: opt_shifts(:,:)                   !< optimal shifts identified
real                             :: TOL_ISO = 1.e-6                   !< LBFGSB tolerance for isotropic search
integer                          :: updateres
logical                          :: didupdateres

! data structures for patch-based motion correction
type(motion_patched)        :: motion_patch
real(dp),       allocatable :: patched_polyn(:)                !< polynomial from patched-based correction

! data structures used by both isotropic & patch-based correction
type(image), allocatable :: movie_frames_shifted_saved(:) !< shifted movie frames
real,        allocatable :: corrs(:)                     !< per-frame correlations
real,        allocatable :: shifts_toplot(:,:)           !< shifts for plotting & parsing
real,        allocatable :: frameweights(:)              !< array of frameweights
real,        allocatable :: acc_doses(:)                 !< accumulated doses
complex,     allocatable :: cmat_sum(:,:,:)              !< complex matrices for OpenMP reduction

! module global variables
integer :: nframes        = 0        !< number of frames
integer :: fixed_frame    = 0        !< fixed frame of reference for isotropic alignment (0,0)
integer :: ldim(3)        = [0,0,0]  !< logical dimension of frame
integer :: ldim_orig(3)   = [0,0,0]  !< logical dimension of frame (original, for use in motion_correct_iter)
integer :: ldim_scaled(3) = [0,0,0]  !< shrunken logical dimension of frame
real    :: hp             = 0.       !< high-pass limit
real    :: lp             = 0.       !< low-pass limit
real    :: resstep        = 0.       !< resolution step size (in angstrom)
real    :: smpd           = 0.       !< sampling distance
real    :: smpd_scaled    = 0.       !< sampling distance
real    :: kV             = 300.     !< acceleration voltage
real    :: dose_rate      = 0.       !< dose rate
logical :: do_scale                    = .false.  !< scale or not
logical :: motion_correct_with_patched = .false.  !< run patch-based aniso or not
character(len=:), allocatable :: patched_shift_fname    !< file name for shift plot for patched-based alignment
character(len=:), allocatable :: mc_starfile_fname      !< file name for starfile rel. to motion correct
type(starfile_table_type)     :: mc_starfile            !< starfile for motion correct output

! module global constants
real,    parameter :: NSIGMAS                       = 5.       !< Number of standard deviations for outliers detection
real,    parameter :: SMALLSHIFT                    = 1.       !< small initial shift to blur out fixed pattern noise
logical, parameter :: FITSHIFTS                     = .true.
logical, parameter :: ISO_POLYN_DIRECT              = .false.  !< use polynomial constraint for isotropic motion correction
logical, parameter :: ISO_UNCONSTR_AFTER            = .false.  !< run a unconstrained (direct) as the second step (at highest resolution)
logical, parameter :: DO_PATCHED_POLYN              = .false.  !< run polynomially constrained motion correction for patch-based motion correction
logical, parameter :: DO_PATCHED_POLYN_DIRECT_AFTER = .false.  !< run a direct polynomial optimization for patch-based motion correction as the second step (at highest resolution)
contains

    ! PUBLIC METHODS, ISOTROPIC MOTION CORRECTION

    subroutine motion_correct_init( movie_stack_fname, ctfvars, nsig, err, gainref_fname )
        character(len=*),           intent(in)  :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(in)  :: ctfvars           !< CTF parameters
        real,                       intent(in)  :: nsig              !< # of std deviations for outliers detection
        logical,                    intent(out) :: err               !< error flag
        character(len=*), optional, intent(in)  :: gainref_fname     !< gain reference filename
        integer,       parameter :: HWINSZ = 6
        type(image), allocatable :: movie_frames(:)
        real,        allocatable :: rmat_pad(:,:), win(:,:), rmat_sum(:,:,:)
        logical,     allocatable :: outliers(:,:)
        real,      pointer :: prmat(:,:,:)
        type(image)        :: tmpmovsum, gainref
        real               :: moldiam, dimo4, time_per_frame, current_time
        integer            :: iframe, ncured, deadhot(2), i, j, winsz, shp(3)
        logical            :: l_gain
        ! get number of frames & dim from stack
        call find_ldim_nptcls(movie_stack_fname, ldim, nframes)
        ldim_orig = ldim
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
            do_scale    = .false.
        endif
        ! set sampling distance
        smpd        = ctfvars%smpd
        smpd_scaled = ctfvars%smpd/params_glob%scale
        ! set fixed frame to central one
        fixed_frame = nint(0.5*real(nframes))
        ! set reslims
        dimo4   = (real(minval(ldim_scaled(1:2))) * smpd_scaled) / 4.
        moldiam = 0.7 * real(params_glob%box) * smpd_scaled
        hp      = min(dimo4,2000.)
        lp      = params_glob%lpstart
        resstep = (params_glob%lpstart - params_glob%lpstop) / 3.
        ! allocate abstract data structures
        allocate( movie_frames(nframes), movie_frames_scaled(nframes),stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 1; simple_motion_correct')
        allocate( shifts_toplot(nframes, 2), source=0., stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 2; simple_motion_correct')
        if ( motion_correct_with_patched ) then
            allocate( movie_frames_shifted_saved(nframes), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('motion_correct_init 3; simple_motion_correct')
        end if
        ! retrieve array shapes
        winsz = 2 * HWINSZ + 1
        ! additional allocations
        allocate(win(winsz,winsz), corrs(nframes), opt_shifts(nframes,2),&
        &rmat_sum(ldim(1),ldim(2),1), rmat_pad(1-HWINSZ:ldim(1)+HWINSZ,1-HWINSZ:ldim(2)+HWINSZ),&
        &frameweights(nframes), source=0.,stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 5; simple_motion_correct')
        frameweights = 1./real(nframes)
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
        !$omp parallel do schedule(static) default(shared) private(iframe,prmat) proc_bind(close) reduction(+:rmat_sum)
        do iframe=1,nframes
            if( l_gain ) call movie_frames(iframe)%mul(gainref) ! gain correction
            call movie_frames(iframe)%get_rmat_ptr(prmat)
            rmat_sum(:,:,1) = rmat_sum(:,:,1) + prmat(1:ldim(1),1:ldim(2),1)
        end do
        !$omp end parallel do
        if( l_gain ) call gainref%kill
        call tmpmovsum%new(ldim, smpd)
        call tmpmovsum%set_rmat(rmat_sum)
        call tmpmovsum%cure_outliers(ncured, nsig, deadhot, outliers)
        call tmpmovsum%kill
        deallocate(rmat_sum)
        write(logfhandle,'(a,1x,i7)') '>>> # DEAD PIXELS:', deadhot(1)
        write(logfhandle,'(a,1x,i7)') '>>> # HOT  PIXELS:', deadhot(2)
        if( any(outliers) )then ! remove the outliers & do the rest
            !$omp parallel do schedule(static) default(shared) private(iframe,prmat,rmat_pad,i,j,win) proc_bind(close)
            do iframe=1,nframes
                call movie_frames(iframe)%get_rmat_ptr(prmat)
                rmat_pad(:,:) = median( reshape(prmat(1:ldim(1),1:ldim(2),1),(/(ldim(1)*ldim(2))/)) )
                rmat_pad(1:ldim(1),1:ldim(2)) = prmat(1:ldim(1),1:ldim(2),1)
                do i=1,ldim(1)
                    do j=1,ldim(2)
                        if( outliers(i,j) )then
                            win = rmat_pad( i-HWINSZ:i+HWINSZ, j-HWINSZ:j+HWINSZ )
                            prmat(i,j,1) = median( reshape(win,(/winsz**2/)) )
                        endif
                    enddo
                enddo
                call movie_frames(iframe)%fft()
                call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
                call movie_frames(iframe)%clip(movie_frames_scaled(iframe))
                call movie_frames(iframe)%kill
            enddo
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
            do iframe=1,nframes
                call movie_frames(iframe)%fft()
                call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
                call movie_frames(iframe)%clip(movie_frames_scaled(iframe))
                call movie_frames(iframe)%kill
            end do
            !$omp end parallel do
        endif
        ! additional allocations
        shp = movie_frames_scaled(1)%get_array_shape()
        allocate(cmat_sum(shp(1),shp(2),shp(3)), source=cmplx(0.,0.))
        ! check if we are doing dose weighting
        if( params_glob%l_dose_weight )then
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
        deallocate(rmat_pad, win, outliers, movie_frames)
    end subroutine motion_correct_init

    subroutine motion_correct_iso_polyn_direct_callback( aPtr, align_iso_polyn_direct, converged )
        class(*), pointer,                    intent(inout) :: aPtr
        class(motion_align_iso_polyn_direct), intent(inout) :: align_iso_polyn_direct
        logical,                              intent(out)   :: converged
        logical :: didupdateres
        didupdateres = .false.
        select case(updateres)
        case(0)
            call update_res( updateres )
        case(1)
            call update_res( updateres )
        case(2)
            call update_res( updateres )
            call align_iso_polyn_direct%set_factr_pgtol(1d+6, 1d-6)
        case DEFAULT
            ! nothing to do
        end select
        if( updateres > 2 .and. .not. didupdateres) then
            call align_iso_polyn_direct%reset_tols
            converged = .true.
        else
            converged = .false.
        end if

    contains

        subroutine update_res( which_update )
            integer, intent(in) :: which_update
            lp = lp - resstep
            call align_iso_polyn_direct%set_hp_lp(hp, lp)
            write(logfhandle,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', lp
            ! need to indicate that we updated resolution limit
            updateres  = updateres + 1
            ! indicate that reslim was updated
            didupdateres = .true.
        end subroutine update_res

    end subroutine motion_correct_iso_polyn_direct_callback

    !> isotropic motion_correction of DDD movie
    subroutine motion_correct_iso( movie_stack_fname, ctfvars, bfactor, shifts, gainref_fname, nsig )
        character(len=*),           intent(in)    :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(inout) :: ctfvars           !< CTF params
        real,                       intent(in)    :: bfactor           !< B-factor
        real,          allocatable, intent(out)   :: shifts(:,:)       !< the nframes shifts identified
        character(len=*), optional, intent(in)    :: gainref_fname     !< gain reference filename
        real,             optional, intent(in)    :: nsig              !< # sigmas (for outlier removal)
        real                                :: ave, sdev, var, minw, maxw, corr, nsig_here
        logical                             :: err, err_stat
        type(motion_align_hybrid)           :: hybrid_srch
        type(motion_align_iso_polyn_direct) :: align_iso_polyn_direct
        class(*), pointer                   :: callback_ptr
        callback_ptr => null()
        ! initialise
        if (ISO_POLYN_DIRECT) then
            call align_iso_polyn_direct%new
        end if
        nsig_here = NSIGMAS
        if( present(nsig) ) nsig_here = nsig
        call motion_correct_init(movie_stack_fname, ctfvars, nsig_here, err, gainref_fname)
        if( err ) return
        call ftexp_transfmat_init(movie_frames_scaled(1))
        if (ISO_POLYN_DIRECT) then
            call align_iso_polyn_direct%set_frames(movie_frames_scaled, nframes)
            call align_iso_polyn_direct%set_hp_lp(hp,lp)
            updateres = 0
            call align_iso_polyn_direct%set_callback( motion_correct_iso_polyn_direct_callback )
            call align_iso_polyn_direct%align_polyn(callback_ptr)
            if (ISO_UNCONSTR_AFTER) then
                call align_iso_polyn_direct%refine_direct
            end if
            corr = align_iso_polyn_direct%get_corr()
            if( corr < 0. )then
                write(logfhandle,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
                THROW_WARN('OPTIMAL CORRELATION < 0.0')
            endif
            call align_iso_polyn_direct%get_opt_shifts(opt_shifts)
            if( allocated(shifts) ) deallocate(shifts)
            allocate(shifts(nframes,2), source=opt_shifts)
            call align_iso_polyn_direct%get_weights(frameweights)
            call align_iso_polyn_direct%get_shifts_toplot(shifts_toplot)
            call align_iso_polyn_direct%kill
        else
            call hybrid_srch%new(movie_frames_scaled)
            call hybrid_srch%set_group_frames(.false.)
            call hybrid_srch%set_reslims(hp, params_glob%lpstart, params_glob%lpstop)
            call hybrid_srch%set_bfactor(bfactor)
            call hybrid_srch%set_trs(params_glob%scale*params_glob%trs)
            call hybrid_srch%set_rand_init_shifts(.true.)
            call hybrid_srch%set_shsrch_tol(TOL_ISO)
            call hybrid_srch%set_fitshifts(FITSHIFTS)
            call hybrid_srch%set_fixed_frame(fixed_frame)
            call hybrid_srch%align
            corr = hybrid_srch%get_corr()
            call hybrid_srch%get_opt_shifts(opt_shifts)
            call hybrid_srch%get_weights(frameweights)
            call hybrid_srch%get_shifts_toplot(shifts_toplot)
            call hybrid_srch%kill
            if( corr < 0. )then
               write(logfhandle,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
               THROW_WARN('OPTIMAL CORRELATION < 0.0')
            endif
            shifts = opt_shifts
        end if
        call moment(frameweights, ave, sdev, var, err_stat)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        write(logfhandle,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT :', ave
        write(logfhandle,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS:', sdev
        write(logfhandle,'(a,7x,f7.4)') '>>> MIN WEIGHT     :', minw
        write(logfhandle,'(a,7x,f7.4)') '>>> MAX WEIGHT     :', maxw
        ! report the sampling distance of the possibly scaled movies
        ctfvars%smpd = smpd_scaled
    end subroutine motion_correct_iso

    ! generates sums & shift frames
    subroutine motion_correct_iso_calc_sums( movie_sum, movie_sum_corrected, movie_sum_ctf)
        type(image), intent(inout) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        type(image),   allocatable :: movie_frames_shifted(:)
        real,          allocatable :: filtarr(:,:)
        real                       :: wvec(1:nframes), scalar_weight
        integer                    :: iframe
        scalar_weight = 1. / real(nframes)
        ! copy unaligned frames
         allocate(movie_frames_shifted(nframes))
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%new(ldim_scaled, smpd_scaled)
            call movie_frames_shifted(iframe)%copy(movie_frames_scaled(iframe))
        end do
        !$omp end parallel do
        ! Unaligned, unweighted, unfiltered sum for comparison
        wvec = scalar_weight
        call sum_frames(movie_sum, wvec)
        ! shift frames
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%shift2Dserial(-opt_shifts(iframe,:))
            if (motion_correct_with_patched) then
                ! Save the isotropically corrected frames for anisotropic alignment
                call movie_frames_shifted_saved(iframe)%new(ldim_scaled, smpd_scaled)
                call movie_frames_shifted_saved(iframe)%copy(movie_frames_shifted(iframe))
                call movie_frames_shifted_saved(iframe)%ifft
            endif
        end do
        !$omp end parallel do
        ! Following sums are performed in case anisotropic correction will be unsuccessfull
        ! Unweighted, unfiltered sum for CTF estimation
        wvec = scalar_weight
        where( frameweights < 1.e-6 ) wvec = 0.
        call sum_frames(movie_sum_ctf, wvec)
        ! Weighted, filtered sum for micrograph
        if( params_glob%l_dose_weight )then
            call gen_dose_weight_filter(filtarr, movie_frames_shifted)
            !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
            do iframe=1,nframes
                call movie_frames_shifted(iframe)%apply_filter_serial(filtarr(iframe,:))
            end do
            !$omp end parallel do
            deallocate(filtarr)
        endif
        call sum_frames(movie_sum_corrected,frameweights)
        ! cleanup
        do iframe=1,nframes
            call movie_frames_shifted(iframe)%kill
        end do
        deallocate(movie_frames_shifted)

        contains

            subroutine sum_frames( img_sum, weights )
                class(image), intent(inout) :: img_sum
                real,         intent(in)    :: weights(1:nframes)
                complex,  pointer :: pcmat(:,:,:)
                integer :: iframe
                call img_sum%new(ldim_scaled, smpd_scaled)
                cmat_sum = cmplx(0.,0.)
                !$omp parallel do default(shared) private(iframe,pcmat) proc_bind(close)&
                !$omp schedule(static) reduction(+:cmat_sum)
                do iframe=1,nframes
                    if( weights(iframe) < TINY ) cycle
                    call movie_frames_shifted(iframe)%get_cmat_ptr(pcmat)
                    cmat_sum = cmat_sum + pcmat * weights(iframe)
                end do
                !$omp end parallel do
                call img_sum%set_ft(.true.)
                call img_sum%set_cmat(cmat_sum)
                call img_sum%ifft()
                nullify(pcmat)
            end subroutine sum_frames
    end subroutine motion_correct_iso_calc_sums

    subroutine motion_correct_iso_calc_sums_tomo( frame_counter, time_per_frame, movie_sum, movie_sum_corrected, movie_sum_ctf )
        integer,     intent(inout) :: frame_counter  !< frame counter
        real,        intent(in)    :: time_per_frame !< time resolution
        type(image), intent(inout) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        real    :: current_time
        integer :: iframe
        ! re-calculates doses
        do iframe=1,nframes
            frame_counter     = frame_counter + 1
            current_time      = real(frame_counter)*time_per_frame ! unit: seconds
            acc_doses(iframe) = dose_rate*current_time             ! unit e/A2
        end do
        call  motion_correct_iso_calc_sums(movie_sum, movie_sum_corrected, movie_sum_ctf)
    end subroutine motion_correct_iso_calc_sums_tomo

    ! open, write isotropic shifts and do not close star file
    subroutine write_iso2star( mc_starfile_fname, moviename, gainref_fname )
        character(len=*),           intent(in) :: mc_starfile_fname, moviename
        character(len=*), optional, intent(in) :: gainref_fname
        real    :: doserateperframe
        integer :: iframe
        ! open as new file
        call starfile_table__new(mc_starfile)
        call starfile_table__open_ofile(mc_starfile, mc_starfile_fname)
        ! global fields
        call starfile_table__addObject(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .true.)
        call starfile_table__setname(mc_starfile, "general")
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_X, ldim_orig(1))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Y, ldim_orig(2))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Z, ldim_orig(3))
        call starfile_table__setValue_string(mc_starfile, EMDL_MICROGRAPH_MOVIE_NAME, simple_abspath(moviename))
        if (present(gainref_fname)) then
            call starfile_table__setValue_string(mc_starfile, EMDL_MICROGRAPH_GAIN_NAME, trim(gainref_fname))
        end if
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_BINNING, real(1./params_glob%scale,dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, real(smpd, dp))
        doserateperframe = 0.
        if( params_glob%l_dose_weight )then
            doserateperframe = params_glob%exp_time*params_glob%dose_rate   ! total dose
            doserateperframe = doserateperframe / real(nframes)             ! per frame
        endif
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_DOSE_RATE, real(doserateperframe, dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_PRE_EXPOSURE, 0.0_dp)
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_VOLTAGE, real(params_glob%kv, dp))
        call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_START_FRAME, 1)
        call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_MOTION_MODEL_VERSION, 1)
        call starfile_table__write_ofile(mc_starfile)
        ! isotropic shifts
        call starfile_table__clear(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .false.)
        call starfile_table__setName(mc_starfile, "global_shift")
        do iframe = 1, nframes
            call starfile_table__addObject(mc_starfile)
            call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, real(shifts_toplot(iframe, 1), dp))
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, real(shifts_toplot(iframe, 2), dp))
        end do
        call starfile_table__write_ofile(mc_starfile)
    end subroutine write_iso2star

    subroutine motion_correct_iso_kill
        integer :: iframe
        if( allocated(movie_frames_scaled) )then
            do iframe=1,size(movie_frames_scaled)
                call movie_frames_scaled(iframe)%kill
            end do
            deallocate(movie_frames_scaled)
        endif
        if (allocated(opt_shifts)) deallocate(opt_shifts)
        call ftexp_transfmat_kill
    end subroutine motion_correct_iso_kill


    ! PUBLIC METHODS, PATCH-BASED MOTION CORRECTION

    !> patch-based motion_correction of DDD movie
    subroutine motion_correct_patched( bfac, chisq )
        real, intent(in)  :: bfac
        real, intent(out) :: chisq(2)     !< whether polynomial fitting was within threshold
        real    :: scale, smpd4scale
        integer :: ldim4scale(3)
        logical :: doscale
        chisq = huge(chisq(1))
        call motion_patch%new(motion_correct_ftol = params_glob%motion_correctftol, &
            motion_correct_gtol = params_glob%motion_correctgtol, trs = params_glob%scale*params_glob%trs)
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
        write(logfhandle,'(A,I2,A3,I2,A1)') '>>> PATCH-BASED REFINEMENT (',&
            &params_glob%nxpatch,' x ',params_glob%nxpatch,')'
        PRINT_NEVALS = .false.
        if (DO_PATCHED_POLYN) then
            call motion_patch%correct_polyn( hp, resstep, movie_frames_shifted_saved,&
                &patched_shift_fname, DO_PATCHED_POLYN_DIRECT_AFTER, shifts_toplot)
            chisq = 0.
        else
            call motion_patch%set_fitshifts(FITSHIFTS)
            call motion_patch%set_frameweights( frameweights )
            call motion_patch%set_fixed_frame(fixed_frame)
            call motion_patch%set_interp_fixed_frame(fixed_frame)
            call motion_patch%set_bfactor(bfac)
            call motion_patch%correct(hp, resstep, movie_frames_shifted_saved,patched_shift_fname, shifts_toplot)
            chisq = motion_patch%get_polyfit_chisq()
        end if
        call motion_patch%get_poly4star(patched_polyn)
    end subroutine motion_correct_patched

    ! weighted & un-weighted sums, dose-weighting
    subroutine motion_correct_patched_calc_sums( movie_sum_corrected, movie_sum_ctf )
        type(image), intent(inout) :: movie_sum_corrected, movie_sum_ctf
        real,            pointer :: prmat(:,:,:)
        type(image), allocatable :: movie_frames_shifted_patched(:)
        real,        allocatable :: rmat_sum(:,:,:), filtarr(:,:)
        real              :: scalar_weight
        integer           :: iframe
        scalar_weight = 1. / real(nframes)
        call movie_sum_ctf%new(ldim_scaled, smpd_scaled)
        allocate(rmat_sum(ldim_scaled(1),ldim_scaled(2),1), source=0.)
        allocate(movie_frames_shifted_patched(nframes))
        ! micrograph for CTF estimation
        !$omp parallel do default(shared) private(iframe,prmat) proc_bind(close) schedule(static) reduction(+:rmat_sum)
        do iframe=1,nframes
            if( frameweights(iframe) < 1.e-6 ) cycle
            call motion_patch%polytransfo(iframe, movie_frames_shifted_saved(iframe), movie_frames_shifted_patched(iframe))
            call movie_frames_shifted_patched(iframe)%get_rmat_ptr(prmat)
            rmat_sum(:,:,1) = rmat_sum(:,:,1) + scalar_weight*prmat(1:ldim_scaled(1),1:ldim_scaled(2),1)
        end do
        !$omp end parallel do
        call movie_sum_ctf%set_rmat(rmat_sum)
        ! micrograph
        call movie_sum_corrected%new(ldim_scaled, smpd_scaled)
        if( params_glob%l_dose_weight ) call gen_dose_weight_filter(filtarr, movie_frames_shifted_patched)
        rmat_sum = 0.
        !$omp parallel do default(shared) private(iframe,prmat) proc_bind(close) schedule(static) reduction(+:rmat_sum)
        do iframe=1,nframes
            if( frameweights(iframe) < 1.e-6 ) cycle
            if( params_glob%l_dose_weight )then
                ! dose weighing
                call movie_frames_shifted_saved(iframe)%fft
                call movie_frames_shifted_saved(iframe)%apply_filter_serial(filtarr(iframe,:))
                ! real space interpolation
                call movie_frames_shifted_saved(iframe)%ifft
                call motion_patch%polytransfo(iframe, movie_frames_shifted_saved(iframe), movie_frames_shifted_patched(iframe))
            endif
            ! sum
            call movie_frames_shifted_patched(iframe)%get_rmat_ptr(prmat)
            rmat_sum(:,:,1) = rmat_sum(:,:,1) + frameweights(iframe) * prmat(1:ldim_scaled(1),1:ldim_scaled(2),1)
        end do
        !$omp end parallel do
        call movie_sum_corrected%set_rmat(rmat_sum)
        ! cleanup
        do iframe=1,nframes
            call movie_frames_shifted_patched(iframe)%kill
        enddo
        deallocate(movie_frames_shifted_patched,rmat_sum)
    end subroutine motion_correct_patched_calc_sums

    ! write anisotropic shifts
    subroutine write_aniso2star
        integer :: iframe
        call starfile_table__clear(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .false.)
        call starfile_table__setName(mc_starfile, "local_motion_model")
        do iframe = 1, size(patched_polyn,1)
            call starfile_table__addObject(mc_starfile)
            call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_MOTION_COEFFS_IDX, iframe-1)
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_MOTION_COEFF, patched_polyn(iframe))
        end do
        call starfile_table__write_ofile(mc_starfile)
    end subroutine write_aniso2star

    subroutine motion_correct_patched_kill
        integer :: iframe
        call motion_patch%kill
        if( allocated(movie_frames_shifted_saved) )then
            do iframe=1,size(movie_frames_shifted_saved)
                call movie_frames_shifted_saved(iframe)%kill
            end do
            deallocate(movie_frames_shifted_saved)
        endif
        call ftexp_transfmat_kill
    end subroutine motion_correct_patched_kill


    ! PUBLIC COMMON

    subroutine close_starfile
        call starfile_table__close_ofile(mc_starfile)
        call starfile_table__delete(mc_starfile)
    end subroutine close_starfile

    subroutine motion_correct_kill_common
        if( allocated(corrs)              ) deallocate(corrs)
        if( allocated(shifts_toplot)      ) deallocate(shifts_toplot)
        if( allocated(frameweights)       ) deallocate(frameweights)
        if( allocated(acc_doses)          ) deallocate(acc_doses)
        if( allocated(cmat_sum)           ) deallocate(cmat_sum)
        call ftexp_transfmat_kill
    end subroutine motion_correct_kill_common


    ! COMMON PRIVATE UTILITY METHODS

    subroutine gen_dose_weight_filter( filtarr, movie_frames )
        real, allocatable,                          intent(inout) :: filtarr(:,:)
        type(image), allocatable, target, optional, intent(in)    :: movie_frames(:)
        type(image), pointer :: movie_frames_here(:)
        real    :: ksum_sq
        integer :: k, sz, filtsz, iframe
        if( allocated(filtarr) )deallocate(filtarr)
        if( present(movie_frames) ) then
            movie_frames_here => movie_frames
        else
            movie_frames_here => movie_frames_scaled
        end if
        sz = movie_frames_here(1)%get_filtsz()
        filtsz = 2*sz ! the filter goes well beyond nyquist so we dont'have to worry about dimensions
        allocate(filtarr(nframes,filtsz))
        do iframe=1,nframes
            filtarr(iframe,:) = acc_dose2filter(movie_frames_here(1), acc_doses(iframe), kV, filtsz)
        end do
        ! frame normalisation
        do k = 1,filtsz
            ksum_sq      = sum(frameweights * filtarr(:,k)**2.)
            filtarr(:,k) = filtarr(:,k) / sqrt(ksum_sq / sum(frameweights))
        enddo
    end subroutine gen_dose_weight_filter

end module simple_motion_correct
