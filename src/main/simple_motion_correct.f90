! motion_correct does motion correction, dose-weighting and frame-weighting of DDD movies
module simple_motion_correct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ft_expanded,                   only: ftexp_transfmat_init, ftexp_transfmat_kill
use simple_motion_patched,                only: motion_patched
use simple_motion_align_hybrid,           only: motion_align_hybrid
use simple_motion_align_iso_polyn_direct, only: motion_align_iso_polyn_direct
use simple_opt_image_weights,             only: opt_image_weights
use simple_image,                         only: image
use simple_parameters,                    only: params_glob
use simple_opt_lbfgsb,                    only: PRINT_NEVALS
implicit none

! Stage drift
public :: motion_correct_iso, motion_correct_iso_calc_sums, motion_correct_iso_calc_sums_tomo, motion_correct_iso_shift_frames
public :: motion_correct_iso_kill
! Beam-induced motion correction
public :: motion_correct_patched, motion_correct_patched_calc_sums, motion_correct_patched_kill
public :: motion_correct_with_patched
! Common & convenience
public :: motion_correct_kill_common, motion_correct_mic2spec, patched_shift_fname
public :: motion_correct_write2star, motion_correct_calc_opt_weights
private
#include "simple_local_flags.inc"

! data structures for isotropic correction
type(image), target, allocatable :: movie_frames_scaled(:) !< scaled movie frames
real,                allocatable :: opt_shifts(:,:)        !< optimal shifts identified
real                             :: TOL_ISO = 1.e-6        !< LBFGSB tolerance for isotropic search
integer                          :: updateres
logical                          :: didupdateres

! data structures for patch-based motion correction
type(motion_patched)     :: motion_patch
real(dp),    allocatable :: patched_polyn(:)              !< polynomial from patched-based correction

! data structures used by both isotropic & patch-based correction
real,        allocatable :: shifts_toplot(:,:)            !< shifts for plotting & parsing
real,        allocatable :: frameweights(:)               !< array of frameweights
real,        allocatable :: acc_doses(:)                  !< accumulated doses
complex,     allocatable :: cmat_sum(:,:,:)               !< complex matrice

! module global variables
integer :: nframes        = 0                             !< number of frames
integer :: fixed_frame    = 0                             !< fixed frame of reference for isotropic alignment (0,0)
integer :: ldim(3)        = [0,0,0]                       !< logical dimension of frame
integer :: ldim_orig(3)   = [0,0,0]                       !< logical dimension of frame (original, for use in motion_correct_iter)
integer :: ldim_scaled(3) = [0,0,0]                       !< shrunken logical dimension of frame
real    :: hp             = 0.                            !< high-pass limit
real    :: lp             = 0.                            !< low-pass limit
real    :: resstep        = 0.                            !< resolution step size (in angstrom)
real    :: smpd           = 0.                            !< sampling distance
real    :: smpd_scaled    = 0.                            !< sampling distance
real    :: kV             = 300.                          !< acceleration voltage
real    :: dose_rate      = 0.                            !< dose rate
logical :: do_scale                    = .false.          !< scale or not
logical :: motion_correct_with_patched = .false.          !< run patch-based aniso or not
character(len=:), allocatable :: patched_shift_fname      !< file name for shift plot for patched-based alignment

! module global constants
real,    parameter :: NSIGMAS                       = 6.       !< Number of standard deviations for outliers detection
logical, parameter :: FITSHIFTS                     = .true.
logical, parameter :: ISO_POLYN_DIRECT              = .false.   !< use polynomial constraint for isotropic motion correction
logical, parameter :: ISO_UNCONSTR_AFTER            = .false.   !< run a unconstrained (direct) as the second step (at highest resolution)
logical, parameter :: DO_PATCHED_POLYN              = .false.   !< run polynomially constrained motion correction for patch-based motion correction
logical, parameter :: DO_PATCHED_POLYN_DIRECT_AFTER = .false.   !< run a direct polynomial optimization for patch-based motion correction as the second step (at highest resolution)
logical, parameter :: DO_OPT_WEIGHTS                = .false.   !< continuously optimize weights after alignment

! benchmarking
logical                 :: L_BENCH = .false.
integer(timer_int_kind) :: t_read, t_cure, t_forctf, t_mic, t_fft_clip, t_patched, t_patched_forctf, t_patched_mic
integer(timer_int_kind) :: t_correct_iso_init, t_correct_iso_transfmat,  t_correct_iso_align, t_dose, t_new, t
real(timer_int_kind)    :: rt_read, rt_cure, rt_forctf, rt_mic, rt_fft_clip, rt_patched, rt_patched_forctf, rt_patched_mic
real(timer_int_kind)    :: rt_correct_iso_init, rt_correct_iso_transfmat, rt_correct_iso_align, rt_dose, rt_new

contains

    ! PUBLIC METHODS, ISOTROPIC MOTION CORRECTION

    subroutine motion_correct_init( movie_stack_fname, ctfvars, err, movie_sum, gainref_fname )
        character(len=*),           intent(in)    :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(in)    :: ctfvars           !< CTF parameters
        logical,                    intent(out)   :: err               !< error flag
        type(image),                intent(inout) :: movie_sum
        character(len=*), optional, intent(in)    :: gainref_fname     !< gain reference filename
        type(image), allocatable :: movie_frames(:)
        complex,     pointer     :: pcmat(:,:,:)
        real     :: moldiam, dimo4, time_per_frame, current_time
        integer  :: iframe, shp(3)
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
        hp      = min(dimo4,1000.)
        lp      = params_glob%lpstart
        resstep = (params_glob%lpstart - params_glob%lpstop) / 3.
        ! allocate abstract data structures
        allocate( movie_frames(nframes), movie_frames_scaled(nframes),stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 1; simple_motion_correct')
        allocate( shifts_toplot(nframes, 2), source=0., stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 2; simple_motion_correct')
        ! additional allocations
        allocate(opt_shifts(nframes,2),frameweights(nframes), source=0.,stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('motion_correct_init 4; simple_motion_correct')
        frameweights = 1./real(nframes)
        ! allocate & read frames
        if( l_BENCH ) t_new = tic()
        !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
        do iframe=1,nframes
            call movie_frames(iframe)%new(ldim, smpd, wthreads=.false.)
        enddo
        !$omp end parallel do
        if( l_BENCH ) rt_new = toc(t_new)
        if( l_BENCH ) t_read = tic()
        do iframe=1,nframes
            call movie_frames(iframe)%read(movie_stack_fname, iframe)
        end do
        if( l_BENCH ) rt_read = toc(t_read)
        ! gain reference & outliers detections
        if( l_BENCH ) t_cure = tic()
        call cure_outliers( movie_frames, gainref_fname )
        if( l_BENCH ) rt_cure = toc(t_cure)
        if( l_BENCH ) t_fft_clip = tic()
        call movie_sum%new(ldim_scaled, smpd_scaled)
        shp = movie_sum%get_array_shape()
        !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%new(ldim_scaled, smpd_scaled, wthreads=.false.)
            call movie_frames(iframe)%fft()
            call movie_frames(iframe)%clip(movie_frames_scaled(iframe))
            call movie_frames(iframe)%kill
        enddo
        !$omp end parallel do
        ! generate sum
        allocate(cmat_sum(shp(1),shp(2),shp(3)), source=cmplx(0.,0.))
        do iframe = 1,nframes
            call movie_frames_scaled(iframe)%get_cmat_ptr(pcmat)
            !$omp parallel workshare proc_bind(close)
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + pcmat(:,:,:)
            !$omp end parallel workshare
        enddo
        call movie_sum%set_cmat(cmat_sum)
        call movie_sum%ifft
        if( l_BENCH ) rt_fft_clip = toc(t_fft_clip)
        ! check if we are doing dose weighting
        if( params_glob%l_dose_weight )then
            if( allocated(acc_doses) ) deallocate(acc_doses)
            allocate( acc_doses(nframes), stat=alloc_stat )
            kV = ctfvars%kv
            time_per_frame = params_glob%exp_time/real(nframes)  ! unit: s
            dose_rate      = params_glob%dose_rate
            do iframe=1,nframes
                current_time      = real(iframe)*time_per_frame ! unit: s
                acc_doses(iframe) = dose_rate*current_time      ! unit: e/A2/s * s = e/A2
            end do
        endif
        deallocate(movie_frames)
        if( L_BENCH )then
            print *,'t_fft_clip: ',rt_fft_clip
            print *,'t_cure:     ',rt_cure
            print *,'t_new:      ',rt_new
            print *,'t_read:     ',rt_read
        endif
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
    subroutine motion_correct_iso( movie_stack_fname, ctfvars, bfactor, movie_sum, gainref_fname )
        character(len=*),           intent(in)    :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(inout) :: ctfvars           !< CTF params
        real,                       intent(in)    :: bfactor           !< B-factor
        type(image),                intent(inout) :: movie_sum
        character(len=*), optional, intent(in)    :: gainref_fname     !< gain reference filename
        real                                :: ave, sdev, var, minw, maxw, corr
        logical                             :: err, err_stat
        type(motion_align_hybrid)           :: hybrid_srch
        type(motion_align_iso_polyn_direct) :: align_iso_polyn_direct
        class(*), pointer                   :: callback_ptr
        callback_ptr => null()
        ! initialise
        if (ISO_POLYN_DIRECT) then
            call align_iso_polyn_direct%new
        end if
        if( l_BENCH ) t_correct_iso_init = tic()
        call motion_correct_init(movie_stack_fname, ctfvars, err, movie_sum, gainref_fname)
        if( l_BENCH ) rt_correct_iso_init = toc(t_correct_iso_init)
        if( err ) return
        if( l_BENCH ) t_correct_iso_transfmat = tic()
        call ftexp_transfmat_init(movie_frames_scaled(1), params_glob%lpstop)
        if( l_BENCH ) rt_correct_iso_transfmat = toc(t_correct_iso_transfmat)
        if( l_BENCH ) t_correct_iso_align = tic()
        if (ISO_POLYN_DIRECT) then
            call align_iso_polyn_direct%set_frames(movie_frames_scaled, nframes)
            call align_iso_polyn_direct%set_hp_lp(hp,lp)
            updateres = 0
            call align_iso_polyn_direct%set_callback( motion_correct_iso_polyn_direct_callback )
            call align_iso_polyn_direct%set_group_frames(trim(params_glob%groupframes).eq.'yes')
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
        end if
        if( l_BENCH ) rt_correct_iso_align = toc(t_correct_iso_align)
        ! deals with reference frame convention
        select case(trim(params_glob%mcconvention))
        case('relion')
            opt_shifts(:,1) = opt_shifts(:,1) - opt_shifts(1,1)
            opt_shifts(:,2) = opt_shifts(:,2) - opt_shifts(1,2)
            shifts_toplot   = opt_shifts
        case DEFAULT
            ! using central frame
        end select
        call moment(frameweights, ave, sdev, var, err_stat)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        write(logfhandle,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT :', ave
        write(logfhandle,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS:', sdev
        write(logfhandle,'(a,7x,f7.4)') '>>> MIN WEIGHT     :', minw
        write(logfhandle,'(a,7x,f7.4)') '>>> MAX WEIGHT     :', maxw
        ! report the sampling distance of the possibly scaled movies
        ctfvars%smpd = smpd_scaled
        if( L_BENCH )then
            print *,'rt_correct_iso_init:      ',rt_correct_iso_init
            print *,'rt_correct_iso_transfmat: ',rt_correct_iso_transfmat
            print *,'rt_correct_iso_align:     ',rt_correct_iso_align
        endif
    end subroutine motion_correct_iso

    ! generates sums & shift frames
    subroutine motion_correct_iso_shift_frames()
        integer :: iframe
        ! copy unaligned frames & shift frames
        if( l_BENCH ) t = tic()
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%shift2Dserial(-opt_shifts(iframe,:))
        end do
        !$omp end parallel do
        if( L_BENCH ) print *,'t_shift:      ', toc(t)
    end subroutine motion_correct_iso_shift_frames

    ! Optimal weights, frames assumed in Fourier space
    subroutine motion_correct_calc_opt_weights()
        type(opt_image_weights)    :: opt_weights
        if (DO_OPT_WEIGHTS) then
            if( l_BENCH ) t = tic()
            call opt_weights%new(movie_frames_scaled, hp, params_glob%lpstop)
            call opt_weights%calc_opt_weights
            frameweights = opt_weights%get_weights()
            call opt_weights%kill
            if( L_BENCH ) print *,'t_opt_weights:',toc(t)
        end if
    end subroutine motion_correct_calc_opt_weights

    !>  Generates sums, movie_frames_scaled are assumed shifted
    subroutine motion_correct_iso_calc_sums( movie_sum_corrected, movie_sum_ctf)
        type(image), intent(inout) :: movie_sum_corrected, movie_sum_ctf
        complex,           pointer :: pcmat(:,:,:)
        real                       :: scalar_weight
        integer                    :: iframe
        ! for CTF estimation
        if( l_BENCH ) t_forctf = tic()
        call movie_sum_ctf%new(ldim_scaled, smpd_scaled)
        scalar_weight = 1. / real(nframes)
        ! in case anisotropic alignement failed
        if( l_BENCH ) t_forctf = tic()
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%fft
        enddo
        !$omp end parallel do
        ! Sum for CTF estimation
        cmat_sum = cmplx(0.,0.)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%get_cmat_ptr(pcmat)
            !$omp parallel workshare proc_bind(close)
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + scalar_weight*pcmat(:,:,:)
            !$omp end parallel workshare
        end do
        call movie_sum_ctf%set_cmat(cmat_sum)
        call movie_sum_ctf%ifft
        if( l_BENCH ) rt_forctf = toc(t_forctf)
        ! dose-weighting
        if( l_BENCH ) t_dose = tic()
        call apply_dose_weighting(movie_frames_scaled(:))
        if( l_BENCH ) rt_dose = toc(t_dose)
        ! Micrograph
        if( l_BENCH ) t_mic = tic()
        call movie_sum_corrected%new(ldim_scaled, smpd_scaled)
        cmat_sum = cmplx(0.,0.)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%get_cmat_ptr(pcmat)
            !$omp parallel workshare proc_bind(close)
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + frameweights(iframe)*pcmat(:,:,:)
            !$omp end parallel workshare
        end do
        call movie_sum_corrected%set_cmat(cmat_sum)
        call movie_sum_corrected%ifft
        if( l_BENCH ) rt_mic = toc(t_mic)
        if( L_BENCH )then
            print *,'t_forctf:     ',rt_forctf
            print *,'t_dose:       ',rt_dose
            print *,'t_mic:        ',rt_mic
        endif
    end subroutine motion_correct_iso_calc_sums

    subroutine motion_correct_iso_calc_sums_tomo( frame_counter, time_per_frame, movie_sum_corrected, movie_sum_ctf )
        integer,     intent(inout) :: frame_counter  !< frame counter
        real,        intent(in)    :: time_per_frame !< time resolution
        type(image), intent(inout) :: movie_sum_corrected, movie_sum_ctf
        real    :: current_time
        integer :: iframe
        ! re-calculates doses
        do iframe=1,nframes
            frame_counter     = frame_counter + 1
            current_time      = real(frame_counter)*time_per_frame ! unit: seconds
            acc_doses(iframe) = dose_rate*current_time             ! unit e/A2
        end do
        call  motion_correct_iso_calc_sums(movie_sum_corrected, movie_sum_ctf)
    end subroutine motion_correct_iso_calc_sums_tomo

    ! Write iso/aniso-tropic shifts
    subroutine motion_correct_write2star( mc_starfile_fname, moviename, writepoly, gainref_fname )
        use simple_starfile_wrappers
        character(len=*),           intent(in) :: mc_starfile_fname, moviename
        logical,                    intent(in) :: writepoly
        character(len=*), optional, intent(in) :: gainref_fname
        type(starfile_table_type) :: mc_starfile
        real(dp) :: u,v,poly_coeffs(size(patched_polyn,1))
        real     :: shift(2), doserateperframe
        integer  :: iframe, n
        call starfile_table__new(mc_starfile)
        call starfile_table__open_ofile(mc_starfile, mc_starfile_fname)
        ! global fields
        call starfile_table__addObject(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .true.)
        call starfile_table__setname(mc_starfile, "general")
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_X, ldim_orig(1))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Y, ldim_orig(2))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Z, nframes)
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
        if( writepoly )then
            call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, 1)
        else
            call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, 0)
        endif
        call starfile_table__write_ofile(mc_starfile)
        ! isotropic shifts
        call starfile_table__clear(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .false.)
        call starfile_table__setName(mc_starfile, "global_shift")
        do iframe = 1, nframes
            shift = shifts_toplot(iframe,:) - shifts_toplot(1,:)
            call starfile_table__addObject(mc_starfile)
            call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, real(shift(1)/params_glob%scale, dp))
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, real(shift(2)/params_glob%scale, dp))
            ! call starfile_table__setValue_double(mc_starfile, EMDL_MOVIE_FRAME_WEIGHT, real(frameweights(iframe), dp))
        end do
        call starfile_table__write_ofile(mc_starfile)
        if( writepoly )then
            ! anisotropic shifts
            n = size(patched_polyn,1)
            poly_coeffs = patched_polyn
            if( do_scale )then
                u = 1.d0 / real(params_glob%scale,dp)
                v = real(params_glob%scale,dp)
                poly_coeffs(    1:n/2) = patched_polyn(    1:n/2) * [u,u,u,1.d0,1.d0,1.d0,v,v,v,1.d0,1.d0,1.d0,v,v,v,v,v,v]
                poly_coeffs(n/2+1:n)   = patched_polyn(n/2+1:n)   * [u,u,u,1.d0,1.d0,1.d0,v,v,v,1.d0,1.d0,1.d0,v,v,v,v,v,v]
            endif
            call starfile_table__clear(mc_starfile)
            call starfile_table__setIsList(mc_starfile, .false.)
            call starfile_table__setName(mc_starfile, "local_motion_model")
            do iframe = 1, n
                call starfile_table__addObject(mc_starfile)
                call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_MOTION_COEFFS_IDX, iframe-1)
                call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_MOTION_COEFF, poly_coeffs(iframe))
            end do
            call starfile_table__write_ofile(mc_starfile)
        endif
        call starfile_table__close_ofile(mc_starfile)
        call starfile_table__delete(mc_starfile)
    end subroutine motion_correct_write2star

    subroutine motion_correct_iso_kill
        if (allocated(opt_shifts)) deallocate(opt_shifts)
        call ftexp_transfmat_kill
    end subroutine motion_correct_iso_kill

    ! PUBLIC METHODS, PATCH-BASED MOTION CORRECTION

    !>  patch-based motion_correction of DDD movie
    !>  movie_frames_scaled assumed shifted, in Fourier domain
    subroutine motion_correct_patched( bfac, rmsd )
        real, intent(in)  :: bfac
        real, intent(out) :: rmsd(2)     !< whether polynomial fitting was within threshold
        integer :: iframe
        if( l_BENCH ) t_patched = tic()
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%ifft
        end do
        !$omp end parallel do
        rmsd = huge(rmsd(1))
        call motion_patch%new(motion_correct_ftol = params_glob%motion_correctftol, &
            motion_correct_gtol = params_glob%motion_correctgtol, trs = params_glob%scale*params_glob%trs)
        write(logfhandle,'(A,I2,A3,I2,A1)') '>>> PATCH-BASED REFINEMENT (',&
            &params_glob%nxpatch,' x ',params_glob%nypatch,')'
        PRINT_NEVALS = .false.
        if (DO_PATCHED_POLYN) then
            call motion_patch%set_bfactor(bfac)
            call motion_patch%correct_polyn( hp, resstep, movie_frames_scaled,&
                &patched_shift_fname, DO_PATCHED_POLYN_DIRECT_AFTER, shifts_toplot)
            rmsd = 0.
        else
            call motion_patch%set_fitshifts(FITSHIFTS)
            call motion_patch%set_frameweights(frameweights)
            call motion_patch%set_fixed_frame(fixed_frame)
            call motion_patch%set_interp_fixed_frame(fixed_frame)
            call motion_patch%set_bfactor(bfac)
            call motion_patch%correct(hp, resstep, movie_frames_scaled, patched_shift_fname, shifts_toplot)
            rmsd = motion_patch%get_polyfit_rmsd()
        end if
        call motion_patch%get_poly4star(patched_polyn)
        if( L_BENCH )then
            rt_patched = toc(t_patched)
            print *,'rt_patched:      ',rt_patched
        endif
    end subroutine motion_correct_patched

    ! weighted & un-weighted sums, dose-weighting
    ! movie_frames_shifted assumed in real-space
    subroutine motion_correct_patched_calc_sums( movie_sum_corrected, movie_sum_ctf )
        type(image), intent(inout) :: movie_sum_corrected, movie_sum_ctf
        real                       :: scalar_weight, weights(nframes)
        integer                    :: iframe
        scalar_weight = 1. / real(nframes)
        ! for CTF estimation
        if( l_BENCH ) t_patched_forctf = tic()
        call movie_sum_ctf%new(ldim_scaled, smpd_scaled)
        weights = scalar_weight
        call motion_patch%polytransfo(movie_frames_scaled(:), weights, movie_sum_ctf)
        if( l_BENCH ) rt_patched_forctf = toc(t_patched_forctf)
        if( params_glob%l_dose_weight )then
            ! dose weighing
            if( l_BENCH ) t_dose = tic()
            !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(guided)
            do iframe=1,nframes
                call movie_frames_scaled(iframe)%fft
            end do
            !$omp end parallel do
            call apply_dose_weighting(movie_frames_scaled(:))
            !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(guided)
            do iframe=1,nframes
                call movie_frames_scaled(iframe)%ifft
            end do
            !$omp end parallel do
            if( l_BENCH ) rt_dose = toc(t_dose)
        endif
        ! micrograph
        if( l_BENCH ) t_patched_mic = tic()
        call movie_sum_corrected%new(ldim_scaled, smpd_scaled)
        weights = frameweights
        call motion_patch%polytransfo(movie_frames_scaled(:), weights, movie_sum_corrected)
        if( l_BENCH ) rt_patched_mic = toc(t_patched_mic)
        if( L_BENCH )then
            rt_patched = toc(t_patched)
            print *,'rt_patched_forctf: ',rt_patched_forctf
            print *,'rt_patched_mic   : ',rt_patched_mic
            print *,'rt_dose   : ',rt_dose
        endif
    end subroutine motion_correct_patched_calc_sums

    subroutine motion_correct_patched_kill
        call motion_patch%kill
        call ftexp_transfmat_kill
    end subroutine motion_correct_patched_kill

    ! PUBLIC COMMON

    !> mic2spec calculates the average powerspectrum over a micrograph
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine motion_correct_mic2spec( img_in, box, speckind, lp_backgr_subtr, img_out )
        class(image),      intent(inout) :: img_in
        integer,           intent(in)    :: box
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        type(image)       :: tiles(nthr_glob), tmp(nthr_glob)
        real, allocatable :: rmat_sum(:,:,:)
        real, pointer     :: prmat(:,:,:)
        real    :: smpd
        integer :: ldim(3), i,j, n, hbox, ithr
        logical :: outside
        hbox = box/2
        call img_in%ifft
        smpd = img_in%get_smpd()
        ldim = img_in%get_ldim()
        n    = 0
        call img_out%new([box,box,1], smpd)
        call img_out%zero_and_unflag_ft
        rmat_sum = img_out%get_rmat()
        ! sum individual spectra
        !$omp parallel do private(i,j,ithr,prmat) default(shared) schedule(static)&
        !$omp proc_bind(close) reduction(+:n,rmat_sum)
        do j = 0,ldim(2)-box,hbox
            do i = 0,ldim(1)-box,hbox
                n    = n+1
                ithr = omp_get_thread_num() + 1
                if( tiles(ithr)%exists() )then
                    call tiles(ithr)%zero_and_unflag_ft
                else
                    call tiles(ithr)%new([box,box,1], smpd, wthreads=.false.)
                    call tmp(ithr)%new([box,box,1], smpd, wthreads=.false.)
                endif
                call img_in%window_slim([i,j], box, tiles(ithr), outside)
                call tiles(ithr)%norm()
                call tiles(ithr)%zero_edgeavg
                call tiles(ithr)%fft()
                call tiles(ithr)%ft2img(speckind, tmp(ithr))
                call tmp(ithr)%get_rmat_ptr(prmat)
                rmat_sum(:box,:box,1) = rmat_sum(:box,:box,1) + prmat(:box,:box,1)
            enddo
        enddo
        !$omp end parallel do
        call img_out%set_rmat(rmat_sum)
        ! average
        call img_out%div(real(n))
        ! post-process
        call img_out%dampen_pspec_central_cross
        call img_out%subtr_backgr(lp_backgr_subtr)
        ! cleanup
        do ithr=1,nthr_glob
            call tiles(ithr)%kill
            call tmp(ithr)%kill
        enddo
    end subroutine motion_correct_mic2spec

    subroutine motion_correct_kill_common
        integer :: iframe
        if( allocated(shifts_toplot)      ) deallocate(shifts_toplot)
        if( allocated(frameweights)       ) deallocate(frameweights)
        if( allocated(acc_doses)          ) deallocate(acc_doses)
        if( allocated(cmat_sum)           ) deallocate(cmat_sum)
        if( allocated(movie_frames_scaled) )then
            do iframe=1,size(movie_frames_scaled)
                call movie_frames_scaled(iframe)%kill
            end do
            deallocate(movie_frames_scaled)
        endif
        call ftexp_transfmat_kill
    end subroutine motion_correct_kill_common

    ! COMMON PRIVATE UTILITY METHODS

    ! Following Grant & Grigorieff; eLife 2015;4:e06980
    ! Frames assumed in fourier space
    subroutine apply_dose_weighting( frames )
        class(image), intent(inout) :: frames(nframes)
        real, parameter   :: A=0.245, B=-1.665, C=2.81
        real, allocatable :: qs(:)
        real    :: frame_dose(nframes), spaFreqk
        real    :: twoNe, smpd, spafreq, limhsq,limksq
        integer :: nrflims(3,2), ldim(3), hphys,kphys, iframe, h,k
        if( .not.params_glob%l_dose_weight ) return
        if( .not.frames(1)%is_ft() ) THROW_HARD('Frames should be in in the Fourier domain')
        nrflims = frames(1)%loop_lims(2)
        smpd    = frames(1)%get_smpd()
        allocate(qs(nframes),source=0.)
        ! doses
        ldim   = frames(1)%get_ldim()
        limhsq = (real(ldim(1))*smpd)**2.
        limksq = (real(ldim(2))*smpd)**2.
        do iframe = 1,nframes
            frame_dose(iframe) = acc_doses(iframe)
            if( is_equal(kV,200.) )then
                frame_dose(iframe) = frame_dose(iframe) / 0.8
            else if( is_equal(kV,100.) )then
                frame_dose(iframe) = frame_dose(iframe) / 0.64
            endif
        enddo
        ! dose normalization
        !$omp parallel do private(h,k,spafreq,spafreqk,twone,kphys,hphys,iframe,qs)&
        !$omp default(shared) schedule(static) proc_bind(close)
        do k = nrflims(2,1),nrflims(2,2)
            kphys    = k + 1 + merge(ldim(2),0,k<0)
            spaFreqk = real(k*k)/limksq
            do h = nrflims(1,1),nrflims(1,2)
                hphys   = h + 1
                spaFreq = sqrt( real(h*h)/limhsq + spaFreqk )
                twoNe   = 2.*(A*spaFreq**B + C)
                qs = exp(-frame_dose/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                do iframe = 1,nframes
                    call frames(iframe)%mul_cmat_at([hphys,kphys,1], qs(iframe))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_dose_weighting

    ! gain correction, calculate image sum and identify outliers
    subroutine cure_outliers( frames, gainref_fname )
        class(image),               intent(inout) :: frames(nframes)
        character(len=*), optional, intent(in)    :: gainref_fname
        integer, parameter   :: hwinsz = 6
        type(image)          :: gainref
        real,        pointer :: prmat(:,:,:)
        real,    allocatable :: rsum(:,:), vals(:), new_vals(:,:)
        integer, allocatable :: pos_outliers(:,:)
        real    :: ave, sdev, var, lthresh, uthresh
        integer :: iframe, noutliers, l,r,b,u,i,j,k, nvals, winsz
        logical :: outliers(ldim(1),ldim(2)), err
        allocate(rsum(ldim(1),ldim(2)),source=0.)
        if( l_BENCH ) t_cure = tic()
        write(logfhandle,'(a)') '>>> REMOVING DEAD/HOT PIXELS & FOURIER TRANSFORMING FRAMES'
        if( present(gainref_fname) )then
            if( file_exists(gainref_fname) )then
                call gainref%new(ldim, smpd)
                call gainref%read(gainref_fname)
                !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
                do iframe = 1,nframes
                    call frames(iframe)%mul(gainref)
                enddo
                !$omp end parallel do
                call gainref%kill
            else
                THROW_HARD('gain reference: '//trim(gainref_fname)//' not found; motion_correct_init')
            endif
        endif
        ! sum
        do iframe = 1,nframes
            call frames(iframe)%get_rmat_ptr(prmat)
            !$omp parallel workshare
            rsum(:,:) = rsum(:,:) + prmat(:ldim(1),:ldim(2),1)
            !$omp end parallel workshare
        enddo
        ! outliers detection
        call moment( rsum, ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh   = ave - NSIGMAS * sdev
        uthresh   = ave + NSIGMAS * sdev
        where( rsum<lthresh .or. rsum>uthresh )
            outliers = .true.
        elsewhere
            outliers = .false.
        end where
        noutliers = count(outliers)
        deallocate(rsum)
        ! cure
        if( noutliers > 0 )then
            write(logfhandle,'(a,1x,i7)') '>>> # DEAD/HOT PIXELS:', noutliers
            winsz = 2*HWINSZ+1
            nvals = winsz*winsz
            allocate(new_vals(noutliers,nframes),pos_outliers(2,noutliers),vals(nvals))
            ! gather positions
            k = 0
            do j = 1,ldim(2)
                do i = 1,ldim(1)
                    if( outliers(i,j) )then
                        k = k + 1
                        pos_outliers(:,k) = [i,j]
                    endif
                enddo
            enddo
            !$omp parallel do default(shared) private(iframe,k,i,j,l,r,b,u,vals,prmat)&
            !$omp proc_bind(close) schedule(static)
            do iframe=1,nframes
                call frames(iframe)%get_rmat_ptr(prmat)
                ! calulate new values
                do k = 1,noutliers
                    i = pos_outliers(1,k)
                    j = pos_outliers(2,k)
                    l = max(i-HWINSZ,1)
                    r = min(ldim(1),l+winsz-1)
                    l = r-winsz+1
                    b = max(j-HWINSZ,1)
                    u = min(ldim(2),b+winsz-1)
                    b = u-winsz+1
                    vals = reshape(prmat(l:r,b:u,1),(/nvals/))
                    new_vals(k,iframe) = median_nocopy(vals)
                enddo
                ! substitute
                do k = 1,noutliers
                    i = pos_outliers(1,k)
                    j = pos_outliers(2,k)
                    prmat(i,j,1) = new_vals(k,iframe)
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine cure_outliers

end module simple_motion_correct
