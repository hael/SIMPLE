! motion_correct does motion correction, dose-weighting and frame-weighting of DDD movies
module simple_motion_correct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ft_expanded,                   only: ftexp_transfmat_init, ftexp_transfmat_kill
use simple_motion_patched,                only: motion_patched
use simple_motion_align_hybrid,           only: motion_align_hybrid
use simple_opt_image_weights,             only: opt_image_weights
use simple_image,                         only: image
use simple_eer_factory,                   only: eer_decoder
use simple_parameters,                    only: params_glob
implicit none

! Stage drift
public :: motion_correct_iso, motion_correct_iso_calc_sums, motion_correct_iso_shift_frames
public :: motion_correct_iso_kill
! Beam-induced motion correction
public :: motion_correct_dev, motion_correct_patched, motion_correct_patched_calc_sums, motion_correct_patched_kill
public :: motion_correct_with_patched
! Common & convenience
public :: motion_correct_kill_common, motion_correct_mic2spec, patched_shift_fname
public :: motion_correct_write2star, motion_correct_calc_opt_weights, motion_correct_calc_msd
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
real(dp),    allocatable :: patched_shifts(:,:,:,:)       !< shifts from patched-based correction
real(dp),    allocatable :: patched_centers(:,:,:)        !< patch centers

! data structures used by both isotropic & patch-based correction
type(eer_decoder)        :: eer                           !< eer object
type(image)              :: gain_img                      !< gain reference image
real,        allocatable :: shifts_toplot(:,:)            !< shifts for plotting & parsing
real,        allocatable :: frameweights(:)               !< array of frameweights
complex,     allocatable :: cmat_sum(:,:,:)               !< complex matrice

! module global variables
character(len=:), allocatable :: patched_shift_fname      !< file name for shift plot for patched-based alignment
integer,          allocatable :: pos_outliers(:,:)        !< positions of defects & hot pixels
integer :: nframes        = 0                             !< number of frames used for alignement
integer :: total_nframes  = 0                             !< total number of frames in movie
integer :: fixed_frame    = 0                             !< fixed frame of reference for isotropic alignment (0,0)
integer :: ldim(3)        = [0,0,0]                       !< logical dimension of frame
integer :: ldim_orig(3)   = [0,0,0]                       !< logical dimension of frame (original, for use in motion_correct_iter)
integer :: ldim_scaled(3) = [0,0,0]                       !< shrunken logical dimension of frame
integer :: ldim_gain(3)   = [0,0,0]                       !< gain reference dimensions
integer :: eer_fraction   = 0
real    :: total_dose     = 0.                            !< total dose in e/A2
real    :: hp             = 0.                            !< high-pass limit
real    :: lp             = 0.                            !< low-pass limit
real    :: resstep        = 0.                            !< resolution step size (in angstrom)
real    :: smpd           = 0.                            !< sampling distance
real    :: smpd_scaled    = 0.                            !< sampling distance
real    :: kV             = 300.                          !< acceleration voltage
logical :: do_scale                    = .false.          !< scale or not
logical :: l_eer                       = .false.          !< whether eer format is in use
logical :: motion_correct_with_patched = .false.          !< run patch-based aniso or not


! module global constants
real,    parameter :: NSIGMAS          = 6.       !< Number of standard deviations for outliers detection
logical, parameter :: FITSHIFTS        = .true.
logical, parameter :: DO_OPT_WEIGHTS   = .false.  !< continuously optimize weights after alignment
! benchmarking
logical                 :: L_BENCH = .false.
integer(timer_int_kind) :: t_read, t_cure, t_forctf, t_mic, t_fft_clip, t_patched, t_patched_forctf, t_patched_mic
integer(timer_int_kind) :: t_correct_iso_init, t_correct_iso_transfmat,  t_correct_iso_align, t_dose, t_new, t_, t_mic2spec
real(timer_int_kind)    :: rt_read, rt_cure, rt_forctf, rt_mic, rt_fft_clip, rt_patched, rt_patched_forctf, rt_patched_mic, rt_aniso_align
real(timer_int_kind)    :: rt_correct_iso_init, rt_correct_iso_transfmat, rt_correct_iso_align, rt_dose, rt_new, rt_align_poly

contains

    ! PUBLIC METHODS, ISOTROPIC MOTION CORRECTION

    subroutine motion_correct_init( movie_stack_fname, ctfvars, err, movie_sum, gainref )
        character(len=*),           intent(in)    :: movie_stack_fname !< input filename of stack
        type(ctfparams),            intent(in)    :: ctfvars           !< CTF parameters
        logical,                    intent(out)   :: err               !< error flag
        type(image),                intent(inout) :: movie_sum
        character(len=*), optional, intent(in)    :: gainref           !< gain reference filename
        type(image), allocatable :: movie_frames(:)
        complex,     pointer     :: pcmat(:,:,:)
        real     :: dimo4, dose_per_frame
        integer  :: shp(3), iframe, n_eer_frames
        smpd       = ctfvars%smpd ! un-scaled pixel size
        total_dose = params_glob%total_dose
        ! get number of frames & dim from stack
        select case(fname2format(movie_stack_fname))
        case('K')
            ! EER case, we need to determine fractionation
            l_eer = .true.
            call find_ldim_nptcls(movie_stack_fname, ldim, n_eer_frames)
            if( params_glob%l_eer_fraction )then
                eer_fraction = params_glob%eer_fraction
                nframes      = floor(real(n_eer_frames)/real(eer_fraction))
                ! because we ignore the last eer frames we need to update the dose
                total_dose   = total_dose * real(nframes) * real(eer_fraction) / real(n_eer_frames)
            else
                if( .not.params_glob%l_dose_weight )then
                    THROW_HARD('Either eer_fraction or total_dose must be defined; motion_correct_init')
                endif
                call calc_eer_fraction(n_eer_frames, params_glob%fraction_dose_target, total_dose, nframes, eer_fraction)
            endif
            write(logfhandle,'(A,2I6)')'>>> NUMBER OF FRACTIONS; EER FRAMES PER FRACTION: ', nframes, eer_fraction
            if( params_glob%l_dose_weight )write(logfhandle,'(A,F8.2)')'>>> EFFECTIVE TOTAL DOSE:', total_dose
            select case(params_glob%eer_upsampling)
            case(1)
                ! 4K, physical pixel size is by convention set on import
            case(2)
                ldim    = 2 * ldim
                ldim(3) = 1
                smpd    = smpd / 2. !! 8K upsampling
            case DEFAULT
                THROW_HARD('Unsupported up-sampling: '//int2str(params_glob%eer_upsampling)//'; motion_correct_init')
            end select
        case DEFAULT
            call find_ldim_nptcls(movie_stack_fname, ldim, nframes)
            l_eer = .false.
        end select
        total_nframes = nframes
        ldim_orig     = ldim
        err = .false.
        if( nframes < 2 )then
            err = .true.
            write(logfhandle,*) 'movie: ', trim(movie_stack_fname)
            THROW_WARN('nframes of movie < 2, aborting motion_correct')
            return
        endif
        ldim(3) = 1
        ! dose weighting prep
        if( params_glob%l_dose_weight )then
            kV = ctfvars%kv
            dose_per_frame = total_dose / real(nframes)
            if( params_glob%max_dose > 0.001 )then
                ! adjusting number of frames for maximum dose
                total_dose = 0.0
                do iframe = 1,nframes
                    total_dose  = total_dose + dose_per_frame ! e/A2
                    if( total_dose > params_glob%max_dose ) exit
                enddo
                nframes = min(iframe,nframes)
            endif
            if( nframes < 2 )then
                err = .true.
                write(logfhandle,*) 'movie: ', trim(movie_stack_fname)
                THROW_WARN('nframes of movie < 2, aborting motion_correct')
                return
            endif
        endif
        ! scaling
        if( params_glob%scale < 0.99 )then
            ldim_scaled(1) = round2even(real(ldim(1))*params_glob%scale)
            ldim_scaled(2) = round2even(real(ldim(2))*params_glob%scale)
            ldim_scaled(3) = 1
            do_scale       = .true.
        else
            ldim_scaled = ldim
            do_scale    = .false.
        endif
        ! set sampling distance
        smpd_scaled = smpd/params_glob%scale
        ! set fixed frame to central one
        fixed_frame = nint(0.5*real(nframes))
        ! set reslims
        dimo4   = (real(minval(ldim_scaled(1:2))) * smpd_scaled) / 4.
        hp      = min(dimo4,1000.)
        lp      = params_glob%lpstart
        resstep = (params_glob%lpstart - params_glob%lpstop) / 3.
        ! allocate abstract data structures
        allocate( movie_frames(nframes), movie_frames_scaled(nframes) )
        allocate( shifts_toplot(nframes, 2), source=0. )
        ! additional allocations
        allocate(opt_shifts(nframes,2),frameweights(nframes), source=0. )
        frameweights = 1./real(nframes)
        ! check gain reference existence
        if( present(gainref) )then
            if( .not.file_exists(gainref) )then
                THROW_HARD('gain reference: '//trim(gainref)//' not found; motion_correct_init')
            endif
        endif
        if( L_BENCH )then
            rt_read = 0.
            rt_cure = 0.
            rt_forctf = 0.
            rt_mic = 0.
            rt_fft_clip = 0.
            rt_patched = 0.
            rt_patched_forctf = 0.
            rt_patched_mic = 0.
            rt_aniso_align = 0.
            rt_correct_iso_init = 0.
            rt_correct_iso_transfmat = 0.
            rt_correct_iso_align = 0.
            rt_dose = 0.
            rt_new = 0.
            rt_align_poly = 0.
        endif
        ! allocate & read frames
        if( l_eer )then
            if( l_BENCH ) t_read = tic()
            call eer%new(movie_stack_fname, ctfvars%smpd, params_glob%eer_upsampling) ! 4K, default
            if( l_BENCH ) rt_read = toc(t_read)
            if( l_BENCH ) t_new = tic()
            call eer%decode(movie_frames, eer_fraction)
            if( l_BENCH ) rt_new = toc(t_new)
        else
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
        endif
        ! gain reference & outliers detections
        if( l_BENCH ) t_cure = tic()
        if( present(gainref) ) call correct_gain(movie_frames, gainref)
        call cure_outliers( movie_frames)
        call gain_img%kill
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
        ! generate sum for display
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
        deallocate(movie_frames)
        if( L_BENCH )then
            print *,'t_fft_clip: ',rt_fft_clip
            print *,'t_cure:     ',rt_cure
            print *,'t_new:      ',rt_new
            print *,'t_read:     ',rt_read
        endif
    end subroutine motion_correct_init

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
        ! type(motion_align_iso_polyn_direct) :: align_iso_polyn_direct
        class(*), pointer                   :: callback_ptr
        callback_ptr => null()
        ! initialise
        if( l_BENCH ) t_correct_iso_init = tic()
        call motion_correct_init(movie_stack_fname, ctfvars, err, movie_sum, gainref_fname)
        if( l_BENCH ) rt_correct_iso_init = toc(t_correct_iso_init)
        if( err ) return
        if( l_BENCH ) t_correct_iso_transfmat = tic()
        call ftexp_transfmat_init(movie_frames_scaled(1), params_glob%lpstop)
        if( l_BENCH ) rt_correct_iso_transfmat = toc(t_correct_iso_transfmat)
        if( l_BENCH ) t_correct_iso_align = tic()
        call hybrid_srch%new(movie_frames_scaled)
        call hybrid_srch%set_group_frames(.false.)
        call hybrid_srch%set_reslims(hp, params_glob%lpstart, params_glob%lpstop)
        call hybrid_srch%set_bfactor(bfactor)
        call hybrid_srch%set_trs(params_glob%scale*params_glob%trs)
        call hybrid_srch%set_rand_init_shifts(.true.)
        call hybrid_srch%set_fitshifts(FITSHIFTS)
        call hybrid_srch%set_fixed_frame(fixed_frame)
        call hybrid_srch%align
        call hybrid_srch%get_weights(frameweights)
        corr = hybrid_srch%get_corr()
        call hybrid_srch%get_opt_shifts(opt_shifts)
        call hybrid_srch%get_shifts_toplot(shifts_toplot)
        call hybrid_srch%kill
        ! end if
        if( corr < 0. )then
           write(logfhandle,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
           THROW_WARN('OPTIMAL CORRELATION < 0.0')
        endif
        if( l_BENCH ) rt_correct_iso_align = toc(t_correct_iso_align)
        ! deals with reference frame convention
        select case(trim(params_glob%mcconvention))
        case('relion','first')
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

    ! shift fft-ed frames
    subroutine motion_correct_iso_shift_frames()
        integer :: iframe
        if( l_BENCH ) t_ = tic()
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%shift2Dserial(-opt_shifts(iframe,:))
        end do
        !$omp end parallel do
        if( L_BENCH ) print *,'t_shift:      ', toc(t_)
    end subroutine motion_correct_iso_shift_frames

    ! Optimal weights, frames assumed in Fourier space
    subroutine motion_correct_calc_opt_weights()
        type(opt_image_weights)    :: opt_weights
        if (DO_OPT_WEIGHTS) then
            if( l_BENCH ) t_ = tic()
            call opt_weights%new(movie_frames_scaled, hp, params_glob%lpstop)
            call opt_weights%calc_opt_weights
            frameweights = opt_weights%get_weights()
            call opt_weights%kill
            if( L_BENCH ) print *,'t_opt_weights:',toc(t_)
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

    subroutine motion_correct_calc_msd( include_patch, msd )
        logical, intent(in)  :: include_patch
        real,    intent(out) :: msd
        real    :: shifts(nframes,2)
        integer :: t
        if( include_patch )then
            do t = 1, nframes
                shifts(t,1) = shifts_toplot(t,1) + real(sum(patched_shifts(1,t,:,:))/real(params_glob%nxpatch*params_glob%nypatch,dp))
                shifts(t,2) = shifts_toplot(t,2) + real(sum(patched_shifts(2,t,:,:))/real(params_glob%nxpatch*params_glob%nypatch,dp))
            enddo
        else
            shifts = shifts_toplot
        endif
        msd = sum( (shifts(:,1)-shifts(1,1))**2 )
        msd = msd + sum( (shifts(:,2)-shifts(1,2))**2 )
        msd = msd/real(nframes)
    end subroutine motion_correct_calc_msd

    ! Write iso/aniso-tropic shifts
    subroutine motion_correct_write2star( mc_starfile_fname, moviename, writepoly, gainref_fname )
        use simple_starfile_wrappers
        character(len=*),           intent(in) :: mc_starfile_fname, moviename
        logical,                    intent(in) :: writepoly
        character(len=*), optional, intent(in) :: gainref_fname
        type(starfile_table_type) :: mc_starfile
        real(dp) :: poly_coeffs(size(patched_polyn,1)),dpscale
        real     :: shift(2), doseperframe
        integer  :: i,iframe, n, ndeadpixels, motion_model
        motion_model = 0
        if( writepoly ) motion_model = 1
        dpscale = real(params_glob%scale,dp)
        call starfile_table__new(mc_starfile)
        call starfile_table__open_ofile(mc_starfile, mc_starfile_fname)
        ! global fields
        call starfile_table__addObject(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .true.)
        call starfile_table__setname(mc_starfile, "general")
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_X, ldim_orig(1))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Y, ldim_orig(2))
        call starfile_table__setValue_int(mc_starfile,    EMDL_IMAGE_SIZE_Z, total_nframes)
        call starfile_table__setValue_string(mc_starfile, EMDL_MICROGRAPH_MOVIE_NAME, simple_abspath(moviename))
        if (present(gainref_fname)) then
            call starfile_table__setValue_string(mc_starfile, EMDL_MICROGRAPH_GAIN_NAME, trim(gainref_fname))
        end if
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_BINNING, 1.d0/dpscale)
        if( l_eer )then
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, real(eer%get_smpd_out(),dp))
        else
            call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, real(smpd,dp))
        endif
        doseperframe = 0.
        if( params_glob%l_dose_weight ) doseperframe = params_glob%total_dose / real(total_nframes)
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_DOSE_RATE, real(doseperframe, dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_PRE_EXPOSURE, 0.0_dp)
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_VOLTAGE, real(params_glob%kv, dp))
        call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_START_FRAME, 1)
        if( l_eer )then
            call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_EER_UPSAMPLING, params_glob%eer_upsampling)
            call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_EER_GROUPING, eer_fraction)
        endif
        call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, motion_model)
        call starfile_table__write_ofile(mc_starfile)
        ! isotropic shifts
        call starfile_table__clear(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .false.)
        call starfile_table__setName(mc_starfile, "global_shift")
        do iframe = 1, total_nframes
            call starfile_table__addObject(mc_starfile)
            call starfile_table__setValue_int(mc_starfile,    EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
            if( iframe <= nframes )then
                shift = shifts_toplot(iframe,:) - shifts_toplot(1,:)
                call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, real(shift(1),dp)/dpscale)
                call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, real(shift(2),dp)/dpscale)
            else
                call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, -9999.0d0)
                call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, -9999.0d0)
            endif
        enddo
        call starfile_table__write_ofile(mc_starfile)
        if( writepoly )then
            ! anisotropic shifts
            n = size(patched_polyn,1)
            poly_coeffs = patched_polyn
            if( do_scale )then
                poly_coeffs = patched_polyn / dpscale
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
        ! Defects & hot pixels
        if( allocated(pos_outliers) )then
            call starfile_table__clear(mc_starfile)
            call starfile_table__setIsList(mc_starfile, .false.)
            call starfile_table__setName(mc_starfile, "hot_pixels")
            ndeadpixels = size(pos_outliers,dim=2)
            do i = 1, ndeadpixels
                call starfile_table__addObject(mc_starfile)
                call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_X, real(pos_outliers(1,i)-1,dp))
                call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_y, real(pos_outliers(2,i)-1,dp))
            end do
            call starfile_table__write_ofile(mc_starfile)
        endif
        !! Patches shifts,  Unused for now
        ! if( writepoly )then
        !     call starfile_table__clear(mc_starfile)
        !     call starfile_table__setIsList(mc_starfile, .false.)
        !     call starfile_table__setName(mc_starfile, "local_shift")
        !     do i = 1, params_glob%nxpatch
        !         do j = 1, params_glob%nypatch
        !             do iframe = 1, nframes
        !                 call starfile_table__addObject(mc_starfile)
        !                 call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_X, patched_centers(i,j,1)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_Y, patched_centers(i,j,2)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, patched_shifts(1,iframe,i,j)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, patched_shifts(2,iframe,i,j)/dpscale)
        !             enddo
        !         enddo
        !     enddo
        !     call starfile_table__write_ofile(mc_starfile)
        ! endif
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
    subroutine motion_correct_patched( bfac, rmsd_threshold, rmsd )
        real, intent(in)  :: bfac, rmsd_threshold
        real, intent(out) :: rmsd(2)     !< whether polynomial fitting was within threshold
        integer :: iframe
        if( l_BENCH ) t_patched = tic()
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe=1,nframes
            call movie_frames_scaled(iframe)%ifft
        end do
        !$omp end parallel do
        rmsd = huge(rmsd(1))
        call motion_patch%new(params_glob%scale*params_glob%trs)
        write(logfhandle,'(A,I2,A3,I2,A1)') '>>> PATCH-BASED REFINEMENT (',params_glob%nxpatch,' x ',params_glob%nypatch,')'
        call motion_patch%set_fitshifts(FITSHIFTS)
        call motion_patch%set_frameweights(frameweights)
        call motion_patch%set_fixed_frame(fixed_frame)
        call motion_patch%set_interp_fixed_frame(fixed_frame)
        call motion_patch%set_bfactor(bfac)
        if( trim(params_glob%algorithm).eq.'poly2' )then
            call motion_patch%correct_poly(hp, resstep, rmsd_threshold, movie_frames_scaled, patched_shift_fname, patched_polyn, shifts_toplot)
            rmsd = motion_patch%get_polyfit_rmsd()
        else
            call motion_patch%correct(hp, resstep, movie_frames_scaled, patched_shift_fname, shifts_toplot)
            rmsd = motion_patch%get_polyfit_rmsd()
            call motion_patch%get_poly4star(patched_polyn, patched_shifts, patched_centers)
        endif
        ! end if
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
            print *,'rt_dose          : ',rt_dose
        endif
    end subroutine motion_correct_patched_calc_sums

    subroutine motion_correct_patched_kill
        call motion_patch%kill
        call ftexp_transfmat_kill
    end subroutine motion_correct_patched_kill

    !!!!!!!!!!!!!!!!!!!!!!!!!
    !> motion_correction of DDD movie
    subroutine motion_correct_dev( movie_stack_fname, ctfvars, movie_sum, movie_sum_corrected,&
            &movie_sum_ctf, aniso_success, poly_rmsd, gainref, boxdata )
        use simple_motion_align_poly, only: motion_align_poly
        character(len=*),           intent(in)    :: movie_stack_fname  !< input filename of stack
        type(ctfparams),            intent(inout) :: ctfvars            !< CTF params
        type(image),                intent(inout) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        logical,                    intent(out)   :: aniso_success
        real,                       intent(out)   :: poly_rmsd
        character(len=*), optional, intent(in)    :: gainref            !< gain reference filename
        real,             optional, intent(in)    :: boxdata(:,:)
        type(motion_align_poly)    :: align_obj
        real, allocatable         :: iso_shifts(:,:)
        real(dp)                  :: star_polyn(36)
        real                      :: ave, sdev, var
        integer                   :: t
        logical                   :: err, err_stat
        ! initialise
        if( l_BENCH ) t_correct_iso_init = tic()
        call motion_correct_init(movie_stack_fname, ctfvars, err, movie_sum, gainref)
        ! deals with reference frame convention
        select case(trim(params_glob%mcconvention))
            case('relion','first')
                fixed_frame = 1
            case DEFAULT
                fixed_frame = nint(real(nframes)/2.)
        end select
        !$omp parallel do schedule(guided) default(shared) private(t) proc_bind(close)
        do t=1,nframes
            call movie_frames_scaled(t)%ifft
        enddo
        !$omp end parallel do
        if( err ) return
        call align_obj%new(movie_frames_scaled, fixed_frame)
        call align_obj%gen_patches_dimensions
        call align_obj%gen_tiles(boxdata)
        if( l_BENCH ) rt_correct_iso_init = toc(t_correct_iso_init)
        if( L_BENCH ) t_ = tic()
        call align_obj%align(params_glob%algorithm, poly_rmsd, aniso_success, patched_polyn, star_polyn)
        if( .not.aniso_success )then
            THROW_WARN('Polynomial fitting to patch-determined shifts was of insufficient quality')
            THROW_WARN('Only isotropic/stage-drift correction will be used')
        endif
        if( L_BENCH ) rt_aniso_align = toc(t_)
        ! solution
        call align_obj%plot_shifts(patched_shift_fname)
        call align_obj%get_iso_shifts(iso_shifts)
        call align_obj%get_weights(frameweights)
        shifts_toplot = iso_shifts
        call moment(frameweights, ave, sdev, var, err_stat)
        write(logfhandle,'(A,F7.4)')    '>>> AVERAGE PATCH & FRAMES CORRELATION: ', align_obj%get_corr()
        write(logfhandle,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT :', ave
        write(logfhandle,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS:', sdev
        write(logfhandle,'(a,7x,f7.4)') '>>> MIN WEIGHT     :', minval(frameweights)
        write(logfhandle,'(a,7x,f7.4)') '>>> MAX WEIGHT     :', maxval(frameweights)
        ! shift frames
        !$omp parallel do schedule(guided) default(shared) private(t) proc_bind(close)
        do t=1,nframes
            call movie_frames_scaled(t)%fft
            call movie_frames_scaled(t)%shift2Dserial(-iso_shifts(t,:))
            if( aniso_success ) call movie_frames_scaled(t)%ifft
        enddo
        !$omp end parallel do
        if( aniso_success )then
            ! dummy init to access interpolation routine
            call motion_patch%new(0.)
            call motion_patch%set_nframes(nframes)
            call motion_patch%set_fixed_frame(fixed_frame)
            call motion_patch%set_interp_fixed_frame(fixed_frame)
            call motion_patch%set_poly_coeffs(patched_polyn)
            ! generate micrograph & CTF micrograph
            call motion_correct_patched_calc_sums(movie_sum_corrected, movie_sum_ctf)
        else
            ! generate micrograph & CTF micrograph
            call motion_correct_iso_calc_sums(movie_sum_corrected, movie_sum_ctf)
        endif
        ! now patched_polyn is with reference to first frame for later star ouput
        patched_polyn = star_polyn
        ! report the sampling distance of the possibly scaled movies
        ctfvars%smpd = smpd_scaled
        call align_obj%kill
        if( L_BENCH )then
            print *,'rt_iso_init   : ',rt_correct_iso_init
            print *,'rt_aniso_align: ',rt_aniso_align
        endif
    end subroutine motion_correct_dev

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
        if( L_BENCH ) t_mic2spec = tic()
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
        call img_out%set_rmat(rmat_sum,.false.)
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
        if( L_BENCH ) print *,'rt_mic2spec:        ',toc(t_mic2spec)
    end subroutine motion_correct_mic2spec

    subroutine motion_correct_kill_common
        integer :: iframe
        if( allocated(shifts_toplot)      ) deallocate(shifts_toplot)
        if( allocated(frameweights)       ) deallocate(frameweights)
        if( allocated(cmat_sum)           ) deallocate(cmat_sum)
        if( allocated(pos_outliers)       ) deallocate(pos_outliers)
        call gain_img%kill
        if( allocated(movie_frames_scaled) )then
            do iframe=1,size(movie_frames_scaled)
                call movie_frames_scaled(iframe)%kill
            end do
            deallocate(movie_frames_scaled)
        endif
        call ftexp_transfmat_kill
        call eer%kill
        l_eer = .false.
    end subroutine motion_correct_kill_common

    ! COMMON PRIVATE UTILITY METHODS

    ! Following Grant & Grigorieff; eLife 2015;4:e06980
    ! Frames assumed in fourier space
    subroutine apply_dose_weighting( frames )
        class(image), intent(inout) :: frames(nframes)
        real, parameter   :: A=0.245, B=-1.665, C=2.81
        real, allocatable :: qs(:)
        real    :: acc_doses(nframes), spaFreqk, dose_per_frame
        real    :: twoNe, smpd, spafreq, limhsq,limksq
        integer :: nrflims(3,2), ldim(3), hphys,kphys, iframe, h,k
        if( .not.params_glob%l_dose_weight ) return
        if( .not.frames(1)%is_ft() ) THROW_HARD('Frames should be in in the Fourier domain')
        nrflims = frames(1)%loop_lims(2)
        smpd    = frames(1)%get_smpd()
        allocate(qs(nframes),source=0.)
        ! doses
        ldim           = frames(1)%get_ldim()
        limhsq         = (real(ldim(1))*smpd)**2.
        limksq         = (real(ldim(2))*smpd)**2.
        dose_per_frame = total_dose / real(nframes)
        do iframe=1,nframes
            acc_doses(iframe) = real(iframe) * dose_per_frame
        end do
        if( is_equal(kV,200.) )then
            acc_doses = acc_doses / 0.8
        else if( is_equal(kV,100.) )then
            acc_doses = acc_doses / 0.64
        endif
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
                qs = exp(-acc_doses/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                do iframe = 1,nframes
                    call frames(iframe)%mul_cmat_at([hphys,kphys,1], qs(iframe))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_dose_weighting

    ! gain correction, calculate image sum and identify outliers
    subroutine correct_gain( frames, gainref_fname )
        class(image),     intent(inout) :: frames(nframes)
        character(len=*), intent(in)    :: gainref_fname
        integer     :: iframe, ifoo
        write(logfhandle,'(a)') '>>> PERFORMING GAIN CORRECTION'
        if( l_eer )then
            call eer%prep_gainref(gainref_fname, gain_img)
        else
            if( fname2format(gainref_fname)=='L' )then
                THROW_HARD('''.gain'' files only for use with EER movies! correct_gain')
            endif
            call find_ldim_nptcls(gainref_fname,ldim_gain,ifoo)
            if( ldim_gain(1).ne.ldim(1) .or. ldim_gain(2).ne.ldim(2) )then
                THROW_HARD('Inconsistent dimensions between movie frames & gain reference! correct_gain')
            endif
            call gain_img%new(ldim_gain, smpd)
            call gain_img%read(gainref_fname)
        endif
        !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
        do iframe = 1,nframes
            call frames(iframe)%mul(gain_img)
        enddo
        !$omp end parallel do
    end subroutine correct_gain

    ! gain correction, calculate image sum and identify outliers
    subroutine cure_outliers( frames )
        use simple_image, only: image_ptr
        class(image), intent(inout) :: frames(nframes)
        integer, parameter   :: hwinsz = 5
        type(image_ptr)      :: prmats(nframes)
        real,        pointer :: prmat(:,:,:)
        real,    allocatable :: rsum(:,:), new_vals(:,:), vals(:)
        integer, allocatable :: pos_outliers_here(:,:)
        real    :: ave, sdev, var, lthresh,uthresh, l,u,localave
        integer :: iframe, noutliers, i,j,k,ii,jj, nvals, winsz, n
        logical :: outliers(ldim(1),ldim(2)), err
        allocate(rsum(ldim(1),ldim(2)),source=0.)
        if( l_BENCH ) t_cure = tic()
        write(logfhandle,'(a)') '>>> REMOVING DEAD/HOT PIXELS & FOURIER TRANSFORMING FRAMES'
        ! sum
        do iframe = 1,nframes
            call frames(iframe)%get_rmat_ptr(prmat)
            !$omp parallel workshare
            rsum(:,:) = rsum(:,:) + prmat(:ldim(1),:ldim(2),1)
            !$omp end parallel workshare
        enddo
        nullify(prmat)
        ! outliers detection
        call moment( rsum, ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh   = ave - NSIGMAS * sdev
        uthresh   = ave + NSIGMAS * sdev
        !$omp workshare
        where( rsum<lthresh .or. rsum>uthresh )
            outliers = .true.
        elsewhere
            outliers = .false.
        end where
        !$omp end workshare
        deallocate(rsum)
        ! cure
        noutliers = count(outliers)
        if( noutliers > 0 )then
            write(logfhandle,'(a,1x,i10)') '>>> # DEAD/HOT PIXELS:', noutliers
            write(logfhandle,'(a,1x,2f10.1)') '>>> AVERAGE (STDEV):  ', ave, sdev
            winsz = 2*HWINSZ+1
            nvals = winsz*winsz
            allocate(pos_outliers(2,noutliers))
            ! gather positions for star output only
            k = 0
            do j = 1,ldim(2)
                do i = 1,ldim(1)
                    if( outliers(i,j) )then
                        k = k + 1
                        pos_outliers(:,k) = [i,j]
                    endif
                enddo
            enddo
            ! add eer gain defects for curation but not star output
            if( l_eer .and. gain_img%exists() )then
                call gain_img%get_rmat_ptr(prmat)
                where( is_zero(prmat(:ldim(1),:ldim(2),1)) )
                    outliers = .true.
                end where
                nullify(prmat)
                noutliers = count(outliers)
                if( noutliers > 0 )then
                    write(logfhandle,'(a,1x,i8)') '>>> # DEAD/HOT PIXELS + EER GAIN DEFFECTS:', noutliers
                    ! gather defect positions again for curation
                    allocate(pos_outliers_here(2,noutliers),source=-1)
                    k = 0
                    do j = 1,ldim(2)
                        do i = 1,ldim(1)
                            if( outliers(i,j) )then
                                k = k + 1
                                pos_outliers_here(:,k) = [i,j]
                            endif
                        enddo
                    enddo
                else
                    ! nothing to do
                    return
                endif
            else
                pos_outliers_here = pos_outliers
            endif
            allocate(new_vals(noutliers,nframes),vals(nvals))
            ave  = ave / real(nframes)
            sdev = sdev / real(nframes)
            uthresh = uthresh / real(nframes)
            lthresh = lthresh / real(nframes)
            !$omp parallel do default(shared) private(iframe,k,i,j,n,ii,jj,vals,l,u,localave)&
            !$omp proc_bind(close) schedule(static)
            do iframe=1,nframes
                call frames(iframe)%get_rmat_ptr(prmats(iframe)%rmat)
                ! calulate new values
                do k = 1,noutliers
                    i = pos_outliers_here(1,k)
                    j = pos_outliers_here(2,k)
                    n = 0
                    do jj = j-HWINSZ,j+HWINSZ
                        if( jj < 1 .or. jj > ldim(2) ) cycle
                        do ii = i-HWINSZ,i+HWINSZ
                            if( ii < 1 .or. ii > ldim(1) ) cycle
                            if( outliers(ii,jj) ) cycle
                            n = n + 1
                            vals(n) = prmats(iframe)%rmat(ii,jj,1)
                        enddo
                    enddo
                    if( n > 1 )then
                        if( real(n)/real(nvals) < 0.85 )then
                            ! high defect area
                            l = minval(vals(:n))
                            u = maxval(vals(:n))
                            if( abs(u-l) < sdev/1000.0 ) u = uthresh
                            localave = sum(vals(:n)) / real(n)
                            new_vals(k,iframe) = gasdev(localave, sdev, [l,u])
                        else
                            new_vals(k,iframe) = median_nocopy(vals(:n))
                        endif
                    else
                        new_vals(k,iframe) = gasdev(ave, sdev, [lthresh,uthresh])
                    endif
                enddo
                ! substitute
                do k = 1,noutliers
                    i = pos_outliers_here(1,k)
                    j = pos_outliers_here(2,k)
                    prmats(iframe)%rmat(i,j,1) = new_vals(k,iframe)
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine cure_outliers

        !>  Utility to calculate the number fractions, # eer frames per fraction and adjusted total_dose
    !   while minimizing the number of leftover frames given total dose, total # of eer frames & a dose target
    subroutine calc_eer_fraction(n_eer_frames, fraction_dose_target, tot_dose, nfractions, eerfraction)
        integer, intent(in)    :: n_eer_frames
        real,    intent(in)    :: fraction_dose_target
        real,    intent(inout) :: tot_dose
        integer, intent(out)   :: nfractions, eerfraction
        real    :: dose_per_eer_frame
        integer :: ffloor, fceil, nffloor, nfceil
        ffloor  = floor(real(n_eer_frames) / (tot_dose / fraction_dose_target))
        fceil   = ffloor+1
        nffloor = floor(real(n_eer_frames) / real(ffloor))
        nfceil  = floor(real(n_eer_frames) / real(fceil))
        if( nffloor*ffloor > nfceil*fceil )then
            nfractions   = nffloor
            eer_fraction = ffloor
        else
            nfractions   = nfceil
            eer_fraction = fceil
        endif
        ! adjusting total dose
        dose_per_eer_frame = tot_dose / real(n_eer_frames)
        tot_dose           = dose_per_eer_frame * real(eerfraction) * real(nfractions)
    end subroutine calc_eer_fraction

end module simple_motion_correct
