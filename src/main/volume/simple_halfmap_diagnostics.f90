!@descr: backend-neutral half-map FSC, cFAR, and resolution diagnostics shared by the gridding and PCG reconstruction paths
module simple_halfmap_diagnostics
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_image_msk,  only: image_msk
use simple_fsc,        only: phase_rand_fsc, fsc_area_score_result
implicit none

public :: halfmap_diagnostics_result, evaluate_halfmap_pair, write_halfmap_diagnostics
private
#include "simple_local_flags.inc"

! cFAR cone construction shared by every half-map diagnostic site
integer, parameter :: CFAR_NDIRS               = 256
real,    parameter :: CFAR_CONE_HALF_ANGLE_DEG = 20.
real,    parameter :: CFAR_FSC_THRESHOLD       = 0.143
integer, parameter :: CFAR_MIN_COUNT           = 1

type :: halfmap_diagnostics_result
    real, allocatable :: fsc(:)
    real              :: res_fsc05   = 0.
    real              :: res_fsc0143 = 0.
    real              :: cfar        = 0.
  contains
    procedure :: kill => kill_halfmap_diagnostics_result
end type halfmap_diagnostics_result

contains

    !> Evaluate the half-map FSC, cFAR, and Nyquist-clamped FSC=0.5/0.143
    !! resolutions for one explicitly prepared real-space pair. The evaluator
    !! is backend-neutral: callers own the input representation (gridding
    !! passes its legacy undeapodized base pair through its adapter, PCG its
    !! restored base solve), the spherical mask radius, refine3D artifact
    !! filenames, and every workflow-level write. The only files produced
    !! here are the fscu/fsct/fscn state arrays phase_rand_fsc persists
    !! internally on the envfsc path, identically for both backends. The
    !! caller-owned inputs are never modified; masking operates on copies.
    !! With envfsc enabled the density-envelope mask is returned through the
    !! optional envmask so the caller can write the automask artifact; the
    !! optional cones argument returns the conical FSC result needed for
    !! directional regularization.
    subroutine evaluate_halfmap_pair( params, state, even, odd, average, spherical_mask_radius, &
        &diagnostics, envmask, cones )
        class(parameters),                      intent(in)    :: params
        integer,                                intent(in)    :: state
        class(image),                           intent(in)    :: even, odd, average
        real,                                   intent(in)    :: spherical_mask_radius
        type(halfmap_diagnostics_result),       intent(out)   :: diagnostics
        class(image),                 optional, intent(inout) :: envmask
        class(fsc_area_score_result), optional, intent(inout) :: cones
        type(image)                 :: work_even, work_odd
        type(image_msk)             :: envmask_work
        type(fsc_area_score_result) :: cones_local
        real, allocatable :: fsc_t(:), fsc_n(:), res(:)
        integer :: nyq
        logical :: l_envfsc_preproc
        nyq = even%get_filtsz()
        ! The envelope-masking + phase-randomization preprocessing is a
        ! GRIDDING-path construction: it approximates, after the fact, an
        ! estimate the solver could not constrain. On the PCG backend the
        ! same envelope is installed as the SOLVE SUPPORT (pcg_priors.md dev
        ! item 5), so the pair handed here is already envelope-constrained
        ! and masking it again would double-count. The envelope itself is
        ! still derived and returned, because the automask artifact has other
        ! consumers (postprocess envfsc, matcher fallback, final rec).
        l_envfsc_preproc = params%l_envfsc .and. trim(params%rec_backend) /= 'pcg'
        if( params%l_envfsc .and. present(envmask) )then
            call envmask_work%automask3D(params, average, .false., lp_override=params%envmsklp)
            call envmask%copy(envmask_work)
        endif
        if( l_envfsc_preproc )then
            ! density-envelope masking with phase-randomized FSC correction
            if( .not. present(envmask) ) &
                &call envmask_work%automask3D(params, average, .false., lp_override=params%envmsklp)
            call phase_rand_fsc(even, odd, envmask_work, state, nyq, diagnostics%fsc, fsc_t, fsc_n)
            call work_even%copy(even)
            call work_odd%copy(odd)
            call work_even%zero_env_background(envmask_work)
            call work_odd%zero_env_background(envmask_work)
            call work_even%mul(envmask_work)
            call work_odd%mul(envmask_work)
            deallocate(fsc_t, fsc_n)
            call envmask_work%kill_bimg
        else
            call envmask_work%kill_bimg
            call work_even%copy(even)
            call work_odd%copy(odd)
            call work_even%mask3D_soft(spherical_mask_radius, backgr=0.)
            call work_odd%mask3D_soft(spherical_mask_radius, backgr=0.)
            allocate(diagnostics%fsc(nyq), source=0.)
        endif
        ! calc_fsc_area_score converts the work maps to Fourier space in place,
        ! so the spherical-mask FSC below reads the same masked representation
        if( present(cones) )then
            call cones%new(work_even, CFAR_NDIRS, CFAR_CONE_HALF_ANGLE_DEG, CFAR_FSC_THRESHOLD, &
                &CFAR_MIN_COUNT)
            call cones%calc_fsc_area_score(work_even, work_odd, state=state)
            diagnostics%cfar = cones%cfar
        else
            call cones_local%new(work_even, CFAR_NDIRS, CFAR_CONE_HALF_ANGLE_DEG, CFAR_FSC_THRESHOLD, &
                &CFAR_MIN_COUNT)
            call cones_local%calc_fsc_area_score(work_even, work_odd, state=state)
            diagnostics%cfar = cones_local%cfar
            call cones_local%kill
        endif
        if( .not. l_envfsc_preproc ) call work_even%fsc(work_odd, diagnostics%fsc)
        res = get_resarr(params%box_crop, params%smpd_crop)
        call get_resolution(diagnostics%fsc, res, diagnostics%res_fsc05, diagnostics%res_fsc0143)
        diagnostics%res_fsc05   = max(diagnostics%res_fsc05,   2. * params%smpd_crop)
        diagnostics%res_fsc0143 = max(diagnostics%res_fsc0143, 2. * params%smpd_crop)
        call work_even%kill
        call work_odd%kill
        deallocate(res)
    end subroutine evaluate_halfmap_pair

    !> Write the half-map resolution text report to an explicit filename. The
    !! writer contains no backend-specific policy; filename selection stays
    !! with the workflow caller.
    subroutine write_halfmap_diagnostics( diagnostics, box, smpd, fname )
        class(halfmap_diagnostics_result), intent(in) :: diagnostics
        integer,                           intent(in) :: box
        real,                              intent(in) :: smpd
        class(string),                     intent(in) :: fname
        real, allocatable :: res(:)
        integer :: k, fnr
        if( .not. allocated(diagnostics%fsc) ) THROW_HARD('No half-map FSC available to write')
        res = get_resarr(box, smpd)
        call fopen(fnr, FILE=fname, STATUS='REPLACE', action='WRITE')
        do k = 1, min(size(res), size(diagnostics%fsc))
            write(fnr,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', diagnostics%fsc(k)
        end do
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', diagnostics%res_fsc05
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', diagnostics%res_fsc0143
        write(fnr,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', diagnostics%cfar
        call fclose(fnr)
        deallocate(res)
    end subroutine write_halfmap_diagnostics

    ! DIAGNOSTIC LIFECYCLE

    subroutine kill_halfmap_diagnostics_result( self )
        class(halfmap_diagnostics_result), intent(inout) :: self
        if( allocated(self%fsc) ) deallocate(self%fsc)
        self%res_fsc05   = 0.
        self%res_fsc0143 = 0.
        self%cfar        = 0.
    end subroutine kill_halfmap_diagnostics_result

end module simple_halfmap_diagnostics
