!@descr: backend-neutral half-map FSC, cFAR, and resolution diagnostics shared by the gridding and PCG reconstruction paths
module simple_halfmap_diagnostics
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_image_msk,  only: image_msk
use simple_fsc,        only: phase_rand_fsc, fsc_area_score_result
implicit none

public :: halfmap_diagnostics_result, evaluate_halfmap_pair, write_halfmap_diagnostics
public :: support_provenance_fname, write_support_provenance, read_support_provenance
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
    !! is backend-neutral and applies NO mask of its own (2026-09-09): both
    !! backends ship halves that already carry the soft spherical support at
    !! msk_crop (the PCG solve support; the gridding restoration applies the
    !! identical mask3D_soft after deapodization), and masking them again
    !! here would square the soft edge. Callers own the input representation,
    !! refine3D artifact filenames, and every workflow-level write. The only
    !! files produced here are the fscu/fsct/fscn state arrays phase_rand_fsc
    !! persists internally on the envfsc path, identically for both backends.
    !! The caller-owned inputs are never modified; the envelope preprocessing
    !! and the in-place Fourier transforms operate on copies. With envfsc
    !! enabled the density-envelope mask is returned through the optional
    !! envmask so the caller can write the automask artifact; the optional
    !! cones argument returns the conical FSC result needed for directional
    !! regularization.
    subroutine evaluate_halfmap_pair( params, state, even, odd, average, diagnostics, envmask, cones, &
        &l_pair_support_constrained )
        class(parameters),                      intent(in)    :: params
        integer,                                intent(in)    :: state
        class(image),                           intent(in)    :: even, odd, average
        type(halfmap_diagnostics_result),       intent(out)   :: diagnostics
        class(image),                 optional, intent(inout) :: envmask
        class(fsc_area_score_result), optional, intent(inout) :: cones
        logical,                      optional, intent(in)    :: l_pair_support_constrained
        type(image)                 :: work_even, work_odd
        type(image_msk)             :: envmask_work
        type(fsc_area_score_result) :: cones_local
        real, allocatable :: fsc_t(:), fsc_n(:), res(:)
        integer :: nyq
        logical :: l_envfsc_preproc
        nyq = even%get_filtsz()
        ! The envelope-masking + phase-randomization preprocessing
        ! approximates, after the fact, an estimate the solver could not
        ! constrain. When the caller states that the pair was ALREADY
        ! density-envelope-constrained in the estimator (the PCG solve
        ! support), masking it again would double-count and the preproc is
        ! skipped. The decision follows the support actually installed for
        ! this pair, never the backend name (code review 2026-09-02 P1):
        ! a PCG run whose base solves fell back to the spherical support
        ! (missing/startvol reference, bootstrap) passes .false. and gets
        ! the ordinary envelope-corrected FSC. The envelope itself is still
        ! derived and returned, because the automask artifact has other
        ! consumers (postprocess envfsc, final rec).
        l_envfsc_preproc = params%l_envfsc
        if( present(l_pair_support_constrained) )then
            if( l_pair_support_constrained ) l_envfsc_preproc = .false.
        endif
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
            allocate(diagnostics%fsc(nyq), source=0.)
        endif
        ! calc_fsc_area_score converts the work maps to Fourier space in place,
        ! so the radial FSC below reads the same representation
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

    !> Support-provenance sidecar of a shipped state volume
    !! (<vol>_pcg_support.txt, the historical name kept for compatibility).
    !! Records whether the shipped half pair was estimated inside the
    !! conservative density envelope (PCG P H P with the density support) or
    !! carries the plain soft spherical support at msk_crop, and what kind of
    !! estimate the primary pair is: a PCG base, regularized or bootstrap-
    !! mixed solve, or a gridding restoration (2026-09-09; the gridding
    !! products carry the same spherical support after deapodization).
    !! Consumers: the PCG trailing bootstrap reads the support field for its
    !! lag-one FSC pair; the PCG base warm-start selector reads the solve
    !! kind (gridding products never seed it); postprocess skips its post-hoc
    !! mask for any volume carrying the sidecar.
    function support_provenance_fname( volname ) result( fname )
        type(string), intent(in) :: volname
        type(string) :: fname
        fname = swap_suffix(add2fbody(volname, MRC_EXT, '_pcg_support'), TXT_EXT, MRC_EXT)
    end function support_provenance_fname

    subroutine write_support_provenance( volname, l_constrained, solve_kind )
        type(string),     intent(in) :: volname
        logical,          intent(in) :: l_constrained
        character(len=*), intent(in) :: solve_kind
        type(string) :: fname
        integer :: funit
        select case( trim(solve_kind) )
            case( 'base', 'regularized', 'mixed', 'gridding' )
            case default
                THROW_HARD('invalid solve kind for the support provenance sidecar')
        end select
        fname = support_provenance_fname(volname)
        call fopen(funit, file=fname, status='replace', action='write')
        write(funit,'(A)') 'solve_support='//merge('density', 'sphere ', l_constrained)
        write(funit,'(A)') 'solve_kind='//trim(solve_kind)
        call fclose(funit)
        call fname%kill
    end subroutine write_support_provenance

    subroutine read_support_provenance( volname, l_constrained, l_found, solve_kind, l_kind_found )
        type(string),               intent(in)  :: volname
        logical,                    intent(out) :: l_constrained, l_found
        character(len=*), optional, intent(out) :: solve_kind
        logical,          optional, intent(out) :: l_kind_found
        type(string) :: fname
        character(len=64) :: line
        integer :: funit, io_stat
        l_constrained = .false.
        l_found       = .false.
        if( present(solve_kind) )   solve_kind   = ''
        if( present(l_kind_found) ) l_kind_found = .false.
        fname = support_provenance_fname(volname)
        if( .not. file_exists(fname) )then
            call fname%kill
            return
        endif
        call fopen(funit, file=fname, status='old', action='read')
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( index(line, 'solve_support=') == 1 )then
                l_found       = .true.
                l_constrained = index(line, 'density') > 0
            else if( index(line, 'solve_kind=') == 1 )then
                if( present(solve_kind) ) solve_kind = adjustl(line(len('solve_kind=')+1:))
                if( present(l_kind_found) ) l_kind_found = .true.
            endif
        enddo
        call fclose(funit)
        call fname%kill
    end subroutine read_support_provenance

    subroutine kill_halfmap_diagnostics_result( self )
        class(halfmap_diagnostics_result), intent(inout) :: self
        if( allocated(self%fsc) ) deallocate(self%fsc)
        self%res_fsc05   = 0.
        self%res_fsc0143 = 0.
        self%cfar        = 0.
    end subroutine kill_halfmap_diagnostics_result

end module simple_halfmap_diagnostics
