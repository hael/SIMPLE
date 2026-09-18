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
public :: copy_support_provenance, rename_support_provenance, remove_support_provenance
private
#include "simple_local_flags.inc"

! cFAR cone construction shared by every half-map diagnostic site
integer, parameter :: CFAR_NDIRS               = 256
real,    parameter :: CFAR_CONE_HALF_ANGLE_DEG = 20.
real,    parameter :: CFAR_FSC_THRESHOLD       = 0.143
integer, parameter :: CFAR_MIN_COUNT           = 1
character(len=*), parameter :: FSC_MODE_UNKNOWN = &
    &'backend=unknown support=unknown posthoc_mask=none phase_randomization=no'

type :: halfmap_diagnostics_result
    real, allocatable :: fsc(:)
    real              :: res_fsc05   = 0.
    real              :: res_fsc0143 = 0.
    real              :: cfar        = 0.
    character(len=200):: fsc_mode    = FSC_MODE_UNKNOWN
  contains
    procedure :: kill => kill_halfmap_diagnostics_result
end type halfmap_diagnostics_result

contains

    !> Evaluate the half-map FSC, cFAR, and Nyquist-clamped FSC=0.5/0.143
    !! resolutions for one explicitly prepared real-space pair, with explicit
    !! backend/support/mask provenance. Gridding may apply the selected density
    !! or NU envelope post hoc and phase-randomize; PCG never does either
    !! because its envelope belongs inside the estimator. Caller-owned inputs
    !! are never modified. Density mode returns the generated mask through the
    !! optional envmask; NU mode requires the caller to supply envmask. The
    !! optional cones result supports directional regularization.
    subroutine evaluate_halfmap_pair( params, state, even, odd, average, diagnostics, backend, envmask, cones, &
        &l_pair_support_constrained, support_kind, mask_kind )
        class(parameters),                      intent(in)    :: params
        integer,                                intent(in)    :: state
        class(image),                           intent(in)    :: even, odd, average
        type(halfmap_diagnostics_result),       intent(out)   :: diagnostics
        character(len=*),                       intent(in)    :: backend
        class(image),                 optional, intent(inout) :: envmask
        class(fsc_area_score_result), optional, intent(inout) :: cones
        logical,                      optional, intent(in)    :: l_pair_support_constrained
        character(len=*),             optional, intent(in)    :: support_kind, mask_kind
        type(image)                 :: work_even, work_odd
        type(image_msk)             :: envmask_work
        type(fsc_area_score_result) :: cones_local
        real, allocatable :: fsc_t(:), fsc_n(:), res(:)
        integer :: nyq
        character(len=16) :: support_kind_here, mask_kind_here, posthoc_kind
        logical :: l_phase_randomization
        nyq = even%get_filtsz()
        support_kind_here = 'sphere'
        mask_kind_here    = 'none'
        if( present(support_kind) ) support_kind_here = trim(support_kind)
        if( present(mask_kind) )    mask_kind_here    = trim(mask_kind)
        if( present(l_pair_support_constrained) )then
            if( l_pair_support_constrained .and. .not. present(support_kind) ) support_kind_here = 'envelope'
        endif
        select case(trim(backend))
            case('gridding','pcg')
            case default
                THROW_HARD('half-map diagnostics backend must be gridding or pcg')
        end select
        select case(trim(mask_kind_here))
            case('none')
            case('density')
                call envmask_work%automask3D(params, average, .false., lp_override=params%envmsklp)
                if( present(envmask) ) call envmask%copy(envmask_work)
            case('nu')
                if( .not. present(envmask) ) THROW_HARD('NU FSC mask was selected but not supplied')
                call envmask_work%copy(envmask)
            case default
                THROW_HARD('unsupported half-map FSC mask kind')
        end select
        ! Phase-randomized solvent correction is deliberately gridding-only.
        ! PCG reports the FSC of the estimate formed on its installed support.
        l_phase_randomization = trim(backend) == 'gridding' .and. trim(mask_kind_here) /= 'none'
        posthoc_kind = 'none'
        if( l_phase_randomization ) posthoc_kind = mask_kind_here
        diagnostics%fsc_mode = 'backend='//trim(backend)//' support='//trim(support_kind_here)//&
            &' posthoc_mask='//trim(posthoc_kind)//&
            &' phase_randomization='//merge('yes', 'no ', l_phase_randomization)
        write(logfhandle,'(A,I0,A)') '>>> FSC MODE: STATE ', state, ', '//trim(diagnostics%fsc_mode)
        if( l_phase_randomization )then
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
        if( .not. l_phase_randomization ) call work_even%fsc(work_odd, diagnostics%fsc)
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
        write(fnr,'(A,1X,A)')    '>>> FSC MODE                             :', trim(diagnostics%fsc_mode)
        call fclose(fnr)
        deallocate(res)
    end subroutine write_halfmap_diagnostics

    ! DIAGNOSTIC LIFECYCLE

    !> Support-provenance sidecar of a shipped state volume
    !! (<vol>_pcg_support.txt, the historical name kept for compatibility).
    !! Records the shipped half pair's support kind (sphere, density, NU,
    !! explicit, or mixed) and what kind of estimate the primary pair is: a
    !! PCG base, regularized or bootstrap-mixed solve, or gridding restoration.
    !! Consumers: the PCG trailing bootstrap reads the support field for its
    !! lag-one FSC pair; postprocess skips its post-hoc mask for any volume
    !! carrying the sidecar. The solve kind is recorded for provenance (the
    !! former PCG base warm-start selector that read it went with the
    !! cross-iteration warm starts, 2026-09-10).
    function support_provenance_fname( volname ) result( fname )
        type(string), intent(in) :: volname
        type(string) :: fname
        fname = swap_suffix(add2fbody(volname, MRC_EXT, '_pcg_support'), TXT_EXT, MRC_EXT)
    end function support_provenance_fname

    subroutine write_support_provenance( volname, l_constrained, solve_kind, support_kind, solvent_prior )
        type(string),     intent(in) :: volname
        logical,          intent(in) :: l_constrained
        character(len=*), intent(in) :: solve_kind
        character(len=*), optional, intent(in) :: support_kind
        character(len=*), optional, intent(in) :: solvent_prior !< e.g. 'soft lambda_rel=1.00' (pcg_solvent=yes)
        type(string) :: fname
        character(len=16) :: support_kind_here
        integer :: funit
        select case( trim(solve_kind) )
            case( 'base', 'regularized', 'mixed', 'gridding' )
            case default
                THROW_HARD('invalid solve kind for the support provenance sidecar')
        end select
        fname = support_provenance_fname(volname)
        support_kind_here = merge('density', 'sphere ', l_constrained)
        if( present(support_kind) ) support_kind_here = trim(support_kind)
        select case(trim(support_kind_here))
            case('sphere','density','nu','explicit','mixed')
            case default
                THROW_HARD('invalid support kind for the support provenance sidecar')
        end select
        call fopen(funit, file=fname, status='replace', action='write')
        write(funit,'(A)') 'solve_support='//trim(support_kind_here)
        write(funit,'(A)') 'solve_kind='//trim(solve_kind)
        if( present(solvent_prior) ) write(funit,'(A)') 'solvent_prior='//trim(solvent_prior)
        call fclose(funit)
        call fname%kill
    end subroutine write_support_provenance

    subroutine read_support_provenance( volname, l_constrained, l_found, solve_kind, l_kind_found, support_kind )
        type(string),               intent(in)  :: volname
        logical,                    intent(out) :: l_constrained, l_found
        character(len=*), optional, intent(out) :: solve_kind
        logical,          optional, intent(out) :: l_kind_found
        character(len=*), optional, intent(out) :: support_kind
        type(string) :: fname
        character(len=64) :: line
        integer :: funit, io_stat
        l_constrained = .false.
        l_found       = .false.
        if( present(solve_kind) )   solve_kind   = ''
        if( present(l_kind_found) ) l_kind_found = .false.
        if( present(support_kind) ) support_kind = ''
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
                if( present(support_kind) ) support_kind = adjustl(line(len('solve_support=')+1:))
                l_constrained = index(line, 'sphere') == 0
            else if( index(line, 'solve_kind=') == 1 )then
                if( present(solve_kind) ) solve_kind = adjustl(line(len('solve_kind=')+1:))
                if( present(l_kind_found) ) l_kind_found = .true.
            endif
        enddo
        call fclose(funit)
        call fname%kill
    end subroutine read_support_provenance

    !> The volume and its sidecar are one artifact: every copy, rename or
    !! fresh write of a state volume goes through these so no valid map loses
    !! its provenance and no stale sidecar survives beside a map that has none
    !! (stage snapshots, final copies, start volumes; review 2026-09-09 P1)
    subroutine copy_support_provenance( src_vol, dest_vol )
        type(string), intent(in) :: src_vol, dest_vol
        type(string) :: src, dest
        src  = support_provenance_fname(src_vol)
        dest = support_provenance_fname(dest_vol)
        if( file_exists(src) )then
            call simple_copy_file(src, dest)
        else if( file_exists(dest) )then
            call del_file(dest)
        endif
        call src%kill
        call dest%kill
    end subroutine copy_support_provenance

    subroutine rename_support_provenance( src_vol, dest_vol )
        type(string), intent(in) :: src_vol, dest_vol
        type(string) :: src, dest
        src  = support_provenance_fname(src_vol)
        dest = support_provenance_fname(dest_vol)
        if( file_exists(src) )then
            call simple_rename(src, dest, overwrite=.true.)
        else if( file_exists(dest) )then
            call del_file(dest)
        endif
        call src%kill
        call dest%kill
    end subroutine rename_support_provenance

    subroutine remove_support_provenance( vol )
        type(string), intent(in) :: vol
        type(string) :: fname
        fname = support_provenance_fname(vol)
        if( file_exists(fname) ) call del_file(fname)
        call fname%kill
    end subroutine remove_support_provenance

    subroutine kill_halfmap_diagnostics_result( self )
        class(halfmap_diagnostics_result), intent(inout) :: self
        if( allocated(self%fsc) ) deallocate(self%fsc)
        self%res_fsc05   = 0.
        self%res_fsc0143 = 0.
        self%cfar        = 0.
        self%fsc_mode    = FSC_MODE_UNKNOWN
    end subroutine kill_halfmap_diagnostics_result

end module simple_halfmap_diagnostics
