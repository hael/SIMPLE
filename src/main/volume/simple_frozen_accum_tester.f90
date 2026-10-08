!@descr: unit tests of the solve3D_addon frozen accumulator sets (simple_frozen_accum)
! Frozen set F + cohort equals the direct union (gridding sums/rho, PCG B/D, restored and solved
! maps); a zero cohort gives F exactly. Plus context round trip and refusals. Box 16, seeded noise.
module simple_frozen_accum_tester
use simple_defs,              only: OSMPL_PAD_FAC, STDLEN
use simple_string,            only: string
use simple_fileio,            only: del_file, fopen, fclose
use simple_image,             only: image
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_parameters,        only: parameters
use simple_sp_project,        only: sp_project
use simple_sym,               only: sym
use simple_memoize_ft_maps,   only: forget_ft_maps, memoize_ft_maps
use simple_type_defs,         only: CTFFLAG_NO, OBJFUN_EUCLID, ctfparams, fplane_type
use simple_reconstructor,     only: reconstructor
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_OP_KERNEL
use simple_refine3D_fnames,   only: refine3D_frozen_rec_fbody, refine3D_frozen_manifest_fname, &
    &refine3D_frozen_pcg_fname
use simple_frozen_accum,      only: frozen_accum
use simple_test_utils
implicit none
private
public :: run_all_frozen_accum_tests

integer, parameter :: BOX      = 16
integer, parameter :: NFROZEN  = 14   ! frozen planes (even/odd alternate)
integer, parameter :: NCOHORT  = 10   ! cohort planes
integer, parameter :: PCG_ITS  = 6
integer, parameter :: SEED     = 20260926
real,    parameter :: SMPD     = 2.0
real,    parameter :: LAMBDA   = 1.e-2
! Single-precision sums of O(10) terms reassociated between F u C and F + C:
! relative differences of a few ulps (~1e-7) of the largest magnitude; the
! restored/solved maps propagate them through linear operators of O(1) gain.
real,    parameter :: RAW_RELTOL = 1.e-5
real,    parameter :: MAP_RELTOL = 1.e-4
character(len=*), parameter :: CTX_FNAME  = 'tmp_frozen_accum_tester_context.txt'
character(len=*), parameter :: CTX_FNAME2 = 'tmp_frozen_accum_tester_context_2.txt'
character(len=*), parameter :: RUN_ID    = 'tester_run_1'

contains

    subroutine run_all_frozen_accum_tests()
        write(*,'(A)') '**** running all frozen accumulator tests ****'
        call test_context_round_trip()
        call test_context_refusals()
        call test_gridding_union_equals_sum()
        call test_pcg_union_equals_sum()
    end subroutine run_all_frozen_accum_tests

    ! ---- fixtures -------------------------------------------------------------

    !> two inherited states of NFROZEN and 1 frozen particles; the working
    !! project has NFROZEN+NCOHORT+3 rows, the frozen project its first
    !! NFROZEN+NCOHORT (three appended rows)
    subroutine make_context( ctx, rid, backend )
        type(frozen_accum), intent(inout) :: ctx
        character(len=*),   intent(in)    :: rid, backend
        call ctx%new(rid, backend, NFROZEN + NCOHORT + 3, NFROZEN + NCOHORT, [NFROZEN, 1])
    end subroutine make_context

    !> NFROZEN + NCOHORT well-spread orientations; the first NFROZEN are frozen
    subroutine make_orientations( all_oris )
        type(oris), intent(inout) :: all_oris
        call all_oris%new(NFROZEN + NCOHORT, .false.)
        call all_oris%spiral()
    end subroutine make_orientations

    !> seeded white-noise particle images, one per orientation
    subroutine make_images( imgs )
        real, allocatable, intent(out) :: imgs(:,:,:)
        integer :: i, j, k
        real    :: r
        allocate(imgs(BOX,BOX,NFROZEN+NCOHORT))
        call set_fixed_seed(SEED)
        do k = 1, NFROZEN + NCOHORT
            do j = 1, BOX
                do i = 1, BOX
                    call random_number(r)
                    imgs(i,j,k) = r - 0.5
                enddo
            enddo
        enddo
    end subroutine make_images

    !> the two text files hold the same lines
    logical function same_text( fname1, fname2 ) result( l_same )
        character(len=*), intent(in) :: fname1, fname2
        character(len=256) :: lines1(64), lines2(64)
        integer :: n1, n2
        call read_text(fname1, lines1, n1)
        call read_text(fname2, lines2, n2)
        l_same = n1 == n2 .and. n1 > 0
        if( l_same ) l_same = all(lines1(1:n1) == lines2(1:n2))
    end function same_text

    subroutine read_text( fname, lines, n )
        character(len=*), intent(in)  :: fname
        character(len=*), intent(out) :: lines(:)
        integer,          intent(out) :: n
        integer :: funit, io_stat
        n = 0
        call fopen(funit, file=string(fname), status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 ) return
        do while( n < size(lines) )
            read(funit,'(A)',iostat=io_stat) lines(n+1)
            if( io_stat /= 0 ) exit
            n = n + 1
        enddo
        call fclose(funit)
    end subroutine read_text

    subroutine write_text( fname, lines )
        character(len=*), intent(in) :: fname
        character(len=*), intent(in) :: lines(:)
        integer :: funit, io_stat, i
        call fopen(funit, file=string(fname), status='REPLACE', action='WRITE', iostat=io_stat)
        do i = 1, size(lines)
            write(funit,'(A)') trim(lines(i))
        enddo
        call fclose(funit)
    end subroutine write_text

    real function max_rel_cdiff( a, b ) result( d )
        complex, intent(in) :: a(:,:,:), b(:,:,:)
        d = maxval(abs(a - b)) / max(maxval(abs(b)), tiny(1.))
    end function max_rel_cdiff

    real function max_rel_rdiff( a, b ) result( d )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        d = maxval(abs(a - b)) / max(maxval(abs(b)), tiny(1.))
    end function max_rel_rdiff

    ! ---- run context ----------------------------------------------------------

    subroutine test_context_round_trip()
        type(frozen_accum)    :: ctx, back
        character(len=STDLEN) :: msg
        integer :: status
        write(*,'(A)') 'test_context_round_trip'
        call make_context(ctx, RUN_ID, 'gridding')
        call ctx%write(string(CTX_FNAME))
        call back%read(string(CTX_FNAME), status, msg)
        call assert_int(0, status, 'a written context reads back: '//trim(msg))
        ! run identifier, backend, layout and counts: the context read back is
        ! written again record for record
        call back%write(string(CTX_FNAME2))
        call assert_true(same_text(CTX_FNAME, CTX_FNAME2), 'every context record round trips')
        call assert_int(NFROZEN, back%get_nfrozen_state(1), 'context state 1 count round trip')
        call assert_int(1,       back%get_nfrozen_state(2), 'context state 2 count round trip')
        call back%validate('gridding', 2, NFROZEN+NCOHORT+3, .false., status, msg)
        call assert_int(0, status, 'the context accepts a consumer of its working project')
        call back%validate('gridding', 2, NFROZEN+NCOHORT, .true., status, msg)
        call assert_int(0, status, 'the context accepts a producer of its frozen project')
        call back%validate('pcg', 2, NFROZEN+NCOHORT+3, .false., status, msg)
        call assert_true(status /= 0, 'the context refuses another backend')
        call back%validate('gridding', 1, NFROZEN+NCOHORT+3, .false., status, msg)
        call assert_true(status /= 0, 'the context refuses another state layout')
        call back%validate('gridding', 2, NFROZEN+NCOHORT+4, .false., status, msg)
        call assert_true(status /= 0, 'the context refuses another particle index space')
        call del_file(CTX_FNAME)
        call del_file(CTX_FNAME2)
        call ctx%kill
        call back%kill
    end subroutine test_context_round_trip

    subroutine test_context_refusals()
        type(frozen_accum)    :: ctx
        character(len=STDLEN) :: msg
        character(len=64)     :: good(10)
        integer :: status
        write(*,'(A)') 'test_context_refusals'
        good = [character(len=64) :: 'solve3D_addon_frozen_context 3', 'run_id r1', 'backend pcg', &
            &'nstates 2', 'nrows 30', 'nrows_frozen 25', 'nfrozen 7', 'nfrozen_state 1 3', 'nfrozen_state 2 4', 'end']
        call write_text(CTX_FNAME, good)
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_int(0, status, 'a hand-written valid context is accepted: '//trim(msg))
        call del_file(CTX_FNAME)
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'a missing context is refused')
        call write_text(CTX_FNAME, [character(len=64) :: 'solve3D_addon_frozen_context 2', good(2:)])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'an unsupported schema version is refused')
        call write_text(CTX_FNAME, [character(len=64) :: 'solve3D_addon_frozen_set 1', good(2:)])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'another schema is refused')
        call write_text(CTX_FNAME, good(1:9))
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'a truncated context (no end marker) is refused')
        call write_text(CTX_FNAME, [character(len=64) :: good(1:7), 'nfrozen_state 1 3', 'end'])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'a context missing a state count is refused')
        call write_text(CTX_FNAME, [character(len=64) :: good(1:6), 'nfrozen 8', good(8:)])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'state counts that do not sum to nfrozen are refused')
        call write_text(CTX_FNAME, [character(len=64) :: good(1:5), good(7:)])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'a context without the frozen project''s row count is refused')
        call write_text(CTX_FNAME, [character(len=64) :: good(1:9), 'frozen_rows 1 2 3', 'end'])
        call ctx%read(string(CTX_FNAME), status, msg)
        call assert_true(status /= 0, 'an unknown context field is refused')
        call del_file(CTX_FNAME)
        call ctx%kill
    end subroutine test_context_refusals

    ! ---- gridding -------------------------------------------------------------

    subroutine test_gridding_union_equals_sum()
        type(parameters), target :: params
        type(sp_project)         :: project
        type(reconstructor)      :: f_even, f_odd, c_even, c_odd, u_even, u_odd, scratch_e, scratch_o
        type(image)              :: map_c, map_u
        type(frozen_accum)       :: ctx, other, weighted, recount
        type(oris)               :: all_oris
        type(ori)                :: o
        type(sym)                :: c1sym
        type(ctfparams)          :: ctfparms
        type(fplane_type)        :: fplane
        type(image)              :: obs, obs_pad
        real,    allocatable     :: imgs(:,:,:), rho_c(:,:,:), rho_u(:,:,:), rho_f(:,:,:), rmat_c(:,:,:), rmat_u(:,:,:)
        complex, allocatable     :: cmat_c(:,:,:), cmat_u(:,:,:), cmat_f(:,:,:)
        character(len=STDLEN)    :: msg
        type(string)             :: fname
        integer :: i, status
        write(*,'(A)') 'test_gridding_union_equals_sum'
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        call f_even%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call f_odd%new_accumulator( params, project, expand=.true., wthreads=.false.)
        call c_even%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call c_odd%new_accumulator( params, project, expand=.true., wthreads=.false.)
        call u_even%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call u_odd%new_accumulator( params, project, expand=.true., wthreads=.false.)
        call scratch_e%new_accumulator(params, project, expand=.false., wthreads=.false.)
        call scratch_o%new_accumulator(params, project, expand=.false., wthreads=.false.)
        call make_orientations(all_oris)
        call make_images(imgs)
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call o%new(.false.)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        ! F planes then C planes into U, the same planes into F or C
        do i = 1, NFROZEN + NCOHORT
            call all_oris%get_ori(i, o)
            call obs%set_rmat(imgs(:,:,i:i), .false.)
            call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
            call obs_pad%fft()
            call obs_pad%gen_fplane4rec([0,BOX/2], SMPD, ctfparms, [0.,0.], fplane)
            if( mod(i,2) == 0 )then
                call u_even%insert_plane_oversamp(c1sym, o, fplane)
                if( i <= NFROZEN )then
                    call f_even%insert_plane_oversamp(c1sym, o, fplane)
                else
                    call c_even%insert_plane_oversamp(c1sym, o, fplane)
                endif
            else
                call u_odd%insert_plane_oversamp(c1sym, o, fplane)
                if( i <= NFROZEN )then
                    call f_odd%insert_plane_oversamp(c1sym, o, fplane)
                else
                    call c_odd%insert_plane_oversamp(c1sym, o, fplane)
                endif
            endif
        enddo
        call f_even%compress_exp()
        call f_odd%compress_exp()
        call c_even%compress_exp()
        call c_odd%compress_exp()
        call u_even%compress_exp()
        call u_odd%compress_exp()
        ! Production reads every partial through the raw MRC accumulator format,
        ! which stores box/2 complex values per row (the h = box/2 column is
        ! not persisted), so the fair comparison passes C and F u C through
        ! the same writer and validated reader (as state 2) and adds F (state 1)
        call make_context(ctx, RUN_ID, 'gridding')
        call ctx%write_gridding_set(1, f_even, f_odd)
        call ctx%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_int(0, status, 'a freshly written gridding set validates: '//trim(msg))
        call round_trip(c_even, c_odd)
        call round_trip(u_even, u_odd)
        cmat_c = c_even%get_cmat()
        cmat_u = u_even%get_cmat()
        call assert_true(max_rel_cdiff(cmat_c, cmat_u) > 100.*RAW_RELTOL, &
            &'negative control: the cohort sums alone differ from F u C')
        call ctx%add_gridding_set(1, c_even, c_odd, scratch_e, scratch_o)
        cmat_c = c_even%get_cmat()
        cmat_u = u_even%get_cmat()
        call c_even%get_rho_copy(rho_c)
        call u_even%get_rho_copy(rho_u)
        call assert_true(max_rel_cdiff(cmat_c, cmat_u) < RAW_RELTOL, 'gridding even sums: C + F equals F u C')
        call assert_true(max_rel_rdiff(rho_c, rho_u)   < RAW_RELTOL, 'gridding even densities: C + F equals F u C')
        cmat_c = c_odd%get_cmat()
        cmat_u = u_odd%get_cmat()
        call c_odd%get_rho_copy(rho_c)
        call u_odd%get_rho_copy(rho_u)
        call assert_true(max_rel_cdiff(cmat_c, cmat_u) < RAW_RELTOL, 'gridding odd sums: C + F equals F u C')
        call assert_true(max_rel_rdiff(rho_c, rho_u)   < RAW_RELTOL, 'gridding odd densities: C + F equals F u C')
        ! restored halves
        call map_c%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call map_u%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call c_even%restore_final(map_c, preserve_numerator=.true.)
        call u_even%restore_final(map_u, preserve_numerator=.true.)
        rmat_c = map_c%get_rmat()
        rmat_u = map_u%get_rmat()
        call assert_true(max_rel_rdiff(rmat_c, rmat_u) < MAP_RELTOL, 'restored gridding half: C + F equals F u C')
        ! a cohort with no contribution is the frozen term alone, exactly, on
        ! every persisted Fourier component and on the density
        call c_even%reset
        call c_odd%reset
        call ctx%add_gridding_set(1, c_even, c_odd, scratch_e, scratch_o)
        cmat_c = c_even%get_cmat()
        cmat_f = f_even%get_cmat()
        call c_even%get_rho_copy(rho_c)
        call f_even%get_rho_copy(rho_f)
        call assert_true(all(cmat_c(:BOX/2,:,:) == cmat_f(:BOX/2,:,:)) .and. all(rho_c == rho_f), &
            &'zero cohort + F is F exactly (coefficient one)')
        call delete_gridding_set(2)
        ! refusals: another run, another grid, a size-mismatched or missing component
        call make_context(other, 'another_run', 'gridding')
        call other%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a set of another add-on run is refused')
        call ctx%gridding_set_status(1, BOX, 1.1*SMPD, status, msg)
        call assert_true(status /= 0, 'a set at another sampling is refused')
        call ctx%gridding_set_status(1, BOX+2, SMPD, status, msg)
        call assert_true(status /= 0, 'no set is found for another box')
        call ctx%gridding_set_status(2, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'no set is found for a state that was not written')
        ! a consumer of another reconstruction weighting, another frozen count
        call ctx%write(string(CTX_FNAME))
        call weighted%load(string(CTX_FNAME), 'gridding', 2, NFROZEN+NCOHORT+3, OBJFUN_EUCLID, producer=.false., wset_id=[0_8, 0_8])
        call weighted%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_true(status /= 0 .and. index(msg, 'weighting') > 0, 'a set of another weighting is refused')
        call recount%new(RUN_ID, 'gridding', NFROZEN+NCOHORT+3, NFROZEN+NCOHORT, [NFROZEN-1, 2])
        call recount%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a set whose frozen count differs from the context is refused')
        call del_file(CTX_FNAME)
        fname = string('rho_')//refine3D_frozen_rec_fbody(1, BOX)//'_odd.mrc'
        call scratch_o%write_rho(fname) ! same size: rewrite is accepted
        call ctx%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_int(0, status, 'a same-size component rewrite keeps the recorded sizes')
        call truncate_file(fname)
        call ctx%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a component of the wrong size is refused')
        call del_file(fname)
        call ctx%gridding_set_status(1, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a set with a missing component is refused')
        ! cleanup
        call delete_gridding_set(1)
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call obs%kill()
        call obs_pad%kill()
        call map_c%kill()
        call map_u%kill()
        call f_even%kill(); call f_odd%kill(); call c_even%kill(); call c_odd%kill()
        call u_even%kill(); call u_odd%kill(); call scratch_e%kill(); call scratch_o%kill()
        call all_oris%kill()
        call project%kill()
        call ctx%kill()
        call other%kill()
        call weighted%kill()
        call recount%kill()
        if( allocated(fplane%cmplx_plane)    ) deallocate(fplane%cmplx_plane)
        if( allocated(fplane%ctfsq_plane)    ) deallocate(fplane%ctfsq_plane)
        if( allocated(fplane%transfer_plane) ) deallocate(fplane%transfer_plane)

    contains

        subroutine round_trip( rec_e, rec_o )
            type(reconstructor), intent(inout) :: rec_e, rec_o
            call ctx%write_gridding_set(2, rec_e, rec_o)
            call rec_e%reset
            call rec_o%reset
            call ctx%add_gridding_set(2, rec_e, rec_o, scratch_e, scratch_o)
        end subroutine round_trip

    end subroutine test_gridding_union_equals_sum

    subroutine truncate_file( fname )
        class(string), intent(in) :: fname
        integer :: funit, io_stat
        call del_file(fname)
        call fopen(funit, file=fname, status='REPLACE', action='WRITE', access='STREAM', iostat=io_stat)
        write(funit) 1.0
        call fclose(funit)
    end subroutine truncate_file

    subroutine delete_gridding_set( state )
        integer, intent(in) :: state
        type(string) :: fbody
        fbody = refine3D_frozen_rec_fbody(state, BOX)
        call del_file(fbody//'_even.mrc')
        call del_file(fbody//'_odd.mrc')
        call del_file(string('rho_')//fbody//'_even.mrc')
        call del_file(string('rho_')//fbody//'_odd.mrc')
        call del_file(refine3D_frozen_manifest_fname(state, BOX))
        call fbody%kill
    end subroutine delete_gridding_set

    ! ---- PCG -------------------------------------------------------------------

    subroutine test_pcg_union_equals_sum()
        type(reconstructor_pcg) :: sampler, pcg_f, pcg_c, pcg_u, pcg_z
        type(frozen_accum)      :: ctx, other, weighted
        type(oris)              :: all_oris, f_oris, c_oris
        type(ori)               :: o
        type(image)             :: obs
        real,    allocatable    :: imgs(:,:,:), x_c(:,:,:), x_u(:,:,:), d_c(:,:,:), d_f(:,:,:)
        complex, allocatable    :: planes(:,:,:), b_c(:,:,:), b_f(:,:,:)
        character(len=STDLEN)   :: msg
        real    :: b_relerr, d_relerr
        integer :: i, lims2(2,2), status, nadded
        write(*,'(A)') 'test_pcg_union_equals_sum'
        call make_orientations(all_oris)
        call make_images(imgs)
        call f_oris%new(NFROZEN, .false.)
        call c_oris%new(NCOHORT, .false.)
        call o%new(.false.)
        do i = 1, NFROZEN
            call all_oris%get_ori(i, o)
            call f_oris%set_ori(i, o)
        enddo
        do i = 1, NCOHORT
            call all_oris%get_ori(NFROZEN+i, o)
            call c_oris%set_ori(i, o)
        enddo
        call sampler%new(BOX, SMPD, LAMBDA)
        lims2 = sampler%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NFROZEN+NCOHORT))
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        do i = 1, NFROZEN + NCOHORT
            call obs%set_rmat(imgs(:,:,i:i), .false.)
            call obs%fft()
            planes(:,:,i) = sampler%extract_native_plane(obs)
            call obs%ifft()
        enddo
        ! F alone, C alone, F u C in one pass (F first, then C)
        call accumulate(pcg_f, f_oris, planes(:,:,1:NFROZEN))
        call accumulate(pcg_c, c_oris, planes(:,:,NFROZEN+1:))
        call accumulate(pcg_u, all_oris, planes)
        call make_context(ctx, RUN_ID, 'pcg')
        call ctx%write_pcg_half(1, 0, BOX, pcg_f, NFROZEN)
        call ctx%pcg_half_status(1, 0, BOX, SMPD, status, msg)
        call assert_int(0, status, 'a freshly written PCG half validates: '//trim(msg))
        call pcg_c%compare_raw_accum(pcg_u, b_relerr, d_relerr)
        call assert_true(b_relerr > 100.*RAW_RELTOL .and. d_relerr > 100.*RAW_RELTOL, &
            &'negative control: the cohort B and D alone differ from F u C')
        call ctx%add_pcg_half(1, 0, BOX, SMPD, pcg_c, nadded)
        call assert_int(NFROZEN, nadded, 'the frozen PCG half reports its frozen particle count')
        call pcg_c%compare_raw_accum(pcg_u, b_relerr, d_relerr)
        call assert_true(b_relerr < RAW_RELTOL, 'PCG raw B: C + F equals F u C')
        call assert_true(d_relerr < RAW_RELTOL, 'PCG raw D: C + F equals F u C')
        ! a cohort with no contribution: an open empty reduction plus F is F exactly
        call pcg_z%new(BOX, SMPD, LAMBDA)
        call pcg_z%begin_reduction
        call ctx%add_pcg_half(1, 0, BOX, SMPD, pcg_z, nadded)
        call pcg_z%get_raw_accum(b_c, d_c)
        call pcg_f%get_raw_accum(b_f, d_f)
        call assert_true(all(b_c == b_f) .and. all(d_c == d_f), 'zero cohort + F is F exactly (coefficient one)')
        ! solved maps at a fixed iteration count from zero
        call pcg_c%end_accum(.true.)
        call pcg_c%set_op_mode(PCG_OP_KERNEL)
        call pcg_u%end_accum(.true.)
        call pcg_u%set_op_mode(PCG_OP_KERNEL)
        allocate(x_c(BOX,BOX,BOX), x_u(BOX,BOX,BOX), source=0.)
        call pcg_c%solve_accum(x_c, maxits=PCG_ITS, rtol=0.)
        call pcg_u%solve_accum(x_u, maxits=PCG_ITS, rtol=0.)
        call assert_true(max_rel_rdiff(x_c, x_u) < MAP_RELTOL, 'solved PCG half: C + F equals F u C')
        ! refusals
        call make_context(other, 'another_run', 'pcg')
        call other%pcg_half_status(1, 0, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a PCG half of another add-on run is refused')
        call ctx%pcg_half_status(1, 1, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a half that was not written is not found')
        call ctx%pcg_half_status(1, 0, BOX+2, SMPD, status, msg)
        call assert_true(status /= 0, 'no PCG half is found for another box')
        call ctx%pcg_half_status(1, 0, BOX, 1.1*SMPD, status, msg)
        call assert_true(status /= 0, 'a PCG half at another sampling is refused')
        call ctx%write(string(CTX_FNAME))
        call weighted%load(string(CTX_FNAME), 'pcg', 2, NFROZEN+NCOHORT+3, OBJFUN_EUCLID, producer=.false., wset_id=[0_8, 0_8])
        call weighted%pcg_half_status(1, 0, BOX, SMPD, status, msg)
        call assert_true(status /= 0, 'a PCG half of another weighting is refused')
        call del_file(CTX_FNAME)
        call del_file(refine3D_frozen_pcg_fname(1, BOX, 'even'))
        ! cleanup
        call sampler%kill(); call pcg_f%kill(); call pcg_c%kill(); call pcg_u%kill(); call pcg_z%kill()
        call obs%kill()
        call o%kill()
        call all_oris%kill(); call f_oris%kill(); call c_oris%kill()
        call ctx%kill(); call other%kill(); call weighted%kill()

    contains

        subroutine accumulate( op, os, pl )
            type(reconstructor_pcg), intent(inout) :: op
            type(oris),              intent(inout) :: os
            complex,                 intent(in)    :: pl(lims2(1,1):,lims2(2,1):,:)
            call op%new(BOX, SMPD, LAMBDA)
            call op%prep_particles(os, use_ctf=.false.)
            call op%begin_accum
            call op%accumulate_batch(pl, size(pl,3), 1)
        end subroutine accumulate

    end subroutine test_pcg_union_equals_sum

end module simple_frozen_accum_tester
