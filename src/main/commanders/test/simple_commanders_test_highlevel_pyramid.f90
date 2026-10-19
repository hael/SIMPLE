!@descr: solve3D_pyramid workflow gate: the three-state phantom of 3 600 particles in 60 exact views, 2D-classified, run end to end and checked against the planted states and maps
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_pyramid
use simple_commanders_mkcavgs,          only: commander_make_cavgs
use simple_commanders_solve3D_pyramid,  only: commander_solve3D_pyramid
use simple_flex_pca_application_tester, only: FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
use simple_test_state_fixture,          only: state_fixture, FIX_VIEWS
use simple_test_truth_metrics,          only: match_states, mapped_accuracy, adjusted_rand_index, state_recalls, &
    &match_maps, dock_both_hands, discrimination, finite_min, offdiag
use simple_state_assignment,            only: REFS_DIR, docked_pair_cc, read_refs_record, fit_mskdiam
use simple_image,                       only: image
use simple_imghead,                     only: find_ldim_nptcls, find_img_smpd
use simple_sp_project,                  only: sp_project
use simple_test_gate,                   only: test_gate
use simple_ui,                          only: make_ui
implicit none
#include "simple_local_flags.inc"

!> the three blob states of the multi-state test fixture in exact views: 60 directions, each an exact
!! 2D class of 60 particles over the three states (equal populations; the state of every particle in
!! the fixture's truth document); one stack, so the first set at whole stacks is the phantom itself
integer, parameter :: PYR_NPTCLS  = 3600
integer, parameter :: PYR_NSTATES = 3
integer, parameter :: PYR_NVIEWS  = 60
integer, parameter :: PYR_SEED    = 202610151
real,    parameter :: PYR_SNR     = 1.0
real,    parameter :: DIFF_FRAC   = 0.1 !< the difference mask of the discrimination, as in the multi-state tests
integer, parameter :: DIFF_DILATE = 1
integer, allocatable :: truth_states(:)
!> provisional detection floors, until a profile gives them (the floor rule): the assignment's
!! accuracy against the planted states under the one-to-one mapping (random labels give about 1/3),
!! and the docked correlation of each reference and each refined map with the truth it maps onto
!! (the hand free)
real,    parameter :: PYR_PURITY_MIN = 0.45
real,    parameter :: PYR_REF_CC_MIN = 0.70
real,    parameter :: PYR_MAP_CC_MIN = 0.70

contains

module subroutine exec_test_solve3D_pyramid_blobs( self, cline )
    class(commander_test_solve3D_pyramid_blobs), intent(inout) :: self
    class(cmdline),                              intent(inout) :: cline
    class(parameters), allocatable :: params
    type(string) :: cwd_saved, fixture_root
    integer :: nthr, status
    logical :: all_ok
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_solve3D_pyramid_blobs_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: could not enter fixture directory')
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    all_ok = .true.
    call run_pyramid_gate(nthr, filepath(cwd_saved, 'metrics.tsv'), all_ok)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: could not restore original directory')
    CWD_GLOB = cwd_saved%to_char()
    if( all_ok )then
        call simple_rmdir(fixture_root)
        write(logfhandle,'(A)') 'PASS: solve3D_pyramid assigns the planted states of the three-state phantom'
        call simple_end('**** SIMPLE_TEST_SOLVE3D_PYRAMID_BLOBS NORMAL STOP ****')
    else
        THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED')
    endif
end subroutine exec_test_solve3D_pyramid_blobs

!> the phantom, its class averages (make_cavgs of its exact classes), solve3D_pyramid on it, and its
!! products against the truth
subroutine run_pyramid_gate( nthr, metrics_file, all_ok )
    integer,       intent(in)    :: nthr
    class(string), intent(in)    :: metrics_file
    logical,       intent(inout) :: all_ok
    type(test_gate) :: gate
    class(sp_project), allocatable :: proj
    type(state_fixture) :: fix
    type(string)     :: root, run_dir, projfile, exec_dir, result_proj, stk
    type(string)     :: truths(PYR_NSTATES)
    real    :: smpd
    integer :: ncls, status, k
    call make_ui
    call simple_getcwd(root)
    CWD_GLOB = root%to_char()
    call gate%new(metrics_file)
    fix%family  = 'blob'
    fix%mode    = FIX_VIEWS
    fix%nstates = PYR_NSTATES
    fix%pops(1:PYR_NSTATES) = 1. / real(PYR_NSTATES)
    fix%nptcls  = PYR_NPTCLS
    fix%nviews  = PYR_NVIEWS
    fix%snr     = PYR_SNR
    fix%nthr    = nthr
    call fix%build
    do k = 1, PYR_NSTATES
        truths(k) = fix%truth_fname(k)
    enddo
    truth_states = read_states(fix%truth_doc())
    run_dir  = root//'/run'
    projfile = run_dir//'/phantom.simple'
    call simple_mkdir(run_dir)
    call simple_copy_file(fix%project_fname(), projfile)
    ! the input of solve3D_pyramid is 2D-classified: the class averages of the exact classes
    call make_phantom_cavgs(run_dir, projfile, nthr)
    allocate(proj)
    call proj%read(projfile)
    ncls = 0
    call proj%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
    call proj%kill
    deallocate(proj)
    call gate%check('pyramid: the phantom has its 60 class averages', ncls == PYR_NVIEWS)
    if( ncls == PYR_NVIEWS )then
        call run_pyramid(run_dir, projfile, nthr, exec_dir)
        result_proj = exec_dir//'/phantom.simple'
        call gate%check('pyramid: the report is written', file_exists(exec_dir//'/solve3D_pyramid_report.txt'))
        call gate%check('pyramid: the references are recorded', file_exists(exec_dir//'/'//REFS_DIR//'/state_refs.txt'))
        call gate%check('pyramid: the result project is written', file_exists(result_proj))
        if( file_exists(result_proj) ) call inspect_pyramid(root, exec_dir, result_proj, truths, gate)
    endif
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()
    all_ok = all_ok .and. gate%passed()
    call gate%kill
end subroutine run_pyramid_gate

!> make_cavgs of the phantom's ptcl2D classes in run_dir: the class averages, their halves and the
!! FRCs registered in projfile
subroutine make_phantom_cavgs( run_dir, projfile, nthr )
    class(string), intent(in) :: run_dir, projfile
    integer,       intent(in) :: nthr
    class(commander_make_cavgs), allocatable :: xmake
    class(cmdline),              allocatable :: cline
    integer :: status
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: could not enter the run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xmake, cline)
    call cline%set('prg',      'make_cavgs')
    call cline%set('projfile', projfile)
    call cline%set('refs',     'cavgs_phantom.mrc')
    call cline%set('ncls',     PYR_NVIEWS)
    call cline%set('mkdir',    'no')
    call cline%set('mskdiam',  FLEX_PHANTOM_MSKDIAM)
    call cline%set('ml_reg',   'no')
    call cline%set('nthr',     nthr)
    call xmake%execute(cline)
    call cline%kill
    deallocate(xmake, cline)
end subroutine make_phantom_cavgs

!> solve3D_pyramid on projfile, in a new execution directory of run_dir (exec_dir), where it stays
subroutine run_pyramid( run_dir, projfile, nthr, exec_dir )
    class(string), intent(in)  :: run_dir, projfile
    integer,       intent(in)  :: nthr
    type(string),  intent(out) :: exec_dir
    class(commander_solve3D_pyramid), allocatable :: xpyramid
    class(cmdline),                   allocatable :: cline
    integer :: status
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: could not enter the run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xpyramid, cline)
    call cline%set('prg',      'solve3D_pyramid')
    call cline%set('projfile', projfile)
    call cline%set('mkdir',    'yes')
    call cline%set('nstates',  PYR_NSTATES)
    call cline%set('mskdiam',  FLEX_PHANTOM_MSKDIAM)
    call cline%set('nthr',     nthr)
    call set_fixed_seed(PYR_SEED, propagate=.true.)
    call xpyramid%execute(cline)
    call simple_getcwd(exec_dir)
    call cline%kill
    deallocate(xpyramid, cline)
end subroutine run_pyramid

!> solve3D_pyramid's products against the truth, under one mapping each: the delivered states onto
!! the planted ones (match_states: the accuracy, the adjusted Rand index, the smallest recall), the
!! references and the refined maps onto the truths by their docked correlations (match_maps), and the
!! discrimination of the refined maps docked onto their mapped truths
subroutine inspect_pyramid( root, exec_dir, result_proj, truths, gate )
    class(string),   intent(in)    :: root, exec_dir, result_proj, truths(:)
    type(test_gate), intent(inout) :: gate
    class(sp_project), allocatable :: proj
    type(string), allocatable :: ref_vols(:)
    type(string)              :: inspect_dir, vol, map_copy
    type(image)               :: docked(PYR_NSTATES), truth_imgs(PYR_NSTATES)
    integer, allocatable :: states(:), map(:), map_refs(:), vol_states(:), shuffled(:), map_shuf(:)
    real,    allocatable :: cc_refs(:,:)
    real    :: smpd, ref_smpd, acc, ccd, ccm, margins(PYR_NSTATES), diffcorr(PYR_NSTATES,PYR_NSTATES)
    integer :: s, t, nassigned, nvols, box, ref_box, nest, nmask, k
    logical :: l_docked
    inspect_dir = root//'/inspect'
    call simple_mkdir(inspect_dir)
    allocate(proj)
    call proj%read(result_proj)
    ! the assignment against the planted states, one-to-one
    states    = proj%os_ptcl3D%get_all_asint('state')
    if( size(states) /= PYR_NPTCLS ) THROW_HARD('SOLVE3D_PYRAMID_BLOBS FAILED: the result project lost particles')
    nassigned = count(states >= 1)
    nest      = max(1, maxval(states))
    map       = match_states(truth_states, states, PYR_NSTATES, nest)
    acc       = mapped_accuracy(truth_states, states, map)
    call gate%report('pyramid: fraction of the particles assigned', real(nassigned) / real(PYR_NPTCLS))
    call gate%metric('pyramid: assignment accuracy under the one-to-one mapping', acc, PYR_PURITY_MIN, &
        &acc >= PYR_PURITY_MIN)
    call gate%report('pyramid: adjusted Rand index', adjusted_rand_index(truth_states, states))
    ! the negative control: the truth labels shuffled across the particles must fail the floor
    shuffled = truth_states
    call set_fixed_seed(PYR_SEED + 1)
    call shuffle(shuffled)
    map_shuf = match_states(shuffled, states, PYR_NSTATES, nest)
    call gate%check('pyramid: shuffled truth labels fail the accuracy floor', &
        &.not. (mapped_accuracy(shuffled, states, map_shuf) >= PYR_PURITY_MIN))
    call gate%report('pyramid: smallest per-state recall', minval(state_recalls(truth_states, states, map, PYR_NSTATES)))
    ! each reference against the truth it maps onto, docked in both hands
    if( read_refs_record(exec_dir//'/'//REFS_DIR, PYR_NSTATES, ref_vols, ref_box, ref_smpd) )then
        call gate%report('pyramid: reference volumes', real(size(ref_vols)))
        allocate(cc_refs(size(ref_vols), size(truths)))
        do s = 1,size(ref_vols)
            map_copy = inspect_dir//'/ref_state'//int2str(s)//'.mrc'
            call simple_copy_file(ref_vols(s), map_copy)
            do t = 1,size(truths)
                cc_refs(s,t) = truth_docked_cc(truths(t), map_copy, inspect_dir, 'ref'//int2str(s)//'_truth'//int2str(t))
            enddo
        enddo
        map_refs = match_maps(cc_refs)
        do s = 1,size(ref_vols)
            if( map_refs(s) < 1 ) cycle
            call gate%metric('pyramid: reference '//int2str(s)//' docked onto its mapped truth', cc_refs(s,map_refs(s)), &
                &PYR_REF_CC_MIN, cc_refs(s,map_refs(s)) >= PYR_REF_CC_MIN)
        enddo
    else
        call gate%check('pyramid: the references can be read', .false.)
    endif
    ! each refined state's map against the truth it maps onto, by the label mapping
    nvols = 0
    allocate(vol_states(0))
    do s = 1,nest
        if( .not. proj%isthere_in_osout('vol', s) ) cycle
        call proj%get_vol('vol', s, vol, smpd, box)
        if( .not. file_exists(vol) ) cycle
        nvols      = nvols + 1
        vol_states = [vol_states, s]
        map_copy   = inspect_dir//'/refined_state'//int2str(s)//'.mrc'
        call simple_copy_file(vol, map_copy)
        if( map(s) < 1 ) cycle
        call gate%metric('pyramid: refined state '//int2str(s)//' docked onto its mapped truth', &
            &truth_docked_cc(truths(map(s)), map_copy, inspect_dir, 'map'//int2str(s)), PYR_MAP_CC_MIN, &
            &truth_docked_cc(truths(map(s)), map_copy, inspect_dir, 'map'//int2str(s)) >= PYR_MAP_CC_MIN)
    enddo
    call gate%check('pyramid: at least two refined states are registered', nvols >= 2)
    ! the discrimination of the refined maps at the phantom's box, docked onto their mapped truths
    l_docked = .false.
    do k = 1, PYR_NSTATES
        call truth_imgs(k)%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD)
        call truth_imgs(k)%read(truths(k))
        call docked(k)%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD)
    enddo
    nmask = 0
    do k = 1, PYR_NSTATES
        s = findloc(map, k, dim=1)
        if( s < 1 ) cycle
        map_copy = inspect_dir//'/refined_state'//int2str(s)//'.mrc'
        if( .not. file_exists(map_copy) ) cycle
        call dock_both_hands(truths(k), map_copy, FLEX_PHANTOM_MSKDIAM, 40., 8., 'disc'//int2str(k), &
            &inspect_dir//'/docked_state'//int2str(k)//'.mrc', ccd, ccm)
        call docked(k)%read(inspect_dir//'/docked_state'//int2str(k)//'.mrc')
        nmask = nmask + 1
    enddo
    if( nmask == PYR_NSTATES )then
        call discrimination(docked, truth_imgs, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
        call gate%report('pyramid: discrimination smallest margin', finite_min(margins))
        call gate%report('pyramid: discrimination smallest difference correlation', finite_min(offdiag(diffcorr)))
    else
        ! a truth without a map discriminates nothing: the lowest values
        call gate%report('pyramid: discrimination smallest margin', -2.)
        call gate%report('pyramid: discrimination smallest difference correlation', -1.)
    endif
    do k = 1, PYR_NSTATES
        call truth_imgs(k)%kill
        call docked(k)%kill
    enddo
    call proj%kill
    deallocate(proj)
end subroutine inspect_pyramid

!> the true state of each particle of the fixture's truth document
function read_states( doc ) result( st )
    use simple_oris, only: oris
    class(string), intent(in) :: doc
    integer, allocatable :: st(:)
    type(oris) :: truth
    call truth%new(PYR_NPTCLS, is_ptcl=.true.)
    call truth%read(doc, [1, PYR_NPTCLS])
    st = nint(truth%get_all('state'))
    call truth%kill
end function read_states

!> The docked correlation of the map map_fname with the truth truth_fname brought to its grid (Fourier
!! cropped from the phantom's box when the map's is smaller), the hand free (docked_pair_cc); -1 when
!! the map's box is larger than the phantom's
real function truth_docked_cc( truth_fname, map_fname, dir, tag ) result( cc )
    class(string),    intent(in) :: truth_fname, map_fname, dir
    character(len=*), intent(in) :: tag
    type(image)  :: truth
    type(string) :: truth_grid
    real    :: smpd, cc_direct, cc_mirror
    integer :: ldim(3), nimgs
    cc = -1.
    call find_ldim_nptcls(map_fname, ldim, nimgs)
    smpd = find_img_smpd(map_fname)
    if( ldim(1) > FLEX_PHANTOM_BOX ) return
    truth_grid = dir//'/'//tag//'_truth_grid.mrc'
    call truth%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
    call truth%read(truth_fname)
    if( ldim(1) < FLEX_PHANTOM_BOX )then
        call truth%fft()
        call truth%clip_inplace([ldim(1),ldim(1),ldim(1)])
        call truth%ifft()
    endif
    call truth%set_smpd(smpd)
    call truth%write(truth_grid)
    call truth%kill
    cc = docked_pair_cc(truth_grid, map_fname, smpd, fit_mskdiam(FLEX_PHANTOM_MSKDIAM, ldim(1), smpd), cc_direct, cc_mirror)
    write(logfhandle,'(A,F7.4,A,F7.4)') '>>> '//tag//' docked correlation: direct ', cc_direct, ', mirrored ', cc_mirror
end function truth_docked_cc

end submodule simple_commanders_test_highlevel_pyramid
