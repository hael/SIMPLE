!@descr: solve3D_pyramid workflow gate: the three-state phantom of 3 600 particles in 60 exact views, 2D-classified, run end to end and checked against the planted states and maps
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_pyramid
use simple_commanders_mkcavgs,          only: commander_make_cavgs
use simple_commanders_solve3D_pyramid,  only: commander_solve3D_pyramid
use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project, &
    &FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
use simple_state_assignment,            only: REFS_DIR, docked_pair_cc, read_refs_record, fit_mskdiam
use simple_image,                       only: image
use simple_imghead,                     only: find_ldim_nptcls, find_img_smpd
use simple_sp_project,                  only: sp_project
use simple_test_gate,                   only: test_gate
use simple_ui,                          only: make_ui
implicit none
#include "simple_local_flags.inc"

!> the phantom of the class-average state refinement's gate (simple_commanders_test_highlevel_flex):
!! particle i in state mod(i-1,3)+1, its ptcl2D class its direction, one of 60, an exact 2D class of
!! 60 particles; one stack, so the first set at whole stacks is the phantom itself (a first set of a
!! third would hold the first 20 directions only, a cap of the sphere)
integer, parameter :: PYR_NPTCLS  = 3600
integer, parameter :: PYR_NSTATES = 3
integer, parameter :: PYR_NVIEWS  = 60
integer, parameter :: PYR_SEED    = 202610151
!> provisional detection floors, to be set from the first run (solve3D_pyramid plan, section 6):
!! the assignment's purity against the planted states (each delivered state recoded to the majority
!! truth of its particles; random labels give about 1/3), and the docked correlation of each
!! reference and each refined map with a planted truth (the hand free)
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
    call run_pyramid_gate(nthr, all_ok)
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
subroutine run_pyramid_gate( nthr, all_ok )
    integer, intent(in)    :: nthr
    logical, intent(inout) :: all_ok
    type(test_gate) :: gate
    class(sp_project), allocatable :: proj
    type(string)     :: root, run_dir, projfile, exec_dir, result_proj, truth_mean, truth_diff, oritab, stack, stk
    type(string)     :: truths(PYR_NSTATES)
    real    :: smpd
    integer :: ncls, status
    call make_ui
    call simple_getcwd(root)
    CWD_GLOB = root%to_char()
    call gate%new(string('metrics.tsv'))
    call create_flex_pca_phantom_fixture(root, nthr, truths(1), truths(2), truth_mean, truth_diff, oritab, stack, &
        &nptcls=PYR_NPTCLS, nstates=PYR_NSTATES, nviews=PYR_NVIEWS, truth_c=truths(3))
    run_dir  = root//'/run'
    projfile = run_dir//'/phantom.simple'
    call simple_mkdir(run_dir)
    call build_flex_pca_phantom_project(projfile, stack, oritab, nptcls=PYR_NPTCLS, nstates=PYR_NSTATES, &
        &nviews=PYR_NVIEWS)
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

!> solve3D_pyramid's products against the truth: the assignment's purity against the planted states
!! (particle i in state mod(i-1,3)+1), each reference docked onto its closest truth, and each
!! refined state's map docked onto the majority truth of its particles
subroutine inspect_pyramid( root, exec_dir, result_proj, truths, gate )
    class(string),   intent(in)    :: root, exec_dir, result_proj, truths(:)
    type(test_gate), intent(inout) :: gate
    class(sp_project), allocatable :: proj
    type(string), allocatable :: ref_vols(:)
    type(string)              :: inspect_dir, vol, map_copy
    integer, allocatable :: states(:)
    integer :: conf(PYR_NSTATES,PYR_NSTATES), major(PYR_NSTATES)
    integer :: i, s, t, nassigned, nvols, box, ref_box, ndistinct
    real    :: purity, smpd, ref_smpd, cc, cc_best
    inspect_dir = root//'/inspect'
    call simple_mkdir(inspect_dir)
    allocate(proj)
    call proj%read(result_proj)
    ! the assignment against the planted states
    states = proj%os_ptcl3D%get_all_asint('state')
    conf   = 0
    do i = 1,min(size(states), PYR_NPTCLS)
        if( states(i) < 1 .or. states(i) > PYR_NSTATES ) cycle
        t = mod(i - 1, PYR_NSTATES) + 1
        conf(states(i),t) = conf(states(i),t) + 1
    enddo
    nassigned = sum(conf)
    major     = maxloc(conf, dim=2)
    purity    = 0.
    if( nassigned > 0 ) purity = real(sum(maxval(conf, dim=2))) / real(nassigned)
    call gate%report('pyramid: fraction of the particles assigned', real(nassigned) / real(PYR_NPTCLS))
    call gate%metric('pyramid: assignment purity against the planted states', purity, PYR_PURITY_MIN, &
        &purity >= PYR_PURITY_MIN)
    ndistinct = 0
    do t = 1,PYR_NSTATES
        if( any(major == t .and. sum(conf, dim=2) > 0) ) ndistinct = ndistinct + 1
    enddo
    call gate%report('pyramid: distinct majority truths of the delivered states', real(ndistinct))
    ! each reference against its closest truth, docked in both hands
    if( read_refs_record(exec_dir//'/'//REFS_DIR, PYR_NSTATES, ref_vols, ref_box, ref_smpd) )then
        call gate%report('pyramid: reference volumes', real(size(ref_vols)))
        do s = 1,size(ref_vols)
            map_copy = inspect_dir//'/ref_state'//int2str(s)//'.mrc'
            call simple_copy_file(ref_vols(s), map_copy)
            cc_best = -1.
            do t = 1,size(truths)
                cc_best = max(cc_best, truth_docked_cc(truths(t), map_copy, inspect_dir, 'ref'//int2str(s)//'_truth'//int2str(t)))
            enddo
            call gate%metric('pyramid: reference '//int2str(s)//' docked onto its closest truth', cc_best, PYR_REF_CC_MIN, &
                &cc_best >= PYR_REF_CC_MIN)
        enddo
    else
        call gate%check('pyramid: the references can be read', .false.)
    endif
    ! each refined state's map against the majority truth of its particles
    nvols = 0
    do s = 1,PYR_NSTATES
        if( .not. proj%isthere_in_osout('vol', s) ) cycle
        call proj%get_vol('vol', s, vol, smpd, box)
        if( .not. file_exists(vol) ) cycle
        nvols    = nvols + 1
        map_copy = inspect_dir//'/refined_state'//int2str(s)//'.mrc'
        call simple_copy_file(vol, map_copy)
        cc = truth_docked_cc(truths(major(s)), map_copy, inspect_dir, 'map'//int2str(s))
        call gate%metric('pyramid: refined state '//int2str(s)//' docked onto its majority truth', cc, PYR_MAP_CC_MIN, &
            &cc >= PYR_MAP_CC_MIN)
    enddo
    call gate%check('pyramid: at least two refined states are registered', nvols >= 2)
    call proj%kill
    deallocate(proj)
end subroutine inspect_pyramid

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
