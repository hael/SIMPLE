!@descr: center_ptcls in miniature on the multi-state fixture's beta-galactosidase particles with planted shifts, a gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_center
implicit none
#include "simple_local_flags.inc"

contains

    !> center_ptcls (doc/algorithms/center_ptcls.md) on 2 000 simulated beta-galactosidase particles
    !! of the multi-state test fixture (one state, SNR 0.5) whose project holds no pose and no shift, the simulator
    !! having planted shifts within SHERR_PX. The gate asserts the program's contract: every active particle carries
    !! the same shift and state of centring in both particle fields (a measured width, or the sentinel where the
    !! pass's draw carried no width and left it unmoved), no other parameter changes, every shift lies within the
    !! disk, and the centring improves on the planted error. It reports the residual against the
    !! planted shifts and the cost; a profile sets floors later. Every row goes to metrics.tsv.
    module subroutine exec_test_center_ptcls( self, cline )
        use simple_test_state_fixture,      only: state_fixture, FIX_SEARCH
        use simple_test_gate,               only: test_gate
        use simple_commanders_center_ptcls, only: commander_center_ptcls
        use simple_math,                    only: hpsort
        use simple_timer,                   only: tic, toc, timer_int_kind
        use, intrinsic :: ieee_arithmetic,  only: ieee_is_finite
        class(commander_test_center_ptcls), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        integer, parameter :: NPTCLS      = 2000
        real,    parameter :: SNR         = 0.5
        real,    parameter :: SHERR_PX    = 6.    !< the planted shifts' range (pixels at 2.2 A): 13 A, the picking error's scale
        real,    parameter :: SHIFT_SIGMA = 10.   !< center_ptcls' default prior (A), whose 3-sigma disk bounds every shift
        class(parameters), allocatable :: params
        type(state_fixture)          :: fix
        type(test_gate)              :: gate
        type(sp_project)             :: proj
        type(oris)                   :: truth
        type(commander_center_ptcls) :: xcenter
        type(cmdline)                :: cl
        type(string)                 :: cwd_saved, root, projfile
        integer(timer_int_kind)      :: t0
        real,    allocatable :: e3_before(:), cls_before(:), dx(:), dy(:), rplant(:), rres(:), sorted(:)
        integer, allocatable :: inds(:)
        real    :: smpd, wall, mx, my, sh(2), sh3(2), tsh(2), bound_px, med_plant, med_res
        integer :: i, status, n3, n5, nsent
        logical :: l_centred, l_same, l_bound
        allocate(params)
        call simple_getcwd(cwd_saved)
        root = filepath(cwd_saved, 'test_center_ptcls_'//int2str(get_process_id()))
        if( dir_exists(root) ) call simple_rmdir(root)
        call simple_mkdir(root)
        call simple_chdir(root, status)
        if( status /= 0 ) THROW_HARD('TEST_CENTER_PTCLS FAILED: could not enter the test directory')
        call params%new(cline)
        call gate%new(filepath(cwd_saved, 'metrics.tsv'))
        ! the fixture: raw picks, the truth's planted shifts beside them
        call fix%new('betagal', FIX_SEARCH, 1, NPTCLS, SNR, params%nthr, kind='identical', sherr=SHERR_PX)
        call fix%build
        smpd     = fix%smpd()
        projfile = simple_abspath(fix%project_fname())
        call truth%new(NPTCLS, is_ptcl=.true.)
        call truth%read(fix%truth_doc())
        call proj%read(projfile)
        if( proj%os_ptcl2D%get_noris() /= NPTCLS ) THROW_HARD('the fixture holds another particle count; test_center_ptcls')
        allocate(e3_before(NPTCLS), cls_before(NPTCLS))
        do i = 1,NPTCLS
            e3_before(i)  = proj%os_ptcl2D%e3get(i)
            cls_before(i) = proj%os_ptcl2D%get(i, 'class')
        enddo
        call proj%kill
        ! the program, in place
        call cl%set('prg',      'center_ptcls')
        call cl%set('projfile', projfile%to_char())
        call cl%set('mskdiam',  fix%mskdiam())
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'no')
        t0 = tic()
        call xcenter%execute(cl)
        wall = real(toc(t0))
        call cl%kill
        ! the contract
        call proj%read(projfile)
        allocate(dx(NPTCLS), dy(NPTCLS), rplant(NPTCLS), rres(NPTCLS))
        bound_px  = 3. * SHIFT_SIGMA / smpd + 2.
        l_centred = .true.
        l_same    = .true.
        l_bound   = .true.
        nsent     = 0
        do i = 1,NPTCLS
            sh  = proj%os_ptcl2D%get_2Dshift(i)
            sh3 = proj%os_ptcl3D%get_2Dshift(i)
            tsh = truth%get_2Dshift(i)
            if( .not. all(ieee_is_finite(sh)) .or. any(abs(sh - sh3) > 1.e-6) ) l_centred = .false.
            if( abs(proj%os_ptcl2D%get(i, 'shsig') - proj%os_ptcl3D%get(i, 'shsig')) > 1.e-6 ) l_centred = .false.
            if( proj%os_ptcl2D%get(i, 'shsig') < 0. )then
                ! not centred by the pass: the sentinel exactly, and the particle left where it was
                nsent = nsent + 1
                if( abs(proj%os_ptcl2D%get(i, 'shsig') + 1.) > 1.e-6 .or. any(abs(sh) > 1.e-6) ) l_centred = .false.
            endif
            if( abs(proj%os_ptcl2D%e3get(i) - e3_before(i)) > 1.e-6 .or. &
                &abs(proj%os_ptcl2D%get(i, 'class') - cls_before(i)) > 1.e-6 ) l_same = .false.
            if( sqrt(sum(sh**2)) > bound_px ) l_bound = .false.
            dx(i)     = sh(1) - tsh(1)
            dy(i)     = sh(2) - tsh(2)
            rplant(i) = sqrt(sum(tsh**2)) * smpd
        enddo
        ! the residual after the common offset (the references' centre is defined up to a translation)
        mx = median_of(dx)
        my = median_of(dy)
        n3 = 0
        n5 = 0
        do i = 1,NPTCLS
            rres(i) = sqrt((dx(i) - mx)**2 + (dy(i) - my)**2) * smpd
            if( rres(i) <= 3. ) n3 = n3 + 1
            if( rres(i) <= 5. ) n5 = n5 + 1
        enddo
        med_plant = median_of(rplant)
        med_res   = median_of(rres)
        call gate%check('every_particle_centred_or_marked_alike_in_both_fields', l_centred)
        call gate%check('no_other_particle_parameter_changed',               l_same)
        call gate%check('every_shift_within_the_disk',                       l_bound)
        call gate%check('centring_improves_on_the_planted_error',            med_res < med_plant)
        call gate%report('planted_shift_median_A',     med_plant)
        call gate%report('residual_median_A',          med_res)
        call gate%report('residual_within_3A_fraction', real(n3) / real(NPTCLS))
        call gate%report('residual_within_5A_fraction', real(n5) / real(NPTCLS))
        call gate%report('common_offset_px',           sqrt(mx**2 + my**2))
        call gate%report('not_centred_fraction',       real(nsent) / real(NPTCLS))
        call gate%report('center_ptcls_wall_s',        wall)
        call proj%kill
        call truth%kill
        call fix%kill
        call simple_chdir(cwd_saved, status)
        if( status /= 0 ) THROW_HARD('TEST_CENTER_PTCLS FAILED: could not restore the original directory')
        if( gate%passed() )then
            call gate%kill
            call simple_rmdir(root)
            write(logfhandle,'(A)') 'PASS: center_ptcls in miniature'
            call simple_end('**** SIMPLE_TEST_CENTER_PTCLS NORMAL STOP ****')
        else
            call gate%kill
            THROW_HARD('TEST_CENTER_PTCLS FAILED')
        endif

      contains

        real function median_of( v )
            real, intent(in) :: v(:)
            integer :: j
            sorted = v
            if( allocated(inds) ) deallocate(inds)
            allocate(inds(size(v)))
            inds = [(j, j = 1,size(v))]
            call hpsort(sorted, inds)
            median_of = sorted(max(1, size(v) / 2))
        end function median_of

    end subroutine exec_test_center_ptcls

end submodule simple_commanders_test_highlevel_center
