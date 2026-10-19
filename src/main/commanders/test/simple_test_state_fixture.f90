!@descr: the multi-state test fixture: truth maps of blob phantoms or beta-galactosidase states, one simulate_particles mixture, a project per mode, cached per build
! Families: blob (box 64, a common body and a lobe per state) and betagal (1JYX without its hetero
! atoms, box 96 at 2.2 A: compositional states tetramer, tetramer less chain D and the A-B dimer, or
! the tetramer and the tetramer with chain D rigidly rotated). Modes: search (no poses, no states in
! the project), known (the true poses, no states), views (exact repeated views, class = view, or
! view and state). The truth (states, poses, shifts, CTF) is the simulator's orientation document
! truth.txt beside it.
! The project is written into the caller's working directory; the cache holds the rest.
module simple_test_state_fixture
use simple_core_module_api
use simple_image,  only: image
use simple_oris,   only: oris
use simple_cmdline, only: cmdline
implicit none

public :: state_fixture, FIX_SEARCH, FIX_KNOWN, FIX_VIEWS, betagal_frame, point_group_check
private
#include "simple_local_flags.inc"

integer, parameter :: FIX_SEARCH = 1, FIX_KNOWN = 2, FIX_VIEWS = 3
integer, parameter :: FIXTURE_VERSION = 2 !< part of the cache key: raise it when the fixture or simulator changes its output
character(len=6), parameter :: MODE_NAMES(3) = ['search', 'known ', 'views ']
integer, parameter :: BLOB_BOX = 64, BETAGAL_BOX = 96
real,    parameter :: BLOB_SMPD = 2.0, BETAGAL_SMPD = 2.2
real,    parameter :: BLOB_MSKDIAM = 100., BETAGAL_MSKDIAM = 200.
real,    parameter :: KV = 300., CS = 2.7, FRACA = 0.1

!> what to build; dir, the truths and the point groups are set by build
type :: state_fixture
    character(len=16)   :: family     = 'blob'     !< blob | betagal
    character(len=16)   :: kind       = 'distinct' !< blob: distinct | identical; betagal: compositional | conformational | identical
    integer             :: mode       = FIX_SEARCH
    integer             :: nstates    = 3
    real                :: pops(8)    = 0.         !< population fractions (0: equal)
    integer             :: nptcls     = 600
    integer             :: nviews     = 0          !< views mode: distinct directions
    logical             :: state_classes = .false. !< views mode: the 2D class is the view and the true state
    real                :: snr        = 1.
    real                :: defocus    = 1.5        !< microns
    real                :: dferr      = 0.5        !< defocus spread (microns)
    real                :: astigerr   = 0.1
    real                :: sherr      = 2.5        !< shift range (pixels) of the search and known modes
    real                :: conf_angle = 10.        !< betagal conformational: rotation of chain D (deg)
    integer             :: seed       = 20261017
    integer             :: nthr       = 1
    character(len=1024) :: dir        = ''         !< the cache directory, set by build
    character(len=4)    :: pgrps(8)   = ''         !< the declared point group of each truth
    logical             :: reused     = .false.    !< build found the dataset in the cache
  contains
    procedure :: build
    procedure :: truth_fname
    procedure :: project_fname
    procedure :: truth_doc
    procedure :: stack_fname
    procedure :: box
    procedure :: smpd
    procedure :: mskdiam
    procedure :: write_lowpassed
    procedure, private :: key_text
end type state_fixture

contains

    integer function box( self )
        class(state_fixture), intent(in) :: self
        box = merge(BETAGAL_BOX, BLOB_BOX, trim(self%family) == 'betagal')
    end function box

    real function smpd( self )
        class(state_fixture), intent(in) :: self
        smpd = merge(BETAGAL_SMPD, BLOB_SMPD, trim(self%family) == 'betagal')
    end function smpd

    real function mskdiam( self )
        class(state_fixture), intent(in) :: self
        mskdiam = merge(BETAGAL_MSKDIAM, BLOB_MSKDIAM, trim(self%family) == 'betagal')
    end function mskdiam

    function truth_fname( self, k ) result( fname )
        class(state_fixture), intent(in) :: self
        integer,              intent(in) :: k
        type(string) :: fname
        fname = trim(self%dir)//'/truth_state'//int2str(k)//'.mrc'
    end function truth_fname

    function project_fname( self ) result( fname )
        class(state_fixture), intent(in) :: self
        type(string) :: fname
        type(string) :: cwd
        call simple_getcwd(cwd)
        fname = cwd//'/fixture_'//trim(self%family)//'_'//trim(MODE_NAMES(self%mode))//'.simple'
    end function project_fname

    function truth_doc( self ) result( fname )
        class(state_fixture), intent(in) :: self
        type(string) :: fname
        fname = trim(self%dir)//'/truth.txt'
    end function truth_doc

    function stack_fname( self ) result( fname )
        class(state_fixture), intent(in) :: self
        type(string) :: fname
        fname = trim(self%dir)//'/particles.mrcs'
    end function stack_fname

    !> The text that keys the cache: the build's commit, every parameter of the dataset, the seed with
    !! its offset and whether the dataset is one of exact views (the search and known modes share a
    !! dataset)
    function key_text( self ) result( key )
        use simple_test_utils, only: offset_seed
        use simple_defs,       only: SIMPLE_GIT_VERSION
        class(state_fixture), intent(in) :: self
        character(len=:), allocatable :: key
        character(len=512) :: buf
        integer :: k
        write(buf,'(a,i0,1x,a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,5(1x,es12.5),1x,a,es12.5,1x,a,i0,1x,a,i0,1x,a)') &
            &'version=', FIXTURE_VERSION, trim(self%family), trim(self%kind), 'nstates=', self%nstates, 'nptcls=', self%nptcls, &
            &'nviews=', merge(self%nviews, 0, self%mode == FIX_VIEWS), 'snr/defocus/dferr/astigerr/sherr=', self%snr, &
            &self%defocus, self%dferr, self%astigerr, merge(self%sherr, 0., self%mode /= FIX_VIEWS), &
            &'conf_angle=', self%conf_angle, 'seed=', offset_seed(self%seed), 'box=', self%box(), &
            &'pops='
        key = 'commit='//trim(SIMPLE_GIT_VERSION)//' '//trim(buf)
        do k = 1, self%nstates
            write(buf,'(es12.5)') self%pops(k)
            key = key//trim(adjustl(buf))//','
        enddo
    end function key_text

    !> Builds the dataset, or finds it in the cache <SIMPLE_PATH>/test_runs/fixture_cache/<family>_<hash>,
    !! where hash is the FNV-1a hash of key_text (the cache directory holds key.txt), in a temporary
    !! directory renamed into place; then writes the project of the mode into the working directory.
    !! The process seed is left at the front door's either way: a test sets its own seeds after build.
    subroutine build( self )
        use simple_test_utils, only: fix_test_seed
        class(state_fixture), intent(inout) :: self
        character(len=:), allocatable :: key, root, tmp
        character(len=4096) :: path
        type(string) :: cwd, keystr, hash
        integer :: status
        logical :: l_exists
        if( self%nstates < 1 .or. self%nstates > 8 ) THROW_HARD('a fixture holds 1 to 8 states')
        if( trim(self%family) == 'blob' .and. self%nstates > 3 ) THROW_HARD('the blob family has three states')
        if( self%mode == FIX_VIEWS .and. (self%nviews < 1 .or. mod(self%nptcls, max(1,self%nviews)) /= 0) )&
            &THROW_HARD('the views mode needs nptcls a multiple of nviews')
        key = self%key_text()
        call get_environment_variable('SIMPLE_PATH', path, status=status)
        if( status /= 0 ) THROW_HARD('the fixture cache needs SIMPLE_PATH')
        root = trim(path)//'/test_runs/fixture_cache'
        call simple_mkdir(root)
        keystr   = string(key)
        hash     = keystr%to_fnv1a_hash64()
        self%dir = root//'/'//trim(self%family)//'_'//hash%to_char()
        inquire(file=trim(self%dir)//'/COMPLETE', exist=l_exists)
        self%reused = l_exists
        if( .not. l_exists )then
            call simple_getcwd(cwd)
            tmp = trim(self%dir)//'.tmp'//int2str(get_process_id())
            if( dir_exists(tmp) ) call simple_rmdir(tmp)
            call simple_mkdir(tmp)
            call simple_chdir(tmp)
            call write_text('key.txt', key)
            call make_truths(self)
            call simulate(self)
            call write_text('COMPLETE', 'complete')
            call simple_chdir(cwd)
            ! another process may have finished the same dataset first
            call rename(tmp, trim(self%dir), status)
            if( status /= 0 ) call simple_rmdir(tmp)
            inquire(file=trim(self%dir)//'/COMPLETE', exist=l_exists)
            if( .not. l_exists ) THROW_HARD('the fixture was not written: '//trim(self%dir))
        endif
        call fix_test_seed
        call read_pgrps(self)
        call write_project(self)
    end subroutine build

    !> the truth of state k low-pass filtered to lp (A), written to fname
    subroutine write_lowpassed( self, k, lp, fname )
        class(state_fixture), intent(in) :: self
        integer,              intent(in) :: k
        real,                 intent(in) :: lp
        class(string),        intent(in) :: fname
        type(image) :: vol
        call vol%new([self%box(), self%box(), self%box()], self%smpd())
        call vol%read(self%truth_fname(k))
        call vol%bp(0., lp)
        call vol%ifft
        call vol%write(fname)
        call vol%kill
    end subroutine write_lowpassed

    ! ---- truths ----------------------------------------------------------------------------------

    !> truth_state<k>.mrc and pgrps.txt in the working directory
    subroutine make_truths( self )
        use simple_flex_pca_application_tester, only: write_phantom_truths
        class(state_fixture), intent(inout) :: self
        type(string), allocatable :: truths(:)
        type(string) :: tmean, tdiff, cwd
        integer :: k
        select case(trim(self%family))
            case('blob')
                call simple_getcwd(cwd)
                allocate(truths(self%nstates))
                call write_phantom_truths(cwd, self%nstates, truths, tmean, tdiff)
                do k = 1, self%nstates
                    if( trim(self%kind) == 'identical' )then
                        call simple_copy_file(truths(1), string('truth_state'//int2str(k)//'.mrc'))
                    else
                        call simple_copy_file(truths(k), string('truth_state'//int2str(k)//'.mrc'))
                    endif
                    self%pgrps(k) = 'c1'
                enddo
                do k = 1, self%nstates
                    call del_file(truths(k))
                enddo
                call del_file(tmean)
                call del_file(tdiff)
            case('betagal')
                call make_betagal_truths(self)
            case default
                THROW_HARD('no fixture family '//trim(self%family))
        end select
        call write_pgrps(self)
    end subroutine make_truths

    !> The D2 frame of 1JYX from the centres of its four chains (C-alpha atoms of the residues all
    !! four chains hold): the centre o and the rotation f whose rows are the 2-fold axes, z relating
    !! chains A and B, x relating A and C (y = z cross x). A 2-fold axis through o carries the
    !! midpoint of a point and its image, so it runs along the sum of two related chain centres.
    subroutine betagal_frame( mol, o, f )
        use simple_molecule_data, only: molecule_data
        type(molecule_data), intent(in)  :: mol
        real,                intent(out) :: o(3), f(3,3)
        character(len=1), parameter :: CHAINS(4) = ['A','B','C','D']
        logical, allocatable :: ca(:,:)
        real    :: c(3,4), x(3), y(3), z(3)
        integer :: ichain, ires
        logical :: common
        allocate(ca(mol%n, 4))
        do ichain = 1, 4
            ca(:,ichain) = adjustl(mol%name) == 'CA' .and. .not. mol%het .and. mol%chain == CHAINS(ichain)
        enddo
        ! residues present in all four chains
        do ires = 1, mol%n
            if( .not. ca(ires,1) ) cycle
            common = all([(any(ca(:,ichain) .and. mol%resnum == mol%resnum(ires)), ichain = 2, 4)])
            if( .not. common ) ca(ires,1) = .false.
        enddo
        do ichain = 2, 4
            do ires = 1, mol%n
                if( ca(ires,ichain) ) ca(ires,ichain) = any(ca(:,1) .and. mol%resnum == mol%resnum(ires))
            enddo
        enddo
        do ichain = 1, 4
            c(:,ichain) = [(sum(mol%xyz(:,ires), mask=ca(:,ichain)), ires = 1, 3)] / real(count(ca(:,ichain)))
        enddo
        o = sum(c, dim=2) / 4.
        if( sqrt(sum((c(:,1) + c(:,2) - 2.*o)**2)) < 1. .or. sqrt(sum((c(:,1) + c(:,3) - 2.*o)**2)) < 1. )&
            &THROW_HARD('betagal_frame: a chain pair has no 2-fold axis along its centres')
        z = unit(c(:,1) + c(:,2) - 2.*o)
        x = c(:,1) + c(:,3) - 2.*o
        x = unit(x - dot_product(x, z) * z)
        y = [z(2)*x(3) - z(3)*x(2), z(3)*x(1) - z(1)*x(3), z(1)*x(2) - z(2)*x(1)]
        f(1,:) = x
        f(2,:) = y
        f(3,:) = z

      contains

        function unit( v ) result( u )
            real, intent(in) :: v(3)
            real :: u(3)
            u = v / sqrt(sum(v**2))
        end function unit

    end subroutine betagal_frame

    !> the beta-galactosidase states in one frame: the D2 axes on x, y, z, the tetramer's centre on
    !! the box centre (its density centre checked and corrected once), then each state's atoms
    subroutine make_betagal_truths( self )
        use simple_molecule_data, only: molecule_data, betagal_1jyx
        use simple_atoms,         only: atoms
        class(state_fixture), intent(inout) :: self
        type(molecule_data) :: mol, sub
        type(atoms)  :: molecule
        type(image)  :: vol
        real    :: o(3), f(3,3), ctr(3), com(3), axis(3), r(3,3), dc(3), extent(3)
        integer :: k, i, ilong
        logical, allocatable :: keep(:), dmask(:)
        mol = betagal_1jyx()
        keep = .not. mol%het
        mol  = subset(mol, keep)
        call betagal_frame(mol, o, f)
        ctr = real(BETAGAL_BOX / 2) * BETAGAL_SMPD
        do i = 1, mol%n
            mol%xyz(i,:) = matmul(f, mol%xyz(i,:) - o) + ctr
        enddo
        ! the density centre of the tetramer, which is its symmetry centre, onto the box centre
        call assert_in_frame(mol)
        call molecule%pdb2mrc(volfile=string('tetramer_probe.mrc'), smpd=BETAGAL_SMPD, mol=mol, &
            &vol_dim=[BETAGAL_BOX,BETAGAL_BOX,BETAGAL_BOX])
        call vol%new([BETAGAL_BOX,BETAGAL_BOX,BETAGAL_BOX], BETAGAL_SMPD)
        call vol%read(string('tetramer_probe.mrc'))
        call vol%masscen(com)
        call vol%kill
        call del_file('tetramer_probe.mrc')
        do i = 1, mol%n
            mol%xyz(i,:) = mol%xyz(i,:) - com * BETAGAL_SMPD
        enddo
        ! the long axis of the tetramer
        do i = 1, 3
            extent(i) = maxval(mol%xyz(:,i)) - minval(mol%xyz(:,i))
        enddo
        ilong = maxloc(extent, dim=1)
        do k = 1, self%nstates
            select case(trim(self%kind))
                case('compositional')
                    select case(k)
                        case(1); sub = mol;                                   self%pgrps(k) = 'd2'
                        case(2); sub = subset(mol, mol%chain /= 'D');         self%pgrps(k) = 'c1'
                        case(3); sub = subset(mol, mol%chain == 'A' .or. mol%chain == 'B'); self%pgrps(k) = 'c2'
                        case default; THROW_HARD('the compositional betagal states are three')
                    end select
                case('conformational')
                    sub = mol
                    self%pgrps(k) = merge('d2', 'c1', k == 1)
                    if( k > 1 )then
                        ! chain D rotated by (k-1) times conf_angle about the axis through its centre
                        ! parallel to the tetramer's long axis
                        dmask = sub%chain == 'D'
                        dc    = [(sum(sub%xyz(:,i), mask=dmask), i = 1, 3)] / real(count(dmask))
                        axis  = 0.
                        axis(ilong) = 1.
                        r = rotation(axis, real(k-1) * self%conf_angle)
                        do i = 1, sub%n
                            if( dmask(i) ) sub%xyz(i,:) = matmul(r, sub%xyz(i,:) - dc) + dc
                        enddo
                    endif
                case('identical')
                    sub = mol
                    self%pgrps(k) = 'd2'
                case default
                    THROW_HARD('no betagal kind '//trim(self%kind))
            end select
            call assert_in_frame(sub)
            call molecule%pdb2mrc(volfile=string('truth_state'//int2str(k)//'.mrc'), smpd=BETAGAL_SMPD, mol=sub, &
                &vol_dim=[BETAGAL_BOX,BETAGAL_BOX,BETAGAL_BOX])
        enddo
        call molecule%kill

      contains

        !> pdb2mrc centres a model with a negative coordinate on its own centre, which would move this
        !! state out of the frame the states share
        subroutine assert_in_frame( m )
            type(molecule_data), intent(in) :: m
            if( any(m%xyz(1:m%n,:) < 0.) ) THROW_HARD('a beta-galactosidase state has a negative coordinate')
        end subroutine assert_in_frame

        function rotation( u, ang ) result( rm )
            real, intent(in) :: u(3), ang
            real :: rm(3,3), c, s, kk(3,3)
            integer :: j
            c  = cos(deg2rad(ang))
            s  = sin(deg2rad(ang))
            kk = reshape([0., u(3), -u(2), -u(3), 0., u(1), u(2), -u(1), 0.], [3,3])
            rm = s * kk + (1. - c) * matmul(kk, kk)
            do j = 1, 3
                rm(j,j) = rm(j,j) + 1.
            enddo
        end function rotation

    end subroutine make_betagal_truths

    !> the atoms of mol where keep holds
    function subset( mol, keep ) result( sub )
        use simple_molecule_data, only: molecule_data
        type(molecule_data), intent(in) :: mol
        logical,             intent(in) :: keep(:)
        type(molecule_data) :: sub
        integer :: j
        sub%n = count(keep)
        sub%name      = pack(mol%name, keep)
        sub%resname   = pack(mol%resname, keep)
        sub%chain     = pack(mol%chain, keep)
        sub%element   = pack(mol%element, keep)
        sub%num       = pack(mol%num, keep)
        sub%resnum    = pack(mol%resnum, keep)
        sub%occupancy = pack(mol%occupancy, keep)
        sub%beta      = pack(mol%beta, keep)
        sub%het       = pack(mol%het, keep)
        allocate(sub%xyz(sub%n, 3))
        do j = 1, 3
            sub%xyz(:,j) = pack(mol%xyz(:,j), keep)
        enddo
    end function subset

    ! ---- simulation and projects -----------------------------------------------------------------

    !> one simulate_particles call with every state: random poses and shifts within sherr, or exact
    !! views; the states drawn by the program with the fixture's populations
    subroutine simulate( self )
        use simple_commanders_sim, only: commander_simulate_particles
        use simple_test_utils,     only: set_fixed_seed
        class(state_fixture), intent(in) :: self
        type(commander_simulate_particles) :: xsim
        class(cmdline), allocatable :: cl
        type(oris)    :: views, poses
        character(len=256) :: pops
        character(len=16)  :: buf
        integer :: k, i, nrep
        allocate(cl)
        do k = 1, self%nstates
            call cl%set('vol'//int2str(k), 'truth_state'//int2str(k)//'.mrc')
        enddo
        if( self%nstates > 1 )then
            call cl%set('nstates', self%nstates)
            if( any(self%pops(1:self%nstates) > 0.) )then
                pops = ''
                do k = 1, self%nstates
                    write(buf,'(f8.5)') self%pops(k)
                    pops = trim(pops)//trim(adjustl(buf))
                    if( k < self%nstates ) pops = trim(pops)//','
                enddo
                call cl%set('state_pops', trim(pops))
            endif
        endif
        call cl%set('prg',      'simulate_particles')
        call cl%set('mkdir',    'no')
        call cl%set('smpd',     self%smpd())
        call cl%set('mskdiam',  self%mskdiam())
        call cl%set('nptcls',   self%nptcls)
        call cl%set('pgrp',     'c1')
        call cl%set('ctf',      'yes')
        call cl%set('snr',      self%snr)
        call cl%set('kv',       KV)
        call cl%set('cs',       CS)
        call cl%set('fraca',    FRACA)
        call cl%set('defocus',  self%defocus)
        call cl%set('dferr',    self%dferr)
        call cl%set('astigerr', self%astigerr)
        call cl%set('bfac',     0.)
        call cl%set('bfacerr',  0.)
        call cl%set('outstk',   'particles.mrcs')
        call cl%set('outfile',  'truth.txt')
        call cl%set('nthr',     self%nthr)
        if( self%mode == FIX_VIEWS )then
            ! nviews spiral directions, each repeated, without a state field: the program draws them
            call views%new(self%nviews, is_ptcl=.false.)
            call views%spiral
            nrep = self%nptcls / self%nviews
            call poses%new(self%nptcls, is_ptcl=.false.)
            do i = 1, self%nptcls
                call poses%set_euler(i, views%get_euler((i - 1) / nrep + 1))
            enddo
            call poses%write(string('views.txt'), [1, self%nptcls])
            call cl%set('oritab', 'views.txt')
            call cl%set('sherr', 0.)
            call views%kill
            call poses%kill
        else
            call cl%set('sherr', self%sherr)
        endif
        call set_fixed_seed(self%seed, propagate=.true.)
        call xsim%execute(cl)
        call cl%kill
        deallocate(cl)
    end subroutine simulate

    !> the project of the mode: the stack with each particle's CTF; known and views: the true poses
    !! (and shifts); views: the 2D class is the view (or, with state_classes, the view and the true
    !! state), scored 1 as a finished 2D solution; every particle in state 1
    subroutine write_project( self )
        use simple_sp_project, only: sp_project
        use simple_type_defs,  only: CTFFLAG_YES, ctfparams
        class(state_fixture), intent(in) :: self
        class(sp_project), allocatable :: proj
        type(oris)      :: truth
        type(ctfparams) :: ctf
        type(cmdline)   :: cline_compenv
        integer :: i, nrep
        allocate(proj)
        call truth%new(self%nptcls, is_ptcl=.true.)
        call truth%read(self%truth_doc(), [1, self%nptcls])
        ctf%smpd    = self%smpd()
        ctf%kv      = KV
        ctf%cs      = CS
        ctf%fraca   = FRACA
        ctf%dfx     = self%defocus
        ctf%dfy     = self%defocus
        ctf%ctfflag = CTFFLAG_YES
        call proj%add_stk(self%stack_fname(), ctf)
        nrep = 1
        if( self%mode == FIX_VIEWS ) nrep = self%nptcls / self%nviews
        do i = 1, self%nptcls
            call proj%os_ptcl2D%set_dfx(i, truth%get_dfx(i))
            call proj%os_ptcl2D%set_dfy(i, truth%get_dfy(i))
            call proj%os_ptcl2D%set(i, 'angast', truth%get(i, 'angast'))
            call proj%os_ptcl3D%set_dfx(i, truth%get_dfx(i))
            call proj%os_ptcl3D%set_dfy(i, truth%get_dfy(i))
            call proj%os_ptcl3D%set(i, 'angast', truth%get(i, 'angast'))
            call proj%os_ptcl2D%set_state(i, 1)
            call proj%os_ptcl3D%set_state(i, 1)
            call proj%os_ptcl2D%set(i, 'eo', real(mod(i, 2)))
            call proj%os_ptcl3D%set(i, 'eo', real(mod(i, 2)))
            if( self%mode /= FIX_SEARCH )then
                call proj%os_ptcl3D%set_euler(i, truth%get_euler(i))
                call proj%os_ptcl3D%set_shift(i, truth%get_2Dshift(i))
                call proj%os_ptcl2D%set(i, 'e3', truth%get(i, 'e3'))
                call proj%os_ptcl2D%set_shift(i, truth%get_2Dshift(i))
                call proj%os_ptcl3D%set(i, 'proj', 1.)
            endif
            if( self%mode == FIX_VIEWS )then
                ! a 2D solution: a particle with neither angles nor a score counts as never searched
                call proj%os_ptcl2D%set(i, 'corr', 1.)
                if( self%state_classes )then
                    ! an oracle 2D classification, pure in state: class (view - 1) * nstates + state
                    call proj%os_ptcl2D%set_class(i, ((i - 1) / nrep) * self%nstates + nint(truth%get(i, 'state')))
                else
                    call proj%os_ptcl2D%set_class(i, (i - 1) / nrep + 1)
                endif
            endif
        enddo
        call proj%update_projinfo(self%project_fname())
        call proj%update_compenv(cline_compenv)
        call proj%write(self%project_fname())
        call proj%kill
        call truth%kill
        deallocate(proj)
    end subroutine write_project

    subroutine write_pgrps( self )
        class(state_fixture), intent(in) :: self
        integer :: funit, k
        open(newunit=funit, file='pgrps.txt', status='replace', action='write')
        do k = 1, self%nstates
            write(funit,'(a)') trim(self%pgrps(k))
        enddo
        close(funit)
    end subroutine write_pgrps

    subroutine read_pgrps( self )
        class(state_fixture), intent(inout) :: self
        integer :: funit, k
        open(newunit=funit, file=trim(self%dir)//'/pgrps.txt', status='old', action='read')
        do k = 1, self%nstates
            read(funit,'(a)') self%pgrps(k)
        enddo
        close(funit)
    end subroutine read_pgrps

    ! ---- checks and helpers ----------------------------------------------------------------------

    !> The self-correlation of the map in fname under the operators ops (rotvol, about the box
    !! centre): the smallest of them in cc_group; and under the extra operator of the next larger
    !! group of its family in cc_larger
    subroutine point_group_check( fname, ops, larger, cc_group, cc_larger )
        use simple_simple_volinterp, only: rotvol
        use simple_ori,              only: ori
        class(string), intent(in)  :: fname
        real,          intent(in)  :: ops(:,:,:), larger(3,3)
        real,          intent(out) :: cc_group, cc_larger
        type(image) :: vol, rot
        type(ori)   :: o
        integer :: ldim(3), n, g
        call find_ldim_nptcls(fname, ldim, n)
        call vol%new(ldim, find_img_smpd(fname))
        call vol%read(fname)
        call o%new(.false.)
        cc_group = 1.
        do g = 1, size(ops, 3)
            call o%set_euler(m2euler(ops(:,:,g)))
            rot = rotvol(vol, o)
            cc_group = min(cc_group, rot%real_corr(vol))
            call rot%kill
        enddo
        call o%set_euler(m2euler(larger))
        rot = rotvol(vol, o)
        cc_larger = rot%real_corr(vol)
        call rot%kill
        call vol%kill
        call o%kill
    end subroutine point_group_check

    subroutine write_text( fname, text )
        character(len=*), intent(in) :: fname, text
        integer :: funit
        open(newunit=funit, file=fname, status='replace', action='write')
        write(funit,'(a)') text
        close(funit)
    end subroutine write_text

end module simple_test_state_fixture
