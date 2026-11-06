!@descr: the two-state flex_pca phantom shared by testers and high-level tests: truth maps, simulated stack, project, command line
! Two blob states over shared poses (optionally a third, or repeated exact views), simulated with CTF at
! SNR 1; the project, the flex_pca command line and the weight-table reader of the library phantom test.
module simple_test_flex_phantom_fixture
use, intrinsic :: iso_fortran_env, only: real32
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_syslib,             only: simple_abspath, simple_chdir
use simple_cmdline,            only: cmdline
use simple_commanders_sim,     only: commander_simulate_particles
use simple_state_weight_set,   only: state_weight_set
use simple_image,              only: image
use simple_oris,               only: oris
use simple_sp_project,         only: sp_project
use simple_test_truth_metrics, only: add_gaussian_blob
use simple_test_utils,         only: set_fixed_seed
use simple_type_defs,          only: CTFFLAG_YES, ctfparams
implicit none
private
#include "simple_local_flags.inc"

public :: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project, read_state_weight_table
public :: configure_flex_pca_phantom, write_phantom_truths

integer, parameter         :: BOX                             = 64
integer, parameter         :: DEFAULT_NPTCLS                  = 2000
integer, parameter         :: MIN_NEFF                        = 100
integer, parameter         :: SEED_A                          = 202610061
integer, parameter         :: SEED_B                          = 202610062
integer, parameter         :: SEED_C                          = 202610091
real,    parameter         :: SMPD                            = 2.0
real,    parameter         :: MSKDIAM                         = 100.0
real,    parameter         :: KV                              = 300.0
real,    parameter         :: CS                              = 2.7
real,    parameter         :: FRACA                           = 0.1
real,    parameter         :: DEFOCUS                         = 1.5
real,    parameter         :: SIM_SNR                         = 1.0
! The delivered state ids are recoded independently to their majority truth class. With at most
! three state ids there are only 2^3 binary recodings; under random truth labels Hoeffding's union
! bound puts P(best accuracy >= 0.60) below 8*exp(-2*DEFAULT_NPTCLS*0.10^2) < 4e-17. This is an independent
! detection floor, not the phase-0 production observation used later as a regression floor.
real,    parameter         :: LABEL_ACCURACY_MIN              = 0.60
integer, parameter, public :: FLEX_PHANTOM_BOX                = BOX
integer, parameter, public :: FLEX_PHANTOM_NPTCLS             = DEFAULT_NPTCLS
integer, parameter, public :: FLEX_PHANTOM_NCOMP              = 2
integer, parameter, public :: FLEX_PHANTOM_MAX_STATES         = 3
real,    parameter, public :: FLEX_PHANTOM_SMPD               = SMPD
real,    parameter, public :: FLEX_PHANTOM_MSKDIAM            = MSKDIAM
real,    parameter, public :: FLEX_PHANTOM_LABEL_ACCURACY_MIN = LABEL_ACCURACY_MIN
!> file tags of the phantom states
character(len=1),  parameter :: STATE_TAGS(3)        = ['a', 'b', 'c']

contains

    !> Truth maps, shared poses and the interleaved stack (row nstates*(i-1)+s holds pose i in state s); the
    !! defaults are the two-state library fixture. nstates=3 adds state C (truth_c); nviews > 0 repeats
    !! nviews directions, each one pose, so a direction is an exact 2D class
    subroutine create_flex_pca_phantom_fixture( root, nthr, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, &
        &nptcls, nstates, nviews, truth_c )
        class(string),           intent(in)  :: root
        integer,                 intent(in)  :: nthr
        type(string),            intent(out) :: truth_a, truth_b, truth_mean, truth_diff, oritab, stack
        integer,       optional, intent(in)  :: nptcls, nstates, nviews
        type(string),  optional, intent(out) :: truth_c
        integer, parameter :: SEEDS(3) = [SEED_A, SEED_B, SEED_C]
        type(string), allocatable :: truths(:), stks(:)
        integer :: nst, nv, nposes, s
        call phantom_layout(nst, nv, nposes, nptcls, nstates, nviews)
        if( nst == 3 .and. .not. present(truth_c) ) THROW_HARD('three phantom states need truth_c')
        allocate(truths(nst), stks(nst))
        call write_phantom_truths(root, nst, truths, truth_mean, truth_diff)
        truth_a = truths(1)
        truth_b = truths(2)
        if( nst == 3 ) truth_c = truths(3)
        oritab = root//'/truth_orientations.txt'
        call write_uniform_orientations(oritab, nposes, nv)
        call simple_chdir(root)
        do s = 1, nst
            stks(s) = simulate_stack(truths(s), oritab, string('particles_state_'//STATE_TAGS(s)//'.mrcs'), &
                &SEEDS(s), nposes, nthr)
        enddo
        if( nst == 2 )then
            stack = interleave_stacks(stks, nposes, string('particles_two_state.mrcs'))
        else
            stack = interleave_stacks(stks, nposes, string('particles_three_state.mrcs'))
        endif
        deallocate(truths, stks)
    end subroutine create_flex_pca_phantom_fixture

    !> The state count, direction count and poses per state of the phantom, from the optional arguments of
    !! the fixture routines; every state shares the nposes poses, every direction nposes/nviews of them
    subroutine phantom_layout( nst, nv, nposes, nptcls, nstates, nviews )
        integer,           intent(out) :: nst, nv, nposes
        integer, optional, intent(in)  :: nptcls, nstates, nviews
        integer :: ntot
        ntot = DEFAULT_NPTCLS
        nst  = 2
        nv   = 0
        if( present(nptcls)  ) ntot = nptcls
        if( present(nstates) ) nst  = nstates
        if( present(nviews)  ) nv   = nviews
        if( nst < 2 .or. nst > 3 ) THROW_HARD('the phantom has two or three states')
        if( mod(ntot, nst) /= 0 ) THROW_HARD('phantom particles must split evenly over the states')
        nposes = ntot / nst
        if( nv > 0 )then
            if( mod(nposes, nv) /= 0 ) THROW_HARD('phantom poses must split evenly over the directions')
        endif
    end subroutine phantom_layout

    !> Common asymmetric body plus one narrow lobe: moved by exactly eight pixels between A and B, and in
    !! state C at the third corner of a near-equilateral triangle with them
    subroutine write_phantom_truths( root, nst, truths, truth_mean, truth_diff )
        class(string), intent(in)  :: root
        integer,       intent(in)  :: nst
        type(string),  intent(out) :: truths(nst), truth_mean, truth_diff
        real,    parameter :: LOBE_POS(3,3) = reshape([-8.0, -2.0, 0.0,   8.0, -2.0, 0.0,   0.0, 12.0, 0.0], [3,3])
        integer, parameter :: NBLOBS        = 6
        real,    parameter :: POS(3,NBLOBS) = reshape([ &
            &-18., -10.,   4.,   14.,   8.,  -8.,   -6.,  18.,  12., &
            &  4., -16., -14.,   20., -12.,  16.,  -20.,  13., -11.], [3,NBLOBS])
        real,    parameter :: SIGMA(NBLOBS) = [7.0, 6.0, 5.5, 6.5, 4.5, 5.0]
        real,    parameter :: AMP(NBLOBS)   = [1.0, 0.85, 0.75, 0.70, 0.55, 0.50]
        type(image)       :: vol
        real, allocatable :: base(:,:,:), a(:,:,:), b(:,:,:), mean(:,:,:)
        real    :: ctr(3), d2
        integer :: blob, ix, iy, iz, s
        do s = 1, nst
            truths(s) = root//'/truth_state_'//STATE_TAGS(s)//'.mrc'
        enddo
        truth_mean = root//'/truth_mean.mrc'
        truth_diff = root//'/truth_a_minus_b.mrc'
        allocate(base(BOX,BOX,BOX), source=0.0)
        ctr = real([BOX,BOX,BOX]) / 2.0 + 1.0
        do iz = 1, BOX
            do iy = 1, BOX
                do ix = 1, BOX
                    do blob = 1, NBLOBS
                        d2 = sum((([real(ix),real(iy),real(iz)] - ctr) * SMPD - POS(:,blob))**2)
                        base(ix,iy,iz) = base(ix,iy,iz) + AMP(blob) * exp(-0.5 * d2 / SIGMA(blob)**2)
                    enddo
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol%set_rmat(base, .false.)
        do s = 1, nst
            call vol%write(truths(s))
            call add_gaussian_blob(truths(s), LOBE_POS(:,s), 4.0, 0.85)
        enddo
        allocate(mean(BOX,BOX,BOX), source=0.0)
        do s = 1, nst
            call vol%read(truths(s))
            if( s == 1 ) a = vol%get_rmat()
            if( s == 2 ) b = vol%get_rmat()
            mean = mean + vol%get_rmat()
        enddo
        call vol%set_rmat(mean / real(nst), .false.)
        call vol%write(truth_mean)
        ! the difference of A and B, when there is a B (a single state has none)
        if( nst >= 2 )then
            call vol%set_rmat(a - b, .false.)
            call vol%write(truth_diff)
            deallocate(b)
        endif
        call vol%kill
        deallocate(base, a, mean)
    end subroutine write_phantom_truths

    !> nposes uniform spiral poses; with nviews > 0, nviews spiral directions, each repeated as one pose
    !! nposes/nviews times in a row
    subroutine write_uniform_orientations( fname, nposes, nviews )
        class(string), intent(in) :: fname
        integer,       intent(in) :: nposes, nviews
        class(oris), allocatable :: poses, views
        integer :: i, nrep
        allocate(poses)
        call poses%new(nposes, is_ptcl=.true.)
        if( nviews > 0 )then
            allocate(views)
            call views%new(nviews, is_ptcl=.true.)
            call views%spiral()
            nrep = nposes / nviews
            do i = 1, nposes
                call poses%transfer_ori(i, views, (i - 1) / nrep + 1)
            enddo
            call views%kill
            deallocate(views)
        else
            call poses%spiral()
        endif
        call poses%write(fname, [1,nposes])
        call poses%kill
        deallocate(poses)
    end subroutine write_uniform_orientations

    function simulate_stack( volume, oritab, outstk, seed, nposes, nthr ) result( stack )
        class(string), intent(in) :: volume, oritab, outstk
        integer,       intent(in) :: seed, nposes, nthr
        type(string) :: stack
        class(commander_simulate_particles), allocatable :: simulate
        class(cmdline), allocatable :: cline
        allocate(simulate, cline)
        call set_fixed_seed(seed, propagate=.true.)
        call cline%set('prg',      'simulate_particles')
        call cline%set('mkdir',    'no')
        call cline%set('vol1',     volume)
        call cline%set('oritab',   oritab)
        call cline%set('outstk',   outstk)
        call cline%set('outfile',  outstk//'.txt')
        call cline%set('nptcls',   nposes)
        call cline%set('nthr',     nthr)
        call cline%set('smpd',     SMPD)
        call cline%set('mskdiam',  MSKDIAM)
        call cline%set('pgrp',     'c1')
        call cline%set('ctf',      'yes')
        call cline%set('snr',      SIM_SNR)
        call cline%set('kv',       KV)
        call cline%set('cs',       CS)
        call cline%set('fraca',    FRACA)
        call cline%set('defocus',  DEFOCUS)
        call cline%set('dferr',    0.0)
        call cline%set('astigerr', 0.0)
        call cline%set('bfac',     0.0)
        call cline%set('bfacerr',  0.0)
        call cline%set('sherr',    0.0)
        call simulate%execute(cline)
        call cline%kill
        stack = simple_abspath(outstk)
        deallocate(simulate, cline)
    end function simulate_stack

    !> image i of state stack s goes to row nstates*(i-1)+s
    function interleave_stacks( stks, nposes, outstk ) result( stack )
        class(string), intent(in) :: stks(:), outstk
        integer,       intent(in) :: nposes
        type(string) :: stack
        type(image)  :: img
        integer :: i, s, nst
        nst = size(stks)
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        do i = 1, nposes
            do s = 1, nst
                call img%read(stks(s), i)
                call img%write(outstk, nst*(i - 1) + s)
            enddo
        enddo
        call img%kill
        stack = simple_abspath(outstk)
    end function interleave_stacks

    !> Every truth state occurs in both even/odd halfsets. The rows of one pose (2i-1 and 2i by default)
    !! share it, so orientation coverage cannot become a state cue. With nviews > 0 the ptcl2D class of a
    !! particle is its direction, an exact 2D class: one pose, no in-plane rotation and no shift.
    subroutine build_flex_pca_phantom_project( projfile, stack, oritab, nptcls, nstates, nviews )
        class(string),     intent(in) :: projfile, stack, oritab
        integer, optional, intent(in) :: nptcls, nstates, nviews
        class(sp_project), allocatable :: project
        class(oris),       allocatable :: poses
        type(ctfparams) :: ctf
        class(cmdline), allocatable:: cline_compenv
        integer :: i, ipose, eo, nst, nv, nposes
        allocate(cline_compenv)
        call phantom_layout(nst, nv, nposes, nptcls, nstates, nviews)
        allocate(project, poses)
        call poses%new(nposes, is_ptcl=.true.)
        call poses%read(oritab, [1,nposes])
        ctf%smpd    = SMPD
        ctf%kv      = KV
        ctf%cs      = CS
        ctf%fraca   = FRACA
        ctf%dfx     = DEFOCUS
        ctf%dfy     = DEFOCUS
        ctf%angast  = 0.0
        ctf%ctfflag = CTFFLAG_YES
        call project%add_stk(stack, ctf)
        do i = 1, nst * nposes
            ipose = (i - 1) / nst + 1
            eo     = mod(ipose, 2)
            call project%os_ptcl3D%set_euler(i, poses%get_euler(ipose))
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo',   real(eo))
            call project%os_ptcl3D%set(i, 'proj', 1.0)
            call project%os_ptcl2D%set_euler(i, poses%get_euler(ipose))
            call project%os_ptcl2D%set_state(i, 1)
            call project%os_ptcl2D%set(i, 'eo', real(eo))
            if( nv > 0 ) call project%os_ptcl2D%set_class(i, (ipose - 1) / (nposes / nv) + 1)
        enddo
        call project%update_projinfo(projfile)
        ! the computing environment new_project writes; the distributed master's queue setup reads it
        call project%update_compenv(cline_compenv)
        call project%write(projfile)
        call poses%kill
        call project%kill
        deallocate(project, poses)
    end subroutine build_flex_pca_phantom_project

    !> The project's state weight set as a table: weights(nstates,nptcls) and the hard labels(nptcls)
    subroutine read_state_weight_table( project, field, nstates, weights, labels, status, message )
        class(sp_project),         intent(in)  :: project
        class(oris),               intent(in)  :: field
        integer,                   intent(out) :: nstates
        real(real32), allocatable, intent(out) :: weights(:,:)
        integer,      allocatable, intent(out) :: labels(:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        type(state_weight_set)    :: weight_set
        real(real32), allocatable :: column(:)
        integer :: s
        nstates = 0
        call weight_set%new(project, field, status, message)
        if( status /= 0 ) return
        nstates = weight_set%get_nstates()
        allocate(weights(nstates, weight_set%get_nptcls()))
        do s = 1, nstates
            call weight_set%get_weights(s, column)
            weights(s,:) = column
        enddo
        call weight_set%get_labels(labels)
        call weight_set%kill
    end subroutine read_state_weight_table

    !> The common production command line. nparts=1 selects shared memory; nparts>1 selects
    !! the distributed master and its real worker processes through the local queue system.
    subroutine configure_flex_pca_phantom( cline, truth_mean, basis_backend, state_backend, nparts, nthr, nbins, box_crop )
        class(cmdline),    intent(inout) :: cline
        class(string),     intent(in)    :: truth_mean
        character(len=*),  intent(in)    :: basis_backend, state_backend
        integer,           intent(in)    :: nparts, nthr
        integer, optional, intent(in)    :: nbins, box_crop
        call cline%set('prg',               'flex_pca')
        call cline%set('projfile',          'phantom.simple')
        call cline%set('mkdir',             'no')
        call cline%set('oritype',           'ptcl3D')
        call cline%set('vol1',              truth_mean)
        call cline%set('pgrp',              'c1')
        call cline%set('mskdiam',           MSKDIAM)
        call cline%set('box_crop',          BOX)
        if( present(box_crop) ) call cline%set('box_crop', box_crop)
        call cline%set('lp',                8.0)
        call cline%set('nstates',           1)
        call cline%set('npreimages',        3)
        call cline%set('preimage_auto',     'yes')
        call cline%set('min_neff',          MIN_NEFF)
        call cline%set('neigs',             FLEX_PHANTOM_NCOMP)
        call cline%set('nkern',             1)
        call cline%set('state_axis',        1)
        call cline%set('n_probe_iters',     2)
        call cline%set('column_separation', 2)
        call cline%set('nbins',             1)
        if( present(nbins) ) call cline%set('nbins', nbins)
        call cline%set('rec_states',        'yes')
        call cline%set('rec_backend',       basis_backend)
        call cline%set('rec_states_backend', state_backend)
        call cline%set('maxits_pcg',        4)
        call cline%set('rtol',              0.0)
        call cline%set('objfun',            'euclid')
        call cline%set('umap',              'no')
        call cline%set('qsys_name',         'local')
        call cline%set('nparts',            nparts)
        call cline%set('nthr',              nthr)
    end subroutine configure_flex_pca_phantom

end module simple_test_flex_phantom_fixture
