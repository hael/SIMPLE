module simple_user_interface
use simple_defs
implicit none

public :: simple_program, make_user_interface, get_prg_ptr, list_distr_prgs_in_ui, list_shmem_prgs_in_ui
private

logical, parameter :: DEBUG = .false.

type simple_input_param
    character(len=:), allocatable :: key
    character(len=:), allocatable :: keytype ! (binary|multi|num|str|file)
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: descr_placeholder
    character(len=:), allocatable :: cval_default
    real                          :: rval_default = 0.
    logical :: required = .true.
end type simple_input_param

type :: simple_program
    private
    character(len=:), allocatable :: name
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: executable
    ! image input/output
    type(simple_input_param), allocatable :: img_ios(:)
    ! parameter input/output
    type(simple_input_param), allocatable :: parm_ios(:)
    ! alternative inputs
    type(simple_input_param), allocatable :: alt_ios(:)
    ! search controls
    type(simple_input_param), allocatable :: srch_ctrls(:)
    ! filter controls
    type(simple_input_param), allocatable :: filt_ctrls(:)
    ! mask controls
    type(simple_input_param), allocatable :: mask_ctrls(:)
    ! computer controls
    type(simple_input_param), allocatable :: comp_ctrls(:)
    ! simple project dependence flag
    logical :: sp_required = .false.
    ! existence flag
    logical :: exists = .false.
  contains
    procedure, private :: new
    procedure, private :: set_input_1
    procedure, private :: set_input_2
    procedure, private :: set_input_3
    generic,   private :: set_input => set_input_1, set_input_2, set_input_3
    procedure          :: print_ui
    procedure          :: print_cmdline
    procedure          :: print_prg_descr_long
    procedure          :: write2json
    procedure          :: get_required_keys
    procedure          :: is_distr
    procedure          :: requires_sp_project
    procedure, private :: kill
end type simple_program

! declare protected program specifications here
type(simple_program), target :: cluster2D
type(simple_program), target :: cluster2D_stream
type(simple_program), target :: make_cavgs
type(simple_program), target :: refine3D
type(simple_program), target :: initial_3Dmodel
type(simple_program), target :: cluster3D
type(simple_program), target :: preprocess
type(simple_program), target :: extract
type(simple_program), target :: motion_correct
type(simple_program), target :: ctf_estimate
type(simple_program), target :: pick
type(simple_program), target :: postprocess
type(simple_program), target :: mask
type(simple_program), target :: scale
type(simple_program), target :: volops

! declare common params here, with name same as flag
type(simple_input_param) :: astigtol
type(simple_input_param) :: bfac
type(simple_input_param) :: cs
type(simple_input_param) :: ctf
type(simple_input_param) :: fraca
type(simple_input_param) :: kv
type(simple_input_param) :: deftab
type(simple_input_param) :: dfmin
type(simple_input_param) :: dfmax
type(simple_input_param) :: dfstep
type(simple_input_param) :: frac
type(simple_input_param) :: hp
type(simple_input_param) :: lp
type(simple_input_param) :: lplim_crit
type(simple_input_param) :: inner
type(simple_input_param) :: maxits
type(simple_input_param) :: mirr
type(simple_input_param) :: msk
type(simple_input_param) :: mskfile
type(simple_input_param) :: mw
type(simple_input_param) :: nspace
type(simple_input_param) :: nparts
type(simple_input_param) :: ncls
type(simple_input_param) :: nthr
type(simple_input_param) :: objfun
type(simple_input_param) :: oritab
type(simple_input_param) :: outer
type(simple_input_param) :: outfile
type(simple_input_param) :: outvol
type(simple_input_param) :: phaseplate
type(simple_input_param) :: pgrp
type(simple_input_param) :: projfile
type(simple_input_param) :: pspecsz
type(simple_input_param) :: remap_cls
type(simple_input_param) :: smpd
type(simple_input_param) :: startit
type(simple_input_param) :: stk
type(simple_input_param) :: stktab
type(simple_input_param) :: trs
type(simple_input_param) :: update_frac
type(simple_input_param) :: weights2D

interface set_param
    module procedure set_param_1
    module procedure set_param_2
end interface set_param

contains

    ! public class methods

    subroutine make_user_interface
        call set_common_params
        call new_cluster2D
        call new_cluster2D_stream
        call new_refine3D
        call new_initial_3Dmodel
        call new_cluster3D
        call new_extract
        call new_make_cavgs
        call new_motion_correct
        call new_preprocess
        call new_pick
        call new_ctf_estimate
        call new_postprocess
        call new_mask
        call new_scale
        call new_volops
        ! ...
    end subroutine make_user_interface

    subroutine get_prg_ptr( which_program, ptr2prg )
        character(len=*), intent(in)   :: which_program
        class(simple_program), pointer :: ptr2prg
        select case(trim(which_program))
            case('cluster2D')
                ptr2prg => cluster2D
            case('cluster2D_stream')
                ptr2prg => cluster2D_stream
            case('make_cavgs')
                ptr2prg => make_cavgs
            case('refine3D')
                ptr2prg => refine3D
            case('initial_3Dmodel')
                ptr2prg => initial_3Dmodel
            case('cluster3D')
                ptr2prg => cluster3D
            case('preprocess')
                ptr2prg => preprocess
            case('extract')
                ptr2prg => extract
            case('motion_correct')
                ptr2prg => motion_correct
            case('ctf_estimate')
                ptr2prg => ctf_estimate
            case('pick')
                ptr2prg => pick
            case('postprocess')
                ptr2prg => postprocess
            case('mask')
                ptr2prg => mask
            case('scale')
                ptr2prg => scale
            case('volops')
                ptr2prg => volops
            case DEFAULT
                ptr2prg => null()
        end select
    end subroutine get_prg_ptr

    subroutine list_distr_prgs_in_ui
        write(*,'(A)') cluster2D%name
        write(*,'(A)') cluster2D_stream%name
        write(*,'(A)') ctf_estimate%name
        write(*,'(A)') initial_3Dmodel%name
        write(*,'(A)') make_cavgs%name
        write(*,'(A)') motion_correct%name
        write(*,'(A)') pick%name
        write(*,'(A)') preprocess%name
        write(*,'(A)') refine3D%name
    end subroutine list_distr_prgs_in_ui

    subroutine list_shmem_prgs_in_ui
        write(*,'(A)') extract%name
        write(*,'(A)') postprocess%name
    end subroutine list_shmem_prgs_in_ui

    ! private class methods

    subroutine set_common_params
        call set_param(projfile,      'projfile',      'file',   'Project file', 'SIMPLE projectfile', 'e.g. myproject.simple', .false., 'myproject.simple')
        call set_param(stk,           'stk',           'file',   'Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., 'stk.mrc')
        call set_param(stktab,        'stktab',        'file',   'List of per-micrograph particle stacks', 'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false., 'stktab.txt')
        call set_param(ctf,           'ctf',           'multi',  'CTF correction', 'Contrast Transfer Function correction; flip indicates that images have been phase-flipped prior(yes|no|flip){no}',&
        &'(yes|no|flip){no}', .true., 'no')
        call set_param(smpd,          'smpd',          'num',    'Sampling distance', 'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms', .true., 1.0)
        call set_param(phaseplate,    'phaseplate',    'binary', 'Phase-plate images', 'Images obtained with Volta phase-plate(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(deftab,        'deftab',        'file',   'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format with dfx, dfy and angast values',&
        &'.simple|.txt parameter file', .false., 'deftab'//trim(METADATA_EXT))
        call set_param(oritab,        'oritab',        'file',   'Orientation and CTF parameter file', 'Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false., 'oritab'//trim(METADATA_EXT))
        call set_param(outfile,       'outfile',       'file',   'Output orientation and CTF parameter file', 'Output Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false., 'outfile'//trim(METADATA_EXT))
        call set_param(startit,       'startit',       'num',    'First iteration', 'Index of first iteration when starting from a previous solution', 'start iterations from here', .false., 1.0)
        call set_param(trs,           'trs',           'num',    'Maximum translational shift', 'Maximum half-width for bund-constrained search of rotational origin shifts',&
        &'max shift per iteration in pixels{5}', .false., 0.0)
        call set_param(maxits,        'maxits',        'num',    'Max iterations', 'Maximum number of iterations', 'Max # iterations', .false., 100.)
        call set_param(hp,            'hp',            'num',    'High-pass limit', 'High-pass resolution limit', 'high-pass limit in Angstroms', .false., 100.)
        call set_param(lp,            'hp',            'num',    'Low-pass limit', 'Low-pass resolution limit', 'high-pass limit in Angstroms', .false., 20.)
        call set_param(msk,           'msk',           'num',    'Mask radius', 'Mask radius in pixels for application of a soft-edged circular mask to remove background noise', 'mask radius in pixels', .true., 0.)
        call set_param(inner,         'inner',         'num',    'Inner mask radius', 'Inner mask radius for omitting unordered cores of particles with high radial symmetry, typically icosahedral viruses',&
        &'inner mask radius in pixels', .false., 0.)
        call set_param(ncls,          'ncls',          'num', 'Number of 2D clusters', 'Number of groups to sort the particles &
        &into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .false., 200.)
        call set_param(nparts,        'nparts',        'num',    'Number of parts', 'Number of partitions for distrbuted memory execution. One part typically corresponds to one CPU socket in the distributed &
        &system. On a single-socket machine there may be speed benfits to dividing the jobs into a few (2-4) partitions, depending on memory capacity', 'divide job into # parts', .true., 1.0)
        call set_param(nthr,          'nthr',          'num',    'Number of threads per part', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .false., 1.0)
        call set_param(outer,         'outer',         'num',    'Outer mask radius', 'Outer mask radius for omitting unordered cores of particles with high radial symmetry, typically icosahedral viruses',&
        &'outer mask radius in pixels', .false., 0.)
        call set_param(update_frac,   'update_frac',   'num',    'Fractional update per iteration', 'Fraction of particles to update per iteration in incremental learning scheme for accelerated convergence &
        &rate(0.1-0.5){1.}', 'update this fraction per iter(0.1-0.5){1.0}', .false., 1.0)
        call set_param(frac,          'frac',          'num',    'Fraction of particles to include', 'Fraction of particles to include based on spectral score (median of FRC between reference and particle)',&
        'fraction of particles used(0.1-0.9){1.0}', .false., 1.0)
        call set_param(mskfile,       'mskfile',       'file',   'Input mask file', 'Input mask file to apply to reference volume(s) before projection', 'e.g. automask.mrc from postprocess', .false., 'mskfile.mrc')
        call set_param(pgrp,          'pgrp',          'str',    'Point-group symmetry', 'Point-group symmetry of particle(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .true., 'c1')
        call set_param(nspace,        'nspace',        'num',    'Number of projection directions', 'Number of projection directions &
        &used in discrete 3D orientation search', '# projections used{2500}', .false., 2500.)
        call set_param(objfun,        'objfun',        'binary', 'Objective function', 'Objective function(cc|ccres){cc}', '(cc|ccres){cc}', .false., 'cc')
        call set_param(weights2D,     'weights2D',     'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(remap_cls,     'remap_cls',     'binary', 'Whether to remap 2D clusters', 'Whether to remap the number of 2D clusters(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(kv,            'kv',            'num',    'Acceleration voltage', 'Acceleration voltage in kV', 'in kV', .false., 300.)
        call set_param(lplim_crit,    'lplim_crit',    'num', 'Low-pass limit FSC criterion', 'FSC criterion for determining the low-pass limit(0.143-0.5){0.3}',&
        &'low-pass FSC criterion(0.143-0.5){0.3}', .false., 0.3)
        call set_param(cs,            'cs',            'num', 'Spherical aberration', 'Spherical aberration constant(in mm){2.7}', 'in nm{2.7}', .false., 2.7)
        call set_param(fraca,         'fraca',         'num', 'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.1}', '{0.1}', .false., 0.1)
        call set_param(pspecsz,       'pspecsz',       'num', 'Size of power spectrum', 'Size of power spectrum in pixels', 'in pixels', .false., 512.)
        call set_param(dfmin,         'dfmin',         'num', 'Expected minimum defocus', 'Expected minimum defocus in microns{0.5}', 'in microns{0.5}', .false., 0.5)
        call set_param(dfmax,         'dfmax',         'num', 'Expected maximum defocus', 'Expected minimum defocus in microns{5.0}', 'in microns{5.0}', .false., 5.0)
        call set_param(dfstep,        'dfstep',        'num', 'Defocus step size', 'Defocus step size for grid search in microns{0.05}', 'in microns{0.05}', .false., 0.05)
        call set_param(astigtol,      'astigtol',      'num', 'Expected astigmatism', 'expected (tolerated) astigmatism(in microns){0.05}', 'in microns',  .false., 0.05)
        call set_param(mw,            'mw',            'num', 'Molecular weight','Molecular weight in kDa', 'in kDa', .false., 0.)
        call set_param(mirr,          'mirr',          'multi', 'Perform mirroring', 'Whether to mirror and along which axis(no|x|y){no}', '(no|x|y){no}', .false., 'no')
        call set_param(bfac,          'bfac',          'num', 'B-factor for sharpening','B-factor for sharpening in Angstroms^2', 'B-factor in Angstroms^2', .false., 200.)
        call set_param(outvol,        'outvol',        'file', 'Output volume name', 'Output volume name', 'e.g. outvol.mrc', .false., '')
    end subroutine set_common_params

    subroutine set_param_1( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        type(simple_input_param), intent(inout) :: self
        character(len=*),         intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                  intent(in)    :: required
        real,                     intent(in)    :: default_value
        allocate(self%key,               source=trim(key))
        allocate(self%keytype,           source=trim(keytype))
        allocate(self%descr_short,       source=trim(descr_short))
        allocate(self%descr_long,        source=trim(descr_long))
        allocate(self%descr_placeholder, source=trim(descr_placeholder))
        self%required = required
        if( .not. self%required ) self%rval_default = default_value
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        type(simple_input_param), intent(inout) :: self
        character(len=*),         intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                  intent(in)    :: required
        character(len=*),         intent(in)    :: default_value
        allocate(self%key,               source=trim(key))
        allocate(self%keytype,           source=trim(keytype))
        allocate(self%descr_short,       source=trim(descr_short))
        allocate(self%descr_long,        source=trim(descr_long))
        allocate(self%descr_placeholder, source=trim(descr_placeholder))
        self%required = required
        if( .not. self%required ) allocate(self%cval_default, source=trim(default_value))
    end subroutine set_param_2

    subroutine new_refine3D
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',& ! name
        &'3D volume refinement',&                                                                          ! descr_short
        &'is a distributed workflow for 3D volume refinement based on probabilistic projection matching',& ! descr_long
        &'simple_distr_exec',&                                                                             ! executable
        &1, 0, 0, 12, 7, 5, 2, .true.)                                                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%set_input('img_ios', 1, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. recvol.mrc', .false., 'vol1.mrc')
        ! parameter input/output
        !<empty>
        ! alternative inputs
        !<empty>
        ! search controls
        call refine3D%set_input('srch_ctrls', 1, nspace)
        call refine3D%set_input('srch_ctrls', 2, startit)
        call refine3D%set_input('srch_ctrls', 3, trs)
        call refine3D%set_input('srch_ctrls', 4, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call refine3D%set_input('srch_ctrls', 5, maxits)
        call refine3D%set_input('srch_ctrls', 6, update_frac)
        call refine3D%set_input('srch_ctrls', 7, frac)
        call refine3D%set_input('srch_ctrls', 8, pgrp)
        call refine3D%set_input('srch_ctrls', 9, 'nnn', 'num', 'Number of nearest neighbours', 'Number of nearest projection direction &
        &neighbours in neigh=yes refinement', '# projection neighbours{10% of search space}', .false., 200.)
        call refine3D%set_input('srch_ctrls', 10, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to reconstruct',&
        '# states to reconstruct', .false., 1.0)
        call refine3D%set_input('srch_ctrls', 11, objfun)
        call refine3D%set_input('srch_ctrls', 12, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc|single|multi|greedy_single|greedy_multi|cluster|&
        &clustersym){no}', '(snhc|single|multi|greedy_single|greedy_multi|cluster|clustersym){single}', .false., 'single')
        ! filter controls
        call refine3D%set_input('filt_ctrls', 1, hp)
        call refine3D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call refine3D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call refine3D%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call refine3D%set_input('filt_ctrls', 5, lplim_crit)
        call refine3D%set_input('filt_ctrls', 6, 'eo', 'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){yes}', '(yes|no){yes}', .false., 'no')
        call refine3D%set_input('filt_ctrls', 7, 'weights3D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! mask controls
        call refine3D%set_input('mask_ctrls', 1, msk)
        call refine3D%set_input('mask_ctrls', 2, inner)
        call refine3D%set_input('mask_ctrls', 3, mskfile)
        call refine3D%set_input('mask_ctrls', 4, 'focusmsk', 'num', 'Mask radius in focused refinement', 'Mask radius in pixels for application of a soft-edged circular &
        &mask to remove background noise in focused refinement', 'focused mask radius in pixels', .false., 0.)
        call refine3D%set_input('mask_ctrls', 5, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge{10}', .false., 10.)
        ! computer controls
        call refine3D%set_input('comp_ctrls', 1, nparts)
        call refine3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_refine3D

    subroutine new_cluster3D
        ! PROGRAM SPECIFICATION
        call cluster3D%new(&
        &'cluster3D',& ! name
        &'3D heterogeneity analysis',&                                             ! descr_short
        &'is a distributed workflow for heterogeneity analysis by 3D clustering',& ! descr_long
        &'simple_distr_exec',&                                                     ! executable
        &0, 1, 0, 7, 6, 5, 2, .true.)                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster3D%set_input('parm_ios', 1, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to separate',&
        '# states to separate', .true., 2.0)
        ! alternative inputs
        !<empty>
        ! search controls
        call cluster3D%set_input('srch_ctrls', 1, nspace)
        call cluster3D%set_input('srch_ctrls', 2, startit)
        call cluster3D%set_input('srch_ctrls', 3, maxits)
        call cluster3D%set_input('srch_ctrls', 4, frac)
        call cluster3D%set_input('srch_ctrls', 5, pgrp)
        call cluster3D%set_input('srch_ctrls', 6, objfun)
        call cluster3D%set_input('srch_ctrls', 7, 'refine', 'binary', 'Refinement mode', 'Refinement mode(cluster|clustersym)&
        &){cluster}', '(cluster|clustersym){cluster}', .false., 'cluster')
        ! filter controls
        call cluster3D%set_input('filt_ctrls', 1, hp)
        call cluster3D%set_input('filt_ctrls', 2, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call cluster3D%set_input('filt_ctrls', 3, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call cluster3D%set_input('filt_ctrls', 4, lplim_crit)
        call cluster3D%set_input('filt_ctrls', 5, 'eo', 'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){yes}', '(yes|no){yes}', .false., 'no')
        call cluster3D%set_input('filt_ctrls', 6, 'weights3D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! mask controls
        call cluster3D%set_input('mask_ctrls', 1, msk)
        call cluster3D%set_input('mask_ctrls', 2, inner)
        call cluster3D%set_input('mask_ctrls', 3, mskfile)
        call cluster3D%set_input('mask_ctrls', 4, 'focusmsk', 'num', 'Mask radius in focused refinement', 'Mask radius in pixels for application of a soft-edged circular &
        &mask to remove background noise in focused refinement', 'focused mask radius in pixels', .false., 0.)
        call cluster3D%set_input('mask_ctrls', 5, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge{10}', .false., 10.)
        ! computer controls
        call cluster3D%set_input('comp_ctrls', 1, nparts)
        call cluster3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster3D

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',&                 ! descr_long
        &'simple_distr_exec',&                                                  ! executable
        &1, 0, 1, 10, 7, 2, 2, .true.)                                          ! # entries in each group

        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster2D%set_input('img_ios', 1, 'refs', 'file', 'Initial references',&
        &'Initial 2D references used to bootstrap the search', 'xxx.mrc file with references', .false., 'refs.mrc')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        call cluster2D%set_input('alt_ios', 1, projfile)
        ! search controls
        call cluster2D%set_input('srch_ctrls', 1, ncls)
        cluster2D%srch_ctrls(1)%required = .true.
        call cluster2D%set_input('srch_ctrls', 2, startit)
        call cluster2D%set_input('srch_ctrls', 3, trs)
        call cluster2D%set_input('srch_ctrls', 4, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 'yes')
        call cluster2D%set_input('srch_ctrls', 5, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D%set_input('srch_ctrls', 6, 'dyncls', 'binary', 'Dynamic reallocation of clusters', 'Dynamic reallocation of clusters &
        &that fall below a minimum population by randomization(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D%set_input('srch_ctrls', 7, maxits)
        call cluster2D%set_input('srch_ctrls', 8, update_frac)
        call cluster2D%set_input('srch_ctrls', 9, frac)
        call cluster2D%set_input('srch_ctrls',10, 'bfac', 'num', 'Correlation B-factor','B-factor for the objective function Angstroms^2', 'B-factor in Angstroms^2 (>0.0)', .false., 1000.)
        ! filter controls
        call cluster2D%set_input('filt_ctrls', 1, hp)
        call cluster2D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false., 20.)
        call cluster2D%set_input('filt_ctrls', 4, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &few iterations of search, before the automatic scheme kicks in. Also controls the degree of downsampling in the first &
        phase', 'initial low-pass limit in Angstroms', .false., 15.)
        call cluster2D%set_input('filt_ctrls', 5, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.)
        call cluster2D%set_input('filt_ctrls', 6, 'match_filt', 'binary', 'Matched filter', 'Filter to maximize the signal-to-noise &
        &ratio (SNR) in the presence of additive stochastic noise. Sometimes causes over-fitting and needs to be turned off(yes|no){yes}',&
        '(yes|no){yes}', .false., 'yes')
        call cluster2D%set_input('filt_ctrls', 7, weights2D)
        ! mask controls
        call cluster2D%set_input('mask_ctrls', 1, msk)
        call cluster2D%set_input('mask_ctrls', 2, inner)
        ! computer controls
        call cluster2D%set_input('comp_ctrls', 1, nparts)
        call cluster2D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster2D

    subroutine new_cluster2D_stream
        ! PROGRAM SPECIFICATION
        call cluster2D_stream%new(&
        &'cluster2D_stream',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images in streaming mode',&                         ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm in streaming mode',&  ! descr_long
        &'simple_distr_exec',&                                                  ! executable
        &1, 0, 0, 6, 5, 2, 2, .true.)                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster2D_stream%set_input('img_ios', 1, 'dir_ptcls', 'file', 'Particles directory',&
        &'Directory where particles and CTF parameters are automatically detected', 'e.g. extract/', .true., '')
        ! parameter input/output
        !<empty>
        ! alternative inputs
        !<empty>
        ! search controls
        call cluster2D_stream%set_input('srch_ctrls', 1, 'ncls_start', 'num', 'Starting number of clusters',&
        &'Minimum number of class averagages to initiate 2D clustering', 'initial # clusters', .true., 50.)
        cluster2D_stream%srch_ctrls(1)%required = .true.
        call cluster2D_stream%set_input('srch_ctrls', 2, 'nptcls_per_cls', 'num', 'Particles per cluster',&
        &'Number of incoming particles for which one new class average is generated', '# particles per cluster', .true., 200.)
        cluster2D_stream%srch_ctrls(2)%required = .true.
        call cluster2D_stream%set_input('srch_ctrls', 3, trs)
        call cluster2D_stream%set_input('srch_ctrls', 4, objfun)
        call cluster2D_stream%set_input('srch_ctrls', 5, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 'yes')
        call cluster2D_stream%set_input('srch_ctrls', 6, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call cluster2D_stream%set_input('filt_ctrls', 1, hp)
        call cluster2D_stream%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D_stream%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false., 20.)
        call cluster2D_stream%set_input('filt_ctrls', 4, 'match_filt', 'binary', 'Matched filter', 'Filter to maximize the signal-to-noise &
        &ratio (SNR) in the presence of additive stochastic noise. Sometimes causes over-fitting and needs to be turned off(yes|no){yes}',&
        '(yes|no){yes}', .false., 'yes')
        call cluster2D_stream%set_input('filt_ctrls', 5, weights2D)
        ! mask controls
        call cluster2D_stream%set_input('mask_ctrls', 1, msk)
        call cluster2D_stream%set_input('mask_ctrls', 2, inner)
        ! computer controls
        call cluster2D_stream%set_input('comp_ctrls', 1, nparts)
        call cluster2D_stream%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster2D_stream

    subroutine new_initial_3Dmodel
        ! PROGRAM SPECIFICATION
        call initial_3Dmodel%new(&
        &'initial_3Dmodel',& ! name
        &'3D abinitio model generation from class averages',&                           ! descr_short
        &'is a distributed workflow for generating an initial 3D model from class'&
        &' averages obtained with cluster2D',&                                          ! descr_long
        &'simple_distr_exec',&                                                          ! executable
        &1, 1, 0, 9, 3, 3, 2, .true.)                                                   ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call initial_3Dmodel%set_input('img_ios', 1, 'stk', 'file', 'Class averages image stack', 'Class averages image stack', 'xxx.mrc file with class averages', .true., '')
        ! parameter input/output
        !<empty>
        ! alternative inputs
        !<empty>
        ! search controls
        call initial_3Dmodel%set_input('srch_ctrls', 1, nspace)
        call initial_3Dmodel%set_input('srch_ctrls', 2, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call initial_3Dmodel%set_input('srch_ctrls', 3, maxits)
        call initial_3Dmodel%set_input('srch_ctrls', 4, update_frac)
        call initial_3Dmodel%set_input('srch_ctrls', 5, frac)
        call initial_3Dmodel%set_input('srch_ctrls', 6, pgrp)
        call initial_3Dmodel%set_input('srch_ctrls', 7, 'pgrp_known', 'binary', 'Point-group applied directly', 'Point-group applied direclty rather than first doing a reconstruction &
        &in c1 and searching for the symmerty axis(yes|no){no}', '(yes|no){no}', .false., 'no')
        call initial_3Dmodel%set_input('srch_ctrls', 8, objfun)
        call initial_3Dmodel%set_input('srch_ctrls', 9, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Final low-pass limit controls the degree of down-scaling(yes|no){yes}','(yes|no){yes}', .false., 'yes')
        ! filter controls
        call initial_3Dmodel%set_input('filt_ctrls', 1, hp)
        call initial_3Dmodel%set_input('filt_ctrls', 2, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms', .false., 0.)
        call initial_3Dmodel%set_input('filt_ctrls', 3, 'lpstop',  'num', 'Final low-pass limit',   'Final low-pass limit',   'low-pass limit in Angstroms', .false., 8.)
        ! mask controls
        call initial_3Dmodel%set_input('mask_ctrls', 1, msk)
        call initial_3Dmodel%set_input('mask_ctrls', 2, inner)
        call initial_3Dmodel%set_input('mask_ctrls', 3, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false., 10.)
        ! computer controls
        call initial_3Dmodel%set_input('comp_ctrls', 1, nparts)
        call initial_3Dmodel%set_input('comp_ctrls', 2, nthr)
    end subroutine new_initial_3Dmodel

    subroutine new_preprocess
        ! PROGRAM SPECIFICATION
        call preprocess%new(&
        &'preprocess', & ! name
        &'is a program that performs for movie alignment',&                                 ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in sequence or streaming mode as the microscope collects the data',&
        &'simple_distr_exec',&                                                              ! executable
        &1, 10, 0, 13, 5, 0, 2, .true.)                                                     ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess%set_input('img_ios', 1, 'dir', 'file', 'Output directory', 'Output directory', 'e.g. preprocess/', .false., 'preprocess')
        ! parameter input/output
        call preprocess%set_input('parm_ios', 1, 'dose_rate', 'num', 'Dose rate', 'Dose rate in e/Ang^2/sec', 'in e/Ang^2/sec', .false., 6.0)
        call preprocess%set_input('parm_ios', 2, 'exp_time', 'num', 'Exposure time', 'Exposure time in seconds', 'in seconds', .false., 10.)
        call preprocess%set_input('parm_ios', 3, 'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.0)
        call preprocess%set_input('parm_ios', 4, 'pcontrast', 'binary', 'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false., 'black')
        call preprocess%set_input('parm_ios', 5, 'stream', 'binary', 'Streaming on/off', 'Whether to activate streaming mode(yes|no){yes}', '(yes|no){no}', .false., 'no')
        call preprocess%set_input('parm_ios', 6, 'box_extract', 'num', 'Box size on extraction', 'Box size on extraction in pixels', 'in pixels', .false., 0.)
        call preprocess%set_input('parm_ios', 7, 'refs', 'file', 'Picking 2D references',&
        &'2D references used for automated picking', 'e.g. pickrefs.mrc file with references', .false., 'pickrefs.mrc')
        call preprocess%set_input('parm_ios', 8, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_')
        call preprocess%set_input('parm_ios', 9, pspecsz)
        call preprocess%set_input('parm_ios',10, 'numlen', 'num', 'Length of number string', 'Length of number string', '...', .false., 5.0)
        ! alternative inputs
        !<empty>
        ! search controls
        call preprocess%set_input('srch_ctrls', 1, trs)
        call preprocess%set_input('srch_ctrls', 2, 'startit', 'num', 'Initial movie alignment iteration', 'Initial movie alignment iteration', '...', .false., 1.0)
        call preprocess%set_input('srch_ctrls', 3, 'nframesgrp', 'num', '# frames to group', '# frames to group before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call preprocess%set_input('srch_ctrls', 4, 'fromf', 'num', 'First frame index', 'First frame to include in the alignment', '...', .false., 1.0)
        call preprocess%set_input('srch_ctrls', 5, 'tof', 'num', 'Last frame index', 'Last frame to include in the alignment', '...', .false., 1.0)
        call preprocess%set_input('srch_ctrls', 6, 'nsig', 'num', '# of sigmas in motion_correct', '# of standard deviation threshold for outlier removal in motion_correct{6}', '{6}', .false., 6.)
        call preprocess%set_input('srch_ctrls', 7, dfmin)
        call preprocess%set_input('srch_ctrls', 8, dfmax)
        call preprocess%set_input('srch_ctrls', 9, dfstep)
        call preprocess%set_input('srch_ctrls',10, astigtol)
        call preprocess%set_input('srch_ctrls',11, 'thres', 'num', 'Picking distance threshold','Picking distance filer (in pixels)', 'in pixels', .false., 0.)
        call preprocess%set_input('srch_ctrls',12, 'rm_outliers', 'binary', 'Remove micrograph image outliers for picking',&
        & 'Remove micrograph image outliers for picking(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call preprocess%set_input('srch_ctrls',13, 'ndev', 'num', '# of sigmas for picking clustering', '# of standard deviations threshold for picking one cluster clustering{2}', '{2}', .false., 2.)
        ! filter controls
        call preprocess%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){15}', 'in Angstroms{15}', .false., 15.)
        call preprocess%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess%set_input('filt_ctrls', 3, 'lp_ctf_estimate', 'num', 'Low-pass limit for CTF parameter estimation',&
        & 'Low-pass limit for CTF parameter estimation in Angstroms{5}', 'in Angstroms{5}', .false., 5.)
        call preprocess%set_input('filt_ctrls', 4, 'hp_ctfestimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30.)
        call preprocess%set_input('filt_ctrls', 5, 'lp_pick', 'num', 'Low-pass limit for picking',&
        & 'Low-pass limit for picking in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        ! mask controls
        ! <empty>
        ! computer controls
        call preprocess%set_input('comp_ctrls', 1, nparts)
        call preprocess%set_input('comp_ctrls', 2, nthr)
    end subroutine new_preprocess

    subroutine new_ctf_estimate
        ! PROGRAM SPECIFICATION
        call ctf_estimate%new(&
        &'ctf_estimate', & ! name
        &'is a distributed SIMPLE workflow for CTF parameter fitting',&        ! descr_short
        &'is a distributed SIMPLE workflow for CTF parameter fitting',&        ! descr_long
        &'simple_distr_exec',&                                                 ! executable
        &0, 2, 0, 4, 2, 0, 2, .true.)                                          ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_estimate%set_input('parm_ios', 1, 'dir', 'file', 'Output directory', 'Output directory', 'e.g. ctf_estimate/', .false., 'ctf_estimate')
        call ctf_estimate%set_input('parm_ios', 2, pspecsz)
        ! alternative inputs
        ! <empty>
        ! search controls
        call ctf_estimate%set_input('srch_ctrls', 1, dfmin)
        call ctf_estimate%set_input('srch_ctrls', 2, dfmax)
        call ctf_estimate%set_input('srch_ctrls', 3, dfstep)
        call ctf_estimate%set_input('srch_ctrls', 4, astigtol)
        ! filter controls
        call ctf_estimate%set_input('filt_ctrls', 1, lp)
        call ctf_estimate%set_input('filt_ctrls', 2, hp)
        ! mask controls
        ! <empty>
        ! computer controls
        call ctf_estimate%set_input('comp_ctrls', 1, nparts)
        call ctf_estimate%set_input('comp_ctrls', 2, nthr)
    end subroutine new_ctf_estimate

    subroutine new_pick
        ! PROGRAM SPECIFICATION
        call pick%new(&
        &'pick', & ! name
        &'is a distributed workflow for template-based particle picking',&      ! descr_short
        &'is a distributed workflow for template-based particle picking',&      ! descr_long
        &'simple_distr_exec',&                                                  ! executable
        &0, 2, 0, 2, 1, 0, 1, .true.)                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call pick%set_input('parm_ios', 1, 'refs', 'file', 'Picking 2D references',&
        &'2D references used for automated picking', 'e.g. pickrefs.mrc file with references', .false., '')
        call pick%set_input('parm_ios', 2, 'dir', 'file', 'Output directory', 'Output directory', 'e.g. pick/', .false., 'pick')
        ! alternative inputs
        ! <empty>
        ! search controls
        call pick%set_input('srch_ctrls',1, 'thres', 'num', 'Distance threshold','Distance filer (in pixels)', 'in pixels', .false., 0.)
        call pick%set_input('srch_ctrls',2, 'ndev', 'num', '# of sigmas for clustering', '# of standard deviations threshold for one cluster clustering{2}', '{2}', .false., 2.)
        ! filter controls
        call pick%set_input('filt_ctrls', 1, 'lp', 'num', 'Low-pass limit','Low-pass limit in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        ! mask controls
        ! <empty>
        ! computer controls
        call pick%set_input('comp_ctrls', 1, nthr)
    end subroutine new_pick

    subroutine new_motion_correct
        ! PROGRAM SPECIFICATION
        call motion_correct%new(&
        &'motion_correct', & ! name
        &'is a program that performs for movie alignment',&                                     ! descr_short
        &'is a distributed workflow for movie alignment or or motion correction based on the same&
        & principal strategy as Grigorieffs program (hence the name). There are two important&
        & differences: automatic weighting of the frames using a correlation-based M-estimator and&
        & continuous optimisation of the shift parameters. Input is a textfile with absolute paths&
        & to movie files in addition to a few input parameters, some of which deserve a comment. If&
        & dose_rate and exp_time are given the individual frames will be low-pass filtered accordingly&
        & (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to&
        & the down-scaling factor (for super-resolution movies). If nframesgrp is given the frames will&
        & be pre-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given, a&
        & contiguous subset of frames will be averaged without any dose-weighting applied.',&   ! descr_long
        &'simple_distr_exec',&                                                                  ! executable
        &0, 7, 0, 6, 2, 0, 2, .true.)                                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call motion_correct%set_input('parm_ios', 1, 'dir', 'file', 'Output directory', 'Output directory', 'e.g. motion_correct/', .false., 'motion_correct')
        call motion_correct%set_input('parm_ios', 2, 'dose_rate', 'num', 'Dose rate', 'Dose rate in e/Ang^2/sec', 'in e/Ang^2/sec', .false., 6.)
        call motion_correct%set_input('parm_ios', 3, 'exp_time', 'num', 'Exposure time', 'Exposure time in seconds', 'in seconds', .false., 10.)
        call motion_correct%set_input('parm_ios', 4, 'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.)
        call motion_correct%set_input('parm_ios', 5, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., '')
        call motion_correct%set_input('parm_ios', 6, pspecsz)
        call motion_correct%set_input('parm_ios', 7, 'numlen', 'num', 'Length of number string', 'Length of number string', '...', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call motion_correct%set_input('srch_ctrls', 1, trs)
        call motion_correct%set_input('srch_ctrls', 2, 'startit', 'num', 'Initial iteration', 'Initial iteration', '...', .false., 1.)
        call motion_correct%set_input('srch_ctrls', 3, 'nframesgrp', 'num', '# frames to group', '# frames to group before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call motion_correct%set_input('srch_ctrls', 4, 'fromf', 'num', 'First frame index', 'First frame to include in the alignment', '...', .false., 1.)
        call motion_correct%set_input('srch_ctrls', 5, 'tof', 'num', 'Last frame index', 'Last frame to include in the alignment', '...', .false., 1.)
        call motion_correct%set_input('srch_ctrls', 6, 'nsig', 'num', '# of sigmas', '# of standard deviation threshold for outlier removal', '{6}', .false., 6.)
        ! filter controls
        call motion_correct%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms)', 'in Angstroms', .false., 15.)
        call motion_correct%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms)', 'in Angstroms', .false., 8.)
        ! mask controls
        ! <empty>
        ! computer controls
        call motion_correct%set_input('comp_ctrls', 1, nparts)
        call motion_correct%set_input('comp_ctrls', 2, nthr)
    end subroutine new_motion_correct

    subroutine new_extract
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', & ! name
        &'is a program that extracts particle images from integrated movies',&                  ! descr_short
        &'Boxfiles are assumed to be in EMAN format but we provide a conversion script'//&      ! descr long
        &' (relion2emanbox.pl) for *.star files containing particle coordinates obtained&
        & with Relion. In addition to one single-particle image stack per micrograph the&
        & program produces a parameter files that should be concatenated for use in&
        & conjunction with other SIMPLE programs.',&
        &'simple_exec',&                                                                        ! executable
        &0, 4, 0, 0, 0, 0, 0, .true.)                                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract%set_input('parm_ios', 1, 'dir', 'string', 'Ouput directory',&
        &'Ouput directory for single-particle images & CTF parameters', 'e.g. extract/', .false., 'extract')
        call extract%set_input('parm_ios', 2, 'box', 'num', 'Box size', 'Square box size in pixels', 'in pixels', .false., 0.)
        call extract%set_input('parm_ios', 3, 'pcontrast', 'binary', 'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false., 'black')
        call extract%set_input('parm_ios', 4, 'outside', 'binary', 'Extract outside boundaries', 'Extract boxes outside the micrograph boundaries(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_extract

    subroutine new_make_cavgs
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &     ! name
        &'is used to produce class averages',&     ! descr_short
        &'is used to produce class averages or initial random references&
        &for cluster2D execution',&                ! descr_long
        &'simple_distr_exec',&                     ! executable
        &1, 4, 0, 0, 0, 0, 2, .true.)              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call make_cavgs%set_input('img_ios', 1, 'refs', 'file', 'Output 2D references',&
        &'Output 2D references', 'xxx.mrc file with references', .false., '')
        ! parameter input/output
        call make_cavgs%set_input('parm_ios', 1, ncls)
        call make_cavgs%set_input('parm_ios', 2, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call make_cavgs%set_input('parm_ios', 3, weights2D)
        call make_cavgs%set_input('parm_ios', 4, remap_cls)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call make_cavgs%set_input('comp_ctrls', 1, nparts)
        call make_cavgs%set_input('comp_ctrls', 2, nthr)
    end subroutine new_make_cavgs

    subroutine new_postprocess
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',& ! name
        &'is a program for post-processing of volumes',&                            ! descr_short
        &'Use program volops to estimate the B-factor with the Guinier plot',&      ! descr_long
        &'simple_exec',&                                                            ! executable
        &1, 1, 0, 0, 7, 9, 1, .false.)                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to post-process', &
        & 'input volume e.g. recvol.mrc', .true., '')
        ! parameter input/output
        call postprocess%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%set_input('filt_ctrls', 1, hp)
        call postprocess%set_input('filt_ctrls', 2, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 15.)
        call postprocess%set_input('filt_ctrls', 3, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false., 20.)
        call postprocess%set_input('filt_ctrls', 4, 'vol_filt', 'file', 'Input filter volume', 'Input filter volume',&
        & 'input filter volume e.g. aniso_optlp_state01.mrc', .false., '')
        call postprocess%set_input('filt_ctrls', 5, 'fsc', 'file', 'Binary file with FSC info', 'Binary file with FSC info&
        & for filtering', 'input binary file e.g. fsc_state01.bin', .false., 'fsc_state01.bin')
        call postprocess%set_input('filt_ctrls', 6, bfac)
        call postprocess%set_input('filt_ctrls', 7, mirr)
        ! mask controls
        call postprocess%set_input('mask_ctrls', 1, msk)
        call postprocess%set_input('mask_ctrls', 2, inner)
        call postprocess%set_input('mask_ctrls', 3, mskfile)
        call postprocess%set_input('mask_ctrls', 4, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false., 1.)
        call postprocess%set_input('mask_ctrls', 5, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for enevloppe mask generation', 'Volume threshold', .false., 0.)
        call postprocess%set_input('mask_ctrls', 6, 'automsk', 'binary', 'Perform envelope masking',&
        &'Whether to generate an envelope mask(yes|no){no}', '(yes|no){no}', .false., 'no')
        call postprocess%set_input('mask_ctrls', 7, mw)
        call postprocess%set_input('mask_ctrls', 8, 'width', 'num', 'Inner mask falloff',&
        &'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false., 10.)
        call postprocess%set_input('mask_ctrls', 9, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6.)
        ! computer controls
        call postprocess%set_input('comp_ctrls', 1, nthr)
    end subroutine new_postprocess

    subroutine new_mask
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',& ! name
        &'is a program for masking images and volumes',&                            ! descr_short
        &'If you want to mask your images with a spherical mask with a soft &
        & falloff, set msk to the radius in pixels',&                               ! descr_long
        &'simple_exec',&                                                            ! executable
        &0, 3, 2, 0, 1,11, 1, .false.)                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mask%set_input('parm_ios', 1, smpd)
        call mask%set_input('parm_ios', 2, oritab)
        call mask%set_input('parm_ios', 3, outfile)
        ! alternative inputs
        call mask%set_input('alt_ios', 1, stk)
        call mask%set_input('alt_ios', 2, 'vol1', 'file', 'Volume', 'Volume to mask', &
        & 'input volume e.g. recvol.mrc', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        call mask%set_input('filt_ctrls', 1, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 15.)
        ! mask controls
        call mask%set_input('mask_ctrls', 1, msk)
        mask%mask_ctrls(1)%required = .false.
        call mask%set_input('mask_ctrls', 2, inner)
        call mask%set_input('mask_ctrls', 3, outer)
        call mask%set_input('mask_ctrls', 4, mskfile)
        call mask%set_input('mask_ctrls', 5, 'msktype', 'binary', 'Mask type',&
        &'Type of mask to use(soft|hard){soft}', '(soft|hard){soft}', .false., 'soft')
        call mask%set_input('mask_ctrls', 6, 'automsk', 'multi', 'Perform 2D envelope masking',&
        &'Whether to perform 2D envelope masking(cavg|no){no}', '(cavg|no){no}', .false., 'no')
        call mask%set_input('mask_ctrls', 7, mw)
        call mask%set_input('mask_ctrls', 8, 'width', 'num', 'Inner mask falloff',&
        &'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false., 10.)
        call mask%set_input('mask_ctrls', 9, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6.)
        call mask%set_input('mask_ctrls',10, 'taper_edges', 'binary', 'Taper edges',&
        &'Whether to taper the edges of image/volume(yes|no){no}', '(yes|no){no}', .false., 'no')
        call mask%set_input('mask_ctrls',11, 'pdbfile', 'file', 'PDB for 3D envelope masking',&
        &'PDB file used to determine the mask', 'e.g. molecule.pdb', .false., '')
        ! computer controls
        call mask%set_input('comp_ctrls', 1, nthr)
    end subroutine new_mask

    subroutine new_scale
        ! PROGRAM SPECIFICATION
        call scale%new(&
        &'scale', & ! name
        &'is a program for re-scaling MRC & SPIDER stacks and volumes',&                        ! descr_short
        &'is a program for re-scaling, clipping and padding MRC & SPIDER stacks and volumes',&  ! descr_long
        &'simple_exec',&                                                                        ! executable
        &0, 8, 3, 0, 0, 0, 1, .false.)                                                          ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call scale%set_input('parm_ios', 1, smpd)
        call scale%set_input('parm_ios', 2, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'in pixels', .false., 0.)
        call scale%set_input('parm_ios', 3, 'scale', 'num', 'Scaling ratio', 'Target box ratio for scaling(0-1)', '(0-1)', .false., 1.)
        call scale%set_input('parm_ios', 4, 'scale2', 'num', 'Second scaling ratio', 'Second target box ratio for scaling(0-1)', '(0-1)', .false., 1.)
        call scale%set_input('parm_ios', 5, 'clip', 'num', 'Clipped box size', 'Target box size for clipping in pixels', 'in pixels', .false., 0.)
        call scale%set_input('parm_ios', 6, outvol)
        call scale%set_input('parm_ios', 7, 'outstk', 'file', 'Output stack name', 'Output images stack name', 'e.g. outstk.mrc', .false., '')
        call scale%set_input('parm_ios', 8, 'outstk2', 'file', 'Second output stack name', 'Second output images stack name', 'e.g. outstk2.mrc', .false., '')
        ! alternative inputs
        call scale%set_input('alt_ios', 1, stk)
        scale%alt_ios(1)%required = .false.
        call scale%set_input('alt_ios', 2, 'vol1', 'file', 'Input volume', 'Input volume to re-scale',&
        &'input volume e.g. recvol.mrc', .false., '')
        call scale%set_input('alt_ios', 3, 'filetab', 'file', 'Stacks list',&
        &'List of stacks of images to rescale', 'list input e.g. stktab.txt', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call scale%set_input('comp_ctrls', 1, nthr)
    end subroutine new_scale

    subroutine new_volops
        ! PROGRAM SPECIFICATION
        call volops%new(&
        &'volops',& ! name
        &'is a program for standard volume editing',&                                             ! descr_short
        &'provides standard single-particle image processing routines to MRC or SPIDER volumes',& ! descr_long
        &'simple_exec',&                                                                          ! executable
        &2, 13, 0, 0, 2, 0, 1, .false.)                                                            ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call volops%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to mask', &
        & 'input volume e.g. recvol.mrc', .false., '')
        call volops%set_input('img_ios', 2, outvol)
        ! ! parameter input/output
        call volops%set_input('parm_ios', 1, smpd)
        volops%parm_ios(1)%required = .false.
        call volops%set_input('parm_ios', 2, 'guinier', 'binary', 'Guinier plot','calculate Guinier plot(yes|no){no}', '(yes|no){no}', .false., 'no')
        call volops%set_input('parm_ios', 3, 'neg', 'binary', 'Invert contrast','Invert volume contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call volops%set_input('parm_ios', 4, 'snr', 'num', 'SNR','Adds noise to the volume', 'signal to noise ratio(0.)', .false., 0.)
        call volops%set_input('parm_ios', 5, mirr)
        call volops%set_input('parm_ios', 6, bfac)
        call volops%set_input('parm_ios', 7, 'e1', 'num', 'Rotation along alpha','Alpha Euler angle', 'in degrees', .false., 0.)
        call volops%set_input('parm_ios', 8, 'e2', 'num', 'Rotation along Theta','Theat Euler angle', 'in degrees', .false., 0.)
        call volops%set_input('parm_ios', 9, 'e3', 'num', 'Rotation along Psi','Psi Euler angle', 'in degrees', .false., 0.)
        call volops%set_input('parm_ios',10, 'xsh', 'num', 'Translation along x-axis','Shift along X in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',11, 'ysh', 'num', 'Translation along y-axis','Shift along Y in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',12, 'zsh', 'num', 'Translation along z-axis','Shift along Z in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',13, outfile)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call volops%set_input('filt_ctrls', 1, lp)
        call volops%set_input('filt_ctrls', 2, hp)
        ! mask controls
        ! <empty>
        ! computer controls
        call volops%set_input('comp_ctrls', 1, nthr)
    end subroutine new_volops

    ! instance methods

    subroutine new( self, name, descr_short, descr_long, executable, n_img_ios, n_parm_ios,&
        &n_alt_ios, n_srch_ctrls, n_filt_ctrls, n_mask_ctrls, n_comp_ctrls, sp_required )
        class(simple_program), intent(inout) :: self
        character(len=*),      intent(in)    :: name, descr_short, descr_long, executable
        integer,               intent(in)    :: n_img_ios, n_parm_ios, n_alt_ios, n_srch_ctrls
        integer,               intent(in)    :: n_filt_ctrls, n_mask_ctrls, n_comp_ctrls
        logical,               intent(in)    :: sp_required
        call self%kill
        allocate(self%name,        source=trim(name)       )
        allocate(self%descr_short, source=trim(descr_short))
        allocate(self%descr_long,  source=trim(descr_long) )
        allocate(self%executable,  source=trim(executable) )
        if( n_img_ios    > 0 ) allocate(self%img_ios(n_img_ios)      )
        if( n_parm_ios   > 0 ) allocate(self%parm_ios(n_parm_ios)    )
        if( n_alt_ios    > 0 ) allocate(self%alt_ios(n_alt_ios)      )
        if( n_srch_ctrls > 0 ) allocate(self%srch_ctrls(n_srch_ctrls))
        if( n_filt_ctrls > 0 ) allocate(self%filt_ctrls(n_filt_ctrls))
        if( n_mask_ctrls > 0 ) allocate(self%mask_ctrls(n_mask_ctrls))
        if( n_comp_ctrls > 0 ) allocate(self%comp_ctrls(n_comp_ctrls))
        self%sp_required = sp_required
        self%exists      = .true.
    end subroutine new

    subroutine set_input_1( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                       intent(in)    :: required
        real,                          intent(in)    :: default_value
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input_1'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                allocate(arr(i)%key,               source=trim(key))
                allocate(arr(i)%keytype,           source=trim(keytype))
                allocate(arr(i)%descr_short,       source=trim(descr_short))
                allocate(arr(i)%descr_long,        source=trim(descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(descr_placeholder))
                arr(i)%required = required
                if( .not. arr(i)%required ) arr(i)%rval_default = default_value
            end subroutine set

    end subroutine set_input_1

    subroutine set_input_2( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                       intent(in)    :: required
        character(len=*),              intent(in)    :: default_value
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input_2'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                allocate(arr(i)%key,               source=trim(key))
                allocate(arr(i)%keytype,           source=trim(keytype))
                allocate(arr(i)%descr_short,       source=trim(descr_short))
                allocate(arr(i)%descr_long,        source=trim(descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(descr_placeholder))
                arr(i)%required = required
                if( .not. arr(i)%required ) allocate(arr(i)%cval_default, source=trim(default_value))
            end subroutine set

    end subroutine set_input_2

    subroutine set_input_3( self, which, i, param )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        type(simple_input_param),      intent(in)    :: param
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input_2'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                allocate(arr(i)%key,               source=trim(param%key))
                allocate(arr(i)%keytype,           source=trim(param%keytype))
                allocate(arr(i)%descr_short,       source=trim(param%descr_short))
                allocate(arr(i)%descr_long,        source=trim(param%descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(param%descr_placeholder))
                arr(i)%required = param%required
            end subroutine set

    end subroutine set_input_3

    subroutine print_ui( self )
        use simple_chash, only: chash
        use simple_ansi_ctrls
        class(simple_program), intent(in) :: self
        type(chash) :: ch
        integer     :: i
        write(*,'(a)') ''
        write(*,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name)
        call ch%push('descr_short', self%descr_short)
        call ch%push('descr_long',  self%descr_long)
        call ch%push('executable',  self%executable)
        call ch%print_key_val_pairs
        call ch%kill
        write(*,'(a)') ''
        write(*,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash(self%img_ios)
        write(*,'(a)') ''
        write(*,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash(self%parm_ios)
        write(*,'(a)') ''
        write(*,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash(self%alt_ios)
        write(*,'(a)') ''
        write(*,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%srch_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%filt_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash(self%mask_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash(self%comp_ctrls)

        contains

            subroutine print_param_hash( arr )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        write(*,'(a,1x,i3)') '>>> PARAMETER #', i
                        call ch%new(6)
                        call ch%push('key',               arr(i)%key)
                        call ch%push('keytype',           arr(i)%keytype)
                        call ch%push('descr_short',       arr(i)%descr_short)
                        call ch%push('descr_long',        arr(i)%descr_long)
                        call ch%push('descr_placeholder', arr(i)%descr_placeholder)
                        if( arr(i)%required )then
                            call ch%push('required', 'T')
                        else
                            call ch%push('required', 'F')
                        endif
                        call ch%print_key_val_pairs
                        call ch%kill
                    end do
                endif
            end subroutine print_param_hash

    end subroutine print_ui

    subroutine print_cmdline( self )
        use simple_chash, only: chash
        use simple_ansi_ctrls
        use simple_strings, only: lexSort
        class(simple_program), intent(in) :: self
        type(chash) :: ch
        integer     :: i
        logical     :: l_distr_exec
        l_distr_exec = self%executable .eq. 'simple_distr_exec'

        write(*,'(a)') format_str('USAGE', C_UNDERLINED)
        if( l_distr_exec )then
            write(*,'(a)') format_str('bash-3.2$ simple_distr_exec prg='//self%name//' key1=val1 key2=val2 ...', C_ITALIC)
        else
            write(*,'(a)') format_str('bash-3.2$ simple_exec prg='//self%name//' key1=val1 key2=val2 ...', C_ITALIC)
        endif
        write(*,'(a)') 'Required input parameters in ' // format_str('bold', C_BOLD) // ' (ensure terminal support)'

        if( allocated(self%img_ios) )    write(*,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash(self%img_ios)

        if( allocated(self%parm_ios) )   write(*,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash(self%parm_ios)

        if( allocated(self%alt_ios) )    write(*,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash(self%alt_ios)

        if( allocated(self%srch_ctrls) ) write(*,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%srch_ctrls)

        if( allocated(self%filt_ctrls) ) write(*,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%filt_ctrls)

        if( allocated(self%mask_ctrls) ) write(*,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash(self%mask_ctrls)

        if( allocated(self%comp_ctrls) ) write(*,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash(self%comp_ctrls)

        contains

            subroutine print_param_hash( arr )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                character(len=KEYLEN),    allocatable :: sorted_keys(:), rearranged_keys(:)
                logical,                  allocatable :: required(:)
                integer,                  allocatable :: inds(:)
                integer :: i, nparams, nreq, iopt
                if( allocated(arr) )then
                    nparams = size(arr)
                    call ch%new(nparams)
                    allocate(sorted_keys(nparams), rearranged_keys(nparams), required(nparams))
                    do i=1,nparams
                        call ch%push(arr(i)%key, arr(i)%descr_short//'; '//arr(i)%descr_placeholder)
                        sorted_keys(i) = arr(i)%key
                        required(i)    = arr(i)%required
                    end do
                    call lexSort(sorted_keys, inds=inds)
                    required = required(inds)
                    if( any(required) )then
                        ! fish out the required ones
                        nreq = 0
                        do i=1,nparams
                            if( required(i) )then
                                nreq = nreq + 1
                                rearranged_keys(nreq) = sorted_keys(i)
                            endif
                        enddo
                        ! fish out the optional ones
                        iopt = nreq
                        do i=1,nparams
                            if( .not. required(i) )then
                                iopt = iopt + 1
                                rearranged_keys(iopt) = sorted_keys(i)
                            endif
                        end do
                        ! replace string array
                        sorted_keys = rearranged_keys
                        ! modify logical mask
                        required(:nreq)     = .true.
                        required(nreq + 1:) = .false.
                    endif
                    call ch%print_key_val_pairs(sorted_keys, mask=required)
                    call ch%kill
                    deallocate(sorted_keys, required)
                endif
            end subroutine print_param_hash

    end subroutine print_cmdline

    subroutine print_prg_descr_long( self )
        class(simple_program), intent(in) :: self
        write(*,'(a)') self%descr_long
    end subroutine print_prg_descr_long

    subroutine write2json( self )
        use json_module
        use simple_strings, only: int2str
        class(simple_program), intent(in) :: self
        type(json_core)           :: json
        type(json_value), pointer :: pjson, program
        ! JSON init
        call json%initialize()
        call json%create_object(pjson,'')
        call json%create_object(program, trim(self%name))
        call json%add(pjson, program)
        ! program section
        call json%add(program, 'name',        self%name)
        call json%add(program, 'descr_short', self%descr_short)
        call json%add(program, 'descr_long',  self%descr_long)
        call json%add(program, 'executable',  self%executable)
        ! all sections
        call create_section( 'image input/output',     self%img_ios )
        call create_section( 'parameter input/output', self%parm_ios )
        call create_section( 'alternative inputs',     self%alt_ios )
        call create_section( 'search controls',        self%srch_ctrls )
        call create_section( 'filter controls',        self%filt_ctrls )
        call create_section( 'mask controls',          self%mask_ctrls )
        call create_section( 'computer controls',      self%comp_ctrls )
        ! write & clean
        call json%print(pjson, trim(adjustl(self%name))//'.json')
        if( json%failed() )then
            write(*,*) 'json input/output error for program: ', trim(self%name)
            stop
        endif
        call json%destroy(pjson)

        contains

            subroutine create_section( name, arr )
                use simple_strings, only: split, parsestr
                character(len=*),          intent(in) :: name
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                type(json_value), pointer :: entry, section, options
                character(len=STDLEN)     :: options_str, before
                character(len=KEYLEN)     :: args(8)
                integer                   :: i, j, sz, nargs
                logical :: found, param_is_multi, param_is_binary, exception
                call json%create_array(section, trim(name))
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        call json%create_object(entry, trim(arr(i)%key))
                        call json%add(entry, 'key', trim(arr(i)%key))
                        call json%add(entry, 'keytype', trim(arr(i)%keytype))
                        call json%add(entry, 'descr_short', trim(arr(i)%descr_short))
                        call json%add(entry, 'descr_long', trim(arr(i)%descr_long))
                        call json%add(entry, 'descr_placeholder', trim(arr(i)%descr_placeholder))
                        call json%add(entry, 'required', arr(i)%required)
                        param_is_multi  = trim(arr(i)%keytype).eq.'multi'
                        param_is_binary = trim(arr(i)%keytype).eq.'binary'
                        if( param_is_multi .or. param_is_binary )then
                            options_str = trim(arr(i)%descr_placeholder)
                            call split( options_str, '(', before )
                            call split( options_str, ')', before )
                            call parsestr(before, '|', args, nargs)
                            exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 3)
                            if( exception )then
                                write(*,*)'Poorly formatted options string for entry ', trim(arr(i)%key)
                                write(*,*)trim(arr(i)%descr_placeholder)
                                stop
                            endif
                            call json%add(entry, 'options', args(1:nargs))
                            do j = 1, nargs
                                call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                            enddo
                        endif
                        call json%add(section, entry)
                    enddo
                endif
                call json%add(pjson, section)
            end subroutine create_section

    end subroutine write2json

    function get_required_keys( self ) result( keys )
        class(simple_program), intent(in) :: self
        type(str4arr), allocatable :: keys(:)
        integer :: nreq, ireq
        ! count # required
        nreq = nreq_counter(self%img_ios) + nreq_counter(self%parm_ios)   + nreq_counter(self%alt_ios)    +&
        &nreq_counter(self%srch_ctrls)    + nreq_counter(self%filt_ctrls) + nreq_counter(self%mask_ctrls) +&
        &nreq_counter(self%comp_ctrls)
        ! extract keys
        if( nreq > 0 )then
            allocate(keys(nreq))
            ireq = 0
            call key_extractor(self%img_ios)
            call key_extractor(self%parm_ios)
            call key_extractor(self%alt_ios)
            call key_extractor(self%srch_ctrls)
            call key_extractor(self%filt_ctrls)
            call key_extractor(self%mask_ctrls)
            call key_extractor(self%comp_ctrls)
        endif

        contains

            function nreq_counter( arr ) result( nreq )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                integer :: nreq, i
                nreq = 0
                if( allocated(arr) )then
                    do i=1,size(arr)
                        if( arr(i)%required ) nreq = nreq + 1
                    end do
                endif
            end function nreq_counter

            subroutine key_extractor( arr )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        if( arr(i)%required )then
                            ireq = ireq + 1
                            allocate(keys(ireq)%str, source=arr(i)%key)
                        endif
                    end do
                endif
            end subroutine key_extractor

    end function get_required_keys

    logical function is_distr( self )
        use simple_strings, only: str_has_substr
        class(simple_program), intent(in) :: self
        is_distr = str_has_substr(self%executable, 'distr')
    end function is_distr

    logical function requires_sp_project( self )
        class(simple_program), intent(in) :: self
        requires_sp_project = self%sp_required
    end function requires_sp_project

    subroutine kill( self )
        class(simple_program), intent(inout) :: self
        integer :: i, sz
        if( self%exists )then
            deallocate(self%name, self%descr_short, self%descr_long, self%executable)
            call dealloc_field(self%img_ios)
            call dealloc_field(self%parm_ios)
            call dealloc_field(self%alt_ios)
            call dealloc_field(self%srch_ctrls)
            call dealloc_field(self%filt_ctrls)
            call dealloc_field(self%mask_ctrls)
            call dealloc_field(self%comp_ctrls)
            self%exists = .false.
        endif

        contains

            subroutine dealloc_field( arr )
                type(simple_input_param), allocatable, intent(inout) :: arr(:)
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        if( allocated(arr(i)%key)               ) deallocate(arr(i)%key              )
                        if( allocated(arr(i)%keytype)           ) deallocate(arr(i)%keytype          )
                        if( allocated(arr(i)%descr_short)       ) deallocate(arr(i)%descr_short      )
                        if( allocated(arr(i)%descr_long)        ) deallocate(arr(i)%descr_long       )
                        if( allocated(arr(i)%descr_placeholder) ) deallocate(arr(i)%descr_placeholder)
                    end do
                    deallocate(arr)
                endif
            end subroutine dealloc_field

    end subroutine kill

end module simple_user_interface
