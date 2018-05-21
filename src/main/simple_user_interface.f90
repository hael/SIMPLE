module simple_user_interface
include 'simple_lib.f08'
implicit none

public :: simple_program, make_user_interface, get_prg_ptr, list_distr_prgs_in_ui, list_shmem_prgs_in_ui
private

logical, parameter :: DEBUG = .false.

type simple_input_param
    character(len=:), allocatable :: key
    character(len=:), allocatable :: keytype ! (binary|multi|num|str|file|dir)
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
    ! sp_project required flag
    logical :: sp_required = .true.
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
    procedure          :: get_name
    procedure          :: get_nrequired_keys
    procedure          :: get_required_keys
    procedure          :: is_distr
    procedure          :: requires_sp_project
    procedure, private :: kill
end type simple_program

! declare protected program specifications here
type(simple_program), target :: center
type(simple_program), target :: cluster2D
type(simple_program), target :: cluster2D_stream
type(simple_program), target :: cluster3D
type(simple_program), target :: cluster3D_refine
type(simple_program), target :: cluster_cavgs
type(simple_program), target :: convert
type(simple_program), target :: ctf_estimate
type(simple_program), target :: ctfops
type(simple_program), target :: extract
type(simple_program), target :: filter
type(simple_program), target :: fsc
type(simple_program), target :: info_image
type(simple_program), target :: info_stktab
type(simple_program), target :: initial_3Dmodel
type(simple_program), target :: import_boxes
type(simple_program), target :: import_cavgs
type(simple_program), target :: import_movies
type(simple_program), target :: import_particles
type(simple_program), target :: make_cavgs
type(simple_program), target :: make_oris
type(simple_program), target :: make_pickrefs
type(simple_program), target :: mask
type(simple_program), target :: motion_correct
type(simple_program), target :: motion_correct_tomo
type(simple_program), target :: new_project
type(simple_program), target :: normalize_
type(simple_program), target :: orisops
type(simple_program), target :: oristats
type(simple_program), target :: pick
type(simple_program), target :: postprocess
type(simple_program), target :: powerspecs
type(simple_program), target :: preprocess
type(simple_program), target :: preprocess_stream
type(simple_program), target :: print_fsc
type(simple_program), target :: print_magic_boxes
type(simple_program), target :: print_project_info
type(simple_program), target :: reproject
type(simple_program), target :: reconstruct3D
type(simple_program), target :: refine3D
type(simple_program), target :: refine3D_init
type(simple_program), target :: scale
type(simple_program), target :: scale_project
type(simple_program), target :: select_
type(simple_program), target :: shift
type(simple_program), target :: simulate_movie
type(simple_program), target :: simulate_noise
type(simple_program), target :: simulate_particles
type(simple_program), target :: simulate_subtomogram
type(simple_program), target :: stack
type(simple_program), target :: stackops
type(simple_program), target :: symsrch
type(simple_program), target :: tseries_track
type(simple_program), target :: update_project
type(simple_program), target :: vizoris
type(simple_program), target :: volops

! declare common params here, with name same as flag
type(simple_input_param) :: angerr
type(simple_input_param) :: astigtol
type(simple_input_param) :: bfac
type(simple_input_param) :: box
type(simple_input_param) :: clip
type(simple_input_param) :: cs
type(simple_input_param) :: ctf
type(simple_input_param) :: dferr
type(simple_input_param) :: e1, e2, e3
type(simple_input_param) :: eo
type(simple_input_param) :: fraca
type(simple_input_param) :: frcs
type(simple_input_param) :: fromf
type(simple_input_param) :: job_memory_per_task
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
type(simple_input_param) :: mkdir_
type(simple_input_param) :: msk
type(simple_input_param) :: mskfile
type(simple_input_param) :: mw
type(simple_input_param) :: ncls
type(simple_input_param) :: neigh
type(simple_input_param) :: neg
type(simple_input_param) :: nsig
type(simple_input_param) :: nspace
type(simple_input_param) :: nparts
type(simple_input_param) :: nptcls
type(simple_input_param) :: numlen
type(simple_input_param) :: nthr
type(simple_input_param) :: objfun
type(simple_input_param) :: oritab
type(simple_input_param) :: oritab2
type(simple_input_param) :: oritype
type(simple_input_param) :: outer
type(simple_input_param) :: outfile
type(simple_input_param) :: outstk
type(simple_input_param) :: outvol
type(simple_input_param) :: pcontrast
type(simple_input_param) :: phaseplate
type(simple_input_param) :: pgrp
type(simple_input_param) :: projfile
type(simple_input_param) :: projname
type(simple_input_param) :: pspecsz
type(simple_input_param) :: qsys_partition
type(simple_input_param) :: qsys_qos
type(simple_input_param) :: qsys_reservation
type(simple_input_param) :: remap_cls
type(simple_input_param) :: shellw
type(simple_input_param) :: sherr
type(simple_input_param) :: smpd
type(simple_input_param) :: startit
type(simple_input_param) :: stk
type(simple_input_param) :: stktab
type(simple_input_param) :: time_per_image
type(simple_input_param) :: tof
type(simple_input_param) :: trs
type(simple_input_param) :: update_frac
type(simple_input_param) :: user_account
type(simple_input_param) :: user_email
type(simple_input_param) :: user_project
type(simple_input_param) :: weights2D
type(simple_input_param) :: weights3D

interface set_param
    module procedure set_param_1
    module procedure set_param_2
end interface set_param

contains

    ! public class methods

    subroutine make_user_interface
        call set_common_params
        call new_center
        call new_cluster2D
        call new_cluster2D_stream
        call new_cluster3D
        call new_cluster3D_refine
        call new_cluster_cavgs
        call new_convert
        call new_ctf_estimate
        call new_ctfops
        call new_extract
        call new_filter
        call new_fsc
        call new_info_image
        call new_info_stktab
        call new_initial_3Dmodel
        call new_import_boxes
        call new_import_cavgs
        call new_import_movies
        call new_import_particles
        call new_make_cavgs
        call new_make_oris
        call new_make_pickrefs
        call new_mask
        call new_motion_correct
        call new_motion_correct_tomo
        call new_new_project
        call new_normalize
        call new_orisops
        call new_oristats
        call new_pick
        call new_postprocess
        call new_powerspecs
        call new_preprocess
        call new_preprocess_stream
        call new_print_fsc
        call new_print_magic_boxes
        call new_print_project_info
        call new_reproject
        call new_reconstruct3D
        call new_refine3D
        call new_refine3D_init
        call new_scale
        call new_scale_project
        call new_select_
        call new_shift
        call new_simulate_movie
        call new_simulate_noise
        call new_simulate_particles
        call new_simulate_subtomogram
        call new_stack
        call new_stackops
        call new_symsrch
        call new_tseries_track
        call new_update_project
        call new_vizoris
        call new_volops
        ! ...
        if( DEBUG ) print *, '***DEBUG::simple_user_interface; make_user_interface, DONE'
    end subroutine make_user_interface

    subroutine get_prg_ptr( which_program, ptr2prg )
        character(len=*), intent(in)  :: which_program
        class(simple_program), pointer :: ptr2prg
        select case(trim(which_program))
            case('center')
                ptr2prg => center
            case('cluster2D')
                ptr2prg => cluster2D
            case('cluster2D_stream')
                ptr2prg => cluster2D_stream
            case('cluster3D')
                ptr2prg => cluster3D
            case('cluster3D_refine')
                ptr2prg => cluster3D_refine
            case('cluster_cavgs')
                ptr2prg => cluster_cavgs
            case('convert')
                ptr2prg => convert
            case('ctf_estimate')
                ptr2prg => ctf_estimate
            case('ctfops')
                ptr2prg => ctfops
            case('extract')
                ptr2prg => extract
            case('filter')
                ptr2prg => filter
            case('fsc')
                ptr2prg => fsc
            case('info_image')
                ptr2prg => info_image
            case('info_stktab')
                ptr2prg => info_stktab
            case('initial_3Dmodel')
                ptr2prg => initial_3Dmodel
            case('import_boxes')
                ptr2prg => import_boxes
            case('import_cavgs')
                ptr2prg => import_cavgs
            case('import_movies')
                ptr2prg => import_movies
            case('import_particles')
                ptr2prg => import_particles
            case('make_cavgs')
                ptr2prg => make_cavgs
            case('make_oris')
                ptr2prg => make_oris
            case('make_pickrefs')
                ptr2prg => make_pickrefs
            case('mask')
                ptr2prg => mask
            case('motion_correct')
                ptr2prg => motion_correct
            case('motion_correct_tomo')
                ptr2prg => motion_correct_tomo
            case('new_project')
                ptr2prg => new_project
            case('normalize')
                ptr2prg => normalize_
            case('orisops')
                ptr2prg => orisops
            case('oristats')
                ptr2prg => oristats
            case('pick')
                ptr2prg => pick
            case('postprocess')
                ptr2prg => postprocess
            case('powerspecs')
                ptr2prg => powerspecs
            case('preprocess')
                ptr2prg => preprocess
            case('preprocess_stream')
                ptr2prg => preprocess_stream
            case('print_fsc')
                ptr2prg => print_fsc
            case('print_magic_boxes')
                ptr2prg => print_magic_boxes
            case('print_project_info')
                ptr2prg => print_project_info
            case('reproject')
                ptr2prg => reproject
            case('reconstruct3D')
                ptr2prg => reconstruct3D
            case('refine3D')
                ptr2prg => refine3D
            case('refine3D_init')
                ptr2prg => refine3D_init
            case('scale')
                ptr2prg => scale
            case('scale_project')
                ptr2prg => scale_project
            case('select')
                ptr2prg => select_
            case('shift')
                ptr2prg => shift
            case('simulate_movie')
                ptr2prg => simulate_movie
            case('simulate_noise')
                ptr2prg => simulate_noise
            case('simulate_particles')
                ptr2prg => simulate_particles
            case('simulate_subtomogram')
                ptr2prg => simulate_subtomogram
            case('stack')
                ptr2prg => stack
            case('stackops')
                ptr2prg => stackops
            case('symsrch')
                ptr2prg => symsrch
            case('tseries_track')
                ptr2prg => tseries_track
            case('update_project')
                ptr2prg => update_project
            case('vizoris')
                ptr2prg => vizoris
            case('volops')
                ptr2prg => volops
            case DEFAULT
                ptr2prg => null()
        end select
    end subroutine get_prg_ptr

    subroutine list_distr_prgs_in_ui
        write(*,'(A)') cluster2D%name
        write(*,'(A)') cluster2D_stream%name
        write(*,'(A)') cluster3D%name
        write(*,'(A)') cluster3D_refine%name
        write(*,'(A)') ctf_estimate%name
        write(*,'(A)') initial_3Dmodel%name
        write(*,'(A)') make_cavgs%name
        write(*,'(A)') motion_correct%name
        write(*,'(A)') motion_correct_tomo%name
        write(*,'(A)') pick%name
        write(*,'(A)') powerspecs%name
        write(*,'(A)') preprocess%name
        write(*,'(A)') preprocess_stream%name
        write(*,'(A)') reconstruct3D%name
        write(*,'(A)') refine3D%name
        write(*,'(A)') refine3D_init%name
        write(*,'(A)') scale_project%name
        write(*,'(A)') symsrch%name
        write(*,'(A)') tseries_track%name
        stop
    end subroutine list_distr_prgs_in_ui

    subroutine list_shmem_prgs_in_ui
        write(*,'(A)') center%name
        write(*,'(A)') cluster_cavgs%name
        write(*,'(A)') convert%name
        write(*,'(A)') ctfops%name
        write(*,'(A)') extract%name
        write(*,'(A)') filter%name
        write(*,'(A)') fsc%name
        write(*,'(A)') info_image%name
        write(*,'(A)') info_stktab%name
        write(*,'(A)') import_boxes%name
        write(*,'(A)') import_cavgs%name
        write(*,'(A)') import_movies%name
        write(*,'(A)') import_particles%name
        write(*,'(A)') make_oris%name
        write(*,'(A)') make_pickrefs%name
        write(*,'(A)') mask%name
        write(*,'(A)') new_project%name
        write(*,'(A)') normalize_%name
        write(*,'(A)') orisops%name
        write(*,'(A)') oristats%name
        write(*,'(A)') postprocess%name
        write(*,'(A)') print_fsc%name
        write(*,'(A)') print_magic_boxes%name
        write(*,'(A)') print_project_info%name
        write(*,'(A)') reproject%name
        write(*,'(A)') select_%name
        write(*,'(A)') shift%name
        write(*,'(A)') simulate_movie%name
        write(*,'(A)') simulate_noise%name
        write(*,'(A)') simulate_particles%name
        write(*,'(A)') simulate_subtomogram%name
        write(*,'(A)') scale%name
        write(*,'(A)') stack%name
        write(*,'(A)') stackops%name
        write(*,'(A)') update_project%name
        write(*,'(A)') vizoris%name
        write(*,'(A)') volops%name
        stop
    end subroutine list_shmem_prgs_in_ui

    ! private class methods

    subroutine set_common_params
        call set_param(projfile,      'projfile',      'file',   'Project file', 'SIMPLE projectfile', 'e.g. myproject.simple', .true., 'myproject.simple')
        call set_param(stk,           'stk',           'file',   'Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., 'stk.mrc')
        call set_param(stktab,        'stktab',        'file',   'List of per-micrograph particle stacks', 'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false., 'stktab.txt')
        call set_param(ctf,           'ctf',           'multi',  'CTF status', 'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){no}',&
        &'(yes|no|flip){no}', .true., 'no')
        call set_param(smpd,          'smpd',          'num',    'Sampling distance', 'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms', .true., 1.0)
        call set_param(phaseplate,    'phaseplate',    'binary', 'Phase-plate images', 'Images obtained with Volta phase-plate(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(deftab,        'deftab',        'file',   'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format with dfx, dfy and angast values',&
        &'.simple|.txt parameter file', .false., 'deftab'//trim(METADATA_EXT))
        call set_param(oritab,        'oritab',        'file',   'Orientation and CTF parameter file', 'Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false., 'oritab'//trim(METADATA_EXT))
        call set_param(oritab2,       'oritab2',        'file',   '2nd orientation and CTF parameter file', '2nd orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false., 'oritab2'//trim(METADATA_EXT))
        call set_param(outfile,       'outfile',       'file',   'Output orientation and CTF parameter file', 'Output Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false., 'outfile'//trim(METADATA_EXT))
        call set_param(startit,       'startit',       'num',    'First iteration', 'Index of first iteration when starting from a previous run', 'start iterations from here', .false., 1.0)
        call set_param(trs,           'trs',           'num',    'Maximum translational shift', 'Maximum half-width for bund-constrained search of rotational origin shifts',&
        &'max shift per iteration in pixels{5}', .false., 0.0)
        call set_param(maxits,        'maxits',        'num',    'Max iterations', 'Maximum number of iterations', 'Max # iterations', .false., 100.)
        call set_param(hp,            'hp',            'num',    'High-pass limit', 'High-pass resolution limit', 'high-pass limit in Angstroms', .false., 100.)
        call set_param(lp,            'lp',            'num',    'Low-pass limit', 'Low-pass resolution limit', 'low-pass limit in Angstroms', .false., 20.)
        call set_param(msk,           'msk',           'num',    'Mask radius', 'Mask radius in pixels for application of a soft-edged circular mask to remove background noise', 'mask radius in pixels', .true., 0.)
        call set_param(inner,         'inner',         'num',    'Inner mask radius', 'Inner mask radius for omitting unordered cores of particles with high radial symmetry, typically icosahedral viruses',&
        &'inner mask radius in pixels', .false., 0.)
        call set_param(ncls,          'ncls',          'num',    'Number of 2D clusters', 'Number of groups to sort the particles &
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
        'fraction of particles(0.1-0.9){1.0}', .false., 1.0)
        call set_param(mskfile,       'mskfile',       'file',   'Input mask file', 'Input mask file to apply to reference volume(s) before projection', 'e.g. automask.mrc from postprocess', .false., 'mskfile.mrc')
        call set_param(pgrp,          'pgrp',          'str',    'Point-group symmetry', 'Point-group symmetry of particle(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .true., 'c1')
        call set_param(nspace,        'nspace',        'num',    'Number of projection directions', 'Number of projection directions &
        &used', '# projections', .false., 2500.)
        call set_param(objfun,        'objfun',        'binary', 'Objective function', 'Objective function(cc|ccres){cc}', '(cc|ccres){cc}', .false., 'cc')
        call set_param(weights2D,     'weights2D',     'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(weights3D,     'weights3D',     'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding re-projection(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(remap_cls,     'remap_cls',     'binary', 'Whether to remap 2D clusters', 'Whether to remap the number of 2D clusters(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(kv,            'kv',            'num',    'Acceleration voltage', 'Acceleration voltage in kV', 'in kV', .false., 300.)
        call set_param(lplim_crit,    'lplim_crit',    'num',    'Low-pass limit FSC criterion', 'FSC criterion for determining the low-pass limit(0.143-0.5){0.3}',&
        &'low-pass FSC criterion(0.143-0.5){0.3}', .false., 0.3)
        call set_param(cs,            'cs',            'num',    'Spherical aberration', 'Spherical aberration constant(in mm){2.7}', 'in mm{2.7}', .false., 2.7)
        call set_param(fraca,         'fraca',         'num',    'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.1}', 'fraction{0.1}', .false., 0.1)
        call set_param(pspecsz,       'pspecsz',       'num',    'Size of power spectrum', 'Size of power spectrum in pixels', 'in pixels', .false., 512.)
        call set_param(dfmin,         'dfmin',         'num',    'Expected minimum defocus', 'Expected minimum defocus in microns{0.5}', 'in microns{0.5}', .false., 0.5)
        call set_param(dfmax,         'dfmax',         'num',    'Expected maximum defocus', 'Expected minimum defocus in microns{5.0}', 'in microns{5.0}', .false., 5.0)
        call set_param(dfstep,        'dfstep',        'num',    'Defocus step size', 'Defocus step size for grid search in microns{0.05}', 'in microns{0.05}', .false., 0.05)
        call set_param(astigtol,      'astigtol',      'num',    'Expected astigmatism', 'expected (tolerated) astigmatism(in microns){0.05}', 'in microns',  .false., 0.05)
        call set_param(mw,            'mw',            'num',    'Molecular weight','Molecular weight in kDa', 'in kDa', .false., 0.)
        call set_param(mirr,          'mirr',          'multi',  'Perform mirroring', 'Whether to mirror and along which axis(no|x|y){no}', '(no|x|y){no}', .false., 'no')
        call set_param(bfac,          'bfac',          'num',    'B-factor for sharpening','B-factor for sharpening in Angstroms^2', 'B-factor in Angstroms^2', .false., 200.)
        call set_param(outvol,        'outvol',        'file',   'Output volume name', 'Output volume name', 'e.g. outvol.mrc', .false., '')
        call set_param(eo,            'eo',            'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){yes}', '(yes|no){yes}', .false., 'no')
        call set_param(job_memory_per_task, 'job_memory_per_task','str', 'Memory per part', 'Memory in MB per part in distributed execution{1600}', 'MB per part{1600}', .false., 1600.)
        call set_param(qsys_partition,'qsys_partition','str',    'Name of SLURM/PBS partition', 'Name of target partition of distributed computer system (SLURM/PBS)', 'give part name', .false., '')
        call set_param(qsys_qos,      'qsys_qos',      'str',    'Schedule priority', 'Job scheduling priority (SLURM/PBS)', 'give priority', .false., '')
        call set_param(qsys_reservation, 'qsys_reservation', 'str', 'Name of reserved partition', 'Name of reserved target partition of distributed computer system (SLURM/PBS)', 'give yourpart', .false., '')
        call set_param(box,            'box',          'num',    'Square image size','Square image size(in pixels)', '# pixels of box', .true., 0.)
        call set_param(nptcls,         'nptcls',       'num',    'Number of particles', 'Number of particle images', '# particles', .true., 0.)
        call set_param(outstk,         'outstk',       'file',   'Output stack name', 'Output images stack name', 'e.g. outstk.mrc', .false., '')
        call set_param(pcontrast,      'pcontrast',    'binary', 'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false., 'black')
        call set_param(clip,           'clip',         'num',    'Clipped box size', 'Target box size for clipping in pixels', 'in pixels', .false., 0.)
        call set_param(neg,            'neg',          'binary', 'Invert contrast','Invert contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(sherr,          'sherr',        'num',    'Shift error half-width', 'Uniform rotational origin shift error half-width(in pixels)', 'shift error in pixels', .false., 0.)
        call set_param(angerr,         'angerr',       'num',    'Rotation angle error half-width', 'Uniform rotation angle shift error half-width(in degrees)', 'rotation error in degrees', .false., 0.)
        call set_param(dferr,          'dferr',        'num',    'Underfocus error half-width',  'Uniform underfoucs error half-width(in microns)',  'defocus error in microns', .false., 1.)
        call set_param(oritype,        'oritype',      'multi',  'Oritype segment in project',  'Oritype segment in project(ptcl2D|ptcl3D|cls3D){ptcl3D}', '(ptcl2D|ptcl3D|cls3D){ptcl3D}', .false., 'ptcl3D')
        call set_param(e1,             'e1',           'num',    'Rotation along Phi',  'Phi Euler angle',   'in degrees', .false., 0.)
        call set_param(e2,             'e2',           'num',    'Rotation along Theta','Theat Euler angle', 'in degrees', .false., 0.)
        call set_param(e3,             'e3',           'num',    'Rotation along Psi',  'Psi Euler angle',   'in degrees', .false., 0.)
        call set_param(numlen,         'numlen',       'num',    'Length of number string', 'Length of number string', '# characters', .false., 5.0)
        call set_param(nsig,           'nsig',         'num',    'Number of sigmas for outlier removal', 'Number of standard deviations threshold for pixel outlier removal{6}', '# standard deviations{6}', .false., 6.)
        call set_param(fromf,          'fromf',        'num',    'First frame to include in subsum', 'First frame index to include in subsum', 'give index', .false., 1.)
        call set_param(tof,            'tof',          'num',    'Last frame to include in subsum', 'Last frame index to include in subsum', 'give index', .false., 1.)
        call set_param(neigh,          'neigh',        'binary', 'Neighbourhood refinement', 'Neighbourhood refinement(yes|no){yes}', '(yes|no){no}', .false., 'no')
        call set_param(projname,       'projname',     'str',    'Project name', 'Name of project to create ./myproject/myproject.simple file for',&
        &'e.g. to create ./myproject/myproject.simple', .true., '')
        call set_param(user_email,     'user_email',   'str',    'Your e-mail address', 'Your e-mail address', 'e.g. myname@uni.edu', .false., '')
        call set_param(time_per_image, 'time_per_image', 'num', 'Time per image', 'Estimated time per image in seconds for forecasting total execution time{100}', 'in seconds{100}', .false., 100.)
        call set_param(user_account,   'user_account', 'str',    'User account name in SLURM/PBS', 'User account name in SLURM/PBS system', 'e.g. Account084', .false., '')
        call set_param(user_project,   'user_project', 'str',    'User project name in SLURM/PBS', 'User project name in SLURM/PBS system', 'e.g. Project001', .false., '')
        call set_param(frcs,           'frcs',         'str',    'Projection FRCs file', 'Projection FRCs file', 'e.g. frcs.bin', .false., '')
        call set_param(mkdir_,         'mkdir',        'binary', 'Make auto-named dir for output', 'Make auto-named consequtively numbered dir for output(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(shellw,         'shellw',       'binary', 'B-factor weighted reconstruction', 'Whether to perform B-factor weighted reconstruction(yes|no){no}',  '(yes|no){no}',  .false., 'no')
        if( DEBUG ) print *, '***DEBUG::simple_user_interface; set_common_params, DONE'
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

    ! TEMPLATE
    ! INPUT PARAMETER SPECIFICATIONS
    ! image input/output
    ! <empty>
    ! parameter input/output
    ! <empty>
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

    subroutine new_center
        ! PROGRAM SPECIFICATION
        call center%new(&
        &'center',&                    ! name
        &'Center volume',&             ! descr_short
        &'is a program for centering a volume and mapping the shift parameters back to the particle images',& ! descr_long
        &'simple_exec',&               ! executable
        &1, 3, 0, 0, 1, 0, 1, .false.) ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call center%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to center', &
        & 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call center%set_input('parm_ios', 1, smpd)
        call center%set_input('parm_ios', 2, oritab)
        call center%set_input('parm_ios', 3, outfile)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call center%set_input('filt_ctrls', 1, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        ! <empty>
        ! computer controls
        call center%set_input('comp_ctrls', 1, nthr)
    end subroutine new_center

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',&                 ! descr_long
        &'simple_distr_exec',&                                                  ! executable
        &1, 1, 0, 11, 8, 2, 2, .true.)                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster2D%set_input('img_ios', 1, 'refs', 'file', 'Initial references',&
        &'Initial 2D references used to bootstrap the search', 'xxx.mrc file with references', .false., 'refs.mrc')
        ! parameter input/output
        call cluster2D%set_input('parm_ios', 1, mkdir_)
        ! alternative inputs
        ! <empty>
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
        call cluster2D%set_input('srch_ctrls',10, 'bfac', 'num', 'Correlation B-factor','B-factor for the objective function in Angstroms^2', 'B-factor in Angstroms^2(>0.0){200}', .false., 200.)
        call cluster2D%set_input('srch_ctrls',11, 'objfun','binary', 'Objective function', 'Objective function(cc|ccres){ccres}', '(cc|ccres){ccres}', .false., 'ccres')
        ! filter controls
        call cluster2D%set_input('filt_ctrls', 1, hp)
        call cluster2D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false., 20.)
        call cluster2D%set_input('filt_ctrls', 4, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &few iterations of search, before the automatic scheme kicks in. Also controls the degree of downsampling in the first &
        &phase', 'initial low-pass limit in Angstroms', .false., 15.)
        call cluster2D%set_input('filt_ctrls', 5, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.)
        call cluster2D%set_input('filt_ctrls', 6, 'match_filt', 'binary', 'Matched filter', 'Filter to maximize the signal-to-noise &
        &ratio (SNR) in the presence of additive stochastic noise. Sometimes causes over-fitting and needs to be turned off(yes|no){yes}',&
        '(yes|no){yes}', .false., 'yes')
        call cluster2D%set_input('filt_ctrls', 7, weights2D)
        call cluster2D%set_input('filt_ctrls', 8, shellw)
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
        &'simple_distr_exec',&                                                                                            ! executable
        &0, 2, 0, 6, 5, 2, 2, .true.)                                                                                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster2D_stream%set_input('parm_ios', 1, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preprocess_stream', .true., '')
        call cluster2D_stream%set_input('parm_ios', 2, mkdir_)
        ! alternative inputs
        ! <empty>
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

    subroutine new_cluster3D
        ! PROGRAM SPECIFICATION
        call cluster3D%new(&
        &'cluster3D',& ! name
        &'3D heterogeneity analysis',&                                             ! descr_short
        &'is a distributed workflow for heterogeneity analysis by 3D clustering',& ! descr_long
        &'simple_distr_exec',&                                                     ! executable
        &0, 2, 0, 7, 6, 5, 2, .true.)                                              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster3D%set_input('parm_ios', 1, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to separate',&
        '# states to separate', .true., 2.0)
        call cluster3D%set_input('parm_ios', 2, mkdir_)
        ! alternative inputs
        ! <empty>
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
        call cluster3D%set_input('filt_ctrls', 5, eo)
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

    subroutine new_cluster3D_refine
        ! PROGRAM SPECIFICATION
        call cluster3D_refine%new(&
        &'cluster3D_refine',&                                                ! name
        &'cluster 3D refinement',&                                           ! descr_short
        &'is a distributed workflow based on probabilistic projection matching &
        &for refinement of 3D heterogeneity analysis by cluster3D ',&        ! descr_long
        &'simple_distr_exec',&                                               ! executable
        &2, 2, 0, 13, 7, 3, 2, .true.)                                       ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster3D_refine%set_input('img_ios', 1, 'msklist', 'file', 'List of mask files', 'List (.txt file) of mask files for the different states', 'e.g. mskfiles.txt', .true., '')
        call cluster3D_refine%set_input('img_ios', 2, 'vollist', 'file', 'List of reference volumes files', 'List (.txt file) of reference volumes for the different states', 'e.g. refvols.txt', .true., '')
        ! parameter input/output
        call cluster3D_refine%set_input('parm_ios', 1,  'state', 'num', 'State to refine', 'Index of state to refine', 'give state index', .false., 1.)
        call cluster3D_refine%set_input('parm_ios', 2, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster3D_refine%set_input('srch_ctrls', 1, nspace)
        call cluster3D_refine%set_input('srch_ctrls', 2, startit)
        call cluster3D_refine%set_input('srch_ctrls', 3, trs)
        call cluster3D_refine%set_input('srch_ctrls', 4, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster3D_refine%set_input('srch_ctrls', 5, maxits)
        call cluster3D_refine%set_input('srch_ctrls', 6, update_frac)
        call cluster3D_refine%set_input('srch_ctrls', 7, frac)
        call cluster3D_refine%set_input('srch_ctrls', 8, pgrp)
        call cluster3D_refine%set_input('srch_ctrls', 9, 'nnn', 'num', 'Number of nearest neighbours', 'Number of nearest projection direction &
        &neighbours in neigh=yes refinement', '# projection neighbours{10% of search space}', .false., 200.)
        call cluster3D_refine%set_input('srch_ctrls', 10, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to reconstruct',&
        '# states to reconstruct', .false., 1.0)
        call cluster3D_refine%set_input('srch_ctrls', 11, objfun)
        call cluster3D_refine%set_input('srch_ctrls', 12, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc|single|multi|greedy_single|greedy_multi|cluster|&
        &clustersym){no}', '(snhc|single|multi|greedy_single|greedy_multi|cluster|clustersym){single}', .false., 'single')
        call cluster3D_refine%set_input('srch_ctrls', 13, neigh)
        ! filter controls
        call cluster3D_refine%set_input('filt_ctrls', 1, hp)
        call cluster3D_refine%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster3D_refine%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call cluster3D_refine%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call cluster3D_refine%set_input('filt_ctrls', 5, lplim_crit)
        call cluster3D_refine%set_input('filt_ctrls', 6, 'eo', 'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){yes}', '(yes|no){yes}', .false., 'no')
        call cluster3D_refine%set_input('filt_ctrls', 7, 'weights3D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! mask controls
        call cluster3D_refine%set_input('mask_ctrls', 1, msk)
        call cluster3D_refine%set_input('mask_ctrls', 2, inner)
        call cluster3D_refine%set_input('mask_ctrls', 3, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge{10}', .false., 10.)
        ! computer controls
        call cluster3D_refine%set_input('comp_ctrls', 1, nparts)
        call cluster3D_refine%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster3D_refine

    subroutine new_cluster_cavgs
        ! PROGRAM SPECIFICATION
        call cluster_cavgs%new(&
        &'cluster_cavgs',&                                                         ! name
        &'Analysis of class averages with affinity propagation',&                  ! descr_short
        &'is a program for analyzing class averages with affinity propagation, &
        &in order to get a better understanding of the view distribution. The balance flag is used &
        &to apply a balancing restraint (on the class population). Adjust balance until you are &
        &satisfied with the shape of the histogram',&                              ! descr_long
        &'simple_exec',&                                                           ! executable
        &1, 2, 0, 2, 2, 1, 1, .true.)                                              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster_cavgs%set_input('img_ios', 1, 'stk', 'file', 'Stack of class averages', 'Stack of class averages', 'e.g. cavgs.mrc', .true., '')
        ! parameter input/output
        call cluster_cavgs%set_input('parm_ios', 1, smpd)
        call cluster_cavgs%set_input('parm_ios', 2, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster_cavgs%set_input('srch_ctrls', 1, 'balance', 'num', 'Max population for balance restraint', 'Max population for balance restraint', 'max # cluster members', .false., 0.)
        call cluster_cavgs%set_input('srch_ctrls', 2, objfun)
        ! filter controls
        call cluster_cavgs%set_input('filt_ctrls', 1, hp)
        call cluster_cavgs%set_input('filt_ctrls', 2, lp)
        lp%required = .true.
        ! mask controls
        call cluster_cavgs%set_input('mask_ctrls', 1, msk)
        ! computer controls
        call cluster_cavgs%set_input('comp_ctrls', 1, nthr)
    end subroutine new_cluster_cavgs

    subroutine new_convert
        ! PROGRAM SPECIFICATION
        call convert%new(&
        &'convert',&                                                    ! name
        &'Convert between SPIDER and MRC formats',&                     ! descr_short
        &'is a program for converting between SPIDER and MRC formats',& ! descr_long
        &'simple_exec',&                                                ! executable
        &2, 0, 2, 0, 0, 0, 0, .false.)                                  ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call convert%set_input('img_ios', 1, outvol)
        call convert%set_input('img_ios', 2, outstk)
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        call convert%set_input('alt_ios', 1, 'vol1', 'file', 'Volume', 'Volume to convert', &
        & 'input volume e.g. vol.spi', .false., '')
        call convert%set_input('alt_ios', 2, 'stk', 'file', 'Stack', 'Stack to convert',&
        & 'input stack e.g. imgs.spi', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_convert

    subroutine new_ctf_estimate
        ! PROGRAM SPECIFICATION
        call ctf_estimate%new(&
        &'ctf_estimate', &                                              ! name
        &'CTF parameter fitting',&                                      ! descr_short
        &'is a distributed SIMPLE workflow for CTF parameter fitting',& ! descr_long
        &'simple_distr_exec',&                                          ! executable
        &0, 1, 0, 4, 2, 0, 2, .true.)                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_estimate%set_input('parm_ios', 1, pspecsz)
        ! alternative inputs
        ! <empty>
        ! search controls
        call ctf_estimate%set_input('srch_ctrls', 1, dfmin)
        call ctf_estimate%set_input('srch_ctrls', 2, dfmax)
        call ctf_estimate%set_input('srch_ctrls', 3, dfstep)
        call ctf_estimate%set_input('srch_ctrls', 4, astigtol)
        ! filter controls
        call ctf_estimate%set_input('filt_ctrls', 1, lp)
        ctf_estimate%filt_ctrls(1)%required     = .false.
        call ctf_estimate%set_input('filt_ctrls', 2, hp)
        ctf_estimate%filt_ctrls(2)%required     = .false.
        ! mask controls
        ! <empty>
        ! computer controls
        call ctf_estimate%set_input('comp_ctrls', 1, nparts)
        call ctf_estimate%set_input('comp_ctrls', 2, nthr)
    end subroutine new_ctf_estimate

    subroutine new_ctfops
        ! PROGRAM SPECIFICATION
        call ctfops%new(&
        &'ctfops', &                                         ! name
        &'Apply CTF to stacked images',&                     ! descr_short
        &'is a program for applying CTF to stacked images',& ! descr long
        &'simple_exec',&                                     ! executable
        &2, 4, 0, 0, 2, 0, 1, .false.)                       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ctfops%set_input('img_ios', 1, stk)
        call ctfops%set_input('img_ios', 2, outstk)
        ! parameter input/output
        call ctfops%set_input('parm_ios', 1, smpd)
        call ctfops%set_input('parm_ios', 2, 'neg', 'binary', 'Invert contrast','Invert contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call ctfops%set_input('parm_ios', 3, oritab)
        call ctfops%set_input('parm_ios', 4, deftab)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ctfops%set_input('filt_ctrls', 1, ctf)
        call ctfops%set_input('filt_ctrls', 2, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        ! mask controls
        ! <empty>
        ! computer controls
        call ctfops%set_input('comp_ctrls', 1, nthr)
    end subroutine new_ctfops

    subroutine new_extract
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', &                                                           ! name
        &'Extract particle images from integrated movies',&                     ! descr_short
        &'is a program for extracting particle images from integrated movies',& ! descr long
        &'simple_exec',&                                                        ! executable
        &0, 3, 0, 0, 0, 0, 0, .true.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract%set_input('parm_ios', 1, 'box', 'num', 'Box size', 'Square box size in pixels', 'in pixels', .false., 0.)
        call extract%set_input('parm_ios', 2, pcontrast)
        call extract%set_input('parm_ios', 3, 'outside', 'binary', 'Extract outside boundaries', 'Extract boxes outside the micrograph boundaries(yes|no){no}', '(yes|no){no}', .false., 'no')
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

    subroutine new_filter
        ! PROGRAM SPECIFICATION
        call filter%new(&
        &'filter',&                                   ! name
        &'Filter stack/volume',&                      ! descr_short
        &'is a program for filtering stack/volume',&  ! descr_long
        &'simple_exec',&                              ! executable
        &2, 1, 2, 0, 8, 0, 1, .false.)                ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call filter%set_input('img_ios', 1, outstk)
        call filter%set_input('img_ios', 2, outvol)
        ! parameter input/output
        call filter%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        call filter%set_input('alt_ios', 1, 'stk',  'file', 'Stack to filter',  'Stack of images to filter', 'e.g. refs.mrc',     .false., '')
        call filter%set_input('alt_ios', 2, 'vol1', 'file', 'Volume to filter', 'Volume to filter',          'e.g. vol.mrc file', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        call filter%set_input('filt_ctrls', 1, lp)
        filter%filt_ctrls(1)%required = .false.
        call filter%set_input('filt_ctrls', 2, hp)
        call filter%set_input('filt_ctrls', 3, 'phrand', 'binary', 'Phase randomization', 'Fouirer phase randomization by white noise substitution(yes|no){no}', '(yes|no){no}', .false., 'no')
        call filter%set_input('filt_ctrls', 4, 'bfac', 'num', 'B-factor of Gaussian low-/high-pass filter','B-factor of Gaussian low-/high-pass filter in Angstroms^2', 'B-factor in Angstroms^2{0}', .false., 0.)
        call filter%set_input('filt_ctrls', 5, 'winsz', 'num', 'Half-window size', 'Half-window size(in pixels)', 'winsz in pixels', .false., 1.0)
        call filter%set_input('filt_ctrls', 6, 'width', 'num', 'Cosine low-pass filter falloff',&
        &'Number of cosine edge pixels of Fourier low-pass filter in pixels', '# pixels cosine edge', .false., 10.)
        call filter%set_input('filt_ctrls', 7, 'real_filter', 'binary', 'Real-space filter',&
        &'Real-space filter(yes|no){no}', '(yes|no){no}', .false., 'no')
        call filter%set_input('filt_ctrls', 8, frcs)
        ! mask controls
        ! <empty>
        ! computer controls
        call filter%set_input('comp_ctrls', 1, nthr)
    end subroutine new_filter

    subroutine new_fsc
        ! PROGRAM SPECIFICATION
        call fsc%new(&
        &'fsc', &                                                               ! name
        &'Calculate FSC between the two input volumes',&                        ! descr_short
        &'is a program for calculating the FSC between the two input volumes',& ! descr_long
        &'simple_exec',&                                                        ! executable
        &2, 1, 0, 0, 0, 2, 1, .false.)                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call fsc%set_input('img_ios', 1, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call fsc%set_input('img_ios', 2, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call fsc%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call fsc%set_input('mask_ctrls', 1, msk)
        call fsc%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call fsc%set_input('comp_ctrls', 1, nthr)
    end subroutine new_fsc

    subroutine new_info_image
        ! PROGRAM SPECIFICATION
        call info_image%new(&
        &'info_image', & ! name
        &'Print header information',&                                                          ! descr_short
        &'is a program for printing header information in MRC and SPIDER stacks and volumes',& ! descr_long
        &'simple_exec',&                                                                       ! executable
        &1, 2, 0, 0, 0, 0, 0, .false.)                                                         ! # entries in each group, requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call info_image%set_input('img_ios', 1, 'fname', 'file', 'Name of image file', 'Name of image file', 'xxx.mrc file', .true., '')
        ! parameter input/output
        call info_image%set_input('parm_ios', 1, 'stats', 'binary', 'Output statistics', 'Output statistics(yes|no){no}',             '(yes|no){no}', .false., 'no')
        call info_image%set_input('parm_ios', 2, 'vis',   'binary', 'Visualize image',   'Visualize image with gnuplot(yes|no){yes}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
    end subroutine new_info_image

    subroutine new_info_stktab
        ! PROGRAM SPECIFICATION
        call info_stktab%new(&
        &'info_stktab', & ! name
        &'Print stktab information',&                                            ! descr_short
        &'is a program for printing information about stktab (list of stacks)',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &1, 0, 0, 0, 0, 0, 0, .false.)                                           ! # entries in each group, requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        stktab%required = .true.
        call info_stktab%set_input('img_ios', 1, stktab)
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
    end subroutine new_info_stktab

    subroutine new_initial_3Dmodel
        ! PROGRAM SPECIFICATION
        call initial_3Dmodel%new(&
        &'initial_3Dmodel',&                                                            ! name
        &'3D ab initio model generation from class averages',&                          ! descr_short
        &'is a distributed workflow for generating an initial 3D model from class'&
        &' averages obtained with cluster2D',&                                          ! descr_long
        &'simple_distr_exec',&                                                          ! executable
        &0, 1, 0, 10, 5, 3, 2, .true.)                                                  ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call initial_3Dmodel%set_input('parm_ios', 1, mkdir_)
        ! alternative inputs
        ! <empty>
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
        call initial_3Dmodel%set_input('srch_ctrls', 10, 'nrestarts', 'num', 'Number of restarts', 'Number of program restarts to execute{1}', '# restarts{1}', .false., 1.0)
        ! filter controls
        call initial_3Dmodel%set_input('filt_ctrls', 1, hp)
        call initial_3Dmodel%set_input('filt_ctrls', 2, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms', .false., 0.)
        call initial_3Dmodel%set_input('filt_ctrls', 3, 'lpstop',  'num', 'Final low-pass limit',   'Final low-pass limit',   'low-pass limit in Angstroms', .false., 8.)
        call initial_3Dmodel%set_input('filt_ctrls', 4, eo)
        initial_3Dmodel%filt_ctrls(4)%cval_default = 'no'
        call initial_3Dmodel%set_input('filt_ctrls', 5, shellw)
        ! mask controls
        call initial_3Dmodel%set_input('mask_ctrls', 1, msk)
        call initial_3Dmodel%set_input('mask_ctrls', 2, inner)
        call initial_3Dmodel%set_input('mask_ctrls', 3, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false., 10.)
        ! computer controls
        call initial_3Dmodel%set_input('comp_ctrls', 1, nparts)
        call initial_3Dmodel%set_input('comp_ctrls', 2, nthr)
    end subroutine new_initial_3Dmodel

    subroutine new_import_boxes
        ! PROGRAM SPECIFICATION
        call import_boxes%new(&
        &'import_boxes',&                                  ! name
        &'Import EMAN box coordinates to SIMPLE project',& ! descr_short
        &'is a program for importing EMAN1.9 box coordinates to the project. The *box (text) files should be listed in boxtab',&
        &'simple_exec',&                                   ! executable
        &0, 1, 0, 0, 0, 0, 0, .true.)                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call import_boxes%set_input('parm_ios', 1, 'boxtab', 'file', 'List of box files', 'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .true., '')
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
    end subroutine new_import_boxes

    subroutine new_import_cavgs
        ! PROGRAM SPECIFICATION
        call import_cavgs%new(&
        &'import_cavgs',&                                        ! name
        &'Import class averages to SIMPLE project',&             ! descr_short
        &'is a program for importing class averages movies to the project',&
        &'simple_exec',&                                         ! executable
        &1, 1, 0, 0, 0, 0, 0, .true.)                            ! # entries in each group, requires sp_project
        call import_cavgs%set_input('img_ios', 1, 'stk', 'file', 'Stack of class averages',&
        &'Stack of class average images to import', 'e.g. cavgs.mrcs', .true., '')
        ! parameter input/output
        call import_cavgs%set_input('parm_ios', 1, smpd)
    end subroutine new_import_cavgs

    subroutine new_import_movies
        ! PROGRAM SPECIFICATION
        call import_movies%new(&
        &'import_movies',&                                       ! name
        &'Import movies to SIMPLE project',&                     ! descr_short
        &'is a program for importing DDD movies to the project. The movies can be located in any read-only location '&
        &'accessible to the project. If the movies contain only a single frame, they will be interpreted as motion-corrected '&
        &'and integrated. Box files (in EMAN format) can be imported along with the movies',&
        &'simple_exec',&                                         ! executable
        &1, 7, 0, 0, 0, 0, 0, .true.)                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call import_movies%set_input('img_ios', 1, 'filetab', 'file', 'List of movie files', 'List of movie files (*.mrcs) to import', 'e.g. movies.txt', .true., '')
        ! parameter input/output
        call import_movies%set_input('parm_ios', 1, smpd)
        call import_movies%set_input('parm_ios', 2, kv)
        import_movies%parm_ios(2)%required = .true.
        call import_movies%set_input('parm_ios', 3, cs)
        import_movies%parm_ios(3)%required = .true.
        call import_movies%set_input('parm_ios', 4, fraca)
        import_movies%parm_ios(4)%required = .true.
        call import_movies%set_input('parm_ios', 5, ctf)
        import_movies%parm_ios(5)%required = .true.
        call import_movies%set_input('parm_ios', 6, phaseplate)
        call import_movies%set_input('parm_ios', 7, 'boxtab', 'file', 'List of box files', 'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .false., '')
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
    end subroutine new_import_movies

    subroutine new_import_particles
        ! PROGRAM SPECIFICATION
        call import_particles%new(&
        &'import_particles',&                                       ! name
        &'Import particles to SIMPLE project',&                     ! descr_short
        &'is a program for importing extracted particle images to the project',&
        &'simple_exec',&                                            ! executable
        &0, 12, 2, 0, 0, 0, 0, .true.)                             ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call import_particles%set_input('parm_ios', 1, smpd)
        call import_particles%set_input('parm_ios', 2, kv)
        call import_particles%set_input('parm_ios', 3, cs)
        call import_particles%set_input('parm_ios', 4, fraca)
        call import_particles%set_input('parm_ios', 5, ctf)
        import_particles%parm_ios(5)%required = .true.
        call import_particles%set_input('parm_ios', 6, phaseplate)
        call import_particles%set_input('parm_ios', 7, oritab)
        call import_particles%set_input('parm_ios', 8, deftab)
        call import_particles%set_input('parm_ios', 9, 'plaintexttab', 'file', 'Plain text file of input parameters',&
        'Plain text file of tabulated per-particle input parameters: dfx, dfy, angast, phshift', 'e.g. params.txt', .false., '')
        call import_particles%set_input('parm_ios', 10, 'dfunit', 'binary', 'Underfocus unit', 'Underfocus unit(A|microns){microns}', '(A|microns){microns}', .false., 'microns')
        call import_particles%set_input('parm_ios', 11, 'angast', 'binary', 'Angle of astigmatism unit', 'Angle of astigmatism unit(radians|degrees){degrees}', '(radians|degrees){degrees}', .false., 'degrees')
        call import_particles%set_input('parm_ios', 12, 'phshiftunit', 'binary', 'Phase-shift unit', 'Phase-shift unit(radians|degrees){degrees}', '(radians|degrees){degrees}', .false., 'degrees')
        ! alternative inputs
        call import_particles%set_input('alt_ios', 1, 'stktab', 'file', 'List of per-micrograph particle stacks',&
        &'List of per-micrograph particle image stacks to import', 'per-micrograph stack list; e.g. stktab.txt', .false., '')
        call import_particles%set_input('alt_ios', 2, 'stk', 'file', 'Stack of particles',&
        &'Stack of particle images to import', 'e.g. stk.mrcs', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_import_particles

    subroutine new_make_cavgs
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &                           ! name
        &'Make class averages',&               ! descr_short
        &'is a distributed workflow for generating class averages or initial random references&
        &for cluster2D execution',&                ! descr_long
        &'simple_distr_exec',&                     ! executable
        &1, 5, 0, 0, 1, 0, 2, .true.)              ! # entries in each group, requires sp_project
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
        call make_cavgs%set_input('parm_ios', 5, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call make_cavgs%set_input('filt_ctrls', 1, shellw)
        ! mask controls
        ! <empty>
        ! computer controls
        call make_cavgs%set_input('comp_ctrls', 1, nparts)
        call make_cavgs%set_input('comp_ctrls', 2, nthr)
    end subroutine new_make_cavgs

    subroutine new_make_oris
        ! PROGRAM SPECIFICATION
        call make_oris%new(&
        &'make_oris',&                       ! name
        &'Make orientations',&           ! descr_short
        &'is a program for making SIMPLE orientation files. Make_oris generates random Euler angles e1.in.[0,360], e2.in.[0,180] '&
        &'and e3.in.[0,360] and random origin shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete is set to an integer number '&
        &'> 0, the shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete is set to an integer number > 0, the orientations '&
        &'produced are randomly sampled from the set of ndiscrete quasi-even projection directions, and the in-plane parameters '&
        &'are assigned randomly. If even=yes, then all nptcls orientations are assigned quasi-even projection directions and '&
        &'random in-plane parameters. If nstates is set to some integer number > 0, then states are assigned randomly .in.[1,nstates]',&
        &'simple_exec',&                     ! executable
        &0, 10, 0, 0, 0, 0, 1, .false.)      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call make_oris%set_input('parm_ios', 1,  'nptcls', 'num', 'Number of per-particle orientations', 'Number of per-particle orientations to produce', '# per-ptcl oris', .true., 1.0)
        call make_oris%set_input('parm_ios', 2,  'ncls', 'num', 'Number of random class labels', 'Number of random class labels to produce', '# classes', .false., 0.)
        call make_oris%set_input('parm_ios', 3,  outfile)
        call make_oris%set_input('parm_ios', 4,  'nstates', 'num', 'Number of random state labels', 'Number of random state labels to produce', '# states', .false., 0.0)
        call make_oris%set_input('parm_ios', 5,  pgrp)
        make_oris%parm_ios(5)%required = .false.
        call make_oris%set_input('parm_ios', 6,  sherr)
        call make_oris%set_input('parm_ios', 7,  angerr)
        call make_oris%set_input('parm_ios', 8,  'even', 'binary', 'Generate even projections', 'Generate quasi-even projection directions(yes|no){no}', '(yes|no){no}', .false., 'no')
        call make_oris%set_input('parm_ios', 9,  'ndiscrete', 'num', 'Number of discrete projection directions', 'Number of discrete projection directions to sample from', '# discrete projs', .false., 0.)
        call make_oris%set_input('parm_ios', 10, oritype)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call make_oris%set_input('comp_ctrls', 1, nthr)
    end subroutine new_make_oris

    subroutine new_make_pickrefs
        ! PROGRAM SPECIFICATION
        call make_pickrefs%new(&
        &'make_pickrefs', &                            ! name
        &'Generate picking references',&               ! descr_short
        &'is a program for generating references for template-based particle picking',& ! descr_long
        &'simple_exec',&                               ! executable
        &0, 3, 2, 1, 0, 0, 1, .false.)                 ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call make_pickrefs%set_input('parm_ios', 1, smpd)
        call make_pickrefs%set_input('parm_ios', 2, pcontrast)
        call make_pickrefs%set_input('parm_ios', 3, mkdir_)
        ! alternative inputs
        call make_pickrefs%set_input('alt_ios', 1, 'stk', 'file', 'Stack of 2D picking references', 'Stack of 2D picking references', 'e.g. refs.mrc', .false., '')
        call make_pickrefs%set_input('alt_ios', 2, 'vol1', 'file', 'Volume', 'Volume to re-project', 'vol.mrc file', .false., '')
        ! search controls
        call make_pickrefs%set_input('srch_ctrls', 1, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call make_pickrefs%set_input('comp_ctrls', 1, nthr)
    end subroutine new_make_pickrefs

    subroutine new_mask
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',& ! name
        &'Mask images/volumes',&                                         ! descr_short
        &'is a program for masking of 2D images and volumes. If you want to mask your images with a spherical mask with a soft &
        & falloff, set msk to the radius in pixels',&                         ! descr_long
        &'simple_exec',&                                                      ! executable
        &0, 3, 2, 1, 1,10, 1, .false.)                                        ! # entries in each group, requires sp_project
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
        & 'input volume e.g. vol.mrc', .false., '')
        ! search controls
        call mask%set_input('srch_ctrls', 1, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
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
        call mask%set_input('mask_ctrls', 6, mw)
        call mask%set_input('mask_ctrls', 7, 'width', 'num', 'Inner mask falloff',&
        &'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false., 10.)
        call mask%set_input('mask_ctrls', 8, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6.)
        call mask%set_input('mask_ctrls', 9, 'taper_edges', 'binary', 'Taper edges',&
        &'Whether to taper the edges of image/volume(yes|no){no}', '(yes|no){no}', .false., 'no')
        call mask%set_input('mask_ctrls',10, 'pdbfile', 'file', 'PDB for 3D envelope masking',&
        &'PDB file used to determine the mask', 'e.g. molecule.pdb', .false., '')
        ! computer controls
        call mask%set_input('comp_ctrls', 1, nthr)
    end subroutine new_mask

    subroutine new_motion_correct
        ! PROGRAM SPECIFICATION
        call motion_correct%new(&
        &'motion_correct', &                                                                  ! name
        &'Motion correction of movies',&                                                      ! descr_short
        &'is a distributed workflow for motion correction of movies based on the same&
        & principal strategy as Grigorieffs program. There are two important&
        & differences: automatic weighting of the frames using a correlation-based M-estimator and&
        & continuous optimisation of the shift parameters. If&
        & dose_rate and exp_time are given the individual frames will be low-pass filtered accordingly&
        & (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to&
        & the down-scaling factor (for super-resolution movies). If nframesgrp is given the frames will&
        & be pre-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given, a&
        & contiguous subset of frames will be averaged without any dose-weighting applied',&   ! descr_long
        &'simple_distr_exec',&                                                                 ! executable
        &0, 5, 0, 6, 2, 0, 2, .true.)                                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call motion_correct%set_input('parm_ios', 1, 'dose_rate', 'num', 'Dose rate', 'Dose rate in e/Ang^2/sec', 'in e/Ang^2/sec', .false., 6.)
        call motion_correct%set_input('parm_ios', 2, 'exp_time', 'num', 'Exposure time', 'Exposure time in seconds', 'in seconds', .false., 10.)
        call motion_correct%set_input('parm_ios', 3, 'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.)
        call motion_correct%set_input('parm_ios', 4, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., '')
        call motion_correct%set_input('parm_ios', 5, pspecsz)
        ! alternative inputs
        ! <empty>
        ! search controls
        call motion_correct%set_input('srch_ctrls', 1, trs)
        call motion_correct%set_input('srch_ctrls', 2, startit)
        call motion_correct%set_input('srch_ctrls', 3, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call motion_correct%set_input('srch_ctrls', 4, fromf)
        call motion_correct%set_input('srch_ctrls', 5, tof)
        call motion_correct%set_input('srch_ctrls', 6, nsig)
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

    subroutine new_motion_correct_tomo
        ! PROGRAM SPECIFICATION
        call motion_correct_tomo%new(&
        &'motion_correct_tomo', &                   ! name
        &'Motion correction of tomography movies',& ! descr_short
        &'is a distributed workflow for motion correction of tomography movies based on the same &
        &principal strategy as Grigorieffs program. There are two important differences: automatic &
        &weighting of the frames using a correlation-based M-estimator and continuous optimisation of &
        &the shift parameters. If dose_rate and exp_time are given, the individual frames will be &
        &low-pass filtered accordingly (dose-weighting strategy). The exp_doc document should contain &
        &per line exp_time=X and dose_rate=Y. It is asssumed that the input list of movies (one per tilt) &
        &are ordered temporally. This is necessary for correct dose-weighting of tomographic tilt series. &
        &If scale is given, the movie will be Fourier cropped according to the down-scaling factor &
        &(for super-resolution movies). If nframesgrp is given the frames will be pre-averaged in the given &
        &chunk size (Falcon 3 movies)',& ! descr_long
        &'simple_distr_exec',&           ! executable
        &0, 8, 0, 4, 3, 0, 1, .false.)   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call motion_correct_tomo%set_input('parm_ios', 1, 'tomoseries', 'file', '.txt filetable of filetables of tomograms',&
        &'.txt filetable of filetables; each line referring to a .txt file listing all movies in the tilt-series', 'e.g. filetab_of_filetabs.txt', .true., '')
        call motion_correct_tomo%set_input('parm_ios', 2, 'exp_doc', 'file', '.txt file with exp_time and dose_rate per tomogram',&
        &'.txt file with exp_time and dose_rate per tomogram', 'e.g. exp_doc.txt', .true., '')
        call motion_correct_tomo%set_input('parm_ios', 3, smpd)
        call motion_correct_tomo%set_input('parm_ios', 4, 'dir', 'dir', 'Output directory', 'Output directory', 'e.g. motion_correct/', .false., 'motion_correct')
        call motion_correct_tomo%set_input('parm_ios', 5, 'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.)
        call motion_correct_tomo%set_input('parm_ios', 6, pspecsz)
        call motion_correct_tomo%set_input('parm_ios', 7, numlen)
        call motion_correct_tomo%set_input('parm_ios', 8, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call motion_correct_tomo%set_input('srch_ctrls', 1, trs)
        call motion_correct_tomo%set_input('srch_ctrls', 2, startit)
        call motion_correct_tomo%set_input('srch_ctrls', 3, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call motion_correct_tomo%set_input('srch_ctrls', 4, nsig)
        ! filter controls
        call motion_correct_tomo%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms)', 'in Angstroms', .false., 15.)
        call motion_correct_tomo%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms)', 'in Angstroms', .false., 8.)
        call motion_correct_tomo%set_input('filt_ctrls', 3, kv)
        ! mask controls
        ! <empty>
        ! computer controls
        call motion_correct_tomo%set_input('comp_ctrls', 1, nthr)
    end subroutine new_motion_correct_tomo

    subroutine new_new_project
        ! PROGRAM SPECIFICATION
        call new_project%new(&
        &'new_project',&                     ! name
        &'Create a new project',&            ! descr_short
        &'is a program for creating a new project. SIMPLE3.0 relies on a monolithic project file for controlling &
        &execution on distributed and shared-memory systems and for unified meta-data management. This program &
        &creates a directory named projname and a file projname.simple inside that directory that contains all &
        &information about the project as well as all meta data generated by the different SIMPLE programs. This &
        &file is mirrored by an abstract data type in the back-end, which manages the parameters and &
        &meta-data I/O required for execution of SIMPLE',& ! descr_longg
        &'simple_exec',&                     ! executable
        &0, 2, 0, 0, 0, 0, 7, .false.)       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call new_project%set_input('parm_ios', 1, projname)
        call new_project%set_input('parm_ios', 2, user_email)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call new_project%set_input('comp_ctrls', 1, time_per_image)
        call new_project%set_input('comp_ctrls', 2, user_account)
        call new_project%set_input('comp_ctrls', 3, user_project)
        call new_project%set_input('comp_ctrls', 4, qsys_partition)
        call new_project%set_input('comp_ctrls', 5, qsys_qos)
        call new_project%set_input('comp_ctrls', 6, qsys_reservation)
        call new_project%set_input('comp_ctrls', 7, job_memory_per_task)
    end subroutine new_new_project

    subroutine new_pick
        ! PROGRAM SPECIFICATION
        call pick%new(&
        &'pick', & ! name
        &'Template-based particle picking',&                               ! descr_short
        &'is a distributed workflow for template-based particle picking',& ! descr_long
        &'simple_distr_exec',&                                             ! executable
        &0, 2, 0, 2, 1, 0, 1, .true.)                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call pick%set_input('parm_ios', 1, 'refs', 'file', 'Picking 2D references',&
        &'2D references used for automated picking', 'e.g. pickrefs.mrc file with references', .true., '')
        call pick%set_input('parm_ios', 2, 'dir', 'dir', 'Output directory', 'Output directory', 'e.g. pick/', .false., 'pick')
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

    subroutine new_postprocess
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',& ! name
        &'Post-processing of volume',&                                        ! descr_short
        &'is a program for map post-processing. Use program volops to estimate the B-factor with the Guinier plot',& ! descr_long
        &'simple_exec',&                                                            ! executable
        &1, 4, 0, 0, 7, 9, 1, .false.)                                              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to post-process', &
        & 'input volume e.g. vol.mrc', .false., '')
        ! parameter input/output
        call postprocess%set_input('parm_ios', 1, smpd)
        call postprocess%set_input('parm_ios', 2, mkdir_)
        call postprocess%set_input('parm_ios', 3, projfile)
        postprocess%parm_ios(3)%required = .false.
        call postprocess%set_input('parm_ios', 4, 'state', 'num', 'State to postprocess', 'State to postprocess{1}', 'Input state{1}', .false., 1.0)
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
        call postprocess%set_input('filt_ctrls', 5, 'fsc', 'file', 'FSC file', 'Binary file with FSC info&
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

    subroutine new_powerspecs
        ! PROGRAM SPECIFICATION
        call powerspecs%new(&
        &'powerspecs', &                                                              ! name
        &'Calculate powerspectra from micrographs',&                                  ! descr_short
        &'is a program for generating powerspectra from a filetable of micrographs',& ! descr_long
        &'simple_distr_exec',&                                                        ! executable
        &0, 6, 2, 0, 1, 0, 2, .false.)                                                ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call powerspecs%set_input('parm_ios', 1, smpd)
        call powerspecs%set_input('parm_ios', 2, 'fbody', 'string', 'Template output powerspec name',&
        &'Template output powerspec name', 'e.g. spec_', .false., '')
        call powerspecs%set_input('parm_ios', 3, pspecsz)
        call powerspecs%set_input('parm_ios', 4, 'speckind', 'multi', 'Power spectrum kind', 'Power spectrum kind(real|power|sqrt|log|phase){sqrt}', '(real|power|sqrt|log|phase){sqrt}', .false., 'sqrt')
        call powerspecs%set_input('parm_ios', 5, clip)
        call powerspecs%set_input('parm_ios', 6, mkdir_)
        ! alternative inputs
        call powerspecs%set_input('alt_ios', 1, stk)
        call powerspecs%set_input('alt_ios', 2, 'filetab', 'file', 'List of micrograph files', 'List of micrograph files to generate powerspectra from', 'e.g. filetab.txt', .true., '')
        ! search controls
        ! <empty>
        ! filter controls
        call powerspecs%set_input('filt_ctrls', 1, 'lp', 'num', 'Low-pass limit of resolution mask', 'Low-pass resolution limit of resolution mask', 'low-pass limit in Angstroms', .false., 6.)
        ! mask controls
        ! <empty>
        ! computer controls
        call powerspecs%set_input('comp_ctrls', 1, nparts)
        call powerspecs%set_input('comp_ctrls', 2, nthr)
    end subroutine new_powerspecs

    subroutine new_preprocess
        ! PROGRAM SPECIFICATION
        call preprocess%new(&
        &'preprocess', & ! name
        &'Preprocessing',&                                                                  ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in sequence',&
        &'simple_distr_exec',&                                                              ! executable
        &0, 9, 0, 13, 5, 0, 2, .true.)                                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call preprocess%set_input('parm_ios', 1,  'dose_rate', 'num', 'Dose rate', 'Dose rate in e/Ang^2/sec', 'in e/Ang^2/sec', .false., 6.0)
        call preprocess%set_input('parm_ios', 2,  'exp_time', 'num', 'Exposure time', 'Exposure time in seconds', 'in seconds', .false., 10.)
        call preprocess%set_input('parm_ios', 3,  'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.0)
        call preprocess%set_input('parm_ios', 4,  pcontrast)
        call preprocess%set_input('parm_ios', 5,  'box_extract', 'num', 'Box size on extraction', 'Box size on extraction in pixels', 'in pixels', .false., 0.)
        call preprocess%set_input('parm_ios', 6,  'refs', 'file', 'Picking 2D references',&
        &'2D references used for automated picking', 'e.g. pickrefs.mrc file with references', .false., 'pickrefs.mrc')
        call preprocess%set_input('parm_ios', 7,  'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_')
        call preprocess%set_input('parm_ios', 8,  pspecsz)
        call preprocess%set_input('parm_ios', 9,  numlen)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess%set_input('srch_ctrls', 1, trs)
        call preprocess%set_input('srch_ctrls', 2, startit)
        call preprocess%set_input('srch_ctrls', 3, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call preprocess%set_input('srch_ctrls', 4, fromf)
        call preprocess%set_input('srch_ctrls', 5, tof)
        call preprocess%set_input('srch_ctrls', 6, nsig)
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
        call preprocess%set_input('filt_ctrls', 4, 'hp_ctf_estimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30.)
        call preprocess%set_input('filt_ctrls', 5, 'lp_pick', 'num', 'Low-pass limit for picking',&
        & 'Low-pass limit for picking in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        ! mask controls
        ! <empty>
        ! computer controls
        call preprocess%set_input('comp_ctrls', 1, nparts)
        call preprocess%set_input('comp_ctrls', 2, nthr)
    end subroutine new_preprocess

    subroutine new_preprocess_stream
        ! PROGRAM SPECIFICATION
        call preprocess_stream%new(&
        &'preprocess_stream', & ! name
        &'Preprocessing in streaming mode',&                                                ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_distr_exec',&                                                              ! executable
        &1,13, 0, 13, 5, 0, 2, .true.)                                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess_stream%set_input('img_ios', 1, 'dir_movies', 'dir', 'Input movies directory', 'Where the movies ot process will squentially appear', 'e.g. data/', .true., 'preprocess/')
        ! parameter input/output
        call preprocess_stream%set_input('parm_ios', 1, 'dose_rate', 'num', 'Dose rate', 'Dose rate in e/Ang^2/sec', 'in e/Ang^2/sec', .false., 6.0)
        call preprocess_stream%set_input('parm_ios', 2, 'exp_time', 'num', 'Exposure time', 'Exposure time in seconds', 'in seconds', .false., 10.)
        call preprocess_stream%set_input('parm_ios', 3, 'scale', 'num', 'Down-scaling factor', 'Down-scaling factor to apply to the movies', '(0-1)', .false., 1.0)
        call preprocess_stream%set_input('parm_ios', 4, pcontrast)
        call preprocess_stream%set_input('parm_ios', 5, 'box_extract', 'num', 'Box size on extraction', 'Box size on extraction in pixels', 'in pixels', .false., 0.)
        call preprocess_stream%set_input('parm_ios', 6, 'refs', 'file', 'Picking 2D references',&
        &'2D references used for automated picking', 'e.g. pickrefs.mrc file with references', .false., 'pickrefs.mrc')
        call preprocess_stream%set_input('parm_ios', 7, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_')
        call preprocess_stream%set_input('parm_ios', 8, pspecsz)
        call preprocess_stream%set_input('parm_ios', 9, kv)
        preprocess_stream%parm_ios(9)%required = .true.
        call preprocess_stream%set_input('parm_ios', 10, cs)
        preprocess_stream%parm_ios(10)%required = .true.
        call preprocess_stream%set_input('parm_ios', 11, fraca)
        preprocess_stream%parm_ios(11)%required = .true.
        call preprocess_stream%set_input('parm_ios', 12, smpd)
        preprocess_stream%parm_ios(12)%required = .true.
        call preprocess_stream%set_input('parm_ios', 13, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess_stream%set_input('srch_ctrls', 1, trs)
        call preprocess_stream%set_input('srch_ctrls', 2, startit)
        call preprocess_stream%set_input('srch_ctrls', 3, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3)', '{0}', .false., 0.)
        call preprocess_stream%set_input('srch_ctrls', 4, fromf)
        call preprocess_stream%set_input('srch_ctrls', 5, tof)
        call preprocess_stream%set_input('srch_ctrls', 6, nsig)
        call preprocess_stream%set_input('srch_ctrls', 7, dfmin)
        call preprocess_stream%set_input('srch_ctrls', 8, dfmax)
        call preprocess_stream%set_input('srch_ctrls', 9, dfstep)
        call preprocess_stream%set_input('srch_ctrls',10, astigtol)
        call preprocess_stream%set_input('srch_ctrls',11, 'thres', 'num', 'Picking distance threshold','Picking distance filer (in pixels)', 'in pixels', .false., 0.)
        call preprocess_stream%set_input('srch_ctrls',12, 'rm_outliers', 'binary', 'Remove micrograph image outliers for picking',&
        & 'Remove micrograph image outliers for picking(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call preprocess_stream%set_input('srch_ctrls',13, 'ndev', 'num', '# of sigmas for picking clustering', '# of standard deviations threshold for picking one cluster clustering{2}', '{2}', .false., 2.)
        ! filter controls
        call preprocess_stream%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){15}', 'in Angstroms{15}', .false., 15.)
        call preprocess_stream%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess_stream%set_input('filt_ctrls', 3, 'lp_ctf_estimate', 'num', 'Low-pass limit for CTF parameter estimation',&
        & 'Low-pass limit for CTF parameter estimation in Angstroms{5}', 'in Angstroms{5}', .false., 5.)
        call preprocess_stream%set_input('filt_ctrls', 4, 'hp_ctf_estimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30.)
        call preprocess_stream%set_input('filt_ctrls', 5, 'lp_pick', 'num', 'Low-pass limit for picking',&
        & 'Low-pass limit for picking in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        ! mask controls
        ! <empty>
        ! computer controls
        call preprocess_stream%set_input('comp_ctrls', 1, nparts)
        call preprocess_stream%set_input('comp_ctrls', 2, nthr)
    end subroutine new_preprocess_stream

    subroutine new_print_fsc
        ! PROGRAM SPECIFICATION
        call print_fsc%new(&
        &'print_fsc', &                                                          ! name
        &'Print FSC file produced by REFINE3D',&                                 ! descr_short
        &'is a program for printing the binary FSC files produced by REFINE3D',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &0, 3, 0, 0, 0, 0, 0, .false.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_fsc%set_input('parm_ios', 1, smpd)
        call print_fsc%set_input('parm_ios', 2, box)
        call print_fsc%set_input('parm_ios', 3, 'fsc', 'file', 'FSC file', 'Binary file with FSC info',&
        'input binary file e.g. fsc_state01.bin', .true., 'fsc_state01.bin')
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
    end subroutine new_print_fsc

    subroutine new_print_magic_boxes
        ! PROGRAM SPECIFICATION
        call print_magic_boxes%new(&
        &'print_magic_boxes', &                                   ! name
        &'Print magic boxes (fast FFT)',&                         ! descr_short
        &'is a program for printing magic box sizes (fast FFT)',& ! descr_long
        &'simple_exec',&                                          ! executable
        &0, 3, 0, 0, 0, 0, 0, .false.)                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_magic_boxes%set_input('parm_ios', 1, smpd)
        call print_magic_boxes%set_input('parm_ios', 2, box)
        call print_magic_boxes%set_input('parm_ios', 3, 'moldiam', 'num', 'Molecular diameter', 'Molecular diameter(in pixels)',&
        'give # pixels of diameter', .false., 140.)
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
    end subroutine new_print_magic_boxes

    subroutine new_print_project_info
        ! PROGRAM SPECIFICATION
        call print_project_info%new(&
        &'print_project_info', &                                             ! name
        &'Print project info',&                                              ! descr_short
        &'is a program prints information abourt a *.simple project file',&  ! descr_long
        &'simple_exec',&                                                     ! executable
        &0, 0, 0, 0, 0, 0, 0, .true.)                                        ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
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
    end subroutine new_print_project_info

    subroutine new_reproject
        ! PROGRAM SPECIFICATION
        call reproject%new(&
        &'reproject',&                           ! name
        &'Re-project volume',&                    ! descr_short
        &'is a program for re-projecting a volume using Fourier interpolation. Input is a SPIDER or &
        &MRC volume. Output is a stack of projection images of the same format as the inputted volume. Projections &
        &are generated by extracting central sections from the Fourier volume and back transforming the 2D FTs. &
        &nspace controls the number of projection directions. The  oritab parameter allows you to input the orientations &
        &that you wish to have your volume projected in. pgrp controls the point-group symmetry &
        &c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahedral). The point-group symmetry is &
        &used to restrict the set of projections to within the asymmetric unit. &
        &neg inverts the contrast of the projections',&
        &'simple_exec',&                       ! executable
        &2, 3, 0, 2, 0, 1, 1, .false.)         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call reproject%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume for creating 2D central &
        & sections', 'input volume e.g. vol.mrc', .true., 'vol1.mrc')
        call reproject%set_input('img_ios', 2,  outstk)
        ! parameter input/output
        call reproject%set_input('parm_ios', 1, smpd)
        call reproject%set_input('parm_ios', 2, oritab)
        call reproject%set_input('parm_ios', 3, 'neg', 'binary', 'Invert contrast','Invert contrast of projections(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        call reproject%set_input('srch_ctrls', 1, nspace)
        call reproject%set_input('srch_ctrls', 2, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call reproject%set_input('mask_ctrls', 1, msk)
        reproject%mask_ctrls(1)%required = .false.
        ! computer controls
        call reproject%set_input('comp_ctrls', 1, nthr)
    end subroutine new_reproject

    subroutine new_normalize
        ! PROGRAM SPECIFICATION
        call normalize_%new(&
        &'normalize',&                         ! name
        &'Normalize volume/stack',&            ! descr_short
        &'is a program for normalization of MRC or SPIDER stacks and volumes',&
        &'simple_exec',&                       ! executable
        &0, 4, 2, 0, 0, 1, 1, .false.)         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call normalize_%set_input('parm_ios', 1, smpd)
        call normalize_%set_input('parm_ios', 2, 'norm',       'binary', 'Normalize',       'Statistical normalization: avg=zero, var=1(yes|no){no}',     '(yes|no){no}', .false., 'no')
        call normalize_%set_input('parm_ios', 3, 'noise_norm', 'binary', 'Noise normalize', 'Statistical normalization based on background(yes|no){no}',  '(yes|no){no}', .false., 'no')
        call normalize_%set_input('parm_ios', 4, 'shell_norm', 'binary', 'Power whitening', 'Normalisation of each Fourier shell to power=1(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        call normalize_%set_input('alt_ios', 1, 'stk',  'file', 'Stack to normalize',  'Stack of images to normalize', 'e.g. imgs.mrc', .false., '')
        call normalize_%set_input('alt_ios', 2, 'vol1', 'file', 'Volume to normalize', 'Volume to normalize',          'e.g. vol.mrc',  .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call normalize_%set_input('mask_ctrls', 1, msk)
        ! computer controls
        call normalize_%set_input('comp_ctrls', 1, nthr)
    end subroutine new_normalize

    subroutine new_orisops
        ! PROGRAM SPECIFICATION
        call orisops%new(&
        &'orisops',&                      ! name
        &'Standard orientation editing',& ! descr_short
        &'is a program for modifying SIMPLE orientation/parameter files. If errify=yes, '&
        &'uniform random angular errors .in.[-angerr,angerr], and uniform origin shift errors .in.[-sherr,sherr], '&
        &'and uniform random defocus errors .in.[-dferr,dferr] are introduced. If nstates > 1 then random states are assigned '&
        &'.in.[1,nstates]. If mirr=2d, then the Euler angles in oritab are mirrored according to the relation '&
        &'e1=e1, e2=180.+e2, e3=-e3. If mirr=3d, then the Euler angles in oritab are mirrored according to the '&
        &'relation R=M(M*R), where R is the rotation matrix calculated from the Euler angle triplet and M is a '&
        &'3D reflection matrix. If e1, e2, or e3 is '&
        &'inputted, the orientations are rotated correspondingly. If you input state, '&
        &'only the orientations assigned to state state are rotated. If mul is defined, the origin shifts are multiplied with mul. '&
        &'If zero=yes, then the shifts are zeroed',&
        &'simple_exec',&                       ! executable
        &0, 19, 0, 0, 0, 0, 0, .false.)        ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call orisops%set_input('parm_ios', 1,  oritab)
        orisops%parm_ios(1)%required = .true.
        call orisops%set_input('parm_ios', 2,  outfile)
        call orisops%set_input('parm_ios', 3,  e1)
        call orisops%set_input('parm_ios', 4,  e2)
        call orisops%set_input('parm_ios', 5,  e3)
        call orisops%set_input('parm_ios', 6,  'nstates', 'num', 'Number of random state labels', 'Number of random state labels to insert', '# states', .false., 0.0)
        call orisops%set_input('parm_ios', 7,  pgrp)
        orisops%parm_ios(7)%required = .false.
        call orisops%set_input('parm_ios', 8,  ctf)
        orisops%parm_ios(8)%required = .false.
        call orisops%set_input('parm_ios', 9,  angerr)
        call orisops%set_input('parm_ios', 10, sherr)
        call orisops%set_input('parm_ios', 11, dferr)
        call orisops%set_input('parm_ios', 12, 'zero', 'binary', 'Zero shifts', 'Zero shifts(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 13, 'ndiscrete', 'num', 'Number of discrete projection directions',&
        &'Number of projection directions to use for discretization of input orientations', '# discrete projs', .false., 0.)
        call orisops%set_input('parm_ios', 14, 'state', 'num', 'State to modify', 'Index of state to modify', 'give state index', .false., 1.)
        call orisops%set_input('parm_ios', 15, 'errify', 'binary', 'Errify input orientations/CTF parameters',&
        &'Introduce uniform random errors in input orientations/CTF parameters(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 16, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call orisops%set_input('parm_ios', 17, 'mirr', 'multi', 'Mirror orientations', 'Mirror orientations(2d|3d|no){no}', '(2d|3d|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 18, 'symrnd', 'binary', 'Randomize over subgroubs of point-group', 'Expand orientations over entire unit sphere by &
        &permutation according to randomly selected subgroup symmetry(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 19, oritype)
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
    end subroutine new_orisops

    subroutine new_oristats
        ! PROGRAM SPECIFICATION
        call oristats%new(&
        &'oristats',&                             ! name
        &'Statistical analyses of orientations',& ! descr_short
        &'is a program for analyzing SIMPLE orientation/parameter files. If two orientation '&
        &'tables (oritab and oritab2) are inputted, statistics of the distances between the orientations '&
        &'in the two documents are provided',&
        &'simple_exec',&                          ! executable
        &0, 11, 0, 0, 0, 0, 1, .false.)           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call oristats%set_input('parm_ios', 1,  oritab)
        oristats%parm_ios(1)%required = .true.
        call oristats%set_input('parm_ios', 2,  oritab2)
        call oristats%set_input('parm_ios', 3,  pgrp)
        oristats%parm_ios(3)%required = .false.
        call oristats%set_input('parm_ios', 4,  nspace)
        call oristats%set_input('parm_ios', 5,  oritype)
        call oristats%set_input('parm_ios', 6,  weights2D)
        call oristats%set_input('parm_ios', 7,  weights3D)
        call oristats%set_input('parm_ios', 8,  'ctfstats', 'binary', 'CTF statistics',         'Provide statistics about CTF(yes|no){no}',                      '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 9,  'classtats', 'binary', 'Class statistics',      'Provide statistics about 2D clusters(yes|no){no}',              '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 10, 'projstats', 'binary', 'Projection statistics', 'Provide statistics about projection directions(yes|no){no}',    '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 11, 'trsstats', 'binary', 'Shift statistics',       'Provide statistics about rotational origin shifts(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call oristats%set_input('comp_ctrls', 1, nthr)
    end subroutine new_oristats

    subroutine new_reconstruct3D
        ! PROGRAM SPECIFICATION
        call reconstruct3D%new(&
        &'reconstruct3D',&                                                     ! name
        &'3D reconstruction from oriented particles',&                         ! descr_long
        &'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks, '&
        &'given input orientations and state assignments. The algorithm is based on direct Fourier inversion '&
        &'with a Kaiser-Bessel (KB) interpolation kernel',&
        &'simple_distr_exec',&                                                 ! executable
        &0, 1, 0, 2, 2, 2, 2, .true.)                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call reconstruct3D%set_input('parm_ios', 1, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call reconstruct3D%set_input('srch_ctrls', 1, pgrp)
        call reconstruct3D%set_input('srch_ctrls', 2, frac)
        ! filter controls
        call reconstruct3D%set_input('filt_ctrls', 1, eo)
        call reconstruct3D%set_input('filt_ctrls', 2, shellw)
        ! mask controls
        call reconstruct3D%set_input('mask_ctrls', 1, msk)
        call reconstruct3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call reconstruct3D%set_input('comp_ctrls', 1, nparts)
        call reconstruct3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_reconstruct3D

    subroutine new_refine3D
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',&                                                                               ! name
        &'3D refinement',&                                                                          ! descr_short
        &'is a distributed workflow for 3D refinement based on probabilistic projection matching',& ! descr_long
        &'simple_distr_exec',&                                                                      ! executable
        &1, 1, 0, 13, 8, 5, 2, .true.)                                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%set_input('img_ios', 1, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. vol.mrc', .false., 'vol1.mrc')
        ! parameter input/output
        call refine3D%set_input('parm_ios', 1, mkdir_)
        ! alternative inputs
        ! <empty>
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
        &clustersym){no}', '(snhc|single|multi|greedy_single|cont_single|greedy_multi|cluster|clustersym){single}', .false., 'single')
        call refine3D%set_input('srch_ctrls', 13, neigh)
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
        &filtering and resolution estimation(yes|no){no}', '(yes|no){no}', .false., 'no')
        call refine3D%set_input('filt_ctrls', 7, 'weights3D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 'no')
        call refine3D%set_input('filt_ctrls', 8, shellw)
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

    subroutine new_refine3D_init
        ! PROGRAM SPECIFICATION
        call refine3D_init%new(&
        &'refine3D_init',& ! name
        &'Random initialisation of 3D refinement',&                                                            ! descr_short
        &'is a distributed workflow for generating a random initial 3D model for initialisation of refine3D',& ! descr_long
        &'simple_distr_exec',&                                                                                 ! executable
        &0, 1, 0, 3, 1, 2, 2, .true.)                                                                          ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call refine3D_init%set_input('parm_ios', 1, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D_init%set_input('srch_ctrls', 1, nspace)
        call refine3D_init%set_input('srch_ctrls', 2, pgrp)
        call refine3D_init%set_input('srch_ctrls', 3, 'nran', 'num', 'Number of random samples', 'Number of images to randomly sample for 3D reconstruction',&
        &'# random samples', .false., 0.)
        ! filter controls
        call refine3D_init%set_input('filt_ctrls', 1, shellw)
        ! mask controls
        call refine3D_init%set_input('mask_ctrls', 1, msk)
        call refine3D_init%set_input('mask_ctrls', 2, inner)
        ! computer controls
        call refine3D_init%set_input('comp_ctrls', 1, nparts)
        call refine3D_init%set_input('comp_ctrls', 2, nthr)
    end subroutine new_refine3D_init

    subroutine new_scale
        ! PROGRAM SPECIFICATION
        call scale%new(&
        &'scale', & ! name
        &'Re-scaling MRC & SPIDER stacks and volumes',&                                         ! descr_short
        &'is a program for re-scaling, clipping and padding MRC & SPIDER stacks and volumes',&  ! descr_long
        &'simple_exec',&                                                                        ! executable
        &0, 8, 3, 0, 0, 0, 1, .false.)                                                          ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call scale%set_input('parm_ios', 1, smpd)
        call scale%set_input('parm_ios', 2, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'new box in pixels', .false., 0.)
        call scale%set_input('parm_ios', 3, 'scale', 'num', 'Scaling ratio', 'Target box ratio for scaling(0-1)', '(0-1)', .false., 1.)
        call scale%set_input('parm_ios', 4, 'scale2', 'num', 'Second scaling ratio', 'Second target box ratio for scaling(0-1)', '(0-1)', .false., 1.)
        call scale%set_input('parm_ios', 5, clip)
        call scale%set_input('parm_ios', 6, outvol)
        call scale%set_input('parm_ios', 7, outstk)
        call scale%set_input('parm_ios', 8, 'outstk2', 'file', 'Second output stack name', 'Second output images stack name', 'e.g. outstk2.mrc', .false., '')
        ! alternative inputs
        call scale%set_input('alt_ios', 1, stk)
        scale%alt_ios(1)%required = .false.
        call scale%set_input('alt_ios', 2, 'vol1', 'file', 'Input volume', 'Input volume to re-scale',&
        &'input volume e.g. vol.mrc', .false., '')
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

    subroutine new_scale_project
        ! PROGRAM SPECIFICATION
        call scale_project%new(&
        &'scale', & ! name
        &'Re-scaling of MRC & SPIDER stacks',&                                             ! descr_short
        &'is a program for re-scaling MRC & SPIDER stacks part of project specification',& ! descr_long
        &'simple_exec',&                                                                   ! executable
        &0, 1, 0, 0, 0, 0, 2, .false.)                                                     ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call scale_project%set_input('parm_ios', 1, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'new box in pixels', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call scale_project%set_input('comp_ctrls', 1, nparts)
        call scale_project%set_input('comp_ctrls', 2, nthr)
    end subroutine new_scale_project

    subroutine new_select_
        ! PROGRAM SPECIFICATION
        call select_%new(&
        &'select',&                                         ! name
        &'Select images',&                                  ! descr_short
        &'is a program for selecting files based on image correlation matching',& ! descr_long
        &'simple_exec',&                                    ! executable
        &8, 0, 0, 0, 0, 0, 1, .false.)                      ! # entries in each group
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call select_%set_input('img_ios', 1, stk )
        call select_%set_input('img_ios', 2, 'stk2',    'file', 'Stack of selected images', 'Stack of selected images', 'e.g. selected(cavgs).mrc', .true., '')
        call select_%set_input('img_ios', 3, 'stk3',    'file', 'Stack of images to select from', 'Stack of images to select from', 'e.g. (cavgs)2selectfrom.mrc', .false., '')
        call select_%set_input('img_ios', 4, 'filetab', 'file', 'List of files to select from', 'List of files to select from', 'e.g. filetab.txt', .false., '')
        call select_%set_input('img_ios', 5,  outfile)
        call select_%set_input('img_ios', 6,  outstk)
        call select_%set_input('img_ios', 7,  'dir_select', 'dir', 'Directory for selected images', 'Move selected files to here{selected}', 'select dir', .false., '')
        call select_%set_input('img_ios', 8,  'dir_reject', 'dir', 'Directory for rejected images', 'Move rejected files to here{rejected}', 'reject dir', .false., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call select_%set_input('comp_ctrls', 1, nthr)
    end subroutine new_select_

    subroutine new_shift
        ! PROGRAM SPECIFICATION
        call shift%new(&
        &'shift',&                                                                  ! name
        &'Shift images to rotational origin',&                                      ! descr_short
        &'is a program for shifting a stack according to origin shifts in oritab',& ! descr_long
        &'simple_exec',&                                                            ! executable
        &3, 3, 0, 0, 0, 0, 1, .false.)                                              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call shift%set_input('img_ios', 1, stk)
        shift%img_ios(1)%required = .true.
        call shift%set_input('img_ios', 2, outstk)
        call shift%set_input('img_ios', 3, oritab)
        ! parameter input/output
        call shift%set_input('parm_ios', 1, smpd)
        call shift%set_input('parm_ios', 2, 'mul',     'num',   'Shift multiplication factor', 'Shift multiplication factor{1.0}', 'multiplier', .false., 1.0)
        call shift%set_input('parm_ios', 3, oritype)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call shift%set_input('comp_ctrls', 1, nthr)
    end subroutine new_shift

    subroutine new_simulate_movie
        ! PROGRAM SPECIFICATION
        call simulate_movie%new(&
        &'simulate_movie',&                                 ! name
        &'Simulate DDD movie',&                             ! descr_short
        &'is a program for crude simulation of a DDD movie. Input is a set of projection images to place. &
        &Movie frames are then generated related by randomly shifting the base image and applying noise',& ! descr_long
        &'simple_exec',&                                    ! executable
        &1, 10, 0, 0, 1, 0, 1, .false.)                     ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_movie%set_input('img_ios', 1, stk)
        ! parameter input/output
        call simulate_movie%set_input('parm_ios',  1, smpd)
        call simulate_movie%set_input('parm_ios',  2, 'snr', 'num', 'SNR', 'Signal-to-noise ratio of movie frame', 'signal-to-noise ratio(0.)', .false., 0.)
        call simulate_movie%set_input('parm_ios',  3, kv)
        call simulate_movie%set_input('parm_ios',  4, cs)
        call simulate_movie%set_input('parm_ios',  5, fraca)
        call simulate_movie%set_input('parm_ios',  6, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2.)
        call simulate_movie%set_input('parm_ios',  7, trs)
        call simulate_movie%set_input('parm_ios',  8, 'nframes',  'num', 'Number of frames', 'Number of movie frames', '# frames', .false., 0.)
        call simulate_movie%set_input('parm_ios',  9, 'xdim',  'num', 'x-dimension', 'Number of pixels in x-direction', '# pixels in x', .false., 0.)
        call simulate_movie%set_input('parm_ios', 10, 'ydim',  'num', 'y-dimension', 'Number of pixels in y-direction', '# pixels in y', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call simulate_movie%set_input('filt_ctrls', 1, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_movie%set_input('comp_ctrls', 1, nthr)
    end subroutine new_simulate_movie

    subroutine new_simulate_noise
        ! PROGRAM SPECIFICATION
        call simulate_noise%new(&
        &'simulate_noise',&                                ! name
        &'White noise simulation',&                        ! descr_short
        &'is a program for generating pure noise images',& ! descr_long
        &'simple_exec',&                                   ! executable
        &0, 2, 0, 0, 0, 0, 0, .false.)                     ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call simulate_noise%set_input('parm_ios', 1, box)
        call simulate_noise%set_input('parm_ios', 2, nptcls)
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
    end subroutine new_simulate_noise

    subroutine new_simulate_particles
        ! PROGRAM SPECIFICATION
        call simulate_particles%new(&
        &'simulate_particles',&                                          ! name
        &'Simulate single-particle images',&                             ! descr_short
        &'is a program for simulating single-particle cryo-EM images. It is not a verysophisticated simulator, but &
        &it is nevertheless useful for testing purposes. It does not do any multi-slice simulation and it cannot be &
        &used for simulating molecules containing heavy atoms. It does not even accept a PDB file as an input. Input &
        &is a cryo-EM map, which we usually generate from a PDB file using EMANs program pdb2mrc. The volume is &
        &projected using Fourier interpolation, 20% of the total noise is added to the images (pink noise), they are &
        &then Fourier transformed and multiplied with astigmatic CTF and B-factor. Next, the they are inverse FTed &
        &before the remaining 80% of the noise (white noise) is added',& ! descr_long
        &'simple_exec',&                                                 ! executable
        &1, 15, 0, 1, 2, 1, 1, .false.)                                   ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_particles%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to project', 'input volume e.g. vol.mrc', .false., '')
        ! parameter input/output
        call simulate_particles%set_input('parm_ios', 1,  smpd)
        call simulate_particles%set_input('parm_ios', 2,  nptcls)
        call simulate_particles%set_input('parm_ios', 3,  'snr', 'num', 'SNR', 'Signal-to-noise ratio of particle images', 'signal-to-noise ratio(0.)', .false., 0.)
        call simulate_particles%set_input('parm_ios', 4,  oritab)
        call simulate_particles%set_input('parm_ios', 5,  outfile)
        call simulate_particles%set_input('parm_ios', 6,  outstk)
        call simulate_particles%set_input('parm_ios', 7,  'ndiscrete', 'num', 'Number of discrete projection directions', 'Number of discrete projection directions used in simulation', '# discrete projs', .false., 0.)
        call simulate_particles%set_input('parm_ios', 8,  sherr)
        call simulate_particles%set_input('parm_ios', 9,  kv)
        call simulate_particles%set_input('parm_ios', 10, cs)
        call simulate_particles%set_input('parm_ios', 11, fraca)
        call simulate_particles%set_input('parm_ios', 12, deftab)
        call simulate_particles%set_input('parm_ios', 13, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2.)
        call simulate_particles%set_input('parm_ios', 14, dferr)
        call simulate_particles%set_input('parm_ios', 15, 'astigerr', 'num', 'Astigmatism error', 'Uniform astigmatism error(in microns)', 'error in microns', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call simulate_particles%set_input('srch_ctrls', 1, pgrp)
        ! filter controls
        call simulate_particles%set_input('filt_ctrls', 1, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        call simulate_particles%set_input('filt_ctrls', 2, 'bfacerr', 'num', 'B-factor error', 'Uniform B-factor error(in Angstroms^2)', 'error(in Angstroms^2)', .false., 50.)
        ! mask controls
        call simulate_particles%set_input('mask_ctrls', 1, msk)
        ! computer controls
        call simulate_particles%set_input('comp_ctrls', 1, nthr)
    end subroutine new_simulate_particles

    subroutine new_simulate_subtomogram
        ! PROGRAM SPECIFICATION
        call simulate_subtomogram%new(&
        &'simulate_subtomogram',&                               ! name
        &'Simulate subtomogram',&                               ! descr_short
        &'is a program for crude simulation of a subtomogram',& ! descr_long
        &'simple_exec',&                                        ! executable
        &1, 3, 0, 0, 0,0, 1, .false.)                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_subtomogram%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to use for simulation', 'input volume e.g. vol.mrc', .false., '')
        ! parameter input/output
        call simulate_subtomogram%set_input('parm_ios', 1,  smpd)
        call simulate_subtomogram%set_input('parm_ios', 2,  nptcls)
        call simulate_subtomogram%set_input('parm_ios', 3,  'snr', 'num', 'SNR', 'Signal-to-noise ratio of particle images', 'signal-to-noise ratio(0.)', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_subtomogram%set_input('comp_ctrls', 1, nthr)
    end subroutine new_simulate_subtomogram

    subroutine new_stack
        ! PROGRAM SPECIFICATION
        call stack%new(&
        &'stack',&                     ! name
        &'Stack images',&              ! descr_short
        &'is a program for stacking individual images (list) or multiple stacks into one',& ! descr_long
        &'simple_exec',&               ! executable
        &2, 2, 0, 0, 0, 0, 0, .false.) ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call stack%set_input('img_ios', 1, 'filetab', 'file', 'Stacks list',&
        &'List of stacks of images to stack into one', 'list input e.g. stktab.txt', .true., '')
        call stack%set_input('img_ios', 2, outstk)
        ! parameter input/output
        call stack%set_input('parm_ios', 1, smpd)
        call stack%set_input('parm_ios', 2, clip)
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
    end subroutine new_stack

    subroutine new_stackops
        ! PROGRAM SPECIFICATION
        call stackops%new(&
        &'stack',&                                   ! name
        &'Standard stack editing',&                  ! descr_short
        &'is a program that provides standard single-particle image processing routines for MRC or SPIDER '&
        &'stacks. To extract a particular state, give oritab and set state. '&
        &'To select the fraction of best particles, give oritab '&
        &'and set frac. State and frac options can be combined. '&
        &'To apply noise, give the desired signal-to-noise ratio via snr. To calculate the autocorrelation '&
        &'function, set acf=yes. To extract a contiguous subset of images from the stack, set '&
        &'fromp and top. To extract a number of particle images at random, set nran to the desired number.'&
        &'With avg=yes the global average of the stack is calculated. '&
        &'If nframesgrp is set to some integer number >1, averages with chunk sizes of nframesgrp are produced, '&
        &'which may be useful for analysis of dose-fractionated image series. neg inverts the contrast of the images',& ! descr_long
        &'simple_exec',&                             ! executable
        &2, 17, 0, 0, 0, 0, 1, .false.)              ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call stackops%set_input('img_ios', 1, stk)
        stackops%img_ios(1)%required = .true.
        call stackops%set_input('img_ios', 2, outstk)
        ! parameter input/output
        call stackops%set_input('parm_ios', 1,  smpd)
        call stackops%set_input('parm_ios', 2,  oritab)
        call stackops%set_input('parm_ios', 3,  mirr)
        call stackops%set_input('parm_ios', 4,  'nran',  'num', 'Number of random samples', 'Number of images to randomly sample', '# random samples', .false., 0.)
        call stackops%set_input('parm_ios', 5,  frac)
        call stackops%set_input('parm_ios', 6,  'state', 'num', 'State index', 'Index of state to extract', 'give state index', .false., 1.)
        call stackops%set_input('parm_ios', 7,  'class', 'num', 'Class index', 'Index of class to extract', 'give class index', .false., 1.)
        call stackops%set_input('parm_ios', 8,  neg)
        call stackops%set_input('parm_ios', 9,  'acf',   'binary', 'Autocorrelation, A * conjg(A)', 'Generate autocorrelation function: A * conjg(A)(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%set_input('parm_ios', 10, 'avg',   'binary', 'Average stack', 'Generate global stack average(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%set_input('parm_ios', 11, 'nframesgrp', 'num', 'Number of stack entries to group & average', 'Number of stack entries to group & average{0}', '# frames', .false., 0.)
        call stackops%set_input('parm_ios', 12, 'vis',   'binary', 'Visualize stack images', 'Visualize stack images with gnuplot(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%set_input('parm_ios', 13, 'snr',   'num', 'Apply noise to give SNR', 'Apply noise to give this signal-to-noise ratio of output', 'signal-to-noise ratio(0.)', .false., 0.)
        call stackops%set_input('parm_ios', 14, 'fromp', 'num', 'From particle index', 'Start index for stack copy', 'start index', .false., 1.0)
        call stackops%set_input('parm_ios', 15, 'top',   'num', 'To particle index', 'Stop index for stack copy', 'stop index', .false., 1.0)
        call stackops%set_input('parm_ios', 16, outfile)
        call stackops%set_input('parm_ios', 17, 'stats', 'binary', 'Provide statistics', 'Provide statistics about images in stack(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call stackops%set_input('comp_ctrls', 1, nthr)
    end subroutine new_stackops

    subroutine new_symsrch
        ! PROGRAM SPECIFICATION
        call symsrch%new(&
        &'symsrch',&                                                                                               ! name
        &'Search for symmetry axis',&                                                                              ! descr_short
        &'is a distributed workflow for searching for the principal symmetry axis of a volume reconstructed in C1. &
        &The 3D reconstruction is projected and common lines-based optimisation is used to identify the principal &
        &symmetry axis. The rotational transformation is applied to the inputted orientations and a new document &
        &is produced. If you are unsure about the point-group, use the compare=yes mode and input the highest &
        &conceviable point-group. The program then calculates probabilities for all lower groups inclusive',&      ! descr_long
        &'simple_exec',&                                                                                           ! executable
        &1, 4, 0, 3, 3, 1, 2, .false.)                                                                             ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symsrch%set_input('img_ios', 1, 'vol1', 'file', 'C1 Volume to identify symmetry axis of', 'C1 Volume to identify symmetry axis of', &
        & 'input volume e.g. vol_C1.mrc', .true., '')
        ! parameter input/output
        call symsrch%set_input('parm_ios', 1, smpd)
        call symsrch%set_input('parm_ios', 2, oritab)
        call symsrch%set_input('parm_ios', 3, outfile)
        call symsrch%set_input('parm_ios', 4, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symsrch%set_input('srch_ctrls', 1, pgrp)
        call symsrch%set_input('srch_ctrls', 2, nspace)
        call symsrch%set_input('srch_ctrls', 3, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symsrch%set_input('filt_ctrls', 1, lp)
        call symsrch%set_input('filt_ctrls', 2, hp)
        call symsrch%set_input('filt_ctrls', 3, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        call symsrch%set_input('mask_ctrls', 1, msk)
        ! computer controls
        call symsrch%set_input('comp_ctrls', 1, nparts)
        call symsrch%set_input('comp_ctrls', 2, nthr)
    end subroutine new_symsrch

    subroutine new_tseries_track
        ! PROGRAM SPECIFICATION
        call tseries_track%new(&
        &'tseries_track',&                                                       ! name
        &'Track particles in time-series',&                                      ! descr_short
        &'is a distributed workflow for particle tracking in time-series data',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &1, 4, 0, 1, 2, 0, 1, .false.)                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tseries_track%set_input('img_ios', 1, 'filetab', 'file', 'List of of movie frames',&
        &'List of movie frames to track in', 'e.g. movie_frames.txt', .true., '')
        ! parameter input/output
        call tseries_track%set_input('parm_ios', 1, 'fbody', 'string', 'Template output tracked series',&
        &'Template output tracked series', 'e.g. tracked_ptcl', .true., '')
        call tseries_track%set_input('parm_ios', 2, smpd)
        call tseries_track%set_input('parm_ios', 3, 'boxfile', 'file', 'List of particle coordiantes',&
        &'.txt file with EMAN particle coordinates', 'e.g. coords.box', .true., '')
        call tseries_track%set_input('parm_ios', 4, mkdir_)
        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_track%set_input('srch_ctrls', 1, 'offset', 'num', 'Shift half-width search bound', 'Shift half-width search bound(in pixels)',&
        'e.g. pixels window halfwidth', .false., 10.)
        ! <empty>
        ! filter controls
        call tseries_track%set_input('filt_ctrls', 1, lp)
        tseries_track%filt_ctrls(1)%required = .false.
        call tseries_track%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the particle and centering', 'centering low-pass limit in Angstroms{5}', .false., 5.)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_track%set_input('comp_ctrls', 1, nparts)
    end subroutine new_tseries_track

    subroutine new_update_project
        ! PROGRAM SPECIFICATION
        call update_project%new(&
        &'new_project',&                     ! name
        &'Update an existing project',&      ! descr_short
        &'is a program for updating an existing project: changing the name/user_email/computer controls',& ! descr_long
        &'simple_exec',&                     ! executable
        &0, 2, 0, 0, 0, 0, 7, .false.)       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call update_project%set_input('parm_ios', 1, projname)
        call update_project%set_input('parm_ios', 2, user_email)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call update_project%set_input('comp_ctrls', 1, time_per_image)
        call update_project%set_input('comp_ctrls', 2, user_account)
        call update_project%set_input('comp_ctrls', 3, user_project)
        call update_project%set_input('comp_ctrls', 4, qsys_partition)
        call update_project%set_input('comp_ctrls', 5, qsys_qos)
        call update_project%set_input('comp_ctrls', 6, qsys_reservation)
        call update_project%set_input('comp_ctrls', 7, job_memory_per_task)
    end subroutine new_update_project

    subroutine new_vizoris
        ! PROGRAM SPECIFICATION
        call vizoris%new(&
        &'vizoris',&                                                                                               ! name
        &'Visualization of orientation distribution',&                                                             ! descr_short
        &'is a program for extracting projection directions from orientations for visualization in UCSF Chimera',& ! descr_long
        &'simple_exec',&                                                                                           ! executable
        &0, 5, 0, 0, 0, 0, 0, .false.)                                                                             ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call vizoris%set_input('parm_ios', 1, oritab)
        vizoris%parm_ios(1)%required = .true.
        call vizoris%set_input('parm_ios', 2, nspace)
        call vizoris%set_input('parm_ios', 3, pgrp)
        call vizoris%set_input('parm_ios', 4, oritype)
        call vizoris%set_input('parm_ios', 5, 'tseries', 'binary', 'Time series anlysis', 'Orientations come from analysis of a time-series(yes|no){no}', '(yes|no){no}', .false., 'no')
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
    end subroutine new_vizoris

    subroutine new_volops
        ! PROGRAM SPECIFICATION
        call volops%new(&
        &'volops',&                                                                               ! name
        &'Standard volume editing',&                                                              ! descr_short
        &'is a program that provides standard single-particle image processing routines for MRC or SPIDER volumes',& ! descr_long
        &'simple_exec',&                                                                          ! executable
        &2, 12, 0, 0, 3, 1, 1, .false.)                                                           ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call volops%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to mask', &
        & 'input volume e.g. vol.mrc', .true., '')
        call volops%set_input('img_ios', 2, outvol)
        ! ! parameter input/output
        call volops%set_input('parm_ios', 1, smpd)
        volops%parm_ios(1)%required = .false.
        call volops%set_input('parm_ios', 2, 'guinier', 'binary', 'Guinier plot','calculate Guinier plot(yes|no){no}', '(yes|no){no}', .false., 'no')
        call volops%set_input('parm_ios', 3, 'neg', 'binary', 'Invert contrast', 'Invert volume contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call volops%set_input('parm_ios', 4, 'snr', 'num', 'SNR','Adds noise to the volume', 'signal-to-noise ratio(0.)', .false., 0.)
        call volops%set_input('parm_ios', 5, mirr)
        call volops%set_input('parm_ios', 6, e1)
        call volops%set_input('parm_ios', 7, e2)
        call volops%set_input('parm_ios', 8, e3)
        call volops%set_input('parm_ios', 9, 'xsh', 'num', 'Translation along x-axis','Shift along X in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',10, 'ysh', 'num', 'Translation along y-axis','Shift along Y in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',11, 'zsh', 'num', 'Translation along z-axis','Shift along Z in pixels', 'in pixels', .false., 0.)
        call volops%set_input('parm_ios',12, outfile)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call volops%set_input('filt_ctrls', 1, lp)
        volops%filt_ctrls(1)%required = .false.
        call volops%set_input('filt_ctrls', 2, hp)
        call volops%set_input('filt_ctrls', 3, bfac)
        ! mask controls
        call volops%set_input('mask_ctrls', 1, msk)
        volops%mask_ctrls(1)%required = .false.
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
                integer,                  intent(in)    :: i
                type(simple_input_param), intent(inout) :: arr(i)
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
        use simple_ansi_ctrls
        class(simple_program), intent(in) :: self
        type(chash) :: ch
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
        use simple_ansi_ctrls
        class(simple_program), intent(in) :: self
        type(chash) :: ch
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
                character(len=*),          intent(in) :: name
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                type(json_value), pointer :: entry, section
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

    character(len=STDLEN) function get_name( self )
        class(simple_program), intent(in) :: self
        get_name = trim(self%name)
    end function get_name

    integer function get_nrequired_keys( self )
        class(simple_program), intent(in) :: self
        get_nrequired_keys = nreq_counter(self%img_ios) + nreq_counter(self%parm_ios) +&
        &nreq_counter(self%srch_ctrls) + nreq_counter(self%filt_ctrls) +&
        &nreq_counter(self%mask_ctrls) + nreq_counter(self%comp_ctrls)
        if( get_nrequired_keys == 0 .and. allocated(self%alt_ios) ) get_nrequired_keys = 1

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

    end function get_nrequired_keys

    function get_required_keys( self ) result( keys )
        class(simple_program), intent(in) :: self
        type(str4arr), allocatable :: keys(:)
        integer :: nreq, ireq
        ! count # required
        nreq = self%get_nrequired_keys()
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
