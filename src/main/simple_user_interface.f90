module simple_user_interface
include 'simple_lib.f08'
implicit none

public :: simple_program, make_user_interface, get_prg_ptr, list_simple_prgs_in_ui
public :: write_ui_json, print_ui_latex, list_single_prgs_in_ui
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG = .false.

type simple_input_param
    character(len=:), allocatable :: key
    character(len=:), allocatable :: keytype ! (binary|multi|num|str|file|dir)
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: descr_placeholder
    character(len=:), allocatable :: gui_submenu
    character(len=:), allocatable :: active_flags
    character(len=:), allocatable :: exclusive_group
    character(len=:), allocatable :: cval_default
    real                          :: rval_default = 0.
    logical :: required = .true.
    logical :: advanced = .true.
    logical :: online   = .false.
end type simple_input_param

type :: simple_program
    private
    character(len=:), allocatable :: name
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: executable
    character(len=:), allocatable :: gui_submenu_list
    logical                       :: advanced = .true.
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
    procedure, private :: set_gui_params
    procedure          :: print_ui
    procedure          :: print_cmdline
    procedure          :: print_cmdline_latex
    procedure          :: print_prg_descr_long
    procedure          :: write2json
    procedure          :: get_name
    procedure          :: get_executable
    procedure          :: get_nrequired_keys
    procedure          :: get_required_keys
    procedure          :: requires_sp_project
    procedure, private :: kill
end type simple_program

! declare simple_exec and single_exec program specifications here
! instances of this class - special

type(simple_program), target :: analysis2D_nano
type(simple_program), target :: assign_optics_groups
type(simple_program), target :: automask
type(simple_program), target :: automask2D
type(simple_program), target :: autorefine3D_nano
type(simple_program), target :: binarize
type(simple_program), target :: calc_pspec
type(simple_program), target :: ced_2D_filter
type(simple_program), target :: center
type(simple_program), target :: cleanup2D
type(simple_program), target :: center2D_nano
type(simple_program), target :: cluster2D
type(simple_program), target :: cluster2D_nano
type(simple_program), target :: cluster2D_stream
type(simple_program), target :: cluster2D_subsets
type(simple_program), target :: cluster3D
type(simple_program), target :: cluster3D_refine
type(simple_program), target :: cluster_cavgs
type(simple_program), target :: comparemc
type(simple_program), target :: convert
type(simple_program), target :: ctf_estimate
type(simple_program), target :: ctfops
type(simple_program), target :: detect_atoms
type(simple_program), target :: dock_volpair
type(simple_program), target :: estimate_diam
type(simple_program), target :: export_relion
type(simple_program), target :: export_starproject
type(simple_program), target :: extract
type(simple_program), target :: filter
type(simple_program), target :: fsc
type(simple_program), target :: gen_pspecs_and_thumbs
type(simple_program), target :: import_boxes
type(simple_program), target :: import_cavgs
type(simple_program), target :: import_movies
type(simple_program), target :: import_particles
type(simple_program), target :: import_starproject
type(simple_program), target :: info_image
type(simple_program), target :: info_stktab
type(simple_program), target :: initial_3Dmodel
type(simple_program), target :: make_cavgs
type(simple_program), target :: make_oris
type(simple_program), target :: map_cavgs_selection
type(simple_program), target :: map_cavgs_states
type(simple_program), target :: mask
type(simple_program), target :: mkdir_
type(simple_program), target :: merge_stream_projects
type(simple_program), target :: motion_correct
type(simple_program), target :: new_project
type(simple_program), target :: nununiform_filter2D
type(simple_program), target :: nununiform_filter3D
type(simple_program), target :: normalize_
type(simple_program), target :: orisops
type(simple_program), target :: oristats
type(simple_program), target :: pick
type(simple_program), target :: postprocess
type(simple_program), target :: preprocess
type(simple_program), target :: preprocess_stream
type(simple_program), target :: preprocess_stream_dev
type(simple_program), target :: print_dose_weights
type(simple_program), target :: print_fsc
type(simple_program), target :: print_magic_boxes
type(simple_program), target :: print_project_field
type(simple_program), target :: print_project_info
type(simple_program), target :: prune_project
type(simple_program), target :: atoms_stats
type(simple_program), target :: reconstruct3D
type(simple_program), target :: reextract
type(simple_program), target :: refine3D
type(simple_program), target :: refine3D_nano
type(simple_program), target :: remoc
type(simple_program), target :: replace_project_field
type(simple_program), target :: selection
type(simple_program), target :: reproject
type(simple_program), target :: scale
type(simple_program), target :: scale_project
type(simple_program), target :: select_
type(simple_program), target :: simulate_atoms
type(simple_program), target :: simulate_movie
type(simple_program), target :: simulate_noise
type(simple_program), target :: simulate_particles
type(simple_program), target :: simulate_subtomogram
type(simple_program), target :: split_
type(simple_program), target :: stack
type(simple_program), target :: stackops
type(simple_program), target :: symaxis_search
type(simple_program), target :: symmetrize_map
type(simple_program), target :: symmetry_test
type(simple_program), target :: tseries_atoms_analysis
type(simple_program), target :: tseries_import
type(simple_program), target :: tseries_import_particles
type(simple_program), target :: tseries_make_pickavg
type(simple_program), target :: tseries_motion_correct
type(simple_program), target :: tseries_swap_stack
type(simple_program), target :: tseries_track_particles
type(simple_program), target :: tseries_reconstruct3D
type(simple_program), target :: graphene_subtr
type(simple_program), target :: uniform_filter2D
type(simple_program), target :: uniform_filter3D
type(simple_program), target :: update_project
type(simple_program), target :: vizoris
type(simple_program), target :: volops
type(simple_program), target :: whiten_and_filter
type(simple_program), target :: write_classes

! declare common params here, with name same as flag
type(simple_input_param) :: algorithm
type(simple_input_param) :: angerr
type(simple_input_param) :: astigtol
type(simple_input_param) :: automsk
type(simple_input_param) :: bfac
type(simple_input_param) :: box
type(simple_input_param) :: box_extract
type(simple_input_param) :: clip
type(simple_input_param) :: clustermode
type(simple_input_param) :: cn
type(simple_input_param) :: cn_min
type(simple_input_param) :: cn_max
type(simple_input_param) :: combine_eo
type(simple_input_param) :: cs
type(simple_input_param) :: ctf
type(simple_input_param) :: ctfpatch
type(simple_input_param) :: ctf_yes
type(simple_input_param) :: ctfresthreshold
type(simple_input_param) :: deftab
type(simple_input_param) :: dferr
type(simple_input_param) :: dfmax
type(simple_input_param) :: dfmin
type(simple_input_param) :: dir_movies
type(simple_input_param) :: e1, e2, e3
type(simple_input_param) :: eer_fraction
type(simple_input_param) :: eer_upsampling
type(simple_input_param) :: element
type(simple_input_param) :: eo
type(simple_input_param) :: focusmskdiam
type(simple_input_param) :: frac
type(simple_input_param) :: fraca
type(simple_input_param) :: fraction_dose_target
type(simple_input_param) :: frcs
type(simple_input_param) :: gainref
type(simple_input_param) :: graphene_filt
type(simple_input_param) :: groupframes
type(simple_input_param) :: hp
type(simple_input_param) :: icefracthreshold
type(simple_input_param) :: job_memory_per_task
type(simple_input_param) :: kv
type(simple_input_param) :: lp
type(simple_input_param) :: lp_backgr
type(simple_input_param) :: lp_lowres
type(simple_input_param) :: lplim_crit
type(simple_input_param) :: lpthres
type(simple_input_param) :: max_dose
type(simple_input_param) :: max_rad
type(simple_input_param) :: maxits
type(simple_input_param) :: mcpatch
type(simple_input_param) :: mcpatch_thres
type(simple_input_param) :: mcconvention
type(simple_input_param) :: min_rad
type(simple_input_param) :: mirr
type(simple_input_param) :: moldiam
type(simple_input_param) :: mskdiam
type(simple_input_param) :: mskfile
type(simple_input_param) :: envfsc
type(simple_input_param) :: ml_reg
type(simple_input_param) :: mul
type(simple_input_param) :: mw
type(simple_input_param) :: nchunks
type(simple_input_param) :: ncls
type(simple_input_param) :: ncls_start
type(simple_input_param) :: neg
type(simple_input_param) :: nparts
type(simple_input_param) :: nparts_chunk
type(simple_input_param) :: nparts_pool
type(simple_input_param) :: nptcls
type(simple_input_param) :: nptcls_per_cls
type(simple_input_param) :: nrestarts
type(simple_input_param) :: nsearch
type(simple_input_param) :: nsig
type(simple_input_param) :: nspace
type(simple_input_param) :: nthr
type(simple_input_param) :: nonuniform
type(simple_input_param) :: numlen
type(simple_input_param) :: nxpatch
type(simple_input_param) :: nypatch
type(simple_input_param) :: objfun
type(simple_input_param) :: oritab
type(simple_input_param) :: oritab2
type(simple_input_param) :: oritype
type(simple_input_param) :: outfile
type(simple_input_param) :: outstk
type(simple_input_param) :: outvol
type(simple_input_param) :: pcontrast
type(simple_input_param) :: pgrp
type(simple_input_param) :: phaseplate
type(simple_input_param) :: projfile
type(simple_input_param) :: projfile_target
type(simple_input_param) :: projname
type(simple_input_param) :: prune
type(simple_input_param) :: pspecsz
type(simple_input_param) :: ptclw
type(simple_input_param) :: qsys_name
type(simple_input_param) :: qsys_partition
type(simple_input_param) :: qsys_qos
type(simple_input_param) :: qsys_reservation
type(simple_input_param) :: remap_cls
type(simple_input_param) :: scale_movies
type(simple_input_param) :: script
type(simple_input_param) :: sherr
type(simple_input_param) :: sigma
type(simple_input_param) :: sigma_est
type(simple_input_param) :: smooth_ext
type(simple_input_param) :: smpd
type(simple_input_param) :: star_datadir
type(simple_input_param) :: starfile
type(simple_input_param) :: star_mic
type(simple_input_param) :: star_model
type(simple_input_param) :: star_ptcl
type(simple_input_param) :: startit
type(simple_input_param) :: startype
type(simple_input_param) :: stk
type(simple_input_param) :: stk2
type(simple_input_param) :: stktab
type(simple_input_param) :: stepsz
type(simple_input_param) :: time_per_image
type(simple_input_param) :: total_dose
type(simple_input_param) :: trs
type(simple_input_param) :: tseries
type(simple_input_param) :: update_frac
type(simple_input_param) :: user_account
type(simple_input_param) :: user_email
type(simple_input_param) :: user_project
type(simple_input_param) :: walltime
type(simple_input_param) :: wcrit
type(simple_input_param) :: width
type(simple_input_param) :: wiener

! this is for making an array of pointers to all programs
type simple_prg_ptr
    type(simple_program), pointer :: ptr2prg => null()
end type simple_prg_ptr
integer, parameter   :: NMAX_PTRS  = 200
integer              :: n_prg_ptrs = 0
type(simple_prg_ptr) :: prg_ptr_array(NMAX_PTRS)

interface set_param
    module procedure set_param_1
    module procedure set_param_2
end interface set_param

contains

    ! public class methods

    subroutine make_user_interface
        call set_common_params
        call set_prg_ptr_array
        call new_analysis2D_nano
        call new_assign_optics_groups
        call new_automask
        call new_automask2D
        call new_autorefine3D_nano
        call new_binarize
        call new_calc_pspec
        call new_ced_2D_filter
        call new_center
        call new_cleanup2D
        call new_center2D_nano
        call new_cluster2D
        call new_cluster2D_nano
        call new_cluster2D_stream
        call new_cluster2D_subsets
        call new_cluster3D
        call new_cluster3D_refine
        call new_cluster_cavgs
        call new_comparemc
        call new_convert
        call new_ctf_estimate
        call new_ctfops
        call new_detect_atoms
        call new_dock_volpair
        call new_estimate_diam
        call new_extract
        call new_export_relion
        call new_export_starproject
        call new_filter
        call new_fsc
        call new_gen_pspecs_and_thumbs
        call new_info_image
        call new_info_stktab
        call new_initial_3Dmodel
        call new_import_boxes
        call new_import_cavgs
        call new_import_movies
        call new_import_particles
        call new_import_starproject
        call new_make_cavgs
        call new_make_oris
        call new_map_cavgs_selection
        call new_map_cavgs_states
        call new_mask
        call new_merge_stream_projects
        call new_mkdir_
        call new_motion_correct
        call new_new_project
        call new_nununiform_filter2D
        call new_nununiform_filter3D
        call new_normalize
        call new_orisops
        call new_oristats
        call new_pick
        call new_postprocess
        call new_preprocess
        call new_preprocess_stream
        call new_preprocess_stream_dev
        call new_print_dose_weights
        call new_print_fsc
        call new_print_magic_boxes
        call new_print_project_info
        call new_print_project_field
        call new_prune_project
        call new_atoms_stats
        call new_reproject
        call new_reconstruct3D
        call new_reextract
        call new_refine3D
        call new_refine3D_nano
        call new_remoc
        call new_replace_project_field
        call new_selection
        call new_scale
        call new_scale_project
        call new_select_
        call new_simulate_atoms
        call new_simulate_movie
        call new_simulate_noise
        call new_simulate_particles
        call new_simulate_subtomogram
        call new_split_
        call new_stack
        call new_stackops
        call new_symaxis_search
        call new_symmetrize_map
        call new_symmetry_test
        call new_tseries_atoms_analysis
        call new_tseries_import
        call new_tseries_import_particles
        call new_tseries_motion_correct
        call new_tseries_make_pickavg
        call new_tseries_swap_stack
        call new_tseries_track_particles
        call new_tseries_reconstruct3D
        call new_graphene_subtr
        call new_uniform_filter2D
        call new_uniform_filter3D
        call new_update_project
        call new_vizoris
        call new_volops
        call new_whiten_and_filter
        call new_write_classes
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; make_user_interface, DONE'
    end subroutine make_user_interface

    subroutine set_prg_ptr_array
        n_prg_ptrs = 0
        call push2prg_ptr_array(analysis2D_nano)
        call push2prg_ptr_array(assign_optics_groups)
        call push2prg_ptr_array(automask)
        call push2prg_ptr_array(automask2D)
        call push2prg_ptr_array(autorefine3D_nano)
        call push2prg_ptr_array(binarize)
        call push2prg_ptr_array(calc_pspec)
        call push2prg_ptr_array(ced_2D_filter)
        call push2prg_ptr_array(center)
        call push2prg_ptr_array(cleanup2D)
        call push2prg_ptr_array(center2D_nano)
        call push2prg_ptr_array(cluster2D)
        call push2prg_ptr_array(cluster2D_nano)
        call push2prg_ptr_array(cluster2D_stream)
        call push2prg_ptr_array(cluster2D_subsets)
        call push2prg_ptr_array(cluster3D)
        call push2prg_ptr_array(cluster3D_refine)
        call push2prg_ptr_array(cluster_cavgs)
        call push2prg_ptr_array(comparemc)
        call push2prg_ptr_array(convert)
        call push2prg_ptr_array(ctf_estimate)
        call push2prg_ptr_array(ctfops)
        call push2prg_ptr_array(detect_atoms)
        call push2prg_ptr_array(dock_volpair)
        call push2prg_ptr_array(extract)
        call push2prg_ptr_array(export_relion)
        call push2prg_ptr_array(export_starproject)
        call push2prg_ptr_array(filter)
        call push2prg_ptr_array(fsc)
        call push2prg_ptr_array(gen_pspecs_and_thumbs)
        call push2prg_ptr_array(info_image)
        call push2prg_ptr_array(info_stktab)
        call push2prg_ptr_array(initial_3Dmodel)
        call push2prg_ptr_array(import_boxes)
        call push2prg_ptr_array(import_cavgs)
        call push2prg_ptr_array(import_movies)
        call push2prg_ptr_array(import_particles)
        call push2prg_ptr_array(import_starproject)
        call push2prg_ptr_array(make_cavgs)
        call push2prg_ptr_array(make_oris)
        call push2prg_ptr_array(map_cavgs_selection)
        call push2prg_ptr_array(map_cavgs_states)
        call push2prg_ptr_array(mask)
        call push2prg_ptr_array(mkdir_)
        call push2prg_ptr_array(motion_correct)
        call push2prg_ptr_array(new_project)
        call push2prg_ptr_array(nununiform_filter2D)
        call push2prg_ptr_array(nununiform_filter3D)
        call push2prg_ptr_array(normalize_)
        call push2prg_ptr_array(orisops)
        call push2prg_ptr_array(oristats)
        call push2prg_ptr_array(pick)
        call push2prg_ptr_array(postprocess)
        call push2prg_ptr_array(preprocess)
        call push2prg_ptr_array(preprocess_stream)
        call push2prg_ptr_array(preprocess_stream_dev)
        call push2prg_ptr_array(print_dose_weights)
        call push2prg_ptr_array(print_fsc)
        call push2prg_ptr_array(print_magic_boxes)
        call push2prg_ptr_array(print_project_info)
        call push2prg_ptr_array(print_project_field)
        call push2prg_ptr_array(prune_project)
        call push2prg_ptr_array(atoms_stats)
        call push2prg_ptr_array(reproject)
        call push2prg_ptr_array(reconstruct3D)
        call push2prg_ptr_array(reextract)
        call push2prg_ptr_array(refine3D)
        call push2prg_ptr_array(refine3D_nano)
        call push2prg_ptr_array(replace_project_field)
        call push2prg_ptr_array(selection)
        call push2prg_ptr_array(scale)
        call push2prg_ptr_array(scale_project)
        call push2prg_ptr_array(select_)
        call push2prg_ptr_array(simulate_atoms)
        call push2prg_ptr_array(simulate_movie)
        call push2prg_ptr_array(simulate_noise)
        call push2prg_ptr_array(simulate_particles)
        call push2prg_ptr_array(simulate_subtomogram)
        call push2prg_ptr_array(split_)
        call push2prg_ptr_array(stack)
        call push2prg_ptr_array(stackops)
        call push2prg_ptr_array(symaxis_search)
        call push2prg_ptr_array(symmetrize_map)
        call push2prg_ptr_array(symmetry_test)
        call push2prg_ptr_array(tseries_atoms_analysis)
        call push2prg_ptr_array(tseries_import)
        call push2prg_ptr_array(tseries_import_particles)
        call push2prg_ptr_array(tseries_make_pickavg)
        call push2prg_ptr_array(tseries_motion_correct)
        call push2prg_ptr_array(tseries_swap_stack)
        call push2prg_ptr_array(tseries_track_particles)
        call push2prg_ptr_array(tseries_reconstruct3D)
        call push2prg_ptr_array(graphene_subtr)
        call push2prg_ptr_array(uniform_filter2D)
        call push2prg_ptr_array(uniform_filter3D)
        call push2prg_ptr_array(update_project)
        call push2prg_ptr_array(vizoris)
        call push2prg_ptr_array(volops)
        call push2prg_ptr_array(whiten_and_filter)
        call push2prg_ptr_array(write_classes)
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; set_prg_ptr_array, DONE'
        contains

            subroutine push2prg_ptr_array( prg )
                type(simple_program), target :: prg
                n_prg_ptrs = n_prg_ptrs + 1
                prg_ptr_array(n_prg_ptrs)%ptr2prg => prg
            end subroutine push2prg_ptr_array

    end subroutine set_prg_ptr_array

    subroutine get_prg_ptr( which_program, ptr2prg )
        character(len=*), intent(in)  :: which_program
        type(simple_program), pointer :: ptr2prg
        select case(trim(which_program))
            case('analysis2D_nano')
                ptr2prg => analysis2D_nano
            case('assign_optics_groups')
                ptr2prg => assign_optics_groups
            case('automask')
                ptr2prg => automask
            case('automask2D')
                ptr2prg => automask2D
            case('autorefine3D_nano')
                ptr2prg => autorefine3D_nano
            case('binarize')
                ptr2prg => binarize
            case('calc_pspec')
                ptr2prg => calc_pspec
            case('ced_2D_filter')
                ptr2prg => ced_2D_filter
            case('center')
                ptr2prg => center
            case('cleanup2D')
                ptr2prg => cleanup2D
            case('center2D_nano')
                ptr2prg => center2D_nano
            case('cluster2D')
                ptr2prg => cluster2D
            case('cluster2D_nano')
                ptr2prg => cluster2D_nano
            case('cluster2D_stream')
                ptr2prg => cluster2D_stream
            case('cluster2D_subsets')
                ptr2prg => cluster2D_subsets
            case('cluster3D')
                ptr2prg => cluster3D
            case('cluster3D_refine')
                ptr2prg => cluster3D_refine
            case('cluster_cavgs')
                ptr2prg => cluster_cavgs
            case('comparemc')
                ptr2prg => comparemc
            case('convert')
                ptr2prg => convert
            case('ctf_estimate')
                ptr2prg => ctf_estimate
            case('ctfops')
                ptr2prg => ctfops
            case('detect_atoms')
                ptr2prg => detect_atoms
            case('dock_volpair')
                ptr2prg => dock_volpair
            case('estimate_diam')
                ptr2prg => estimate_diam
            case('extract')
                ptr2prg => extract
            case('export_relion')
                ptr2prg => export_relion
            case('export_starproject')
                ptr2prg => export_starproject
            case('filter')
                ptr2prg => filter
            case('fsc')
                ptr2prg => fsc
            case('gen_pspecs_and_thumbs')
                ptr2prg => gen_pspecs_and_thumbs
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
            case('import_starproject')
                ptr2prg => import_starproject
            case('make_cavgs')
                ptr2prg => make_cavgs
            case('make_oris')
                ptr2prg => make_oris
            case('map_cavgs_selection')
                ptr2prg => map_cavgs_selection
            case('map_cavgs_states')
                ptr2prg => map_cavgs_states
            case('mask')
                ptr2prg => mask
            case('merge_stream_projects')
                ptr2prg => merge_stream_projects
            case('mkdir')
                ptr2prg => mkdir_
            case('motion_correct')
                ptr2prg => motion_correct
            case('new_project')
                ptr2prg => new_project
            case('nununiform_filter2D')
                ptr2prg => nununiform_filter2D
            case('nununiform_filter3D')
                ptr2prg => nununiform_filter3D
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
            case('preprocess')
                ptr2prg => preprocess
            case('preprocess_stream')
                ptr2prg => preprocess_stream
            case('preprocess_stream_dev')
                ptr2prg => preprocess_stream_dev
            case('print_dose_weights')
                ptr2prg => print_dose_weights
            case('print_fsc')
                ptr2prg => print_fsc
            case('print_magic_boxes')
                ptr2prg => print_magic_boxes
            case('print_project_info')
                ptr2prg => print_project_info
            case('print_project_field')
                ptr2prg => print_project_field
            case('prune_project')
                ptr2prg => prune_project
            case('atoms_stats')
                ptr2prg => atoms_stats
            case('reproject')
                ptr2prg => reproject
            case('reconstruct3D')
                ptr2prg => reconstruct3D
            case('reextract')
                ptr2prg => reextract
            case('refine3D')
                ptr2prg => refine3D
            case('refine3D_nano')
                ptr2prg => refine3D_nano
            case('remoc')
                ptr2prg => remoc
            case('replace_project_field')
                ptr2prg => replace_project_field
            case('selection')
                ptr2prg => selection
            case('scale')
                ptr2prg => scale
            case('scale_project')
                ptr2prg => scale_project
            case('select')
                ptr2prg => select_
            case('simulate_atoms')
                ptr2prg => simulate_atoms
            case('simulate_movie')
                ptr2prg => simulate_movie
            case('simulate_noise')
                ptr2prg => simulate_noise
            case('simulate_particles')
                ptr2prg => simulate_particles
            case('simulate_subtomogram')
                ptr2prg => simulate_subtomogram
            case('split')
                ptr2prg => split_
            case('stack')
                ptr2prg => stack
            case('stackops')
                ptr2prg => stackops
            case('symaxis_search')
                ptr2prg => symaxis_search
            case('symmetrize_map')
                ptr2prg => symmetrize_map
            case('symmetry_test')
                ptr2prg => symmetry_test
            case('tseries_atoms_analysis')
                ptr2prg => tseries_atoms_analysis
            case('tseries_import')
                ptr2prg => tseries_import
            case('tseries_import_particles')
                ptr2prg => tseries_import_particles
            case('tseries_make_pickavg')
                ptr2prg => tseries_make_pickavg
            case('tseries_motion_correct')
                ptr2prg => tseries_motion_correct
            case('tseries_swap_stack')
                ptr2prg => tseries_swap_stack
            case('tseries_track_particles')
                ptr2prg => tseries_track_particles
            case('tseries_reconstruct3D')
                ptr2prg => tseries_reconstruct3D
            case('graphene_subtr')
                ptr2prg => graphene_subtr
            case('uniform_filter2D')
                ptr2prg => uniform_filter2D
            case('uniform_filter3D')
                ptr2prg => uniform_filter3D
            case('update_project')
                ptr2prg => update_project
            case('vizoris')
                ptr2prg => vizoris
            case('volops')
                ptr2prg => volops
            case('whiten_and_filter')
                ptr2prg => whiten_and_filter
            case('write_classes')
                ptr2prg => write_classes
            case DEFAULT
                ptr2prg => null()
        end select
    end subroutine get_prg_ptr

    subroutine list_simple_prgs_in_ui
        write(logfhandle,'(A)') assign_optics_groups%name
        write(logfhandle,'(A)') automask%name
        write(logfhandle,'(A)') automask2D%name
        write(logfhandle,'(A)') binarize%name
        write(logfhandle,'(A)') calc_pspec%name
        write(logfhandle,'(A)') ced_2D_filter%name
        write(logfhandle,'(A)') center%name
        write(logfhandle,'(A)') cleanup2D%name
        write(logfhandle,'(A)') cluster_cavgs%name
        write(logfhandle,'(A)') cluster2D%name
        write(logfhandle,'(A)') cluster2D_stream%name
        write(logfhandle,'(A)') cluster2D_subsets%name
        write(logfhandle,'(A)') cluster3D%name
        write(logfhandle,'(A)') cluster3D_refine%name
        write(logfhandle,'(A)') comparemc%name
        write(logfhandle,'(A)') convert%name
        write(logfhandle,'(A)') ctf_estimate%name
        write(logfhandle,'(A)') ctfops%name
        write(logfhandle,'(A)') dock_volpair%name
        write(logfhandle,'(A)') extract%name
        write(logfhandle,'(A)') export_relion%name
        write(logfhandle,'(A)') export_starproject%name
        write(logfhandle,'(A)') filter%name
        write(logfhandle,'(A)') fsc%name
        write(logfhandle,'(A)') gen_pspecs_and_thumbs%name
        write(logfhandle,'(A)') initial_3Dmodel%name
        write(logfhandle,'(A)') info_image%name
        write(logfhandle,'(A)') info_stktab%name
        write(logfhandle,'(A)') import_boxes%name
        write(logfhandle,'(A)') import_cavgs%name
        write(logfhandle,'(A)') import_movies%name
        write(logfhandle,'(A)') import_particles%name
        write(logfhandle,'(A)') import_starproject%name
        write(logfhandle,'(A)') make_cavgs%name
        write(logfhandle,'(A)') make_oris%name
        write(logfhandle,'(A)') map_cavgs_selection%name
        write(logfhandle,'(A)') map_cavgs_states%name
        write(logfhandle,'(A)') mask%name
        write(logfhandle,'(A)') merge_stream_projects%name
        write(logfhandle,'(A)') mkdir_%name
        write(logfhandle,'(A)') motion_correct%name
        write(logfhandle,'(A)') new_project%name
        write(logfhandle,'(A)') nununiform_filter2D%name
        write(logfhandle,'(A)') nununiform_filter3D%name
        write(logfhandle,'(A)') normalize_%name
        write(logfhandle,'(A)') orisops%name
        write(logfhandle,'(A)') oristats%name
        write(logfhandle,'(A)') pick%name
        write(logfhandle,'(A)') postprocess%name
        write(logfhandle,'(A)') preprocess%name
        write(logfhandle,'(A)') preprocess_stream%name
        write(logfhandle,'(A)') print_dose_weights%name
        write(logfhandle,'(A)') print_fsc%name
        write(logfhandle,'(A)') print_magic_boxes%name
        write(logfhandle,'(A)') print_project_info%name
        write(logfhandle,'(A)') print_project_field%name
        write(logfhandle,'(A)') prune_project%name
        write(logfhandle,'(A)') reconstruct3D%name
        write(logfhandle,'(A)') reextract%name
        write(logfhandle,'(A)') refine3D%name
        write(logfhandle,'(A)') refine3D_nano%name
        write(logfhandle,'(A)') remoc%name
        write(logfhandle,'(A)') replace_project_field%name
        write(logfhandle,'(A)') selection%name
        write(logfhandle,'(A)') reproject%name
        write(logfhandle,'(A)') select_%name
        write(logfhandle,'(A)') simulate_movie%name
        write(logfhandle,'(A)') simulate_noise%name
        write(logfhandle,'(A)') simulate_particles%name
        write(logfhandle,'(A)') simulate_subtomogram%name
        write(logfhandle,'(A)') scale%name
        write(logfhandle,'(A)') scale_project%name
        write(logfhandle,'(A)') split_%name
        write(logfhandle,'(A)') stack%name
        write(logfhandle,'(A)') stackops%name
        write(logfhandle,'(A)') symaxis_search%name
        write(logfhandle,'(A)') symmetrize_map%name
        write(logfhandle,'(A)') symmetry_test%name
        write(logfhandle,'(A)') uniform_filter2D%name
        write(logfhandle,'(A)') uniform_filter3D%name
        write(logfhandle,'(A)') update_project%name
        write(logfhandle,'(A)') vizoris%name
        write(logfhandle,'(A)') volops%name
        write(logfhandle,'(A)') whiten_and_filter%name
        write(logfhandle,'(A)') write_classes%name
    end subroutine list_simple_prgs_in_ui

    subroutine list_single_prgs_in_ui
        write(logfhandle,'(A)') format_str('PROJECT MANAGEMENT PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') new_project%name
        write(logfhandle,'(A)') update_project%name
        write(logfhandle,'(A)') print_project_info%name
        write(logfhandle,'(A)') print_project_field%name
        write(logfhandle,'(A)') tseries_import%name
        write(logfhandle,'(A)') import_particles%name
        write(logfhandle,'(A)') tseries_import_particles%name
        write(logfhandle,'(A)') prune_project%name
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('TIME-SERIES PRE-PROCESSING PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') tseries_make_pickavg%name
        write(logfhandle,'(A)') tseries_motion_correct%name
        write(logfhandle,'(A)') tseries_track_particles%name
        write(logfhandle,'(A)') graphene_subtr%name
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('PARTICLE 3D RECONSTRUCTION PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') analysis2D_nano%name
        write(logfhandle,'(A)') center2D_nano%name
        write(logfhandle,'(A)') cluster2D_nano%name
        write(logfhandle,'(A)') map_cavgs_selection%name
        write(logfhandle,'(A)') estimate_diam%name
        write(logfhandle,'(A)') simulate_atoms%name
        write(logfhandle,'(A)') refine3D_nano%name
        write(logfhandle,'(A)') autorefine3D_nano%name
        write(logfhandle,'(A)') tseries_reconstruct3D%name
        write(logfhandle,'(A)') tseries_swap_stack%name
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('VALIDATION PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') vizoris%name
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('MODEL BULDING/ANALYSIS PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') detect_atoms%name
        write(logfhandle,'(A)') atoms_stats%name
        write(logfhandle,'(A)') tseries_atoms_analysis%name
    end subroutine list_single_prgs_in_ui

    ! private class methods

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

    subroutine set_common_params
        call set_param(projfile,      'projfile',      'file',   'Project file', 'SIMPLE projectfile', 'e.g. myproject.simple', .true., '')
        call set_param(projfile_target,'projfile_target','file', 'Another project file', 'SIMPLE projectfile', 'e.g. myproject2.simple', .true., '')
        call set_param(stk,           'stk',           'file',   'Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., '')
        call set_param(stk2,          'stk2',          'file',   'Second Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., 'stk2.mrc')
        call set_param(stktab,        'stktab',        'file',   'List of per-micrograph particle stacks', 'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false., 'stktab.txt')
        call set_param(ctf,           'ctf',           'multi',  'CTF status', 'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){no}',&
        &'(yes|no|flip){no}', .true., 'no')
        call set_param(ctf_yes,       'ctf',           'multi',  'CTF status', 'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){yes}', '(yes|no|flip){yes}', .false., 'yes')
        call set_param(ctfpatch,      'ctfpatch',      'binary', 'Patch CTF estimation', 'Whether to perform patch CTF estimation(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(ctfresthreshold,'ctfresthreshold','num',  'CTF Resolution rejection threshold', 'Micrographs with a CTF resolution above the threshold (in Angs) will be ignored from further processing{50}', 'CTF resolution threshold(in Angstroms){50}', .false., 50.0)
        call set_param(icefracthreshold,'icefracthreshold','num','Ice Fraction rejection threshold', 'Micrographs with an ice ring/1st pspec maxima fraction above the threshold will be ignored from further processing{1.0}', 'Ice fraction threshold{1.0}', .false., 1.0)
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
        call set_param(dir_movies,    'dir_movies',    'dir',    'Input movies directory', 'Where the movies to process are located or will squentially appear', 'e.g. /cryodata/', .true., 'preprocess/')
        call set_param(startit,       'startit',       'num',    'First iteration', 'Index of first iteration when starting from a previous run', 'start iterations from here', .false., 1.0)
        call set_param(trs,           'trs',           'num',    'Maximum translational shift', 'Maximum half-width for bund-constrained search of rotational origin shifts',&
        &'max shift per iteration in pixels{5}', .false., 5.0)
        call set_param(maxits,        'maxits',        'num',    'Max iterations', 'Maximum number of iterations', 'Max # iterations', .false., 100.)
        call set_param(hp,            'hp',            'num',    'High-pass limit', 'High-pass resolution limit', 'high-pass limit in Angstroms', .false., 100.)
        call set_param(lp,            'lp',            'num',    'Low-pass limit', 'Low-pass resolution limit', 'low-pass limit in Angstroms', .false., 20.)
        call set_param(lp_backgr,     'lp_backgr',     'num',    'Background low-pass resolution', 'Low-pass resolution for solvent blurring', 'low-pass limit in Angstroms', .false., 20.)
        call set_param(mskdiam,       'mskdiam',       'num',    'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to remove background noise', 'mask diameter in A', .true., 0.)
        call set_param(ncls,          'ncls',          'num',    'Number of 2D clusters', 'Number of groups to sort the particles &
        &into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .true., 200.)
        call set_param(ncls_start,    'ncls_start',    'num',    'Number of 2D clusters per subset of particles', 'Number of class averages used in the independent 2D classification of each subset of particles', '# 2D clusters / subset', .true., 50.)
        call set_param(nparts,        'nparts',        'num',    'Number of computing nodes', 'Number of partitions for distributed memory execution. One part typically corresponds to one CPU socket in the distributed &
        &system. On a single-socket machine there may be speed benefits to dividing the jobs into a few (2-4) partitions, depending on memory capacity', 'divide job into # parts', .true., 1.0)
        call set_param(nparts_chunk,  'nparts_chunk',  'num',    'Number of computing nodes per subset', 'Number of computing nodes allocated to 2D classification of each particles subset (chunk){1}', '# of nodes per subset{1}', .false., 1.0)
        call set_param(nparts_pool,   'nparts_pool',   'num',    'Number of computing nodes for the pooled subsets', 'Number of computing nodes allocated to 2D classification of the pooled particles subsets', '# of nodes for the pooled subsets', .false., 2.0)
        call set_param(nthr,          'nthr',          'num',    'Number of threads per computing node, give 0 if unsure', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .true., 0.)
        call set_param(nonuniform,    'nonuniform',    'binary', 'Nonuniform filter', 'Apply nonuniform filter(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(update_frac,   'update_frac',   'num',    'Fractional update per iteration', 'Fraction of particles to update per iteration in incremental learning scheme for accelerated convergence &
        &rate(0.1-0.5){1.}', 'update this fraction per iter(0.1-0.5){1.0}', .false., 1.0)
        call set_param(frac,          'frac',          'num',    'Fraction of particles to include', 'Fraction of particles to include based on spectral score (median of FRC between reference and particle)',&
        'fraction of particles(0.1-0.9){1.0}', .false., 1.0)
        call set_param(mskfile,       'mskfile',       'file',   'Input mask file', 'Input mask file to apply to reference volume(s) before projection', 'e.g. automask.mrc from postprocess', .false., 'mskfile.mrc')
        call set_param(pgrp,          'pgrp',          'str',    'Point-group symmetry', 'Point-group symmetry of particle(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .true., 'c1')
        call set_param(prune,         'prune',         'binary', 'Automated particles pruning', 'Whether to prune deselected particles(yes|no){no}', 'Automated particles pruning(yes|no){no}', .false., 'no')
        call set_param(nspace,        'nspace',        'num',    'Number of projection directions', 'Number of projection directions &
        &used', '# projections', .false., 2500.)
        call set_param(objfun,        'objfun',        'multi',  'Objective function', 'Objective function(euclid|cc|prob){euclid}', '(euclid|cc|prob){euclid}', .false., 'euclid')
        call set_param(remap_cls,     'remap_cls',     'binary', 'Whether to remap 2D clusters', 'Whether to remap the number of 2D clusters(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(kv,            'kv',            'num',    'Acceleration voltage', 'Acceleration voltage in kV{300}', 'in kV{300}', .false., 300.)
        call set_param(lplim_crit,    'lplim_crit',    'num',    'Low-pass limit FSC criterion', 'FSC criterion for determining the low-pass limit(0.143-0.5){0.143}',&
        &'low-pass FSC criterion(0.143-0.5){0.143}', .false., 0.143)
        call set_param(cs,            'cs',            'num',    'Spherical aberration', 'Spherical aberration constant(in mm){2.7}', 'in mm{2.7}', .false., 2.7)
        call set_param(total_dose,    'total_dose',     'num',    'Total exposure dose (e/Ang^2)', 'Total exposure dose (e/Ang^2)', 'in e/Ang^2', .false., 50.)
        call set_param(fraction_dose_target,'fraction_dose_target','num','EER fraction dose target (e/Ang^2)', 'EER fraction dose target, used to determine how many EER frames are included in each movie fraction(e/Ang^2)', 'in e/Ang^2', .false., 1.)
        call set_param(fraca,         'fraca',         'num',    'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.1}', 'fraction{0.1}', .false., 0.1)
        call set_param(pspecsz,       'pspecsz',       'num',    'Size of power spectrum', 'Size of power spectrum in pixels{512}', 'give # pixels{512}', .false., 512.)
        call set_param(dfmin,         'dfmin',         'num',    'Expected minimum defocus', 'Expected minimum defocus in microns{0.2}', 'in microns{0.2}', .false., DFMIN_DEFAULT)
        call set_param(dfmax,         'dfmax',         'num',    'Expected maximum defocus', 'Expected maximum defocus in microns{5.0}', 'in microns{5.0}', .false., DFMAX_DEFAULT)
        call set_param(astigtol,      'astigtol',      'num',    'Expected astigmatism', 'expected (tolerated) astigmatism(in microns){0.05}', 'in microns{0.05}',  .false., 0.05)
        call set_param(mw,            'mw',            'num',    'Molecular weight','Molecular weight in kDa', 'in kDa', .false., 0.)
        call set_param(mirr,          'mirr',          'multi',  'Perform mirroring', 'Whether to mirror and along which axis(no|x|y){no}', '(no|x|y){no}', .false., 'no')
        call set_param(bfac,          'bfac',          'num',    'B-factor for sharpening','B-factor for sharpening in Angstroms^2', 'B-factor in Angstroms^2', .false., 200.)
        call set_param(outvol,        'outvol',        'file',   'Output volume name', 'Output volume name', 'e.g. outvol.mrc', .false., '')
        call set_param(eo,            'eo',            'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(job_memory_per_task, 'job_memory_per_task','str', 'Memory per computing node', 'Memory in MB per part/computing node in distributed execution{16000}', 'MB per part{16000}', .false., 16000.)
        call set_param(qsys_name,     'qsys_name',     'multi',  'Queue system kind', 'Queue system kind(local|slurm|pbs)', '(local|slurm|pbs)', .false., 'local')
        call set_param(qsys_partition,'qsys_partition','str',    'Name of SLURM/PBS partition', 'Name of target partition of distributed computer system (SLURM/PBS)', 'give partition name', .false., '')
        call set_param(qsys_qos,      'qsys_qos',      'str',    'Schedule priority', 'Job scheduling priority (SLURM/PBS)', 'give priority', .false., '')
        call set_param(qsys_reservation, 'qsys_reservation', 'str', 'Name of reserved partition', 'Name of reserved target partition of distributed computer system (SLURM/PBS)', 'give your part', .false., '')
        call set_param(box,            'box',          'num',    'Particle box size', 'Particle box size(in pixels)', '# pixels of box', .true., 0.)
        call set_param(box_extract,    'box_extract',  'num',    'Extracted particle image size', 'Extracted particle image size(in pixels)', 'Extracted particle image size', .false., 0.)
        call set_param(nchunks,        'nchunks',      'num',    'Number of subsets to classify simultaneously', 'Maximum number of particles subsets (chunks) to classify simultaneously', '# of substsets', .true., 2.)
        call set_param(nptcls,         'nptcls',       'num',    'Number of particles', 'Number of particle images', '# particles', .true., 0.)
        call set_param(nptcls_per_cls, 'nptcls_per_cls','num',   'Number of particles per cluster', 'Initial number of particles per cluster{35}', '# initial particles per cluster{35}', .false., 35.)
        call set_param(outstk,         'outstk',       'file',   'Output stack name', 'Output images stack name', 'e.g. outstk.mrc', .false., '')
        call set_param(pcontrast,      'pcontrast',    'multi',  'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false., 'black')
        call set_param(clip,           'clip',         'num',    'Clipped box size', 'Target box size for clipping in pixels', 'in pixels', .false., 0.)
        call set_param(clustermode,    'clustermode',  'multi',  'Feature used for clustering', 'Feature used for clustering (ar|dist|ang|maxint|intint){ar}', '(ar|dist|ang|maxint|intint){ar}', .true., 'ar')
        call set_param(neg,            'neg',          'binary', 'Invert contrast','Invert contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(sherr,          'sherr',        'num',    'Shift error half-width', 'Uniform rotational origin shift error half-width(in pixels)', 'shift error in pixels', .false., 0.)
        call set_param(angerr,         'angerr',       'num',    'Rotation angle error half-width', 'Uniform rotation angle shift error half-width(in degrees)', 'rotation error in degrees', .false., 0.)
        call set_param(dferr,          'dferr',        'num',    'Underfocus error half-width',  'Uniform underfoucs error half-width(in microns)',  'defocus error in microns', .false., 1.)
        call set_param(oritype,        'oritype',      'multi',  'Oritype segment in project',  'Oritype segment in project(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}',&
        &'(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}', .false., 'ptcl3D')
        call set_param(e1,             'e1',           'num',    'Rotation along Phi',  'Phi Euler angle',   'in degrees', .false., 0.)
        call set_param(e2,             'e2',           'num',    'Rotation along Theta','Theat Euler angle', 'in degrees', .false., 0.)
        call set_param(e3,             'e3',           'num',    'Rotation along Psi',  'Psi Euler angle',   'in degrees', .false., 0.)
        call set_param(eer_fraction,   'eer_fraction', 'num',    '# of EER frames to fraction together', 'Number of raw EER frames to fraction together', '# EER frames{20}', .false., 20.)
        call set_param(eer_upsampling, 'eer_upsampling','multi', 'EER up-sampling', 'EER up-sampling(1=4K|2=8K){1}', '(1|2){1}', .false., 1.)
        call set_param(gainref,        'gainref',      'file',   'Gain reference', 'Gain reference image', 'input image e.g. gainref.mrc', .false., '')
        call set_param(groupframes,    'groupframes',  'binary', 'Patch motion correction frames averaging', 'Whether to perform frames averaging during motion correction - for patchesonly(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(mcpatch,        'mcpatch',      'binary', 'Patch-based motion correction', 'Whether to perform Patch-based motion correction(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(mcpatch_thres,'mcpatch_thres','binary','Use motion correction patch threshold', 'Whether to use the threshold for motion correction patch solution(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(mcconvention,   'mcconvention', 'str',    'Frame of reference during movie alignment', 'Frame of reference during movie alignment; simple/unblur:central; relion/motioncorr:first(simple|unblur|relion|motioncorr){simple}', '(simple|unblur|relion|motioncorr){simple}', .false., 'simple')
        call set_param(nxpatch,        'nxpatch',      'num',    '# of patches along x-axis', 'Motion correction # of patches along x-axis', '# x-patches{5}', .false., 5.)
        call set_param(nypatch,        'nypatch',      'num',    '# of patches along y-axis', 'Motion correction # of patches along y-axis', '# y-patches{5}', .false., 5.)
        call set_param(numlen,         'numlen',       'num',    'Length of number string', 'Length of number string', '# characters', .false., 5.0)
        call set_param(nsig,           'nsig',         'num',    'Number of sigmas for outlier removal', 'Number of standard deviations threshold for pixel outlier removal{6}', '# standard deviations{6}', .false., 6.)
        call set_param(projname,       'projname',     'str',    'Project name', 'Name of project to create ./myproject/myproject.simple file for',&
        &'e.g. to create ./myproject/myproject.simple', .true., '')
        call set_param(user_email,     'user_email',   'str',    'Your e-mail address', 'Your e-mail address', 'e.g. myname@uni.edu', .false., '')
        call set_param(time_per_image, 'time_per_image','num',   'Time per image', 'Estimated time per image in seconds for forecasting total execution time{100}', 'in seconds{100}', .false., 100.)
        call set_param(walltime,       'walltime',     'num',    'Walltime', 'Maximum execution time for job scheduling and management(23h59mins){86340}', 'in seconds(23h59mins){86340}', .false., 86340.)
        call set_param(user_account,   'user_account', 'str',    'User account name in SLURM/PBS', 'User account name in SLURM/PBS system', 'e.g. Account084', .false., '')
        call set_param(user_project,   'user_project', 'str',    'User project name in SLURM/PBS', 'User project name in SLURM/PBS system', 'e.g. Project001', .false., '')
        call set_param(frcs,           'frcs',         'str',    'Projection FRCs file', 'Projection FRCs file', 'e.g. frcs.bin', .false., '')
        call set_param(focusmskdiam,   'focusmskdiam', 'num',    'Mask diameter in focused refinement', 'Mask diameter in A for application of a soft-edged circular mask to remove background noise in focused refinement', 'focused mask diameter in A', .false., 0.)
        call set_param(nrestarts,      'nrestarts',    'num',    'Number of restarts', 'Number of program restarts to execute{1}', '# restarts{1}', .false., 1.0)
        call set_param(star_datadir,   'star_datadir', 'file',   'STAR project data directory', 'Pathname of STAR image/data files', 'e.g. Micrographs', .false., '')
        call set_param(starfile,       'starfile',     'file',   'STAR-format file name', 'STAR-formatted filename', 'e.g. proj.star', .false., '')
        call set_param(star_mic,       'star_mic',     'file',   'Micrographs STAR file name', 'Micrographs STAR-formatted filename', 'e.g. micrographs.star', .true., '')
        call set_param(star_model,     'star_model',   'file',   'Model STAR file name', 'Model STAR-formatted filename', 'e.g. model.star', .true., '')
        call set_param(star_ptcl,      'star_ptcl',    'file',   'Particles STAR file name', 'Particles STAR-formatted filename', 'e.g. particles.star', .true., '')
        call set_param(startype,       'startype',     'str',    'STAR-format export type', 'STAR experiment type used to define variables in export file', 'e.g. micrographs or class2d or refine3d', .false., '')
        call set_param(scale_movies,   'scale',        'num',    'Down-scaling factor(0-1)', 'Down-scaling factor to apply to the movies(0-1){1.}', '{1.}', .false., 1.0)
        call set_param(ptclw,          'ptclw',        'binary', 'Particle weights', 'Particle weights(yes|no){no}',  '(yes|no){no}',  .false., 'no')
        call set_param(envfsc,         'envfsc',       'binary', 'Envelope mask e/o maps for FSC', 'Envelope mask even/odd pairs prior to FSC calculation(yes|no){yes}',  '(yes|no){yes}',  .false., 'yes')
        call set_param(graphene_filt,  'graphene_filt','binary', 'Omit graphene bands from corr calc', 'Omit graphene bands from corr calc(yes|no){no}',  '(yes|no){no}',  .false., 'no')
        call set_param(wcrit,          'wcrit',        'multi',  'Correlation to weights conversion scheme', 'Correlation to weights conversion scheme(softmax|zscore|sum|cen|exp|inv|no){softmax}',  '(softmax|zscore|sum|cen|exp|inv|no){softmax}',  .false., 'softmax')
        call set_param(element,        'element',      'str',    'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '  ')
        call set_param(tseries,        'tseries',      'binary', 'Stack is time-series', 'Stack is time-series(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(max_rad,        'max_rad',      'num',    'Maximum radius in A', 'Maximum radius in A {300.}', '{300.}', .true., 300.)
        call set_param(min_rad,        'min_rad',      'num',    'Minimum radius in A', 'Minimum radius in A {50.} ', '{50.}',  .true., 50.)
        call set_param(cn,             'cn',           'num',    'Fixed std coordination number', 'Minimum std cn to consider for dipole calc ', '8',  .false., 8.)
        call set_param(cn_min,         'cn_min',       'num',    'Minimum std coordination number', 'Minimum std cn to consider ', '4',  .false., 4.)
        call set_param(cn_max,         'cn_max',       'num',    'Maximum std coordination number', 'Maximum std cn to consider ', '12', .false., 12.)
        call set_param(stepsz,         'stepsz',       'num',    'Steps size in A', 'Step size in A {10.} ', '{10.}',  .true., 10.)
        call set_param(moldiam,        'moldiam',      'num',    'Molecular diameter', 'Molecular diameter(in Angstroms)','In Angstroms',.false., 0.)
        call set_param(mul,            'mul',          'num',    'Multiplication factor', 'Multiplication factor{1.}','{1.}',.false., 1.)
        call set_param(algorithm,      'algorithm',    'multi',  'Algorithm for motion correction','Algorithm for motion correction(patch|patch_refine){patch}','(patch|patch_refine){patch}', .false.,'patch')
        call set_param(width,          'width',        'num',    'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge{10}', .false., 10.)
        call set_param(automsk,        'automsk',      'multi',  'Perform envelope masking', 'Whether to generate/apply an envelope mask(yes|tight|no){no}', '(yes|tight|no){no}', .false., 'no')
        call set_param(wiener,         'wiener',       'multi',  'Wiener restoration', 'Wiener restoration, full or partial (full|partial){full}','(full|partial){full}', .false., 'full')
        call set_param(max_dose,       'max_dose',     'num',    'Maximum dose threshold(e/A2)', 'Threshold for maximum dose and number of frames used during movie alignment(e/A2), if <=0 all frames are used{0.0}','{0.0}',.false., 0.0)
        call set_param(script,         'script',       'binary', 'Generate script for shared-mem exec on cluster', 'Generate script for shared-mem exec on cluster(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(lp_lowres,      'lp_lowres',    'num',    'Low resolution limit, nonuniform filter', 'Low-pass limit in discrete nonuniform filter search{30 A}', 'input low-pass limit{30 A}', .false., 30.)
        call set_param(nsearch,        'nsearch',      'num',    'Number of points to search in nonuniform filter', 'Number of points to search in discrete nonuniform filter{40}', '# points to search{40}', .false., 40.)
        call set_param(smooth_ext,     'smooth_ext',   'num'   , 'Smoothing window extension', 'Smoothing window extension for nonuniform filter optimization in pixels{20}', 'give # pixels{2D=20,3D=8}', .false., 20.)
        call set_param(lpthres,        'lpthres',      'num',    'Resolution rejection threshold', 'Classes with lower resolution are iteratively rejected in Angstroms{30}', 'give rejection threshold in angstroms{30}', .false., 30.)
        call set_param(ml_reg,         'ml_reg',       'binary', 'ML regularization', 'Regularization (ML-style) based on the signal power(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(sigma_est,      'sigma_est',    'multi',  'Sigma estimation method', 'Sigma estimation method(group|global){group}', '(group|global){group}', .false., 'group')
        call set_param(combine_eo,     'combine_eo',   'binary', 'whether e/o references are combined for final alignment(yes|no){no}', 'whether e/o references are combined for final alignment(yes|no){no}', '(yes|no){no}', .false., 'no')
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; set_common_params, DONE'
    end subroutine set_common_params

    ! CONSTRUCTOR TEMPLATE
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

    subroutine new_analysis2D_nano
        ! PROGRAM SPECIFICATION
        call analysis2D_nano%new(&
        &'analysis2D_nano', &                                         ! name
        &'2D analysis (centering, diameter estimation & clustering) for nanycrystal time-series',& ! descr_short
        &'is a program for 2D analysis for nanycrystal time-series',& ! descr long
        &'single_exec',&                                              ! executable
        &0, 1, 0, 1, 0, 1, 3, .true.)                                 ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call analysis2D_nano%set_input('parm_ios', 1, element)
        ! alternative inputs
        ! <empty>
        ! search controls
        call analysis2D_nano%set_input('srch_ctrls', 1, nptcls_per_cls)
        ! filter controls
        ! <empty>
        ! mask controls
        call analysis2D_nano%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call analysis2D_nano%set_input('comp_ctrls', 1, nparts)
        analysis2D_nano%comp_ctrls(1)%required = .false.
        call analysis2D_nano%set_input('comp_ctrls', 2, nthr)
        call analysis2D_nano%set_input('comp_ctrls', 3, script)
    end subroutine new_analysis2D_nano

    subroutine new_assign_optics_groups
        ! PROGRAM SPECIFICATION
        call assign_optics_groups%new(&
        &'assign_optics_groups', &                                              ! name
        &'Assign optics groups',&                                               ! descr_short
        &'is a program to assign optics groups',&                               ! descr long
        &'simple_exec',&                                                        ! executable
        &0, 4, 0, 0, 0, 0, 0, .true.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics_groups%set_input('parm_ios', 1, 'xmldir', 'dir', 'Directory containing per movie EPU XML files',&
        & 'Directory containing per movie EPY XML files', 'e.g. /data/datasetid/xml', .false., '')
        call assign_optics_groups%set_input('parm_ios', 2, 'maxpop', 'num', 'Maximum number of movies/micrographs/stacks in each optics group',&
        & 'Maximum number of movies/micrographs/stacks in each optics group', 'e.g. 100', .false., '')
        call assign_optics_groups%set_input('parm_ios', 3, 'optics_offset', 'num', 'Numbering offset to apply to optics groups',&
        & 'Numbering offset to apply to optics groups. Aids with combining datasets', 'e.g. 10', .false., '')
        call assign_optics_groups%set_input('parm_ios', 4, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
    end subroutine new_assign_optics_groups

    subroutine new_automask
        ! PROGRAM SPECIFICATION
        call automask%new(&
        &'automask',&                                    ! name
        &'envelope masking',&                            ! descr_short
        &'is a program for automated envelope masking',& ! descr_long
        &'simple_exec',&                                 ! executable
        &1, 1, 0, 0, 1, 5, 1, .false.)                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume subjected to envelope masking', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call automask%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask%set_input('filt_ctrls', 1, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 12.)
        ! mask controls
        call automask%set_input('mask_ctrls', 1, mskdiam)
        call automask%set_input('mask_ctrls', 2, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false., 1.)
        call automask%set_input('mask_ctrls', 3, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for enevloppe mask generation', 'Volume threshold', .false., 0.)
        call automask%set_input('mask_ctrls', 4, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6.)
        call automask%set_input('mask_ctrls', 5, mw)
        ! computer controls
        call automask%set_input('comp_ctrls', 1, nthr)
    end subroutine new_automask

    subroutine new_automask2D
        ! PROGRAM SPECIFICATION
        call automask2D%new(&
        &'automask2D',&                                        ! name
        &'2D envelope masking',&                               ! descr_short
        &'is a program for automated envelope masking in 2D',& ! descr_long
        &'simple_exec',&                                       ! executable
        &1, 1, 0, 0, 2, 4, 1, .false.)                         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask2D%set_input('img_ios', 1, stk)
        automask2D%img_ios(1)%required = .true.
        ! parameter input/output
        call automask2D%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask2D%set_input('filt_ctrls', 1, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms{20 A}', 'low-pass limit in Angstroms{20 A}', .false., 20.)
        call automask2D%set_input('filt_ctrls', 2, 'winsz', 'num', 'Window size for median filter',&
        &'Window size for median filter(in pixels)', 'winsz in pixels', .false., 5.0)
        ! mask controls
        call automask2D%set_input('mask_ctrls', 1, mskdiam)
        call automask2D%set_input('mask_ctrls', 2, 'ngrow', 'num', '# layers to grow',&
        &'Binary layers grown for molecular envelope in pixels{3}', 'width of binary layers grown in pixels{3}', .false., 3.)
        call automask2D%set_input('mask_ctrls', 3, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6.)
        call automask2D%set_input('mask_ctrls', 4, 'positive', 'binary', 'Consider only positive pixels',&
        &'Consider only positive pixels for threshold determination(yes|no){no}', 'only positive(yes|no){no}', .false., 'no')
        ! computer controls
        call automask2D%set_input('comp_ctrls', 1, nthr)
    end subroutine new_automask2D

    subroutine new_autorefine3D_nano
        ! PROGRAM SPECIFICATION
        call autorefine3D_nano%new(&
        &'autorefine3D_nano',&                                                            ! name
        &'auto 3D refinement of metallic nanoparticles',&                                 ! descr_short
        &'is a distributed workflow for automated 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! descr_long
        &'single_exec',&                                                                  ! executable
        &1, 2, 0, 5, 3, 1, 2, .true.)                                                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call autorefine3D_nano%set_input('img_ios', 1, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call autorefine3D_nano%set_input('parm_ios', 1, smpd)
        call autorefine3D_nano%set_input('parm_ios', 2, element)
        ! alternative inputs
        ! <empty>
        ! search controls
        call autorefine3D_nano%set_input('srch_ctrls', 1, nspace)
        call autorefine3D_nano%set_input('srch_ctrls', 2, trs)
        call autorefine3D_nano%set_input('srch_ctrls', 3, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call autorefine3D_nano%set_input('srch_ctrls', 4, 'maxits', 'num', 'Max iterations', 'Maximum number of iterations', 'Max # iterations{5}', .false., 5.)
        call autorefine3D_nano%set_input('srch_ctrls', 5, pgrp)
        ! filter controls
        call autorefine3D_nano%set_input('filt_ctrls', 1, hp)
        call autorefine3D_nano%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5.)
        call autorefine3D_nano%set_input('filt_ctrls', 3, 'lp', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms{1.5}', .true., 1.5)
        ! mask controls
        call autorefine3D_nano%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call autorefine3D_nano%set_input('comp_ctrls', 1, nthr)
        call autorefine3D_nano%set_input('comp_ctrls', 2, script)
    end subroutine new_autorefine3D_nano

    subroutine new_binarize
        call binarize%new(&
        &'binarize',&                                     ! name
        &'Binarization routines for volumes and stacks',& ! descr_long
        &'Binarization routines for volumes and stacks',& ! descr_long
        &'simple_exec',&                                  ! executable
        &0, 3, 2, 0, 0, 0, 1, .false.)                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call binarize%set_input('parm_ios', 1, 'fill_holes', 'binary', 'Fill holes', 'Fill holes(yes|no){no}', '(yes|no){no}', .false., 'no')
        call binarize%set_input('parm_ios', 2, 'ndev', 'num', 'Binarization threshold', 'Binarization threshold in # sigmas', '# sigmas', .false., 0.)
        call binarize%set_input('parm_ios', 3, 'winsz', 'num', 'Half-window size', 'Half-window size(in pixels)', 'winsz in pixels', .false., 15.0)
        ! alternative inputs
        call binarize%set_input('alt_ios', 1, 'vol1', 'file', 'Volume', 'Volume to binarize',&
        & 'input volume e.g. vol.mrc', .false., '')
        call binarize%set_input('alt_ios', 2, 'stk', 'file', 'Stack', 'Stack to binarize',&
        & 'input stack e.g. imgs.mrcs', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        call binarize%set_input('comp_ctrls', 1, nthr)
    end subroutine new_binarize

    subroutine new_calc_pspec
        ! PROGRAM SPECIFICATION
        call calc_pspec%new(&
        &'calc_pspec',&                                                     ! name
        &'Calculate individual particles power spectra; internal use oly',& ! descr_long
        &'Calculate individual particles power spectra; internal use oly',& ! descr_long
        &'simple_exec',&                                                    ! executable
        &0, 0, 0, 0, 0, 0, 2, .true.)                                       ! # entries in each group, requires sp_project
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
        call calc_pspec%set_input('comp_ctrls', 1, nparts)
        call calc_pspec%set_input('comp_ctrls', 2, nthr)
    end subroutine new_calc_pspec

    subroutine new_ced_2D_filter
        ! PROGRAM SPECIFICATION
        call ced_2D_filter%new(&
        &'ced_2D_filter',&             ! name
        &'Coherence-Enhancing Diffution 2D Filter',&             ! descr_short
        &'is a filter that smoothes images in a 2D stack using the diffusion equation.',& ! descr_long
        &'simple_exec',&               ! executable
        &1, 1, 0, 0, 1, 1, 1, .false.)                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ced_2D_filter%set_input('img_ios', 1, 'stk',  'file', 'Stack',  'Stack',  'stack.mrc file', .true., '')
        ! parameter input/output
        call ced_2D_filter%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        !sigma.required = .true.
        call ced_2D_filter%set_input('filt_ctrls',1, 'sigma', 'num', 'Sigma, for Gaussian generation', 'Sigma, for Gaussian generation(in pixels)', &
        & '{1.}', .false., 1.0)
        ! mask controls
        call ced_2D_filter%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call ced_2D_filter%set_input('comp_ctrls', 1, nthr)
    end subroutine new_ced_2D_filter

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

    subroutine new_cleanup2D
        ! PROGRAM SPECIFICATION
        call cleanup2D%new(&
        &'cleanup2D',&                                                          ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm&
        & suitable for the first pass of cleanup after picking',&               ! descr_long
        &'simple_exec',&                                                        ! executable
        &0, 0, 0, 6, 5, 1, 2, .true.)                                           ! # entries in each group, requires sp_project
        cleanup2D%gui_submenu_list = "search,mask,filter"
        cleanup2D%advanced = .false.
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call cleanup2D%set_input('srch_ctrls', 1, ncls)
        call cleanup2D%set_gui_params('srch_ctrls', 1, submenu="search", advanced=.false.)
        call cleanup2D%set_input('srch_ctrls', 2, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cleanup2D%set_gui_params('srch_ctrls', 2, submenu="search")
        call cleanup2D%set_input('srch_ctrls', 3, maxits)
        call cleanup2D%set_gui_params('srch_ctrls', 3, submenu="search")
        call cleanup2D%set_input('srch_ctrls', 4, update_frac)
        call cleanup2D%set_gui_params('srch_ctrls', 4, submenu="search")
        call cleanup2D%set_input('srch_ctrls', 5, objfun)
        call cleanup2D%set_gui_params('srch_ctrls', 5, submenu="search")
        call cleanup2D%set_input('srch_ctrls', 6, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated execution(yes|no){yes}','(yes|no){yes}', .false., 'yes')
        call cleanup2D%set_gui_params('srch_ctrls', 6, submenu="search")
        ! filter controls
        call cleanup2D%set_input('filt_ctrls', 1, hp)
        call cleanup2D%set_gui_params('filt_ctrls', 1, submenu="filter")
        call cleanup2D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cleanup2D%set_gui_params('filt_ctrls', 2, submenu="filter")
        call cleanup2D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 15.)
        call cleanup2D%set_gui_params('filt_ctrls', 3, submenu="filter")
        call cleanup2D%set_input('filt_ctrls', 4, wiener)
        call cleanup2D%set_gui_params('filt_ctrls', 4, submenu="filter")
        call cleanup2D%set_input('filt_ctrls', 5, ml_reg)
        cleanup2D%filt_ctrls(5)%descr_long        = 'Regularization (ML-style) based on the signal power(yes|no){no}'
        cleanup2D%filt_ctrls(5)%descr_placeholder = '(yes|no){no}'
        cleanup2D%filt_ctrls(5)%cval_default      = 'no'
        call cleanup2D%set_gui_params('filt_ctrls', 5, submenu="filter")
        ! mask controls
        call cleanup2D%set_input('mask_ctrls', 1, mskdiam)
        call cleanup2D%set_gui_params('mask_ctrls', 1, submenu="mask", advanced=.false.)
        cleanup2D%mask_ctrls(1)%required = .false.
        ! computer controls
        call cleanup2D%set_input('comp_ctrls', 1, nparts)
        cleanup2D%comp_ctrls(1)%required = .false.
        call cleanup2D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cleanup2D

    subroutine new_center2D_nano
        ! PROGRAM SPECIFICATION
        call center2D_nano%new(&
        &'center2D_nano',&                                                      ! name
        &'Simultaneous 2D alignment and clustering of nanoparticle images',&    ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm&
        & suitable for the first pass of cleanup after time-series tracking',&  ! descr_long
        &'single_exec',&                                                        ! executable
        &0, 0, 0, 4, 3, 1, 3, .true.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call center2D_nano%set_input('srch_ctrls', 1, ncls)
        center2D_nano%srch_ctrls(1)%required = .false.
        call center2D_nano%set_input('srch_ctrls', 2, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call center2D_nano%set_input('srch_ctrls', 3, maxits)
        call center2D_nano%set_input('srch_ctrls', 4, trs)
        ! filter controls
        call center2D_nano%set_input('filt_ctrls', 1, hp)
        call center2D_nano%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5.)
        call center2D_nano%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms{1.0}', .false., 1.)
        ! mask controls
        call center2D_nano%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call center2D_nano%set_input('comp_ctrls', 1, nparts)
        center2D_nano%comp_ctrls(1)%required = .false.
        call center2D_nano%set_input('comp_ctrls', 2, nthr)
        call center2D_nano%set_input('comp_ctrls', 3, script)
    end subroutine new_center2D_nano

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',&                                                          ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',&                 ! descr_long
        &'simple_exec',&                                                        ! executable
        &1, 0, 0, 10, 7, 1, 2, .true.)                                          ! # entries in each group, requires sp_project
        cluster2D%gui_submenu_list = "search,mask,filter"
        cluster2D%advanced = .false.
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster2D%set_input('img_ios', 1, 'refs', 'file', 'Initial references',&
        &'Initial 2D references used to bootstrap the search', 'xxx.mrc file with references', .false., 'refs.mrc')
        call cluster2D%set_gui_params('img_ios', 1, submenu="search")
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D%set_input('srch_ctrls', 1, ncls)
        call cluster2D%set_gui_params('srch_ctrls', 1, submenu="search", advanced=.false.)
        call cluster2D%set_input('srch_ctrls', 2, startit)
        call cluster2D%set_gui_params('srch_ctrls', 2, submenu="search")
        call cluster2D%set_input('srch_ctrls', 3, trs)
        call cluster2D%set_gui_params('srch_ctrls', 3, submenu="search")
        call cluster2D%set_input('srch_ctrls', 4, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 'yes')
        call cluster2D%set_gui_params('srch_ctrls', 4, submenu="search")
        call cluster2D%set_input('srch_ctrls', 5, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D%set_gui_params('srch_ctrls', 5, submenu="search")
        call cluster2D%set_input('srch_ctrls', 6, maxits)
        call cluster2D%set_gui_params('srch_ctrls', 6, submenu="search")
        call cluster2D%set_input('srch_ctrls', 7, update_frac)
        call cluster2D%set_gui_params('srch_ctrls', 7, submenu="search")
        call cluster2D%set_input('srch_ctrls', 8, objfun)
        call cluster2D%set_gui_params('srch_ctrls', 8, submenu="search")
        call cluster2D%set_input('srch_ctrls', 9, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc|greedy){snhc}', '(snhc|greedy){snhc}', .false., 'snhc')
        call cluster2D%set_gui_params('srch_ctrls', 9, submenu="search")
        call cluster2D%set_input('srch_ctrls', 10, sigma_est)
        call cluster2D%set_gui_params('srch_ctrls', 10, submenu="search")
        ! filter controls
        call cluster2D%set_input('filt_ctrls', 1, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D%set_gui_params('filt_ctrls', 1, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 2, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false., 20.)
        call cluster2D%set_gui_params('filt_ctrls', 2, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 3, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &few iterations of search, before the automatic scheme kicks in. Also controls the degree of downsampling in the first &
        &phase', 'initial low-pass limit in Angstroms', .false., 15.)
        call cluster2D%set_gui_params('filt_ctrls', 3, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.)
        call cluster2D%set_gui_params('filt_ctrls', 4, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 5,  wiener)
        call cluster2D%set_gui_params('filt_ctrls', 5, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 6,  graphene_filt)
        call cluster2D%set_gui_params('filt_ctrls', 6, submenu="filter")
        call cluster2D%set_input('filt_ctrls', 7, ml_reg)
        cluster2D%filt_ctrls(7)%descr_long        = 'Regularization (ML-style) based on the signal power(yes|no){no}'
        cluster2D%filt_ctrls(7)%descr_placeholder = '(yes|no){no}'
        cluster2D%filt_ctrls(7)%cval_default      = 'no'
        call cluster2D%set_gui_params('filt_ctrls', 7, submenu="filter")
        ! mask controls
        call cluster2D%set_input('mask_ctrls', 1, mskdiam)
        call cluster2D%set_gui_params('mask_ctrls', 1, submenu="mask")
        ! computer controls
        call cluster2D%set_input('comp_ctrls', 1, nparts)
        cluster2D%comp_ctrls(1)%required = .false.
        call cluster2D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster2D

    subroutine new_cluster2D_nano
        ! PROGRAM SPECIFICATION
        call cluster2D_nano%new(&
        &'cluster2D_nano',&                                                                 ! name
        &'Simultaneous 2D alignment and clustering of time-series of nanoparticle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm for time-series of nanoparticle images',& ! descr_long
        &'single_exec',&                                                                    ! executable
        &0, 1, 0, 5, 5, 1, 3, .true.)                                                       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster2D_nano%set_input('parm_ios', 1, moldiam)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_nano%set_input('srch_ctrls', 1, nptcls_per_cls)
        call cluster2D_nano%set_input('srch_ctrls', 2, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D_nano%set_input('srch_ctrls', 3, 'winsz', 'num', 'Half-window size', 'Half-window size(frames)', 'winsz in # frames', .false., 3.0)
        call cluster2D_nano%set_input('srch_ctrls', 4, maxits)
        call cluster2D_nano%set_input('srch_ctrls', 5, trs)
        ! filter controls
        call cluster2D_nano%set_input('filt_ctrls', 1, hp)
        call cluster2D_nano%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{5.0}', .false., 5.)
        call cluster2D_nano%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit{1.0}', 'low-pass limit in Angstroms', .false., 1.)
        call cluster2D_nano%set_input('filt_ctrls', 4, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'initial low-pass limit in Angstroms', .false., 1.)
        call cluster2D_nano%set_input('filt_ctrls', 5, 'lpstop', 'num', 'Final low-pass limit', 'Final low-pass limit{1.0}', 'final low-pass limit in Angstroms', .false., 1.)
        ! mask controls
        call cluster2D_nano%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call cluster2D_nano%set_input('comp_ctrls', 1, nparts)
        cluster2D_nano%comp_ctrls(1)%required = .false.
        call cluster2D_nano%set_input('comp_ctrls', 2, nthr)
        call cluster2D_nano%set_input('comp_ctrls', 3, script)
    end subroutine new_cluster2D_nano

    subroutine new_cluster2D_stream
        ! PROGRAM SPECIFICATION
        call cluster2D_stream%new(&
        &'cluster2D_stream',&                                                                     ! name
        &'Simultaneous 2D alignment and clustering of single-particle images in streaming mode',& ! descr_short
        &'is a distributed workflow implementing cluster2D in streaming mode',&                   ! descr_long
        &'simple_exec',&                                                                          ! executable
        &0, 1, 0, 7, 6, 1, 5, .true.)                                                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster2D_stream%set_input('parm_ios', 1, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preprocess_stream', .true., '')
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_stream%set_input('srch_ctrls', 1, ncls_start)
        call cluster2D_stream%set_input('srch_ctrls', 2, nptcls_per_cls)
        cluster2D_stream%srch_ctrls(2)%required = .true.
        cluster2D_stream%srch_ctrls(2)%rval_default = 300.
        call cluster2D_stream%set_input('srch_ctrls', 3, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 'yes')
        call cluster2D_stream%set_input('srch_ctrls', 4, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
            &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D_stream%set_input('srch_ctrls', 5, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of groups to sort the particles into prior to averaging to create 2D class averages with improved SNR', 'Maximum # 2D clusters', .true., 200.)
        call cluster2D_stream%set_input('srch_ctrls', 6, lpthres)
        call cluster2D_stream%set_input('srch_ctrls', 7, 'refine', 'multi', 'Refinement mode', '2D Refinement mode(no|greedy){no}', '(no|greedy){no}', .false., 'no')
        ! filter controls
        call cluster2D_stream%set_input('filt_ctrls', 1, hp)
        call cluster2D_stream%set_input('filt_ctrls', 2, 'cenlp',      'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D_stream%set_input('filt_ctrls', 3, 'lp',         'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 15.)
        call cluster2D_stream%set_input('filt_ctrls', 4, 'lpstop',     'num', 'Resolution limit', 'Resolution limit', 'Resolution limit in Angstroms', .false., 2.0*MAX_SMPD)
        call cluster2D_stream%set_input('filt_ctrls', 5,  wiener)
        call cluster2D_stream%set_input('filt_ctrls', 6, lplim_crit)
        ! mask controls
        call cluster2D_stream%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call cluster2D_stream%set_input('comp_ctrls', 1, nchunks)
        call cluster2D_stream%set_input('comp_ctrls', 2, nparts_chunk)
        call cluster2D_stream%set_input('comp_ctrls', 3, nparts)
        call cluster2D_stream%set_input('comp_ctrls', 4, nthr)
        call cluster2D_stream%set_input('comp_ctrls', 5, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740.)
    end subroutine new_cluster2D_stream

    subroutine new_cluster2D_subsets
        ! PROGRAM SPECIFICATION
        call cluster2D_subsets%new(&
        &'cluster2D_subsets',&                                                                     ! name
        &'Simultaneous 2D alignment and clustering of single-particle images in streaming mode',& ! descr_short
        &'is a distributed workflow implementing cluster2D in streaming mode',&                   ! descr_long
        &'simple_exec',&                                                                          ! executable
        &0, 0, 0, 7, 5, 1, 5, .true.)                                                             ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_subsets%set_input('srch_ctrls', 1, ncls_start)
        call cluster2D_subsets%set_input('srch_ctrls', 2, nptcls_per_cls)
        cluster2D_subsets%srch_ctrls(2)%required = .true.
        cluster2D_subsets%srch_ctrls(2)%rval_default = 300.
        call cluster2D_subsets%set_input('srch_ctrls', 3, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 'yes')
        call cluster2D_subsets%set_input('srch_ctrls', 4, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
            &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D_subsets%set_input('srch_ctrls', 5, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .true., 200.)
        call cluster2D_subsets%set_input('srch_ctrls', 6, lpthres)
        call cluster2D_subsets%set_input('srch_ctrls', 7, 'refine', 'multi', 'Refinement mode', '2D Refinement mode(no|greedy){no}', '(no|greedy){no}', .false., 'no')
        ! filter controls
        call cluster2D_subsets%set_input('filt_ctrls', 1, hp)
        call cluster2D_subsets%set_input('filt_ctrls', 2, 'cenlp',      'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster2D_subsets%set_input('filt_ctrls', 3, 'lp',         'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 15.)
        call cluster2D_subsets%set_input('filt_ctrls', 4, 'lpstop',     'num', 'Resolution limit', 'Resolution limit', 'Resolution limit in Angstroms', .false., 2.0*MAX_SMPD)
        call cluster2D_subsets%set_input('filt_ctrls', 5,  wiener)
        ! mask controls
        call cluster2D_subsets%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call cluster2D_subsets%set_input('comp_ctrls', 1, nchunks)
        call cluster2D_subsets%set_input('comp_ctrls', 2, nparts_chunk)
        call cluster2D_subsets%set_input('comp_ctrls', 3, nparts_pool)
        cluster2D_subsets%comp_ctrls(3)%required = .true.
        call cluster2D_subsets%set_input('comp_ctrls', 4, nthr)
        call cluster2D_subsets%set_input('comp_ctrls', 5, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740.)
    end subroutine new_cluster2D_subsets

    subroutine new_cluster3D
        ! PROGRAM SPECIFICATION
        call cluster3D%new(&
        &'cluster3D',&                                                             ! name
        &'3D heterogeneity analysis',&                                             ! descr_short
        &'is a distributed workflow for heterogeneity analysis by 3D clustering',& ! descr_long
        &'simple_exec',&                                                     ! executable
        &0, 1, 0, 7, 6, 3, 2, .true.)                                              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster3D%set_input('parm_ios', 1, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to separate',&
        '# states to separate', .true., 2.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster3D%set_input('srch_ctrls', 1, nspace)
        call cluster3D%set_input('srch_ctrls', 2, startit)
        call cluster3D%set_input('srch_ctrls', 3, maxits)
        call cluster3D%set_input('srch_ctrls', 4, frac)
        call cluster3D%set_input('srch_ctrls', 5, pgrp)
        call cluster3D%set_input('srch_ctrls', 6, objfun)
        call cluster3D%set_input('srch_ctrls', 7, 'refine', 'multi', 'Refinement mode', 'Refinement mode(cluster|clustersym)&
        &){cluster}', '(cluster|clustersym){cluster}', .false., 'cluster')
        ! filter controls
        call cluster3D%set_input('filt_ctrls', 1, hp)
        call cluster3D%set_input('filt_ctrls', 2, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call cluster3D%set_input('filt_ctrls', 3, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call cluster3D%set_input('filt_ctrls', 4, lplim_crit)
        call cluster3D%set_input('filt_ctrls', 5, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass resolution limit','low-pass limit in Angstroms', .false., 0.)
        call cluster3D%set_input('filt_ctrls', 6, envfsc)
        ! mask controls
        call cluster3D%set_input('mask_ctrls', 1, mskdiam)
        call cluster3D%set_input('mask_ctrls', 2, mskfile)
        call cluster3D%set_input('mask_ctrls', 3, focusmskdiam)
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
        &'simple_exec',&                                                     ! executable
        &2, 1, 0, 9, 6, 1, 2, .true.)                                        ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster3D_refine%set_input('img_ios', 1, 'msklist', 'file', 'List of mask files', 'List (.txt file) of mask files for the different states', 'e.g. mskfiles.txt', .false., '')
        call cluster3D_refine%set_input('img_ios', 2, 'vollist', 'file', 'List of reference volumes files', 'List (.txt file) of reference volumes for the different states', 'e.g. refvols.txt', .false., '')
        ! parameter input/output
        call cluster3D_refine%set_input('parm_ios', 1,  'state', 'num', 'State to refine', 'Index of state to refine', 'give state index', .false., 1.)
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
        call cluster3D_refine%set_input('srch_ctrls', 9, objfun)
        ! filter controls
        call cluster3D_refine%set_input('filt_ctrls', 1, hp)
        call cluster3D_refine%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call cluster3D_refine%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call cluster3D_refine%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call cluster3D_refine%set_input('filt_ctrls', 5, lplim_crit)
        call cluster3D_refine%set_input('filt_ctrls', 6, envfsc)
        ! mask controls
        call cluster3D_refine%set_input('mask_ctrls', 1, mskdiam)
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
        &to apply a balancing constraint (on the class population). Adjust balance until you are &
        &satisfied with the shape of the histogram',&                              ! descr_long
        &'simple_exec',&                                                           ! executable
        &0, 2, 0, 0, 3, 1, 1, .true.)                                              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster_cavgs%set_input('parm_ios', 1, 'bin_cls', 'multi', 'Perform good/bad classification based on common lines',&
        &'Classes with lower common line correlation to the rest are rejected by Otsu(yes|no|only){yes}', '(yes|no|only){yes}', .false., 'yes')
        call cluster_cavgs%set_input('parm_ios', 2,  'ncls', 'num', 'Number of highest populated AP clusters', 'Number of highest populated AP clusters to map the AP solution onto', '# classes', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call cluster_cavgs%set_input('filt_ctrls', 1, hp)
        call cluster_cavgs%set_input('filt_ctrls', 2, lp)
        cluster_cavgs%filt_ctrls(2)%required = .true.
        call cluster_cavgs%set_input('filt_ctrls', 3, lpthres)
        ! mask controls
        call cluster_cavgs%set_input('mask_ctrls', 1, mskdiam)
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

    subroutine new_comparemc
        ! PROGRAM SPECIFICATION
        call comparemc%new(&
        &'comparemc',&                                                    ! name
        &'utility',&                     ! descr_short
        &'is a utility program',& ! descr_long
        &'simple_exec',&                                                ! executable
        &0, 1, 0, 0, 0, 0, 0, .false.)                                  ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call comparemc%set_input('parm_ios', 1, 'infile', 'file', 'List of project files', 'List of project files', 'e.g. projects.txt', .true., '')
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
    end subroutine new_comparemc

    subroutine new_ctf_estimate
        ! PROGRAM SPECIFICATION
        call ctf_estimate%new(&
        &'ctf_estimate', &                                              ! name
        &'CTF parameter fitting',&                                      ! descr_short
        &'is a distributed SIMPLE workflow for CTF parameter fitting',& ! descr_long
        &'simple_exec',&                                                ! executable
        &0, 2, 0, 3, 2, 0, 2, .true.)                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_estimate%set_input('parm_ios', 1, pspecsz)
        call ctf_estimate%set_input('parm_ios', 2, ctfpatch)
        ! alternative inputs
        ! <empty>
        ! search controls
        call ctf_estimate%set_input('srch_ctrls', 1, dfmin)
        call ctf_estimate%set_input('srch_ctrls', 2, dfmax)
        call ctf_estimate%set_input('srch_ctrls', 3, astigtol)
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
        call ctfops%set_input('parm_ios', 2, 'neg', 'binary', 'Invert contrast','Invert contrast(yes|no){no}',&
            '(yes|no){no}', .false., 'no')
        call ctfops%set_input('parm_ios', 3, oritab)
        call ctfops%set_input('parm_ios', 4, deftab)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ctfops%set_input('filt_ctrls', 1, ctf)
        call ctfops%set_input('filt_ctrls', 2, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', &
            'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        ! mask controls
        ! <empty>
        ! computer controls
        call ctfops%set_input('comp_ctrls', 1, nthr)
    end subroutine new_ctfops

    subroutine new_detect_atoms
        ! PROGRAM SPECIFICATION
        call detect_atoms%new(&
        &'detect_atoms', &                                      ! name
        &'Detect atoms in atomic-resolution nanoparticle map',& ! descr_short
        &'is a program for identifying atoms in atomic-resolution nanoparticle maps and generating bin and connected-comp map',& ! descr long
        &'single_exec',&                                        ! executable
        &2, 4, 0, 0, 1, 1, 1, .false.)                         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call detect_atoms%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Nanoparticle volume to analyse', &
        & 'input volume e.g. vol.mrc', .true., '')
        call detect_atoms%set_input('img_ios', 2, 'vol2', 'file', 'Volume', 'Nanoparticle volume to use for lattice fitting', &
        & 'input volume 4 lattice fit e.g. vol2.mrc', .false., '')
        ! parameter input/output
        call detect_atoms%set_input('parm_ios', 1, smpd)
        call detect_atoms%set_input('parm_ios', 2, 'corr_thres', 'num', 'Per-atom corr threshold', 'Per-atom validation correlation threshold for discarding atoms(0.3-0.5){0.5}', &
        & 'Corr threshold for discarding atoms(0.3-0.5){0.5}', .false., 0.5)
        call detect_atoms%set_input('parm_ios', 3, 'use_thres', 'binary', 'Use contact-based thresholding', 'Use contact-based thresholding(yes|no){yes}', &
        & '(yes|no){yes}', .false., 'yes')
        call detect_atoms%set_input('parm_ios', 4, 'cs_thres', 'num', 'Contact score thresholding', 'Use contact score threshold during autorefine3D_nano(2-3){2}', '(2-3){2}', .false., 2.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call detect_atoms%set_input('filt_ctrls', 1, element)
        ! mask controls
        call detect_atoms%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call detect_atoms%set_input('comp_ctrls', 1, nthr)
    end subroutine new_detect_atoms

    subroutine new_dock_volpair
        ! PROGRAM SPECIFICATION
        call dock_volpair%new(&
        &'dock_volpair', &                              ! name
        &'Dock a pair of volumes',&                     ! descr_short
        &'is a program for docking a pair of volumes',& ! descr long
        &'simple_exec',&                                ! executable
        &3, 2, 0, 2, 3, 1, 1, .false.)                  ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call dock_volpair%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Reference volume', &
        & 'input reference volume e.g. vol1.mrc', .true., '')
        call dock_volpair%set_input('img_ios', 2, 'vol2', 'file', 'Volume', 'Target volume', &
        & 'input target volume e.g. vol2.mrc', .true., '')
        call dock_volpair%set_input('img_ios', 3, outvol)
        ! parameter input/output
        call dock_volpair%set_input('parm_ios', 1, smpd)
        call dock_volpair%set_input('parm_ios', 2,  outfile)
        ! alternative inputs
        ! <empty>
        ! search controls
        call dock_volpair%set_input('srch_ctrls', 1, trs)
        dock_volpair%srch_ctrls(1)%rval_default = 5.
        call dock_volpair%set_input('srch_ctrls', 2, 'dockmode', 'multi', 'Docking mode', 'Docking mode(rot|shift|rotshift|refine){rotshift}', '(rot|shift|rotshift|refine){rotshift}', .false., 'rotshift')
        ! filter controls
        call dock_volpair%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass resolution limit', 'low-pass limit in Angstroms', .true., 0.)
        call dock_volpair%set_input('filt_ctrls', 2, 'lpstop',   'num', 'Final low-pass limit',   'Final low-pass resolution limit',   'low-pass limit in Angstroms', .true., 0.)
        call dock_volpair%set_input('filt_ctrls', 3, hp)
        ! mask controls
        call dock_volpair%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call dock_volpair%set_input('comp_ctrls', 1, nthr)
    end subroutine new_dock_volpair

    subroutine new_estimate_diam
        ! PROGRAM SPECIFICATION
        call estimate_diam%new(&
        &'estimate_diam',&                                                                                            ! name
        &'Estimation of a suitable mask diameter for nanoparticle time-series',&                                        ! descr_short
        &'is a program for estimation of a suitable mask diameter for spherical masking of nanoparticle time-series ',& ! descr_long
        &'single_exec',&                                                                                              ! executable
        &1, 2, 0, 0, 1, 1, 1, .false.)                                               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call estimate_diam%set_input('img_ios', 1, stk)
        estimate_diam%img_ios(1)%required = .true.
        ! parameter input/output
        call estimate_diam%set_input('parm_ios', 1, smpd)
        call estimate_diam%set_input('parm_ios', 2, 'roavg', 'binary', 'Rotationally average', 'Rotationally average before diam estimate(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call estimate_diam%set_input('filt_ctrls', 1, lp)
        estimate_diam%filt_ctrls(1)%descr_short  = 'low-pass limit in Angstroms{7.}'
        estimate_diam%filt_ctrls(1)%rval_default = 7.
        ! mask controls
        call estimate_diam%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call estimate_diam%set_input('comp_ctrls', 1, nthr)
    end subroutine new_estimate_diam

    subroutine new_extract
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', &                                                           ! name
        &'Extract particle images from integrated movies',&                     ! descr_short
        &'is a program for extracting particle images from integrated movies',& ! descr long
        &'simple_exec',&                                                  ! executable
        &1, 4, 0, 0, 0, 0, 1, .true.)                                           ! # entries in each group, requires sp_project
        extract%gui_submenu_list = "extract"
        extract%advanced = .false.
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call extract%set_input('img_ios', 1, 'dir_box', 'file', 'Box files directory', 'Directory to read the box files from', 'e.g. boxes/', .false., '')
        call extract%set_gui_params('img_ios', 1, submenu="extract")
        ! parameter input/output
        call extract%set_input('parm_ios', 1, box)
        call extract%set_gui_params('parm_ios', 1, submenu="extract", advanced=.false.)
        extract%parm_ios(1)%required = .false.
        call extract%set_input('parm_ios', 2, pcontrast)
        call extract%set_gui_params('parm_ios', 2, submenu="extract", advanced=.false.)
        call extract%set_input('parm_ios', 3, 'outside', 'binary', 'Extract outside boundaries', 'Extract boxes outside the micrograph boundaries(yes|no){no}', '(yes|no){no}', .false., 'no')
        call extract%set_gui_params('parm_ios', 3, submenu="extract")
        call extract%set_input('parm_ios', 4, 'ctf', 'multi', 'Whether to extract particles with phases flipped', 'Whether to extract particles with phases flipped(flip|no){no}', '(flip|no){no}', .false., 'no')
        call extract%set_gui_params('parm_ios', 4, submenu="extract")
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call extract%set_input('comp_ctrls', 1, nparts)
    end subroutine new_extract

    subroutine new_export_starproject
        ! PROGRAM SPECIFICATION
        call export_starproject%new(&
        &'export_starproject', &                                                ! name
        &'Export projectfile in star format',&                                  ! descr_short
        &'is a program to export a SIMPLE projectfile in star format',&         ! descr long
        &'simple_exec',&                                                        ! executable
        &0, 0, 0, 0, 0, 0, 0, .true.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
       ! call export_starproject%set_input('parm_ios', 1, projfile)
        ! parameter input/output
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
    end subroutine new_export_starproject

    subroutine new_filter
        ! PROGRAM SPECIFICATION
        call filter%new(&
        &'filter',&                                   ! name
        &'Filter stack/volume',&                      ! descr_short
        &'is a program for filtering stack/volume',&  ! descr_long
        &'simple_exec',&                              ! executable
        &2, 1, 2, 0, 12, 0, 1, .false.)               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call filter%set_input('img_ios', 1, outstk)
        call filter%set_input('img_ios', 2, outvol)
        ! parameter input/output
        call filter%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        call filter%set_input('alt_ios', 1, 'stk',  'file', 'Stack to filter',  'Stack of images to filter', 'e.g. stk.mrcs',     .false., '')
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
        call filter%set_input('filt_ctrls', 7, 'real_filter', 'multi', 'Real-space filter',&
        &'Real-space filter(median|average|stdev|bman|NLmean|no){no}', '(median|average|stdev|bman|NLmean|no){no}', .false., 'no')
        call filter%set_input('filt_ctrls', 8, 'fsc', 'file', 'FSC file', 'FSC file',          'e.g. fsc_state01.bin file', .false., '')
        call filter%set_input('filt_ctrls', 9, frcs)
        call filter%set_input('filt_ctrls',10, 'filter', 'multi', 'Filter type(tv|nlmean|no){no}', 'Filter type(tv|nlmean|corr|no){no}', '(tv|nlmean|no){no}', .false., 'no')
        call filter%set_input('filt_ctrls',11, 'lambda', 'num', 'Tv filter lambda', 'Strength of noise reduction', '(0.5-3.0){1.0}', .false., 1.0)
        call filter%set_input('filt_ctrls',12, 'sigma', 'num', 'sigma, for Gaussian generation', 'sigma, for Gaussian generation(in pixels)', &
        & '{1.}', .false., 1.0)
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
        &2, 1, 0, 0, 2, 2, 1, .false.)                                          ! # entries in each group, requires sp_project
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
        hp%required = .false.
        lp%required = .false.
        call fsc%set_input('filt_ctrls', 1, hp)
        call fsc%set_input('filt_ctrls', 2, lp)
        ! mask controls
        call fsc%set_input('mask_ctrls', 1, mskdiam)
        call fsc%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call fsc%set_input('comp_ctrls', 1, nthr)
    end subroutine new_fsc

    subroutine new_gen_pspecs_and_thumbs
        ! PROGRAM SPECIFICATION
        call gen_pspecs_and_thumbs%new(&
        &'gen_pspecs_and_thumbs', &                                              ! name
        &'Motion correction of movies',&                                         ! descr_short
        &'is a distributed workflow for generating power spectra and thumbnails&
        & for imported integrated movies',&                                      ! descr_long
        &'simple_exec',&                                                   ! executable
        &0, 1, 0, 0, 0, 0, 2, .true.)                                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call gen_pspecs_and_thumbs%set_input('parm_ios', 1, pspecsz)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call gen_pspecs_and_thumbs%set_input('comp_ctrls', 1, nparts)
        call gen_pspecs_and_thumbs%set_input('comp_ctrls', 2, nthr)
    end subroutine new_gen_pspecs_and_thumbs

    subroutine new_import_starproject
        ! PROGRAM SPECIFICATION
        call import_starproject%new(&
        &'import_starproject', &                                                ! name
        &'Import project in in star format',&                                   ! descr_short
        &'is a program to import a SIMPLE projectfile from star format',&       ! descr long
        &'simple_exec',&                                                        ! executable
        &0, 3, 0, 0, 0, 0, 0, .false.)                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        import_starproject%gui_submenu_list = "data"
        import_starproject%advanced = .false.
        ! image input/output
        ! parameter input/output
        call import_starproject%set_input('parm_ios', 1, 'import_dir',  'dir',   'Import directory', 'Directory to import data from. In auto mode should be output from an external job e.g. relion', 'e.g. MotionCorr/job001', .true., '')
        call import_starproject%set_gui_params('parm_ios', 1, submenu="data", advanced=.false.)
        call import_starproject%set_input('parm_ios', 2, 'starfile',    'file',  'Metadata starfile', 'Path to starfile containing micrograph/particle metadata. Only required when not using auto mode', 'e.g. micrographs.star', .false., '')
        call import_starproject%set_gui_params('parm_ios', 2, submenu="data", advanced=.false.)
        call import_starproject%set_input('parm_ios', 3, 'import_type', 'multi', 'Import type', 'Type of data contained in starfile (auto|mic|ptcl2D|ptcl3D){auto}. Auto mode (default) will attempt to determine this automatically', '(auto|mic|ptcl2D|ptcl3D){auto}', .false., 'auto')
        call import_starproject%set_gui_params('parm_ios', 3, submenu="data", advanced=.false.)
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
    end subroutine new_import_starproject

    subroutine new_info_image
        ! PROGRAM SPECIFICATION
        call info_image%new(&
        &'info_image', &                                                                       ! name
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
        &'info_stktab', &                                                        ! name
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
        &'initial_3Dmodel',&                                                          ! name
        &'3D ab initio model generation from class averages',&                        ! descr_short
        &'is a distributed workflow for generating an initial 3D model from class&
        & averages obtained with cluster2D',&                                         ! descr_long
        &'simple_exec',&                                                              ! executable
        &0, 0, 0, 5, 6, 2, 2, .true.)                                                 ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call initial_3Dmodel%set_input('srch_ctrls', 1, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call initial_3Dmodel%set_input('srch_ctrls', 2, frac)
        call initial_3Dmodel%set_input('srch_ctrls', 3, pgrp)
        call initial_3Dmodel%set_input('srch_ctrls', 4, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Final low-pass limit controls the degree of down-scaling(yes|no){yes}','(yes|no){yes}', .false., 'yes')
        call initial_3Dmodel%set_input('srch_ctrls', 5, 'pgrp_start','str', 'Initial point-group symmetry',&
        &'Initial point-group symmetry(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .false., 'c1')
        ! filter controls
        call initial_3Dmodel%set_input('filt_ctrls', 1, hp)
        call initial_3Dmodel%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call initial_3Dmodel%set_input('filt_ctrls', 3, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass resolution limit for the first stage of ab-initio model generation',&
            &'low-pass limit in Angstroms', .false., 0.)
        call initial_3Dmodel%set_input('filt_ctrls', 4, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the second stage (no e/o cavgs refinement) in Angstroms', .false., 8.)
        call initial_3Dmodel%set_input('filt_ctrls', 5, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
            & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 15.)
        call initial_3Dmodel%set_input('filt_ctrls', 6, nonuniform)
        ! mask controls
        call initial_3Dmodel%set_input('mask_ctrls', 1, mskdiam)
        call initial_3Dmodel%set_input('mask_ctrls', 2, automsk)
        ! computer controls
        call initial_3Dmodel%set_input('comp_ctrls', 1, nparts)
        initial_3Dmodel%comp_ctrls(1)%required = .false.
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
        call import_boxes%set_input('parm_ios', 1, 'boxtab', 'file', 'List of box files', &
            'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .true., '')
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
        &'is a program for importing DDD movies to the project. The movies can be located in any read-only location&
        & accessible to the project. If the movies contain only a single frame, they will be interpreted as motion-corrected&
        & and integrated. Box files (in EMAN format) can be imported along with the movies',&
        &'simple_exec',&                                         ! executable
        &1, 8, 0, 0, 0, 0, 0, .true.)                            ! # entries in each group, requires sp_project
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
        call import_movies%set_input('parm_ios', 5, ctf_yes)
        call import_movies%set_input('parm_ios', 6, phaseplate)
        call import_movies%set_input('parm_ios', 7, 'boxtab', 'file', 'List of box files', 'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .false., '')
        call import_movies%set_input('parm_ios', 8, 'deftab', 'file','Pre-determined per-micrograph CTF parameters',&
        &'List of CTF parmeters for micrographs import only', 'e.g. deftab.txt', .false., '')
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
        &'all',&                                                   ! executable
        &0, 12, 2, 0, 0, 0, 0, .true.)                             ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call import_particles%set_input('parm_ios', 1, smpd)
        call import_particles%set_input('parm_ios', 2, kv)
        import_particles%parm_ios(2)%required = .true.
        call import_particles%set_input('parm_ios', 3, cs)
        import_particles%parm_ios(3)%required = .true.
        call import_particles%set_input('parm_ios', 4, fraca)
        import_particles%parm_ios(4)%required = .true.
        call import_particles%set_input('parm_ios', 5, ctf_yes)
        call import_particles%set_input('parm_ios', 6, phaseplate)
        call import_particles%set_input('parm_ios', 7, oritab)
        call import_particles%set_input('parm_ios', 8, deftab)
        call import_particles%set_input('parm_ios', 9, 'plaintexttab', 'file', 'Plain text file of input parameters',&
        'Plain text file of tabulated per-particle input parameters: dfx, dfy, angast, phshift', 'e.g. params.txt', .false., '')
        call import_particles%set_input('parm_ios', 10, 'dfunit', 'multi', 'Underfocus unit', 'Underfocus unit(A|microns){microns}', '(A|microns){microns}', .false., 'microns')
        call import_particles%set_input('parm_ios', 11, 'angastunit', 'multi', 'Angle of astigmatism unit', 'Angle of astigmatism unit(radians|degrees){degrees}', '(radians|degrees){degrees}', .false., 'degrees')
        call import_particles%set_input('parm_ios', 12, 'phshiftunit', 'multi', 'Phase-shift unit', 'Phase-shift unit(radians|degrees){radians}', '(radians|degrees){radians}', .false., 'degrees')
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

    subroutine new_export_relion
        ! PROGRAM SPECIFICATION
        call export_relion%new(&
        &'export_relion',&                                              ! name
        &'Export project to relion ',&                                  ! descr_short
        &'is a program to export simple project to relion',&
        &'simple_exec',&                                                ! executable
        &0, 5, 0, 0, 0, 0, 0, .true.)                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call export_relion%set_input('parm_ios', 1, 'tiltgroupmax', 'num', 'Max movies in a tilt/optics group', &
            'Sub-divide beamtilt/optics groups', '0', .false., 0.0)
        call export_relion%set_input('parm_ios', 2, 'reliongroups', 'num', 'Number of Relion groups based on defocus', &
            'Divide particles into X groups based on defocus for relion', '# micrographs', .false., 0.0)
        call export_relion%set_input('parm_ios', 3, 'xmlloc', 'file', 'Pathname of EPU XML files',&
            'Pathname of EPU XML files ', 'e.g. /data/xml', .false., 'NONE')
        call export_relion%set_input('parm_ios', 4, 'tilt_thres', 'num', 'Distance threshold',&
            'Distance threshold for hierarchical clustering of beamtilt/shift groups ', '{0.05}', .false., 0.05)
        call export_relion%set_input('parm_ios', 5, 'optics_offset', 'num', 'Offset to apply to optics group numbering', &
            'Offset to apply to optics group numbering', '{0}', .false., 0.0)
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_export_relion

    subroutine new_make_cavgs
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &                           ! name
        &'Make class averages',&                   ! descr_short
        &'is a distributed workflow for generating class averages or initial random references&
        & for cluster2D execution',&               ! descr_long
        &'simple_exec',&                           ! executable
        &1, 3, 0, 0, 1, 0, 2, .true.)              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call make_cavgs%set_input('img_ios', 1, 'refs', 'file', 'Output 2D references',&
        &'Output 2D references', 'xxx.mrc file with references', .false., '')
        ! parameter input/output
        call make_cavgs%set_input('parm_ios', 1, ncls)
        call make_cavgs%set_input('parm_ios', 2, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call make_cavgs%set_input('parm_ios', 3, remap_cls)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call make_cavgs%set_input('filt_ctrls', 1, wiener)
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
        &'Make orientations',&               ! descr_short
        &'is a program for making SIMPLE orientation files. Random Euler angles and random origin shifts are generated.&
        & If ndiscrete is set to an integer number > 0, the orientations produced are randomly sampled from the set of&
        & ndiscrete quasi-even projection directions, and the in-plane parameters are assigned randomly. If even is set&
        & to yes, then all nptcls orientations are assigned quasi-even projection directions and  random in-plane parameters.&
        & If nstates is set to some integer number > 0, then states are assigned randomly',&
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

    subroutine new_map_cavgs_selection
        ! PROGRAM SPECIFICATION
        call map_cavgs_selection%new(&
        &'map_cavgs_selection',&                                         ! name
        &'Map class average selection to particles in project file',&    ! descr_short
        &'is a program for mapping selection based on class averages to the individual particles using correlation matching',& ! descr_long
        &'all',&                                                         ! executable
        &2, 1, 0, 0, 0, 0, 0, .true.)                                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call  map_cavgs_selection%set_input('img_ios', 1, 'stk', 'file', 'Stack of cavgs to select from', 'Stack of cavgs to select from', 'e.g. cavgs_iter0XX.mrc', .false., '')
        call  map_cavgs_selection%set_input('img_ios', 2, 'stk2', 'file', 'Stack of selected cavgs', 'Stack of selected cavgs', 'e.g. selected.spi', .true., '')
        ! parameter input/output
        call  map_cavgs_selection%set_input('parm_ios',1, prune)
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
    end subroutine new_map_cavgs_selection

    subroutine new_map_cavgs_states
        ! PROGRAM SPECIFICATION
        call map_cavgs_states%new(&
        &'map_cavgs_states',&                                            ! name
        &'Map class average state selection by common lines-based clustering to particles in project file',& ! descr_short
        &'is a program for mapping state selection by common lines-based clustering of class averages to the individual particles using correlation matching',& ! descr_long
        &'all',&                                                         ! executable
        &2, 0, 0, 0, 0, 0, 0, .true.)                                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call map_cavgs_states%set_input('img_ios', 1, 'stk',    'file', 'Stack of cavgs to select from', 'Stack of cavgs to select from', 'e.g. cavgs_iter0XX.mrc', .false., '')
        call map_cavgs_states%set_input('img_ios', 2, 'stktab', 'file', 'Stacks of class averages list',&
        &'List of stacks of class averages to use for mapping states', 'list input e.g. stktab.txt', .true., '')
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
    end subroutine new_map_cavgs_states

    subroutine new_mask
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',&                                                        ! name
        &'Mask images/volumes',&                                         ! descr_short
        &'is a program for masking of 2D images and volumes. If you want to mask your images with a spherical mask with a soft &
        & falloff, set mskdiam to the diameter in A',&                   ! descr_long
        &'simple_exec',&                                                 ! executable
        &0, 3, 2, 1, 1, 8, 1, .false.)                                   ! # entries in each group, requires sp_project
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
        call mask%set_input('filt_ctrls', 1, lp_backgr)
        ! mask controls
        call mask%set_input('mask_ctrls', 1, mskdiam)
        mask%mask_ctrls(1)%required = .false.
        call mask%set_input('mask_ctrls', 2, mskfile)
        call mask%set_input('mask_ctrls', 3, 'msktype', 'multi', 'Mask type',&
        &'Type of mask to use(soft|hard){soft}', '(soft|hard){soft}', .false., 'soft')
        call mask%set_input('mask_ctrls', 4, mw)
        call mask%set_input('mask_ctrls', 5, width)
        call mask%set_input('mask_ctrls', 6, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6.)
        call mask%set_input('mask_ctrls', 7, 'taper_edges', 'binary', 'Taper edges',&
        &'Whether to taper the edges of image/volume(yes|no){no}', '(yes|no){no}', .false., 'no')
        call mask%set_input('mask_ctrls', 8, 'pdbfile', 'file', 'PDB for 3D envelope masking',&
        &'PDB file used to determine the mask', 'e.g. molecule.pdb', .false., '')
        ! computer controls
        call mask%set_input('comp_ctrls', 1, nthr)
    end subroutine new_mask

    subroutine new_merge_stream_projects
        ! PROGRAM SPECIFICATION
        call merge_stream_projects%new(&
        &'merge_stream_projects',&                   ! name
        &'merge stream two projects',&               ! descr_short
        &'is a program for discarding deselected data (particles,stacks) from a project',& ! descr_long
        &'simple_exec',&                             ! executable
        &0, 2, 0, 0, 0, 0, 0, .false.)               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call merge_stream_projects%set_input('parm_ios', 1, projfile)
        call merge_stream_projects%set_input('parm_ios', 2, projfile_target)
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
    end subroutine new_merge_stream_projects

    subroutine new_mkdir_
        ! PROGRAM SPECIFICATION
        call mkdir_%new(&
        &'mkdir',&                                                       ! name
        &'Make directory',&                                              ! descr_short
        &'is a program for making an automatically numbered directory',& ! descr_long
        &'simple_exec',&                                                 ! executable
        &0, 1, 0, 0, 0, 0, 0, .false.)                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mkdir_%set_input('parm_ios', 1, 'dir', 'dir', 'Name of directory to create', 'Name of directory name to create(e.g. X_dir/)', 'e.g. X_dir/', .true., ' ')
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
    end subroutine new_mkdir_

    subroutine new_motion_correct
        ! PROGRAM SPECIFICATION
        call motion_correct%new(&
        &'motion_correct', &                                                      ! name
        &'Anisotropic motion correction of movies',&                              ! descr_short
        &'is a distributed workflow for anisotropic motion correction of movies.&
        & If then total dose is given the individual frames will be filtered accordingly&
        & (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to&
        & the down-scaling factor (for super-resolution movies). If nframesgrp is given the frames will&
        & be pre-averaged in the given chunk size (Falcon 3 movies).',&     ! descr_long
        &'simple_exec',&                                                    ! executable
        &1, 8, 0, 9, 3, 0, 2, .true.)                                       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call motion_correct%set_input('img_ios', 1, gainref)
        ! parameter input/output
        call motion_correct%set_input('parm_ios', 1, total_dose)
        call motion_correct%set_input('parm_ios', 2, fraction_dose_target)
        call motion_correct%set_input('parm_ios', 3, max_dose)
        call motion_correct%set_input('parm_ios', 4, scale_movies)
        call motion_correct%set_input('parm_ios', 5, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., '')
        call motion_correct%set_input('parm_ios', 6, pspecsz)
        call motion_correct%set_input('parm_ios', 7, eer_fraction)
        call motion_correct%set_input('parm_ios', 8, eer_upsampling)
        ! alternative inputs
        ! <empty>
        ! search controls
        call motion_correct%set_input('srch_ctrls', 1, trs)
        motion_correct%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{20}'
        motion_correct%srch_ctrls(1)%rval_default      = 20.
        call motion_correct%set_input('srch_ctrls', 2, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3){0}', '{0}', .false., 0.)
        call motion_correct%set_input('srch_ctrls', 3, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50.)
        call motion_correct%set_input('srch_ctrls', 4, mcpatch)
        call motion_correct%set_input('srch_ctrls', 5, nxpatch)
        call motion_correct%set_input('srch_ctrls', 6, nypatch)
        call motion_correct%set_input('srch_ctrls', 7, mcconvention)
        call motion_correct%set_input('srch_ctrls', 8, algorithm)
        call motion_correct%set_input('srch_ctrls', 9, mcpatch_thres)
        ! filter controls
        call motion_correct%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call motion_correct%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call motion_correct%set_input('filt_ctrls', 3, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call motion_correct%set_input('comp_ctrls', 1, nparts)
        call motion_correct%set_input('comp_ctrls', 2, nthr)
    end subroutine new_motion_correct

    subroutine new_nununiform_filter2D
        ! PROGRAM SPECIFICATION
        call nununiform_filter2D%new(&
        &'nununiform_filter2D',&                                                  ! name
        &'Optimization (search) based 2D filter (uniform/nonuniform)',&     ! descr_short
        &'is a program for 2D uniform/nonuniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                                                    ! executable
        &3, 1, 0, 0, 4, 0, 1, .false.)                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call nununiform_filter2D%set_input('img_ios', 1, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '')
        call nununiform_filter2D%set_input('img_ios', 2, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '')
        call nununiform_filter2D%set_input('img_ios', 3, 'stk3', 'file', 'Mask stack', 'Mask stack', 'stack_mask.mrc file',  .false., '')
        ! parameter input/output
        call nununiform_filter2D%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call nununiform_filter2D%set_input('filt_ctrls', 1, smooth_ext)
        call nununiform_filter2D%set_input('filt_ctrls', 2, lp_lowres)
        call nununiform_filter2D%set_input('filt_ctrls', 3, nsearch)
        frcs%required = .true.
        call nununiform_filter2D%set_input('filt_ctrls', 4, frcs)
        ! mask controls
        ! <empty>
        ! computer controls
        call nununiform_filter2D%set_input('comp_ctrls', 1, nthr)
    end subroutine new_nununiform_filter2D

    subroutine new_nununiform_filter3D
        ! PROGRAM SPECIFICATION
        call nununiform_filter3D%new(&
        &'nununiform_filter3D',&                                      ! name
        &'Butterworth 3D filter (uniform/nonuniform)',&         ! descr_short
        &'is a program for 3D uniform/nonuniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                                        ! executable
        &2, 1, 0, 0, 5, 2, 1, .false.)                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call nununiform_filter3D%set_input('img_ios', 1, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call nununiform_filter3D%set_input('img_ios', 2, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call nununiform_filter3D%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call nununiform_filter3D%set_input('filt_ctrls', 1, smooth_ext)
        call nununiform_filter3D%set_input('filt_ctrls', 2, lp_lowres)
        call nununiform_filter3D%set_input('filt_ctrls', 3, nsearch)
        call nununiform_filter3D%set_input('filt_ctrls', 4, 'fsc', 'file', 'FSC file', 'FSC file', 'e.g. fsc_state01.bin file', .true., '')
        call nununiform_filter3D%set_input('filt_ctrls', 5, 'lpstop', 'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .false., -1.)
        ! mask controls
        call nununiform_filter3D%set_input('mask_ctrls', 1, mskdiam)
        call nununiform_filter3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call nununiform_filter3D%set_input('comp_ctrls', 1, nthr)
    end subroutine new_nununiform_filter3D

    subroutine new_uniform_filter3D
        ! PROGRAM SPECIFICATION
        call uniform_filter3D%new(&
        &'uniform_filter3D',&                                   ! name
        &'Uniform Butterworth 3D filter',&                      ! descr_short
        &'is a program for 3D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                                        ! executable
        &2, 1, 0, 0, 4, 2, 1, .false.)                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter3D%set_input('img_ios', 1, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call uniform_filter3D%set_input('img_ios', 2, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call uniform_filter3D%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter3D%set_input('filt_ctrls', 1, lp_lowres)
        call uniform_filter3D%set_input('filt_ctrls', 2, nsearch)
        call uniform_filter3D%set_input('filt_ctrls', 3, 'fsc', 'file', 'FSC file', 'FSC file', 'e.g. fsc_state01.bin file', .true., '')
        call uniform_filter3D%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .false., -1.)
        ! mask controls
        call uniform_filter3D%set_input('mask_ctrls', 1, mskdiam)
        call uniform_filter3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call uniform_filter3D%set_input('comp_ctrls', 1, nthr)
    end subroutine new_uniform_filter3D

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
        &'all',&                             ! executable
        &0, 1, 2, 0, 0, 0, 9, .false.)       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call new_project%set_input('parm_ios', 1, user_email)
        ! alternative inputs
        call new_project%set_input('alt_ios', 1, projname)
        new_project%alt_ios(1)%required = .false.
        call new_project%set_input('alt_ios', 2, 'dir', 'dir', 'Project directory', 'Project directory', 'give dir', .false., '')
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
        call new_project%set_input('comp_ctrls', 8, qsys_name)
        call new_project%set_input('comp_ctrls', 9, walltime)
    end subroutine new_new_project

    subroutine new_pick
        ! PROGRAM SPECIFICATION
        call pick%new(&
        &'pick', &                                                         ! name
        &'Template-based particle picking',&                               ! descr_short
        &'is a distributed workflow for template-based particle picking',& ! descr_long
        &'simple_exec',&                                                   ! executable
        &1, 4, 0, 3, 1, 0, 2, .true.)                                      ! # entries in each group, requires sp_project
        pick%gui_submenu_list = "picking"
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call pick%set_input('img_ios', 1, 'pickrefs', 'file', 'Stack of class-averages/reprojections for picking', 'Stack of class-averages/reprojections for picking', 'e.g. pickrefs.mrc', .false., '')
        call pick%set_gui_params('img_ios', 1, submenu="picking", advanced=.false.)
        ! parameter input/output
        call pick%set_input('parm_ios', 1, 'dir', 'dir', 'Output directory', 'Output directory', 'e.g. pick/', .false., 'pick')
        call pick%set_gui_params('parm_ios', 1, submenu="picking")
        call pick%set_input('parm_ios', 2, pcontrast)
        call pick%set_gui_params('parm_ios', 2, submenu="picking")
        call pick%set_input('parm_ios', 3, moldiam)
        call pick%set_gui_params('parm_ios', 3, submenu="picking")
        call pick%set_input('parm_ios', 4, 'picker', 'picker', 'Which picker to use', 'Which picker to use(old|new){old}', '(old|new){old}', .false., 'old')
        call pick%set_gui_params('parm_ios', 4, submenu="picking")
        ! alternative inputs
        ! <empty>
        ! search controls
        call pick%set_input('srch_ctrls', 1, 'thres', 'num', 'Distance threshold in Angs','Distance filter in Angs{24}', '{24}', .false., 24.)
        call pick%set_gui_params('srch_ctrls', 1, submenu="picking", advanced=.false.)
        call pick%set_input('srch_ctrls', 2, 'ndev', 'num', '# of sigmas for outlier detection', '# of standard deviations threshold for outlier detection{2.5}', '{2.5}', .false., 2.5)
        call pick%set_gui_params('srch_ctrls', 2, submenu="picking", advanced=.false.)
        call pick%set_input('srch_ctrls', 3, pgrp)
        call pick%set_gui_params('srch_ctrls', 3, submenu="picking", advanced=.false.)
        pick%srch_ctrls(3)%required = .false.
        ! filter controls
        call pick%set_input('filt_ctrls', 1, 'lp', 'num', 'Low-pass limit','Low-pass limit in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        call pick%set_gui_params('filt_ctrls', 1, submenu="picking")
        ! mask controls
        ! <empty>
        ! computer controls
        call pick%set_input('comp_ctrls', 1, nparts)
        call pick%set_input('comp_ctrls', 2, nthr)
    end subroutine new_pick

    subroutine new_postprocess
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',&                                                      ! name
        &'Post-processing of volume',&                                        ! descr_short
        &'is a program for map post-processing. Use program volops to estimate the B-factor with the Guinier plot',& ! descr_long
        &'simple_exec',&                                                      ! executable
        &0, 1, 0, 0, 6, 7, 1, .true.)                                         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call postprocess%set_input('parm_ios', 1, 'state', 'num', 'State to postprocess', 'State to postprocess{1}', 'Input state{1}', .false., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%set_input('filt_ctrls', 1, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 12.)
        call postprocess%set_input('filt_ctrls', 2, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false., 20.)
        call postprocess%set_input('filt_ctrls', 3, bfac)
        call postprocess%set_input('filt_ctrls', 4, mirr)
        call postprocess%set_input('filt_ctrls', 5, lp_backgr)
        call postprocess%set_input('filt_ctrls', 6, nonuniform)
        ! mask controls
        call postprocess%set_input('mask_ctrls', 1, mskdiam)
        call postprocess%set_input('mask_ctrls', 2, mskfile)
        call postprocess%set_input('mask_ctrls', 3, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false., 1.)
        call postprocess%set_input('mask_ctrls', 4, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for enevloppe mask generation', 'Volume threshold', .false., 0.)
        call postprocess%set_input('mask_ctrls', 5, automsk)
        call postprocess%set_input('mask_ctrls', 6, mw)
        call postprocess%set_input('mask_ctrls', 7, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6.)
        ! computer controls
        call postprocess%set_input('comp_ctrls', 1, nthr)
    end subroutine new_postprocess

    subroutine new_preprocess
        ! PROGRAM SPECIFICATION
        call preprocess%new(&
        &'preprocess', &                                                                    ! name
        &'Preprocessing',&                                                                  ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in sequence',&
        &'simple_exec',&                                                                    ! executable
        &2, 11, 0, 15, 5, 0, 2, .true.)                                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess%set_input('img_ios', 1, gainref)
        call preprocess%set_input('img_ios', 2, 'pickrefs', 'file', 'Reference images for picking', 'Stack of images for picking', 'e.g. cavgs.mrc', .false., '')
        ! parameter input/output
        call preprocess%set_input('parm_ios', 1,  total_dose)
        call preprocess%set_input('parm_ios', 2,  fraction_dose_target)
        call preprocess%set_input('parm_ios', 3,  max_dose)
        call preprocess%set_input('parm_ios', 4,  scale_movies)
        call preprocess%set_input('parm_ios', 5,  eer_fraction)
        call preprocess%set_input('parm_ios', 6,  eer_upsampling)
        call preprocess%set_input('parm_ios', 7,  pcontrast)
        call preprocess%set_input('parm_ios', 8,  'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_')
        call preprocess%set_input('parm_ios', 9,  pspecsz)
        call preprocess%set_input('parm_ios',10,  numlen)
        call preprocess%set_input('parm_ios',11,  ctfpatch)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess%set_input('srch_ctrls', 1, trs)
        preprocess%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{20}'
        preprocess%srch_ctrls(1)%rval_default      = 20.
        call preprocess%set_input('srch_ctrls', 2, 'nframesgrp', 'num', 'Number of contigous frames to sum', '# contigous frames to sum before motion_correct(Falcon 3){0}', '{0}', .false., 0.)
        call preprocess%set_input('srch_ctrls', 3, dfmin)
        call preprocess%set_input('srch_ctrls', 4, dfmax)
        call preprocess%set_input('srch_ctrls', 5, astigtol)
        call preprocess%set_input('srch_ctrls', 6, 'thres', 'num', 'Picking distance threshold','Picking distance filter (in Angs)', 'in Angs{24.}', .false., 24.)
        call preprocess%set_input('srch_ctrls', 7, 'ndev', 'num', '# of sigmas for picking outlier detection', '# of standard deviations threshold for picking oulier detection{2.5}', '{2.5}', .false., 2.5)
        call preprocess%set_input('srch_ctrls', 8, pgrp)
        preprocess%srch_ctrls(8)%required = .false.
        call preprocess%set_input('srch_ctrls', 9, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50.)
        call preprocess%set_input('srch_ctrls',10, mcpatch)
        call preprocess%set_input('srch_ctrls',11, nxpatch)
        call preprocess%set_input('srch_ctrls',12, nypatch)
        call preprocess%set_input('srch_ctrls',13, mcconvention)
        call preprocess%set_input('srch_ctrls',14, algorithm)
        call preprocess%set_input('srch_ctrls',15, mcpatch_thres)
        ! filter controls
        call preprocess%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
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
        &'preprocess_stream', &                                                             ! name
        &'Preprocessing in streaming mode',&                                                ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_exec',&                                                                    ! executable
        &4, 14, 0, 14, 5, 0, 2, .true.)                                                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess_stream%set_input('img_ios', 1, dir_movies)
        call preprocess_stream%set_input('img_ios', 2, gainref)
        call preprocess_stream%set_input('img_ios', 3, 'pickrefs', 'file', 'References images for picking', 'Stack of class-averages for picking', 'e.g. cavgs.mrc', .false., '')
        call preprocess_stream%set_input('img_ios', 4, 'dir_prev', 'file', 'Previous run directory',&
            &'Directory where a previous preprocess_stream application was run', 'e.g. 2_preprocess_stream', .false., '')
        ! parameter input/output
        call preprocess_stream%set_input('parm_ios', 1, total_dose)
        call preprocess_stream%set_input('parm_ios', 2, fraction_dose_target)
        call preprocess_stream%set_input('parm_ios', 3, max_dose)
        call preprocess_stream%set_input('parm_ios', 4, scale_movies)
        call preprocess_stream%set_input('parm_ios', 5, eer_fraction)
        call preprocess_stream%set_input('parm_ios', 6, eer_upsampling)
        call preprocess_stream%set_input('parm_ios', 7, pcontrast)
        call preprocess_stream%set_input('parm_ios', 8, box_extract)
        call preprocess_stream%set_input('parm_ios', 9, pspecsz)
        call preprocess_stream%set_input('parm_ios',10, kv)
        preprocess_stream%parm_ios(10)%required = .true.
        call preprocess_stream%set_input('parm_ios',11, cs)
        preprocess_stream%parm_ios(11)%required = .true.
        call preprocess_stream%set_input('parm_ios',12, fraca)
        preprocess_stream%parm_ios(12)%required = .true.
        call preprocess_stream%set_input('parm_ios',13, smpd)
        preprocess_stream%parm_ios(13)%required = .true.
        call preprocess_stream%set_input('parm_ios',14, ctfpatch)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess_stream%set_input('srch_ctrls', 1, trs)
        preprocess_stream%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{20}'
        preprocess_stream%srch_ctrls(1)%rval_default      = 20.
        call preprocess_stream%set_input('srch_ctrls', 2, dfmin)
        call preprocess_stream%set_input('srch_ctrls', 3, dfmax)
        call preprocess_stream%set_input('srch_ctrls', 4, astigtol)
        call preprocess_stream%set_input('srch_ctrls', 5, 'thres', 'num', 'Picking distance threshold','Picking distance filter (in Angs)', 'in Angs{24.}', .false., 24.)
        call preprocess_stream%set_input('srch_ctrls', 6, 'ndev', 'num', '# of sigmas for picking clustering', '# of standard deviations threshold for picking one cluster clustering{2}', '{2}', .false., 2.)
        call preprocess_stream%set_input('srch_ctrls', 7, pgrp)
        preprocess_stream%srch_ctrls(7)%required = .false.
        call preprocess_stream%set_input('srch_ctrls', 8, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50.)
        call preprocess_stream%set_input('srch_ctrls', 9, mcpatch)
        call preprocess_stream%set_input('srch_ctrls',10, nxpatch)
        call preprocess_stream%set_input('srch_ctrls',11, nypatch)
        call preprocess_stream%set_input('srch_ctrls',12, mcconvention)
        call preprocess_stream%set_input('srch_ctrls',13, algorithm)
        call preprocess_stream%set_input('srch_ctrls',14, mcpatch_thres)
        ! filter controls
        call preprocess_stream%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess_stream%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
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

    subroutine new_preprocess_stream_dev
        ! PROGRAM SPECIFICATION
        call preprocess_stream_dev%new(&
        &'preprocess_stream_dev', &                                                         ! name
        &'Preprocessing in streaming mode',&                                                ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_exec',&                                                                    ! executable
        &5, 15, 0, 22, 10, 1, 9, .true.)                                                    ! # entries in each group, requires sp_project
        preprocess_stream_dev%gui_submenu_list = "data,motion correction,CTF estimation,picking,cluster 2D"
        preprocess_stream_dev%advanced = .false.
        ! image input/output
        call preprocess_stream_dev%set_input('img_ios', 1, dir_movies)
        call preprocess_stream_dev%set_gui_params('img_ios', 1, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('img_ios', 2, gainref)
        call preprocess_stream_dev%set_gui_params('img_ios', 2, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('img_ios', 3, 'pickrefs', 'file', 'References images for picking', 'Stack of class-averages for picking', 'e.g. cavgs.mrc', .false., '')
        call preprocess_stream_dev%set_gui_params('img_ios', 3, submenu="picking", advanced=.false., exclusive_group="pickrefs" )
        call preprocess_stream_dev%set_input('img_ios', 4, 'dir_prev', 'file', 'Previous run directory',&
            &'Directory where a previous preprocess_stream application was run', 'e.g. 2_preprocess_stream', .false., '')
        call preprocess_stream_dev%set_gui_params('img_ios', 4, submenu="data")
        call preprocess_stream_dev%set_input('img_ios', 5, 'dir_meta', 'dir', 'Directory containing per-movie metadata in XML format',&
            &'Directory containing per-movie metadata XML files from EPU', 'e.g. /dataset/metadata', .false., '')
        call preprocess_stream_dev%set_gui_params('img_ios', 5, submenu="data", advanced=.false.)
        ! parameter input/output
        call preprocess_stream_dev%set_input('parm_ios', 1, total_dose)
        call preprocess_stream_dev%set_gui_params('parm_ios', 1, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios', 2, fraction_dose_target)
        call preprocess_stream_dev%set_gui_params('parm_ios', 2, submenu="data")
        call preprocess_stream_dev%set_input('parm_ios', 3, scale_movies)
        call preprocess_stream_dev%set_gui_params('parm_ios', 3, submenu="motion correction", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios', 4, eer_fraction)
        call preprocess_stream_dev%set_gui_params('parm_ios', 4, submenu="motion correction", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios', 5, eer_upsampling)
        call preprocess_stream_dev%set_gui_params('parm_ios', 5, submenu="motion correction", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios', 6, pcontrast)
        call preprocess_stream_dev%set_gui_params('parm_ios', 6, submenu="picking")
        call preprocess_stream_dev%set_input('parm_ios', 7, box_extract)
        call preprocess_stream_dev%set_gui_params('parm_ios', 7, submenu="picking")
        call preprocess_stream_dev%set_input('parm_ios', 8, pspecsz)
        call preprocess_stream_dev%set_gui_params('parm_ios', 8, submenu="picking")
        call preprocess_stream_dev%set_input('parm_ios', 9, max_dose)
        call preprocess_stream_dev%set_gui_params('parm_ios', 9, submenu="motion correction")
        call preprocess_stream_dev%set_input('parm_ios',10, kv)
        preprocess_stream_dev%parm_ios(10)%required = .true.
        call preprocess_stream_dev%set_gui_params('parm_ios', 10, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios',11, cs)
        preprocess_stream_dev%parm_ios(11)%required = .true.
        call preprocess_stream_dev%set_gui_params('parm_ios', 11, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios',12, fraca)
        preprocess_stream_dev%parm_ios(12)%required = .true.
        call preprocess_stream_dev%set_gui_params('parm_ios', 12, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios',13, smpd)
        preprocess_stream_dev%parm_ios(13)%required = .true.
        call preprocess_stream_dev%set_gui_params('parm_ios', 13, submenu="data", advanced=.false.)
        call preprocess_stream_dev%set_input('parm_ios',14, ctfpatch)
        call preprocess_stream_dev%set_gui_params('parm_ios', 14, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('parm_ios',15, prune)
        call preprocess_stream_dev%set_gui_params('parm_ios', 15, submenu="cluster 2D")
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess_stream_dev%set_input('srch_ctrls', 1, trs)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 1, submenu="motion correction")
        preprocess_stream_dev%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{20}'
        preprocess_stream_dev%srch_ctrls(1)%rval_default      = 20.
        call preprocess_stream_dev%set_input('srch_ctrls', 2, dfmin)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 2, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('srch_ctrls', 3, dfmax)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 3, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('srch_ctrls', 4, astigtol)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 4, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('srch_ctrls', 5, 'thres', 'num', 'Picking distance threshold','Picking distance filter (in Angs)', 'in Angs{24.}', .false., 24.)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 5, submenu="picking", advanced=.false.)
        call preprocess_stream_dev%set_input('srch_ctrls', 6, 'ndev',  'num', '# of sigmas for picking clustering', '# of standard deviations threshold for picking one cluster clustering{2}', '{2}', .false., 2.)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 6, submenu="picking", advanced=.false.)
        call preprocess_stream_dev%set_input('srch_ctrls', 7, pgrp)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 7, submenu="picking", advanced=.false.)
        preprocess_stream_dev%srch_ctrls(7)%required = .false.
        call preprocess_stream_dev%set_input('srch_ctrls', 8, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50.)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 8, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls', 9, mcpatch)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 9, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls',10, nxpatch)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 10, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls',11, nypatch)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 11, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls',12, mcconvention)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 12, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls',13, algorithm)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 13, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls',14, ncls_start)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 14, submenu="cluster 2D", advanced=.false.)
        preprocess_stream_dev%srch_ctrls(14)%required = .false.
        call preprocess_stream_dev%set_input('srch_ctrls',15, nptcls_per_cls)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 15, submenu="cluster 2D", advanced=.false.)
        preprocess_stream_dev%srch_ctrls(15)%rval_default = 300.
        call preprocess_stream_dev%set_input('srch_ctrls',16, 'autoscale', 'binary', 'Automatic down-scaling for 2D classification', 'Automatic down-scaling of images &
        &for accelerated convergence rate. (yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 16, submenu="cluster 2D")
        call preprocess_stream_dev%set_input('srch_ctrls',17, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 17, submenu="cluster 2D")
        call preprocess_stream_dev%set_input('srch_ctrls',18, 'ncls', 'num', 'Maximum number of 2D clusters*',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .false., 200.)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 18, submenu="cluster 2D", advanced=.false.)
        call preprocess_stream_dev%set_input('srch_ctrls',19, lpthres)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 19, submenu="cluster 2D", online=.true.)
        call preprocess_stream_dev%set_input('srch_ctrls',20, 'refine', 'multi', 'Refinement mode', '2D Refinement mode(no|greedy){no}', '(no|greedy){no}', .false., 'no')
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 20, submenu="cluster 2D")
        call preprocess_stream_dev%set_input('srch_ctrls',21, mcpatch_thres)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 21, submenu="motion correction")
        call preprocess_stream_dev%set_input('srch_ctrls', 22, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05)
        call preprocess_stream_dev%set_gui_params('srch_ctrls', 22, submenu="motion correction", online=.true.)
        ! filter controls
        call preprocess_stream_dev%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 1, submenu="motion correction")
        call preprocess_stream_dev%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 2, submenu="motion correction")
        call preprocess_stream_dev%set_input('filt_ctrls', 3, 'lp_ctf_estimate', 'num', 'Low-pass limit for CTF parameter estimation',&
        & 'Low-pass limit for CTF parameter estimation in Angstroms{5}', 'in Angstroms{5}', .false., 5.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 3, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('filt_ctrls', 4, 'hp_ctf_estimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 4, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('filt_ctrls', 5, 'lp_pick', 'num', 'Low-pass limit for picking',&
        & 'Low-pass limit for picking in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 5, submenu="picking")
        call preprocess_stream_dev%set_input('filt_ctrls', 6, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 6, submenu="picking")
        call preprocess_stream_dev%set_input('filt_ctrls', 7, 'lp2D', 'num', 'Static low-pass limit for 2D classification', 'Static low-pass limit for 2D classification',&
        &'low-pass limit in Angstroms', .false., 15.)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 7, submenu="cluster 2D")
        call preprocess_stream_dev%set_input('filt_ctrls', 8,  wiener)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 8, submenu="cluster 2D")
        call preprocess_stream_dev%set_input('filt_ctrls', 9, ctfresthreshold)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 9, submenu="CTF estimation")
        call preprocess_stream_dev%set_input('filt_ctrls', 10, icefracthreshold)
        call preprocess_stream_dev%set_gui_params('filt_ctrls', 10, submenu="CTF estimation")
        ! mask controls
        call preprocess_stream_dev%set_input('mask_ctrls', 1, mskdiam)
        call preprocess_stream_dev%set_gui_params('mask_ctrls', 1, submenu="cluster 2D", advanced=.false.)
        preprocess_stream_dev%mask_ctrls(1)%required = .false.
        preprocess_stream_dev%mask_ctrls(1)%descr_short = 'Mask Diameter*'
        ! computer controls
        call preprocess_stream_dev%set_input('comp_ctrls', 1, nchunks)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 1, advanced=.false.)
        preprocess_stream_dev%comp_ctrls(1)%required = .false.
        call preprocess_stream_dev%set_input('comp_ctrls', 2, nparts_chunk)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 2, advanced=.false.)
        call preprocess_stream_dev%set_input('comp_ctrls', 3, nparts_pool)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 3, advanced=.false.)
        call preprocess_stream_dev%set_input('comp_ctrls', 4, nparts)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 4, advanced=.false.)
        preprocess_stream_dev%comp_ctrls(4)%descr_short = 'Number of computing nodes allocated to preprocessing'
        call preprocess_stream_dev%set_input('comp_ctrls', 5, nthr)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 5, advanced=.false.)
        preprocess_stream_dev%comp_ctrls(5)%descr_short = 'Number of threads/node for preprocessing'
        preprocess_stream_dev%comp_ctrls(5)%descr_long  = 'Number of threads per node allocated to preprocessing steps (motion correction, CTF estimation, picking)'
        call preprocess_stream_dev%set_input('comp_ctrls', 6, 'nthr2D', 'num', 'Number of threads/node for 2D classification', 'Number of threads per node allocated to 2D classification',&
        &'# of threads for per', .true., 1.)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 6)
        call preprocess_stream_dev%set_input('comp_ctrls', 7, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740.)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 7)
        call preprocess_stream_dev%set_input('comp_ctrls', 8, 'job_memory_per_task2D','str', 'Memory per 2D computing node',&
        &'Memory dedicated to 2D classification per computing node (in MB){16000}', 'MB per part{16000}', .false., 16000.)
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 8)
        call preprocess_stream_dev%set_input('comp_ctrls', 9, 'qsys_partition2D','str', 'Name of SLURM/PBS partition for 2D classification',&
        &'Name of target partition of distributed computer system (SLURM/PBS) dedicated to 2D classification', 'parttion name', .false., '')
        call preprocess_stream_dev%set_gui_params('comp_ctrls', 9)
    end subroutine new_preprocess_stream_dev

    subroutine new_print_dose_weights
        ! PROGRAM SPECIFICATION
        call print_dose_weights%new(&
        &'print_dose_weights', &                                                  ! name
        &'Print dose weights used in motion correction',&                         ! descr_short
        &'is a program for printing the dose weights used in motion correction',& ! descr_long
        &'simple_exec',&                                                          ! executable
        &0, 5, 0, 0, 0, 0, 0, .false.)                                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_dose_weights%set_input('parm_ios', 1, smpd)
        call print_dose_weights%set_input('parm_ios', 2, box)
        call print_dose_weights%set_input('parm_ios', 3, 'nframes',   'num', 'Number of frames', 'Number of movie frames', '# frames', .true., 0.)
        call print_dose_weights%set_input('parm_ios', 4, kv)
        call print_dose_weights%set_input('parm_ios', 5, total_dose)
        print_dose_weights%parm_ios(5)%required = .true.
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
    end subroutine new_print_dose_weights

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
        call print_magic_boxes%set_input('parm_ios', 3, moldiam)
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

    subroutine new_print_project_field
        ! PROGRAM SPECIFICATION
        call print_project_field%new(&
        &'print_project_field', &                                             ! name
        &'Print project field',&                                              ! descr_short
        &'is a program for printing an orientation field in the project data structure (segment in *.simple project file)',&  ! descr_long
        &'all',&                                                          ! executable
        &0, 1, 0, 0, 0, 0, 0, .true.)                                        ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_project_field%set_input('parm_ios', 1, oritype)
        print_project_field%parm_ios(1)%required = .true.
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
    end subroutine new_print_project_field

    subroutine new_print_project_info
        ! PROGRAM SPECIFICATION
        call print_project_info%new(&
        &'print_project_info', &                                           ! name
        &'Print project info',&                                            ! descr_short
        &'is a program prints information about a *.simple project file',& ! descr_long
        &'all',&                                                        ! executable
        &0, 0, 0, 0, 0, 0, 0, .true.)                                      ! # entries in each group, requires sp_project
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

    subroutine new_reextract
        ! PROGRAM SPECIFICATION
        call reextract%new(&
        &'reextract', &                                                         ! name
        &'Re-extract particle images from integrated movies',&                  ! descr_short
        &'is a program for re-extracting particle images from integrated movies based on determined 2D/3D shifts',& ! descr long
        &'simple_exec',&                                                  ! executable
        &0, 4, 0, 0, 0, 0, 1, .true.)                                           ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call reextract%set_input('parm_ios', 1, box)
        reextract%parm_ios(1)%required = .false.
        call reextract%set_input('parm_ios', 2, oritype)
        reextract%parm_ios(2)%descr_placeholder = '(ptcl2D|ptcl3D){ptcl3D}'
        call reextract%set_input('parm_ios', 3, pcontrast)
        call reextract%set_input('parm_ios', 4, 'ctf', 'multi', 'Whether to extract particles with phases flipped', 'Whether to extract particles with phases flipped(flip|no){no}', '(flip|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call reextract%set_input('comp_ctrls', 1, nparts)
    end subroutine new_reextract

    subroutine new_reproject
        ! PROGRAM SPECIFICATION
        call reproject%new(&
        &'reproject',&                         ! name
        &'Re-project volume',&                 ! descr_short
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
        call reproject%set_input('mask_ctrls', 1, mskdiam)
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
        &0, 3, 2, 0, 0, 1, 1, .false.)         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call normalize_%set_input('parm_ios', 1, smpd)
        call normalize_%set_input('parm_ios', 2, 'norm',       'binary', 'Normalize',       'Statistical normalization: avg=zero, var=1(yes|no){no}',     '(yes|no){no}', .false., 'no')
        call normalize_%set_input('parm_ios', 3, 'noise_norm', 'binary', 'Noise normalize', 'Statistical normalization based on background(yes|no){no}',  '(yes|no){no}', .false., 'no')
        ! alternative inputs
        call normalize_%set_input('alt_ios', 1, 'stk',  'file', 'Stack to normalize',  'Stack of images to normalize', 'e.g. imgs.mrc', .false., '')
        call normalize_%set_input('alt_ios', 2, 'vol1', 'file', 'Volume to normalize', 'Volume to normalize',          'e.g. vol.mrc',  .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call normalize_%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call normalize_%set_input('comp_ctrls', 1, nthr)
    end subroutine new_normalize

    subroutine new_orisops
        ! PROGRAM SPECIFICATION
        call orisops%new(&
        &'orisops',&                      ! name
        &'Standard orientation editing',& ! descr_short
        &'is a program for modifying SIMPLE orientation/parameter files. If errify=yes,&
        & uniform random angular errors, and uniform origin shift errors, &
        & and uniform random defocus errors are introduced. If nstates > 1 then random states are assigned.&
        & If mirr=2d, then the Euler angles in oritab are mirrored according to the relation&
        & e1=e1, e2=180+e2, e3=-e3. If mirr=3d, then the Euler angles in oritab are mirrored according to the&
        & relation R=M(M*R), where R is the rotation matrix calculated from the Euler angle triplet and M is a&
        & 3D reflection matrix. If e1, e2, or e3 is&
        & inputted, the orientations are rotated correspondingly. If you input state,&
        & only the orientations assigned to state state are rotated. If mul is defined, the origin shifts are multiplied with mul.&
        & If zero=yes, then the shifts are zeroed',&
        &'simple_exec',&                  ! executable
        &0, 18, 0, 0, 0, 0, 0, .false.)   ! # entries in each group, requires sp_project
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
        call orisops%set_input('parm_ios', 15, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call orisops%set_input('parm_ios', 16, 'mirr', 'multi', 'Mirror orientations', 'Mirror orientations(2d|3d|no){no}', '(2d|3d|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 17, 'symrnd', 'binary', 'Randomize over subgroubs of point-group', 'Expand orientations over entire unit sphere by &
        &permutation according to randomly selected subgroup symmetry(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%set_input('parm_ios', 18, oritype)
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
        &'is a program for analyzing SIMPLE orientation/parameter files. If two orientation&
        & tables (oritab and oritab2) are inputted, statistics of the distances between the orientations&
        & in the two documents are provided',&
        &'simple_exec',&                          ! executable
        &0, 10, 0, 0, 0, 0, 1, .false.)            ! # entries in each group, requires sp_project
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
        call oristats%set_input('parm_ios', 6,  'ctfstats',  'binary', 'CTF statistics',        'Provide statistics about CTF(yes|no){no}',                      '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 7,  'classtats', 'binary', 'Class statistics',      'Provide statistics about 2D clusters(yes|no){no}',              '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 8,  'projstats', 'binary', 'Projection statistics', 'Provide statistics about projection directions(yes|no){no}',    '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 9,  'trsstats',  'binary', 'Shift statistics',      'Provide statistics about rotational origin shifts(yes|no){no}', '(yes|no){no}', .false., 'no')
        call oristats%set_input('parm_ios', 10, 'specstats', 'binary', 'Specscore statistics',  'Provide statistics about spectral score(yes|no){no}', '(yes|no){no}', .false., 'no')
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

    subroutine new_prune_project
        ! PROGRAM SPECIFICATION
        call prune_project%new(&
        &'prune_project',&                            ! name
        &'discards deselected data from a project',&  ! descr_short
        &'is a program for discarding deselected data (particles,stacks) from a project',& ! descr_long
        &'all',&                                      ! executable
        &0, 1, 0, 0, 0, 0, 1, .true.)                 ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call prune_project%set_input('parm_ios', 1, 'state', 'num', 'State index', 'Index of state to extract', 'give state index', .false., 1.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call prune_project%set_input('comp_ctrls', 1, nparts)
    end subroutine new_prune_project

    subroutine new_reconstruct3D
        ! PROGRAM SPECIFICATION
        call reconstruct3D%new(&
        &'reconstruct3D',&                                               ! name
        &'3D reconstruction from oriented particles',&                   ! descr_short
        &'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,&
        & given input orientations and state assignments. The algorithm is based on direct Fourier inversion&
        & with a Kaiser-Bessel (KB) interpolation kernel',&
        &'simple_exec',&                                                 ! executable
        &0, 0, 0, 2, 3, 2, 2, .true.)                                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call reconstruct3D%set_input('srch_ctrls', 1, pgrp)
        call reconstruct3D%set_input('srch_ctrls', 2, frac)
        ! filter controls
        call reconstruct3D%set_input('filt_ctrls', 1, ptclw)
        call reconstruct3D%set_input('filt_ctrls', 2, envfsc)
        call reconstruct3D%set_input('filt_ctrls', 3, wiener)
        ! mask controls
        call reconstruct3D%set_input('mask_ctrls', 1, mskdiam)
        call reconstruct3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call reconstruct3D%set_input('comp_ctrls', 1, nparts)
        reconstruct3D%comp_ctrls(1)%required = .false.
        call reconstruct3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_reconstruct3D

    subroutine new_refine3D
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',&                                                                               ! name
        &'3D refinement',&                                                                          ! descr_short
        &'is a distributed workflow for 3D refinement based on probabilistic projection matching',& ! descr_long
        &'simple_exec',&                                                                            ! executable
        &1, 0, 0, 13, 13, 4, 2, .true.)                                                             ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%set_input('img_ios', 1, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. vol.mrc', .false., 'vol1.mrc')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D%set_input('srch_ctrls', 1, nspace)
        call refine3D%set_input('srch_ctrls', 2, trs)
        call refine3D%set_input('srch_ctrls', 3, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call refine3D%set_input('srch_ctrls', 4, maxits)
        call refine3D%set_input('srch_ctrls', 5, update_frac)
        call refine3D%set_input('srch_ctrls', 6, frac)
        call refine3D%set_input('srch_ctrls', 7, pgrp)
        call refine3D%set_input('srch_ctrls', 8, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to reconstruct',&
        '# states to reconstruct', .false., 1.0)
        call refine3D%set_input('srch_ctrls', 9, objfun)
        call refine3D%set_input('srch_ctrls', 10, 'refine', 'multi', 'Refinement mode', 'Refinement mode(shc|neigh|cont|inpl|proj|cluster|clustersym){shc}', '(snhc|shc|neigh|cont|inpl|proj|cluster|clustersym){shc}', .false., 'shc')
        call refine3D%set_input('srch_ctrls', 11, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}', '(yes|no){no}', .false., 'no')
        call refine3D%set_input('srch_ctrls', 12, 'lp_iters', 'num', '# iterations prior to e/o refinement', '# of iterations after which low-pass limited alignment is switched to e/o(0:never){1}',&
        &'# iterations prior to e/o refinement{1}', .false., 1.)
        call refine3D%set_input('srch_ctrls', 13, sigma_est)
        ! filter controls
        call refine3D%set_input('filt_ctrls', 1, hp)
        call refine3D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        call refine3D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20.)
        call refine3D%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0)
        call refine3D%set_input('filt_ctrls', 5, lplim_crit)
        call refine3D%set_input('filt_ctrls', 6, lp_backgr)
        call refine3D%set_input('filt_ctrls', 7, ptclw)
        call refine3D%set_input('filt_ctrls', 8, envfsc)
        call refine3D%set_input('filt_ctrls', 9, nonuniform)
        call refine3D%set_input('filt_ctrls', 10, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 12.)
        call refine3D%set_input('filt_ctrls', 11, wiener)
        call refine3D%set_input('filt_ctrls', 12, ml_reg)
        call refine3D%set_input('filt_ctrls', 13, combine_eo)
        ! mask controls
        call refine3D%set_input('mask_ctrls', 1, mskdiam)
        call refine3D%set_input('mask_ctrls', 2, mskfile)
        call refine3D%set_input('mask_ctrls', 3, focusmskdiam)
        call refine3D%set_input('mask_ctrls', 4, automsk)
        ! computer controls
        call refine3D%set_input('comp_ctrls', 1, nparts)
        refine3D%comp_ctrls(1)%required = .false.
        call refine3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_refine3D

    subroutine new_refine3D_nano
        ! PROGRAM SPECIFICATION
        call refine3D_nano%new(&
        &'refine3D_nano',&                                                                                                    ! name
        &'3D refinement of metallic nanoparticles',&                                                                          ! descr_short
        &'is a distributed workflow for 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! descr_long
        &'single_exec',&                                                                                                      ! executable
        &3, 0, 0, 8, 4, 2, 2, .true.)                                                                                         ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D_nano%set_input('img_ios', 1, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .false., '')
        call refine3D_nano%set_input('img_ios', 2, 'vol_odd',  'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .false., '')
        call refine3D_nano%set_input('img_ios', 3, 'vol_even', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .false., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D_nano%set_input('srch_ctrls', 1, nspace)
        call refine3D_nano%set_input('srch_ctrls', 2, trs)
        call refine3D_nano%set_input('srch_ctrls', 3, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call refine3D_nano%set_input('srch_ctrls', 4, maxits)
        call refine3D_nano%set_input('srch_ctrls', 5, update_frac)
        call refine3D_nano%set_input('srch_ctrls', 6, frac)
        call refine3D_nano%set_input('srch_ctrls', 7, pgrp)
        call refine3D_nano%set_input('srch_ctrls', 8, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! filter controls
        call refine3D_nano%set_input('filt_ctrls', 1, hp)
        call refine3D_nano%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5.)
        call refine3D_nano%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms{1.0}', .false., 1.)
        call refine3D_nano%set_input('filt_ctrls', 4, ptclw)
        ! mask controls
        call refine3D_nano%set_input('mask_ctrls', 1, mskdiam)
        call refine3D_nano%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call refine3D_nano%set_input('comp_ctrls', 1, nparts)
        refine3D_nano%comp_ctrls(1)%required = .false.
        call refine3D_nano%set_input('comp_ctrls', 2, nthr)
    end subroutine new_refine3D_nano

    subroutine new_remoc
        ! PROGRAM SPECIFICATION
        call remoc%new(&
        &'remoc', &                             ! name
        &'Regularized motion correction of movies',&     ! descr_short
        &'Regularized motion correction of movies',&     ! descr_long
        &'simple_exec',&                                                    ! executable
        &0, 4, 0, 1, 2, 0, 2, .true.)                                       ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call remoc%set_input('parm_ios', 1, star_ptcl)
        call remoc%set_input('parm_ios', 2, star_mic)
        call remoc%set_input('parm_ios', 3, star_model)
        call remoc%set_input('parm_ios', 4, star_datadir)
        remoc%parm_ios(4)%required = .true.
        ! alternative inputs
        ! <empty>
        ! search controls
        call remoc%set_input('srch_ctrls', 1, trs)
        remoc%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{20}'
        remoc%srch_ctrls(1)%rval_default      = 20.
        ! call motion_correct%set_input('srch_ctrls', 2, mcconvention)
        ! filter controls
        call remoc%set_input('filt_ctrls', 1, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call remoc%set_input('filt_ctrls', 2, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call remoc%set_input('comp_ctrls', 1, nparts)
        call remoc%set_input('comp_ctrls', 2, nthr)
    end subroutine new_remoc

    subroutine new_replace_project_field
        ! PROGRAM SPECIFICATION
        call replace_project_field%new(&
        &'replace_project_field',&                    ! name
        &'hard substitution of project field',&       ! descr_short
        &'is a program for hard substitution of project field, for development purposes',& ! descr_long
        &'simple_exec',&                             ! executable
        &0, 3, 0, 0, 0, 0, 0, .false.)              ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call replace_project_field%set_input('parm_ios', 1, projfile)
        call replace_project_field%set_input('parm_ios', 2, projfile_target)
        call replace_project_field%set_input('parm_ios', 3, oritype)
        replace_project_field%parm_ios(3)%required = .true.
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
    end subroutine new_replace_project_field

    subroutine new_selection
        ! PROGRAM SPECIFICATION
        call selection%new(&
        &'selection',&                                                           ! name
        &'Reports external selection through state 0/1 tags to project',&               ! descr_short
        &'is a program for reporting external (GUI) selections to the SIMPLE project',& ! descr_long
        &'simple_exec',&                                                                ! executable
        &0, 4, 0, 0, 0, 0, 0, .true.)                                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call selection%set_input('parm_ios', 1, 'infile', 'file', 'File with selection state (0/1) flags', 'Plain text file (.txt) with selection state (0/1) flags',&
        &'give .txt selection file', .true., '')
        call selection%set_input('parm_ios', 2, oritype)
        call selection%set_input('parm_ios', 3, 'state', 'num', 'State number', 'Map selection to oris with this state only', '{1}', .false., 1.0)
        call selection%set_input('parm_ios', 4, prune)
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
    end subroutine new_selection

    subroutine new_scale
        ! PROGRAM SPECIFICATION
        call scale%new(&
        &'scale', &                                                                             ! name
        &'Re-scaling MRC and SPIDER stacks and volumes',&                                         ! descr_short
        &'is a program for re-scaling, clipping and padding MRC and SPIDER stacks and volumes',&  ! descr_long
        &'simple_exec',&                                                                        ! executable
        &0, 6, 3, 0, 0, 0, 1, .false.)                                                          ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call scale%set_input('parm_ios', 1, smpd)
        call scale%set_input('parm_ios', 2, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'new box in pixels', .false., 0.)
        call scale%set_input('parm_ios', 3, 'scale', 'num', 'Scaling ratio', 'Target box ratio for scaling(0-1+)', '(0-1+)', .false., 1.)
        call scale%set_input('parm_ios', 4, clip)
        call scale%set_input('parm_ios', 5, outvol)
        call scale%set_input('parm_ios', 6, outstk)
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
        &'scale_project', &                                                                ! name
        &'Re-scaling of MRC and SPIDER stacks',&                                           ! descr_short
        &'is a distributed workflow for re-scaling MRC and SPIDER stacks part of project specification',& ! descr_long
        &'simple_exec',&                                                             ! executable
        &0, 1, 0, 0, 0, 0, 2, .true.)                                                      ! # entries in each group, requires sp_project
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
        &8, 0, 0, 0, 0, 0, 1, .false.)                      ! # entries in each group, requires sp_project
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

    subroutine new_simulate_atoms
        ! PROGRAM SPECIFICATION
        call simulate_atoms%new(&
        &'simulate_atoms',&                                 ! name
        &'Simulate atoms or FCC lattice density',&          ! descr_short
        &'is a program for simulation of atoms or FCC lattice density',& ! descr_long
        &'single_exec',&                                    ! executable
        &2, 4, 0, 0, 0, 0, 1, .false.)                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_atoms%set_input('img_ios', 1, 'pdbfile', 'file', 'PDB', 'Input coordinates file in PDB format', 'Input coordinates file', .false., '')
        call simulate_atoms%set_input('img_ios', 2, outvol)
        ! parameter input/output
        call simulate_atoms%set_input('parm_ios',  1, smpd)
        call simulate_atoms%set_input('parm_ios',  2, box)
        call simulate_atoms%set_input('parm_ios',  3, element)
        simulate_atoms%parm_ios(3)%required = .false.
        call simulate_atoms%set_input('parm_ios',  4, moldiam)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_atoms%set_input('comp_ctrls', 1, nthr)
    end subroutine new_simulate_atoms

    subroutine new_simulate_movie
        ! PROGRAM SPECIFICATION
        call simulate_movie%new(&
        &'simulate_movie',&                                 ! name
        &'Simulate DDD movie',&                             ! descr_short
        &'is a program for crude simulation of a DDD movie. Input is a set of projection images to place. &
        &Movie frames are then generated related by randomly shifting the base image and applying noise',& ! descr_long
        &'simple_exec',&                                    ! executable
        &1, 10, 0, 0, 1, 0, 1, .false.)                     ! # entries in each group, requires sp_project
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
        &0, 2, 0, 0, 0, 0, 0, .false.)                     ! # entries in each group, requires sp_project
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
        &'simulate_particles',&                                           ! name
        &'Simulate single-particle images',&                              ! descr_short
        &'is a program for simulating single-particle cryo-EM images. It is not a very sophisticated simulator, but&
        & it is nevertheless useful for testing purposes. It does not do any multi-slice simulation and it cannot be&
        & used for simulating molecules containing heavy atoms. It does not even accept a PDB file as an input. Input&
        & is a cryo-EM map, which we usually generate from a PDB file using EMANs program pdb2mrc. The volume is&
        & projected using Fourier interpolation, 20% of the total noise is added to the images (pink noise), they are&
        & then Fourier transformed and multiplied with astigmatic CTF and B-factor. Next, the they are inverse FTed&
        & before the remaining 80% of the noise (white noise) is added',& ! descr_long
        &'simple_exec',&                                                  ! executable
        &1, 16, 0, 1, 2, 1, 1, .false.)                                   ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_particles%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to project', 'input volume e.g. vol.mrc', .false., '')
        simulate_particles%img_ios(1)%required = .true.
        ! parameter input/output
        call simulate_particles%set_input('parm_ios', 1,  smpd)
        call simulate_particles%set_input('parm_ios', 2,  nptcls)
        call simulate_particles%set_input('parm_ios', 3,  'snr', 'num', 'SNR', 'Signal-to-noise ratio of particle images', 'signal-to-noise ratio(0.)', .true., 0.)
        call simulate_particles%set_input('parm_ios', 4,  oritab)
        call simulate_particles%set_input('parm_ios', 5,  outfile)
        call simulate_particles%set_input('parm_ios', 6,  outstk)
        call simulate_particles%set_input('parm_ios', 7,  'even', 'binary', 'Generate even projections', 'Generate quasi-even projection directions(yes|no){no}', '(yes|no){no}', .false., 'no')
        call simulate_particles%set_input('parm_ios', 8,  sherr)
        call simulate_particles%set_input('parm_ios', 9,  kv)
        call simulate_particles%set_input('parm_ios', 10, cs)
        call simulate_particles%set_input('parm_ios', 11, fraca)
        call simulate_particles%set_input('parm_ios', 12, deftab)
        call simulate_particles%set_input('parm_ios', 13, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2.)
        call simulate_particles%set_input('parm_ios', 14, dferr)
        call simulate_particles%set_input('parm_ios', 15, 'astigerr', 'num', 'Astigmatism error', 'Uniform astigmatism error(in microns)', 'error in microns', .false., 0.)
        call simulate_particles%set_input('parm_ios', 16, ctf)
        ! alternative inputs
        ! <empty>
        ! search controls
        call simulate_particles%set_input('srch_ctrls', 1, pgrp)
        ! filter controls
        call simulate_particles%set_input('filt_ctrls', 1, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        call simulate_particles%set_input('filt_ctrls', 2, 'bfacerr', 'num', 'B-factor error', 'Uniform B-factor error(in Angstroms^2)', 'error(in Angstroms^2)', .false., 50.)
        ! mask controls
        call simulate_particles%set_input('mask_ctrls', 1, mskdiam)
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
        &1, 3, 0, 0, 0,0, 1, .false.)                           ! # entries in each group, requires sp_project
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

    subroutine new_split_
        call split_%new(&
        &'split',&                                   ! name
        &'Split stack into substacks',&              ! descr_short
        &'is a program for splitting a stack into evenly partitioned substacks',& ! descr_long
        &'simple_exec',&               ! executable
        &1, 1, 0, 0, 0, 0, 1, .false.) ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call split_%set_input('img_ios', 1, stk)
        split_%img_ios(1)%required = .true.
        ! parameter input/output
        call split_%set_input('parm_ios', 1, smpd)
        ! computer controls
        call split_%set_input('comp_ctrls', 1, nparts)
    end subroutine new_split_

    subroutine new_stack
        ! PROGRAM SPECIFICATION
        call stack%new(&
        &'stack',&                     ! name
        &'Stack images',&              ! descr_short
        &'is a program for stacking individual images (list) or multiple stacks into one',& ! descr_long
        &'simple_exec',&               ! executable
        &2, 2, 0, 0, 0, 0, 0, .false.) ! # entries in each group, requires sp_project
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
        &'stackops',&                                ! name
        &'Standard stack editing',&                  ! descr_short
        &'is a program that provides standard single-particle image processing routines for MRC or SPIDER&
        & stacks. To extract a particular state, give oritab and set state.&
        & To select the fraction of best particles, give oritab&
        & and set frac. State and frac options can be combined.&
        & To apply noise, give the desired signal-to-noise ratio via snr. To calculate the autocorrelation&
        & function, set acf=yes. To extract a contiguous subset of images from the stack, set&
        & fromp and top. To extract a number of particle images at random, set nran to the desired number.&
        & With avg=yes the global average of the stack is calculated.&
        & If nframesgrp is set to some integer number >1, averages with chunk sizes of nframesgrp are produced,&
        & which may be useful for analysis of dose-fractionated image series. neg inverts the contrast of the images',& ! descr_long
        &'simple_exec',&                             ! executable
        &2, 20, 0, 0, 0, 0, 1, .false.)              ! # entries in each group, requires sp_project
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
        call stackops%set_input('parm_ios', 18, 'roavg', 'binary', 'Rotationally average', 'Rotationally average images in stack(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%set_input('parm_ios', 19, 'angstep', 'num', 'Angular stepsize', 'Angular stepsize for rotational averaging(in degrees)', 'give degrees', .false., 5.)
        call stackops%set_input('parm_ios', 20, 'makemovie', 'binary', 'Whether to make a movie', 'Generates images and script to make a movie with FFmpeg(yes|no){no}', '(yes|no){no}', .false., 'no')
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

    subroutine new_symaxis_search
        ! PROGRAM SPECIFICATION
        call symaxis_search%new(&
        &'symaxis_search',&                                                                                 ! name
        &'Search for symmetry axis',&                                                                       ! descr_short
        &'is a program for searching for the principal symmetry axis of a volume reconstructed in C1. &
        &The rotational transformation is applied to the oritype field in the project and the project &
        &file is updated. If you are unsure about the point-group, use the symmetry_test program instead',& ! descr_long
        &'simple_exec',&                                                                                    ! executable
        &1, 1, 0, 2, 3, 1, 1, .false.)                                                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symaxis_search%set_input('img_ios', 1, 'vol1', 'file', 'C1 Volume to identify symmetry axis of', 'C1 Volume to identify symmetry axis of', &
        & 'input volume e.g. vol_C1.mrc', .true., '')
        ! parameter input/output
        call symaxis_search%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symaxis_search%set_input('srch_ctrls', 1, pgrp)
        call symaxis_search%set_input('srch_ctrls', 2, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symaxis_search%set_input('filt_ctrls', 1, lp)
        call symaxis_search%set_input('filt_ctrls', 2, hp)
        call symaxis_search%set_input('filt_ctrls', 3, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        call symaxis_search%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call symaxis_search%set_input('comp_ctrls', 1, nthr)
    end subroutine new_symaxis_search

    subroutine new_symmetrize_map
        ! PROGRAM SPECIFICATION
        call symmetrize_map%new(&
        &'symmetrize_map',&                                                                                           ! name
        &'Symmetrization of density map',&                                                                            ! descr_short
        &'is a program that implements symmetrization of the input density map. &
        &Input is a volume and point-group symmetry, output is the volume aligned to the principal symmetry axis and averaged over the symmetry operations',& ! descr long
        &'simple_exec',&                                                                                             ! executable
        &2, 1, 0, 2, 3, 1, 1, .false.)                                                                               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symmetrize_map%set_input('img_ios', 1, 'vol1', 'file', 'Volume to symmetrize', 'Volume to symmetrize', &
        & 'input volume e.g. vol.mrc', .true., '')
        call symmetrize_map%set_input('img_ios', 2, outvol)
        ! parameter input/output
        call symmetrize_map%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symmetrize_map%set_input('srch_ctrls', 1, pgrp)
        call symmetrize_map%set_input('srch_ctrls', 2, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symmetrize_map%set_input('filt_ctrls', 1, lp)
        call symmetrize_map%set_input('filt_ctrls', 2, hp)
        call symmetrize_map%set_input('filt_ctrls', 3, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 20.)
        ! mask controls
        call symmetrize_map%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call symmetrize_map%set_input('comp_ctrls', 1, nthr)
    end subroutine new_symmetrize_map

    subroutine new_symmetry_test
        ! PROGRAM SPECIFICATION
        call symmetry_test%new(&
        &'symmetry_test',&                                                                                           ! name
        &'Statistical test for symmetry',&                                                                           ! descr_short
        &'is a program that implements a statistical test for point-group symmetry. &
        &Input is a volume reconstructed without symmetry (c1) and output is the most likely point-group symmetry',& ! descr long
        &'simple_exec',&                                                                                             ! executable
        &1, 1, 0, 3, 3, 1, 1, .false.)                                                                               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symmetry_test%set_input('img_ios', 1, 'vol1', 'file', 'C1 Volume to identify symmetry of', 'C1 Volume to identify symmetry of', &
        & 'input volume e.g. vol_C1.mrc', .true., '')
        ! parameter input/output
        call symmetry_test%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symmetry_test%set_input('srch_ctrls', 1, 'cn_stop',  'num', 'Rotational symmetry order stop index',  'Rotational symmetry order stop index',  'give stop index',  .false., 10.)
        call symmetry_test%set_input('srch_ctrls', 2, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call symmetry_test%set_input('srch_ctrls', 3, 'platonic', 'binary', 'Search for Platonic symmetries', 'Search for Platonic symmetries(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symmetry_test%set_input('filt_ctrls', 1, lp)
        call symmetry_test%set_input('filt_ctrls', 2, hp)
        call symmetry_test%set_input('filt_ctrls', 3, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        call symmetry_test%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call symmetry_test%set_input('comp_ctrls', 1, nthr)
    end subroutine new_symmetry_test

    subroutine new_atoms_stats
        ! PROGRAM SPECIFICATION
        call atoms_stats%new(&
        &'atoms_stats',&                                                                              ! name
        &'Statistical test for radial dependent symmetry',&                                           ! descr_short
        &'is a program that generates statistics at different radii and across the whold nano map.',& ! descr long
        &'single_exec',&                                                                              ! executable
        &3, 4, 0, 0, 1, 1, 1, .false.)                                                                ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call atoms_stats%set_input('img_ios', 1, 'vol1', 'file', 'Raw volume', 'Raw volume of grey valued pixel intensities', &
        & 'input volume e.g. vol.mrc', .true., '')
        call atoms_stats%set_input('img_ios', 2, 'vol2', 'file', 'Connected components volume', 'Connected components volume produced by detect atoms', &
        & 'input volume e.g. *CC.mrc', .true., '')
        call atoms_stats%set_input('img_ios', 3, 'vol3', 'file', 'Volume', 'Nanoparticle volume to use for lattice fitting', &
        & 'input volume 4 lattice fit e.g. vol3.mrc', .false., '')
        ! parameter input/output
        call atoms_stats%set_input('parm_ios', 1, smpd)
        call atoms_stats%set_input('parm_ios', 2, 'pdbfile', 'file', 'PDB', 'Input coords file in PDB format', 'Input coords file in PDB format', .true., '')
        call atoms_stats%set_input('parm_ios', 3, 'pdbfile2', 'file', 'PDB', 'subset coords for stats calc', 'subset coords file in PDB format for stats calc', .false., '')
        call atoms_stats%set_input('parm_ios', 4, 'rmsd_file','file', 'bin', 'per-atom e/o rmsd:s', 'per-atom e/o rmsd:s from CS model building', .false., '')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call atoms_stats%set_input('filt_ctrls', 1, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        ! mask controls
        call atoms_stats%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call atoms_stats%set_input('comp_ctrls', 1, nthr)
    end subroutine new_atoms_stats

    subroutine new_tseries_atoms_analysis
        ! PROGRAM SPECIFICATION
        call tseries_atoms_analysis%new(&
        &'tseries_atoms_analysis',&                                                   ! name
        &'Analysis of results obtianed with tseries_reconstruct3D and detect_atoms',& ! descr_short
        &'is a program that analysis atomic time-series coordinates',&                ! descr long
        &'single_exec',&                                                              ! executable
        &0, 3, 0, 0, 1, 0, 0, .false.)                                                ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_atoms_analysis%set_input('parm_ios', 1, smpd)
        call tseries_atoms_analysis%set_input('parm_ios', 2, 'pdbfiles',  'file', 'txt', 'List of PDB format coords files',  'List of input coords files in PDB format', .true., '')
        call tseries_atoms_analysis%set_input('parm_ios', 3, 'frac_diam', 'num',  'Fraction of atomic diameter', 'Fraction of atomic diameter used for thresholding{0.8}', '{0.8}', .false., 0.8)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call tseries_atoms_analysis%set_input('filt_ctrls', 1, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_tseries_atoms_analysis

    subroutine new_tseries_import
        ! PROGRAM SPECIFICATION
        call tseries_import%new(&
        &'tseries_import',&                               ! name
        &'Imports time-series datasets',&                 ! descr_short
        &'is a workflow for importing time-series data',& ! descr_long
        &'single_exec',&                                  ! executable
        &1, 4, 0, 0, 0, 0, 0, .true.)                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tseries_import%set_input('img_ios', 1, 'filetab', 'file', 'List of individual movie frame files', 'List of frame files (*.mrcs) to import', 'e.g. movie_frames.txt', .true., '')
        ! parameter input/output
        call tseries_import%set_input('parm_ios', 1, smpd)
        call tseries_import%set_input('parm_ios', 2, kv)
        tseries_import%parm_ios(2)%required = .true.
        call tseries_import%set_input('parm_ios', 3, 'cs', 'num', 'Spherical aberration', 'Spherical aberration constant(in mm){0.0}', 'in mm{0.0}', .true., 0.0)
        call tseries_import%set_input('parm_ios', 4, 'fraca', 'num', 'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.4}', 'fraction{0.4}', .true., 0.4)
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
    end subroutine new_tseries_import

    subroutine new_tseries_import_particles
        ! PROGRAM SPECIFICATION
        call tseries_import_particles%new(&
        &'tseries_import_particles',&                     ! name
        &'Imports time-series particles stack',&          ! descr_short
        &'is a workflow for importing time-series data',& ! descr_long
        &'single_exec',&                                  ! executable
        &1, 2, 0, 0, 0, 0, 0, .true.)                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tseries_import_particles%set_input('img_ios', 1, stk)
        tseries_import_particles%img_ios(1)%required = .true.
        ! parameter input/output
        call tseries_import_particles%set_input('parm_ios', 1, smpd)
        call tseries_import_particles%set_input('parm_ios', 2, deftab)
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
    end subroutine new_tseries_import_particles

    subroutine new_tseries_motion_correct
        ! PROGRAM SPECIFICATION
        call tseries_motion_correct%new(&
        &'tseries_motion_correct', &                                       ! name
        &'Anisotropic motion correction of time-series of nanoparticles',& ! descr_short
        &'is a distributed workflow for anisotropic motion correction of time-series (movies) of nanoparticles.',& ! descr_long
        &'single_exec',&                                                   ! executable
        &0, 1, 0, 7, 3, 0, 2, .true.)                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_motion_correct%set_input('parm_ios', 1, 'boxfile', 'file', 'List of particle coordinates',&
        &'.txt file with EMAN-convention particle coordinates', 'e.g. coords.box', .false., '')        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_motion_correct%set_input('srch_ctrls', 1, trs)
        tseries_motion_correct%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{10}'
        tseries_motion_correct%srch_ctrls(1)%rval_default      = 10.
        call tseries_motion_correct%set_input('srch_ctrls', 2, 'nframesgrp', 'num', '# frames in time moving time window', '# frames in time moving time window subjected to correction', '{5}', .false., 5.)
        call tseries_motion_correct%set_input('srch_ctrls', 3, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{5}', .false., 5.)
        call tseries_motion_correct%set_input('srch_ctrls', 4, mcpatch)
        call tseries_motion_correct%set_input('srch_ctrls', 5, nxpatch)
        tseries_motion_correct%srch_ctrls(5)%descr_placeholder = '# x-patches{3}'
        tseries_motion_correct%srch_ctrls(5)%rval_default = 3.
        call tseries_motion_correct%set_input('srch_ctrls', 6, nypatch)
        tseries_motion_correct%srch_ctrls(6)%descr_placeholder = '# y-patches{3}'
        tseries_motion_correct%srch_ctrls(6)%rval_default = 3.
        call tseries_motion_correct%set_input('srch_ctrls', 7, algorithm)
        ! filter controls
        call tseries_motion_correct%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call tseries_motion_correct%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){3}', 'in Angstroms{3}', .false., 3.)
        call tseries_motion_correct%set_input('filt_ctrls', 3, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_motion_correct%set_input('comp_ctrls', 1, nparts)
        call tseries_motion_correct%set_input('comp_ctrls', 2, nthr)
    end subroutine new_tseries_motion_correct

    subroutine new_tseries_make_pickavg
        ! PROGRAM SPECIFICATION
        call tseries_make_pickavg%new(&
        &'tseries_make_pickavg',&                                                        ! name
        &'Align & average the first few frames of the time-series',&                     ! descr_short
        &'is a program for aligning & averaging the first few frames of the time-series&
        & to accomplish SNR enhancement for particle identification',&                   ! descr_long
        &'single_exec',&                                                                 ! executable
        &0, 1, 0, 5, 3, 0, 1, .true.)                                                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_make_pickavg%set_input('parm_ios', 1, 'nframesgrp', 'num', '# contigous frames to average', 'Number of contigous frames to average using correlation-based weights{10}', '{10}', .false., 10.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_make_pickavg%set_input('srch_ctrls', 1, trs)
        tseries_make_pickavg%srch_ctrls(1)%descr_placeholder = 'max shift per iteration in pixels{10}'
        tseries_make_pickavg%srch_ctrls(1)%rval_default      = 10.
        call tseries_make_pickavg%set_input('srch_ctrls', 2, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{5}', .false., 5.)
        call tseries_make_pickavg%set_input('srch_ctrls', 3, mcpatch)
        call tseries_make_pickavg%set_input('srch_ctrls', 4, nxpatch)
        tseries_make_pickavg%srch_ctrls(4)%rval_default = 3.
        tseries_make_pickavg%srch_ctrls(4)%descr_placeholder = '# x-patches{3}'
        call tseries_make_pickavg%set_input('srch_ctrls', 5, nypatch)
        tseries_make_pickavg%srch_ctrls(5)%rval_default = 3.
        tseries_make_pickavg%srch_ctrls(5)%descr_placeholder = '# x-patches{3}'
        ! filter controls
        call tseries_make_pickavg%set_input('filt_ctrls', 1, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call tseries_make_pickavg%set_input('filt_ctrls', 2, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){3}', 'in Angstroms{3}', .false., 3.)
        call tseries_make_pickavg%set_input('filt_ctrls', 3, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_make_pickavg%set_input('comp_ctrls', 1, nthr)
    end subroutine new_tseries_make_pickavg

    subroutine new_tseries_swap_stack
        ! PROGRAM SPECIFICATION
        call tseries_swap_stack%new(&
        &'tseries_swap_stack',&                                           ! name
        &'Substitutes stack into an existing project',&                   ! descr_short
        &'is a program for substituting stack into an existing project',& ! descr_long
        &'single_exec',&                                                  ! executable
        &1, 0, 0, 0, 0, 0, 0, .true.)                                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tseries_swap_stack%set_input('img_ios', 1, stk)
        tseries_swap_stack%img_ios(1)%required = .true.
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
    end subroutine new_tseries_swap_stack

    subroutine new_tseries_track_particles
        ! PROGRAM SPECIFICATION
        call tseries_track_particles%new(&
        &'tseries_track_particles',&                                             ! name
        &'Track particles in time-series',&                                      ! descr_short
        &'is a distributed workflow for particle tracking in time-series data',& ! descr_long
        &'single_exec',&                                                         ! executable
        &0, 3, 0, 2, 4, 0, 1, .true.)                                            ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_track_particles%set_input('parm_ios', 1, 'fbody', 'string', 'Template output tracked series',&
        &'Template output tracked series', 'e.g. tracked_ptcl', .true., '')
        call tseries_track_particles%set_input('parm_ios', 2, 'boxfile', 'file', 'List of particle coordinates',&
        &'.txt file with EMAN particle coordinates', 'e.g. coords.box', .true., '')
        call tseries_track_particles%set_input('parm_ios', 3, 'neg', 'binary', 'Invert contrast', 'Invert image contrast(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_track_particles%set_input('srch_ctrls', 1, 'offset', 'num', 'Shift half-width search bound', 'Shift half-width search bound(in pixels)',&
        'e.g. pixels window halfwidth', .false., 10.)
        call tseries_track_particles%set_input('srch_ctrls', 2, 'nframesgrp', 'num', 'Number of contigous frames to average', '# contigous frames to average before tracking{30}', '{30}', .false., 30.)
        ! <empty>
        ! filter controls
        call tseries_track_particles%set_input('filt_ctrls', 1, lp)
        tseries_track_particles%filt_ctrls(1)%required     = .false.
        tseries_track_particles%filt_ctrls(1)%rval_default = 2.3
        tseries_track_particles%filt_ctrls(1)%descr_placeholder = 'Low-pass limit in Angstroms{2.3}'
        call tseries_track_particles%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the particle and centering', 'centering low-pass limit in Angstroms{5}', .false., 5.)
        call tseries_track_particles%set_input('filt_ctrls', 3, 'filter', 'multi','Alternative filter for particle tracking',&
            &'Alternative filter for particle tracking(no|tv|nlmean){tv}', '(no|tv|nlmean){tv}', .false., 'tv')
        call tseries_track_particles%set_input('filt_ctrls', 4, hp)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_track_particles%set_input('comp_ctrls', 1, nthr)
    end subroutine new_tseries_track_particles

    subroutine new_tseries_reconstruct3D
        ! PROGRAM SPECIFICATION
        call tseries_reconstruct3D%new(&
        &'tseries_reconstruct3D',&                                       ! name
        &'Time windowed 3D reconstruction from oriented particles',&     ! descr_long
        &'Time windowed 3D reconstruction from oriented particles',&
        &'single_exec',&                                                 ! executable
        &0, 3, 0, 1, 0, 2, 2, .true.)                                    ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_reconstruct3D%set_input('parm_ios', 1, 'stepsz',  'num', 'Time window size (# frames){500}', 'Time window size (# frames) for windowed 3D rec{500}', 'give # frames',  .false., 500.)
        call tseries_reconstruct3D%set_input('parm_ios', 2, 'fromp', 'num', 'From particle index', 'Start index for 3D reconstruction', 'start index', .false., 1.0)
        call tseries_reconstruct3D%set_input('parm_ios', 3, 'top',   'num', 'To particle index', 'Stop index for 3D reconstruction', 'stop index', .false., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_reconstruct3D%set_input('srch_ctrls', 1, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call tseries_reconstruct3D%set_input('mask_ctrls', 1, mskdiam)
        call tseries_reconstruct3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call tseries_reconstruct3D%set_input('comp_ctrls', 1, nparts)
        tseries_reconstruct3D%comp_ctrls(1)%required = .false.
        call tseries_reconstruct3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_tseries_reconstruct3D

    subroutine new_graphene_subtr
        ! PROGRAM SPECIFICATION
        call graphene_subtr%new(&
        &'graphene_subtr',&                        ! name
        &'Removes graphene Fourier peaks in time-series',& ! descr_short
        &'Removes graphene Fourier peaks in time-series',& ! descr_long
        &'single_exec',&                                   ! executable
        &3, 1, 0, 0, 0, 0, 1, .false.)                     ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call graphene_subtr%set_input('img_ios', 1, stk)
        graphene_subtr%img_ios(1)%required = .true.
        graphene_subtr%img_ios(1)%descr_placeholder = 'Input tracked particles, eg NP_X.mrc'
        call graphene_subtr%set_input('img_ios', 2, stk2)
        graphene_subtr%img_ios(2)%required = .true.
        graphene_subtr%img_ios(2)%descr_placeholder = 'Input background power spectra stack, eg NP_X_background_pspec.mrc'
        call graphene_subtr%set_input('img_ios', 3, outstk)
        ! parameter input/output
        call graphene_subtr%set_input('parm_ios', 1, smpd)
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
        call graphene_subtr%set_input('comp_ctrls', 1, nthr)
    end subroutine new_graphene_subtr

    subroutine new_uniform_filter2D
        ! PROGRAM SPECIFICATION
        call uniform_filter2D%new(&
        &'uniform_filter2D',&           ! name
        &'Uniform 2D filter',&           ! descr_short
        &'is a program for 2D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                 ! executable
        &3, 1, 0, 0, 3, 0, 1, .false.)                                      ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter2D%set_input('img_ios', 1, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '')
        call uniform_filter2D%set_input('img_ios', 2, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '')
        call uniform_filter2D%set_input('img_ios', 3, 'stk3', 'file', 'Mask stack', 'Mask stack', 'stack_mask.mrc file',  .false., '')
        ! parameter input/output
        call uniform_filter2D%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter2D%set_input('filt_ctrls', 1, lp_lowres)
        call uniform_filter2D%set_input('filt_ctrls', 2, nsearch)
        frcs%required = .true.
        call uniform_filter2D%set_input('filt_ctrls', 3, frcs)
        ! mask controls
        ! <empty>
        ! computer controls
        call uniform_filter2D%set_input('comp_ctrls', 1, nthr)
    end subroutine new_uniform_filter2D

    subroutine new_update_project
        ! PROGRAM SPECIFICATION
        call update_project%new(&
        &'update_project',&                  ! name
        &'Update an existing project',&      ! descr_short
        &'is a program for updating an existing project: changing the name/user_email/computer controls',& ! descr_long
        &'all',&                          ! executable
        &0, 2, 0, 0, 0, 0, 9, .true.)        ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call update_project%set_input('parm_ios', 1, projfile)
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
        call update_project%set_input('comp_ctrls', 8, qsys_name)
        call update_project%set_input('comp_ctrls', 9, walltime)
    end subroutine new_update_project

    subroutine new_vizoris
        ! PROGRAM SPECIFICATION
        call vizoris%new(&
        &'vizoris',&                                                                                               ! name
        &'Visualization of orientation distribution',&                                                             ! descr_short
        &'is a program for extracting projection directions from orientations for visualization in UCSF Chimera',& ! descr_long
        &'all',&                                                                                                   ! executable
        &0, 5, 0, 0, 0, 0, 0, .false.)                                                                             ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call vizoris%set_input('parm_ios', 1, oritab)
        vizoris%parm_ios(1)%required = .true.
        call vizoris%set_input('parm_ios', 2, nspace)
        call vizoris%set_input('parm_ios', 3, pgrp)
        call vizoris%set_input('parm_ios', 4, oritype)
        call vizoris%set_input('parm_ios', 5, 'tseries', 'binary', 'Time series analysis', 'Orientations originate from analysis of a time-series(yes|no){no}', '(yes|no){no}', .false., 'no')
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
        &2, 13, 0, 0, 3, 1, 1, .false.)                                                           ! # entries in each group, requires sp_project
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
        call volops%set_input('parm_ios',13, mul)
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
        call volops%set_input('mask_ctrls', 1, mskdiam)
        volops%mask_ctrls(1)%required = .false.
        ! computer controls
        call volops%set_input('comp_ctrls', 1, nthr)
    end subroutine new_volops

    subroutine new_whiten_and_filter
        ! PROGRAM SPECIFICATION
        call whiten_and_filter%new(&
        &'whiten_and_filter', &                                                      ! name
        &'Whiten the noise power and apply optimal FSC filter',&                     ! descr_short
        &'is a program for whitening the noise power filter the two input volumes',& ! descr_long
        &'simple_exec',&                                                             ! executable
        &2, 1, 0, 0, 0, 1, 1, .false.)                                               ! # entries in each group, requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call whiten_and_filter%set_input('img_ios', 1, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call whiten_and_filter%set_input('img_ios', 2, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call whiten_and_filter%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call whiten_and_filter%set_input('mask_ctrls', 1, mskdiam)
        ! computer controls
        call whiten_and_filter%set_input('comp_ctrls', 1, nthr)
    end subroutine new_whiten_and_filter

    subroutine new_write_classes
        ! PROGRAM SPECIFICATION
        call write_classes%new(&
        &'write_classes',&                                                                                  ! name
        &'Writes the class averages and the individual (rotated and shifted) particles part of the class',& ! descr_short
        &'is a program for the class averages and the individual (rotated and shifted) particles part of the classto to individual stacks',& ! descr_long
        &'simple_exec',&                                                                                    ! executable
        &0, 0, 0, 0, 0, 0, 0, .true.)                                                                       ! # entries in each group, requires sp_project
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
    end subroutine new_write_classes

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
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; simple_program :: set_input_1')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)    :: i
                type(simple_input_param), intent(inout) :: arr(:)
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
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; simple_program :: set_input_2')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)  :: i
                type(simple_input_param), intent(inout) :: arr(:)
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
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; simple_program :: set_input_3')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)  :: i
                type(simple_input_param), intent(inout) :: arr(:)
                allocate(arr(i)%key,               source=trim(param%key))
                allocate(arr(i)%keytype,           source=trim(param%keytype))
                allocate(arr(i)%descr_short,       source=trim(param%descr_short))
                allocate(arr(i)%descr_long,        source=trim(param%descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(param%descr_placeholder))
                arr(i)%required = param%required
            end subroutine set

    end subroutine set_input_3

    subroutine set_gui_params( self, which, i, submenu, exclusive_group, advanced, active_flags, online)
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        character(len=*),              intent(in),optional    :: submenu
        character(len=*),              intent(in),optional    :: exclusive_group
        character(len=*),              intent(in),optional    :: active_flags
        logical,                       intent(in),optional    :: advanced
        logical,                       intent(in),optional    :: online
        integer,                       intent(in)    :: i
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
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; simple_program :: set_input_3')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)  :: i
                type(simple_input_param), intent(inout) :: arr(:)
                if( present(submenu) ) then
                  allocate(arr(i)%gui_submenu, source=trim(submenu))
                endif
                if( present (exclusive_group) ) then
                  allocate(arr(i)%exclusive_group, source=trim(exclusive_group))
                endif
                if( present(advanced) ) then
                  arr(i)%advanced = advanced
                endif
                if( present(active_flags) ) then
                  allocate(arr(i)%active_flags, source=trim(active_flags))
                endif
                if( present(online) ) then
                  arr(i)%online = online
                endif
            end subroutine set

    end subroutine set_gui_params

    subroutine print_ui( self )
        class(simple_program), intent(in) :: self
        type(chash) :: ch
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name)
        call ch%push('descr_short', self%descr_short)
        call ch%push('descr_long',  self%descr_long)
        call ch%push('executable',  self%executable)
        call ch%print_key_val_pairs(logfhandle)
        call ch%kill
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash(self%img_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash(self%parm_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash(self%alt_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%srch_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%filt_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash(self%mask_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash(self%comp_ctrls)

        contains

            subroutine print_param_hash( arr )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        write(logfhandle,'(a,1x,i3)') '>>> PARAMETER #', i
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
                        call ch%print_key_val_pairs(logfhandle)
                        call ch%kill
                    end do
                endif
            end subroutine print_param_hash

    end subroutine print_ui

    subroutine print_ui_latex
        integer :: i
        type(simple_program), pointer :: ptr => null()
        do i=1,n_prg_ptrs
            ptr => prg_ptr_array(i)%ptr2prg
            write(logfhandle, '(a)') '\subsection{' // str2latex(ptr%name) // '}'
            write(logfhandle, '(a)') '\label{'      // ptr%name // '}'
            write(logfhandle, '(a)') '\prgname{' // str2latex(ptr%name) // '} ' // str2latex(ptr%descr_long // '. Executed by ' // ptr%executable)
            call ptr%print_cmdline_latex
            write(logfhandle, '(a)') ''
        end do
    end subroutine print_ui_latex

    subroutine print_cmdline( self )
        class(simple_program), intent(in) :: self
        write(logfhandle,'(a)') format_str('USAGE', C_UNDERLINED)
        write(logfhandle,'(a)') format_str('bash-3.2$ simple_exec prg=' //self%name // ' key1=val1 key2=val2 ...', C_ITALIC)
        write(logfhandle,'(a)') 'Required input parameters in ' // format_str('bold', C_BOLD) // ' (ensure terminal support)'
        if( allocated(self%img_ios) )    write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash(self%img_ios)
        if( allocated(self%parm_ios) )   write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash(self%parm_ios)
        if( allocated(self%alt_ios) )    write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash(self%alt_ios)
        if( allocated(self%srch_ctrls) ) write(logfhandle,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%srch_ctrls)
        if( allocated(self%filt_ctrls) ) write(logfhandle,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%filt_ctrls)
        if( allocated(self%mask_ctrls) ) write(logfhandle,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash(self%mask_ctrls)
        if( allocated(self%comp_ctrls) ) write(logfhandle,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash(self%comp_ctrls)
    end subroutine print_cmdline

    subroutine print_cmdline_latex( self )
        class(simple_program), intent(in) :: self
        logical     :: l_distr_exec
        l_distr_exec = self%executable .eq. 'simple_exec'
        write(logfhandle,'(a)') '\begin{Verbatim}[commandchars=+\[\],fontsize=\small,breaklines=true]'
        write(logfhandle,'(a)') '+underline[USAGE]'
        if( l_distr_exec )then
            write(logfhandle,'(a)') '+textit[bash-3.2$ simple_exec prg=' // self%name // ' key1=val1 key2=val2 ...]'
        else
            write(logfhandle,'(a)') '+textit[bash-3.2$ simple_exec prg='       // self%name // ' key1=val1 key2=val2 ...]'
        endif
        write(logfhandle,'(a)') 'Required input parameters in ' // '+textbf[bold]' // ' (ensure terminal support)'
        if( allocated(self%img_ios) )    write(logfhandle,'(a)') '+underline[IMAGE INPUT/OUTPUT]'
        call print_param_hash(self%img_ios,   latex=.true.)
        if( allocated(self%parm_ios) )   write(logfhandle,'(a)') '+underline[PARAMETER INPUT/OUTPUT]'
        call print_param_hash(self%parm_ios,  latex=.true.)
        if( allocated(self%alt_ios) )    write(logfhandle,'(a)') '+underline[ALTERNATIVE INPUTS]'
        call print_param_hash(self%alt_ios,   latex=.true.)
        if( allocated(self%srch_ctrls) ) write(logfhandle,'(a)') '+underline[SEARCH CONTROLS]'
        call print_param_hash(self%srch_ctrls, latex=.true.)
        if( allocated(self%filt_ctrls) ) write(logfhandle,'(a)') '+underline[FILTER CONTROLS]'
        call print_param_hash(self%filt_ctrls, latex=.true.)
        if( allocated(self%mask_ctrls) ) write(logfhandle,'(a)') '+underline[MASK CONTROLS]'
        call print_param_hash(self%mask_ctrls, latex=.true.)
        if( allocated(self%comp_ctrls) ) write(logfhandle,'(a)') '+underline[COMPUTER CONTROLS]'
        call print_param_hash(self%comp_ctrls, latex=.true.)
        write(logfhandle,'(a)') '\end{Verbatim}'
    end subroutine print_cmdline_latex

    ! supporting the print_cmdline routines (above)
    subroutine print_param_hash( arr, latex )
        type(simple_input_param), allocatable, intent(in) :: arr(:)
        logical, optional,                     intent(in) :: latex
        character(len=KEYLEN),    allocatable :: sorted_keys(:), rearranged_keys(:)
        logical,                  allocatable :: required(:)
        integer,                  allocatable :: inds(:)
        type(chash) :: ch
        integer     :: i, nparams, nreq, iopt
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
            call ch%print_key_val_pairs(logfhandle, sorted_keys, mask=required, latex=latex)
            call ch%kill
            deallocate(sorted_keys, required)
        endif
    end subroutine print_param_hash

    subroutine print_prg_descr_long( self )
        class(simple_program), intent(in) :: self
        write(logfhandle,'(a)') self%descr_long
    end subroutine print_prg_descr_long

    subroutine write_ui_json
        use json_module
        type(json_core)           :: json
        type(json_value), pointer :: program_entry, program, all_programs
        integer :: iprg
        ! JSON init
        call json%initialize()
        ! create array of program entries
        call json%create_array(all_programs, 'SIMPLE User Interface')
        do iprg=1,n_prg_ptrs
            call create_program_entry
            call json%add(all_programs, program_entry)
        end do
        ! write & clean
        call json%print(all_programs, 'simple_user_interface.json')
        if( json%failed() )then
            write(logfhandle,*) 'json input/output error for simple_user_interface'
            stop
        endif
        call json%destroy(all_programs)

        contains

            subroutine create_program_entry
                call json%create_object(program_entry,'')
                call json%create_object(program, trim(prg_ptr_array(iprg)%ptr2prg%name))
                call json%add(program_entry, program)
                ! program section
                call json%add(program, 'name',        prg_ptr_array(iprg)%ptr2prg%name)
                call json%add(program, 'descr_short', prg_ptr_array(iprg)%ptr2prg%descr_short)
                call json%add(program, 'descr_long',  prg_ptr_array(iprg)%ptr2prg%descr_long)
                call json%add(program, 'executable',  prg_ptr_array(iprg)%ptr2prg%executable)
                call json%add(program, 'advanced', prg_ptr_array(iprg)%ptr2prg%advanced)
                if( allocated(prg_ptr_array(iprg)%ptr2prg%gui_submenu_list)) then
                    call json%add(program, 'gui_submenu_list', prg_ptr_array(iprg)%ptr2prg%gui_submenu_list)
                endif
                ! all sections
                call create_section( 'image input/output',     prg_ptr_array(iprg)%ptr2prg%img_ios )
                call create_section( 'parameter input/output', prg_ptr_array(iprg)%ptr2prg%parm_ios )
                call create_section( 'alternative inputs',     prg_ptr_array(iprg)%ptr2prg%alt_ios )
                call create_section( 'search controls',        prg_ptr_array(iprg)%ptr2prg%srch_ctrls )
                call create_section( 'filter controls',        prg_ptr_array(iprg)%ptr2prg%filt_ctrls )
                call create_section( 'mask controls',          prg_ptr_array(iprg)%ptr2prg%mask_ctrls )
                call create_section( 'computer controls',      prg_ptr_array(iprg)%ptr2prg%comp_ctrls )
            end subroutine create_program_entry

            subroutine create_section( name, arr )
                character(len=*),          intent(in) :: name
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                type(json_value), pointer :: entry, section
                character(len=STDLEN)     :: options_str, before
                character(len=KEYLEN)     :: args(10)
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
                        if( allocated(arr(i)%gui_submenu) ) then
                          call json%add(entry, 'gui_submenu', trim(arr(i)%gui_submenu))
                        endif
                        if( allocated(arr(i)%exclusive_group) ) then
                          call json%add(entry, 'exclusive_group', trim(arr(i)%exclusive_group))
                        endif
                        if( allocated(arr(i)%active_flags) ) then
                          call json%add(entry, 'active_flags', trim(arr(i)%active_flags))
                        endif
                        call json%add(entry, 'advanced', arr(i)%advanced)
                        call json%add(entry, 'online', arr(i)%online)
                        param_is_multi  = trim(arr(i)%keytype).eq.'multi'
                        param_is_binary = trim(arr(i)%keytype).eq.'binary'
                        if( param_is_multi .or. param_is_binary )then
                            options_str = trim(arr(i)%descr_placeholder)
                            call split( options_str, '(', before )
                            call split( options_str, ')', before )
                            call parsestr(before, '|', args, nargs)
                            exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 2)
                            if( exception )then
                                write(logfhandle,*)'Poorly formatted options string for entry ', trim(arr(i)%key)
                                write(logfhandle,*)trim(arr(i)%descr_placeholder)
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
                call json%add(program_entry, section)
            end subroutine create_section

    end subroutine write_ui_json

    subroutine write2json( self )
        use json_module
        class(simple_program), intent(in) :: self
        type(json_core)           :: json
        type(json_value), pointer :: program_entry, program
        ! JSON init
        call json%initialize()
        call json%create_object(program_entry,'')
        call json%create_object(program, trim(self%name))
        call json%add(program_entry, program)
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
        call json%print(program_entry, trim(adjustl(self%name))//'.json')
        if( json%failed() )then
            write(logfhandle,*) 'json input/output error for program: ', trim(self%name)
            stop
        endif
        call json%destroy(program_entry)

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
                            exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 2)
                            if( exception )then
                                write(logfhandle,*)'Poorly formatted options string for entry ', trim(arr(i)%key)
                                write(logfhandle,*)trim(arr(i)%descr_placeholder)
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
                call json%add(program_entry, section)
            end subroutine create_section

    end subroutine write2json

    function get_name( self ) result( name )
        class(simple_program), intent(in) :: self
        character(len=:), allocatable :: name
        allocate(name, source=trim(self%name))
    end function get_name

    function get_executable( self ) result( name )
        class(simple_program), intent(in) :: self
        character(len=:), allocatable :: name
        allocate(name, source=trim(self%executable))
    end function get_executable

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
