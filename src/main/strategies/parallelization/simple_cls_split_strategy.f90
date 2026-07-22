module simple_cls_split_strategy
use simple_core_module_api
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder,                only: builder
use simple_parameters,             only: parameters
use simple_cmdline,                only: cmdline
use simple_defs_fname,             only: SRCHSPACE_MAP_FNAME
use simple_qsys_env,               only: qsys_env
use simple_sp_project,             only: sp_project
use simple_image,                  only: image
use simple_image_msk,              only: automask2D
use simple_classaverager,          only: transform_ptcls
use simple_srchspace_map2D_io,     only: write_srchspace_map2D
use simple_srch_sort_loc,          only: hpsort
implicit none

public :: cls_split_strategy
public :: cls_split_shmem_strategy
public :: cls_split_worker_strategy
public :: cls_split_master_strategy
public :: create_cls_split_strategy

private
#include "simple_local_flags.inc"

integer, parameter :: CLS_SPLIT_ICM_RANK_MAXITS = 16 
real,    parameter :: CLS_SPLIT_ICM_RANK_BETA_FRAC = 0.35       ! neighbor penalty: larger values discourage fragmented spectra such as keep/drop/keep.
real,    parameter :: CLS_SPLIT_ICM_RANK_COMPLEXITY_FRAC = 0.10 ! Late-rank complexity penalty: larger values make high-index eigenfeatures harder to retain.
real,    parameter :: CLS_SPLIT_ICM_RANK_LOWER_SEED_FRAC = 0.50 ! Lower-seed fraction of the eigengap upper bound; smaller values probe more compact latent spaces.

type, abstract :: cls_split_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type cls_split_strategy

type, extends(cls_split_strategy) :: cls_split_shmem_strategy
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type cls_split_shmem_strategy

type, extends(cls_split_strategy) :: cls_split_worker_strategy
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
end type cls_split_worker_strategy

type, extends(cls_split_strategy) :: cls_split_master_strategy
    type(qsys_env)                :: qenv
    type(chash)                   :: job_descr
    type(chash), allocatable      :: part_params(:)
    integer                       :: nparts_run = 1
    integer                       :: nthr_master = 1
contains
    procedure :: initialize   => master_initialize
    procedure :: execute      => master_execute
    procedure :: finalize_run => master_finalize_run
    procedure :: cleanup      => master_cleanup
end type cls_split_master_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: cls_split_strategy, parameters, builder, cmdline
        class(cls_split_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: cls_split_strategy, parameters, builder, cmdline
        class(cls_split_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: cls_split_strategy, parameters, builder, cmdline
        class(cls_split_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params)
        import :: cls_split_strategy, parameters
        class(cls_split_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    function create_cls_split_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(cls_split_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker, is_master
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        is_worker = cline%defined('part')
        is_master = (nparts > 1) .and. (.not. is_worker)
        if( is_master )then
            allocate(cls_split_master_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED CLS_SPLIT (MASTER)'
        else if( is_worker )then
            allocate(cls_split_worker_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> CLS_SPLIT (WORKER)'
        else
            allocate(cls_split_shmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> CLS_SPLIT (SHARED-MEMORY)'
        endif
    end function create_cls_split_strategy

    subroutine apply_defaults(cline)
        use simple_default_clines, only: set_automask2D_defaults
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')    ) call cline%set('neigs',    DIFFMAP_NEIGS_SCAN_DEFAULT)
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'diffusion_maps')
        if( .not. cline%defined('k_nn')     ) call cline%set('k_nn',      DIFFMAP_GRAPH_KNN_DEFAULT)
        if( .not. cline%defined('graph') )then
            if( cline%get_carg('oritype') == 'ptcl3D' )then
                call cline%set('graph', 'ori')
            else
                call cline%set('graph', 'euc')
            endif
        endif
        call set_automask2D_defaults(cline)
    end subroutine apply_defaults

    subroutine init_common(params, build, cline)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        logical :: do3d
        call apply_defaults(cline)
        do3d = .false.
        if( cline%defined('oritype') ) do3d = (cline%get_carg('oritype') .eq. 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=do3d)
        call validate_cls_split(params, cline)
    end subroutine init_common

    subroutine validate_cls_split(params, cline)
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        logical :: l_fixed_nsubcls
        l_fixed_nsubcls = cline%defined('ncls') .and. params%ncls > 1
        select case(trim(params%oritype))
            case('ptcl2D','ptcl3D')
            case DEFAULT
                THROW_HARD('cls_split supports oritype=ptcl2D|ptcl3D only')
        end select
        select case(trim(params%pca_mode))
            case('diffusion_maps','kpca')
            case DEFAULT
                THROW_HARD('cls_split supports pca_mode=diffusion_maps|kpca only')
        end select
        select case(trim(lowercase(params%graph)))
            case('euc','ori')
            case DEFAULT
                THROW_HARD('cls_split supports graph=euc|ori only')
        end select
        if( params%k_nn < 1 ) THROW_HARD('cls_split requires k_nn >= 1')
        if( params%neigs < 0 ) THROW_HARD('cls_split requires neigs >= 0; use neigs=0 for auto scan')
        if( cline%defined('ncls') )then
            if( params%ncls < 0 ) THROW_HARD('cls_split requires ncls >= 0; use ncls=0 for auto')
            if( params%ncls == 1 ) THROW_HARD('cls_split ncls=1 is invalid; use ncls=0 for auto or ncls>=2')
        endif
        if( .not. l_fixed_nsubcls )then
            if( params%nsubcls_min < 2 ) THROW_HARD('cls_split requires nsubcls_min >= 2')
            if( params%nsubcls_max < params%nsubcls_min ) THROW_HARD('cls_split requires nsubcls_max >= nsubcls_min')
        endif
    end subroutine validate_cls_split

    subroutine shmem_initialize(self, params, build, cline)
        class(cls_split_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call init_common(params, build, cline)
    end subroutine shmem_initialize

    subroutine shmem_execute(self, params, build, cline)
        class(cls_split_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_local_split(params, build, cline, spproj, part=0, l_write_project=.true.)
        call spproj%kill
    end subroutine shmem_execute

    subroutine shmem_finalize_run(self, params, build, cline)
        class(cls_split_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup(self, params)
        class(cls_split_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
    end subroutine shmem_cleanup

    subroutine worker_initialize(self, params, build, cline)
        class(cls_split_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call init_common(params, build, cline)
        if( .not. cline%defined('part') ) THROW_HARD('PART must be defined for cls_split worker execution')
        if( .not. cline%defined('class_assignment') ) THROW_HARD('CLASS_ASSIGNMENT must be defined for cls_split worker execution')
    end subroutine worker_initialize

    subroutine worker_execute(self, params, build, cline)
        use simple_qsys_funs, only: qsys_job_finished
        class(cls_split_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_local_split(params, build, cline, spproj, part=params%part, l_write_project=.false.)
        call qsys_job_finished(params, string('simple_cls_split_strategy :: worker_execute'))
        call spproj%kill
    end subroutine worker_execute

    subroutine worker_finalize_run(self, params, build, cline)
        use simple_qsys_funs, only: qsys_job_finished
        class(cls_split_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_cluster2D :: exec_cls_split'))
    end subroutine worker_finalize_run

    subroutine worker_cleanup(self, params)
        class(cls_split_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
    end subroutine worker_cleanup

    subroutine master_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(cls_split_master_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        integer, allocatable :: cls_inds(:), cls_pops(:)
        call init_common(params, build, cline)
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        self%nparts_run = min(max(1, params%nparts), size(cls_inds))
        write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls split worker planning: requested_nparts=', params%nparts, &
            ' eligible_classes=', size(cls_inds), ' effective_nparts=', self%nparts_run
        call flush(logfhandle)
        if( self%nparts_run < params%nparts )then
            write(logfhandle,'(A)') 'Cls split note: reducing worker count because classes are the scheduling unit.'
            call flush(logfhandle)
        endif
        call set_master_num_threads(self%nthr_master, string('CLS_SPLIT'))
        call prepare_class_partitions(self, params, cline, cls_inds, cls_pops)
        call self%qenv%new(params, self%nparts_run, numlen=params%numlen)
        ! Keep distributed workers inside the master's execution directory so
        ! JOB_FINISHED flags and part outputs land where the master is watching.
        call cline%set('mkdir', 'no')
        call cline%gen_job_descr(self%job_descr)
        call self%job_descr%set('mkdir', 'no')
        call self%job_descr%set('nparts', int2str(self%nparts_run))
        call self%job_descr%set('numlen', int2str(params%numlen))
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
    end subroutine master_initialize

    subroutine master_execute(self, params, build, cline)
        class(cls_split_master_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, part_params=self%part_params, &
                                                     array=L_USE_SLURM_ARR, extra_params=params)
        call merge_worker_outputs(params, self%nparts_run)
    end subroutine master_execute

    subroutine master_finalize_run(self, params, build, cline)
        class(cls_split_master_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
    end subroutine master_finalize_run

    subroutine master_cleanup(self, params)
        use simple_qsys_funs, only: qsys_cleanup
        class(cls_split_master_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        integer :: ipart
        call self%qenv%kill
        call qsys_cleanup(params)
        if( allocated(self%part_params) )then
            do ipart = 1, size(self%part_params)
                call self%part_params(ipart)%kill
            end do
            deallocate(self%part_params)
        endif
        call self%job_descr%kill
    end subroutine master_cleanup

    subroutine prepare_class_partitions(self, params, cline, cls_inds, cls_pops)
        class(cls_split_master_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        integer,                          intent(in)    :: cls_inds(:), cls_pops(:)
        integer :: order(size(cls_inds))
        integer, allocatable :: part_counts(:), part_cls(:,:)
        integer(kind=8), allocatable :: part_weights(:)
        integer :: ncls, ipart, iord, icls, lightest
        type(string) :: fname
        ncls = size(cls_inds)
        allocate(part_counts(self%nparts_run), part_cls(ncls, self%nparts_run), part_weights(self%nparts_run))
        order = [(icls, icls=1,ncls)]
        call sort_order_by_weight_desc(order, cls_pops)
        part_counts  = 0
        part_weights = 0_8
        part_cls     = 0
        do iord = 1, ncls
            icls = order(iord)
            lightest = minloc(part_weights, dim=1)
            part_counts(lightest) = part_counts(lightest) + 1
            part_cls(part_counts(lightest), lightest) = cls_inds(icls)
            part_weights(lightest) = part_weights(lightest) + int(cls_pops(icls), kind=8)
        end do
        allocate(self%part_params(self%nparts_run))
        do ipart = 1, self%nparts_run
            if( part_counts(ipart) > 1 ) call hpsort(part_cls(1:part_counts(ipart), ipart))
            fname = string('cls_split_classes_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call arr2txtfile(part_cls(1:part_counts(ipart), ipart), fname)
            call self%part_params(ipart)%new(1)
            call self%part_params(ipart)%set('class_assignment', fname%to_char())
            call fname%kill
        end do
        deallocate(part_counts, part_cls, part_weights)
    end subroutine prepare_class_partitions

    subroutine collect_split_classes(cline, params, build, cls_inds, cls_pops)
        class(cmdline),       intent(inout) :: cline
        type(parameters),     intent(inout) :: params
        type(builder),        intent(inout) :: build
        integer, allocatable, intent(out)   :: cls_inds(:), cls_pops(:)
        integer, allocatable :: pinds(:), assigned_classes(:)
        logical, allocatable :: keep_mask(:)
        type(string) :: label
        integer :: i
        call determine_split_label(params, build, label)
        cls_inds = build%spproj_field%get_label_inds(label%to_char())
        if( cline%defined('class') ) cls_inds = pack(cls_inds, mask=(cls_inds == params%class))
        if( cline%defined('class_assignment') )then
            call read_int_file(cline%get_carg('class_assignment'), assigned_classes)
            allocate(keep_mask(size(cls_inds)), source=.false.)
            do i = 1, size(assigned_classes)
                keep_mask = keep_mask .or. (cls_inds == assigned_classes(i))
            end do
            cls_inds = pack(cls_inds, mask=keep_mask)
            deallocate(keep_mask)
            if( allocated(assigned_classes) ) deallocate(assigned_classes)
        endif
        if( size(cls_inds) < 1 ) THROW_HARD('No classes selected for cls_split')
        allocate(cls_pops(size(cls_inds)), source=0)
        do i = 1, size(cls_inds)
            call build%spproj_field%get_pinds(cls_inds(i), label%to_char(), pinds)
            if( allocated(pinds) )then
                cls_pops(i) = size(pinds)
                deallocate(pinds)
            endif
        end do
        cls_inds = pack(cls_inds, mask=cls_pops > 2)
        cls_pops = pack(cls_pops, mask=cls_pops > 2)
        if( size(cls_inds) < 1 ) THROW_HARD('No classes with enough particles to split')
        call label%kill
    end subroutine collect_split_classes

    subroutine determine_split_label(params, build, label)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(string),     intent(out)   :: label
        select case(trim(params%oritype))
            case('ptcl2D')
                label = 'class'
            case('ptcl3D')
                label = 'proj'
                call build%spproj_field%proj2class
            case DEFAULT
                THROW_HARD('cls_split only supports ORITYPE=ptcl2D or ptcl3D')
        end select
    end subroutine determine_split_label

    subroutine determine_phase_flip(spproj, params, l_phflip)
        type(sp_project), intent(inout) :: spproj
        type(parameters), intent(inout) :: params
        logical,          intent(out)   :: l_phflip
        l_phflip = .false.
        select case( spproj%get_ctfflag_type(params%oritype) )
            case(CTFFLAG_NO)
                THROW_WARN('No CTF information could be found, phase flipping is deactivated')
            case(CTFFLAG_FLIP)
                THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
            case(CTFFLAG_YES)
                l_phflip = .true.
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG')
        end select
    end subroutine determine_phase_flip

    subroutine run_local_split(params, build, cline, spproj, part, l_write_project)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: part
        logical,          intent(in)    :: l_write_project
        type(string) :: label, map_fname, assign_fname
        integer, allocatable :: cls_inds(:), cls_pops(:), pinds(:), labels(:), new_class(:), new_parent(:), parent_of_subcls(:), pop_of_subcls(:)
        integer :: i, j, k, iglob, nsplit, funit_map, funit_assign
        logical :: l_phflip, l_fixed_nsubcls, l_mskdiam_override
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        call determine_split_label(params, build, label)
        l_fixed_nsubcls = cline%defined('ncls') .and. params%ncls > 1
        l_mskdiam_override = cline%defined('mskdiam')
        if( trim(params%gen_model) == 'yes' )then
            THROW_WARN('gen_model is ignored by cls_split; use denoise_project for generated/denoised particles')
        endif
        call determine_phase_flip(spproj, params, l_phflip)
        if( l_write_project )then
            map_fname = string('cls_split_class_map.txt')
            open(newunit=funit_map, file=map_fname%to_char(), status='replace', action='write')
            write(funit_map,'(A)') '# global_subclass parent_class local_subclass pop'
            select case(trim(params%oritype))
                case('ptcl2D')
                    allocate(new_class(spproj%os_ptcl2D%get_noris()), new_parent(spproj%os_ptcl2D%get_noris()), source=0)
                case('ptcl3D')
                    allocate(new_class(spproj%os_ptcl3D%get_noris()), new_parent(spproj%os_ptcl3D%get_noris()), source=0)
            end select
            allocate(parent_of_subcls(sum(cls_pops)), pop_of_subcls(sum(cls_pops)), source=0)
        else
            map_fname    = string('cls_split_class_map_part')//int2str_pad(part, params%numlen)//TXT_EXT
            assign_fname = string('cls_split_assignments_part')//int2str_pad(part, params%numlen)//TXT_EXT
            open(newunit=funit_map,    file=map_fname%to_char(),    status='replace', action='write')
            open(newunit=funit_assign, file=assign_fname%to_char(), status='replace', action='write')
            write(funit_map,'(A)')    '# local_subclass_row parent_class local_subclass pop'
            write(funit_assign,'(A)') '# particle_index parent_class local_subclass'
        endif
        iglob = 0
        do i = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8)') 'Cls split starting: class=', cls_inds(i), ' nptcls=', cls_pops(i)
            call flush(logfhandle)
            call split_one_parent_class(params, build, spproj, cls_inds(i), l_phflip, l_fixed_nsubcls, &
                                        l_mskdiam_override, nsplit, pinds, labels)
            if( nsplit < 1 ) cycle
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls split summary: class=', cls_inds(i), ' nptcls=', size(labels), ' nsubcls=', nsplit
            call flush(logfhandle)
            do j = 1, nsplit
                iglob = iglob + 1
                if( l_write_project )then
                    parent_of_subcls(iglob) = cls_inds(i)
                    pop_of_subcls(iglob)    = count(labels == j)
                    write(funit_map,'(I8,1X,I8,1X,I8,1X,I8)') iglob, cls_inds(i), j, count(labels == j)
                else
                    write(funit_map,'(I8,1X,I8,1X,I8,1X,I8)') iglob, cls_inds(i), j, count(labels == j)
                endif
            end do
            if( l_write_project )then
                do k = 1, size(labels)
                    if( labels(k) <= 0 ) cycle
                    new_class(pinds(k))  = iglob - nsplit + labels(k)
                    new_parent(pinds(k)) = cls_inds(i)
                end do
            else
                do k = 1, size(labels)
                    if( labels(k) <= 0 ) cycle
                    write(funit_assign,'(I10,1X,I10,1X,I10)') pinds(k), cls_inds(i), labels(k)
                end do
            endif
            if( allocated(pinds) ) deallocate(pinds)
            if( allocated(labels) ) deallocate(labels)
        end do
        close(funit_map)
        if( .not. l_write_project ) close(funit_assign)
        if( l_write_project )then
            call apply_split_project_updates(spproj, params, iglob, new_class, new_parent, parent_of_subcls(1:iglob), pop_of_subcls(1:iglob))
            deallocate(new_class, new_parent, parent_of_subcls, pop_of_subcls)
        endif
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
        call label%kill
        call map_fname%kill
        if( assign_fname%is_allocated() ) call assign_fname%kill
    end subroutine run_local_split

    subroutine split_one_parent_class(params, build, spproj, cls_id, l_phflip, l_fixed_nsubcls, &
                                      l_mskdiam_override, nsplit, pinds, labels)
        use simple_diff_map_graphs,  only: diffmap_graph
        use simple_imgarr_utils,     only: dealloc_imgarr, copy_imgarr
        use simple_imgproc,          only: make_pcavecs
        use simple_clustering_utils, only: cluster_dmat, silhouette_score
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id
        logical,                  intent(in)    :: l_phflip, l_fixed_nsubcls, l_mskdiam_override
        integer,                  intent(out)   :: nsplit
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        type(parameters) :: params_mask
        type(image), allocatable :: imgs(:), imgs_ppca(:), class_mask(:)
        type(image) :: cavg_raw
        real, allocatable :: avg(:), pcavecs(:,:), coords(:,:), eigvals(:), dmat(:,:)
        real, allocatable :: class_diams(:), class_shifts(:,:)
        type(diffmap_graph) :: split_graph
        integer, allocatable :: i_medoids(:)
        integer :: nptcls, npix, j, k, class_ldim(3)
        real    :: class_moldiam, class_mskdiam, class_mskrad, dval, sdev_noise
        nsplit = 0
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(labels) ) deallocate(labels)
        call transform_ptcls(params, build, spproj, params%oritype, cls_id, imgs, pinds, phflip=l_phflip, cavg=cavg_raw)
        if( .not. allocated(imgs) )then
            write(logfhandle,'(A,I8)') 'Diffusion class splitting warning: no transformed images returned for parent class ', cls_id
            call flush(logfhandle)
            return
        endif
        nptcls = size(imgs)
        if( nptcls < 3 )then
            call dealloc_imgarr(imgs)
            if( allocated(pinds) ) deallocate(pinds)
            return
        endif
        imgs_ppca = copy_imgarr(imgs)
        call imgs_ppca(1)%memoize_mask_coords()
        allocate(class_mask(1))
        call class_mask(1)%copy(cavg_raw)
        params_mask%ngrow   = params%ngrow
        params_mask%winsz   = params%winsz
        params_mask%edge    = params%edge
        params_mask%amsklp  = params%amsklp
        params_mask%automsk = params%automsk
        params_mask%part    = params%part
        class_ldim        = cavg_raw%get_ldim()
        params_mask%box   = class_ldim(1)
        params_mask%smpd  = cavg_raw%get_smpd()
        params_mask%msk   = real(class_ldim(1) / 2) - COSMSKHALFWIDTH
        if( l_mskdiam_override )then
            class_mskdiam = params%mskdiam
        else
            call automask2D(params_mask, class_mask, params_mask%ngrow, nint(params_mask%winsz), params_mask%edge, class_diams, class_shifts)
            class_moldiam = params_mask%smpd * real(min(round2even(class_diams(1) / params_mask%smpd + 2. * COSMSKHALFWIDTH), class_ldim(1)))
            class_mskdiam = class_moldiam * MSK_EXP_FAC
        endif
        class_mskrad  = min(real(class_ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * class_mskdiam / params_mask%smpd)
        if( l_mskdiam_override )then
            write(logfhandle,'(A,I8,A,F8.2,A,F8.2)') 'Cls split mask update: class=', cls_id, &
                ' input_mskdiam=', class_mskdiam, ' mskrad=', class_mskrad
        else
            write(logfhandle,'(A,I8,A,F8.2,A,F8.2,A,F8.2)') 'Cls split mask update: class=', cls_id, &
                ' automask_diam=', class_diams(1), ' mskdiam=', class_mskdiam, ' mskrad=', class_mskrad
        endif
        call flush(logfhandle)
        do j = 1, nptcls
            call imgs_ppca(j)%norm_noise(build%lmsk, sdev_noise)
            call imgs_ppca(j)%mask2D_softavg(class_mskrad)
        end do
        call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
        call make_split_embedding(params, spproj, pinds, cls_id, nptcls, npix, pcavecs, imgs_ppca, coords, eigvals, &
                                  split_graph)
        call sanitize_embedding_coords(trim(params%pca_mode), cls_id, coords)
        allocate(dmat(nptcls,nptcls), source=0.)
        do j = 1, nptcls - 1
            do k = j + 1, nptcls
                dval = euclid(coords(:,j), coords(:,k))
                dmat(j,k) = dval
                dmat(k,j) = dval
            end do
        end do
        call sanitize_distance_matrix(trim(params%pca_mode), cls_id, dmat)
        if( l_fixed_nsubcls )then
            nsplit = params%ncls
            call cluster_dmat(dmat, 'kmed', nsplit, i_medoids, labels)
            write(logfhandle,'(A,I8,A,I8,A,F10.5)') 'Cls split fixed ncls silhouette: class=', cls_id, &
                ' selected=', nsplit, ' score=', silhouette_score(labels, dmat)
            call flush(logfhandle)
        else
            call select_kmedoids_by_silhouette(dmat, cls_id, params%nsubcls_min, params%nsubcls_max, &
                                               nsplit, i_medoids, labels)
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8)') 'Cls split auto ncls: class=', cls_id, &
                ' size=', nptcls, ' trial_min=', params%nsubcls_min, &
                ' trial_max=', params%nsubcls_max, ' selected=', nsplit
            call flush(logfhandle)
        endif
        call cavg_raw%kill
        if( allocated(coords)       ) deallocate(coords)
        if( allocated(eigvals)      ) deallocate(eigvals)
        if( allocated(dmat)         ) deallocate(dmat)
        call split_graph%kill()
        if( allocated(i_medoids)    ) deallocate(i_medoids)
        if( allocated(avg)          ) deallocate(avg)
        if( allocated(pcavecs)      ) deallocate(pcavecs)
        if( allocated(class_diams)  ) deallocate(class_diams)
        if( allocated(class_shifts) ) deallocate(class_shifts)
        if( allocated(class_mask)   ) call dealloc_imgarr(class_mask)
        if( allocated(imgs_ppca)    ) call dealloc_imgarr(imgs_ppca)
        if( allocated(imgs)         ) call dealloc_imgarr(imgs)
    end subroutine split_one_parent_class

    subroutine select_kmedoids_by_silhouette(dmat, cls_id, k_min_in, k_max_in, nsplit, i_medoids, labels)
        use simple_clustering_utils, only: cluster_dmat, silhouette_score
        real,                 intent(in)    :: dmat(:,:)
        integer,              intent(in)    :: cls_id, k_min_in, k_max_in
        integer,              intent(out)   :: nsplit
        integer, allocatable, intent(inout) :: i_medoids(:), labels(:)
        integer, allocatable :: trial_medoids(:), trial_labels(:), best_medoids(:), best_labels(:)
        integer :: nptcls, k_min, k_max, k_trial, ntrial, best_k, icls
        real    :: score, best_score
        logical :: l_all_populated
        nptcls = size(dmat, dim=1)
        k_min  = max(2, min(k_min_in, max(2, nptcls - 1)))
        k_max  = max(k_min, min(k_max_in, max(2, nptcls - 1)))
        best_score = -huge(best_score)
        best_k     = k_min
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8)') 'Cls split silhouette k search: class=', cls_id, &
            ' size=', nptcls, ' k_min=', k_min, ' k_max=', k_max
        call flush(logfhandle)
        do k_trial = k_min, k_max
            ntrial = k_trial
            call cluster_dmat(dmat, 'kmed', ntrial, trial_medoids, trial_labels)
            l_all_populated = .true.
            do icls = 1, ntrial
                if( count(trial_labels == icls) == 0 )then
                    l_all_populated = .false.
                    exit
                endif
            end do
            if( l_all_populated )then
                score = silhouette_score(trial_labels, dmat)
            else
                score = -huge(score)
            endif
            write(logfhandle,'(A,I8,A,I8,A,F10.5)') 'Cls split silhouette k trial: class=', cls_id, &
                ' k=', k_trial, ' score=', score
            call flush(logfhandle)
            if( ieee_is_finite(score) .and. score > best_score + 1.e-6 )then
                best_score = score
                best_k     = ntrial
                if( allocated(best_medoids) ) deallocate(best_medoids)
                if( allocated(best_labels)  ) deallocate(best_labels)
                allocate(best_medoids(size(trial_medoids)), source=trial_medoids)
                allocate(best_labels(size(trial_labels)), source=trial_labels)
                deallocate(trial_medoids, trial_labels)
            else
                if( allocated(trial_medoids) ) deallocate(trial_medoids)
                if( allocated(trial_labels)  ) deallocate(trial_labels)
            endif
        end do
        if( .not. allocated(best_labels) )then
            best_k = k_min
            call cluster_dmat(dmat, 'kmed', best_k, best_medoids, best_labels)
            best_score = silhouette_score(best_labels, dmat) / real(best_k)
        endif
        if( allocated(i_medoids) ) deallocate(i_medoids)
        if( allocated(labels)    ) deallocate(labels)
        allocate(i_medoids(size(best_medoids)), source=best_medoids)
        allocate(labels(size(best_labels)), source=best_labels)
        deallocate(best_medoids, best_labels)
        nsplit = best_k
        write(logfhandle,'(A,I8,A,I8,A,F10.5)') 'Cls split silhouette k selected: class=', cls_id, &
            ' k=', nsplit, ' score=', best_score
        call flush(logfhandle)
    end subroutine select_kmedoids_by_silhouette

    subroutine make_split_embedding(params, spproj, pinds, cls_id, nptcls, npix, pcavecs, imgs_ppca, coords, eigvals, &
                                    graph)
        use simple_diff_map_graphs, only: diffmap_graph, build_cls_split_graph
        use simple_diffusion_maps,  only: embed_graph
        use simple_kpca_svd,        only: kpca_svd
        type(parameters),  intent(inout) :: params
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: pinds(:), cls_id, nptcls, npix
        real,              intent(in)    :: pcavecs(npix,nptcls)
        type(image),       intent(inout) :: imgs_ppca(:)
        real, allocatable, intent(out)   :: coords(:,:)
        real, allocatable, intent(out)   :: eigvals(:)
        type(diffmap_graph), intent(inout) :: graph
        type(kpca_svd) :: kpca_model
        real, allocatable :: feat(:)
        integer :: neigs, neigs_scan, neigs_used, i
        call graph%kill()
        if( params%neigs > 0 )then
            neigs_scan = params%neigs
        else
            neigs_scan = cls_split_auto_neigs_scan(trim(params%pca_mode), nptcls)
        endif
        select case(trim(params%pca_mode))
            case('diffusion_maps')
                neigs = min(max(neigs_scan, 0), max(nptcls-2, 1))
            case DEFAULT
                neigs = min(max(neigs_scan, 1), max(nptcls-1, 1))
        end select
        select case(trim(params%pca_mode))
            case('diffusion_maps')
                call build_cls_split_graph(params, spproj, pinds, pcavecs, graph)
                call embed_graph(graph, neigs, coords, eigvals)
            case('kpca')
                call kpca_model%new(nptcls, npix, neigs)
                call kpca_model%set_params(params%nthr, params%kpca_ker, params%kpca_backend, params%kpca_nystrom_npts, &
                    params%kpca_rbf_gamma, params%kpca_nystrom_local_nbrs, params%kpca_cosine_weight_power)
                call kpca_model%master(pcavecs, 15)
                allocate(coords(neigs,nptcls), source=0.)
                do i = 1, nptcls
                    feat = kpca_model%get_feat(i)
                    coords(:,i) = feat
                    if( allocated(feat) ) deallocate(feat)
                end do
                eigvals = kpca_model%get_eigvals()
                call kpca_model%kill
            case DEFAULT
                THROW_HARD('unsupported pca_mode for class splitting; use diffusion_maps or kpca')
        end select
        if( allocated(eigvals) )then
            neigs_used = select_neigs_icm(eigvals, trim(params%pca_mode), size(coords,1), cls_id)
            call trim_split_embedding(coords, eigvals, neigs_used)
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8)') 'Cls split selected neigs: mode=', trim(params%pca_mode), &
                ' class=', cls_id, ' size=', nptcls, ' scan=', neigs, ' selected=', size(coords,1)
            call flush(logfhandle)
        endif
        if( graph%n > 0 )then
            write(logfhandle,'(A,A,A,A,A,I8,A,I8,A,I8,A,I8)') 'Cls split embedding: mode=', &
                trim(params%pca_mode), ' metric=', trim(graph%metric), &
                ' class=', cls_id, ' size=', nptcls, ' dims=', size(coords,1), ' nnz=', graph%nnz
        else
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8)') 'Cls split embedding: mode=', trim(params%pca_mode), &
                ' class=', cls_id, ' size=', nptcls, ' dims=', size(coords,1)
        endif
        call flush(logfhandle)
    end subroutine make_split_embedding

    integer function cls_split_auto_neigs_scan(mode, nptcls) result(neigs_scan)
        character(len=*), intent(in) :: mode
        integer,          intent(in) :: nptcls
        select case(trim(mode))
            case('diffusion_maps')
                neigs_scan = min(DIFFMAP_NEIGS_AUTO_SCAN_MAX, max(1, nptcls - 2))
                if( nptcls > 3 ) neigs_scan = max(2, neigs_scan)
            case DEFAULT
                neigs_scan = min(DIFFMAP_NEIGS_AUTO_SCAN_MAX, max(1, nptcls - 1))
        end select
    end function cls_split_auto_neigs_scan

    integer function select_neigs_icm(eigvals, mode, max_neigs, cls_id) result(nkeep)
        real,             intent(in) :: eigvals(:)
        character(len=*), intent(in) :: mode
        integer,          intent(in) :: max_neigs, cls_id
        real,    allocatable :: spec(:)
        real :: smin, smax, delta, beta, complexity, score, best_score
        integer :: i, n, nmin_rank, upper_rank, lower_rank, best_seed, k_trial
        integer :: seeds(2)
        character(len=12) :: seed_names(2)
        n = min(size(eigvals), max(1, max_neigs))
        nmin_rank = cls_split_min_neigs(trim(mode), n)
        if( n <= 0 )then
            nkeep = 1
            return
        endif
        allocate(spec(n), source=0.)
        do i = 1, n
            if( ieee_is_finite(eigvals(i)) .and. eigvals(i) > DTINY )then
                spec(i) = log(eigvals(i))
            else
                spec(i) = log(DTINY)
            endif
        end do
        smin  = minval(spec)
        smax  = maxval(spec)
        delta = smax - smin
        if( .not. ieee_is_finite(delta) .or. delta <= 1.e-6 )then
            nkeep = nmin_rank
            deallocate(spec)
            return
        endif
        spec = (spec - smin) / delta
        upper_rank = cls_split_initial_neigs_from_gap(spec, nmin_rank)
        upper_rank = max(nmin_rank, min(n, upper_rank))
        lower_rank = nint(CLS_SPLIT_ICM_RANK_LOWER_SEED_FRAC * real(upper_rank))
        lower_rank = max(nmin_rank, min(upper_rank, lower_rank))
        seeds(1) = upper_rank
        seeds(2) = lower_rank
        seed_names = [character(len=12) :: 'upper_gap', 'fraction']
        beta       = CLS_SPLIT_ICM_RANK_BETA_FRAC
        complexity = CLS_SPLIT_ICM_RANK_COMPLEXITY_FRAC
        write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8,A,I8)') 'Cls split ICM neigs seeds: mode=', trim(mode), &
            ' class=', cls_id, ' scan=', n, ' min=', nmin_rank, ' upper=', upper_rank, &
            ' lower_fraction=', lower_rank
        call flush(logfhandle)
        best_score = huge(best_score)
        nkeep = upper_rank
        best_seed = 1
        do i = 1, size(seeds)
            call run_icm_rank_seed(spec, trim(mode), cls_id, trim(seed_names(i)), seeds(i), nmin_rank, upper_rank, beta, complexity, &
                k_trial, score)
            if( score < best_score )then
                best_score = score
                nkeep = k_trial
                best_seed = i
            endif
        end do
        write(logfhandle,'(A,A,A,I8,A,A,A,I8,A,F12.5)') 'Cls split ICM neigs selected: mode=', trim(mode), &
            ' class=', cls_id, ' seed=', trim(seed_names(best_seed)), ' keep=', nkeep, ' score=', best_score
        call flush(logfhandle)
        deallocate(spec)
    end function select_neigs_icm

    subroutine run_icm_rank_seed(spec, mode, cls_id, seed_name, seed_rank, nmin_rank, nmax_rank, beta, alpha, nkeep, score)
        real,             intent(in)  :: spec(:), beta, alpha
        character(len=*), intent(in)  :: mode, seed_name
        integer,          intent(in)  :: cls_id, seed_rank, nmin_rank, nmax_rank
        integer,          intent(out) :: nkeep
        real,             intent(out) :: score
        integer, allocatable :: labels(:), prev_labels(:)
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: iter, i, n, nchanged, maxits
        n = size(spec)
        allocate(labels(n), prev_labels(n), source=0)
        labels = 0
        labels(1:max(nmin_rank, min(nmax_rank, seed_rank))) = 1
        if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
        write(logfhandle,'(A,A,A,I8,A,A,A,I8,A,I8,A,I8)') 'Cls split ICM neigs init: mode=', trim(mode), &
            ' class=', cls_id, ' seed=', trim(seed_name), ' init=', cls_split_rank_prefix(labels, nmin_rank, nmax_rank), &
            ' min=', nmin_rank, ' max=', nmax_rank
        call flush(logfhandle)
        maxits = max(CLS_SPLIT_ICM_RANK_MAXITS, n)
        do iter = 1, maxits
            prev_labels = labels
            call estimate_icm_rank_stats(spec, prev_labels, mu_drop, mu_keep, var_drop, var_keep)
            do i = 1, nmax_rank
                call update_icm_rank_site(spec, prev_labels, i, beta, alpha, nmin_rank, nmax_rank, &
                    mu_drop, mu_keep, var_drop, var_keep, labels(i))
            end do
            labels(1:nmin_rank) = 1
            if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
            nchanged = count(labels /= prev_labels)
            score = score_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
            write(logfhandle,'(A,A,A,I8,A,A,A,I8,A,I8,A,I8,A,F12.5)') 'Cls split ICM neigs iter: mode=', trim(mode), &
                ' class=', cls_id, ' seed=', trim(seed_name), ' iter=', iter, ' changed=', nchanged, &
                ' keep=', cls_split_rank_prefix(labels, nmin_rank, nmax_rank), ' score=', score
            call flush(logfhandle)
            if( nchanged == 0 ) exit
        end do
        nkeep = cls_split_rank_prefix(labels, nmin_rank, nmax_rank)
        score = score_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
        deallocate(labels, prev_labels)
    end subroutine run_icm_rank_seed

    integer function cls_split_rank_prefix(labels, nmin_rank, nmax_rank) result(nkeep)
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        integer :: i
        nkeep = 0
        do i = 1, min(size(labels), nmax_rank)
            if( labels(i) == 1 ) nkeep = i
        end do
        nkeep = max(nmin_rank, min(nmax_rank, nkeep))
    end function cls_split_rank_prefix

    integer function cls_split_min_neigs(mode, max_neigs) result(nmin_rank)
        character(len=*), intent(in) :: mode
        integer,          intent(in) :: max_neigs
        nmin_rank = 1
        if( trim(mode) .eq. 'diffusion_maps' .and. max_neigs >= 2 ) nmin_rank = 2
    end function cls_split_min_neigs

    integer function cls_split_initial_neigs_from_gap(spec, nmin_rank) result(nkeep)
        real,    intent(in) :: spec(:)
        integer, intent(in) :: nmin_rank
        real :: best_gap, gap
        integer :: i
        if( size(spec) <= 1 )then
            nkeep = 1
            return
        endif
        best_gap = -huge(best_gap)
        nkeep = min(max(nmin_rank, 2), size(spec))
        do i = 1, size(spec) - 1
            gap = spec(i) - spec(i+1)
            if( gap > best_gap )then
                best_gap = gap
                nkeep = i
            endif
        end do
        nkeep = max(nmin_rank, min(size(spec), nkeep))
    end function cls_split_initial_neigs_from_gap

    subroutine estimate_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
        real,    intent(in)  :: spec(:)
        integer, intent(in)  :: labels(:)
        real,    intent(out) :: mu_drop, mu_keep, var_drop, var_keep
        integer :: ndrop, nkeep, i
        ndrop = count(labels == 0)
        nkeep = count(labels == 1)
        if( ndrop > 0 )then
            mu_drop = sum(spec, mask=(labels == 0)) / real(ndrop)
        else
            mu_drop = spec(size(spec))
        endif
        if( nkeep > 0 )then
            mu_keep = sum(spec, mask=(labels == 1)) / real(nkeep)
        else
            mu_keep = spec(1)
        endif
        var_drop = 0.
        var_keep = 0.
        do i = 1, size(spec)
            if( labels(i) == 0 )then
                var_drop = var_drop + (spec(i) - mu_drop)**2
            else
                var_keep = var_keep + (spec(i) - mu_keep)**2
            endif
        end do
        var_drop = max(var_drop / real(max(1, ndrop)), 1.e-3)
        var_keep = max(var_keep / real(max(1, nkeep)), 1.e-3)
    end subroutine estimate_icm_rank_stats

    subroutine update_icm_rank_site(spec, labels, ind, beta, alpha, nmin_rank, nmax_rank, &
        mu_drop, mu_keep, var_drop, var_keep, label_new)
        real,    intent(in)  :: spec(:), beta, alpha, mu_drop, mu_keep, var_drop, var_keep
        integer, intent(in)  :: labels(:), ind, nmin_rank, nmax_rank
        integer, intent(out) :: label_new
        real :: cost_drop, cost_keep
        integer :: kfree
        cost_drop = (spec(ind) - mu_drop)**2 / var_drop
        cost_keep = (spec(ind) - mu_keep)**2 / var_keep
        kfree = min(max(nmin_rank, 4), nmax_rank)
        cost_keep = cost_keep + alpha * real(max(0, ind - kfree)) / real(max(1, nmax_rank - kfree))
        if( ind > 1 )then
            if( labels(ind-1) /= 0 ) cost_drop = cost_drop + beta
            if( labels(ind-1) /= 1 ) cost_keep = cost_keep + beta
        endif
        if( ind < size(labels) )then
            if( labels(ind+1) /= 0 ) cost_drop = cost_drop + beta
            if( labels(ind+1) /= 1 ) cost_keep = cost_keep + beta
        endif
        if( ind > 1 )then
            if( labels(ind-1) == 0 ) cost_keep = cost_keep + 2. * beta
        endif
        if( ind < size(labels) )then
            if( labels(ind+1) == 1 ) cost_drop = cost_drop + 2. * beta
        endif
        if( cost_keep < cost_drop )then
            label_new = 1
        else
            label_new = 0
        endif
    end subroutine update_icm_rank_site

    real function score_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank) result(score)
        real,    intent(in) :: spec(:), beta, alpha
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: i, kfree
        call estimate_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
        kfree = min(max(nmin_rank, 4), nmax_rank)
        score = 0.
        do i = 1, size(spec)
            if( labels(i) == 1 )then
                score = score + (spec(i) - mu_keep)**2 / var_keep
                score = score + alpha * real(max(0, i - kfree)) / real(max(1, nmax_rank - kfree))
            else
                score = score + (spec(i) - mu_drop)**2 / var_drop
            endif
        end do
        do i = 1, size(spec) - 1
            if( labels(i) /= labels(i+1) ) score = score + beta
            if( labels(i) == 0 .and. labels(i+1) == 1 ) score = score + 2. * beta
        end do
    end function score_icm_rank_solution

    subroutine trim_split_embedding(coords, eigvals, nkeep)
        real, allocatable, intent(inout) :: coords(:,:), eigvals(:)
        integer,           intent(in)    :: nkeep
        real, allocatable :: coords_tmp(:,:), eigvals_tmp(:)
        integer :: ndims, nvals
        ndims = max(1, min(nkeep, size(coords,1)))
        if( ndims < size(coords,1) )then
            allocate(coords_tmp(ndims,size(coords,2)), source=coords(1:ndims,:))
            deallocate(coords)
            allocate(coords(ndims,size(coords_tmp,2)), source=coords_tmp)
            deallocate(coords_tmp)
        endif
        if( allocated(eigvals) )then
            nvals = max(1, min(ndims, size(eigvals)))
            if( nvals < size(eigvals) )then
                allocate(eigvals_tmp(nvals), source=eigvals(1:nvals))
                deallocate(eigvals)
                allocate(eigvals(nvals), source=eigvals_tmp)
                deallocate(eigvals_tmp)
            endif
        endif
    end subroutine trim_split_embedding

    subroutine sanitize_embedding_coords(mode, cls_id, coords)
        character(len=*), intent(in)    :: mode
        integer,          intent(in)    :: cls_id
        real,             intent(inout) :: coords(:,:)
        integer :: i, j, nbad
        nbad = 0
        do j = 1, size(coords,2)
            do i = 1, size(coords,1)
                if( ieee_is_finite(coords(i,j)) ) cycle
                coords(i,j) = 0.
                nbad = nbad + 1
            end do
        end do
        if( nbad > 0 )then
            write(logfhandle,'(A,A,A,I8,A,I8)') 'Cls split embedding warning: mode=', trim(mode), &
                ' class=', cls_id, ' nonfinite_coords_reset=', nbad
            call flush(logfhandle)
        endif
    end subroutine sanitize_embedding_coords

    subroutine sanitize_distance_matrix(mode, cls_id, dmat)
        character(len=*), intent(in)    :: mode
        integer,          intent(in)    :: cls_id
        real,             intent(inout) :: dmat(:,:)
        real :: finite_max, smin, smax, delta
        integer :: i, j, nbad
        logical :: have_finite
        nbad = 0
        finite_max = 0.
        have_finite = .false.
        do i = 1, size(dmat,1)
            dmat(i,i) = 0.
            do j = i + 1, size(dmat,2)
                if( .not. ieee_is_finite(dmat(i,j)) )then
                    nbad = nbad + 1
                else
                    finite_max = max(finite_max, dmat(i,j))
                    have_finite = .true.
                endif
            end do
        end do
        if( .not. have_finite ) finite_max = 1.
        do i = 1, size(dmat,1)
            do j = i + 1, size(dmat,2)
                if( ieee_is_finite(dmat(i,j)) ) cycle
                dmat(i,j) = finite_max
                dmat(j,i) = finite_max
            end do
        end do
        if( nbad > 0 )then
            write(logfhandle,'(A,A,A,I8,A,I8,A,F9.4)') 'Cls split distance warning: mode=', trim(mode), &
                ' class=', cls_id, ' nonfinite_pairs_reset=', nbad, ' fill=', finite_max
            call flush(logfhandle)
        endif
        smin  = minval(dmat)
        smax  = maxval(dmat)
        delta = smax - smin
        if( .not. ieee_is_finite(delta) .or. delta <= DTINY )then
            dmat = 0.
            write(logfhandle,'(A,A,A,I8,A,F9.4,A,F9.4)') 'Cls split distance warning: mode=', trim(mode), &
                ' class=', cls_id, ' degenerate distance matrix; min=', smin, ' max=', smax
            call flush(logfhandle)
        else
            dmat = (dmat - smin) / delta
        endif
        do i = 1, size(dmat,1)
            dmat(i,i) = 0.
        end do
    end subroutine sanitize_distance_matrix

    subroutine apply_split_project_updates(spproj, params, nsplit, new_class, new_parent, parent_of_subcls, pop_of_subcls)
        type(sp_project), intent(inout) :: spproj
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nsplit
        integer,          intent(in)    :: new_class(:), new_parent(:)
        integer,          intent(in)    :: parent_of_subcls(:), pop_of_subcls(:)
        integer :: i
        if( trim(params%oritype) .eq. 'ptcl2D' )then
            call spproj%os_ptcl2D%set_all2single('class',   0)
            call spproj%os_ptcl2D%set_all2single('cluster', 0)
            do i = 1, size(new_class)
                if( new_class(i) <= 0 ) cycle
                call spproj%os_ptcl2D%set(i, 'class',   new_class(i))
                call spproj%os_ptcl2D%set(i, 'cluster', new_parent(i))
            end do
        else
            call spproj%os_ptcl3D%set_all2single('class',   0)
            call spproj%os_ptcl3D%set_all2single('cluster', 0)
            do i = 1, size(new_class)
                if( new_class(i) <= 0 ) cycle
                call spproj%os_ptcl3D%set(i, 'class',   new_class(i))
                call spproj%os_ptcl3D%set(i, 'cluster', new_parent(i))
            end do
        endif
        call spproj%os_cls2D%new(nsplit, is_ptcl=.false.)
        call spproj%os_cls3D%new(nsplit, is_ptcl=.false.)
        do i = 1, nsplit
            call spproj%os_cls2D%set(i, 'cluster', parent_of_subcls(i))
            call spproj%os_cls2D%set(i, 'pop',     pop_of_subcls(i))
            call spproj%os_cls2D%set(i, 'accept',  1)
            call spproj%os_cls2D%set(i, 'state',   1)
            call spproj%os_cls3D%set(i, 'cluster', parent_of_subcls(i))
            call spproj%os_cls3D%set(i, 'accept',  1)
            call spproj%os_cls3D%set(i, 'state',   1)
        end do
        if( nsplit > 0 )then
            call write_srchspace_map2D(parent_of_subcls(1:nsplit), string(SRCHSPACE_MAP_FNAME))
        endif
        call spproj%write(params%projfile)
    end subroutine apply_split_project_updates

    subroutine merge_worker_outputs(params, nparts_run)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        type(sp_project) :: spproj
        integer, allocatable :: part_counts(:), part_localstack(:,:), part_parent(:,:), part_local(:,:)
        integer, allocatable :: part_pop(:,:), part_global(:,:)
        integer, allocatable :: comb_part(:), comb_row(:), comb_parent(:), comb_local(:), comb_pop(:)
        integer, allocatable :: new_class(:), new_parent(:)
        type(string) :: map_fname, assign_fname
        integer :: ipart, nlocal, total, max_count, idx, i, funit, ios, pind, parent_cls, local_cls, global_cls
        character(len=XLONGSTRLEN) :: line
        call spproj%read(params%projfile)
        write(logfhandle,'(A,I8)') 'Cls split merge: start nparts=', nparts_run
        call flush(logfhandle)
        allocate(part_counts(nparts_run), source=0)
        total = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_split_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            write(logfhandle,'(A,I8,A,A)') 'Cls split merge: counting map part=', ipart, ' file=', trim(map_fname%to_char())
            call flush(logfhandle)
            call count_data_lines(map_fname, nlocal)
            part_counts(ipart) = nlocal
            total = total + nlocal
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls split merge: counted map part=', ipart, ' rows=', nlocal, &
                ' running_total=', total
            call flush(logfhandle)
            call map_fname%kill
        end do
        if( total < 1 ) THROW_HARD('No subclass outputs produced by distributed cls_split workers')
        max_count = maxval(part_counts)
        write(logfhandle,'(A,I8,A,I8)') 'Cls split merge: maps counted total=', total, ' max_count=', max_count
        call flush(logfhandle)
        allocate(part_localstack(max_count, nparts_run), part_parent(max_count, nparts_run), part_local(max_count, nparts_run), &
                 part_pop(max_count, nparts_run), part_global(max_count, nparts_run), source=0)
        allocate(comb_part(total), comb_row(total), comb_parent(total), comb_local(total), comb_pop(total), source=0)
        idx = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_split_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            write(logfhandle,'(A,I8,A,A)') 'Cls split merge: reading map part=', ipart, ' file=', trim(map_fname%to_char())
            call flush(logfhandle)
            call read_part_map(map_fname, part_counts(ipart), part_localstack(1:part_counts(ipart), ipart), &
                               part_parent(1:part_counts(ipart), ipart), part_local(1:part_counts(ipart), ipart), &
                               part_pop(1:part_counts(ipart), ipart))
            write(logfhandle,'(A,I8,A,I8)') 'Cls split merge: read map part=', ipart, ' rows=', part_counts(ipart)
            call flush(logfhandle)
            do i = 1, part_counts(ipart)
                idx = idx + 1
                comb_part(idx)   = ipart
                comb_row(idx)    = i
                comb_parent(idx) = part_parent(i, ipart)
                comb_local(idx)  = part_local(i, ipart)
                comb_pop(idx)    = part_pop(i, ipart)
            end do
            call map_fname%kill
        end do
        write(logfhandle,'(A,I8)') 'Cls split merge: sorting subclass map rows=', total
        call flush(logfhandle)
        call sort_combined_maps(comb_part, comb_row, comb_parent, comb_local, comb_pop)
        map_fname = string('cls_split_class_map.txt')
        open(newunit=funit, file=map_fname%to_char(), status='replace', action='write')
        write(funit,'(A)') '# global_subclass parent_class local_subclass pop'
        do idx = 1, total
            if( idx > 1 )then
                if( comb_parent(idx) == comb_parent(idx-1) .and. comb_local(idx) == comb_local(idx-1) )then
                    THROW_HARD('Duplicate parent/local subclass pair detected while merging cls_split outputs')
                endif
            endif
            part_global(comb_row(idx), comb_part(idx)) = idx
            write(funit,'(I8,1X,I8,1X,I8,1X,I8)') idx, comb_parent(idx), comb_local(idx), comb_pop(idx)
        end do
        close(funit)
        write(logfhandle,'(A,I8)') 'Cls split merge: subclass map merged total=', total
        call flush(logfhandle)
        if( trim(params%oritype) .eq. 'ptcl2D' )then
            allocate(new_class(spproj%os_ptcl2D%get_noris()), new_parent(spproj%os_ptcl2D%get_noris()), source=0)
        else
            allocate(new_class(spproj%os_ptcl3D%get_noris()), new_parent(spproj%os_ptcl3D%get_noris()), source=0)
        endif
        do ipart = 1, nparts_run
            assign_fname = string('cls_split_assignments_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            open(newunit=funit, file=assign_fname%to_char(), status='old', action='read', iostat=ios)
            call fileiochk('merge_worker_outputs opening '//assign_fname%to_char(), ios)
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*) pind, parent_cls, local_cls
                global_cls = lookup_part_global(part_counts(ipart), part_parent(1:part_counts(ipart), ipart), &
                                                part_local(1:part_counts(ipart), ipart), part_global(1:part_counts(ipart), ipart), &
                                                parent_cls, local_cls)
                if( global_cls <= 0 ) THROW_HARD('Could not map worker-local subclass to global subclass during merge')
                new_class(pind)  = global_cls
                new_parent(pind) = parent_cls
            end do
            close(funit)
            call assign_fname%kill
        end do
        call apply_split_project_updates(spproj, params, total, new_class, new_parent, comb_parent, comb_pop)
        call spproj%kill
        deallocate(new_class, new_parent, part_counts, part_localstack, part_parent, part_local, part_pop, part_global, &
                   comb_part, comb_row, comb_parent, comb_local, comb_pop)
        call map_fname%kill
    end subroutine merge_worker_outputs

    integer function lookup_part_global(nrows, parents, locals, globals, parent_cls, local_cls) result(global_cls)
        integer, intent(in) :: nrows, parents(:), locals(:), globals(:), parent_cls, local_cls
        integer :: i
        global_cls = 0
        do i = 1, nrows
            if( parents(i) == parent_cls .and. locals(i) == local_cls )then
                global_cls = globals(i)
                return
            endif
        end do
    end function lookup_part_global

    subroutine count_data_lines(fname, nlines)
        type(string), intent(in)  :: fname
        integer,      intent(out) :: nlines
        integer :: funit, ios
        character(len=XLONGSTRLEN) :: line
        nlines = 0
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=ios)
        call fileiochk('count_data_lines opening '//fname%to_char(), ios)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            nlines = nlines + 1
        end do
        close(funit)
    end subroutine count_data_lines

    subroutine read_part_map(fname, nrows, localstack, parents, locals, pops)
        type(string), intent(in)  :: fname
        integer,      intent(in)  :: nrows
        integer,      intent(out) :: localstack(:), parents(:), locals(:), pops(:)
        integer :: funit, ios, irow
        character(len=XLONGSTRLEN) :: line
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=ios)
        call fileiochk('read_part_map opening '//fname%to_char(), ios)
        irow = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            irow = irow + 1
            if( irow > nrows ) exit
            read(line,*) localstack(irow), parents(irow), locals(irow), pops(irow)
        end do
        close(funit)
    end subroutine read_part_map

    subroutine read_int_file(fname, vals)
        type(string),         intent(in)  :: fname
        integer, allocatable, intent(out) :: vals(:)
        integer :: nvals, funit, ios, i
        character(len=XLONGSTRLEN) :: line
        call count_data_lines(fname, nvals)
        allocate(vals(nvals))
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=ios)
        call fileiochk('read_int_file opening '//fname%to_char(), ios)
        i = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            i = i + 1
            read(line,*) vals(i)
        end do
        close(funit)
    end subroutine read_int_file

    subroutine sort_order_by_weight_desc(order, weights)
        integer, intent(inout) :: order(:)
        integer, intent(in)    :: weights(:)
        integer :: idx(size(order)), perm(size(order)), sortable(size(order))
        integer :: i, n, tmp
        n = size(order)
        if( n <= 1 ) return
        idx = order
        perm = [(i, i=1, n)]
        do i = 1, n
            sortable(i) = weights(idx(i))
        end do
        call hpsort(sortable, perm)
        order = idx(perm)
        do i = 1, n / 2
            tmp = order(i)
            order(i) = order(n - i + 1)
            order(n - i + 1) = tmp
        end do
    end subroutine sort_order_by_weight_desc

    subroutine sort_combined_maps(parts, rows, parents, locals, pops)
        integer, intent(inout) :: parts(:), rows(:), parents(:), locals(:), pops(:)
        integer :: keys(size(parts)), perm(size(parts))
        integer :: parts_in(size(parts)), rows_in(size(parts)), parents_in(size(parts)), locals_in(size(parts)), pops_in(size(parts))
        integer :: n, max_local, scale, max_parent, i
        n = size(parts)
        if( n <= 1 ) return
        perm = [(i, i=1, n)]
        max_local  = maxval(locals)
        max_parent = maxval(parents)
        scale = max_local + 1
        if( scale <= 0 ) THROW_HARD('sort_combined_maps: invalid local subclass labels')
        if( max_parent > 0 )then
            if( max_parent > huge(scale) / scale )then
                THROW_HARD('sort_combined_maps: key overflow risk for parent/local lexicographic sort')
            endif
        endif
        parts_in   = parts
        rows_in    = rows
        parents_in = parents
        locals_in  = locals
        pops_in    = pops
        keys = (parents_in - 1) * scale + locals_in
        call hpsort(keys, perm)
        parts   = parts_in(perm)
        rows    = rows_in(perm)
        parents = parents_in(perm)
        locals  = locals_in(perm)
        pops    = pops_in(perm)
    end subroutine sort_combined_maps

end module simple_cls_split_strategy
