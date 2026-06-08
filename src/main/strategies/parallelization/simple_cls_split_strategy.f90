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
        if( .not. cline%defined('neigs')    ) call cline%set('neigs',    0)
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'diffusion_maps')
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
    end subroutine init_common

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
        use simple_imgarr_utils,     only: dealloc_imgarr
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: part
        logical,          intent(in)    :: l_write_project
        type(string) :: label, map_fname, assign_fname, raw_fname, den_fname, coeff_fname
        type(string) :: den_ptcl_fname, ptcl_map_fname
        type(image), allocatable :: raw_subavgs(:), den_subavgs(:), coeff_subavgs(:)
        type(image), allocatable :: den_ptcls(:), coeff_ptcls(:)
        integer, allocatable :: cls_inds(:), cls_pops(:), pinds(:), labels(:), new_class(:), new_parent(:), parent_of_subcls(:), pop_of_subcls(:)
        integer :: i, j, k, iglob, iptcl_glob, nsplit, subcls_offset, funit_map, funit_assign, funit_ptcl_map
        logical :: l_phflip, l_pre_norm, l_fixed_nsubcls, l_write_ptcls
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        call determine_split_label(params, build, label)
        l_pre_norm = (trim(params%pre_norm) .eq. 'yes')
        l_fixed_nsubcls = cline%defined('ncls') .and. params%ncls > 1
        select case(trim(params%gen_model))
            case('yes')
                l_write_ptcls = .true.
            case('no')
                l_write_ptcls = .false.
            case DEFAULT
                THROW_HARD('gen_model must be yes or no; cls_split')
        end select
        call determine_phase_flip(spproj, params, l_phflip)
        if( l_write_project )then
            map_fname = string('cls_split_class_map.txt')
            raw_fname = string('cls_split_subclass_avgs.mrcs')
            den_fname = string('cls_split_subclass_avgs_conditioned.mrcs')
            coeff_fname = string('cls_split_subclass_avgs_steer_coeffproj.mrcs')
            if( l_write_ptcls )then
                den_ptcl_fname    = string('cls_split_particles_denoised.mrcs')
                ptcl_map_fname    = string('cls_split_particles_map.txt')
            endif
            call del_file(raw_fname)
            call del_file(den_fname)
            call del_file(coeff_fname)
            if( l_write_ptcls )then
                call del_file(string('cls_split_particles_original.mrcs'))
                call del_file(den_ptcl_fname)
                call del_file(string('cls_split_particles.mrcs'))
                call del_file(string('cls_split_particles_steer_coeffproj.mrcs'))
            endif
            open(newunit=funit_map, file=map_fname%to_char(), status='replace', action='write')
            write(funit_map,'(A)') '# global_subclass parent_class local_subclass pop'
            if( l_write_ptcls )then
                open(newunit=funit_ptcl_map, file=ptcl_map_fname%to_char(), status='replace', action='write')
                write(funit_ptcl_map,'(A)') '# stack_index particle_index parent_class local_subclass global_subclass'
            endif
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
            raw_fname    = string('cls_split_subclass_avgs_part')//int2str_pad(part, params%numlen)//'.mrcs'
            den_fname    = string('cls_split_subclass_avgs_conditioned_part')//int2str_pad(part, params%numlen)//'.mrcs'
            coeff_fname  = string('cls_split_subclass_avgs_steer_coeffproj_part')//int2str_pad(part, params%numlen)//'.mrcs'
            if( l_write_ptcls )then
                den_ptcl_fname    = string('cls_split_particles_denoised_part')//int2str_pad(part, params%numlen)//'.mrcs'
                ptcl_map_fname    = string('cls_split_particles_map_part')//int2str_pad(part, params%numlen)//TXT_EXT
            endif
            call del_file(raw_fname)
            call del_file(den_fname)
            call del_file(coeff_fname)
            if( l_write_ptcls )then
                call del_file(string('cls_split_particles_original_part')//int2str_pad(part, params%numlen)//'.mrcs')
                call del_file(den_ptcl_fname)
                call del_file(string('cls_split_particles_part')//int2str_pad(part, params%numlen)//'.mrcs')
                call del_file(string('cls_split_particles_steer_coeffproj_part')//int2str_pad(part, params%numlen)//'.mrcs')
            endif
            open(newunit=funit_map,    file=map_fname%to_char(),    status='replace', action='write')
            open(newunit=funit_assign, file=assign_fname%to_char(), status='replace', action='write')
            write(funit_map,'(A)')    '# local_stack_index parent_class local_subclass pop'
            write(funit_assign,'(A)') '# particle_index parent_class local_subclass'
            if( l_write_ptcls )then
                open(newunit=funit_ptcl_map, file=ptcl_map_fname%to_char(), status='replace', action='write')
                write(funit_ptcl_map,'(A)') '# stack_index particle_index parent_class local_subclass global_subclass'
            endif
        endif
        iglob = 0
        iptcl_glob = 0
        do i = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8)') 'Cls split starting: class=', cls_inds(i), ' nptcls=', cls_pops(i)
            call flush(logfhandle)
            call split_one_parent_class(params, build, spproj, cls_inds(i), l_phflip, l_pre_norm, l_fixed_nsubcls, &
                                        l_write_ptcls, &
                                        nsplit, pinds, labels, raw_subavgs, den_subavgs, coeff_subavgs, &
                                        den_ptcls, coeff_ptcls)
            if( nsplit < 1 ) cycle
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls split summary: class=', cls_inds(i), ' nptcls=', size(labels), ' nsubcls=', nsplit
            call flush(logfhandle)
            do j = 1, nsplit
                iglob = iglob + 1
                if( l_write_project )then
                    parent_of_subcls(iglob) = cls_inds(i)
                    pop_of_subcls(iglob)    = count(labels == j)
                    write(funit_map,'(I8,1X,I8,1X,I8,1X,I8)') iglob, cls_inds(i), j, count(labels == j)
                    call raw_subavgs(j)%write(raw_fname, iglob)
                    call den_subavgs(j)%write(den_fname, iglob)
                    if( allocated(coeff_subavgs) ) call coeff_subavgs(j)%write(coeff_fname, iglob)
                else
                    write(funit_map,'(I8,1X,I8,1X,I8,1X,I8)') iglob, cls_inds(i), j, count(labels == j)
                    call raw_subavgs(j)%write(raw_fname, iglob)
                    call den_subavgs(j)%write(den_fname, iglob)
                    if( allocated(coeff_subavgs) ) call coeff_subavgs(j)%write(coeff_fname, iglob)
                endif
            end do
            subcls_offset = iglob - nsplit
            if( l_write_ptcls .and. allocated(den_ptcls) )then
                do k = 1, size(labels)
                    if( labels(k) <= 0 ) cycle
                    iptcl_glob = iptcl_glob + 1
                    call den_ptcls(k)%write(den_ptcl_fname, iptcl_glob)
                    write(funit_ptcl_map,'(I10,1X,I10,1X,I10,1X,I10,1X,I10)') &
                        iptcl_glob, pinds(k), cls_inds(i), labels(k), subcls_offset + labels(k)
                end do
            endif
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
            if( allocated(raw_subavgs) ) call dealloc_imgarr(raw_subavgs)
            if( allocated(den_subavgs) ) call dealloc_imgarr(den_subavgs)
            if( allocated(coeff_subavgs) ) call dealloc_imgarr(coeff_subavgs)
            if( allocated(den_ptcls)    ) call dealloc_imgarr(den_ptcls)
            if( allocated(coeff_ptcls)  ) call dealloc_imgarr(coeff_ptcls)
            if( allocated(pinds) ) deallocate(pinds)
            if( allocated(labels) ) deallocate(labels)
        end do
        close(funit_map)
        if( l_write_ptcls ) close(funit_ptcl_map)
        if( .not. l_write_project ) close(funit_assign)
        if( l_write_project .and. l_write_ptcls ) call sort_cls_split_particle_outputs_by_pind(params%smpd_crop)
        if( l_write_project )then
            call apply_split_project_updates(spproj, params, iglob, new_class, new_parent, parent_of_subcls(1:iglob), pop_of_subcls(1:iglob))
            deallocate(new_class, new_parent, parent_of_subcls, pop_of_subcls)
        endif
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
        call label%kill
        call map_fname%kill
        if( assign_fname%is_allocated() ) call assign_fname%kill
        if( raw_fname%is_allocated() ) call raw_fname%kill
        if( den_fname%is_allocated() ) call den_fname%kill
        if( coeff_fname%is_allocated() ) call coeff_fname%kill
        if( den_ptcl_fname%is_allocated()    ) call den_ptcl_fname%kill
        if( ptcl_map_fname%is_allocated()    ) call ptcl_map_fname%kill
    end subroutine run_local_split

    subroutine split_one_parent_class(params, build, spproj, cls_id, l_phflip, l_pre_norm, l_fixed_nsubcls, &
                                      l_write_ptcls, &
                                      nsplit, pinds, labels, raw_subavgs, den_subavgs, coeff_subavgs, &
                                      den_ptcls, coeff_ptcls)
        use simple_diffusion_maps,   only: steerable_transport_denoise, steerable_coeffproj_denoise, graph_coeffproj_denoise
        use simple_imgarr_utils,     only: dealloc_imgarr, copy_imgarr
        use simple_imgproc,          only: make_pcavecs
        use simple_clustering_utils, only: cluster_dmat, silhouette_score
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id
        logical,                  intent(in)    :: l_phflip, l_pre_norm, l_fixed_nsubcls, l_write_ptcls
        integer,                  intent(out)   :: nsplit
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        type(image), allocatable, intent(out)   :: raw_subavgs(:), den_subavgs(:), coeff_subavgs(:)
        type(image), allocatable, intent(out)   :: den_ptcls(:), coeff_ptcls(:)
        type(parameters) :: params_mask
        type(image), allocatable :: imgs(:), imgs_ppca(:), class_mask(:)
        type(image) :: cavg_raw, cavg_den, cavg_coeff
        real, allocatable :: avg(:), pcavecs(:,:), den_pcavecs(:,:), coords(:,:), eigvals(:), dmat(:,:), steer_aff(:,:), steer_theta(:,:)
        real, allocatable :: steer_shift_x(:,:), steer_shift_y(:,:)
        real, allocatable :: class_diams(:), class_shifts(:,:)
        integer, allocatable :: i_medoids(:)
        integer :: nptcls, npix, nsplit_count, j, k, class_ldim(3)
        real    :: class_moldiam, class_mskdiam, class_mskrad, dval, sdev_noise
        logical :: l_so3_split
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
        call automask2D(params_mask, class_mask, params_mask%ngrow, nint(params_mask%winsz), params_mask%edge, class_diams, class_shifts)
        class_moldiam = params_mask%smpd * real(min(round2even(class_diams(1) / params_mask%smpd + 2. * COSMSKHALFWIDTH), class_ldim(1)))
        class_mskdiam = class_moldiam * MSK_EXP_FAC
        class_mskrad  = min(real(class_ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * class_mskdiam / params_mask%smpd)
        write(logfhandle,'(A,I8,A,F8.2,A,F8.2,A,F8.2)') 'Cls split mask update: class=', cls_id, &
            ' automask_diam=', class_diams(1), ' mskdiam=', class_mskdiam, ' mskrad=', class_mskrad
        call flush(logfhandle)
        do j = 1, nptcls
            call imgs_ppca(j)%norm_noise(build%lmsk, sdev_noise)
            call imgs_ppca(j)%mask2D_softavg(class_mskrad)
        end do
        if( l_pre_norm )then
            do j = 1, nptcls
                call imgs_ppca(j)%norm
            end do
        endif
        l_so3_split = trim(params%pca_mode) == 'diff_map_so3'
        if( .not. l_so3_split )then
            call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
        endif
        if( l_so3_split )then
            call make_so3_split_embedding(params, spproj, pinds, cls_id, nptcls, imgs_ppca, coords, eigvals, &
                                          steer_aff, steer_theta, steer_shift_x, steer_shift_y)
        else
            if( l_write_ptcls )then
                call make_split_embedding(params, cls_id, nptcls, npix, pcavecs, imgs_ppca, coords, eigvals, &
                                          steer_aff, steer_theta, den_pcavecs)
            else
                call make_split_embedding(params, cls_id, nptcls, npix, pcavecs, imgs_ppca, coords, eigvals, &
                                          steer_aff, steer_theta)
            endif
        endif
        if( l_write_ptcls .and. trim(params%pca_mode) == 'diffusion_maps' .and. allocated(den_pcavecs) )then
            den_ptcls = copy_imgarr(imgs_ppca)
            do j = 1, nptcls
                call den_ptcls(j)%unserialize(avg + den_pcavecs(:,j))
            end do
            write(logfhandle,'(A,I8,A,I8,A)') 'Cls split diffusion-map decoder denoising: class=', cls_id, &
                ' n=', nptcls, ' method=joint'
            call flush(logfhandle)
        endif
        if( l_write_ptcls .and. trim(params%pca_mode) == 'steerable_diff_map' )then
            call steerable_transport_denoise(params, imgs_ppca, avg, steer_aff, steer_theta, den_ptcls)
            call steerable_coeffproj_denoise(params, imgs_ppca, avg, steer_aff, steer_theta, coeff_ptcls)
            if( .not. allocated(coeff_ptcls) )then
                write(logfhandle,'(A,I8,A)') 'Cls split steerable coefficient denoising fallback: class=', cls_id, &
                    ' using transported/normal particles'
                call flush(logfhandle)
                if( allocated(den_ptcls) )then
                    coeff_ptcls = copy_imgarr(den_ptcls)
                else
                    coeff_ptcls = copy_imgarr(imgs_ppca)
                endif
            endif
        endif
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
            nsplit_count = particle_count_nsplit(nptcls, params%nptcls_per_subcls, params%nsubcls_min, params%nsubcls_max)
            call select_kmedoids_by_silhouette(dmat, cls_id, params%nsubcls_min, params%nsubcls_max, &
                                               nsplit, i_medoids, labels)
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,I8)') 'Cls split auto ncls: class=', cls_id, &
                ' size=', nptcls, ' count_suggest=', nsplit_count, ' min=', params%nsubcls_min, &
                ' max=', params%nsubcls_max, ' selected=', nsplit
            call flush(logfhandle)
        endif
        if( l_write_ptcls .and. trim(params%pca_mode) == 'diff_map_so3' )then
            call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
            call graph_coeffproj_denoise(params, imgs_ppca, avg, steer_aff, steer_theta, steer_shift_x, steer_shift_y, &
                                         trim(params%so3_steering), coeff_ptcls)
            if( .not. allocated(coeff_ptcls) )then
                write(logfhandle,'(A,I8,A,A,A)') 'Cls split graph coefficient denoising fallback: class=', cls_id, &
                    ' steering=', trim(params%so3_steering), ' no coeff-projected averages written'
                call flush(logfhandle)
            endif
        endif
        call cavg_den%new(cavg_raw%get_ldim(), cavg_raw%get_smpd())
        if( allocated(coeff_ptcls) ) call cavg_coeff%new(cavg_raw%get_ldim(), cavg_raw%get_smpd())
        allocate(raw_subavgs(nsplit), den_subavgs(nsplit))
        if( allocated(coeff_ptcls) ) allocate(coeff_subavgs(nsplit))
        do j = 1, nsplit
            call cavg_raw%zero_and_unflag_ft
            call cavg_den%zero_and_unflag_ft
            if( allocated(coeff_ptcls) ) call cavg_coeff%zero_and_unflag_ft
            do k = 1, size(labels)
                if( labels(k) /= j ) cycle
                call cavg_raw%add(imgs(k))
                if( allocated(den_ptcls) )then
                    call cavg_den%add(den_ptcls(k))
                else
                    call cavg_den%add(imgs_ppca(k))
                endif
                if( allocated(coeff_ptcls) ) call cavg_coeff%add(coeff_ptcls(k))
            end do
            call cavg_raw%div(real(max(count(labels == j), 1)))
            call cavg_den%div(real(max(count(labels == j), 1)))
            if( allocated(coeff_ptcls) ) call cavg_coeff%div(real(max(count(labels == j), 1)))
            call raw_subavgs(j)%copy(cavg_raw)
            call den_subavgs(j)%copy(cavg_den)
            if( allocated(coeff_ptcls) ) call coeff_subavgs(j)%copy(cavg_coeff)
        end do
        call cavg_raw%kill
        call cavg_den%kill
        if( cavg_coeff%exists() ) call cavg_coeff%kill
        if( allocated(coords)       ) deallocate(coords)
        if( allocated(eigvals)      ) deallocate(eigvals)
        if( allocated(dmat)         ) deallocate(dmat)
        if( allocated(steer_aff)    ) deallocate(steer_aff)
        if( allocated(steer_theta)  ) deallocate(steer_theta)
        if( allocated(steer_shift_x)) deallocate(steer_shift_x)
        if( allocated(steer_shift_y)) deallocate(steer_shift_y)
        if( allocated(i_medoids)    ) deallocate(i_medoids)
        if( allocated(avg)          ) deallocate(avg)
        if( allocated(pcavecs)      ) deallocate(pcavecs)
        if( allocated(den_pcavecs)  ) deallocate(den_pcavecs)
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
                score = silhouette_score(trial_labels, dmat) / real(k_trial)
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

    subroutine make_split_embedding(params, cls_id, nptcls, npix, pcavecs, imgs_ppca, coords, eigvals, steer_aff, steer_theta, den_pcavecs)
        use simple_diffusion_maps,           only: diffusion_map_embedder, steerable_diffusion_map_embedder
        use simple_kpca_svd,                 only: kpca_svd
        use simple_ppca,                     only: ppca
        type(parameters),  intent(inout) :: params
        integer,           intent(in)    :: cls_id, nptcls, npix
        real,              intent(in)    :: pcavecs(npix,nptcls)
        type(image),        intent(in)    :: imgs_ppca(:)
        real, allocatable, intent(out)   :: coords(:,:)
        real, allocatable, intent(out)   :: eigvals(:)
        real, allocatable, optional, intent(out) :: steer_aff(:,:), steer_theta(:,:)
        real, allocatable, optional, intent(out) :: den_pcavecs(:,:)
        type(diffusion_map_embedder) :: diffmap
        type(steerable_diffusion_map_embedder) :: steerable_diffmap
        type(kpca_svd) :: kpca_model
        type(ppca)     :: ppca_model
        real, allocatable :: feat(:), tmpvec(:)
        integer :: neigs, neigs_scan, neigs_used, i
        logical :: l_auto_neigs, l_decode
        l_auto_neigs = params%neigs <= 0
        l_decode = present(den_pcavecs) .and. trim(params%pca_mode) == 'diffusion_maps'
        if( l_auto_neigs )then
            neigs_scan = cls_split_auto_neigs_scan(trim(params%pca_mode), nptcls)
        else
            neigs_scan = params%neigs
        endif
        if( trim(params%pca_mode) .eq. 'diffusion_maps' )then
            neigs = min(max(neigs_scan, 0), max(nptcls-2, 1))
        else
            neigs = min(max(neigs_scan, 1), max(nptcls-1, 1))
        endif
        select case(trim(params%pca_mode))
            case('diffusion_maps')
                call diffmap%set_params(neigs, min(max(2, params%k_nn), max(2, nptcls-1)))
                if( l_decode )then
                    call diffmap%set_preimage_params(.true., decoder_ridge=1.e-2, &
                        joint_decoder=.true., &
                        joint_decoder_iters=80, joint_decoder_weight=0.5, joint_decoder_tol=1.e-4)
                    call diffmap%embed(pcavecs, coords, eigvals, pcavecs)
                else
                    call diffmap%embed(pcavecs, coords, eigvals)
                endif
            case('steerable_diff_map')
                call steerable_diffmap%set_params(neigs, min(max(2, params%k_nn), max(2, nptcls-1)), params%steerable_nmodes)
                if( present(steer_aff) .and. present(steer_theta) )then
                    call steerable_diffmap%embed(params, imgs_ppca, coords, eigvals, steer_aff, steer_theta)
                else
                    call steerable_diffmap%embed(params, imgs_ppca, coords, eigvals)
                endif
            case('ppca')
                call ppca_model%new(nptcls, npix, neigs)
                call ppca_model%set_verbose(.false.)
                call ppca_model%master(pcavecs, 15)
                allocate(coords(neigs,nptcls), source=0.)
                do i = 1, nptcls
                    feat = ppca_model%get_feat(i)
                    coords(:,i) = feat
                    if( allocated(feat) ) deallocate(feat)
                end do
                eigvals = ppca_model%get_signal_eigvals()
                call ppca_model%kill
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
        end select
        if( l_auto_neigs .and. allocated(eigvals) )then
            neigs_used = select_neigs_icm(eigvals, trim(params%pca_mode), size(coords,1), cls_id)
            if( l_decode ) call diffmap%trim_preimage_model(neigs_used)
            call trim_split_embedding(coords, eigvals, neigs_used)
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8)') 'Cls split auto neigs: mode=', trim(params%pca_mode), &
                ' class=', cls_id, ' size=', nptcls, ' scan=', neigs, ' selected=', size(coords,1)
            call flush(logfhandle)
        endif
        if( l_decode )then
            allocate(den_pcavecs(npix,nptcls), source=0.)
            do i = 1, nptcls
                call diffmap%joint_decode(coords(:,i), tmpvec)
                den_pcavecs(:,i) = tmpvec
                if( allocated(tmpvec) ) deallocate(tmpvec)
            end do
            call diffmap%kill_preimage_model()
        endif
        write(logfhandle,'(A,A,A,I8,A,I8,A,I8)') 'Cls split embedding: mode=', trim(params%pca_mode), &
            ' class=', cls_id, ' size=', nptcls, ' dims=', size(coords,1)
        call flush(logfhandle)
    end subroutine make_split_embedding

    subroutine make_so3_split_embedding(params, spproj, pinds, cls_id, nptcls, imgs_ppca, coords, eigvals, &
                                        graph_aff, graph_theta, graph_shift_x, graph_shift_y)
        use simple_diff_map_graphs, only: build_so3_split_affinity
        use simple_diffusion_maps,   only: embed_affinity, embed_so2_steerable_affinity, embed_se2_steerable_affinity
        type(parameters),  intent(inout) :: params
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: pinds(:), cls_id, nptcls
        type(image),       intent(inout) :: imgs_ppca(:)
        real, allocatable, intent(out)   :: coords(:,:), eigvals(:)
        real, allocatable, optional, intent(out) :: graph_aff(:,:), graph_theta(:,:), graph_shift_x(:,:), graph_shift_y(:,:)
        real, allocatable :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        real :: se2_shift_scale
        character(len=STDLEN) :: so3_graph, so3_steering
        integer :: neigs, neigs_scan, neigs_used
        logical :: l_auto_neigs

        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('pca_mode=diff_map_so3 requires oritype=ptcl3D')
        so3_graph = lowercase(trim(params%so3_graph))
        select case(trim(so3_graph))
            case('cluster2d', 'projection_registration')
            case DEFAULT
                THROW_HARD('so3_graph must be cluster2d or projection_registration')
        end select
        so3_steering = lowercase(trim(params%so3_steering))
        select case(trim(so3_steering))
            case('none', 'so2', 'se2')
            case DEFAULT
                THROW_HARD('so3_steering must be none, so2, or se2')
        end select
        l_auto_neigs = params%neigs <= 0
        if( l_auto_neigs )then
            neigs_scan = cls_split_auto_neigs_scan(trim(params%pca_mode), nptcls)
        else
            neigs_scan = params%neigs
        endif
        neigs = min(max(neigs_scan, 1), max(nptcls - 2, 1))
        call build_so3_split_affinity(params, spproj, pinds, imgs_ppca, aff, theta, shift_x, shift_y)
        select case(trim(so3_steering))
            case('none')
                call embed_affinity(aff, neigs, coords, eigvals)
            case('se2')
                se2_shift_scale = estimate_se2_shift_scale(params, aff, shift_x, shift_y)
                write(logfhandle,'(A,I8,A,F10.4)') 'SO3 SE2 steering: class=', cls_id, &
                    ' shift_scale=', se2_shift_scale
                call flush(logfhandle)
                call embed_se2_steerable_affinity(aff, theta, shift_x, shift_y, neigs, params%steerable_nmodes, &
                    se2_shift_scale, coords, eigvals)
            case('so2')
                call embed_so2_steerable_affinity(aff, theta, neigs, params%steerable_nmodes, coords, eigvals)
        end select
        if( l_auto_neigs .and. allocated(eigvals) )then
            neigs_used = select_neigs_icm(eigvals, trim(params%pca_mode), size(coords,1), cls_id)
            call trim_split_embedding(coords, eigvals, neigs_used)
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8)') 'Cls split auto neigs: mode=', trim(params%pca_mode), &
                ' class=', cls_id, ' size=', nptcls, ' scan=', neigs, ' selected=', size(coords,1)
            call flush(logfhandle)
        endif
        write(logfhandle,'(A,A,A,A,A,A,A,I8,A,I8,A,I8)') 'Cls split embedding: mode=', trim(params%pca_mode), &
            ' so3_graph=', trim(so3_graph), ' so3_steering=', trim(so3_steering), &
            ' class=', cls_id, ' size=', nptcls, ' dims=', size(coords,1)
        call flush(logfhandle)
        if( present(graph_aff)     ) allocate(graph_aff(size(aff,1),size(aff,2)), source=aff)
        if( present(graph_theta)   ) allocate(graph_theta(size(theta,1),size(theta,2)), source=theta)
        if( present(graph_shift_x) ) allocate(graph_shift_x(size(shift_x,1),size(shift_x,2)), source=shift_x)
        if( present(graph_shift_y) ) allocate(graph_shift_y(size(shift_y,1),size(shift_y,2)), source=shift_y)
        if( allocated(aff)   ) deallocate(aff)
        if( allocated(theta) ) deallocate(theta)
        if( allocated(shift_x) ) deallocate(shift_x)
        if( allocated(shift_y) ) deallocate(shift_y)
    end subroutine make_so3_split_embedding

    real function estimate_se2_shift_scale(params, aff, shift_x, shift_y) result(shift_scale)
        type(parameters), intent(in) :: params
        real,             intent(in) :: aff(:,:), shift_x(:,:), shift_y(:,:)
        real :: observed
        if( any(aff > DTINY) )then
            observed = max(maxval(abs(shift_x), mask=(aff > DTINY)), maxval(abs(shift_y), mask=(aff > DTINY)))
        else
            observed = 0.
        endif
        shift_scale = max(1., observed)
        if( params%trs > 0. ) shift_scale = max(shift_scale, params%trs)
    end function estimate_se2_shift_scale

    integer function cls_split_auto_neigs_scan(mode, nptcls) result(neigs_scan)
        character(len=*), intent(in) :: mode
        integer,          intent(in) :: nptcls
        select case(trim(mode))
            case('diffusion_maps', 'diff_map_so3')
                neigs_scan = min(24, max(1, nptcls - 2))
                if( nptcls > 3 ) neigs_scan = max(2, neigs_scan)
            case DEFAULT
                neigs_scan = min(24, max(1, nptcls - 1))
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
        if( (trim(mode) .eq. 'diffusion_maps' .or. trim(mode) .eq. 'diff_map_so3') .and. max_neigs >= 2 ) nmin_rank = 2
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

    integer function particle_count_nsplit(nptcls, target_pop, nmin, nmax) result(nsplit)
        integer, intent(in) :: nptcls, target_pop, nmin, nmax
        integer :: nlo, nhi, target
        nlo    = max(1, nmin)
        nhi    = max(nlo, nmax)
        target = max(1, target_pop)
        nsplit = ceiling(real(nptcls) / real(target))
        nsplit = max(nlo, min(nhi, nsplit))
    end function particle_count_nsplit

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
        use simple_stack_io, only: stack_io
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        type(sp_project) :: spproj
        type(image)      :: img
        type(stack_io), allocatable :: raw_ios(:), den_ios(:), coeff_ios(:)
        type(stack_io) :: raw_out_io, den_out_io, coeff_out_io
        integer, allocatable :: part_counts(:), part_localstack(:,:), part_parent(:,:), part_local(:,:)
        integer, allocatable :: part_pop(:,:), part_global(:,:)
        integer, allocatable :: comb_part(:), comb_row(:), comb_parent(:), comb_local(:), comb_pop(:), last_localstack(:)
        integer, allocatable :: new_class(:), new_parent(:)
        type(string) :: map_fname, assign_fname, raw_fname, den_fname, coeff_fname
        integer, parameter :: MERGE_STACK_BUFSZ = 64
        integer :: ipart, nlocal, total, max_count, idx, i, funit, ios, pind, parent_cls, local_cls, global_cls, ldim(3), nstk
        real    :: smpd
        logical :: have_coeff_avgs
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
        write(logfhandle,'(A)') 'Cls split merge: deleting previous global class-average stacks'
        call flush(logfhandle)
        call del_file(string('cls_split_subclass_avgs.mrcs'))
        call del_file(string('cls_split_subclass_avgs_conditioned.mrcs'))
        call del_file(string('cls_split_subclass_avgs_steer_coeffproj.mrcs'))
        raw_fname = string('cls_split_subclass_avgs_part')//int2str_pad(comb_part(1), params%numlen)//'.mrcs'
        coeff_fname = string('cls_split_subclass_avgs_steer_coeffproj_part')//int2str_pad(comb_part(1), params%numlen)//'.mrcs'
        have_coeff_avgs = trim(params%pca_mode) == 'steerable_diff_map' .and. file_exists(coeff_fname%to_char())
        call coeff_fname%kill
        if( have_coeff_avgs )then
            do ipart = 1, nparts_run
                coeff_fname = string('cls_split_subclass_avgs_steer_coeffproj_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
                if( .not. file_exists(coeff_fname%to_char()) )then
                    write(logfhandle,'(A,I8,A,A)') 'Cls split merge: disabling coeff-stack merge; missing part=', ipart, &
                        ' file=', trim(coeff_fname%to_char())
                    call flush(logfhandle)
                    have_coeff_avgs = .false.
                    exit
                endif
                call coeff_fname%kill
            end do
        endif
        if( .not. file_exists(raw_fname%to_char()) )then
            write(logfhandle,'(A,A)') 'Cls split merge error: missing first raw class-average stack file=', &
                trim(raw_fname%to_char())
            call flush(logfhandle)
            THROW_HARD('Missing raw class-average stack while merging cls_split outputs')
        endif
        write(logfhandle,'(A,A)') 'Cls split merge: probing class-average stack file=', trim(raw_fname%to_char())
        call flush(logfhandle)
        call find_ldim_nptcls(raw_fname, ldim, nstk)
        write(logfhandle,'(A,3(I8,1X),A,I8,A,L1)') 'Cls split merge: class-average stack dims=', ldim, &
            ' nstk=', nstk, ' have_coeff=', have_coeff_avgs
        call flush(logfhandle)
        smpd = params%smpd_crop
        ldim(3) = 1
        call img%new(ldim, smpd)
        allocate(raw_ios(nparts_run), den_ios(nparts_run))
        if( have_coeff_avgs ) allocate(coeff_ios(nparts_run))
        do ipart = 1, nparts_run
            if( part_counts(ipart) < 1 ) cycle
            raw_fname = string('cls_split_subclass_avgs_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            den_fname = string('cls_split_subclass_avgs_conditioned_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            if( .not. file_exists(raw_fname%to_char()) )then
                write(logfhandle,'(A,A)') 'Cls split merge error: missing raw class-average stack file=', trim(raw_fname%to_char())
                call flush(logfhandle)
                THROW_HARD('Missing raw class-average stack while merging cls_split outputs')
            endif
            if( .not. file_exists(den_fname%to_char()) )then
                write(logfhandle,'(A,A)') 'Cls split merge error: missing conditioned class-average stack file=', &
                    trim(den_fname%to_char())
                call flush(logfhandle)
                THROW_HARD('Missing conditioned class-average stack while merging cls_split outputs')
            endif
            write(logfhandle,'(A,I8,A,I8)') 'Cls split merge: opening class-average part=', ipart, ' rows=', part_counts(ipart)
            call flush(logfhandle)
            call raw_ios(ipart)%open(raw_fname, smpd, 'READ', bufsz=MERGE_STACK_BUFSZ)
            call den_ios(ipart)%open(den_fname, smpd, 'READ', bufsz=MERGE_STACK_BUFSZ)
            if( have_coeff_avgs )then
                coeff_fname = string('cls_split_subclass_avgs_steer_coeffproj_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
                call coeff_ios(ipart)%open(coeff_fname, smpd, 'READ', bufsz=MERGE_STACK_BUFSZ)
                call coeff_fname%kill
            endif
            call raw_fname%kill
            call den_fname%kill
        end do
        call raw_out_io%open(string('cls_split_subclass_avgs.mrcs'), smpd, 'WRITE', box=ldim(1), bufsz=MERGE_STACK_BUFSZ)
        call den_out_io%open(string('cls_split_subclass_avgs_conditioned.mrcs'), smpd, 'WRITE', box=ldim(1), &
            bufsz=MERGE_STACK_BUFSZ)
        if( have_coeff_avgs ) call coeff_out_io%open(string('cls_split_subclass_avgs_steer_coeffproj.mrcs'), &
            smpd, 'WRITE', box=ldim(1), bufsz=MERGE_STACK_BUFSZ)
        allocate(last_localstack(nparts_run), source=0)
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
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8)') 'Cls split merge avg: idx=', idx, ' total=', total, &
                ' part=', comb_part(idx), ' row=', comb_row(idx), ' stack=', part_localstack(comb_row(idx), comb_part(idx))
            call flush(logfhandle)
            if( part_localstack(comb_row(idx), comb_part(idx)) <= last_localstack(comb_part(idx)) )then
                THROW_HARD('Non-monotonic class-average stack merge order; merge_worker_outputs')
            endif
            last_localstack(comb_part(idx)) = part_localstack(comb_row(idx), comb_part(idx))
            call raw_ios(comb_part(idx))%read(part_localstack(comb_row(idx), comb_part(idx)), img)
            call raw_out_io%write(idx, img)
            call den_ios(comb_part(idx))%read(part_localstack(comb_row(idx), comb_part(idx)), img)
            call den_out_io%write(idx, img)
            if( have_coeff_avgs )then
                call coeff_ios(comb_part(idx))%read(part_localstack(comb_row(idx), comb_part(idx)), img)
                call coeff_out_io%write(idx, img)
            endif
            write(logfhandle,'(A,I8)') 'Cls split merge avg: done idx=', idx
            call flush(logfhandle)
        end do
        close(funit)
        call raw_out_io%close
        call den_out_io%close
        if( have_coeff_avgs ) call coeff_out_io%close
        do ipart = 1, nparts_run
            if( part_counts(ipart) < 1 ) cycle
            call raw_ios(ipart)%close
            call den_ios(ipart)%close
            if( have_coeff_avgs ) call coeff_ios(ipart)%close
        end do
        deallocate(raw_ios, den_ios, last_localstack)
        if( allocated(coeff_ios) ) deallocate(coeff_ios)
        write(logfhandle,'(A,I8)') 'Cls split merge: class-average stacks merged total=', total
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
        if( trim(params%gen_model) == 'yes' )then
            call merge_worker_particle_stacks(params, nparts_run, part_counts, part_parent, part_local, part_global)
        endif
        call spproj%kill
        call img%kill
        deallocate(new_class, new_parent, part_counts, part_localstack, part_parent, part_local, part_pop, part_global, &
                   comb_part, comb_row, comb_parent, comb_local, comb_pop)
        call map_fname%kill
    end subroutine merge_worker_outputs

    subroutine merge_worker_particle_stacks(params, nparts_run, part_counts, part_parent, part_local, part_global)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        integer,          intent(in)    :: part_counts(:), part_parent(:,:), part_local(:,:), part_global(:,:)
        type(image) :: img
        type(string) :: part_map_fname, den_fname, out_map_fname
        integer :: ipart, funit_in, funit_out, ios, stack_idx, pind, parent_cls, local_cls, local_global
        integer :: global_cls, global_stack, first_part, ldim(3), nstk
        real    :: smpd
        logical :: have_ptcls
        character(len=XLONGSTRLEN) :: line
        have_ptcls = .false.
        first_part = 0
        do ipart = 1, nparts_run
            den_fname = string('cls_split_particles_denoised_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            if( file_exists(den_fname%to_char()) )then
                have_ptcls = .true.
                first_part = ipart
                exit
            endif
            call den_fname%kill
        end do
        if( .not. have_ptcls ) return
        if( den_fname%is_allocated() ) call den_fname%kill
        call del_file(string('cls_split_particles_original.mrcs'))
        call del_file(string('cls_split_particles.mrcs'))
        call del_file(string('cls_split_particles_denoised.mrcs'))
        call del_file(string('cls_split_particles_steer_coeffproj.mrcs'))
        den_fname = string('cls_split_particles_denoised_part')//int2str_pad(first_part, params%numlen)//'.mrcs'
        call find_ldim_nptcls(den_fname, ldim, nstk)
        smpd = params%smpd_crop
        ldim(3) = 1
        call img%new(ldim, smpd)
        out_map_fname = string('cls_split_particles_map.txt')
        open(newunit=funit_out, file=out_map_fname%to_char(), status='replace', action='write', iostat=ios)
        if( ios /= 0 )then
            write(logfhandle,'(A,A)') 'Cls split particle-stack merge skipped: could not open ', out_map_fname%to_char()
            call flush(logfhandle)
            call img%kill
            call out_map_fname%kill
            if( den_fname%is_allocated() ) call den_fname%kill
            return
        endif
        write(funit_out,'(A)') '# stack_index particle_index parent_class local_subclass global_subclass'
        global_stack = 0
        do ipart = 1, nparts_run
            part_map_fname = string('cls_split_particles_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            den_fname      = string('cls_split_particles_denoised_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            if( .not. file_exists(part_map_fname%to_char()) .or. .not. file_exists(den_fname%to_char()) )then
                call part_map_fname%kill
                call den_fname%kill
                cycle
            endif
            open(newunit=funit_in, file=part_map_fname%to_char(), status='old', action='read', iostat=ios)
            if( ios /= 0 )then
                write(logfhandle,'(A,A)') 'Cls split particle-stack merge skipping unreadable map ', part_map_fname%to_char()
                call flush(logfhandle)
                call part_map_fname%kill
                call den_fname%kill
                cycle
            endif
            do
                read(funit_in,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*,iostat=ios) stack_idx, pind, parent_cls, local_cls, local_global
                if( ios /= 0 )then
                    write(logfhandle,'(A,A)') 'Cls split particle-stack merge skipping malformed row in ', part_map_fname%to_char()
                    call flush(logfhandle)
                    cycle
                endif
                global_cls = lookup_part_global(part_counts(ipart), part_parent(1:part_counts(ipart), ipart), &
                                                part_local(1:part_counts(ipart), ipart), part_global(1:part_counts(ipart), ipart), &
                                                parent_cls, local_cls)
                if( global_cls <= 0 )then
                    write(logfhandle,'(A,I8,A,I8)') 'Cls split particle-stack merge skipping unmapped row: parent=', parent_cls, &
                        ' local=', local_cls
                    call flush(logfhandle)
                    cycle
                endif
                global_stack = global_stack + 1
                call img%read(den_fname, stack_idx)
                call img%write(string('cls_split_particles_denoised.mrcs'), global_stack)
                write(funit_out,'(I10,1X,I10,1X,I10,1X,I10,1X,I10)') &
                    global_stack, pind, parent_cls, local_cls, global_cls
            end do
            close(funit_in)
            call part_map_fname%kill
            call den_fname%kill
        end do
        close(funit_out)
        call sort_cls_split_particle_outputs_by_pind(params%smpd_crop)
        write(logfhandle,'(A,I10)') 'Cls split merged particle stacks: n=', global_stack
        call flush(logfhandle)
        call img%kill
        call out_map_fname%kill
        if( den_fname%is_allocated() ) call den_fname%kill
    end subroutine merge_worker_particle_stacks

    subroutine sort_cls_split_particle_outputs_by_pind(smpd)
        real, intent(in) :: smpd
        type(image) :: img
        type(string) :: map_fname, den_fname, tmp_map_fname, tmp_den_fname
        integer, allocatable :: stack_idxs(:), pinds(:), parent_cls(:), local_cls(:), global_cls(:)
        integer, allocatable :: keys(:), perm(:)
        integer :: nrows, funit_in, funit_out, ios, irow, iout, src_row, ldim(3), nstk
        character(len=XLONGSTRLEN) :: line
        map_fname    = string('cls_split_particles_map.txt')
        den_fname    = string('cls_split_particles_denoised.mrcs')
        if( .not. file_exists(map_fname%to_char()) .or. .not. file_exists(den_fname%to_char()) )then
            call map_fname%kill
            call den_fname%kill
            return
        endif
        call count_data_lines(map_fname, nrows)
        if( nrows < 2 )then
            call map_fname%kill
            call den_fname%kill
            return
        endif
        allocate(stack_idxs(nrows), pinds(nrows), parent_cls(nrows), local_cls(nrows), global_cls(nrows), &
                 keys(nrows), perm(nrows), source=0)
        open(newunit=funit_in, file=map_fname%to_char(), status='old', action='read', iostat=ios)
        call fileiochk('sort_cls_split_particle_outputs_by_pind opening '//map_fname%to_char(), ios)
        irow = 0
        do
            read(funit_in,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            irow = irow + 1
            if( irow > nrows ) exit
            read(line,*) stack_idxs(irow), pinds(irow), parent_cls(irow), local_cls(irow), global_cls(irow)
        end do
        close(funit_in)
        if( irow /= nrows ) THROW_HARD('Unexpected row count while sorting cls_split particle map')
        keys = pinds
        perm = [(irow, irow=1,nrows)]
        call hpsort(keys, perm)
        do irow = 2, nrows
            if( keys(irow) == keys(irow-1) ) THROW_HARD('Duplicate particle index while sorting cls_split particle stack')
        end do
        if( all(perm == [(irow, irow=1,nrows)]) )then
            call map_fname%kill
            call den_fname%kill
            deallocate(stack_idxs, pinds, parent_cls, local_cls, global_cls, keys, perm)
            return
        endif
        tmp_map_fname    = string('cls_split_particles_map_sorted_tmp.txt')
        tmp_den_fname    = string('cls_split_particles_denoised_sorted_tmp.mrcs')
        call del_file(tmp_map_fname)
        call del_file(tmp_den_fname)
        call find_ldim_nptcls(den_fname, ldim, nstk)
        if( nstk < maxval(stack_idxs) ) THROW_HARD('Particle stack shorter than cls_split particle map while sorting')
        ldim(3) = 1
        call img%new(ldim, smpd)
        do iout = 1, nrows
            src_row = perm(iout)
            call img%read(den_fname, stack_idxs(src_row))
            call img%write(tmp_den_fname, iout)
        end do
        if( img%exists() ) call img%kill
        open(newunit=funit_out, file=tmp_map_fname%to_char(), status='replace', action='write', iostat=ios)
        call fileiochk('sort_cls_split_particle_outputs_by_pind opening '//tmp_map_fname%to_char(), ios)
        write(funit_out,'(A)') '# stack_index particle_index parent_class local_subclass global_subclass'
        do iout = 1, nrows
            src_row = perm(iout)
            write(funit_out,'(I10,1X,I10,1X,I10,1X,I10,1X,I10)') &
                iout, pinds(src_row), parent_cls(src_row), local_cls(src_row), global_cls(src_row)
        end do
        close(funit_out)
        call simple_rename(tmp_den_fname, den_fname)
        call simple_rename(tmp_map_fname, map_fname)
        write(logfhandle,'(A,I10)') 'Cls split denoised particle stack sorted by particle index: n=', nrows
        call flush(logfhandle)
        call map_fname%kill
        call den_fname%kill
        call tmp_map_fname%kill
        call tmp_den_fname%kill
        deallocate(stack_idxs, pinds, parent_cls, local_cls, global_cls, keys, perm)
    end subroutine sort_cls_split_particle_outputs_by_pind

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
