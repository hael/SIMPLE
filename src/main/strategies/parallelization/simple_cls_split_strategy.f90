module simple_cls_split_strategy
use simple_core_module_api
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_cmdline,            only: cmdline
use simple_defs_fname,         only: SRCHSPACE_MAP_FNAME
use simple_qsys_env,           only: qsys_env
use simple_sp_project,         only: sp_project
use simple_image,              only: image
use simple_image_msk,          only: automask2D
use simple_classaverager,      only: transform_ptcls
use simple_srchspace_map2D_io, only: write_srchspace_map2D
use simple_srch_sort_loc,      only: hpsort
implicit none

public :: cls_split_strategy
public :: cls_split_shmem_strategy
public :: cls_split_worker_strategy
public :: cls_split_master_strategy
public :: create_cls_split_strategy

private
#include "simple_local_flags.inc"

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
        class(cls_split_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_local_split(params, build, cline, spproj, part=params%part, l_write_project=.false.)
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
        type(string) :: label, map_fname, assign_fname, raw_fname, den_fname
        type(image), allocatable :: raw_subavgs(:), den_subavgs(:)
        integer, allocatable :: cls_inds(:), cls_pops(:), pinds(:), labels(:), new_class(:), new_parent(:), parent_of_subcls(:), pop_of_subcls(:)
        integer :: i, j, k, iglob, nsplit, funit_map, funit_assign
        logical :: l_phflip, l_pre_norm, l_fixed_nsubcls
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        call determine_split_label(params, build, label)
        l_pre_norm = (trim(params%pre_norm) .eq. 'yes')
        l_fixed_nsubcls = cline%defined('ncls') .and. params%ncls > 1
        call determine_phase_flip(spproj, params, l_phflip)
        if( l_write_project )then
            map_fname = string('cls_split_class_map.txt')
            raw_fname = string('cls_split_subclass_avgs.mrcs')
            den_fname = string('cls_split_subclass_avgs_conditioned.mrcs')
            call del_file(raw_fname)
            call del_file(den_fname)
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
            raw_fname    = string('cls_split_subclass_avgs_part')//int2str_pad(part, params%numlen)//'.mrcs'
            den_fname    = string('cls_split_subclass_avgs_conditioned_part')//int2str_pad(part, params%numlen)//'.mrcs'
            call del_file(raw_fname)
            call del_file(den_fname)
            open(newunit=funit_map,    file=map_fname%to_char(),    status='replace', action='write')
            open(newunit=funit_assign, file=assign_fname%to_char(), status='replace', action='write')
            write(funit_map,'(A)')    '# local_stack_index parent_class local_subclass pop'
            write(funit_assign,'(A)') '# particle_index parent_class local_subclass'
        endif
        iglob = 0
        do i = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8)') 'Cls split starting: class=', cls_inds(i), ' nptcls=', cls_pops(i)
            call flush(logfhandle)
            call split_one_parent_class(params, build, spproj, cls_inds(i), l_phflip, l_pre_norm, l_fixed_nsubcls, &
                                        nsplit, pinds, labels, raw_subavgs, den_subavgs)
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
                else
                    write(funit_map,'(I8,1X,I8,1X,I8,1X,I8)') iglob, cls_inds(i), j, count(labels == j)
                    call raw_subavgs(j)%write(raw_fname, iglob)
                    call den_subavgs(j)%write(den_fname, iglob)
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
            if( allocated(raw_subavgs) ) call dealloc_imgarr(raw_subavgs)
            if( allocated(den_subavgs) ) call dealloc_imgarr(den_subavgs)
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
        if( raw_fname%is_allocated() ) call raw_fname%kill
        if( den_fname%is_allocated() ) call den_fname%kill
    end subroutine run_local_split

    subroutine split_one_parent_class(params, build, spproj, cls_id, l_phflip, l_pre_norm, l_fixed_nsubcls, nsplit, pinds, labels, raw_subavgs, den_subavgs)
        use simple_imgarr_utils,     only: dealloc_imgarr, copy_imgarr
        use simple_imgproc,          only: make_pcavecs
        use simple_clustering_utils, only: cluster_dmat
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id
        logical,                  intent(in)    :: l_phflip, l_pre_norm, l_fixed_nsubcls
        integer,                  intent(out)   :: nsplit
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        type(image), allocatable, intent(out)   :: raw_subavgs(:), den_subavgs(:)
        type(parameters) :: params_mask
        type(image), allocatable :: imgs(:), imgs_ppca(:), class_mask(:)
        type(image) :: cavg_raw, cavg_den
        real, allocatable :: avg(:), pcavecs(:,:), coords(:,:), eigvals(:), dmat(:,:)
        real, allocatable :: class_diams(:), class_shifts(:,:)
        integer, allocatable :: i_medoids(:)
        integer :: nptcls, npix, nsplit_count, nsplit_spec, j, k, class_ldim(3)
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
        call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
        call make_split_embedding(params, cls_id, nptcls, npix, pcavecs, coords, eigvals)
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
            call cluster_dmat(dmat, 'hclust', nsplit, i_medoids, labels)
        else
            nsplit_count = particle_count_nsplit(nptcls, params%nptcls_per_subcls, params%nsubcls_min, params%nsubcls_max)
            if( allocated(eigvals) )then
                nsplit_spec = spectral_nsplit(eigvals, params%nsubcls_min, params%nsubcls_max)
                nsplit = nsplit_spec
                print *, 'SUGGESTED NSPLIT FROM SPECTRUM: ', nsplit_spec
            else
                nsplit_spec = params%nsubcls_max
                nsplit      = max(params%nsubcls_min, min(params%nsubcls_max, nsplit_count))
            endif
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8)') 'Cls split auto ncls: class=', cls_id, &
                ' size=', nptcls, ' count_suggest=', nsplit_count, ' spectral_suggest=', nsplit_spec, ' selected=', nsplit
            call flush(logfhandle)
            call cluster_dmat(dmat, 'hclust', nsplit, i_medoids, labels)
        endif
        call cavg_den%new(cavg_raw%get_ldim(), cavg_raw%get_smpd())
        allocate(raw_subavgs(nsplit), den_subavgs(nsplit))
        do j = 1, nsplit
            call cavg_raw%zero_and_unflag_ft
            call cavg_den%zero_and_unflag_ft
            do k = 1, size(labels)
                if( labels(k) /= j ) cycle
                call cavg_raw%add(imgs(k))
                call cavg_den%add(imgs_ppca(k))
            end do
            call cavg_raw%div(real(max(count(labels == j), 1)))
            call cavg_den%div(real(max(count(labels == j), 1)))
            call raw_subavgs(j)%copy(cavg_raw)
            call den_subavgs(j)%copy(cavg_den)
        end do
        call cavg_raw%kill
        call cavg_den%kill
        if( allocated(coords)       ) deallocate(coords)
        if( allocated(eigvals)      ) deallocate(eigvals)
        if( allocated(dmat)         ) deallocate(dmat)
        if( allocated(i_medoids)    ) deallocate(i_medoids)
        if( allocated(avg)          ) deallocate(avg)
        if( allocated(pcavecs)      ) deallocate(pcavecs)
        if( allocated(class_diams)  ) deallocate(class_diams)
        if( allocated(class_shifts) ) deallocate(class_shifts)
        if( allocated(class_mask)   ) call dealloc_imgarr(class_mask)
        if( allocated(imgs_ppca)    ) call dealloc_imgarr(imgs_ppca)
        if( allocated(imgs)         ) call dealloc_imgarr(imgs)
    end subroutine split_one_parent_class

    subroutine make_split_embedding(params, cls_id, nptcls, npix, pcavecs, coords, eigvals)
        use simple_diffusion_maps, only: diffusion_map_embedder
        use simple_kpca_svd,       only: kpca_svd
        use simple_ppca,           only: ppca
        type(parameters),  intent(inout) :: params
        integer,           intent(in)    :: cls_id, nptcls, npix
        real,              intent(in)    :: pcavecs(npix,nptcls)
        real, allocatable, intent(out)   :: coords(:,:)
        real, allocatable, intent(out)   :: eigvals(:)
        type(diffusion_map_embedder) :: diffmap
        type(kpca_svd) :: kpca_model
        type(ppca)     :: ppca_model
        real, allocatable :: feat(:)
        integer :: neigs, i
        neigs = params%neigs
        if( neigs <= 0 )then
            select case(trim(params%pca_mode))
                case('diffusion_maps')
                    neigs = 0
                case DEFAULT
                    neigs = min(5, max(1, nptcls - 1))
                    write(logfhandle,'(A,A,A,I8,A,I8,A,I8)') 'Cls split auto-selected embedding dims: mode=', trim(params%pca_mode), &
                        ' class=', cls_id, ' size=', nptcls, ' neigs=', neigs
                    call flush(logfhandle)
            end select
        endif
        if( trim(params%pca_mode) .eq. 'diffusion_maps' )then
            neigs = min(max(neigs, 0), max(nptcls-2, 1))
        else
            neigs = min(max(neigs, 1), max(nptcls-1, 1))
        endif
        select case(trim(params%pca_mode))
            case('diffusion_maps')
                call diffmap%set_params(neigs, min(max(2, params%k_nn), max(2, nptcls-1)))
                call diffmap%embed(pcavecs, coords, eigvals)
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
            case DEFAULT
                THROW_HARD('cls_split pca_mode must be ppca, kpca, or diffusion_maps')
        end select
        write(logfhandle,'(A,A,A,I8,A,I8,A,I8)') 'Cls split embedding: mode=', trim(params%pca_mode), &
            ' class=', cls_id, ' size=', nptcls, ' dims=', size(coords,1)
        call flush(logfhandle)
    end subroutine make_split_embedding

    integer function particle_count_nsplit(nptcls, target_pop, nmin, nmax) result(nsplit)
        integer, intent(in) :: nptcls, target_pop, nmin, nmax
        integer :: nlo, nhi, target
        nlo    = max(1, nmin)
        nhi    = max(nlo, nmax)
        target = max(1, target_pop)
        nsplit = ceiling(real(nptcls) / real(target))
        nsplit = max(nlo, min(nhi, nsplit))
    end function particle_count_nsplit

    integer function spectral_nsplit(eigvals, nmin, nmax) result(nsplit)
        real,    intent(in) :: eigvals(:)
        integer, intent(in) :: nmin, nmax
        real :: best_gap, gap
        integer :: i, nlo, nhi, maxcand
        nlo = max(1, nmin)
        nhi = max(nlo, nmax)
        nsplit = nlo
        maxcand = min(size(eigvals), nhi)
        if( maxcand <= 1 ) return
        best_gap = -huge(best_gap)
        do i = 1, maxcand - 1
            gap = eigvals(i) - eigvals(i+1)
            if( gap > best_gap )then
                best_gap = gap
                nsplit = i + 1
            endif
        end do
        nsplit = max(nlo, min(nhi, nsplit))
    end function spectral_nsplit

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
        type(image)      :: img
        integer, allocatable :: part_counts(:), part_localstack(:,:), part_parent(:,:), part_local(:,:), part_pop(:,:), part_global(:,:)
        integer, allocatable :: comb_part(:), comb_row(:), comb_parent(:), comb_local(:), comb_pop(:)
        integer, allocatable :: new_class(:), new_parent(:)
        type(string) :: map_fname, assign_fname, raw_fname, den_fname
        integer :: ipart, nlocal, total, max_count, idx, i, funit, ios, pind, parent_cls, local_cls, global_cls, ldim(3), nstk
        real    :: smpd
        character(len=XLONGSTRLEN) :: line
        call spproj%read(params%projfile)
        allocate(part_counts(nparts_run), source=0)
        total = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_split_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call count_data_lines(map_fname, nlocal)
            part_counts(ipart) = nlocal
            total = total + nlocal
            call map_fname%kill
        end do
        if( total < 1 ) THROW_HARD('No subclass outputs produced by distributed cls_split workers')
        max_count = maxval(part_counts)
        allocate(part_localstack(max_count, nparts_run), part_parent(max_count, nparts_run), part_local(max_count, nparts_run), &
                 part_pop(max_count, nparts_run), part_global(max_count, nparts_run), source=0)
        allocate(comb_part(total), comb_row(total), comb_parent(total), comb_local(total), comb_pop(total), source=0)
        idx = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_split_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call read_part_map(map_fname, part_counts(ipart), part_localstack(1:part_counts(ipart), ipart), &
                               part_parent(1:part_counts(ipart), ipart), part_local(1:part_counts(ipart), ipart), &
                               part_pop(1:part_counts(ipart), ipart))
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
        call sort_combined_maps(comb_part, comb_row, comb_parent, comb_local, comb_pop)
        call del_file(string('cls_split_subclass_avgs.mrcs'))
        call del_file(string('cls_split_subclass_avgs_conditioned.mrcs'))
        raw_fname = string('cls_split_subclass_avgs_part')//int2str_pad(comb_part(1), params%numlen)//'.mrcs'
        call find_ldim_nptcls(raw_fname, ldim, nstk, smpd=smpd)
        ldim(3) = 1
        call img%new(ldim, smpd)
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
            raw_fname = string('cls_split_subclass_avgs_part')//int2str_pad(comb_part(idx), params%numlen)//'.mrcs'
            den_fname = string('cls_split_subclass_avgs_conditioned_part')//int2str_pad(comb_part(idx), params%numlen)//'.mrcs'
            call img%read(raw_fname, part_localstack(comb_row(idx), comb_part(idx)))
            call img%write(string('cls_split_subclass_avgs.mrcs'), idx)
            call img%read(den_fname, part_localstack(comb_row(idx), comb_part(idx)))
            call img%write(string('cls_split_subclass_avgs_conditioned.mrcs'), idx)
        end do
        close(funit)
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
        call img%kill
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
