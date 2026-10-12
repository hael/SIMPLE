!@descr: cls_expansion execution strategies: shared memory, distributed master and distributed worker
! ncls is the total number of subclasses: the master (or the shared-memory run) groups the classes by
! class-average correlation and budgets the subclasses per class; the restoration pass reuses the split's
! soft weights and class map.
module simple_cls_expansion_strategy
use simple_core_module_api
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder,                only: builder
use simple_parameters,             only: parameters
use simple_cmdline,                only: cmdline
use simple_qsys_env,               only: qsys_env
use simple_sp_project,             only: sp_project
use simple_image,                  only: image
use simple_classaverager,          only: transform_ptcls
use simple_imgarr_utils,           only: alloc_imgarr, dealloc_imgarr
use simple_srch_sort_loc,          only: hpsort
implicit none

public :: cls_expansion_strategy
public :: FLEX_CLS_CAVGS_FILE, FLEX_CLS_RESTORE_FLAG
public :: cls_expansion_shmem_strategy
public :: cls_expansion_worker_strategy
public :: cls_expansion_master_strategy
public :: create_cls_expansion_strategy

private
#include "simple_local_flags.inc"

integer, parameter :: FLEX_CLS_NEIGS_DEFAULT  = 8    !< rank of the per-class covariance model
integer, parameter :: FLEX_CLS_SEED_BASE      = 20260930 !< per-class random seed = base + class index
real,    parameter :: FLEX_CLS_LP_FIT_DEFAULT = 30.  !< the fit band stops here (A) unless lp= is given;
                                                     !! the sub-class averages always use the full band
character(len=*), parameter :: FLEX_CLS_WEIGHTS_FILE = 'cls_expansion_weights.txt'    !< pind parent subclass nsub w(1:nsub)
character(len=*), parameter :: FLEX_CLS_MAP_FILE     = 'cls_expansion_class_map.txt'  !< subclass parent local pop neff sep repro
character(len=*), parameter :: FLEX_CLS_ALLOC_FILE   = 'cls_expansion_allocation.txt' !< class group pop nsub: the split's budget
character(len=*), parameter :: FLEX_CLS_GROUPS_FBODY = 'cls_expansion_groups'         !< class-average groups: stacks and table
character(len=*), parameter :: FLEX_CLS_CAVGS_FILE   = 'cls_expansion_cavgs.mrc'      !< weighted sub-class averages
character(len=*), parameter :: FLEX_CLS_CAVGS_EVEN   = 'cls_expansion_cavgs_even.mrc' !< the same from the even members
character(len=*), parameter :: FLEX_CLS_CAVGS_ODD    = 'cls_expansion_cavgs_odd.mrc'  !< the same from the odd members
character(len=*), parameter :: FLEX_CLS_RESTORE_FLAG = 'cls_expansion_restore.txt'    !< present: restore from the split's weights

!> the split's soft weights of every particle (restoration input), indexed by particle
type split_weights
    integer, allocatable :: off(:)     !< first weight of the particle in w, 0 without a row
    integer, allocatable :: nsub(:)    !< subclasses of the particle's parent
    integer, allocatable :: sub(:)     !< the particle's hard subclass (global index)
    integer, allocatable :: parent(:)  !< the particle's parent class
    real,    allocatable :: w(:)
end type split_weights

type, abstract :: cls_expansion_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type cls_expansion_strategy

type, extends(cls_expansion_strategy) :: cls_expansion_shmem_strategy
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type cls_expansion_shmem_strategy

type, extends(cls_expansion_strategy) :: cls_expansion_worker_strategy
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
end type cls_expansion_worker_strategy

type, extends(cls_expansion_strategy) :: cls_expansion_master_strategy
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
end type cls_expansion_master_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: cls_expansion_strategy, parameters, builder, cmdline
        class(cls_expansion_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: cls_expansion_strategy, parameters, builder, cmdline
        class(cls_expansion_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: cls_expansion_strategy, parameters, builder, cmdline
        class(cls_expansion_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params)
        import :: cls_expansion_strategy, parameters
        class(cls_expansion_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    function create_cls_expansion_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(cls_expansion_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker, is_master
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        is_worker = cline%defined('part')
        is_master = (nparts > 1) .and. (.not. is_worker)
        if( is_master )then
            allocate(cls_expansion_master_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED CLS_EXPANSION (MASTER)'
        else if( is_worker )then
            allocate(cls_expansion_worker_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> CLS_EXPANSION (WORKER)'
        else
            allocate(cls_expansion_shmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> CLS_EXPANSION (SHARED-MEMORY)'
        endif
    end function create_cls_expansion_strategy

    subroutine apply_defaults(cline)
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')   ) call cline%set('neigs',   FLEX_CLS_NEIGS_DEFAULT)
        if( .not. cline%defined('objfun')  ) call cline%set('objfun',  'euclid')
    end subroutine apply_defaults

    subroutine init_common(params, build, cline)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call apply_defaults(cline)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call validate_cls_expansion(params, cline)
    end subroutine init_common

    subroutine validate_cls_expansion(params, cline)
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        if( trim(params%oritype) /= 'ptcl2D' ) THROW_HARD('cls_expansion supports oritype=ptcl2D only')
        if( .not. cline%defined('ncls') .or. params%ncls < 1 ) THROW_HARD('cls_expansion requires ncls, the total number of subclasses')
        if( params%neigs  < 1 ) THROW_HARD('cls_expansion requires neigs>=1')
        if( params%minpop < 2 ) THROW_HARD('cls_expansion requires minpop>=2')
        if( params%nclust < 1 ) THROW_HARD('cls_expansion requires nclust>=1')
    end subroutine validate_cls_expansion

    subroutine shmem_initialize(self, params, build, cline)
        class(cls_expansion_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call init_common(params, build, cline)
    end subroutine shmem_initialize

    subroutine shmem_execute(self, params, build, cline)
        class(cls_expansion_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_local_split(params, build, cline, spproj, part=0, l_write_project=.true.)
        call spproj%kill
    end subroutine shmem_execute

    subroutine shmem_finalize_run(self, params, build, cline)
        class(cls_expansion_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup(self, params)
        class(cls_expansion_shmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
    end subroutine shmem_cleanup

    subroutine worker_initialize(self, params, build, cline)
        class(cls_expansion_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call init_common(params, build, cline)
        if( .not. cline%defined('part') ) THROW_HARD('PART must be defined for cls_expansion worker execution')
        if( .not. cline%defined('class_assignment') ) THROW_HARD('CLASS_ASSIGNMENT must be defined for cls_expansion worker execution')
    end subroutine worker_initialize

    subroutine worker_execute(self, params, build, cline)
        use simple_qsys_funs, only: qsys_declare_part_finished
        class(cls_expansion_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_local_split(params, build, cline, spproj, part=params%part, l_write_project=.false.)
        call qsys_declare_part_finished(params, string('simple_cls_expansion_strategy :: worker_execute'))
        call spproj%kill
    end subroutine worker_execute

    subroutine worker_finalize_run(self, params, build, cline)
        use simple_qsys_funs, only: qsys_declare_part_finished
        class(cls_expansion_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call qsys_declare_part_finished(params, string('simple_commanders_denoise :: exec_cls_expansion'))
    end subroutine worker_finalize_run

    subroutine worker_cleanup(self, params)
        class(cls_expansion_worker_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
    end subroutine worker_cleanup

    subroutine master_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(cls_expansion_master_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(sp_project)     :: spproj
        integer, allocatable :: cls_inds(:), cls_pops(:), nsub(:), first_sub(:)
        real,    allocatable :: sep(:)
        call init_common(params, build, cline)
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        self%nparts_run = min(max(1, params%nparts), size(cls_inds))
        write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls expansion worker planning: requested_nparts=', params%nparts, &
            ' eligible_classes=', size(cls_inds), ' effective_nparts=', self%nparts_run
        call flush(logfhandle)
        if( self%nparts_run < params%nparts )then
            write(logfhandle,'(A)') 'Cls expansion note: reducing worker count because classes are the scheduling unit.'
            call flush(logfhandle)
        endif
        call set_master_num_threads(self%nthr_master, string('CLS_EXPANSION'))
        if( flex_restore_mode() )then
            ! the workers restore subsets: the whole set of parents is checked against the split here
            call read_split_map(cls_inds, .true., nsub, first_sub, sep)
            deallocate(nsub, first_sub, sep)
        else
            ! the subclass budget is set once, here, and read by the workers
            call spproj%read(params%projfile)
            call allocate_subclasses(params, spproj, cls_inds, cls_pops, self%nthr_master, nsub)
            call spproj%kill
            deallocate(nsub)
        endif
        call prepare_class_partitions(self, params, cline, cls_inds, cls_pops)
        call self%qenv%new(params, self%nparts_run, numlen=params%numlen)
        ! workers run inside the master's execution directory
        call cline%set('mkdir', 'no')
        call cline%gen_job_descr(self%job_descr)
        call self%job_descr%set('mkdir', 'no')
        call self%job_descr%set('nparts', int2str(self%nparts_run))
        call self%job_descr%set('numlen', int2str(params%numlen))
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
    end subroutine master_initialize

    subroutine master_execute(self, params, build, cline)
        class(cls_expansion_master_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, part_params=self%part_params, &
                                                     array=L_USE_SLURM_ARR, extra_params=params)
        call merge_worker_outputs(params, self%nparts_run)
    end subroutine master_execute

    subroutine master_finalize_run(self, params, build, cline)
        class(cls_expansion_master_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
    end subroutine master_finalize_run

    subroutine master_cleanup(self, params)
        use simple_qsys_funs, only: qsys_cleanup
        class(cls_expansion_master_strategy), intent(inout) :: self
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
        class(cls_expansion_master_strategy), intent(inout) :: self
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
            fname = string('cls_expansion_classes_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
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
        ! the distinct values of the label among the active particles
        cls_inds = split_label_inds(build, label)
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
        if( size(cls_inds) < 1 ) THROW_HARD('No classes selected for cls_expansion')
        allocate(cls_pops(size(cls_inds)), source=0)
        do i = 1, size(cls_inds)
            call build%spproj_field%get_pinds(cls_inds(i), label%to_char(), pinds)
            if( allocated(pinds) )then
                cls_pops(i) = size(pinds)
                deallocate(pinds)
            endif
        end do
        ! every labelled class is a parent: one too small to split is kept whole
        cls_inds = pack(cls_inds, mask=cls_pops > 0)
        cls_pops = pack(cls_pops, mask=cls_pops > 0)
        if( size(cls_inds) < 1 ) THROW_HARD('No populated classes for cls_expansion')
        call label%kill
    end subroutine collect_split_classes

    function split_label_inds( build, label ) result( inds )
        type(builder), intent(inout) :: build
        type(string),  intent(in)    :: label
        integer, allocatable :: inds(:)
        logical, allocatable :: isthere(:)
        integer :: i, v, vmax
        vmax = 0
        do i = 1, build%spproj_field%get_noris()
            if( build%spproj_field%get_state(i) <= 0 ) cycle
            vmax = max(vmax, build%spproj_field%get_int(i, label%to_char()))
        end do
        allocate(isthere(max(vmax,1)), source=.false.)
        do i = 1, build%spproj_field%get_noris()
            if( build%spproj_field%get_state(i) <= 0 ) cycle
            v = build%spproj_field%get_int(i, label%to_char())
            if( v >= 1 ) isthere(v) = .true.
        end do
        inds = pack([(i, i=1,size(isthere))], mask=isthere)
    end function split_label_inds

    subroutine determine_split_label(params, build, label)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(string),     intent(out)   :: label
        label = 'class'
        ! restoration: the parents are the clusters the split wrote
        if( flex_restore_mode() ) label = 'cluster'
    end subroutine determine_split_label


    !> the classes of this run (all, or a worker's part): split under the budget, or restored from the split's
    !! weights; writes the class map, weights, stacks and, outside a worker, the project
    subroutine run_local_split(params, build, cline, spproj, part, l_write_project)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: part
        logical,          intent(in)    :: l_write_project
        type(split_weights)      :: sw
        type(string)             :: map_fname, assign_fname, wts_fname, stk_fname, stk_even, stk_odd
        type(image), allocatable :: cavgs(:), cavgs_even(:), cavgs_odd(:)
        integer,     allocatable :: cls_inds(:), cls_pops(:), nsub(:), first_sub(:), pinds(:), labels(:)
        integer,     allocatable :: new_class(:), new_parent(:), parent_of_subcls(:), pop_of_subcls(:)
        real,        allocatable :: weights(:,:), repro(:), sep_of_cls(:), neff_of_subcls(:), sep_of_subcls(:), repro_of_subcls(:)
        real    :: separation
        integer :: i, iglob, ntot_sub, funit_map, funit_assign, funit_wts, kfit
        logical :: l_sigma, l_restore
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        l_restore = flex_restore_mode()
        if( l_restore )then
            ! the split's record is read before any output of this pass replaces it
            call read_split_map(cls_inds, l_write_project, nsub, first_sub, sep_of_cls)
            call read_split_weights(spproj%os_ptcl2D%get_noris(), sw)
        else if( l_write_project )then
            call allocate_subclasses(params, spproj, cls_inds, cls_pops, params%nthr, nsub)
        else
            call read_allocation(cls_inds, nsub)
        endif
        ntot_sub = sum(nsub)
        call prepare_flex_sigma2(params, build, spproj, cline, l_sigma)
        kfit = flex_fit_band(params, cline)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,L1,A,L1)') 'Cls expansion flex: classes=', size(cls_inds), &
            &' subclasses=', ntot_sub, ' neigs=', params%neigs, ' fit band kto=', kfit, ' sigma2 loaded=', l_sigma, &
            &' restoration=', l_restore
        call flush(logfhandle)
        if( l_write_project )then
            map_fname = string(FLEX_CLS_MAP_FILE)
            wts_fname = string(FLEX_CLS_WEIGHTS_FILE)
            stk_fname = string(FLEX_CLS_CAVGS_FILE)
            stk_even  = string(FLEX_CLS_CAVGS_EVEN)
            stk_odd   = string(FLEX_CLS_CAVGS_ODD)
            allocate(new_class(spproj%os_ptcl2D%get_noris()), new_parent(spproj%os_ptcl2D%get_noris()), source=0)
            allocate(parent_of_subcls(ntot_sub), pop_of_subcls(ntot_sub), source=0)
            allocate(neff_of_subcls(ntot_sub), sep_of_subcls(ntot_sub), repro_of_subcls(ntot_sub), source=0.)
        else
            map_fname    = string('cls_expansion_class_map_part')//int2str_pad(part, params%numlen)//TXT_EXT
            assign_fname = string('cls_expansion_assignments_part')//int2str_pad(part, params%numlen)//TXT_EXT
            wts_fname    = flex_weights_part_fname(part, params%numlen)
            stk_fname    = flex_cavgs_part_fname(part, params%numlen)
            stk_even     = flex_cavgs_part_fname(part, params%numlen, 'even')
            stk_odd      = flex_cavgs_part_fname(part, params%numlen, 'odd')
            open(newunit=funit_assign, file=assign_fname%to_char(), status='replace', action='write')
            write(funit_assign,'(A)') '# particle_index parent_class local_subclass'
        endif
        ! in a part the subclass column is the part's row (map) or the local subclass (weights)
        open(newunit=funit_map, file=map_fname%to_char(), status='replace', action='write')
        write(funit_map,'(A)') '# subclass parent_class local_subclass pop neff separation repro'
        open(newunit=funit_wts, file=wts_fname%to_char(), status='replace', action='write')
        write(funit_wts,'(A)') '# particle_index parent_class subclass nsub weights(1:nsub)'
        iglob = 0
        do i = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8,A,I6)') 'Cls expansion starting: class=', cls_inds(i), ' nptcls=', cls_pops(i), &
                &' nsub=', nsub(i)
            call flush(logfhandle)
            if( l_restore )then
                if( l_write_project .and. iglob + 1 /= first_sub(i) ) THROW_HARD('cls_expansion restoration: subclass numbering differs from the split')
                separation = sep_of_cls(i)
                call restore_class_flex(params, build, spproj, cls_inds(i), nsub(i), first_sub(i), separation, kfit, l_sigma, &
                    &sw, pinds, labels, weights, cavgs, repro, cavgs_even, cavgs_odd)
            else
                call split_class_flex(params, build, spproj, cls_inds(i), nsub(i), params%neigs, kfit, l_sigma, &
                    &pinds, labels, weights, cavgs, separation, repro, cavgs_even, cavgs_odd)
            endif
            call record_class(cls_inds(i), nsub(i), pinds, labels, weights, cavgs, separation, repro, cavgs_even, cavgs_odd)
            if( allocated(pinds)      ) deallocate(pinds)
            if( allocated(repro)      ) deallocate(repro)
            if( allocated(labels)     ) deallocate(labels)
            if( allocated(weights)    ) deallocate(weights)
            if( allocated(cavgs)      ) call dealloc_imgarr(cavgs)
            if( allocated(cavgs_even) ) call dealloc_imgarr(cavgs_even)
            if( allocated(cavgs_odd)  ) call dealloc_imgarr(cavgs_odd)
        end do
        close(funit_map)
        close(funit_wts)
        if( .not. l_write_project ) close(funit_assign)
        if( l_write_project )then
            if( iglob /= ntot_sub ) THROW_HARD('cls_expansion: the delivered subclasses differ from the budget')
            call apply_split_project_updates(spproj, params, iglob, new_class, new_parent, parent_of_subcls, pop_of_subcls, &
                &neff_of_subcls, sep_of_subcls, repro_of_subcls, stk_fname)
            deallocate(new_class, new_parent, parent_of_subcls, pop_of_subcls, neff_of_subcls, sep_of_subcls, repro_of_subcls)
        endif
        call kill_split_weights(sw)
        call map_fname%kill
        call wts_fname%kill
        call stk_fname%kill
        call stk_even%kill
        call stk_odd%kill
        if( assign_fname%is_allocated() ) call assign_fname%kill

    contains

        !> one class into the map, weights, stacks and assignments (part) or project arrays (whole run)
        subroutine record_class( cls_id, nsplit, pinds, labels, weights, cavgs, separation, repro, cavgs_even, cavgs_odd )
            integer,                  intent(in)    :: cls_id, nsplit
            integer,                  intent(in)    :: pinds(:), labels(:)
            real,                     intent(in)    :: weights(:,:)
            type(image), allocatable, intent(inout) :: cavgs(:)
            real,                     intent(in)    :: separation
            real,        allocatable, intent(in)    :: repro(:)
            type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
            integer :: j, k, isub
            real    :: neff
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls expansion summary: class=', cls_id, ' nptcls=', size(labels), ' nsubcls=', nsplit
            call flush(logfhandle)
            do j = 1, nsplit
                iglob = iglob + 1
                neff  = sum(weights(:,j))
                if( l_write_project )then
                    parent_of_subcls(iglob) = cls_id
                    pop_of_subcls(iglob)    = count(labels == j)
                    neff_of_subcls(iglob)   = neff
                    sep_of_subcls(iglob)    = separation
                    repro_of_subcls(iglob)  = repro(j)
                endif
                write(funit_map,'(I8,1X,I8,1X,I8,1X,I8,1X,F12.3,1X,F8.3,1X,F8.3)') iglob, cls_id, j, count(labels == j), &
                    &neff, separation, repro(j)
                ! the sub-class average lands at its row: global index here, local row in a part
                call cavgs(j)%write(stk_fname, iglob)
                call cavgs_even(j)%write(stk_even, iglob)
                call cavgs_odd(j)%write(stk_odd, iglob)
            end do
            do k = 1, size(labels)
                if( labels(k) <= 0 ) cycle
                if( l_write_project )then
                    isub = iglob - nsplit + labels(k)
                    new_class(pinds(k))  = isub
                    new_parent(pinds(k)) = cls_id
                else
                    isub = labels(k)
                    write(funit_assign,'(I10,1X,I10,1X,I10)') pinds(k), cls_id, labels(k)
                endif
                write(funit_wts,'(I10,1X,I10,1X,I10,1X,I6,*(1X,F8.6))') pinds(k), cls_id, isub, nsplit, weights(k,:)
            end do
        end subroutine record_class

    end subroutine run_local_split

    !> The split's subclass budget: the classes grouped by class-average correlation (one group per class without
    !! class averages), then ncls subclasses spread over groups and classes; written to FLEX_CLS_ALLOC_FILE
    subroutine allocate_subclasses( params, spproj, cls_inds, cls_pops, nthr, nsub )
        use simple_view_partition_sampling, only: has_clusterable_cavgs, cluster_cavgs_into_groups
        use simple_flex_cls_expansion,      only: flex_cls_allocate_subclasses, flex_cls_split_floor
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: cls_inds(:), cls_pops(:), nthr
        integer, allocatable, intent(out)   :: nsub(:)
        integer, allocatable  :: groups(:)
        logical, allocatable  :: l_sel(:)
        character(len=STDLEN) :: msg
        integer :: m, ncls2D, minsplit, ngot, ngroups, funit, i, g
        m = size(cls_inds)
        if( params%ncls < m )then
            write(msg,'(A,I0,A,I0,A)') 'cls_expansion: ncls=', params%ncls, ' (the total number of subclasses) is below the ', &
                &m, ' classes'
            THROW_HARD(trim(msg))
        endif
        if( m > 1 )then
            if( any(cls_inds(2:) <= cls_inds(:m-1)) ) THROW_HARD('cls_expansion: the classes must be in increasing order')
        endif
        ncls2D = spproj%os_cls2D%get_noris()
        if( has_clusterable_cavgs(spproj) .and. maxval(cls_inds) <= ncls2D )then
            allocate(l_sel(ncls2D), source=.false.)
            l_sel(cls_inds) = .true.
            call cluster_cavgs_into_groups(params, spproj, nthr, l_sel, params%nclust, FLEX_CLS_GROUPS_FBODY, groups)
        else
            write(logfhandle,'(A)') '>>> CLS_EXPANSION: no class average for every class (or no frc2D); every class is its own group'
            groups = [(i, i=1,m)]
        endif
        ngroups  = maxval(groups)
        minsplit = max(2 * params%minpop, flex_cls_split_floor(2, params%neigs))
        allocate(nsub(m))
        call flex_cls_allocate_subclasses(cls_pops, groups, params%ncls, params%minpop, minsplit, nsub, ngot)
        if( ngot < params%ncls )then
            write(msg,'(A,I0,A,I0,A,I0)') 'cls_expansion: the class populations hold ', ngot, ' subclasses at minpop=', &
                &params%minpop, ' instead of ncls=', params%ncls
            THROW_WARN(trim(msg))
        endif
        open(newunit=funit, file=FLEX_CLS_ALLOC_FILE, status='replace', action='write')
        write(funit,'(A)') '# class group pop nsub'
        do i = 1, m
            write(funit,'(I8,1X,I6,1X,I10,1X,I6)') cls_inds(i), groups(i), cls_pops(i), nsub(i)
        end do
        close(funit)
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> CLS_EXPANSION SUBCLASS BUDGET: ', sum(nsub), ' subclasses of ', &
            &m, ' classes in ', ngroups, ' groups (ncls=', params%ncls, ', minpop=', params%minpop, ')'
        write(logfhandle,'(A,I0,A,I0)') '    Classes kept whole: ', count(nsub == 1), '; split threshold (particles): ', minsplit
        write(logfhandle,'(A)') '    Group   Classes    Particles   Subclasses   Particles/subclass'
        do g = 1, ngroups
            if( .not. any(groups == g) ) cycle
            write(logfhandle,'(4X,I5,2X,I8,2X,I11,2X,I11,2X,F19.1)') g, count(groups == g), sum(cls_pops, mask=(groups == g)), &
                &sum(nsub, mask=(groups == g)), real(sum(cls_pops, mask=(groups == g))) / real(sum(nsub, mask=(groups == g)))
        end do
        write(logfhandle,'(A)') ''
        call flush(logfhandle)
    end subroutine allocate_subclasses

    !> a worker's subclass counts, from the master's budget file
    subroutine read_allocation( cls_inds, nsub )
        integer,              intent(in)  :: cls_inds(:)
        integer, allocatable, intent(out) :: nsub(:)
        character(len=XLONGSTRLEN) :: line
        integer :: funit, ios, cls, grp, pop, n, i
        allocate(nsub(size(cls_inds)), source=0)
        open(newunit=funit, file=FLEX_CLS_ALLOC_FILE, status='old', action='read', iostat=ios)
        call fileiochk('read_allocation opening '//FLEX_CLS_ALLOC_FILE, ios)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            read(line,*) cls, grp, pop, n
            i = findloc(cls_inds, cls, dim=1)
            if( i > 0 ) nsub(i) = n
        end do
        close(funit)
        if( any(nsub < 1) ) THROW_HARD('a class of this part has no subclass count in '//FLEX_CLS_ALLOC_FILE)
    end subroutine read_allocation

    !> the split's class map: per class of cls_inds its subclass count, first global subclass and separation;
    !! l_all requires cls_inds to hold every parent of the map (a whole run, not a worker's part)
    subroutine read_split_map( cls_inds, l_all, nsub, first_sub, sep )
        integer,              intent(in)  :: cls_inds(:)
        logical,              intent(in)  :: l_all
        integer, allocatable, intent(out) :: nsub(:), first_sub(:)
        real,    allocatable, intent(out) :: sep(:)
        character(len=XLONGSTRLEN) :: line
        real    :: neff, sepv, repro
        integer :: funit, ios, isub, parent, loc, pop, i
        allocate(nsub(size(cls_inds)), first_sub(size(cls_inds)), source=0)
        allocate(sep(size(cls_inds)), source=0.)
        open(newunit=funit, file=FLEX_CLS_MAP_FILE, status='old', action='read', iostat=ios)
        call fileiochk('read_split_map opening '//FLEX_CLS_MAP_FILE, ios)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            read(line,*) isub, parent, loc, pop, neff, sepv, repro
            i = findloc(cls_inds, parent, dim=1)
            if( i == 0 .and. l_all ) THROW_HARD('a class of the split has no members to restore; its subclasses would be renumbered')
            if( i == 0 ) cycle
            ! the subclasses of a class are consecutive: global = first + local - 1
            if( nsub(i) == 0 ) first_sub(i) = isub - loc + 1
            if( isub - loc + 1 /= first_sub(i) ) THROW_HARD('the subclasses of a class are not consecutive in '//FLEX_CLS_MAP_FILE)
            nsub(i) = nsub(i) + 1
            sep(i)  = sepv
        end do
        close(funit)
        if( any(nsub < 1) ) THROW_HARD('a class to restore has no subclasses in '//FLEX_CLS_MAP_FILE)
    end subroutine read_split_map

    !> the split's soft weights of every particle with a row (nptcls = the particle field)
    subroutine read_split_weights( nptcls, sw )
        integer,             intent(in)    :: nptcls
        type(split_weights), intent(inout) :: sw
        character(len=:), allocatable :: line
        real,             allocatable :: wbuf(:)
        integer :: funit, ios, pind, parent, isub, n, ntot, nmax, j
        call kill_split_weights(sw)
        allocate(sw%off(nptcls), sw%nsub(nptcls), sw%sub(nptcls), sw%parent(nptcls), source=0)
        open(newunit=funit, file=FLEX_CLS_WEIGHTS_FILE, status='old', action='read', iostat=ios)
        call fileiochk('read_split_weights opening '//FLEX_CLS_WEIGHTS_FILE, ios)
        ! rows are as wide as their parent has subclasses: count first, then read
        ntot = 0
        nmax = 0
        do
            call read_record(funit, line, ios)
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            read(line,*,iostat=ios) pind, parent, isub, n
            if( ios /= 0 ) THROW_HARD('read error in '//FLEX_CLS_WEIGHTS_FILE)
            if( pind < 1 .or. pind > nptcls .or. n < 1 ) THROW_HARD('invalid row in '//FLEX_CLS_WEIGHTS_FILE)
            ntot = ntot + n
            nmax = max(nmax, n)
        end do
        allocate(sw%w(max(1,ntot)), source=0.)
        allocate(wbuf(max(1,nmax)))
        rewind(funit)
        ntot = 0
        do
            call read_record(funit, line, ios)
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            read(line,*,iostat=ios) pind, parent, isub, n, (wbuf(j), j=1,n)
            if( ios /= 0 ) THROW_HARD('short or invalid row in '//FLEX_CLS_WEIGHTS_FILE)
            sw%off(pind)        = ntot + 1
            sw%nsub(pind)       = n
            sw%sub(pind)        = isub
            sw%parent(pind)     = parent
            sw%w(ntot+1:ntot+n) = wbuf(1:n)
            ntot = ntot + n
        end do
        close(funit)
        deallocate(wbuf)
    end subroutine read_split_weights

    !> one whole record of a formatted sequential file, however long (ios /= 0 at the end of the file)
    subroutine read_record( funit, line, ios )
        integer,                       intent(in)  :: funit
        character(len=:), allocatable, intent(out) :: line
        integer,                       intent(out) :: ios
        character(len=1024) :: buf
        integer :: nread
        line = ''
        do
            read(funit,'(A)',advance='no',iostat=ios,size=nread) buf
            line = line//buf(1:nread)
            if( is_iostat_eor(ios) )then
                ios = 0
                return
            endif
            if( ios /= 0 ) return
        end do
    end subroutine read_record

    subroutine kill_split_weights( sw )
        type(split_weights), intent(inout) :: sw
        if( allocated(sw%off)    ) deallocate(sw%off)
        if( allocated(sw%nsub)   ) deallocate(sw%nsub)
        if( allocated(sw%sub)    ) deallocate(sw%sub)
        if( allocated(sw%parent) ) deallocate(sw%parent)
        if( allocated(sw%w)      ) deallocate(sw%w)
    end subroutine kill_split_weights


    ! ===== the flex model

    function flex_weights_part_fname( part, numlen ) result( fname )
        integer, intent(in) :: part, numlen
        type(string) :: fname
        fname = string('cls_expansion_weights_part')//int2str_pad(part, numlen)//TXT_EXT
    end function flex_weights_part_fname

    function flex_cavgs_part_fname( part, numlen, which ) result( fname )
        integer,                    intent(in) :: part, numlen
        character(len=*), optional, intent(in) :: which   !< 'even' or 'odd' for the half stacks
        type(string) :: fname
        if( present(which) )then
            fname = string('cls_expansion_cavgs_'//trim(which)//'_part')//int2str_pad(part, numlen)//MRC_EXT
        else
            fname = string('cls_expansion_cavgs_part')//int2str_pad(part, numlen)//MRC_EXT
        endif
    end function flex_cavgs_part_fname

    !> the last Fourier shell of the covariance fit (lp= or the default, capped at Nyquist)
    integer function flex_fit_band( params, cline ) result( kfit )
        type(parameters), intent(in) :: params
        class(cmdline),   intent(in) :: cline
        real    :: lp_fit
        integer :: nyq
        nyq    = fdim(params%box) - 1
        lp_fit = FLEX_CLS_LP_FIT_DEFAULT
        if( cline%defined('lp') ) lp_fit = params%lp
        lp_fit = max(lp_fit, 2. * params%smpd)
        kfit   = max(2, min(nyq, calc_fourier_index(lp_fit, params%box, params%smpd)))
    end function flex_fit_band

    !> the canonical sigma2 spectra of every particle, loaded over the whole particle field
    subroutine prepare_flex_sigma2( params, build, spproj, cline, loaded )
        use simple_sigma2_files, only: canonical_sigma2_consumable, load_sigma2_groups
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        class(cmdline),   intent(inout) :: cline
        logical,          intent(out)   :: loaded
        character(len=STDLEN) :: message
        integer :: fromp_bak, top_bak
        loaded = .false.
        if( params%cc_objfun /= OBJFUN_EUCLID )then
            write(logfhandle,'(A)') 'Cls expansion flex: objfun=cc, the covariance fit uses unit noise weights'
            return
        endif
        if( .not. canonical_sigma2_consumable(spproj, spproj%os_ptcl2D, params%box, params%smpd, params%l_sigma_glob, message) )then
            ! either grouping policy's spectrum serves as the noise weight
            if( canonical_sigma2_consumable(spproj, spproj%os_ptcl2D, params%box, params%smpd, .not. params%l_sigma_glob, message) )then
                params%l_sigma_glob = .not. params%l_sigma_glob
                if( params%l_sigma_glob )then
                    params%sigma_est = 'global'
                else
                    params%sigma_est = 'group'
                endif
                write(logfhandle,'(A)') 'Cls expansion flex: canonical sigma2 state found under sigma_est='//trim(params%sigma_est)
            else
                write(logfhandle,'(A)') 'Cls expansion flex WARNING: no consumable canonical sigma2 state ('//trim(message)// &
                    &'); the covariance fit uses the class residual spectrum as noise weights. Run calc_pspec &
                    &(or a euclid 2D/3D pass) on the project first.'
                return
            endif
        endif
        fromp_bak = params%fromp
        top_bak   = params%top
        params%fromp = 1
        params%top   = spproj%os_ptcl2D%get_noris()
        call load_sigma2_groups(params, build%pftc, build%esig, spproj, spproj%os_ptcl2D, loaded)
        params%fromp = fromp_bak
        params%top   = top_bak
    end subroutine prepare_flex_sigma2

    !> the members of one class as Fourier planes in the class frame with their CTF, half-plane lattice, shells 1..nyq;
    !! noise weights from the canonical sigma2 or the class's own residual spectrum; nptcls = 0 without members
    subroutine flex_class_planes( params, build, spproj, cls_id, kfit, l_sigma, pinds, y, c, w, wq, &
                                  hidx, kidx, shell, fitmask, gridcorr_img, ldim, ncoeff, nptcls )
        use simple_flex_cls_expansion,  only: flex_cls_shell_noise
        use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps, ft_map_lims, ft_map_phys_addrh, ft_map_phys_addrk
        use simple_ctf,             only: ctf
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id, kfit
        logical,                  intent(in)    :: l_sigma
        integer,     allocatable, intent(out)   :: pinds(:), hidx(:), kidx(:), shell(:)
        complex(sp), allocatable, intent(out)   :: y(:,:)
        real(sp),    allocatable, intent(out)   :: c(:,:), w(:,:), wq(:)
        logical,     allocatable, intent(out)   :: fitmask(:)
        type(image),              intent(inout) :: gridcorr_img
        integer,                  intent(out)   :: ldim(3), ncoeff, nptcls
        type(image), allocatable :: imgs(:)
        type(ctfparams)      :: ctfparms
        type(ctf)            :: tfun
        type(ctfvars)        :: ctfvals
        real(sp),    allocatable :: s2(:)
        real    :: mat(2,2), loc(2), sfsq, ang, df, tval, e3, sum_df, diff_df, angast, phshift, amp_contr, wl, half_wl2_cs
        integer :: lims(3,2), nyq, i, j, h, k, sh, iptcl
        logical :: l_ctf, l_flip
        type(string) :: label
        integer, allocatable :: members(:)
        nptcls = 0
        ncoeff = 0
        if( allocated(pinds) ) deallocate(pinds)
        ! the members by the split label (class, or cluster when the labels come from the project)
        call determine_split_label(params, build, label)
        call spproj%os_ptcl2D%get_pinds(cls_id, label%to_char(), members)
        call label%kill
        if( .not. allocated(members) ) return
        call transform_ptcls(params, build, spproj, params%oritype, cls_id, imgs, pinds, phflip=.false., &
            &keep_ft=.true., gridcorr=gridcorr_img, pinds_in=members)
        deallocate(members)
        if( .not. allocated(imgs) )then
            write(logfhandle,'(A,I8)') 'Cls expansion flex warning: no transformed images returned for parent class ', cls_id
            call flush(logfhandle)
            return
        endif
        nptcls = size(imgs)
        ldim = [params%box, params%box, 1]
        call memoize_ft_maps(ldim, params%smpd)
        lims = ft_map_lims
        nyq  = fdim(params%box) - 1
        do h = 0, lims(1,2)
            do k = lims(2,1), lims(2,2)
                if( h == 0 .and. k < 0 ) cycle
                sh = nint(hyp(h,k))
                if( sh < 1 .or. sh > nyq ) cycle
                ncoeff = ncoeff + 1
            end do
        end do
        allocate(hidx(ncoeff), kidx(ncoeff), shell(ncoeff), wq(ncoeff), fitmask(ncoeff))
        j = 0
        do h = 0, lims(1,2)
            do k = lims(2,1), lims(2,2)
                if( h == 0 .and. k < 0 ) cycle
                sh = nint(hyp(h,k))
                if( sh < 1 .or. sh > nyq ) cycle
                j = j + 1
                hidx(j)  = h
                kidx(j)  = k
                shell(j) = sh
                wq(j)    = 2.
                if( h == 0 .and. k == 0 ) wq(j) = 1.
                fitmask(j) = sh <= kfit
            end do
        end do
        allocate(y(ncoeff,nptcls), c(ncoeff,nptcls), w(ncoeff,nptcls))
        select case( spproj%get_ctfflag_type(params%oritype, pinds(1)) )
            case(CTFFLAG_NO)
                l_ctf = .false.; l_flip = .false.
            case(CTFFLAG_FLIP)
                l_ctf = .true.;  l_flip = .true.
            case(CTFFLAG_YES)
                l_ctf = .true.;  l_flip = .false.
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG; flex_class_planes')
        end select
        !$omp parallel do default(shared) schedule(static) proc_bind(close) &
        !$omp private(i,j,iptcl,e3,mat,loc,sfsq,ang,df,tval,ctfparms,tfun,ctfvals,sum_df,diff_df,angast,phshift,amp_contr,wl,half_wl2_cs)
        do i = 1, nptcls
            iptcl = pinds(i)
            e3    = spproj%os_ptcl2D%e3get(iptcl)
            ! the CTF at the particle-frame location of the class-frame coefficient
            call rotmat2d(-e3, mat)
            if( l_ctf )then
                ctfparms = spproj%get_ctfparams(params%oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
                ctfvals     = tfun%get_ctfvars(ctfparms%phshift)
                wl          = ctfvals%wl
                half_wl2_cs = 0.5 * wl * wl * ctfvals%cs
                sum_df      = ctfvals%dfx + ctfvals%dfy
                diff_df     = ctfvals%dfx - ctfvals%dfy
                angast      = ctfvals%angast
                phshift     = ctfvals%phshift
                amp_contr   = ctfvals%amp_contr_const
            endif
            do j = 1, ncoeff
                y(j,i) = imgs(i)%get_cmat_at(ft_map_phys_addrh(hidx(j),kidx(j)), ft_map_phys_addrk(hidx(j),kidx(j)), 1)
                if( l_ctf )then
                    loc  = matmul(real([hidx(j),kidx(j)]), mat)
                    sfsq = (loc(1)/real(ldim(1)))**2 + (loc(2)/real(ldim(2)))**2
                    ang  = atan2(loc(2), loc(1))
                    df   = 0.5 * (sum_df + cos(2.0 * (ang - angast)) * diff_df)
                    tval = sin(PI * wl * sfsq * (df - half_wl2_cs * sfsq) + phshift + amp_contr)
                    if( l_flip ) tval = abs(tval)
                    c(j,i) = tval
                else
                    c(j,i) = 1.
                endif
                if( l_sigma )then
                    w(j,i) = 1. / max(build%esig%sigma2_noise(shell(j), iptcl), TINY)
                else
                    w(j,i) = 1.
                endif
            end do
        end do
        !$omp end parallel do
        call dealloc_imgarr(imgs)
        call forget_ft_maps
        if( .not. l_sigma )then
            ! no canonical sigma2: the class's own residual spectrum sets the noise weights
            allocate(s2(nyq))
            call flex_cls_shell_noise(y, c, shell, nyq, s2)
            do i = 1, nptcls
                do j = 1, ncoeff
                    w(j,i) = 1. / s2(shell(j))
                end do
            end do
            write(logfhandle,'(A,I8,A,ES10.3,A,ES10.3)') 'Cls expansion flex noise from residuals: class=', cls_id, &
                &' sigma2 shell 1=', s2(1), ' shell kfit=', s2(kfit)
            deallocate(s2)
        endif
    end subroutine flex_class_planes

    !> placement into ncls subclasses and the weighted averages of one class, from its fitted model
    subroutine flex_class_deliver( params, cls_id, model, ncls, y, c, w, wq, fitmask, shell, nsh, hidx, kidx, gridcorr_img, ldim, &
                                   ncoeff, nptcls, labels, weights, cavgs, separation, repro, cavgs_even, cavgs_odd )
        use simple_flex_cls_expansion,  only: flex_cls_model, flex_cls_place_states
        type(parameters),         intent(inout) :: params
        integer,                  intent(in)    :: cls_id, ncls, ldim(3), ncoeff, nptcls, nsh
        type(flex_cls_model),     intent(in)    :: model
        complex(sp),              intent(in)    :: y(:,:)
        real(sp),                 intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,                  intent(in)    :: fitmask(:)
        integer,                  intent(in)    :: shell(:), hidx(:), kidx(:)
        type(image),              intent(inout) :: gridcorr_img
        integer,     allocatable, intent(out)   :: labels(:)
        real,        allocatable, intent(out)   :: weights(:,:)
        type(image), allocatable, intent(inout) :: cavgs(:)
        real,                     intent(out)   :: separation
        real,        allocatable, intent(out)   :: repro(:)
        type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
        real(sp), allocatable :: neff(:)
        allocate(weights(nptcls,ncls), labels(nptcls), neff(ncls))
        call flex_cls_place_states(model, ncls, weights, labels, neff, separation)
        call flex_class_restore_deliver(params, cls_id, ncls, y, c, w, wq, fitmask, shell, nsh, hidx, kidx, gridcorr_img, ldim, &
            &ncoeff, nptcls, labels, weights, cavgs, repro, cavgs_even, cavgs_odd, neff=neff, separation=separation)
        deallocate(neff)
    end subroutine flex_class_deliver

    !> the restoration of one class from its weights: regularised sub-averages, even/odd versions, reproducibility
    subroutine flex_class_restore_deliver( params, cls_id, ncls, y, c, w, wq, fitmask, shell, nsh, hidx, kidx, gridcorr_img, &
                                           ldim, ncoeff, nptcls, labels, weights, cavgs, repro, cavgs_even, cavgs_odd, neff, separation )
        use simple_flex_cls_expansion,  only: flex_cls_restore_states, flex_cls_half_reproducibility, flex_cls_signal_power
        use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps, ft_map_phys_addrh, ft_map_phys_addrk
        type(parameters),         intent(inout) :: params
        integer,                  intent(in)    :: cls_id, ncls, ldim(3), ncoeff, nptcls, nsh
        complex(sp),              intent(in)    :: y(:,:)
        real(sp),                 intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,                  intent(in)    :: fitmask(:)
        integer,                  intent(in)    :: shell(:), hidx(:), kidx(:)
        type(image),              intent(inout) :: gridcorr_img
        integer,                  intent(in)    :: labels(:)
        real,                     intent(in)    :: weights(:,:)
        type(image), allocatable, intent(inout) :: cavgs(:)
        real,        allocatable, intent(out)   :: repro(:)
        type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
        real,        optional,    intent(in)    :: neff(:), separation
        complex(sp), allocatable :: avgs(:,:), avgs_e(:,:), avgs_o(:,:)
        real(dp),    allocatable :: tau2(:,:), den_e(:,:), den_o(:,:)
        real    :: nf, sep
        integer :: s
        allocate(repro(ncls))
        allocate(avgs(ncoeff,ncls), avgs_e(ncoeff,ncls), avgs_o(ncoeff,ncls), tau2(nsh,ncls), den_e(nsh,ncls), den_o(nsh,ncls))
        ! unregularised even/odd sub-averages -> FRC and noise variance -> prior signal power
        call flex_cls_half_reproducibility(y, c, w, wq, fitmask, weights, repro, avgs_e, avgs_o, shell=shell, &
            &den_even=den_e, den_odd=den_o)
        call flex_cls_signal_power(avgs_e, avgs_o, den_e, den_o, shell, nsh, tau2)
        call flex_cls_restore_states(y, c, w, weights, avgs, shell=shell, tau2=tau2)
        call flex_cls_half_reproducibility(y, c, w, wq, fitmask, weights, repro, avgs_e, avgs_o, shell=shell, tau2=tau2)
        call memoize_ft_maps(ldim, params%smpd)
        call planes_to_images(avgs,   cavgs)
        call planes_to_images(avgs_e, cavgs_even)
        call planes_to_images(avgs_o, cavgs_odd)
        call forget_ft_maps
        sep = 0.
        if( present(separation) ) sep = separation
        do s = 1, ncls
            nf = real(count(labels == s))
            if( present(neff) ) nf = neff(s)
            write(logfhandle,'(A,I8,A,I4,A,I8,A,F10.1,A,F7.3,A,F7.3)') 'Cls expansion flex subclass: class=', cls_id, ' sub=', s, &
                &' pop=', count(labels == s), ' neff=', nf, ' sep=', sep, ' repro=', repro(s)
        end do
        call flush(logfhandle)
        deallocate(avgs, avgs_e, avgs_o, tau2, den_e, den_o)

      contains

        subroutine planes_to_images( planes, imgs )
            complex(sp),              intent(in)    :: planes(:,:)
            type(image), allocatable, intent(inout) :: imgs(:)
            integer :: j, s
            call alloc_imgarr(ncls, ldim, params%smpd, imgs)
            do s = 1, ncls
                call imgs(s)%zero_and_flag_ft
                do j = 1, ncoeff
                    call imgs(s)%set_cmat_at(ft_map_phys_addrh(hidx(j),kidx(j)), ft_map_phys_addrk(hidx(j),kidx(j)), 1, planes(j,s))
                    ! the half-plane keeps the h=0 column for k>=0 only; the physical layout stores
                    ! both halves of that column, so the k<0 mate is set to the conjugate
                    if( hidx(j) == 0 .and. kidx(j) > 0 ) call imgs(s)%set_cmat_at(ft_map_phys_addrh(0,-kidx(j)), &
                        &ft_map_phys_addrk(0,-kidx(j)), 1, conjg(planes(j,s)))
                end do
                call imgs(s)%ifft
                call imgs(s)%mul(gridcorr_img)
            end do
        end subroutine planes_to_images

    end subroutine flex_class_restore_deliver

    !> the restoration pass: the commander's flag file in the execution directory, seen by the workers too
    logical function flex_restore_mode()
        flex_restore_mode = file_exists(FLEX_CLS_RESTORE_FLAG)
    end function flex_restore_mode

    !> one class into nsub subclasses: planes, cross-fitted fit with nuisance, placement, delivery; nsub = 1 keeps the
    !! class whole (one weighted average of all members)
    subroutine split_class_flex( params, build, spproj, cls_id, nsub, ncomp, kfit, l_sigma, pinds, labels, weights, cavgs, &
                                 separation, repro, cavgs_even, cavgs_odd )
        use simple_flex_cls_expansion,  only: flex_cls_model, flex_cls_fit_crossed, flex_cls_weighted_mean, flex_cls_pose_tangents, &
            &flex_cls_split_floor
        use simple_rnd,             only: seed_rnd_fixed
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id, nsub, ncomp, kfit
        logical,                  intent(in)    :: l_sigma
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        real,        allocatable, intent(out)   :: weights(:,:)
        type(image), allocatable, intent(inout) :: cavgs(:)
        real,                     intent(out)   :: separation
        real,        allocatable, intent(out)   :: repro(:)
        type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
        type(image)          :: gridcorr_img
        type(flex_cls_model) :: model
        complex(sp), allocatable :: y(:,:), nuis(:,:)
        complex(dp), allocatable :: mu_all(:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        logical,     allocatable :: fitmask(:)
        integer :: ldim(3), ncoeff, nptcls, j
        separation = 0.
        if( allocated(cavgs) ) call dealloc_imgarr(cavgs)
        call flex_class_planes(params, build, spproj, cls_id, kfit, l_sigma, pinds, y, c, w, wq, &
            &hidx, kidx, shell, fitmask, gridcorr_img, ldim, ncoeff, nptcls)
        if( nptcls < 1 ) THROW_HARD('cls_expansion: a class has no members')
        if( nsub == 1 )then
            allocate(weights(nptcls,1), source=1.)
            allocate(labels(nptcls), source=1)
            call flex_class_restore_deliver(params, cls_id, 1, y, c, w, wq, fitmask, shell, maxval(shell), hidx, kidx, &
                &gridcorr_img, ldim, ncoeff, nptcls, labels, weights, cavgs, repro, cavgs_even, cavgs_odd)
        else
            if( nptcls < flex_cls_split_floor(nsub, ncomp) ) THROW_HARD('cls_expansion: a class is too small for its subclass count')
            call seed_rnd_fixed(FLEX_CLS_SEED_BASE + cls_id)
            ! nuisance columns: the pose tangents and the class mean (per-particle contrast)
            allocate(mu_all(ncoeff), nuis(ncoeff,4))
            call flex_cls_weighted_mean(y, c, w, mu_all)
            call flex_cls_pose_tangents(mu_all, hidx, kidx, params%box, nuis(:,1:3))
            nuis(:,4) = cmplx(mu_all, kind=sp)
            call flex_cls_fit_crossed(model, y, c, w, wq, fitmask, ncomp, verbose=.true., nuis=nuis)
            write(logfhandle,'(A,4F9.4)') 'Cls expansion flex nuisance rms (x, y, rot, contrast; latent units): ', &
                &[(sqrt(sum(model%nu(:,j)**2)/real(nptcls,dp)), j=1,4)]
            deallocate(mu_all, nuis)
            call flex_class_deliver(params, cls_id, model, nsub, y, c, w, wq, fitmask, shell, maxval(shell), hidx, kidx, &
                &gridcorr_img, ldim, ncoeff, nptcls, labels, weights, cavgs, separation, repro, cavgs_even, cavgs_odd)
            call model%kill
        endif
        call gridcorr_img%kill
        deallocate(y, c, w, wq, hidx, kidx, shell, fitmask)
    end subroutine split_class_flex

    !> one class restored on the current alignment from the split's soft weights; the hard labels are the split's
    subroutine restore_class_flex( params, build, spproj, cls_id, nsub, first_sub, separation, kfit, l_sigma, sw, pinds, &
                                   labels, weights, cavgs, repro, cavgs_even, cavgs_odd )
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id, nsub, first_sub, kfit
        real,                     intent(in)    :: separation   !< the split's, for the log
        logical,                  intent(in)    :: l_sigma
        type(split_weights),      intent(in)    :: sw
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        real,        allocatable, intent(out)   :: weights(:,:)
        type(image), allocatable, intent(inout) :: cavgs(:)
        real,        allocatable, intent(out)   :: repro(:)
        type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
        type(image) :: gridcorr_img
        complex(sp), allocatable :: y(:,:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        logical,     allocatable :: fitmask(:)
        integer :: ldim(3), ncoeff, nptcls, k, ip, o
        if( allocated(cavgs) ) call dealloc_imgarr(cavgs)
        call flex_class_planes(params, build, spproj, cls_id, kfit, l_sigma, pinds, y, c, w, wq, &
            &hidx, kidx, shell, fitmask, gridcorr_img, ldim, ncoeff, nptcls)
        if( nptcls < 1 ) THROW_HARD('cls_expansion restoration: a class has no members')
        allocate(weights(nptcls,nsub), source=0.)
        allocate(labels(nptcls), source=0)
        do k = 1, nptcls
            ip = pinds(k)
            o  = sw%off(ip)
            if( o < 1 .or. sw%parent(ip) /= cls_id .or. sw%nsub(ip) /= nsub ) THROW_HARD('cls_expansion restoration: a member has no split weights for its class')
            weights(k,:) = sw%w(o:o+nsub-1)
            labels(k)    = sw%sub(ip) - first_sub + 1
            if( labels(k) < 1 .or. labels(k) > nsub ) THROW_HARD('cls_expansion restoration: a hard label outside its class')
        end do
        call flex_class_restore_deliver(params, cls_id, nsub, y, c, w, wq, fitmask, shell, maxval(shell), hidx, kidx, &
            &gridcorr_img, ldim, ncoeff, nptcls, labels, weights, cavgs, repro, cavgs_even, cavgs_odd, neff=sum(weights, dim=1), &
            &separation=separation)
        call gridcorr_img%kill
        deallocate(y, c, w, wq, hidx, kidx, shell, fitmask)
    end subroutine restore_class_flex

















    !> the split into the project: hard subclass and parent per particle, one cls2D/cls3D row per subclass with its
    !! readouts, and the sub-class average stack
    subroutine apply_split_project_updates(spproj, params, nsplit, new_class, new_parent, parent_of_subcls, pop_of_subcls, &
                                           neff_of_subcls, sep_of_subcls, repro_of_subcls, cavgs_stk)
        type(sp_project), intent(inout) :: spproj
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nsplit
        integer,          intent(in)    :: new_class(:), new_parent(:)
        integer,          intent(in)    :: parent_of_subcls(:), pop_of_subcls(:)
        real,             intent(in)    :: neff_of_subcls(:)  !< soft population: the sum of the subclass's weights
        real,             intent(in)    :: sep_of_subcls(:)   !< Ashman's D of the parent's split (same for its subclasses)
        real,             intent(in)    :: repro_of_subcls(:) !< cross-half reproducibility of the subclass difference
        type(string),     intent(in)    :: cavgs_stk          !< the sub-class average stack to register
        integer :: i
        call spproj%os_ptcl2D%set_all2single('class',   0)
        call spproj%os_ptcl2D%set_all2single('cluster', 0)
        do i = 1, size(new_class)
            if( new_class(i) <= 0 ) cycle
            call spproj%os_ptcl2D%set(i, 'class',   new_class(i))
            call spproj%os_ptcl2D%set(i, 'cluster', new_parent(i))
        end do
        call spproj%os_cls2D%new(nsplit, is_ptcl=.false.)
        call spproj%os_cls3D%new(nsplit, is_ptcl=.false.)
        do i = 1, nsplit
            call spproj%os_cls2D%set(i, 'cluster', parent_of_subcls(i))
            call spproj%os_cls2D%set(i, 'pop',     pop_of_subcls(i))
            call spproj%os_cls2D%set(i, 'accept',  1)
            call spproj%os_cls2D%set(i, 'state',   1)
            call spproj%os_cls2D%set(i, 'neff',    neff_of_subcls(i))
            call spproj%os_cls2D%set(i, 'sep',     sep_of_subcls(i))
            call spproj%os_cls2D%set(i, 'repro',   repro_of_subcls(i))
            call spproj%os_cls3D%set(i, 'cluster', parent_of_subcls(i))
            call spproj%os_cls3D%set(i, 'accept',  1)
            call spproj%os_cls3D%set(i, 'state',   1)
        end do
        call spproj%add_cavgs2os_out(cavgs_stk, params%smpd, imgkind='cavg')
        call spproj%write(params%projfile)
    end subroutine apply_split_project_updates

    !> the part outputs in global subclass order (by parent, then local subclass): class map, assignments into the
    !! project, weights table and the average stacks
    subroutine merge_worker_outputs(params, nparts_run)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        type(sp_project) :: spproj
        integer, allocatable :: part_counts(:), part_parent(:,:), part_local(:,:), part_pop(:,:), part_global(:,:)
        real,    allocatable :: part_neff(:,:), part_sep(:,:), part_repro(:,:)
        integer, allocatable :: comb_part(:), comb_row(:), comb_parent(:), comb_local(:), comb_pop(:)
        integer, allocatable :: new_class(:), new_parent(:), first_glob(:), nsub_of(:)
        real,    allocatable :: neff(:), sep(:), repro(:)
        type(string) :: map_fname, assign_fname
        integer :: ipart, nlocal, total, max_count, idx, i, funit, ios, pind, parent_cls, local_cls, p
        character(len=XLONGSTRLEN) :: line
        call spproj%read(params%projfile)
        write(logfhandle,'(A,I8)') 'Cls expansion merge: start nparts=', nparts_run
        call flush(logfhandle)
        allocate(part_counts(nparts_run), source=0)
        total = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_expansion_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call count_data_lines(map_fname, nlocal)
            part_counts(ipart) = nlocal
            total = total + nlocal
            call map_fname%kill
        end do
        if( total < 1 ) THROW_HARD('No subclass outputs produced by distributed cls_expansion workers')
        max_count = maxval(part_counts)
        allocate(part_parent(max_count, nparts_run), part_local(max_count, nparts_run), part_pop(max_count, nparts_run), &
                 part_global(max_count, nparts_run), source=0)
        allocate(part_neff(max_count, nparts_run), part_sep(max_count, nparts_run), part_repro(max_count, nparts_run), source=0.)
        allocate(comb_part(total), comb_row(total), comb_parent(total), comb_local(total), comb_pop(total), source=0)
        idx = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_expansion_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call read_part_map(map_fname, part_counts(ipart), part_parent(:,ipart), part_local(:,ipart), part_pop(:,ipart), &
                &part_neff(:,ipart), part_sep(:,ipart), part_repro(:,ipart))
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
        ! a parent's subclasses are consecutive: global = first_glob(parent) + local - 1
        allocate(first_glob(maxval(comb_parent)), nsub_of(maxval(comb_parent)), source=0)
        allocate(neff(total), sep(total), repro(total), source=0.)
        do idx = 1, total
            if( idx > 1 )then
                if( comb_parent(idx) == comb_parent(idx-1) .and. comb_local(idx) == comb_local(idx-1) )then
                    THROW_HARD('Duplicate parent/local subclass pair detected while merging cls_expansion outputs')
                endif
            endif
            p = comb_parent(idx)
            if( comb_local(idx) == 1 ) first_glob(p) = idx
            nsub_of(p) = nsub_of(p) + 1
            if( first_glob(p) < 1 .or. idx /= first_glob(p) + comb_local(idx) - 1 )then
                THROW_HARD('the subclasses of a parent are not consecutive while merging cls_expansion outputs')
            endif
            part_global(comb_row(idx), comb_part(idx)) = idx
            neff(idx)  = part_neff(comb_row(idx), comb_part(idx))
            sep(idx)   = part_sep(comb_row(idx), comb_part(idx))
            repro(idx) = part_repro(comb_row(idx), comb_part(idx))
        end do
        map_fname = string(FLEX_CLS_MAP_FILE)
        open(newunit=funit, file=map_fname%to_char(), status='replace', action='write')
        write(funit,'(A)') '# subclass parent_class local_subclass pop neff separation repro'
        do idx = 1, total
            write(funit,'(I8,1X,I8,1X,I8,1X,I8,1X,F12.3,1X,F8.3,1X,F8.3)') idx, comb_parent(idx), comb_local(idx), &
                &comb_pop(idx), neff(idx), sep(idx), repro(idx)
        end do
        close(funit)
        call map_fname%kill
        write(logfhandle,'(A,I8)') 'Cls expansion merge: subclass map merged total=', total
        call flush(logfhandle)
        allocate(new_class(spproj%os_ptcl2D%get_noris()), new_parent(spproj%os_ptcl2D%get_noris()), source=0)
        do ipart = 1, nparts_run
            assign_fname = string('cls_expansion_assignments_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            open(newunit=funit, file=assign_fname%to_char(), status='old', action='read', iostat=ios)
            call fileiochk('merge_worker_outputs opening '//assign_fname%to_char(), ios)
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*) pind, parent_cls, local_cls
                new_class(pind)  = global_subclass(parent_cls, local_cls)
                new_parent(pind) = parent_cls
            end do
            close(funit)
            call assign_fname%kill
        end do
        call merge_part_weights_and_stacks(params, nparts_run, part_counts, part_global, first_glob, nsub_of)
        call apply_split_project_updates(spproj, params, total, new_class, new_parent, comb_parent, comb_pop, &
            &neff, sep, repro, string(FLEX_CLS_CAVGS_FILE))
        call spproj%kill
        deallocate(new_class, new_parent, part_counts, part_parent, part_local, part_pop, part_global, part_neff, part_sep, &
                   part_repro, comb_part, comb_row, comb_parent, comb_local, comb_pop, first_glob, nsub_of, neff, sep, repro)

    contains

        integer function global_subclass( parent, loc )
            integer, intent(in) :: parent, loc
            global_subclass = 0
            if( parent >= 1 .and. parent <= size(first_glob) )then
                if( loc >= 1 .and. loc <= nsub_of(parent) ) global_subclass = first_glob(parent) + loc - 1
            endif
            if( global_subclass < 1 ) THROW_HARD('Could not map worker-local subclass to global subclass during merge')
        end function global_subclass

    end subroutine merge_worker_outputs

    !> the part weight tables and average stacks in global subclass order
    subroutine merge_part_weights_and_stacks(params, nparts_run, part_counts, part_global, first_glob, nsub_of)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run, part_counts(:), part_global(:,:), first_glob(:), nsub_of(:)
        type(string) :: wts_fname, stk_fname, out_wts, out_stk
        type(image)  :: img
        character(len=:), allocatable :: line
        real,             allocatable :: w(:)
        integer :: ipart, funit_in, funit_out, ios, pind, parent_cls, local_cls, n, irow, ldim(3), nimgs, j, ieo
        allocate(w(max(1, maxval(nsub_of))))
        out_wts = string(FLEX_CLS_WEIGHTS_FILE)
        open(newunit=funit_out, file=out_wts%to_char(), status='replace', action='write')
        write(funit_out,'(A)') '# particle_index parent_class subclass nsub weights(1:nsub)'
        do ipart = 1, nparts_run
            wts_fname = flex_weights_part_fname(ipart, params%numlen)
            open(newunit=funit_in, file=wts_fname%to_char(), status='old', action='read', iostat=ios)
            call fileiochk('merge_part_weights_and_stacks opening '//wts_fname%to_char(), ios)
            do
                ! rows are as wide as their parent has subclasses
                call read_record(funit_in, line, ios)
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*,iostat=ios) pind, parent_cls, local_cls, n
                if( ios /= 0 ) THROW_HARD('read error in a cls_expansion part weight table')
                if( n < 1 .or. n > size(w) ) THROW_HARD('inconsistent row in a cls_expansion part weight table')
                read(line,*,iostat=ios) pind, parent_cls, local_cls, n, (w(j), j=1,n)
                if( ios /= 0 ) THROW_HARD('short row in a cls_expansion part weight table')
                if( parent_cls < 1 .or. parent_cls > size(nsub_of) ) THROW_HARD('unknown parent in a cls_expansion part weight table')
                if( n /= nsub_of(parent_cls) .or. local_cls < 1 .or. local_cls > n ) THROW_HARD('inconsistent row in a cls_expansion part weight table')
                write(funit_out,'(I10,1X,I10,1X,I10,1X,I6,*(1X,F8.6))') pind, parent_cls, first_glob(parent_cls) + local_cls - 1, &
                    &n, w(1:n)
            end do
            close(funit_in)
            call wts_fname%kill
        end do
        close(funit_out)
        ! the stacks (full, even, odd): part row r -> global index part_global(r, ipart)
        do ieo = 0, 2
            select case(ieo)
                case(0); out_stk = string(FLEX_CLS_CAVGS_FILE)
                case(1); out_stk = string(FLEX_CLS_CAVGS_EVEN)
                case(2); out_stk = string(FLEX_CLS_CAVGS_ODD)
            end select
            do ipart = 1, nparts_run
                select case(ieo)
                    case(0); stk_fname = flex_cavgs_part_fname(ipart, params%numlen)
                    case(1); stk_fname = flex_cavgs_part_fname(ipart, params%numlen, 'even')
                    case(2); stk_fname = flex_cavgs_part_fname(ipart, params%numlen, 'odd')
                end select
                call find_ldim_nptcls(stk_fname, ldim, nimgs)
                if( nimgs /= part_counts(ipart) ) THROW_HARD('flex part average stack does not match its class map: '//stk_fname%to_char())
                ldim(3) = 1
                call img%new(ldim, params%smpd)
                do irow = 1, nimgs
                    call img%read(stk_fname, irow)
                    call img%write(out_stk, part_global(irow, ipart))
                end do
                call img%kill
                call stk_fname%kill
            end do
        end do
        deallocate(w)
        call out_wts%kill
        call out_stk%kill
    end subroutine merge_part_weights_and_stacks

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

    subroutine read_part_map(fname, nrows, parents, locals, pops, neffs, seps, repros)
        type(string), intent(in)  :: fname
        integer,      intent(in)  :: nrows
        integer,      intent(out) :: parents(:), locals(:), pops(:)
        real,         intent(out) :: neffs(:), seps(:), repros(:)
        integer :: funit, ios, irow, row
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
            read(line,*) row, parents(irow), locals(irow), pops(irow), neffs(irow), seps(irow), repros(irow)
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

end module simple_cls_expansion_strategy
