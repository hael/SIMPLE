!@descr: cls_expansion execution strategies: shared memory, distributed master and distributed worker
! ncls is the total number of subclasses: the master (or the shared-memory run) groups the classes by
! class-average correlation and budgets the subclasses per class. Every class is split once into its subclasses
! (hard labels) and delivers their conditional averages with independent even/odd halves and per-subclass FRCs.
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
public :: FLEX_CLS_CAVGS_FILE
public :: cls_expansion_shmem_strategy
public :: cls_expansion_worker_strategy
public :: cls_expansion_master_strategy
public :: create_cls_expansion_strategy

private
#include "simple_local_flags.inc"

integer, parameter :: FLEX_CLS_NEIGS_DEFAULT  = 16   !< rank of the per-class covariance model
integer, parameter :: FLEX_CLS_SEED_BASE      = 20260930 !< per-class random seed = base + class index
real,    parameter :: FLEX_CLS_BAND_FIRST_LP  = 20.  !< the band-selection fit stops here at most (A); lp= sets the band
character(len=*), parameter :: FLEX_CLS_MAP_FILE     = 'cls_expansion_class_map.txt'  !< subclass parent local pop neff sep repro
character(len=*), parameter :: FLEX_CLS_ALLOC_FILE   = 'cls_expansion_allocation.txt' !< class group pop nsub: the split's budget
character(len=*), parameter :: FLEX_CLS_GROUPS_FBODY = 'cls_expansion_groups'         !< class-average groups: stacks and table
character(len=*), parameter :: FLEX_CLS_CAVGS_FILE   = 'cls_expansion_cavgs.mrc'      !< conditional sub-class averages, filtered
character(len=*), parameter :: FLEX_CLS_CAVGS_EVEN   = 'cls_expansion_cavgs_even.mrc' !< the same from the even half model, unfiltered
character(len=*), parameter :: FLEX_CLS_CAVGS_ODD    = 'cls_expansion_cavgs_odd.mrc'  !< the same from the odd half model, unfiltered
character(len=*), parameter :: FLEX_CLS_FRCS_TXT     = 'cls_expansion_frcs.txt'       !< even/odd FRC per subclass (row frc(1:nyq))
character(len=*), parameter :: FLEX_CLS_FRCS         = 'cls_expansion_frcs.bin'       !< the same as class FRCs (frc2D)

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
        integer, allocatable :: cls_inds(:), cls_pops(:), nsub(:)
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
        ! the subclass budget is set once, here, and read by the workers
        call spproj%read(params%projfile)
        call allocate_subclasses(params, spproj, cls_inds, cls_pops, self%nthr_master, nsub)
        call spproj%kill
        deallocate(nsub)
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
        label = 'class'
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

    !> the classes of this run (all, or a worker's part): each split under the budget into its subclasses; writes the
    !! class map, the per-subclass FRCs, the stacks and, outside a worker, the project (a worker writes its assignments)
    subroutine run_local_split(params, build, cline, spproj, part, l_write_project)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: part
        logical,          intent(in)    :: l_write_project
        type(string)             :: map_fname, assign_fname, frc_fname, stk_fname, stk_even, stk_odd
        type(image), allocatable :: cavgs(:), cavgs_even(:), cavgs_odd(:)
        integer,     allocatable :: cls_inds(:), cls_pops(:), nsub(:), pinds(:), labels(:), half_pops(:,:)
        integer,     allocatable :: new_class(:), new_parent(:), parent_of_subcls(:), pop_of_subcls(:)
        integer,     allocatable :: pop_even_of_subcls(:), pop_odd_of_subcls(:)
        real,        allocatable :: repro(:), frcs(:,:), neff_of_subcls(:), sep_of_subcls(:), repro_of_subcls(:)
        real    :: separation
        integer :: i, iglob, ntot_sub, funit_map, funit_assign, funit_frc, kfit
        logical :: l_sigma
        call collect_split_classes(cline, params, build, cls_inds, cls_pops)
        if( l_write_project )then
            call allocate_subclasses(params, spproj, cls_inds, cls_pops, params%nthr, nsub)
        else
            call read_allocation(cls_inds, nsub)
        endif
        ntot_sub = sum(nsub)
        call prepare_flex_sigma2(params, build, spproj, cline, l_sigma)
        kfit = flex_fit_band(params, cline)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,L1)') 'Cls expansion flex: classes=', size(cls_inds), &
            &' subclasses=', ntot_sub, ' neigs=', params%neigs, ' fit band kto=', kfit, ' (0: data-derived per class) sigma2 loaded=', l_sigma
        call flush(logfhandle)
        if( l_write_project )then
            map_fname = string(FLEX_CLS_MAP_FILE)
            frc_fname = string(FLEX_CLS_FRCS_TXT)
            stk_fname = string(FLEX_CLS_CAVGS_FILE)
            stk_even  = string(FLEX_CLS_CAVGS_EVEN)
            stk_odd   = string(FLEX_CLS_CAVGS_ODD)
            allocate(new_class(spproj%os_ptcl2D%get_noris()), new_parent(spproj%os_ptcl2D%get_noris()), source=0)
            allocate(parent_of_subcls(ntot_sub), pop_of_subcls(ntot_sub), source=0)
            allocate(pop_even_of_subcls(ntot_sub), pop_odd_of_subcls(ntot_sub), source=0)
            allocate(neff_of_subcls(ntot_sub), sep_of_subcls(ntot_sub), repro_of_subcls(ntot_sub), source=0.)
        else
            map_fname    = string('cls_expansion_class_map_part')//int2str_pad(part, params%numlen)//TXT_EXT
            assign_fname = string('cls_expansion_assignments_part')//int2str_pad(part, params%numlen)//TXT_EXT
            frc_fname    = flex_frcs_part_fname(part, params%numlen)
            stk_fname    = flex_cavgs_part_fname(part, params%numlen)
            stk_even     = flex_cavgs_part_fname(part, params%numlen, 'even')
            stk_odd      = flex_cavgs_part_fname(part, params%numlen, 'odd')
            open(newunit=funit_assign, file=assign_fname%to_char(), status='replace', action='write')
            write(funit_assign,'(A)') '# particle_index parent_class local_subclass'
        endif
        ! in a part the subclass column is the part's row (map, FRCs)
        open(newunit=funit_map, file=map_fname%to_char(), status='replace', action='write')
        write(funit_map,'(A)') '# subclass parent_class local_subclass pop neff separation repro pop_even pop_odd'
        open(newunit=funit_frc, file=frc_fname%to_char(), status='replace', action='write')
        if( l_write_project )then
            write(funit_frc,'(A)') '# row frc(1:nyq): even/odd FRC of the conditional average (row: global)'
        else
            write(funit_frc,'(A)') '# row frc(1:nyq): even/odd FRC of the conditional average (row: the part''s)'
        endif
        iglob = 0
        do i = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8,A,I6)') 'Cls expansion starting: class=', cls_inds(i), ' nptcls=', cls_pops(i), &
                &' nsub=', nsub(i)
            call flush(logfhandle)
            call split_class_flex(params, build, spproj, cls_inds(i), nsub(i), params%neigs, kfit, l_sigma, &
                &pinds, labels, cavgs, separation, repro, cavgs_even, cavgs_odd, frcs, half_pops)
            call record_class(cls_inds(i), nsub(i), pinds, labels, cavgs, separation, repro, cavgs_even, cavgs_odd)
            if( allocated(pinds)      ) deallocate(pinds)
            if( allocated(half_pops)  ) deallocate(half_pops)
            if( allocated(repro)      ) deallocate(repro)
            if( allocated(labels)     ) deallocate(labels)
            if( allocated(frcs)       ) deallocate(frcs)
            if( allocated(cavgs)      ) call dealloc_imgarr(cavgs)
            if( allocated(cavgs_even) ) call dealloc_imgarr(cavgs_even)
            if( allocated(cavgs_odd)  ) call dealloc_imgarr(cavgs_odd)
        end do
        close(funit_map)
        close(funit_frc)
        if( .not. l_write_project ) close(funit_assign)
        if( l_write_project )then
            if( iglob /= ntot_sub ) THROW_HARD('cls_expansion: the delivered subclasses differ from the budget')
            call flex_frcs_txt2bin(params, [frc_fname], iglob)
            call apply_split_project_updates(spproj, params, iglob, new_class, new_parent, parent_of_subcls, pop_of_subcls, &
                &neff_of_subcls, sep_of_subcls, repro_of_subcls, pop_even_of_subcls, pop_odd_of_subcls, stk_fname)
            deallocate(new_class, new_parent, parent_of_subcls, pop_of_subcls, neff_of_subcls, sep_of_subcls, repro_of_subcls, &
                &pop_even_of_subcls, pop_odd_of_subcls)
        endif
        call map_fname%kill
        call frc_fname%kill
        call stk_fname%kill
        call stk_even%kill
        call stk_odd%kill
        if( assign_fname%is_allocated() ) call assign_fname%kill

    contains

        !> one class into the map, FRCs, stacks and assignments (part) or project arrays (whole run); neff is the hard
        !! population
        subroutine record_class( cls_id, nsplit, pinds, labels, cavgs, separation, repro, cavgs_even, cavgs_odd )
            integer,                  intent(in)    :: cls_id, nsplit
            integer,                  intent(in)    :: pinds(:), labels(:)
            type(image), allocatable, intent(inout) :: cavgs(:)
            real,                     intent(in)    :: separation
            real,        allocatable, intent(in)    :: repro(:)
            type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
            integer :: j, k
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Cls expansion summary: class=', cls_id, ' nptcls=', size(labels), ' nsubcls=', nsplit
            call flush(logfhandle)
            do j = 1, nsplit
                iglob = iglob + 1
                if( l_write_project )then
                    parent_of_subcls(iglob) = cls_id
                    pop_of_subcls(iglob)    = count(labels == j)
                    neff_of_subcls(iglob)   = real(count(labels == j))
                    sep_of_subcls(iglob)    = separation
                    repro_of_subcls(iglob)  = repro(j)
                    pop_even_of_subcls(iglob) = half_pops(1,j)
                    pop_odd_of_subcls(iglob)  = half_pops(2,j)
                endif
                write(funit_map,'(I8,1X,I8,1X,I8,1X,I8,1X,F12.3,1X,F8.3,1X,F8.3,2(1X,I8))') iglob, cls_id, j, &
                    &count(labels == j), real(count(labels == j)), separation, repro(j), half_pops(1,j), half_pops(2,j)
                write(funit_frc,'(I8,*(1X,F8.5))') iglob, frcs(:,j)
                ! the sub-class average lands at its row: global index here, local row in a part
                call cavgs(j)%write(stk_fname, iglob)
                call cavgs_even(j)%write(stk_even, iglob)
                call cavgs_odd(j)%write(stk_odd, iglob)
            end do
            do k = 1, size(labels)
                if( labels(k) <= 0 ) cycle
                if( l_write_project )then
                    new_class(pinds(k))  = iglob - nsplit + labels(k)
                    new_parent(pinds(k)) = cls_id
                else
                    write(funit_assign,'(I10,1X,I10,1X,I10)') pinds(k), cls_id, labels(k)
                endif
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

    ! ===== the flex model

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

    !> the global fit band (shell) from lp=; 0 without lp= (the bands are then data-derived per class, split_class_flex)
    integer function flex_fit_band( params, cline ) result( kfit )
        type(parameters), intent(in) :: params
        class(cmdline),   intent(in) :: cline
        real    :: lp_fit
        integer :: nyq
        kfit = 0
        if( .not. cline%defined('lp') ) return
        nyq    = fdim(params%box) - 1
        lp_fit = max(params%lp, 2. * params%smpd)
        kfit   = max(2, min(nyq, calc_fourier_index(lp_fit, params%box, params%smpd)))
    end function flex_fit_band

    !> the two canonical sigma2 spectra (one per even/odd half)
    subroutine prepare_flex_sigma2( params, build, spproj, cline, loaded )
        use simple_sigma2_files, only: canonical_sigma2_consumable, load_sigma2
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        class(cmdline),   intent(inout) :: cline
        logical,          intent(out)   :: loaded
        character(len=STDLEN) :: message
        loaded = .false.
        if( params%cc_objfun /= OBJFUN_EUCLID )then
            write(logfhandle,'(A)') 'Cls expansion flex: objfun=cc, the covariance fit uses unit noise weights'
            return
        endif
        if( .not. canonical_sigma2_consumable(spproj, spproj%os_ptcl2D, params%box, params%smpd, message) )then
            write(logfhandle,'(A)') 'Cls expansion flex WARNING: no consumable canonical sigma2 state ('//trim(message)// &
                &'); the covariance fit uses the class residual spectrum as noise weights. Run calc_pspec &
                &(or a euclid 2D/3D pass) on the project first.'
            return
        endif
        call load_sigma2(params, build%pftc, build%esig, spproj, spproj%os_ptcl2D, loaded)
    end subroutine prepare_flex_sigma2

    !> the members of one class as Fourier planes in the class frame (inscribed-circle mask) with their CTF, half-plane
    !! lattice, shells 1..nyq; noise weights from the canonical sigma2 or the class's own residual spectrum
    subroutine flex_class_planes( params, build, spproj, cls_id, kfit, l_sigma, pinds, y, c, w, wq, &
                                  hidx, kidx, shell, fitmask, gridcorr_img, ldim, ncoeff, nptcls )
        use simple_flex_cls_expansion,  only: flex_cls_shell_noise
        use simple_euclid_sigma2,       only: sigma2_half
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
        integer :: lims(3,2), nyq, i, j, h, k, sh, iptcl, ih
        logical :: l_ctf, l_flip
        integer, allocatable :: members(:)
        nptcls = 0
        ncoeff = 0
        if( allocated(pinds) ) deallocate(pinds)
        call spproj%os_ptcl2D%get_pinds(cls_id, 'class', members)
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
        ! the box's inscribed circle, soft edge: drops the box corners the in-plane rotation moves, which would
        ! otherwise become directions of variation, and keeps everything beyond mskdiam
        call imgs(1)%memoize_mask_coords
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nptcls
            call imgs(i)%ifft
            call imgs(i)%mask2D_soft(real(params%box) / 2. - 2., backgr=0.)
            call imgs(i)%fft
        end do
        !$omp end parallel do
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
        !$omp private(i,j,iptcl,ih,e3,mat,loc,sfsq,ang,df,tval,ctfparms,tfun,ctfvals,sum_df,diff_df,angast,phshift,amp_contr,wl,&
        !$omp half_wl2_cs)
        do i = 1, nptcls
            iptcl = pinds(i)
            e3    = spproj%os_ptcl2D%e3get(iptcl)
            ! the column of the member's even/odd half in the two global spectra
            ih    = 1
            if( l_sigma ) ih = sigma2_half(spproj%os_ptcl2D%get_eo(iptcl))
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
                    w(j,i) = 1. / max(build%esig%sigma2_noise(shell(j), ih), TINY)
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
            !$omp parallel do default(shared) private(i,j) schedule(static) proc_bind(close)
            do i = 1, nptcls
                do j = 1, ncoeff
                    w(j,i) = 1. / s2(shell(j))
                end do
            end do
            !$omp end parallel do
            write(logfhandle,'(A,I8,A,ES10.3,A,ES10.3)') 'Cls expansion flex noise from residuals: class=', cls_id, &
                &' sigma2 shell 1=', s2(1), ' shell nyq/4=', s2(max(1,nyq/4))
            deallocate(s2)
        endif
    end subroutine flex_class_planes

    !> one fitted class: placement (hard labels) and the conditional averages of its subclasses
    subroutine flex_class_deliver( params, cls_id, model, ncls, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, &
                                   gridcorr_img, ldim, nptcls, labels, cavgs, cavgs_even, cavgs_odd, frcs, separation, repro )
        use simple_flex_cls_expansion,  only: flex_cls_model, flex_cls_place_states
        type(parameters),         intent(in)    :: params
        integer,                  intent(in)    :: cls_id, ncls, ncomp, ldim(3), nptcls
        type(flex_cls_model),     intent(in)    :: model
        complex(sp),              intent(in)    :: y(:,:)
        real(sp),                 intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,                  intent(in)    :: fitmask(:), even(:)
        integer,                  intent(in)    :: shell(:), hidx(:), kidx(:)
        type(image),              intent(inout) :: gridcorr_img
        integer,     allocatable, intent(out)   :: labels(:)
        type(image), allocatable, intent(inout) :: cavgs(:), cavgs_even(:), cavgs_odd(:)
        real,        allocatable, intent(inout) :: frcs(:,:)
        real,                     intent(out)   :: separation
        real,        allocatable, intent(out)   :: repro(:)
        integer(timer_int_kind) :: t0
        integer :: s
        t0 = tic()
        allocate(labels(nptcls))
        call flex_cls_place_states(model, ncls, labels, separation)
        call flex_time(cls_id, 'placement', t0)
        call flex_class_conditional_averages(params, cls_id, ncls, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, labels, &
            &gridcorr_img, ldim, cavgs, cavgs_even, cavgs_odd, frcs, repro, model=model)
        do s = 1, ncls
            write(logfhandle,'(A,I8,A,I4,A,I8,A,F7.3,A,F7.3)') 'Cls expansion flex subclass: class=', cls_id, ' sub=', s, &
                &' pop=', count(labels == s), ' sep=', separation, ' repro=', repro(s)
        end do
        call flush(logfhandle)
    end subroutine flex_class_deliver

    !> the subclass averages P + D_s of one class from the joint model, the unfiltered even/odd ones from the half models
    !! (too few members: half parents + joint deviations), their FRCs, and the full one filtered by 2F/(1+F) of its FRC.
    !! Without a model (a class kept whole): the CTF-corrected parent and half parents
    subroutine flex_class_conditional_averages( params, cls_id, ncls, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, labels, &
                                  gridcorr_img, ldim, cavgs, cavgs_even, cavgs_odd, frcs, repro, model )
        use simple_flex_cls_expansion, only: flex_cls_model, flex_cls_ctf_mean, flex_cls_extend_basis, &
                                             flex_cls_conditional_averages, flex_cls_shell_frc
        use simple_estimate_ssnr,      only: fsc2optlp_sub
        use simple_memoize_ft_maps,    only: memoize_ft_maps, forget_ft_maps, ft_map_phys_addrh, ft_map_phys_addrk
        type(parameters),         intent(in)    :: params
        integer,                  intent(in)    :: cls_id, ncls, ncomp, ldim(3)
        complex(sp),              intent(in)    :: y(:,:)
        real(sp),                 intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,                  intent(in)    :: fitmask(:), even(:)
        integer,                  intent(in)    :: shell(:), hidx(:), kidx(:), labels(:)
        type(image),              intent(inout) :: gridcorr_img
        type(image), allocatable, intent(inout) :: cavgs(:), cavgs_even(:), cavgs_odd(:)
        real,        allocatable, intent(inout) :: frcs(:,:)
        real,        allocatable, intent(out)   :: repro(:)
        type(flex_cls_model), optional, intent(in) :: model
        complex(sp), allocatable :: par(:,:), dev(:,:), parh(:,:), devh(:,:,:), planes(:,:)
        complex(dp), allocatable :: uext(:,:)
        real(sp),    allocatable :: fpar(:,:)
        real(dp),    allocatable :: num(:), pe(:), po(:)
        real    :: filt(maxval(shell)), fdev(maxval(shell))
        integer :: ncoeff, nsh, s, j, sh
        logical :: l_halves
        integer(timer_int_kind) :: t0
        t0 = tic()
        ncoeff = size(y,1); nsh = maxval(shell)
        l_halves = present(model) .and. min(count(even), count(.not. even)) >= max(2*ncls, 3*(ncomp+1))
        allocate(par(ncoeff,1), dev(ncoeff,ncls))
        dev = CMPLX_ZERO
        ! the joint model: its parent (CTF-corrected, no signal prior) and the subclass deviations on every coefficient
        call flex_cls_ctf_mean(y, c, w, par(:,1))
        if( present(model) )then
            allocate(uext(ncoeff,model%ncomp))
            call flex_cls_extend_basis(model, y, c, w, shell, uext, even=even)
            call flex_cls_conditional_averages(model, labels, uext, dev)
            deallocate(uext)
        endif
        call flex_time(cls_id, 'conditional', t0)
        ! the halves
        allocate(parh(ncoeff,2), devh(ncoeff,ncls,2))
        call flex_class_half_models(params, cls_id, ncls, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, labels, &
            &l_halves, parh, devh)
        if( .not. l_halves )then
            devh(:,:,1) = dev
            devh(:,:,2) = dev
        endif
        call flex_time(cls_id, 'halves', t0)
        ! the agreement of the half parents and of the half deviations (pooled over the subclasses), logged
        allocate(fpar(nsh,1))
        call flex_cls_shell_frc(parh(:,1:1), parh(:,2:2), shell, nsh, fpar)
        allocate(num(nsh), pe(nsh), po(nsh), source=0.d0)
        do s = 1, ncls
            do j = 1, ncoeff
                sh = shell(j)
                if( sh < 1 .or. sh > nsh ) cycle
                num(sh) = num(sh) + real(cmplx(devh(j,s,1),kind=dp) * conjg(cmplx(devh(j,s,2),kind=dp)), dp)
                pe(sh)  = pe(sh)  + abs(cmplx(devh(j,s,1),kind=dp))**2
                po(sh)  = po(sh)  + abs(cmplx(devh(j,s,2),kind=dp))**2
            end do
        end do
        fdev = 0.
        where( pe * po > 0.d0 ) fdev = real(num / sqrt(pe * po))
        ! the full (joint), the even and the odd conditional averages
        allocate(planes(ncoeff,3*ncls))
        do s = 1, ncls
            planes(:,s)        = par(:,1)  + dev(:,s)
            planes(:,ncls+s)   = parh(:,1) + devh(:,s,1)
            planes(:,2*ncls+s) = parh(:,2) + devh(:,s,2)
        end do
        if( allocated(frcs) ) deallocate(frcs)
        allocate(frcs(nsh,ncls), repro(ncls))
        call flex_cls_shell_frc(planes(:,ncls+1:2*ncls), planes(:,2*ncls+1:3*ncls), shell, nsh, frcs)
        ! the full conditional average filtered by its subclass's even/odd FRC
        do s = 1, ncls
            call fsc2optlp_sub(nsh, frcs(:,s), filt)
            do j = 1, ncoeff
                sh = shell(j)
                if( sh >= 1 .and. sh <= nsh ) planes(j,s) = planes(j,s) * filt(sh)
            end do
        end do
        do s = 1, ncls
            repro(s) = 0.
            if( sum(abs(devh(:,s,1))**2) * sum(abs(devh(:,s,2))**2) > 0. ) repro(s) = &
                &real(sum(devh(:,s,1) * conjg(devh(:,s,2)))) / sqrt(sum(abs(devh(:,s,1))**2) * sum(abs(devh(:,s,2))**2))
        end do
        if( present(model) )then
            write(logfhandle,'(A,I8,A,L1,A,F7.2,A,F7.2,A,F7.2,A)') 'Cls expansion conditional averages: class=', cls_id, ' halves=', l_halves, &
                &' FRC 0.143: parent', res143(fpar(:,1)), ' A, deviations', res143(fdev), ' A, subclasses (mean)', &
                &sum([(res143(frcs(:,s)), s=1,ncls)]) / real(ncls), ' A'
        else
            write(logfhandle,'(A,I8,A,F7.2,A)') 'Cls expansion class kept whole: class=', cls_id, ' FRC 0.143 of its half means', &
                &res143(fpar(:,1)), ' A'
        endif
        ! images
        call memoize_ft_maps(ldim, params%smpd)
        call to_images(planes(:,1:ncls),          cavgs)
        call to_images(planes(:,ncls+1:2*ncls),   cavgs_even)
        call to_images(planes(:,2*ncls+1:3*ncls), cavgs_odd)
        call forget_ft_maps
        call flex_time(cls_id, 'images', t0)
        deallocate(par, dev, parh, devh, fpar, num, pe, po, planes)

      contains

        !> the resolution (A) at which an FRC first falls below 0.143 (from shell 3 on)
        real function res143( frc )
            real, intent(in) :: frc(:)
            integer :: k
            res143 = real(params%box) * params%smpd / real(size(frc))
            do k = 3, size(frc)
                if( frc(k) < 0.143 )then
                    res143 = real(params%box) * params%smpd / real(k)
                    return
                endif
            end do
        end function res143

        subroutine to_images( pl, imgs )
            complex(sp),              intent(in)    :: pl(:,:)
            type(image), allocatable, intent(inout) :: imgs(:)
            integer :: jj, ss
            if( allocated(imgs) ) call dealloc_imgarr(imgs)
            call alloc_imgarr(size(pl,2), ldim, params%smpd, imgs)
            do ss = 1, size(pl,2)
                call imgs(ss)%zero_and_flag_ft
                do jj = 1, ncoeff
                    call imgs(ss)%set_cmat_at(ft_map_phys_addrh(hidx(jj),kidx(jj)), ft_map_phys_addrk(hidx(jj),kidx(jj)), 1, pl(jj,ss))
                    ! the half-plane keeps the h=0 column for k>=0 only; its k<0 mate is the conjugate
                    if( hidx(jj) == 0 .and. kidx(jj) > 0 ) call imgs(ss)%set_cmat_at(ft_map_phys_addrh(0,-kidx(jj)), &
                        &ft_map_phys_addrk(0,-kidx(jj)), 1, conjg(pl(jj,ss)))
                end do
                call imgs(ss)%ifft
                call imgs(ss)%mul(gridcorr_img)
            end do
        end subroutine to_images

    end subroutine flex_class_conditional_averages

    !> wall time of one step of one class, logged
    subroutine flex_time( cls_id, step, t0 )
        integer,                 intent(in)    :: cls_id
        character(len=*),        intent(in)    :: step
        integer(timer_int_kind), intent(inout) :: t0
        write(logfhandle,'(A,I8,A,A,A,F9.2)') 'Cls expansion time: class=', cls_id, ' step=', step, ' s=', toc(t0)
        call flush(logfhandle)
        t0 = tic()
    end subroutine flex_time

    !> the even and the odd half of one class: its CTF-corrected parent and, with l_models, its own model on the joint
    !! band (fitmask) and the subclass deviations U z_s of the joint labels; the halves share no particles
    subroutine flex_class_half_models( params, cls_id, ncls, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, labels, &
                                       l_models, parh, devh )
        use simple_flex_cls_expansion, only: flex_cls_model, flex_cls_fit, flex_cls_extend_basis, flex_cls_conditional_averages, &
                                             flex_cls_weighted_mean, flex_cls_pose_tangents, flex_cls_ctf_mean
        use simple_rnd,                only: seed_rnd_fixed
        type(parameters), intent(in)  :: params
        integer,          intent(in)  :: cls_id, ncls, ncomp
        complex(sp),      intent(in)  :: y(:,:)
        real(sp),         intent(in)  :: c(:,:), w(:,:), wq(:)
        logical,          intent(in)  :: fitmask(:), even(:), l_models
        integer,          intent(in)  :: shell(:), hidx(:), kidx(:), labels(:)
        complex(sp),      intent(out) :: parh(:,:), devh(:,:,:)   !< (ncoeff,2), (ncoeff,ncls,2)
        type(flex_cls_model)     :: mh
        complex(sp), allocatable :: yh(:,:), nuis(:,:)
        real(sp),    allocatable :: ch(:,:), wh(:,:)
        complex(dp), allocatable :: uext(:,:), mu_h(:)
        integer,     allocatable :: idx(:)
        integer :: ncoeff, n, h, i
        ncoeff = size(y,1); n = size(y,2)
        allocate(nuis(ncoeff,4), mu_h(ncoeff))
        devh = CMPLX_ZERO
        do h = 1, 2
            idx = pack([(i, i=1,n)], even .eqv. (h == 1))
            allocate(yh(ncoeff,size(idx)), ch(ncoeff,size(idx)), wh(ncoeff,size(idx)))
            do i = 1, size(idx)
                yh(:,i) = y(:,idx(i)); ch(:,i) = c(:,idx(i)); wh(:,i) = w(:,idx(i))
            end do
            call flex_cls_ctf_mean(yh, ch, wh, parh(:,h))
            if( l_models )then
                call flex_cls_weighted_mean(yh, ch, wh, mu_h)
                call flex_cls_pose_tangents(mu_h, hidx, kidx, params%box, nuis(:,1:3))
                nuis(:,4) = cmplx(mu_h, kind=sp)
                call seed_rnd_fixed(FLEX_CLS_SEED_BASE + cls_id + 1000 * h)
                call flex_cls_fit(mh, yh, ch, wh, wq, fitmask, ncomp, verbose=.false., nuis=nuis)
                allocate(uext(ncoeff,mh%ncomp))
                call flex_cls_extend_basis(mh, yh, ch, wh, shell, uext)
                call flex_cls_conditional_averages(mh, labels(idx), uext, devh(:,:,h))
                call mh%kill
                deallocate(uext)
            endif
            deallocate(yh, ch, wh)
        end do
        deallocate(nuis, mu_h)
    end subroutine flex_class_half_models

    !> one class: its members as Fourier planes in the class frame; a class kept whole (nsub = 1) gives its CTF-corrected
    !! parent (with the half parents and the same filter); otherwise the bands (data-derived without lp=), the model,
    !! placement and the conditional averages of its nsub subclasses
    subroutine split_class_flex( params, build, spproj, cls_id, nsub, ncomp, kfit, l_sigma, pinds, labels, cavgs, &
                                 separation, repro, cavgs_even, cavgs_odd, frcs, half_pops )
        use simple_flex_cls_expansion,  only: flex_cls_model, flex_cls_fit, flex_cls_half_frc, flex_cls_frc_cutoff, &
            &flex_cls_extend_basis, flex_cls_weighted_mean, flex_cls_pose_tangents, flex_cls_split_floor
        use simple_rnd,             only: seed_rnd_fixed
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        integer,                  intent(in)    :: cls_id, nsub, ncomp, kfit
        logical,                  intent(in)    :: l_sigma
        integer,     allocatable, intent(out)   :: pinds(:), labels(:)
        type(image), allocatable, intent(inout) :: cavgs(:)
        real,                     intent(out)   :: separation
        real,        allocatable, intent(out)   :: repro(:)
        type(image), allocatable, intent(inout) :: cavgs_even(:), cavgs_odd(:)
        real,        allocatable, intent(inout) :: frcs(:,:)
        integer,     allocatable, intent(out)   :: half_pops(:,:) !< (2,nsub) members of each subclass in the even, odd half
        type(image)          :: gridcorr_img
        type(flex_cls_model) :: model
        complex(sp), allocatable :: y(:,:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        logical,     allocatable :: fitmask(:), even(:)
        integer :: ldim(3), ncoeff, nptcls, j
        separation = 0.
        if( allocated(cavgs) ) call dealloc_imgarr(cavgs)
        call flex_class_planes(params, build, spproj, cls_id, kfit, l_sigma, pinds, y, c, w, wq, &
            &hidx, kidx, shell, fitmask, gridcorr_img, ldim, ncoeff, nptcls)
        if( nptcls < 1 ) THROW_HARD('cls_expansion: a class has no members')
        ! the halves: the project's eo of the members (the particles' gold-standard halves), else members alternately
        allocate(even(nptcls))
        do j = 1, nptcls
            even(j) = nint(spproj%os_ptcl2D%get(pinds(j), 'eo')) == 0
        end do
        if( all(even) .or. .not. any(even) )then
            even = [(mod(j,2) == 0, j=1,nptcls)]
            write(logfhandle,'(A,I8,A)') 'Cls expansion: class=', cls_id, ' has no eo halves; members alternately'
        endif
        if( nsub == 1 )then
            allocate(labels(nptcls), source=1)
            call flex_class_conditional_averages(params, cls_id, 1, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, &
                &labels, gridcorr_img, ldim, cavgs, cavgs_even, cavgs_odd, frcs, repro)
        else
            if( nptcls < flex_cls_split_floor(nsub, ncomp) ) THROW_HARD('cls_expansion: a class is too small for its subclass count')
            call fit_class_model
            call flex_class_deliver(params, cls_id, model, nsub, ncomp, y, c, w, wq, fitmask, shell, hidx, kidx, even, &
                &gridcorr_img, ldim, nptcls, labels, cavgs, cavgs_even, cavgs_odd, frcs, separation, repro)
            call model%kill
        endif
        ! a subclass half without members is delivered inactive (state 0 in cls2D)
        allocate(half_pops(2,nsub))
        do j = 1, nsub
            half_pops(1,j) = count(labels == j .and. even)
            half_pops(2,j) = count(labels == j .and. .not. even)
        end do
        call gridcorr_img%kill
        deallocate(y, c, w, wq, hidx, kidx, shell, fitmask, even)

      contains

        !> the class model: the bands, then one in-sample fit on the global band with the pose tangents and the class
        !! mean (per-particle contrast) as nuisance columns
        subroutine fit_class_model
            type(flex_cls_model)     :: model0
            complex(dp), allocatable :: uext(:,:), mu_all(:)
            complex(sp), allocatable :: nuis(:,:)
            real(dp),    allocatable :: frc(:)
            integer(timer_int_kind)  :: t0
            integer :: nsh, k0, kres, kglob
            t0  = tic()
            nsh = maxval(shell)
            ! nuisance columns: the in-plane pose tangents of the class mean (d/dx, d/dy, d/dtheta) and the class mean
            ! (per-particle contrast), fitted jointly and never placed on or averaged over
            allocate(mu_all(ncoeff), nuis(ncoeff,4))
            call flex_cls_weighted_mean(y, c, w, mu_all)
            call flex_cls_pose_tangents(mu_all, hidx, kidx, params%box, nuis(:,1:3))
            nuis(:,4) = cmplx(mu_all, kind=sp)
            deallocate(mu_all)
            if( kfit > 0 )then
                ! lp= given: the global band
                kglob = kfit
            else
                ! data-derived: a first fit to half the class resolution (at most FLEX_CLS_BAND_FIRST_LP); the even/odd
                ! FRC of its basis on every shell sets the global band at 0.5
                allocate(frc(nsh))
                call flex_cls_half_frc(y, c, w, shell, nsh, frc, even)
                kres = flex_cls_frc_cutoff(frc, 0.143d0)
                k0   = max(4, min(kres / 2, calc_fourier_index(FLEX_CLS_BAND_FIRST_LP, params%box, params%smpd)))
                fitmask = shell >= 1 .and. shell <= k0
                call seed_rnd_fixed(FLEX_CLS_SEED_BASE + cls_id)
                call flex_cls_fit(model0, y, c, w, wq, fitmask, ncomp, verbose=.false., nuis=nuis)
                allocate(uext(ncoeff,ncomp))
                call flex_cls_extend_basis(model0, y, c, w, shell, uext, frc_basis=frc, even=even)
                deallocate(uext)
                call model0%kill
                ! never finer than the first band: reproducible fine variation would crowd out the coarse motion
                kglob = max(max(4, k0 / 2), min(k0, flex_cls_frc_cutoff(frc, 0.5d0)))
                write(logfhandle,'(A,I8,A,F7.1,A,F7.1,A,F7.1,A)') 'Cls expansion bands: class=', cls_id, &
                    &' class resolution ', params%box * params%smpd / real(max(kres,1)), ' A; first fit to ', &
                    &params%box * params%smpd / real(k0), ' A; global latent to ', &
                    &params%box * params%smpd / real(kglob), ' A'
                deallocate(frc)
                call flex_time(cls_id, 'bands', t0)
            endif
            fitmask = shell >= 1 .and. shell <= kglob
            call seed_rnd_fixed(FLEX_CLS_SEED_BASE + cls_id)
            call flex_cls_fit(model, y, c, w, wq, fitmask, ncomp, verbose=.true., nuis=nuis)
            write(logfhandle,'(A,4F9.4)') 'Cls expansion flex nuisance rms (x, y, rot, contrast; latent units): ', &
                &[(sqrt(sum(model%nu(:,j)**2)/real(nptcls,dp)), j=1,4)]
            deallocate(nuis)
            call flex_time(cls_id, 'fit', t0)
        end subroutine fit_class_model

    end subroutine split_class_flex

    function flex_frcs_part_fname( part, numlen ) result( fname )
        integer, intent(in) :: part, numlen
        type(string) :: fname
        fname = string('cls_expansion_frcs_part')//int2str_pad(part, numlen)//TXT_EXT
    end function flex_frcs_part_fname

    !> the conditional averages' even/odd FRCs (text rows 'row frc(1:nyq)') as class FRCs (FLEX_CLS_FRCS, frc2D);
    !! rows map to global subclasses directly, or through row_global(:,ifile) for part files
    subroutine flex_frcs_txt2bin( params, fnames, ntot, row_global )
        use simple_class_frcs, only: class_frcs
        type(parameters),  intent(in) :: params
        type(string),      intent(in) :: fnames(:)
        integer,           intent(in) :: ntot
        integer, optional, intent(in) :: row_global(:,:)
        type(class_frcs)  :: frcs
        real, allocatable :: frc(:)
        character(len=XLONGSTRLEN) :: line
        integer :: ifile, funit, ios, row, iglob, nyq
        if( ntot < 1 ) return
        nyq = fdim(params%box) - 1
        call frcs%new(ntot, params%box, params%smpd)
        allocate(frc(nyq))
        do ifile = 1, size(fnames)
            if( .not. file_exists(fnames(ifile)%to_char()) ) cycle
            open(newunit=funit, file=fnames(ifile)%to_char(), status='old', action='read', iostat=ios)
            if( ios /= 0 ) cycle
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*,iostat=ios) row, frc
                if( ios /= 0 ) cycle
                iglob = row
                if( present(row_global) ) iglob = row_global(row, ifile)
                if( iglob >= 1 .and. iglob <= ntot ) call frcs%set_frc(iglob, frc)
            end do
            close(funit)
        end do
        call frcs%write(string(FLEX_CLS_FRCS))
        ! and as text in global subclass order (the part files map their rows)
        if( present(row_global) )then
            open(newunit=funit, file=FLEX_CLS_FRCS_TXT, status='replace', action='write')
            write(funit,'(A)') '# row frc(1:nyq): even/odd FRC of the conditional average (row: global)'
            do iglob = 1, ntot
                write(funit,'(I8,*(1X,F8.5))') iglob, frcs%get_frc(iglob, params%box)
            end do
            close(funit)
        endif
        call frcs%kill
        deallocate(frc)
    end subroutine flex_frcs_txt2bin

    !> the split into the project: hard subclass and parent per particle, one cls2D/cls3D row per subclass with its
    !! readouts, the sub-class average stack and its FRCs (frc2D)
    subroutine apply_split_project_updates(spproj, params, nsplit, new_class, new_parent, parent_of_subcls, pop_of_subcls, &
                                           neff_of_subcls, sep_of_subcls, repro_of_subcls, pop_even_of_subcls, &
                                           pop_odd_of_subcls, cavgs_stk)
        type(sp_project), intent(inout) :: spproj
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nsplit
        integer,          intent(in)    :: new_class(:), new_parent(:)
        integer,          intent(in)    :: parent_of_subcls(:), pop_of_subcls(:)
        real,             intent(in)    :: neff_of_subcls(:)  !< the hard population
        real,             intent(in)    :: sep_of_subcls(:)   !< Ashman's D of the parent's split (same for its subclasses)
        real,             intent(in)    :: repro_of_subcls(:) !< cross-half reproducibility of the subclass difference
        integer,          intent(in)    :: pop_even_of_subcls(:), pop_odd_of_subcls(:) !< members in the even, odd half
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
            ! a half without members is inactive: its average carries no member of the subclass
            call spproj%os_cls2D%set(i, 'state_even', merge(1, 0, pop_even_of_subcls(i) > 0))
            call spproj%os_cls2D%set(i, 'state_odd',  merge(1, 0, pop_odd_of_subcls(i)  > 0))
            call spproj%os_cls3D%set(i, 'cluster', parent_of_subcls(i))
            call spproj%os_cls3D%set(i, 'accept',  1)
            call spproj%os_cls3D%set(i, 'state',   1)
        end do
        call spproj%add_frcs2os_out(string(FLEX_CLS_FRCS), 'frc2D')
        call spproj%add_cavgs2os_out(cavgs_stk, params%smpd, imgkind='cavg')
        call spproj%write(params%projfile)
    end subroutine apply_split_project_updates

    !> the part outputs in global subclass order (by parent, then local subclass): class map, assignments into the
    !! project, the average stacks and the per-subclass FRCs
    subroutine merge_worker_outputs(params, nparts_run)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        type(sp_project) :: spproj
        integer, allocatable :: part_counts(:), part_parent(:,:), part_local(:,:), part_pop(:,:), part_global(:,:)
        real,    allocatable :: part_neff(:,:), part_sep(:,:), part_repro(:,:)
        integer, allocatable :: part_pop_even(:,:), part_pop_odd(:,:), pop_even(:), pop_odd(:)
        integer, allocatable :: comb_part(:), comb_row(:), comb_parent(:), comb_local(:), comb_pop(:)
        integer, allocatable :: new_class(:), new_parent(:), first_glob(:), nsub_of(:)
        real,    allocatable :: neff(:), sep(:), repro(:)
        type(string) :: map_fname, assign_fname
        type(string), allocatable :: frc_fnames(:)
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
        allocate(part_pop_even(max_count, nparts_run), part_pop_odd(max_count, nparts_run), source=0)
        allocate(comb_part(total), comb_row(total), comb_parent(total), comb_local(total), comb_pop(total), source=0)
        idx = 0
        do ipart = 1, nparts_run
            map_fname = string('cls_expansion_class_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call read_part_map(map_fname, part_counts(ipart), part_parent(:,ipart), part_local(:,ipart), part_pop(:,ipart), &
                &part_neff(:,ipart), part_sep(:,ipart), part_repro(:,ipart), part_pop_even(:,ipart), part_pop_odd(:,ipart))
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
        allocate(pop_even(total), pop_odd(total), source=0)
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
            pop_even(idx) = part_pop_even(comb_row(idx), comb_part(idx))
            pop_odd(idx)  = part_pop_odd(comb_row(idx), comb_part(idx))
        end do
        map_fname = string(FLEX_CLS_MAP_FILE)
        open(newunit=funit, file=map_fname%to_char(), status='replace', action='write')
        write(funit,'(A)') '# subclass parent_class local_subclass pop neff separation repro pop_even pop_odd'
        do idx = 1, total
            write(funit,'(I8,1X,I8,1X,I8,1X,I8,1X,F12.3,1X,F8.3,1X,F8.3,2(1X,I8))') idx, comb_parent(idx), comb_local(idx), &
                &comb_pop(idx), neff(idx), sep(idx), repro(idx), pop_even(idx), pop_odd(idx)
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
        call merge_part_stacks(params, nparts_run, part_counts, part_global)
        ! the per-subclass FRCs, part rows -> global subclasses
        allocate(frc_fnames(nparts_run))
        do ipart = 1, nparts_run
            frc_fnames(ipart) = flex_frcs_part_fname(ipart, params%numlen)
        end do
        call flex_frcs_txt2bin(params, frc_fnames, total, row_global=part_global(:,1:nparts_run))
        call frc_fnames(:)%kill
        deallocate(frc_fnames)
        call apply_split_project_updates(spproj, params, total, new_class, new_parent, comb_parent, comb_pop, &
            &neff, sep, repro, pop_even, pop_odd, string(FLEX_CLS_CAVGS_FILE))
        call spproj%kill
        deallocate(new_class, new_parent, part_counts, part_parent, part_local, part_pop, part_global, part_neff, part_sep, &
                   part_repro, comb_part, comb_row, comb_parent, comb_local, comb_pop, first_glob, nsub_of, neff, sep, repro, &
                   part_pop_even, part_pop_odd, pop_even, pop_odd)

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

    !> the part average stacks (full, even, odd) in global subclass order: part row r -> global index part_global(r, ipart)
    subroutine merge_part_stacks(params, nparts_run, part_counts, part_global)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run, part_counts(:), part_global(:,:)
        type(string) :: stk_fname, out_stk
        type(image)  :: img
        integer :: ipart, irow, ldim(3), nimgs, ieo
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
        call out_stk%kill
    end subroutine merge_part_stacks

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

    subroutine read_part_map(fname, nrows, parents, locals, pops, neffs, seps, repros, pops_even, pops_odd)
        type(string), intent(in)  :: fname
        integer,      intent(in)  :: nrows
        integer,      intent(out) :: parents(:), locals(:), pops(:), pops_even(:), pops_odd(:)
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
            read(line,*) row, parents(irow), locals(irow), pops(irow), neffs(irow), seps(irow), repros(irow), &
                &pops_even(irow), pops_odd(irow)
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
