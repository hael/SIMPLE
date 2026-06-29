module simple_denoise_project_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_cmdline,            only: cmdline
use simple_qsys_env,           only: qsys_env
use simple_sp_project,         only: sp_project
use simple_image,              only: image
use simple_image_msk,          only: automask2D
use simple_diff_map_graphs,    only: diffmap_graph, graph_matvec
use simple_linalg,             only: sparse_eigh
use simple_eulspace_neigh_map, only: eulspace_neigh_map
use simple_classaverager,      only: transform_ptcls
use simple_imgfile,            only: imgfile
use simple_srch_sort_loc,      only: hpsort
implicit none

public :: denoise_project_strategy
public :: denoise_project_shmem_strategy
public :: denoise_project_worker_strategy
public :: denoise_project_master_strategy
public :: create_denoise_project_strategy

private
#include "simple_local_flags.inc"

integer, parameter :: DIFFMAP_DENOISE_ICM_RANK_MAXITS = 16
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_BETA_FRAC = 0.35
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_COMPLEXITY_FRAC = 0.10
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_LOWER_SEED_FRAC = 0.50

type, abstract :: denoise_project_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type denoise_project_strategy

type, extends(denoise_project_strategy) :: denoise_project_shmem_strategy
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type denoise_project_shmem_strategy

type, extends(denoise_project_strategy) :: denoise_project_worker_strategy
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
end type denoise_project_worker_strategy

type, extends(denoise_project_strategy) :: denoise_project_master_strategy
    type(qsys_env)           :: qenv
    type(chash)              :: job_descr
    type(chash), allocatable :: part_params(:)
    integer                  :: nparts_run = 1
    integer                  :: nthr_master = 1
contains
    procedure :: initialize   => master_initialize
    procedure :: execute      => master_execute
    procedure :: finalize_run => master_finalize_run
    procedure :: cleanup      => master_cleanup
end type denoise_project_master_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: denoise_project_strategy, parameters, builder, cmdline
        class(denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(inout) :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: denoise_project_strategy, parameters, builder, cmdline
        class(denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(inout) :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: denoise_project_strategy, parameters, builder, cmdline
        class(denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(in)    :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params)
        import :: denoise_project_strategy, parameters
        class(denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    function create_denoise_project_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(denoise_project_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker, is_master
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        is_worker = cline%defined('part')
        is_master = (nparts > 1) .and. (.not. is_worker)
        if( is_master )then
            allocate(denoise_project_master_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED DENOISE_PROJECT (MASTER)'
        else if( is_worker )then
            allocate(denoise_project_worker_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DENOISE_PROJECT (WORKER)'
        else
            allocate(denoise_project_shmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DENOISE_PROJECT (SHARED-MEMORY)'
        endif
    end function create_denoise_project_strategy

    subroutine apply_defaults(cline)
        use simple_default_clines, only: set_automask2D_defaults
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype',  'ptcl2D')
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'diffusion_maps')
        if( .not. cline%defined('graph')    ) call cline%set('graph',    'euc')
        if( .not. cline%defined('steering') ) call cline%set('steering', 'none')
        if( .not. cline%defined('k_nn')     ) call cline%set('k_nn',      10)
        if( .not. cline%defined('neigs')    ) call cline%set('neigs',     200)
        call set_automask2D_defaults(cline)
    end subroutine apply_defaults

    subroutine init_common(params, build, cline)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call apply_defaults(cline)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
    end subroutine init_common

    subroutine shmem_initialize(self, params, build, cline)
        class(denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                      intent(inout) :: params
        type(builder),                         intent(inout) :: build
        class(cmdline),                        intent(inout) :: cline
        call init_common(params, build, cline)
    end subroutine shmem_initialize

    subroutine shmem_execute(self, params, build, cline)
        class(denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                      intent(inout) :: params
        type(builder),                         intent(inout) :: build
        class(cmdline),                        intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_denoise_project(params, build, cline, spproj, 0, .true.)
        call spproj%kill
    end subroutine shmem_execute

    subroutine shmem_finalize_run(self, params, build, cline)
        class(denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                      intent(in)    :: params
        type(builder),                         intent(inout) :: build
        class(cmdline),                        intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup(self, params)
        class(denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                      intent(in)    :: params
    end subroutine shmem_cleanup

    subroutine worker_initialize(self, params, build, cline)
        class(denoise_project_worker_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        call init_common(params, build, cline)
        if( .not. cline%defined('part') ) THROW_HARD('PART must be defined for denoise_project worker execution')
        if( .not. cline%defined('class_assignment') )then
            THROW_HARD('CLASS_ASSIGNMENT must be defined for denoise_project worker execution')
        endif
    end subroutine worker_initialize

    subroutine worker_execute(self, params, build, cline)
        use simple_qsys_funs, only: qsys_job_finished
        class(denoise_project_worker_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_denoise_project(params, build, cline, spproj, params%part, .false.)
        call qsys_job_finished(params, string('simple_denoise_project_strategy :: worker_execute'))
        call spproj%kill
    end subroutine worker_execute

    subroutine worker_finalize_run(self, params, build, cline)
        class(denoise_project_worker_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
    end subroutine worker_finalize_run

    subroutine worker_cleanup(self, params)
        class(denoise_project_worker_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
    end subroutine worker_cleanup

    subroutine master_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(denoise_project_master_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        type(sp_project) :: spproj
        integer, allocatable :: cls_inds(:), cls_pops(:)
        logical :: l_phflip, l_ctf_no
        call init_common(params, build, cline)
        call spproj%read(params%projfile)
        call validate_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip, l_ctf_no)
        if( trim(lowercase(params%graph)) == 'ori' )then
            self%nparts_run = 1
        else
            self%nparts_run = min(max(1, params%nparts), size(cls_inds))
        endif
        write(logfhandle,'(A,I8,A,I8,A,I8)') 'Denoise worker planning: requested_nparts=', params%nparts, &
            ' eligible_classes=', size(cls_inds), ' effective_nparts=', self%nparts_run
        call flush(logfhandle)
        if( self%nparts_run < params%nparts )then
            write(logfhandle,'(A)') 'Denoise note: reducing worker count because classes are the scheduling unit.'
            call flush(logfhandle)
        endif
        call set_master_num_threads(self%nthr_master, string('DENOISE_PROJECT'))
        call prepare_class_partitions(self, params, cls_inds, cls_pops)
        call self%qenv%new(params, self%nparts_run, numlen=params%numlen)
        call cline%set('mkdir', 'no')
        call cline%gen_job_descr(self%job_descr)
        call self%job_descr%set('mkdir', 'no')
        call self%job_descr%set('nparts', int2str(self%nparts_run))
        call self%job_descr%set('numlen', int2str(params%numlen))
        call spproj%kill
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
    end subroutine master_initialize

    subroutine master_execute(self, params, build, cline)
        class(denoise_project_master_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, part_params=self%part_params, &
                                                     array=L_USE_SLURM_ARR, extra_params=params)
        call finalize_diffmap_worker_outputs(params, self%nparts_run)
    end subroutine master_execute

    subroutine master_finalize_run(self, params, build, cline)
        class(denoise_project_master_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
    end subroutine master_finalize_run

    subroutine master_cleanup(self, params)
        use simple_qsys_funs, only: qsys_cleanup
        class(denoise_project_master_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
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

    subroutine prepare_class_partitions(self, params, cls_inds, cls_pops)
        class(denoise_project_master_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
        integer,                                intent(in)    :: cls_inds(:), cls_pops(:)
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
            fname = string('denoise_classes_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            call arr2txtfile(part_cls(1:part_counts(ipart), ipart), fname)
            call self%part_params(ipart)%new(1)
            call self%part_params(ipart)%set('class_assignment', fname%to_char())
            call fname%kill
        end do
        deallocate(part_counts, part_cls, part_weights)
    end subroutine prepare_class_partitions

    subroutine run_denoise_project(params, build, cline, spproj, part, l_write_project)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: part
        logical,          intent(in)    :: l_write_project
        type(sp_project) :: outproj
        type(string), allocatable :: raw_stks(:), den_stks(:)
        type(string) :: map_fname
        type(image)  :: img_blank
        integer, allocatable :: cls_inds(:), cls_pops(:)
        logical, allocatable :: processed(:)
        logical :: l_phflip, l_ctf_no, l_img_blank_init, l_mskdiam_override
        integer :: icls, iptcl, nptcls, nprocessed, funit_map, ldim_blank(3)
        call validate_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip, l_ctf_no)
        call filter_classes_by_assignment(cline, cls_inds, cls_pops)
        if( trim(params%pca_mode) == 'diffusion_maps' .and. trim(lowercase(params%graph)) == 'ori' )then
            call run_diffmap_so3_mixture_graphs(params, build, spproj, cls_inds)
            if( allocated(cls_inds) ) deallocate(cls_inds)
            if( allocated(cls_pops) ) deallocate(cls_pops)
            return
        endif
        l_mskdiam_override = cline%defined('mskdiam')
        nptcls = spproj%os_ptcl2D%get_noris()
        allocate(processed(nptcls), source=.false.)
        l_img_blank_init = .false.
        funit_map = -1
        allocate(raw_stks(1), den_stks(1))
        if( l_write_project )then
            raw_stks(1) = string('denoise_raw_particles.mrcs')
            den_stks(1) = string('denoise_particles.mrcs')
            map_fname   = string('denoise_particle_map.txt')
        else
            raw_stks(1) = string('denoise_raw_particles_part')//int2str_pad(part, params%numlen)//'.mrcs'
            den_stks(1) = string('denoise_particles_part')//int2str_pad(part, params%numlen)//'.mrcs'
            map_fname = string('denoise_particle_map_part')//int2str_pad(part, params%numlen)//TXT_EXT
        endif
        call del_file(raw_stks(1))
        call del_file(den_stks(1))
        call del_file(map_fname)
        open(newunit=funit_map, file=map_fname%to_char(), status='replace', action='write')
        write(funit_map,'(A)') '# particle_index stack_index stack_local_index output_stack_index'
        nprocessed = 0
        do icls = 1, size(cls_inds)
            call write_diffmap_denoised_class(params, build, spproj, cls_inds(icls), l_phflip, &
                l_mskdiam_override, raw_stks, den_stks, processed, nprocessed, .true., funit_map, l_write_project, icls)
        end do
        if( l_write_project )then
            ldim_blank = [params%box_crop, params%box_crop, 1]
            call img_blank%new(ldim_blank, params%smpd_crop, wthreads=.false.)
            l_img_blank_init = .true.
            call img_blank%zero_and_unflag_ft
            do iptcl = 1, nptcls
                if( spproj%os_ptcl2D%get_state(iptcl) > 0 ) cycle
                call write_diffmap_stack_image(raw_stks(1), iptcl, img_blank)
                call write_diffmap_stack_image(den_stks(1), iptcl, img_blank)
                processed(iptcl) = .true.
            end do
        endif
        close(funit_map)
        if( l_write_project .and. any(.not. processed) )then
            write(logfhandle,'(A,I10)') 'unprocessed particles: ', count(.not. processed)
            THROW_HARD('denoise_project requires every active ptcl2D particle to have a valid class assignment')
        endif
        if( l_write_project )then
            call write_diffmap_project(params, spproj, raw_stks(1), den_stks(1), l_ctf_no, outproj)
            write(logfhandle,'(A,A)') 'Denoise project written: ', params%projfile%to_char()
        endif
        write(logfhandle,'(A,I10)') 'Denoise project particles: ', nprocessed
        call flush(logfhandle)
        call outproj%kill
        if( allocated(raw_stks) ) deallocate(raw_stks)
        if( allocated(den_stks) ) deallocate(den_stks)
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
        if( allocated(processed) ) deallocate(processed)
        if( l_img_blank_init ) call img_blank%kill
        if( map_fname%is_allocated() ) call map_fname%kill
    end subroutine run_denoise_project

    subroutine validate_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip, l_ctf_no)
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer, allocatable, intent(out)   :: cls_inds(:), cls_pops(:)
        logical,              intent(out)   :: l_phflip, l_ctf_no
        type(string) :: ctfstr
        integer, allocatable :: pinds(:), stk_offsets(:), stk_nptcls(:)
        logical, allocatable :: seen_slots(:)
        integer :: iptcl, istk, icls, nptcls, nstks, cls, stkind, indstk, slot, ctf_mode, nptcls_slots, nptcls_active
        if( trim(params%oritype) /= 'ptcl2D' )then
            THROW_HARD('denoise_project supports oritype=ptcl2D input only')
        endif
        if( trim(params%pca_mode) /= 'diffusion_maps' )then
            THROW_HARD('denoise_project supports pca_mode=diffusion_maps only')
        endif
        select case(trim(lowercase(params%graph)))
            case('euc','ori')
            case DEFAULT
                THROW_HARD('denoise_project supports graph=euc|ori for pca_mode=diffusion_maps')
        end select
        if( trim(params%steering) /= 'none' )then
            THROW_HARD('denoise_project supports non-steerable diffusion maps only; use steering=none')
        endif
        if( trim(params%pca_mode) == 'diffusion_maps' )then
            if( trim(lowercase(params%graph)) == 'ori' )then
                if( params%nspace < 2 ) THROW_HARD('denoise_project graph=ori requires nspace >= 2')
                if( params%nspace_sub < 1 ) THROW_HARD('denoise_project graph=ori requires nspace_sub >= 1')
                if( params%nspace_sub > params%nspace )then
                    THROW_HARD('denoise_project graph=ori requires nspace_sub <= nspace')
                endif
                if( mod(params%nspace, 2) /= 0 .or. mod(params%nspace_sub, 2) /= 0 )then
                    THROW_HARD('denoise_project graph=ori requires even nspace and nspace_sub')
                endif
            endif
        endif
        if( params%k_nn < 1 ) THROW_HARD('denoise_project requires k_nn >= 1')
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('empty ptcl2D field; denoise_project')
        if( spproj%os_ptcl3D%get_noris() /= nptcls )then
            THROW_HARD('ptcl2D/ptcl3D particle counts differ; denoise_project')
        endif
        nstks = spproj%os_stk%get_noris()
        if( nstks < 1 ) THROW_HARD('empty stack field; denoise_project')
        allocate(stk_offsets(nstks), stk_nptcls(nstks), source=0)
        nptcls_slots = 0
        ctf_mode = 0
        do istk = 1, nstks
            stk_offsets(istk) = nptcls_slots
            stk_nptcls(istk) = get_diffmap_stack_nptcls(spproj, istk)
            nptcls_slots = nptcls_slots + stk_nptcls(istk)
            if( .not. spproj%os_stk%isthere(istk, 'ctf') ) THROW_HARD('stack ctf key missing; denoise_project')
            call spproj%os_stk%getter(istk, 'ctf', ctfstr)
            select case(trim(ctfstr%to_char()))
                case('yes')
                    if( ctf_mode == 0 ) ctf_mode = CTFFLAG_YES
                    if( ctf_mode /= CTFFLAG_YES ) THROW_HARD('mixed ctf=yes/ctf=flip stacks are not supported; denoise_project')
                case('flip')
                    if( ctf_mode == 0 ) ctf_mode = CTFFLAG_FLIP
                    if( ctf_mode /= CTFFLAG_FLIP ) THROW_HARD('mixed ctf=yes/ctf=flip stacks are not supported; denoise_project')
                case('no')
                    if( ctf_mode == 0 ) ctf_mode = CTFFLAG_NO
                    if( ctf_mode /= CTFFLAG_NO ) THROW_HARD('mixed ctf=no with ctf=yes/ctf=flip stacks are not supported; denoise_project')
                case DEFAULT
                    THROW_HARD('denoise_project requires ctf=yes, ctf=flip, or ctf=no input stacks')
            end select
        enddo
        l_phflip = ctf_mode == CTFFLAG_YES
        l_ctf_no = ctf_mode == CTFFLAG_NO
        allocate(seen_slots(nptcls_slots), source=.false.)
        nptcls_active = 0
        do iptcl = 1, nptcls
            if( spproj%os_ptcl2D%get_state(iptcl) <= 0 ) cycle
            nptcls_active = nptcls_active + 1
            if( .not. spproj%os_ptcl2D%isthere(iptcl, 'class') )then
                THROW_HARD('ptcl2D class assignment missing; denoise_project')
            endif
            cls = spproj%os_ptcl2D%get_int(iptcl, 'class')
            if( cls <= 0 ) THROW_HARD('ptcl2D class assignment must be positive; denoise_project')
            call spproj%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, indstk)
            if( stkind < 1 .or. stkind > nstks ) THROW_HARD('invalid particle stkind; denoise_project')
            if( indstk < 1 .or. indstk > stk_nptcls(stkind) )then
                THROW_HARD('invalid particle indstk; denoise_project')
            endif
            slot = stk_offsets(stkind) + indstk
            if( slot < 1 .or. slot > nptcls_slots ) THROW_HARD('invalid particle stack slot; denoise_project')
            if( seen_slots(slot) ) THROW_HARD('duplicate particle stack slot; denoise_project')
            seen_slots(slot) = .true.
        enddo
        if( nptcls_active < 1 ) THROW_HARD('No active ptcl2D particles found; denoise_project')
        cls_inds = spproj%os_ptcl2D%get_label_inds('class')
        cls_inds = pack(cls_inds, mask=cls_inds > 0)
        if( size(cls_inds) < 1 ) THROW_HARD('No positive ptcl2D classes found; denoise_project')
        allocate(cls_pops(size(cls_inds)), source=0)
        do icls = 1, size(cls_inds)
            call spproj%os_ptcl2D%get_pinds(cls_inds(icls), 'class', pinds)
            if( allocated(pinds) )then
                cls_pops(icls) = size(pinds)
                deallocate(pinds)
            endif
            if( cls_pops(icls) == 0 ) cycle
            if( cls_pops(icls) < 3 )then
                write(logfhandle,'(A,I8,A,I8)') 'class=', cls_inds(icls), ' pop=', cls_pops(icls)
                THROW_HARD('denoise_project requires at least 3 particles per class')
            endif
        end do
        cls_inds = pack(cls_inds, mask=cls_pops > 0)
        cls_pops = pack(cls_pops, mask=cls_pops > 0)
        if( size(cls_inds) < 1 ) THROW_HARD('No active ptcl2D classes found; denoise_project')
        if( sum(cls_pops) /= nptcls_active ) THROW_HARD('class populations do not cover all active particles; denoise_project')
        if( allocated(seen_slots) ) deallocate(seen_slots)
        if( allocated(stk_nptcls) ) deallocate(stk_nptcls)
        if( allocated(stk_offsets) ) deallocate(stk_offsets)
        call ctfstr%kill
    end subroutine validate_denoise_project

    integer function get_diffmap_stack_nptcls(spproj, istk) result(nptcls_stk)
        type(sp_project), intent(in) :: spproj
        integer,          intent(in) :: istk
        integer :: fromp, top
        nptcls_stk = 0
        if( spproj%os_stk%isthere(istk, 'nptcls_stk') )then
            nptcls_stk = spproj%os_stk%get_int(istk, 'nptcls_stk')
        else if( spproj%os_stk%isthere(istk, 'fromp') .and. spproj%os_stk%isthere(istk, 'top') )then
            fromp = spproj%os_stk%get_fromp(istk)
            top   = spproj%os_stk%get_top(istk)
            nptcls_stk = top - fromp + 1
        endif
        if( nptcls_stk < 1 )then
            write(logfhandle,'(A,I8)') 'invalid stack particle count at stkind=', istk
            THROW_HARD('denoise_project requires stack nptcls_stk or valid fromp/top range')
        endif
    end function get_diffmap_stack_nptcls

    integer function get_n_active_ptcl2D(spproj) result(n_active)
        type(sp_project), intent(in) :: spproj
        integer :: iptcl
        n_active = 0
        do iptcl = 1, spproj%os_ptcl2D%get_noris()
            if( spproj%os_ptcl2D%get_state(iptcl) > 0 ) n_active = n_active + 1
        end do
    end function get_n_active_ptcl2D

    subroutine filter_classes_by_assignment(cline, cls_inds, cls_pops)
        class(cmdline),       intent(inout) :: cline
        integer, allocatable, intent(inout) :: cls_inds(:), cls_pops(:)
        integer, allocatable :: assigned_classes(:)
        logical, allocatable :: keep_mask(:)
        integer :: i
        if( .not. cline%defined('class_assignment') ) return
        call read_int_file(cline%get_carg('class_assignment'), assigned_classes)
        allocate(keep_mask(size(cls_inds)), source=.false.)
        do i = 1, size(assigned_classes)
            keep_mask = keep_mask .or. (cls_inds == assigned_classes(i))
        end do
        cls_inds = pack(cls_inds, mask=keep_mask)
        cls_pops = pack(cls_pops, mask=keep_mask)
        deallocate(keep_mask)
        if( allocated(assigned_classes) ) deallocate(assigned_classes)
        if( size(cls_inds) < 1 ) THROW_HARD('No classes selected for denoise_project')
    end subroutine filter_classes_by_assignment

    subroutine write_diffmap_denoised_class(params, build, spproj, cls_id, l_phflip, l_mskdiam_override, raw_stks, den_stks, &
                                            processed, nprocessed, l_dense_output, map_unit, l_index_by_pind, class_index)
        use simple_diff_map_graphs,  only: diffmap_graph, build_cls_split_graph
        use simple_imgarr_utils,     only: dealloc_imgarr, copy_imgarr
        use simple_imgproc,          only: make_pcavecs
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: cls_id
        logical,          intent(in)    :: l_phflip, l_mskdiam_override
        type(string),     intent(in)    :: raw_stks(:), den_stks(:)
        logical,          intent(inout) :: processed(:)
        integer,          intent(inout) :: nprocessed
        logical,          intent(in)    :: l_dense_output
        integer,          intent(in)    :: map_unit
        logical,          intent(in)    :: l_index_by_pind
        integer,          intent(in)    :: class_index
        type(parameters)         :: params_mask
        type(image), allocatable :: imgs(:), imgs_ppca(:), class_mask(:), den_ptcls(:)
        real,        allocatable :: avg(:), pcavecs(:,:)
        real,        allocatable :: class_diams(:), class_shifts(:,:)
        integer,     allocatable :: pinds(:)
        type(image)              :: cavg_raw
        type(diffmap_graph)      :: graph
        character(len=STDLEN)    :: recon_mode
        integer :: nptcls, npix, j, class_ldim(3), stkind, indstk, local_ind
        integer :: den_rank, icm_iters, k_nn_eff
        real    :: class_moldiam, class_mskdiam, class_mskrad, sdev_noise
        real    :: recon_rmse, recon_rel_rmse, icm_score, resid_ratio
        logical :: icm_converged, icm_more_iters
        call transform_ptcls(params, build, spproj, 'ptcl2D', cls_id, imgs, pinds, phflip=l_phflip, cavg=cavg_raw)
        if( .not. allocated(imgs) ) THROW_HARD('transform_ptcls returned no images; denoise_project')
        nptcls = size(imgs)
        if( nptcls < 3 ) THROW_HARD('at least 3 particles per class required; denoise_project')
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
        class_ldim          = cavg_raw%get_ldim()
        params_mask%box     = class_ldim(1)
        params_mask%smpd    = cavg_raw%get_smpd()
        params_mask%msk     = real(class_ldim(1) / 2) - COSMSKHALFWIDTH
        if( l_mskdiam_override )then
            class_mskdiam = params%mskdiam
        else
            call automask2D(params_mask, class_mask, params_mask%ngrow, nint(params_mask%winsz), params_mask%edge, &
                            class_diams, class_shifts, verbose=.false.)
            class_moldiam = params_mask%smpd * real(min(round2even(class_diams(1) / params_mask%smpd + &
                2. * COSMSKHALFWIDTH), class_ldim(1)))
            class_mskdiam = class_moldiam * MSK_EXP_FAC
        endif
        class_mskrad  = min(real(class_ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * class_mskdiam / params_mask%smpd)
        do j = 1, nptcls
            call imgs_ppca(j)%norm_noise(build%lmsk, sdev_noise)
            call imgs_ppca(j)%mask2D_softavg(class_mskrad)
        end do
        call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
        call build_cls_split_graph(params, spproj, pinds, pcavecs, graph)
        if( graph%n /= nptcls ) THROW_HARD('diffusion-map graph size mismatch; denoise_project')
        if( trim(graph%metric) /= 'euc' .or. trim(graph%steering) /= 'none' )then
            THROW_HARD('denoise_project expected a non-steerable Euclidean graph')
        endif
        k_nn_eff = graph%k_nn
        call estimate_diffmap_denoise_rank(params, graph, cls_id, nptcls, den_rank, icm_converged, icm_iters, icm_score)
        call graph_nystrom_residual_preimage(imgs, cavg_raw, graph, den_ptcls, den_rank)
        recon_mode = 'nystrom'
        if( .not. allocated(den_ptcls) ) THROW_HARD('diffusion-map denoising failed; denoise_project')
        call calc_diffmap_reconstruction_error(imgs, den_ptcls, recon_rmse, recon_rel_rmse)
        call calc_diffmap_residual_energy_ratio(imgs, den_ptcls, cavg_raw, resid_ratio)
        icm_more_iters = (.not. icm_converged) .and. icm_iters > 0
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,A,A,ES12.4,A,ES12.4,A,ES12.4,A,I8,A,A,A,A)') &
            '>>> CLASS class=', cls_id, ' index=', class_index, ' nptcls=', nptcls, ' features=', den_rank, &
            ' k_nn=', k_nn_eff, ' recon=', trim(recon_mode), &
            ' recon_rmse=', recon_rmse, ' recon_rel_rmse=', recon_rel_rmse, &
            ' resid_energy=', resid_ratio, ' icm_iters=', icm_iters, &
            ' icm_converged=', merge('yes','no ',icm_converged), ' icm_more_iters=', merge('yes','no ',icm_more_iters)
        call flush(logfhandle)
        do j = 1, nptcls
            if( processed(pinds(j)) ) THROW_HARD('particle was processed twice; denoise_project')
            call spproj%map_ptcl_ind2stk_ind('ptcl2D', pinds(j), stkind, indstk)
            if( l_dense_output )then
                if( l_index_by_pind )then
                    local_ind = pinds(j)
                else
                    local_ind = nprocessed + 1
                endif
                call write_diffmap_stack_image(raw_stks(1), local_ind, imgs(j))
                call write_diffmap_stack_image(den_stks(1), local_ind, den_ptcls(j))
                write(map_unit,'(I10,1X,I8,1X,I8,1X,I10)') pinds(j), stkind, indstk, local_ind
            else
                if( stkind < 1 .or. stkind > size(raw_stks) ) THROW_HARD('invalid output stack index; denoise_project')
                call write_diffmap_stack_image(raw_stks(stkind), indstk, imgs(j))
                call write_diffmap_stack_image(den_stks(stkind), indstk, den_ptcls(j))
            endif
            processed(pinds(j)) = .true.
            nprocessed = nprocessed + 1
        end do
        call cavg_raw%kill
        call graph%kill()
        if( allocated(imgs)         ) call dealloc_imgarr(imgs)
        if( allocated(imgs_ppca)    ) call dealloc_imgarr(imgs_ppca)
        if( allocated(class_mask)   ) call dealloc_imgarr(class_mask)
        if( allocated(den_ptcls)    ) call dealloc_imgarr(den_ptcls)
        if( allocated(pinds)        ) deallocate(pinds)
        if( allocated(avg)          ) deallocate(avg)
        if( allocated(pcavecs)      ) deallocate(pcavecs)
        if( allocated(class_diams)  ) deallocate(class_diams)
        if( allocated(class_shifts) ) deallocate(class_shifts)
    end subroutine write_diffmap_denoised_class

    subroutine setup_diffmap_so3_mixture_map(params, build, eulspace, neigh_map)
        type(parameters),         intent(in)    :: params
        type(builder),            intent(inout) :: build
        type(oris),               intent(inout) :: eulspace
        type(eulspace_neigh_map), intent(inout) :: neigh_map
        type(oris) :: eulspace_sub
        call eulspace%kill
        call neigh_map%kill
        call eulspace%new(params%nspace, is_ptcl=.false.)
        call build%pgrpsyms%build_refspiral(eulspace)
        call eulspace_sub%new(params%nspace_sub, is_ptcl=.false.)
        call build%pgrpsyms%build_refspiral(eulspace_sub)
        call neigh_map%new(eulspace, eulspace_sub, build%pgrpsyms)
        call eulspace_sub%kill
    end subroutine setup_diffmap_so3_mixture_map

    subroutine assign_global_so3_mixtures(build, spproj, eulspace, neigh_map, cls_inds, particle_mix, mix_counts)
        type(builder),            intent(inout) :: build
        type(sp_project),         intent(inout) :: spproj
        type(oris),               intent(in)    :: eulspace
        type(eulspace_neigh_map), intent(in)    :: neigh_map
        integer,                  intent(in)    :: cls_inds(:)
        integer, allocatable,     intent(out)   :: particle_mix(:), mix_counts(:)
        integer, allocatable :: full2sub(:)
        type(ori) :: o
        integer :: iptcl, full_proj, isub, nptcls, cls
        nptcls = spproj%os_ptcl2D%get_noris()
        full2sub = neigh_map%get_full2sub_map()
        allocate(particle_mix(nptcls), source=0)
        allocate(mix_counts(neigh_map%get_nsub()), source=0)
        do iptcl = 1, nptcls
            if( spproj%os_ptcl2D%get_state(iptcl) <= 0 ) cycle
            cls = spproj%os_ptcl2D%get_int(iptcl, 'class')
            if( .not. any(cls_inds == cls) ) cycle
            call spproj%os_ptcl3D%get_ori(iptcl, o)
            full_proj = build%pgrpsyms%find_closest_proj(eulspace, o)
            if( full_proj < 1 .or. full_proj > size(full2sub) )then
                THROW_HARD('SO3 closest projection outside mixture map; denoise_project')
            endif
            isub = full2sub(full_proj)
            if( isub < 1 .or. isub > size(mix_counts) )then
                THROW_HARD('SO3 mixture id outside subspace range; denoise_project')
            endif
            particle_mix(iptcl) = isub
            mix_counts(isub) = mix_counts(isub) + 1
            call o%kill
        end do
        if( allocated(full2sub) ) deallocate(full2sub)
    end subroutine assign_global_so3_mixtures

    subroutine pack_diffmap_so3_mixture_pinds(particle_mix, mix_counts, mix_rowptr, mix_pinds)
        integer,              intent(in)  :: particle_mix(:), mix_counts(:)
        integer, allocatable, intent(out) :: mix_rowptr(:), mix_pinds(:)
        integer, allocatable :: cursor(:)
        integer :: iptcl, isub, pos, nsub, total_count
        nsub = size(mix_counts)
        allocate(mix_rowptr(nsub + 1), source=1)
        do isub = 1, nsub
            mix_rowptr(isub + 1) = mix_rowptr(isub) + mix_counts(isub)
        end do
        total_count = mix_rowptr(nsub + 1) - 1
        allocate(mix_pinds(total_count), source=0)
        allocate(cursor(nsub), source=mix_rowptr(1:nsub))
        do iptcl = 1, size(particle_mix)
            isub = particle_mix(iptcl)
            if( isub == 0 ) cycle
            if( isub < 1 .or. isub > nsub ) THROW_HARD('missing global SO3 mixture assignment')
            pos = cursor(isub)
            if( pos < mix_rowptr(isub) .or. pos >= mix_rowptr(isub + 1) )then
                THROW_HARD('SO3 mixture particle list overflow; denoise_project')
            endif
            mix_pinds(pos) = iptcl
            cursor(isub) = pos + 1
        end do
        do isub = 1, nsub
            if( cursor(isub) /= mix_rowptr(isub + 1) )then
                THROW_HARD('SO3 mixture particle list population mismatch; denoise_project')
            endif
        end do
        deallocate(cursor)
    end subroutine pack_diffmap_so3_mixture_pinds

    subroutine build_diffmap_so3_mixture_graphs(params, spproj, mix_counts, mix_rowptr, mix_pinds)
        use simple_diff_map_graphs, only: diffmap_graph, build_cls_split_graph
        type(parameters), intent(inout) :: params
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: mix_counts(:), mix_rowptr(:), mix_pinds(:)
        type(diffmap_graph) :: graph
        integer :: isub, nptcls, fromp, top
        if( size(mix_rowptr) /= size(mix_counts) + 1 )then
            THROW_HARD('SO3 mixture row pointer/count size mismatch; denoise_project')
        endif
        do isub = 1, size(mix_counts)
            nptcls = mix_counts(isub)
            if( nptcls < 2 ) cycle
            fromp = mix_rowptr(isub)
            top   = mix_rowptr(isub + 1) - 1
            if( fromp < 1 .or. top > size(mix_pinds) .or. top - fromp + 1 /= nptcls )then
                THROW_HARD('SO3 mixture particle list range mismatch; denoise_project')
            endif
            call build_cls_split_graph(params=params, spproj=spproj, pinds=mix_pinds(fromp:top), graph=graph)
            if( graph%n /= nptcls ) THROW_HARD('SO3 mixture graph size mismatch; denoise_project')
            if( trim(graph%metric) /= 'ori' .or. trim(graph%steering) /= 'none' )then
                THROW_HARD('denoise_project expected a non-steerable orientation graph')
            endif
            write(logfhandle,'(A,I8,A,I10,A,A,A,I8,A,I10)') &
                '>>> SO3_MIXTURE mixture=', isub, ' nptcls=', nptcls, ' graph=', trim(graph%metric), &
                ' k_nn=', graph%k_nn, ' nnz=', graph%nnz
            call flush(logfhandle)
            call graph%kill()
        end do
    end subroutine build_diffmap_so3_mixture_graphs

    subroutine run_diffmap_so3_mixture_graphs(params, build, spproj, cls_inds)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: cls_inds(:)
        type(oris) :: eulspace
        type(eulspace_neigh_map) :: neigh_map
        integer, allocatable :: particle_mix(:), mix_counts(:), mix_rowptr(:), mix_pinds(:)
        call setup_diffmap_so3_mixture_map(params, build, eulspace, neigh_map)
        call assign_global_so3_mixtures(build, spproj, eulspace, neigh_map, cls_inds, particle_mix, mix_counts)
        call pack_diffmap_so3_mixture_pinds(particle_mix, mix_counts, mix_rowptr, mix_pinds)
        call build_diffmap_so3_mixture_graphs(params, spproj, mix_counts, mix_rowptr, mix_pinds)
        if( allocated(particle_mix) ) deallocate(particle_mix)
        if( allocated(mix_counts) ) deallocate(mix_counts)
        if( allocated(mix_rowptr) ) deallocate(mix_rowptr)
        if( allocated(mix_pinds) ) deallocate(mix_pinds)
        call neigh_map%kill
        call eulspace%kill
    end subroutine run_diffmap_so3_mixture_graphs

    subroutine estimate_diffmap_denoise_rank(params, graph, cls_id, nptcls, den_rank, icm_converged, icm_iters, icm_score)
        use simple_diff_map_graphs, only: diffmap_graph
        use simple_diffusion_maps,  only: embed_graph
        type(parameters),    intent(in)  :: params
        type(diffmap_graph), intent(in)  :: graph
        integer,             intent(in)  :: cls_id, nptcls
        integer,             intent(out) :: den_rank, icm_iters
        logical,             intent(out) :: icm_converged
        real,                intent(out) :: icm_score
        real, allocatable :: coords(:,:), eigvals(:)
        integer :: rank_scan, max_rank
        max_rank = max(1, nptcls - 2)
        if( params%neigs > 0 )then
            rank_scan = min(max(1, params%neigs), max_rank)
        else
            rank_scan = diffmap_denoise_auto_neigs_scan(nptcls)
        endif
        rank_scan = min(max(1, rank_scan), max_rank)
        call embed_graph(graph, rank_scan, coords, eigvals)
        if( allocated(eigvals) .and. size(eigvals) > 0 )then
            call select_diffmap_denoise_rank_icm(eigvals, size(eigvals), cls_id, den_rank, &
                                                 icm_converged, icm_iters, icm_score)
        else
            den_rank      = rank_scan
            icm_converged = .false.
            icm_iters     = 0
            icm_score     = huge(icm_score)
            write(logfhandle,'(A,I8,A,I8)') 'Diffmap denoise rank warning: no eigenspectrum; class=', cls_id, &
                ' fallback_features=', den_rank
            call flush(logfhandle)
        endif
        den_rank = min(max(1, den_rank), nptcls)
        if( allocated(coords) ) deallocate(coords)
        if( allocated(eigvals) ) deallocate(eigvals)
    end subroutine estimate_diffmap_denoise_rank

    integer function diffmap_denoise_auto_neigs_scan(nptcls) result(neigs_scan)
        integer, intent(in) :: nptcls
        neigs_scan = min(50, max(1, nptcls - 2))
        if( nptcls > 3 ) neigs_scan = max(2, neigs_scan)
    end function diffmap_denoise_auto_neigs_scan

    subroutine select_diffmap_denoise_rank_icm(eigvals, max_neigs, cls_id, nkeep, converged, niter, score)
        real,    intent(in)  :: eigvals(:)
        integer, intent(in)  :: max_neigs, cls_id
        integer, intent(out) :: nkeep, niter
        logical, intent(out) :: converged
        real,    intent(out) :: score
        real, allocatable :: spec(:)
        real :: smin, smax, delta, beta, complexity, trial_score, best_score
        integer :: i, n, nmin_rank, upper_rank, lower_rank, best_seed, k_trial, trial_iter
        integer :: seeds(2)
        character(len=12) :: seed_names(2)
        logical :: trial_converged, best_converged
        n = min(size(eigvals), max(1, max_neigs))
        nmin_rank = diffmap_denoise_min_neigs(n)
        if( n <= 0 )then
            nkeep     = 1
            niter     = 0
            score     = 0.
            converged = .false.
            return
        endif
        allocate(spec(n), source=0.)
        do i = 1, n
            if( ieee_is_finite(eigvals(i)) .and. eigvals(i) > real(DTINY) )then
                spec(i) = log(eigvals(i))
            else
                spec(i) = log(real(DTINY))
            endif
        end do
        smin  = minval(spec)
        smax  = maxval(spec)
        delta = smax - smin
        if( .not. ieee_is_finite(delta) .or. delta <= 1.e-6 )then
            nkeep     = nmin_rank
            niter     = 0
            score     = 0.
            converged = .true.
            deallocate(spec)
            return
        endif
        spec = (spec - smin) / delta
        upper_rank = diffmap_denoise_initial_neigs_from_gap(spec, nmin_rank)
        upper_rank = max(nmin_rank, min(n, upper_rank))
        lower_rank = nint(DIFFMAP_DENOISE_ICM_RANK_LOWER_SEED_FRAC * real(upper_rank))
        lower_rank = max(nmin_rank, min(upper_rank, lower_rank))
        seeds(1) = upper_rank
        seeds(2) = lower_rank
        seed_names = [character(len=12) :: 'upper_gap', 'fraction']
        beta       = DIFFMAP_DENOISE_ICM_RANK_BETA_FRAC
        complexity = DIFFMAP_DENOISE_ICM_RANK_COMPLEXITY_FRAC
        best_score     = huge(best_score)
        nkeep          = upper_rank
        niter          = 0
        best_seed      = 1
        best_converged = .false.
        do i = 1, size(seeds)
            call run_diffmap_denoise_icm_rank_seed(spec, cls_id, trim(seed_names(i)), seeds(i), nmin_rank, upper_rank, &
                beta, complexity, k_trial, trial_score, trial_iter, trial_converged)
            if( trial_score < best_score )then
                best_score     = trial_score
                nkeep          = k_trial
                niter          = trial_iter
                best_seed      = i
                best_converged = trial_converged
            endif
        end do
        score     = best_score
        converged = best_converged
        deallocate(spec)
    end subroutine select_diffmap_denoise_rank_icm

    subroutine run_diffmap_denoise_icm_rank_seed(spec, cls_id, seed_name, seed_rank, nmin_rank, nmax_rank, beta, alpha, &
                                                 nkeep, score, niter, converged)
        real,             intent(in)  :: spec(:), beta, alpha
        character(len=*), intent(in)  :: seed_name
        integer,          intent(in)  :: cls_id, seed_rank, nmin_rank, nmax_rank
        integer,          intent(out) :: nkeep, niter
        real,             intent(out) :: score
        logical,          intent(out) :: converged
        integer, allocatable :: labels(:), prev_labels(:)
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: iter, i, n, nchanged, maxits
        n = size(spec)
        allocate(labels(n), prev_labels(n), source=0)
        labels = 0
        labels(1:max(nmin_rank, min(nmax_rank, seed_rank))) = 1
        if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
        converged = .false.
        niter     = 0
        maxits    = max(DIFFMAP_DENOISE_ICM_RANK_MAXITS, n)
        do iter = 1, maxits
            prev_labels = labels
            call estimate_diffmap_denoise_icm_rank_stats(spec, prev_labels, mu_drop, mu_keep, var_drop, var_keep)
            do i = 1, nmax_rank
                call update_diffmap_denoise_icm_rank_site(spec, prev_labels, i, beta, alpha, nmin_rank, nmax_rank, &
                    mu_drop, mu_keep, var_drop, var_keep, labels(i))
            end do
            labels(1:nmin_rank) = 1
            if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
            nchanged = count(labels /= prev_labels)
            score = score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
            niter = iter
            if( nchanged == 0 )then
                converged = .true.
                exit
            endif
        end do
        nkeep = diffmap_denoise_rank_prefix(labels, nmin_rank, nmax_rank)
        score = score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
        deallocate(labels, prev_labels)
    end subroutine run_diffmap_denoise_icm_rank_seed

    integer function diffmap_denoise_rank_prefix(labels, nmin_rank, nmax_rank) result(nkeep)
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        integer :: i
        nkeep = 0
        do i = 1, min(size(labels), nmax_rank)
            if( labels(i) == 1 ) nkeep = i
        end do
        nkeep = max(nmin_rank, min(nmax_rank, nkeep))
    end function diffmap_denoise_rank_prefix

    integer function diffmap_denoise_min_neigs(max_neigs) result(nmin_rank)
        integer, intent(in) :: max_neigs
        nmin_rank = 1
        if( max_neigs >= 2 ) nmin_rank = 2
    end function diffmap_denoise_min_neigs

    integer function diffmap_denoise_initial_neigs_from_gap(spec, nmin_rank) result(nkeep)
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
    end function diffmap_denoise_initial_neigs_from_gap

    subroutine estimate_diffmap_denoise_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
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
    end subroutine estimate_diffmap_denoise_icm_rank_stats

    subroutine update_diffmap_denoise_icm_rank_site(spec, labels, ind, beta, alpha, nmin_rank, nmax_rank, &
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
    end subroutine update_diffmap_denoise_icm_rank_site

    real function score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank) result(score)
        real,    intent(in) :: spec(:), beta, alpha
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: i, kfree
        call estimate_diffmap_denoise_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
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
    end function score_diffmap_denoise_icm_rank_solution

    subroutine graph_nystrom_residual_preimage(raw_imgs, avg_img, graph, den_imgs, rank_keep)
        type(image),       intent(inout) :: raw_imgs(:)
        type(image),       intent(inout) :: avg_img
        type(diffmap_graph), intent(in)  :: graph
        type(image), allocatable, intent(out) :: den_imgs(:)
        integer,            intent(in)  :: rank_keep
        type(image), allocatable :: mode_imgs(:)
        type(image) :: resid_img
        real, allocatable :: evals(:), evecs(:,:), phi_ext(:,:)
        real :: coeff, op_w
        real :: smpd
        integer :: i, j, p, k, n, ldim(3), rank_used, nev, eig_idx, eig_info, max_basis
        n = size(raw_imgs)
        if( graph%n /= n ) THROW_HARD('nystrom preimage graph/image count mismatch; denoise_project')
        if( n < 3 ) THROW_HARD('nystrom preimage requires at least three graph nodes; denoise_project')
        ldim = avg_img%get_ldim()
        smpd = avg_img%get_smpd()
        rank_used = min(max(1, rank_keep), max(1, n - 2))
        nev = rank_used + 1
        allocate(evals(nev), evecs(n,nev), phi_ext(n,rank_used), source=0.)
        max_basis = min(n, max(160, 8 * nev + 80))
        call sparse_eigh(graph_matvec, graph, n, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=eig_info)
        if( eig_info /= 0 ) THROW_HARD('sparse eigensolve failed in nystrom preimage; denoise_project')
        do k = 1, rank_used
            eig_idx = nev - k
            if( abs(evals(eig_idx)) <= real(DTINY) ) cycle
            do i = 1, n
                coeff = 0.
                do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                    j = graph%colind(p)
                    if( j < 1 .or. j > n ) cycle
                    if( allocated(graph%wnorm) )then
                        op_w = graph%wnorm(p)
                    else
                        op_w = graph%w(p)
                    endif
                    coeff = coeff + op_w * evecs(j,eig_idx)
                end do
                phi_ext(i,k) = coeff / evals(eig_idx)
            end do
        end do
        allocate(den_imgs(n))
        allocate(mode_imgs(rank_used))
        do i = 1, n
            call den_imgs(i)%new(ldim, smpd, wthreads=.false.)
        end do
        do k = 1, rank_used
            call mode_imgs(k)%new(ldim, smpd, wthreads=.false.)
            call mode_imgs(k)%zero()
        end do
        call resid_img%new(ldim, smpd, wthreads=.false.)
        do k = 1, rank_used
            eig_idx = nev - k
            do i = 1, n
                coeff = evecs(i,eig_idx)
                if( abs(coeff) <= real(DTINY) ) cycle
                call resid_img%copy_fast(raw_imgs(i))
                call resid_img%subtr(avg_img)
                call mode_imgs(k)%add(resid_img, coeff)
            end do
        end do
        !$omp parallel do default(shared) private(i,k,coeff) schedule(dynamic) proc_bind(close)
        do i = 1, n
            call den_imgs(i)%copy_fast(avg_img)
            do k = 1, rank_used
                coeff = phi_ext(i,k)
                if( abs(coeff) <= real(DTINY) ) cycle
                call den_imgs(i)%add(mode_imgs(k), coeff)
            end do
        end do
        !$omp end parallel do
        write(logfhandle,'(A,I8,A,I8,A,I8,A,F8.4)') 'Diffmap Nyström preimage: n=', n, ' rank=', rank_used, &
            ' nnz=', graph%nnz, ' lambda_min=', minval(evals(max(1,nev-rank_used):nev-1))
        call flush(logfhandle)
        do k = 1, rank_used
            if( mode_imgs(k)%exists() ) call mode_imgs(k)%kill
        end do
        if( resid_img%exists() ) call resid_img%kill
        deallocate(mode_imgs, evals, evecs, phi_ext)
    end subroutine graph_nystrom_residual_preimage

    subroutine calc_diffmap_residual_energy_ratio(raw_imgs, den_imgs, avg_img, ratio)
        type(image), intent(inout) :: raw_imgs(:), den_imgs(:), avg_img
        real,        intent(out)   :: ratio
        real, allocatable :: raw_vec(:), den_vec(:), avg_vec(:)
        real(kind=8) :: raw_ssq, den_ssq
        integer :: i
        if( size(raw_imgs) /= size(den_imgs) ) THROW_HARD('residual energy image count mismatch; denoise_project')
        avg_vec = avg_img%serialize()
        raw_ssq = 0.d0
        den_ssq = 0.d0
        !$omp parallel do default(shared) private(i,raw_vec,den_vec) schedule(static) proc_bind(close) &
        !$omp& reduction(+:raw_ssq,den_ssq)
        do i = 1, size(raw_imgs)
            raw_vec = raw_imgs(i)%serialize()
            den_vec = den_imgs(i)%serialize()
            if( size(raw_vec) /= size(avg_vec) .or. size(den_vec) /= size(avg_vec) )then
                THROW_HARD('residual energy image size mismatch; denoise_project')
            endif
            raw_ssq = raw_ssq + sum(real(raw_vec - avg_vec, kind=8) * real(raw_vec - avg_vec, kind=8))
            den_ssq = den_ssq + sum(real(den_vec - avg_vec, kind=8) * real(den_vec - avg_vec, kind=8))
            deallocate(raw_vec, den_vec)
        end do
        !$omp end parallel do
        ratio = real(sqrt(den_ssq / max(raw_ssq, real(DTINY, kind=8))))
        deallocate(avg_vec)
    end subroutine calc_diffmap_residual_energy_ratio

    subroutine calc_diffmap_reconstruction_error(raw_imgs, den_imgs, rmse, rel_rmse)
        type(image), intent(inout) :: raw_imgs(:), den_imgs(:)
        real,        intent(out)   :: rmse, rel_rmse
        real, allocatable :: raw_vec(:), den_vec(:)
        real(kind=8) :: ssq, raw_ssq
        integer(kind=8) :: npix_total
        integer :: i
        if( size(raw_imgs) /= size(den_imgs) ) THROW_HARD('reconstruction error image count mismatch; denoise_project')
        ssq        = 0.d0
        raw_ssq    = 0.d0
        npix_total = 0_8
        !$omp parallel do default(shared) private(i,raw_vec,den_vec) schedule(static) proc_bind(close) &
        !$omp& reduction(+:ssq,raw_ssq,npix_total)
        do i = 1, size(raw_imgs)
            raw_vec = raw_imgs(i)%serialize()
            den_vec = den_imgs(i)%serialize()
            if( size(raw_vec) /= size(den_vec) ) THROW_HARD('reconstruction error image size mismatch; denoise_project')
            ssq        = ssq + sum(real(raw_vec - den_vec, kind=8) * real(raw_vec - den_vec, kind=8))
            raw_ssq    = raw_ssq + sum(real(raw_vec, kind=8) * real(raw_vec, kind=8))
            npix_total = npix_total + int(size(raw_vec), kind=8)
            deallocate(raw_vec, den_vec)
        end do
        !$omp end parallel do
        if( npix_total > 0_8 )then
            rmse = real(sqrt(ssq / real(npix_total, kind=8)))
        else
            rmse = 0.
        endif
        rel_rmse = real(sqrt(ssq / max(raw_ssq, real(DTINY, kind=8))))
    end subroutine calc_diffmap_reconstruction_error

    subroutine write_diffmap_stack_image(fname, indstk, img)
        type(string), intent(in)    :: fname
        integer,      intent(in)    :: indstk
        type(image),  intent(inout) :: img
        type(imgfile) :: ioimg
        real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
        integer :: ldim(3)
        if( indstk < 1 ) THROW_HARD('invalid stack-local output index; denoise_project')
        if( .not. img%is_2d() ) THROW_HARD('2D particle image expected; denoise_project')
        ldim = img%get_ldim()
        ldim(3) = 1
        call ioimg%open(fname, ldim, img%get_smpd(), formatchar='M', readhead=file_exists(fname))
        call img%get_rmat_ptr(rmat_ptr)
        call ioimg%wmrcSlices(indstk, indstk, rmat_ptr, ldim, img%is_ft())
        call ioimg%close
        rmat_ptr => null()
    end subroutine write_diffmap_stack_image

    subroutine write_diffmap_project(params, spproj, raw_stk, den_stk, l_ctf_no, outproj)
        type(parameters), intent(in)    :: params
        type(sp_project), intent(inout) :: spproj
        type(string),     intent(in)    :: raw_stk, den_stk
        logical,          intent(in)    :: l_ctf_no
        type(sp_project), intent(inout) :: outproj
        integer :: iptcl, nptcls, nraw, nden, ldim_raw(3), ldim_den(3)
        nptcls = spproj%os_ptcl2D%get_noris()
        call find_ldim_nptcls(raw_stk, ldim_raw, nraw)
        call find_ldim_nptcls(den_stk, ldim_den, nden)
        if( nraw /= nptcls .or. nden /= nptcls ) THROW_HARD('compact output stack count mismatch; denoise_project')
        if( .not. all(ldim_raw(1:2) == ldim_den(1:2)) )then
            THROW_HARD('raw/denoised output stack dimensions differ; denoise_project')
        endif
        call outproj%copy(spproj)
        call outproj%update_projinfo(params%projfile)
        if( outproj%os_stk%get_noris() > 0 ) call outproj%os_stk%kill()
        call outproj%os_stk%new(1, is_ptcl=.false.)
        call outproj%os_stk%set(1, 'stk',        simple_abspath(raw_stk))
        call outproj%os_stk%set(1, 'stk_den',    simple_abspath(den_stk))
        call outproj%os_stk%set(1, 'box',        ldim_raw(1))
        call outproj%os_stk%set(1, 'nptcls',     nraw)
        call outproj%os_stk%set(1, 'nptcls_stk', nraw)
        call outproj%os_stk%set(1, 'fromp',      1)
        call outproj%os_stk%set(1, 'top',        nraw)
        call outproj%os_stk%set(1, 'stkkind',    'single')
        call outproj%os_stk%set(1, 'imgkind',    'ptcl')
        call outproj%os_stk%set(1, 'smpd',       params%smpd_crop)
        call outproj%os_stk%set(1, 'kv',         params%kv)
        call outproj%os_stk%set(1, 'cs',         params%cs)
        call outproj%os_stk%set(1, 'fraca',      params%fraca)
        if( l_ctf_no )then
            call outproj%os_stk%set(1, 'ctf',    'no')
        else
            call outproj%os_stk%set(1, 'ctf',    'flip')
        endif
        call outproj%os_stk%set(1, 'phaseplate', 'no')
        call outproj%os_stk%set(1, 'state',      1)
        do iptcl = 1, nptcls
            call outproj%os_ptcl2D%set(iptcl, 'stkind', 1)
            call outproj%os_ptcl3D%set(iptcl, 'stkind', 1)
            call outproj%os_ptcl2D%set(iptcl, 'indstk', iptcl)
            call outproj%os_ptcl3D%set(iptcl, 'indstk', iptcl)
        end do
        call outproj%os_ptcl2D%zero_inpl()
        call outproj%os_ptcl3D%zero_inpl()
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call preserve_diffmap_frc2D(spproj, outproj)
        call outproj%write(params%projfile)
    end subroutine write_diffmap_project

    subroutine write_diffmap_partial_project(params, spproj, nparts_run, part_counts, row_part, row_local, l_ctf_no, outproj)
        type(parameters), intent(in)    :: params
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: nparts_run
        integer,          intent(in)    :: part_counts(nparts_run), row_part(:), row_local(:)
        logical,          intent(in)    :: l_ctf_no
        type(sp_project), intent(inout) :: outproj
        type(string) :: raw_part, den_part
        integer :: ipart, iptcl, nptcls, nraw, nden, ldim_raw(3), ldim_den(3), fromp, top, stkind, indstk
        nptcls = spproj%os_ptcl2D%get_noris()
        call outproj%copy(spproj)
        call outproj%update_projinfo(params%projfile)
        if( outproj%os_stk%get_noris() > 0 ) call outproj%os_stk%kill()
        call outproj%os_stk%new(nparts_run, is_ptcl=.false.)
        fromp = 1
        do ipart = 1, nparts_run
            raw_part = string('denoise_raw_particles_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            den_part = string('denoise_particles_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            call find_ldim_nptcls(raw_part, ldim_raw, nraw)
            call find_ldim_nptcls(den_part, ldim_den, nden)
            if( nraw /= part_counts(ipart) .or. nden /= part_counts(ipart) )then
                THROW_HARD('denoise partial project stack count mismatch')
            endif
            if( .not. all(ldim_raw(1:2) == ldim_den(1:2)) )then
                THROW_HARD('denoise partial project stack dimensions differ')
            endif
            top = fromp + part_counts(ipart) - 1
            call outproj%os_stk%set(ipart, 'stk',        simple_abspath(raw_part))
            call outproj%os_stk%set(ipart, 'stk_den',    simple_abspath(den_part))
            call outproj%os_stk%set(ipart, 'box',        ldim_raw(1))
            call outproj%os_stk%set(ipart, 'nptcls',     part_counts(ipart))
            call outproj%os_stk%set(ipart, 'nptcls_stk', part_counts(ipart))
            call outproj%os_stk%set(ipart, 'fromp',      fromp)
            call outproj%os_stk%set(ipart, 'top',        top)
            call outproj%os_stk%set(ipart, 'stkkind',    'single')
            call outproj%os_stk%set(ipart, 'imgkind',    'ptcl')
            call outproj%os_stk%set(ipart, 'smpd',       params%smpd_crop)
            call outproj%os_stk%set(ipart, 'kv',         params%kv)
            call outproj%os_stk%set(ipart, 'cs',         params%cs)
            call outproj%os_stk%set(ipart, 'fraca',      params%fraca)
            if( l_ctf_no )then
                call outproj%os_stk%set(ipart, 'ctf',    'no')
            else
                call outproj%os_stk%set(ipart, 'ctf',    'flip')
            endif
            call outproj%os_stk%set(ipart, 'phaseplate', 'no')
            call outproj%os_stk%set(ipart, 'state',      1)
            fromp = top + 1
            call raw_part%kill
            call den_part%kill
        end do
        do iptcl = 1, nptcls
            if( spproj%os_ptcl2D%get_state(iptcl) > 0 )then
                if( row_part(iptcl) < 1 .or. row_part(iptcl) > nparts_run )then
                    THROW_HARD('invalid particle part while writing denoise partial project')
                endif
                if( row_local(iptcl) < 1 .or. row_local(iptcl) > part_counts(row_part(iptcl)) )then
                    THROW_HARD('invalid particle local stack index while writing denoise partial project')
                endif
                call outproj%os_ptcl2D%set(iptcl, 'stkind', row_part(iptcl))
                call outproj%os_ptcl3D%set(iptcl, 'stkind', row_part(iptcl))
                call outproj%os_ptcl2D%set(iptcl, 'indstk', row_local(iptcl))
                call outproj%os_ptcl3D%set(iptcl, 'indstk', row_local(iptcl))
            else
                call outproj%os_ptcl2D%set(iptcl, 'stkind', 1)
                call outproj%os_ptcl3D%set(iptcl, 'stkind', 1)
                call outproj%os_ptcl2D%set(iptcl, 'indstk', 1)
                call outproj%os_ptcl3D%set(iptcl, 'indstk', 1)
            endif
        end do
        do iptcl = 1, nptcls
            if( spproj%os_ptcl2D%get_state(iptcl) <= 0 ) cycle
            call outproj%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, indstk)
            if( stkind /= row_part(iptcl) .or. indstk /= row_local(iptcl) )then
                THROW_HARD('ptcl2D stack mapping mismatch in denoise partial project')
            endif
            call outproj%map_ptcl_ind2stk_ind('ptcl3D', iptcl, stkind, indstk)
            if( stkind /= row_part(iptcl) .or. indstk /= row_local(iptcl) )then
                THROW_HARD('ptcl3D stack mapping mismatch in denoise partial project')
            endif
        end do
        call outproj%os_ptcl2D%zero_inpl()
        call outproj%os_ptcl3D%zero_inpl()
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call preserve_diffmap_frc2D(spproj, outproj)
        call outproj%write(params%projfile)
        write(logfhandle,'(A,I8,A,I10)') 'Denoise partial project stacks: ', nparts_run, ' particles=', nptcls
        call flush(logfhandle)
    end subroutine write_diffmap_partial_project

    subroutine preserve_diffmap_frc2D(spproj, outproj)
        type(sp_project), intent(in)    :: spproj
        type(sp_project), intent(inout) :: outproj
        type(string) :: frcs_fname
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( file_exists(frcs_fname%to_char()) ) call outproj%add_frcs2os_out(frcs_fname, 'frc2D')
        if( frcs_fname%is_allocated() ) call frcs_fname%kill
    end subroutine preserve_diffmap_frc2D

    subroutine finalize_diffmap_worker_outputs(params, nparts_run)
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts_run
        type(sp_project)     :: spproj, outproj
        type(string)         :: raw_part, den_part
        integer, allocatable :: row_part(:), row_local(:), row_stkind(:), row_indstk(:)
        integer, allocatable :: stk_offsets(:), slot_pind(:), stk_nptcls(:)
        integer, allocatable :: part_counts(:), part_nimgs(:), part_seen(:)
        integer, allocatable :: cls_inds(:), cls_pops(:)
        integer :: nptcls, nptcls_active, nstks, istk, pind, ipart, total_count
        integer :: ldim_raw(3), ldim_den(3), nraw, nden
        integer :: nptcls_slots, slot, local_slot
        logical :: l_phflip, l_ctf_no
        call spproj%read(params%projfile)
        call validate_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip, l_ctf_no)
        nptcls = spproj%os_ptcl2D%get_noris()
        nptcls_active = get_n_active_ptcl2D(spproj)
        nstks  = spproj%os_stk%get_noris()
        allocate(row_part(nptcls), row_local(nptcls), row_stkind(nptcls), row_indstk(nptcls), source=0)
        allocate(stk_offsets(nstks), stk_nptcls(nstks), source=0)
        allocate(part_counts(nparts_run), part_nimgs(nparts_run), source=0)
        nptcls_slots = 0
        do istk = 1, nstks
            stk_offsets(istk) = nptcls_slots
            stk_nptcls(istk) = get_diffmap_stack_nptcls(spproj, istk)
            nptcls_slots = nptcls_slots + stk_nptcls(istk)
        end do
        allocate(slot_pind(nptcls_slots), source=0)
        call read_diffmap_worker_maps(params, nparts_run, row_part, row_local, row_stkind, row_indstk, total_count)
        if( total_count /= nptcls_active )then
            write(logfhandle,'(A,I10,A,I10)') 'Denoise merge: mapped particles=', total_count, ' expected_active=', nptcls_active
            THROW_HARD('distributed denoise_project did not produce exactly one row per active particle')
        endif
        do pind = 1, nptcls
            if( spproj%os_ptcl2D%get_state(pind) <= 0 ) cycle
            if( row_part(pind) < 1 .or. row_local(pind) < 1 ) THROW_HARD('missing particle in denoise worker maps')
            if( row_stkind(pind) < 1 .or. row_stkind(pind) > nstks ) THROW_HARD('invalid stack index in denoise worker maps')
            if( row_indstk(pind) < 1 .or. row_indstk(pind) > stk_nptcls(row_stkind(pind)) )then
                THROW_HARD('invalid stack-local index in denoise worker maps')
            endif
            slot = stk_offsets(row_stkind(pind)) + row_indstk(pind)
            if( slot < 1 .or. slot > nptcls_slots ) THROW_HARD('invalid stack slot in denoise worker maps')
            if( slot_pind(slot) /= 0 ) THROW_HARD('duplicate stack slot in denoise worker maps')
            slot_pind(slot) = pind
            if( row_part(pind) > nparts_run ) THROW_HARD('invalid part index in denoise worker maps')
            part_counts(row_part(pind)) = part_counts(row_part(pind)) + 1
        end do
        write(logfhandle,'(A,I8,A,I10)') 'Denoise finalize: partial stacks=', nparts_run, ' particles=', nptcls
        call flush(logfhandle)
        do ipart = 1, nparts_run
            if( part_counts(ipart) < 1 ) THROW_HARD('empty denoise worker part')
            raw_part = string('denoise_raw_particles_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            den_part = string('denoise_particles_part')//int2str_pad(ipart, params%numlen)//'.mrcs'
            if( .not. file_exists(raw_part%to_char()) ) THROW_HARD('missing denoise worker raw stack')
            if( .not. file_exists(den_part%to_char()) ) THROW_HARD('missing denoise worker denoised stack')
            call find_ldim_nptcls(raw_part, ldim_raw, nraw)
            call find_ldim_nptcls(den_part, ldim_den, nden)
            if( nraw /= part_counts(ipart) .or. nden /= part_counts(ipart) )then
                write(logfhandle,'(A,I8,A,I10,A,I10,A,I10)') 'Denoise part count mismatch: part=', ipart, &
                    ' map_count=', part_counts(ipart), ' raw=', nraw, ' denoised=', nden
                THROW_HARD('denoise partial stack count mismatch')
            endif
            if( .not. all(ldim_raw(1:2) == ldim_den(1:2)) )then
                THROW_HARD('raw/denoised partial stack dimensions differ; denoise_project')
            endif
            if( ldim_raw(1) /= params%box_crop .or. ldim_raw(2) /= params%box_crop )then
                THROW_HARD('unexpected denoise partial stack dimensions')
            endif
            part_nimgs(ipart) = nraw
            call raw_part%kill
            call den_part%kill
        end do
        allocate(part_seen(maxval(part_nimgs)), source=0)
        do ipart = 1, nparts_run
            part_seen = 0
            do pind = 1, nptcls
                if( spproj%os_ptcl2D%get_state(pind) <= 0 ) cycle
                if( row_part(pind) /= ipart ) cycle
                local_slot = row_local(pind)
                if( local_slot < 1 .or. local_slot > part_nimgs(ipart) )then
                    THROW_HARD('particle local index out of partial stack range; denoise_project')
                endif
                if( part_seen(local_slot) /= 0 ) THROW_HARD('duplicate particle local index in partial stack')
                part_seen(local_slot) = pind
            end do
            if( count(part_seen(1:part_nimgs(ipart)) > 0) /= part_nimgs(ipart) )then
                THROW_HARD('partial stack local index coverage is incomplete; denoise_project')
            endif
        end do
        call write_diffmap_partial_project(params, spproj, nparts_run, part_counts, row_part, row_local, l_ctf_no, outproj)
        write(logfhandle,'(A,A)') 'Denoise project written: ', params%projfile%to_char()
        call flush(logfhandle)
        call outproj%kill
        call spproj%kill
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
        deallocate(row_part, row_local, row_stkind, row_indstk, stk_offsets, slot_pind, stk_nptcls)
        deallocate(part_counts, part_nimgs, part_seen)
        if( raw_part%is_allocated() ) call raw_part%kill
        if( den_part%is_allocated() ) call den_part%kill
    end subroutine finalize_diffmap_worker_outputs

    subroutine read_diffmap_worker_maps(params, nparts_run, row_part, row_local, row_stkind, row_indstk, total_count)
        type(parameters), intent(in)    :: params
        integer,          intent(in)    :: nparts_run
        integer,          intent(inout) :: row_part(:), row_local(:), row_stkind(:), row_indstk(:)
        integer,          intent(out)   :: total_count
        type(string) :: map_fname
        integer :: ipart, funit, ios, pind, stkind, indstk, local_ind, nptcls
        character(len=XLONGSTRLEN) :: line
        nptcls = size(row_part)
        total_count = 0
        do ipart = 1, nparts_run
            map_fname = string('denoise_particle_map_part')//int2str_pad(ipart, params%numlen)//TXT_EXT
            open(newunit=funit, file=map_fname%to_char(), status='old', action='read', iostat=ios)
            call fileiochk('read_diffmap_worker_maps opening '//map_fname%to_char(), ios)
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                read(line,*) pind, stkind, indstk, local_ind
                if( pind < 1 .or. pind > nptcls ) THROW_HARD('invalid particle index in denoise worker map')
                if( row_part(pind) /= 0 ) THROW_HARD('duplicate particle index in denoise worker maps')
                row_part(pind)   = ipart
                row_local(pind)  = local_ind
                row_stkind(pind) = stkind
                row_indstk(pind) = indstk
                total_count = total_count + 1
            end do
            close(funit)
            call map_fname%kill
        end do
    end subroutine read_diffmap_worker_maps

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

end module simple_denoise_project_strategy
