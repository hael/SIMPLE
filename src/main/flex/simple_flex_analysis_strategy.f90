!@descr: shmem/worker/master sparse diffusion-map 3D variability analysis
module simple_flex_analysis_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,               only: builder
use simple_cmdline,               only: cmdline
use simple_diff_map_denoise,      only: select_spectral_rank_icm
use simple_diff_map_graphs,       only: diffmap_graph, build_gated_euclidean_knn_graph, &
    &find_gated_euclidean_neighbors_rows, build_gated_euclidean_graph_from_neighbors
use simple_diffusion_maps,        only: embed_graph
use simple_flex_diffmap_features, only: prepare_flex_diffmap_features, prepare_flex_diffmap_feature_part, &
    &assemble_flex_diffmap_feature_parts, read_flex_diffmap_feature_parts, flex_projection_directions, &
    &write_flex_mean_projection_stack
use simple_flex_diffmap_preimage, only: select_flex_diffmap_preimages, build_flex_preimage_kernel_weights
use simple_flex_diffmap_rec3D,    only: reconstruct_flex_diffmap_states, reconstruct_flex_diffmap_weighted_states, &
    &write_flex_diffmap_rec_parts, &
    &reduce_flex_diffmap_rec_parts, cleanup_flex_diffmap_rec_parts, canonicalize_flex_preimage_coordinates
use simple_image,                 only: image
use simple_euclid_sigma2,         only: sigma2_star_from_iter
use simple_parameters,            only: parameters
use simple_qsys_env,              only: qsys_env
use simple_sp_project,            only: sp_project
implicit none
private
#include "simple_local_flags.inc"

public :: flex_analysis_strategy, flex_analysis_shmem_strategy
public :: flex_analysis_worker_strategy, flex_analysis_master_strategy
public :: create_flex_analysis_strategy, fit_flex_analysis_embedding, flex_embedding_result
public :: run_flex_preimage_basis_ab_test

! The embedding API is produced by this strategy and consumed by nano3D.
! Keeping the small value type here avoids a one-type transport module.
type :: flex_embedding_result
    character(len=32) :: method = 'none'
    integer :: nptcls = 0
    integer :: ncomp  = 0
    integer, allocatable :: pinds(:)
    real(dp), allocatable :: z(:,:)
    real(dp), allocatable :: eigvals(:)
contains
    procedure :: kill => kill_flex_embedding_result
end type flex_embedding_result

type, abstract :: flex_analysis_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type flex_analysis_strategy

type, extends(flex_analysis_strategy) :: flex_analysis_shmem_strategy
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type flex_analysis_shmem_strategy

type, extends(flex_analysis_strategy) :: flex_analysis_worker_strategy
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
end type flex_analysis_worker_strategy

type, extends(flex_analysis_strategy) :: flex_analysis_master_strategy
    type(qsys_env)           :: qenv
    type(chash)              :: job_descr
    type(chash), allocatable :: part_params(:)
    integer                  :: nparts_run = 1
contains
    procedure :: initialize   => master_initialize
    procedure :: execute      => master_execute
    procedure :: finalize_run => master_finalize_run
    procedure :: cleanup      => master_cleanup
end type flex_analysis_master_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: flex_analysis_strategy, parameters, builder, cmdline
        class(flex_analysis_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        type(builder),                 intent(inout) :: build
        class(cmdline),                intent(inout) :: cline
    end subroutine init_interface
    subroutine exec_interface(self, params, build, cline)
        import :: flex_analysis_strategy, parameters, builder, cmdline
        class(flex_analysis_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        type(builder),                 intent(inout) :: build
        class(cmdline),                intent(inout) :: cline
    end subroutine exec_interface
    subroutine finalize_interface(self, params, build, cline)
        import :: flex_analysis_strategy, parameters, builder, cmdline
        class(flex_analysis_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        type(builder),                 intent(inout) :: build
        class(cmdline),                intent(inout) :: cline
    end subroutine finalize_interface
    subroutine cleanup_interface(self, params)
        import :: flex_analysis_strategy, parameters
        class(flex_analysis_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    subroutine kill_flex_embedding_result( self )
        class(flex_embedding_result), intent(inout) :: self
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( allocated(self%z) ) deallocate(self%z)
        if( allocated(self%eigvals) ) deallocate(self%eigvals)
        self%method='none'; self%nptcls=0; self%ncomp=0
    end subroutine kill_flex_embedding_result

    function create_flex_analysis_strategy( cline ) result(strategy)
        class(cmdline), intent(in) :: cline
        class(flex_analysis_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker,is_master
        nparts=1
        if( cline%defined('nparts') ) nparts=max(1,cline%get_iarg('nparts'))
        is_worker=cline%defined('part')
        is_master=(nparts>1).and.(.not.is_worker)
        if( is_master )then
            allocate(flex_analysis_master_strategy::strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED FLEX_ANALYSIS (MASTER)'
        else if( is_worker )then
            allocate(flex_analysis_worker_strategy::strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> FLEX_ANALYSIS (WORKER)'
        else
            allocate(flex_analysis_shmem_strategy::strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> FLEX_ANALYSIS (SHARED-MEMORY)'
        endif
    end function create_flex_analysis_strategy

    subroutine apply_defaults( cline )
        class(cmdline), intent(inout) :: cline
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir','yes')
        if( .not.cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nstates') ) call cline%set('nstates',1)
        if( .not.cline%defined('npreimages') ) call cline%set('npreimages',8)
        if( .not.cline%defined('neigs') ) call cline%set('neigs',15)
        if( .not.cline%defined('icm') ) call cline%set('icm','yes')
        if( .not.cline%defined('k_nn') ) call cline%set('k_nn',100)
        if( .not.cline%defined('nang_nbrs') ) call cline%set('nang_nbrs',1000)
        if( .not.cline%defined('lp') ) call cline%set('lp',6.0)
        if( .not.cline%defined('bandwidth_mode') ) call cline%set('bandwidth_mode','ferguson')
        if( .not.cline%defined('bandwidth_tune') ) call cline%set('bandwidth_tune',1.0)
        if( .not.cline%defined('outvol') ) call cline%set('outvol','flex_state_001.mrc')
    end subroutine apply_defaults

    subroutine init_common( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        if( .not.cline%defined('nspace') ) &
            &THROW_HARD('flex_analysis requires nspace=<projection grid size>; it cannot be inferred from populated project rows')
        call apply_defaults(cline)
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
    end subroutine init_common

    subroutine shmem_initialize( self, params, build, cline )
        class(flex_analysis_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer :: max_modes
        call init_common(params,build,cline)
        call validate_inputs(params,cline,max_modes)
    end subroutine shmem_initialize

    subroutine shmem_execute( self, params, build, cline )
        class(flex_analysis_shmem_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, allocatable :: pinds(:)
        real, allocatable :: coords(:,:),raw_coords(:,:),spectral_z(:,:),nystrom_coords(:,:),state_weights(:,:), &
            &state_bandwidths(:),state_neff(:)
        integer, allocatable :: medoids(:),labels(:)
        type(string) :: registered_stack,registered_project
        type(builder) :: model_build
        type(parameters) :: model_params
        type(cmdline) :: model_cline
        integer :: nmodes
        call ensure_flex_sigma_support(params,build,cline,for_graph=.true.)
        call run_flex_analysis(params,build,cline,pinds,coords,raw_coords,spectral_z,nystrom_coords,nmodes, &
            &registered_stack,registered_project)
        call select_flex_diffmap_preimages(pinds,raw_coords,params%npreimages,medoids,labels)
        call build_flex_preimage_kernel_weights(raw_coords,medoids,state_weights,state_bandwidths,state_neff)
        call write_preimage_particle_map(build%spproj,pinds,coords,labels,state_weights)
        call init_model_context(cline,registered_project,model_cline,model_params,model_build)
        call ensure_flex_sigma_support(model_params,model_build,model_cline,for_graph=.false.)
        write(logfhandle,'(A,A)') '>>> FLEX PRE-IMAGE DIAGNOSTIC PROJECT (REGISTERED): ',model_params%projfile%to_char()
        call reconstruct_flex_diffmap_weighted_states(model_params,model_build,pinds,state_weights,size(medoids), &
            &fsc_projfile=params%projfile)
        call model_build%kill_general_tbox
        call model_cline%kill
        call finish_analysis_outputs(registered_stack,registered_project)
        deallocate(pinds,coords,raw_coords,spectral_z,nystrom_coords,state_weights,state_bandwidths,state_neff, &
            &medoids,labels)
    end subroutine shmem_execute

    subroutine shmem_finalize_run( self, params, build, cline )
        class(flex_analysis_shmem_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup( self, params )
        class(flex_analysis_shmem_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
    end subroutine shmem_cleanup

    subroutine worker_initialize( self, params, build, cline )
        class(flex_analysis_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        call init_common(params,build,cline)
        if( .not.cline%defined('part') ) THROW_HARD('PART must be defined for flex_analysis worker execution')
        if( .not.cline%defined('infile') ) THROW_HARD('INFILE assignment must be defined for flex_analysis worker execution')
        if( params%stage<1 .or. params%stage>3 ) THROW_HARD('invalid flex_analysis worker stage')
        if( params%stage==3 .and. params%npreimages<2 ) THROW_HARD('invalid flex_analysis worker pre-image count')
    end subroutine worker_initialize

    subroutine worker_execute( self, params, build, cline )
        use simple_qsys_funs, only: qsys_job_finished
        class(flex_analysis_worker_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, allocatable :: pinds(:)
        real, allocatable :: coords(:,:),features(:,:),proj_dirs(:,:),d2s(:,:)
        integer, allocatable :: proj_ids(:),rows(:),nbrs(:,:),ncandidates(:)
        integer :: i
        if( params%stage == 1 .or. params%stage == 3 )then
            call ensure_flex_sigma_support(params,build,cline,for_graph=params%stage==1)
        endif
        select case(params%stage)
            case(1)
                call read_int_file(params%infile,pinds)
                call prepare_flex_diffmap_feature_part(params,build,pinds,params%part)
            case(2)
                call read_int_file(params%infile,pinds)
                call read_flex_diffmap_feature_parts(params,params%nparts,features,proj_ids)
                call flex_projection_directions(build,proj_dirs)
                if( size(pinds)/=params%top-params%fromp+1 ) THROW_HARD('flex graph row assignment size mismatch')
                allocate(rows(params%top-params%fromp+1),source=[(i,i=params%fromp,params%top)])
                call find_gated_euclidean_neighbors_rows(features,proj_ids,proj_dirs,params%k_nn, &
                    &params%nang_nbrs,rows,nbrs,d2s,ncandidates)
                call write_graph_part(params,rows,nbrs,d2s,ncandidates)
                deallocate(features,proj_ids,proj_dirs,rows,nbrs,d2s,ncandidates)
            case(3)
                call read_flex_reconstruction_assignment(params%infile,params%neigs,pinds,coords)
                if( any(pinds>build%spproj%os_ptcl3D%get_noris()) ) &
                    &THROW_HARD('flex reconstruction assignment particle outside native project')
                call write_flex_diffmap_rec_parts(params,build,pinds,coords,params%neigs,params%part)
                deallocate(coords)
        end select
        call qsys_job_finished(params,string('simple_flex_analysis_strategy :: worker_execute'))
        deallocate(pinds)
    end subroutine worker_execute

    subroutine worker_finalize_run( self, params, build, cline )
        class(flex_analysis_worker_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
    end subroutine worker_finalize_run

    subroutine worker_cleanup( self, params )
        class(flex_analysis_worker_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
    end subroutine worker_cleanup

    subroutine master_initialize( self, params, build, cline )
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer :: max_modes
        call init_common(params,build,cline)
        call validate_inputs(params,cline,max_modes)
    end subroutine master_initialize

    subroutine master_execute( self, params, build, cline )
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        type(string) :: registered_stack,registered_project,worker_program
        integer, allocatable :: pinds(:),medoids(:),labels(:),proj_ids(:)
        real, allocatable :: coords(:,:),raw_coords(:,:),spectral_z(:,:),nystrom_coords(:,:),state_weights(:,:), &
            &state_bandwidths(:),state_neff(:),features(:,:),proj_dirs(:,:)
        type(diffmap_graph) :: graph
        type(builder) :: model_build
        type(parameters) :: model_params
        type(cmdline) :: model_cline
        integer :: nmodes,nptcls,max_modes,cand_min,cand_max
        real :: cand_mean
        integer(timer_int_kind) :: t_step
        call ensure_flex_sigma_support(params,build,cline,for_graph=.true.)
        call validate_inputs(params,cline,max_modes)
        call select_particles(params,build,pinds,nptcls)
        call validate_selected_particles(build%spproj,pinds)
        self%nparts_run=min(max(1,params%nparts),size(pinds))
        write(logfhandle,'(A,I0,A,I0,A,I0)') 'Flex worker planning: requested_nparts=',params%nparts, &
            &' selected_particles=',size(pinds),' effective_nparts=',self%nparts_run
        write(logfhandle,'(A)') 'Flex partition unit: contiguous selected-particle rows for features, graph rows, and reconstruction.'
        call self%qenv%new(params,self%nparts_run,numlen=params%numlen,nptcls=size(pinds))
        call prepare_particle_partitions(self,params,pinds)
        worker_program='flex_analysis'
        call cline%gen_job_descr(self%job_descr,prg=worker_program)
        call self%job_descr%set('mkdir','no')
        call self%job_descr%set('nparts',int2str(self%nparts_run))
        call self%job_descr%set('numlen',int2str(params%numlen))

        call self%job_descr%set('stage','1')
        t_step=tic()
        call write_flex_mean_projection_stack(params,build)
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr,part_params=self%part_params, &
            &array=L_USE_SLURM_ARR,extra_params=params)
        call assemble_flex_diffmap_feature_parts(params,build%spproj,pinds,self%nparts_run,self%qenv%parts,registered_project)
        registered_stack='flex_registered_particles_part*.mrcs'
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP distributed_registration_seconds=',toc(t_step)

        call self%job_descr%set('stage','2')
        t_step=tic()
        ! The gated graph requires all residual features as candidates.  Keep
        ! the registration stage distributed, but assemble and evaluate this
        ! global table once in the master: this avoids both duplicated worker
        ! allocations and a second qsys launch over the same part assignments.
        write(logfhandle,'(A)') '>>> FLEX DIFFMAP graph evaluation: master (one global feature table)'
        call flush(logfhandle)
        call read_flex_diffmap_feature_parts(params,self%nparts_run,features,proj_ids,build,pinds)
        call flex_projection_directions(build,proj_dirs)
        call build_gated_euclidean_knn_graph(features,proj_ids,proj_dirs,params%k_nn,params%nang_nbrs,graph, &
            &cand_min,cand_max,cand_mean,params%bandwidth_mode,params%bandwidth_tune)
        deallocate(features,proj_ids,proj_dirs)
        call cleanup_distributed_analysis_parts(params,self%nparts_run)
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP distributed_graph_seconds=',toc(t_step)
        call embed_flex_graph(params,pinds,graph,max_modes,cand_min,cand_max,cand_mean,coords,raw_coords, &
            &spectral_z,nystrom_coords,nmodes)
        call select_flex_diffmap_preimages(pinds,raw_coords,params%npreimages,medoids,labels)
        call build_flex_preimage_kernel_weights(raw_coords,medoids,state_weights,state_bandwidths,state_neff)
        call write_preimage_particle_map(build%spproj,pinds,coords,labels,state_weights)
        call graph%kill()
        call init_model_context(cline,registered_project,model_cline,model_params,model_build)
        call ensure_flex_sigma_support(model_params,model_build,model_cline,for_graph=.false.)
        write(logfhandle,'(A,A)') '>>> FLEX PRE-IMAGE DIAGNOSTIC PROJECT (REGISTERED): ',model_params%projfile%to_char()
        t_step=tic()
        call reconstruct_flex_diffmap_weighted_states(model_params,model_build,pinds,state_weights,size(medoids), &
            &fsc_projfile=params%projfile)
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP distributed_reconstruction_seconds=',toc(t_step)
        call model_build%kill_general_tbox
        call model_cline%kill
        call worker_program%kill
        call finish_analysis_outputs(registered_stack,registered_project)
        deallocate(pinds,coords,raw_coords,spectral_z,nystrom_coords,state_weights,state_bandwidths,state_neff, &
            &medoids,labels)
    end subroutine master_execute

    subroutine init_model_context( source_cline, model_project, model_cline, model_params, model_build )
        class(cmdline), intent(in) :: source_cline
        type(string), intent(in) :: model_project
        type(cmdline), intent(inout) :: model_cline
        type(parameters), intent(inout) :: model_params
        type(builder), intent(inout) :: model_build
        model_cline=source_cline
        call model_cline%set('projfile',model_project%to_char())
        call model_cline%set('ptcl_src','raw')
        call model_cline%set('mkdir','no')
        call init_common(model_params,model_build,model_cline)
    end subroutine init_model_context

    subroutine ensure_flex_sigma_support( params, build, cline, for_graph )
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        logical, intent(in) :: for_graph
        type(string) :: sigma_part_fname, sigma_group_fname
        integer :: fromp_saved, top_saved, noris
        logical :: has_group, loaded
        call carry_over_prior_sigma_files(params)
        call pick_sigma_group_file(params,sigma_group_fname,has_group)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            sigma_part_fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
        endif
        loaded=.false.
        if( has_group .and. params%cc_objfun == OBJFUN_EUCLID )then
            fromp_saved=params%fromp
            top_saved=params%top
            noris=build%spproj_field%get_noris()
            if( noris > 0 )then
                ! Reused sigma-star files from upstream refine/abinitio runs are
                ! often global (single-group). Force global sigma dispatch to
                ! avoid indexing group-by-stkind against a 1-group table.
                params%sigma_est='global'
                params%l_sigma_glob=.true.
                call cline%set('sigma_est','global')
                params%fromp=1
                params%top=noris
                call build%esig%new(params,build%pftc,sigma_part_fname,params%box)
                call build%esig%read_groups(build%spproj_field,fname=sigma_group_fname)
                call build%esig%allocate_ptcls()
                loaded=allocated(build%esig%sigma2_noise)
                params%fromp=fromp_saved
                params%top=top_saved
            endif
        endif
        if( loaded .and. params%cc_objfun == OBJFUN_EUCLID )then
            params%ml_reg='yes'
            params%l_ml_reg=.true.
            call cline%set('ml_reg','yes')
            write(logfhandle,'(A)') '>>> FLEX SIGMA SUPPORT: enabling ml_reg (reused sigma files found)'
        else
            params%ml_reg='no'
            params%l_ml_reg=.false.
            call cline%set('ml_reg','no')
            write(logfhandle,'(A)') '>>> FLEX SIGMA SUPPORT: disabling ml_reg (no reusable sigma files found)'
        endif
        if( for_graph )then
            if( loaded )then
                write(logfhandle,'(A)') '>>> FLEX SIGMA SUPPORT: graph distances will use sigma-scaled residuals'
            else
                write(logfhandle,'(A)') '>>> FLEX SIGMA SUPPORT: graph distances use unscaled residuals (current behavior fallback)'
            endif
        endif
        call sigma_part_fname%kill
        call sigma_group_fname%kill
    end subroutine ensure_flex_sigma_support

    subroutine carry_over_prior_sigma_files( params )
        type(parameters), intent(in) :: params
        type(sp_project) :: inproj
        type(string), allocatable :: siblings(:)
        type(string) :: proj_dir,parent_dir,proj_cwd,proj_base,sibling_dir
        logical :: have_proj_cwd
        integer :: i
        proj_dir=get_fpath(params%projfile)
        parent_dir=string(PATH_PARENT)
        proj_base=basename(params%projfile)
        proj_cwd=''
        have_proj_cwd=.false.
        ! Prefer sibling-run discovery from ../, e.g. ../5_abinitio3D/<projfile>.
        if( parent_dir /= '' .and. proj_base /= '' )then
            call simple_list_files(parent_dir%to_char()//'*/'//proj_base%to_char(), siblings)
            do i=1,size(siblings)
                sibling_dir=get_fpath(siblings(i))
                if( sibling_dir /= '' ) call copy_sigma_files_from_dir(sibling_dir)
                call sibling_dir%kill
            end do
            if( allocated(siblings) ) deallocate(siblings)
        endif
        if( file_exists(params%projfile) )then
            call inproj%read_non_data_segments(params%projfile)
            if( inproj%projinfo%isthere('cwd') )then
                call inproj%projinfo%getter(1,'cwd',proj_cwd)
                have_proj_cwd = proj_cwd /= ''
            endif
            call inproj%kill
        endif
        if( have_proj_cwd ) call copy_sigma_files_from_dir(proj_cwd)
        if( proj_dir /= '' ) call copy_sigma_files_from_dir(proj_dir)
        if( parent_dir /= '' ) call copy_sigma_files_from_dir(parent_dir)
        call proj_base%kill
        call proj_cwd%kill
        call proj_dir%kill
        call parent_dir%kill
    end subroutine carry_over_prior_sigma_files

    subroutine copy_sigma_files_from_dir( src_dir )
        type(string), intent(in) :: src_dir
        type(string), allocatable :: list(:)
        type(string) :: target_name
        character(len=:), allocatable :: src_prefix
        integer :: i
        src_prefix=trim(src_dir%to_char())
        if( len_trim(src_prefix) < 1 ) return
        if( src_prefix(len_trim(src_prefix):len_trim(src_prefix)) /= '/' ) src_prefix=src_prefix//'/'
        call simple_list_files(src_prefix//SIGMA2_FBODY//'*', list)
        do i=1,size(list)
            target_name=string(PATH_HERE)//basename(list(i))
            if( .not. file_exists(target_name) ) call simple_copy_file(list(i),target_name)
            call target_name%kill
        end do
        if( allocated(list) ) deallocate(list)
        call simple_list_files(src_prefix//SIGMA2_GROUP_FBODY//'*'//STAR_EXT, list)
        do i=1,size(list)
            target_name=string(PATH_HERE)//basename(list(i))
            if( .not. file_exists(target_name) ) call simple_copy_file(list(i),target_name)
            call target_name%kill
        end do
        if( allocated(list) ) deallocate(list)
    end subroutine copy_sigma_files_from_dir

    subroutine pick_sigma_group_file( params, sigma_group_fname, found )
        type(parameters), intent(in) :: params
        type(string), intent(out) :: sigma_group_fname
        logical, intent(out) :: found
        type(string), allocatable :: list(:)
        integer :: i,best_iter,iter_here
        sigma_group_fname=''
        found=.false.
        sigma_group_fname=sigma2_star_from_iter(params%which_iter)
        if( file_exists(sigma_group_fname) )then
            found=.true.
            return
        endif
        call sigma_group_fname%kill
        call simple_list_files(SIGMA2_GROUP_FBODY//'*'//STAR_EXT, list)
        if( .not. allocated(list) ) return
        if( size(list) < 1 )then
            deallocate(list)
            return
        endif
        best_iter=-1
        do i=1,size(list)
            iter_here = trailing_iter_from_sigma_group_name(basename(list(i)))
            if( iter_here >= best_iter )then
                best_iter=iter_here
                sigma_group_fname=basename(list(i))
                found=.true.
            endif
        end do
        deallocate(list)
    end subroutine pick_sigma_group_file

    integer function trailing_iter_from_sigma_group_name( fname ) result(iter)
        type(string), intent(in) :: fname
        character(len=:), allocatable :: raw
        integer :: i,scale,d
        raw=fname%to_char()
        iter=0
        scale=1
        do i=len_trim(raw),1,-1
            if( raw(i:i) >= '0' .and. raw(i:i) <= '9' )then
                d=iachar(raw(i:i))-iachar('0')
                iter=iter+d*scale
                scale=scale*10
            else if( scale > 1 )then
                exit
            endif
        end do
    end function trailing_iter_from_sigma_group_name

    subroutine master_finalize_run( self, params, build, cline )
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
    end subroutine master_finalize_run

    subroutine master_cleanup( self, params )
        use simple_qsys_funs, only: qsys_cleanup
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        type(string) :: fname
        integer :: ipart
        call self%qenv%kill
        call qsys_cleanup(params)
        if( allocated(self%part_params) )then
            do ipart=1,size(self%part_params)
                fname=string('flex_particles_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
                call del_file(fname); call fname%kill
                call self%part_params(ipart)%kill
            end do
            deallocate(self%part_params)
        endif
        call self%job_descr%kill
    end subroutine master_cleanup

    subroutine prepare_particle_partitions( self, params, pinds )
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        integer, intent(in) :: pinds(:)
        type(string) :: fname
        integer :: ipart,first,last
        allocate(self%part_params(self%nparts_run))
        do ipart=1,self%nparts_run
            first=self%qenv%parts(ipart,1); last=self%qenv%parts(ipart,2)
            fname=string('flex_particles_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            call arr2txtfile(pinds(first:last),fname)
            call self%part_params(ipart)%new(1)
            call self%part_params(ipart)%set('infile',fname%to_char())
            call fname%kill
        end do
    end subroutine prepare_particle_partitions

    ! Stage 3 uses the normal qsys assignment file, extended with the
    ! per-particle spectral vector.  Particle auxiliary metadata is not a
    ! persistent project transport: SIMPLE serializes only the fixed particle
    ! parameters.  Keeping the vector in the assignment makes the worker
    ! handoff explicit and preserves the selected-row ordering.
    subroutine prepare_reconstruction_partitions( self, params, model_pinds, spectral_z )
        class(flex_analysis_master_strategy), intent(inout) :: self
        type(parameters), intent(in) :: params
        integer, intent(in) :: model_pinds(:)
        real, intent(in) :: spectral_z(:,:)
        type(string) :: fname
        integer :: ipart,first,last,i,u
        if( size(model_pinds)<1 .or. size(spectral_z,1)/=size(model_pinds) .or. size(spectral_z,2)<1 ) &
            &THROW_HARD('invalid flex reconstruction partition table')
        if( .not.allocated(self%part_params) .or. size(self%part_params)/=self%nparts_run ) &
            &THROW_HARD('flex reconstruction partitions were not initialized')
        do ipart=1,self%nparts_run
            first=self%qenv%parts(ipart,1); last=self%qenv%parts(ipart,2)
            if( first<1 .or. last<first .or. last>size(model_pinds) ) &
                &THROW_HARD('invalid flex reconstruction partition bounds')
            fname=string('flex_particles_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            open(newunit=u,file=fname%to_char(),status='replace',action='write')
            write(u,'(A)') '# registered_particle_index spectral_coordinates'
            do i=first,last
                write(u,*) model_pinds(i),spectral_z(i,:)
            end do
            close(u)
            call fname%kill
        end do
    end subroutine prepare_reconstruction_partitions

    subroutine fit_flex_analysis_embedding( params, build, cline, fit_result )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        type(flex_embedding_result), intent(inout) :: fit_result
        integer, allocatable :: pinds(:)
        real, allocatable :: coords(:,:),raw_coords(:,:),spectral_z(:,:),nystrom_coords(:,:)
        type(string) :: registered_stack,registered_project
        integer :: nmodes
        call run_flex_analysis(params,build,cline,pinds,coords,raw_coords,spectral_z,nystrom_coords,nmodes, &
            &registered_stack,registered_project,fit_result)
        write(logfhandle,'(A)') '>>> FLEX DIFFMAP reconstruction=skipped (embedding-only caller)'
        call finish_analysis_outputs(registered_stack,registered_project)
        deallocate(pinds,coords,raw_coords,spectral_z,nystrom_coords)
    end subroutine fit_flex_analysis_embedding

    !> Strict numerical A/B diagnostic for the residual Nyström pre-image.
    !! A uses the raw graph eigenfunctions.  B whitens their full-rank span
    !! and transforms each Nyström target by the inverse basis transform, so
    !! A and B denote the same representative states.  They may differ only
    !! through finite-precision conditioning or an implementation defect.
    subroutine run_flex_preimage_basis_ab_test( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, allocatable :: pinds(:), medoids(:), labels(:)
        real, allocatable :: coords(:,:), raw_coords(:,:), spectral_z(:,:), nystrom_coords(:,:)
        real, allocatable :: target_raw(:,:), z_canonical(:,:), target_canonical(:,:)
        type(string) :: registered_stack, registered_project
        type(builder) :: model_build
        type(parameters) :: model_params, raw_params, canonical_params
        type(cmdline) :: model_cline
        real(dp), allocatable :: transform(:,:)
        real(dp) :: transform_condition, gram_error, target_error
        integer :: nmodes
        call run_flex_analysis(params,build,cline,pinds,coords,raw_coords,spectral_z,nystrom_coords,nmodes, &
            &registered_stack,registered_project)
        call select_flex_diffmap_preimages(pinds,raw_coords,params%npreimages,medoids,labels)
        call collect_target_coefficients(spectral_z,nystrom_coords,medoids,target_raw)
        allocate(z_canonical(size(spectral_z,1),nmodes),target_canonical(size(target_raw,1),nmodes), &
            &transform(nmodes,nmodes))
        call canonicalize_flex_preimage_coordinates(spectral_z,target_raw,z_canonical,target_canonical,transform, &
            &transform_condition,gram_error,target_error)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX PRE-IMAGE A/B canonical_transform_condition=', &
            &transform_condition,' gram_error=',gram_error,' target_error=',target_error
        call flush(logfhandle)
        call init_model_context(cline,registered_project,model_cline,model_params,model_build)
        raw_params=model_params
        raw_params%outvol='flex_preimage_ab_raw_001.mrc'
        canonical_params=model_params
        canonical_params%outvol='flex_preimage_ab_canonical_001.mrc'
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE A/B: A=RAW GRAPH EIGENFUNCTIONS'
        call reconstruct_flex_diffmap_states(raw_params,model_build,pinds,spectral_z,target_raw,size(medoids))
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE A/B: B=WHITENED GRAPH EIGENFUNCTIONS, INVERSE TARGETS'
        call reconstruct_flex_diffmap_states(canonical_params,model_build,pinds,z_canonical,target_canonical,size(medoids))
        call write_preimage_basis_ab_metrics(model_params,size(medoids),transform_condition,gram_error,target_error)
        call model_build%kill_general_tbox
        call model_cline%kill
        call finish_analysis_outputs(registered_stack,registered_project)
        deallocate(pinds,coords,raw_coords,spectral_z,nystrom_coords,target_raw,z_canonical,target_canonical, &
            &transform,medoids,labels)
    end subroutine run_flex_preimage_basis_ab_test

    subroutine write_preimage_basis_ab_metrics( params, nstates, transform_condition, gram_error, target_error )
        type(parameters), intent(in) :: params
        integer, intent(in) :: nstates
        real(dp), intent(in) :: transform_condition, gram_error, target_error
        type(image) :: raw_img, canonical_img
        type(string) :: raw_fname, canonical_fname
        integer :: state, funit
        real(dp) :: map_cc
        open(newunit=funit,file='flex_preimage_basis_ab_metrics.txt',status='replace',action='write')
        write(funit,'(A,ES16.8)') 'canonical_transform_condition=',transform_condition
        write(funit,'(A,ES16.8)') 'canonical_gram_error=',gram_error
        write(funit,'(A,ES16.8)') 'inverse_target_error=',target_error
        write(funit,'(A)') '# state raw_vs_canonical_cc'
        do state=1,nstates
            raw_fname=preimage_ab_state_fname('raw',state)
            canonical_fname=preimage_ab_state_fname('canonical',state)
            call raw_img%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call canonical_img%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call raw_img%read(raw_fname)
            call canonical_img%read(canonical_fname)
            map_cc=raw_img%corr(canonical_img)
            write(funit,'(I0,1X,ES16.8)') state,map_cc
            write(logfhandle,'(A,I0,A,ES12.4)') '>>> FLEX PRE-IMAGE A/B STATE=',state,' raw_vs_canonical_cc=',map_cc
            call raw_img%kill
            call canonical_img%kill
            call raw_fname%kill
            call canonical_fname%kill
        end do
        close(funit)
    end subroutine write_preimage_basis_ab_metrics

    function preimage_ab_state_fname( label, state ) result(fname)
        character(len=*), intent(in) :: label
        integer, intent(in) :: state
        type(string) :: fname
        character(len=3) :: tag
        write(tag,'(I3.3)') state
        fname='flex_preimage_ab_'//label//'_'//tag//'.mrc'
    end function preimage_ab_state_fname

    subroutine run_flex_analysis( params, build, cline, pinds, coords, raw_coords, spectral_z, nystrom_coords, &
        &nmodes, registered_stack, &
        &registered_project, fit_result )
        type(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, allocatable, intent(out) :: pinds(:)
        real, allocatable, intent(out) :: coords(:,:)
        real, allocatable, intent(out) :: raw_coords(:,:)
        real, allocatable, intent(out) :: spectral_z(:,:),nystrom_coords(:,:)
        integer, intent(out) :: nmodes
        type(string), intent(out) :: registered_stack,registered_project
        type(flex_embedding_result), optional, intent(inout) :: fit_result
        type(diffmap_graph) :: graph
        integer, allocatable :: proj_ids(:)
        real, allocatable :: features(:,:),proj_dirs(:,:),coords_mode_major(:,:),raw_mode_major(:,:),eigvals(:)
        real, allocatable :: spectral_mode_major(:,:),nystrom_mode_major(:,:)
        integer :: nptcls,max_modes,icm_iters,cand_min,cand_max
        real :: cand_mean,icm_score
        logical :: icm_converged
        integer(timer_int_kind) :: t_step
        call validate_inputs(params,cline,max_modes)
        call select_particles(params,build,pinds,nptcls)
        call validate_selected_particles(build%spproj,pinds)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,F8.3)') &
            '>>> FLEX DIFFMAP particles=',nptcls,' k_nn=',params%k_nn,' nang_nbrs=',params%nang_nbrs, &
            &' max_modes=',max_modes,' lp=',params%lp
        call flush(logfhandle)
        t_step=tic()
        call prepare_flex_diffmap_features(params,build,pinds,features,proj_ids,proj_dirs,registered_stack,registered_project)
        write(logfhandle,'(A,F10.3,A,I0)') '>>> FLEX DIFFMAP registration_seconds=',toc(t_step), &
            &' feature_values=',size(features,kind=8)
        t_step=tic()
        call build_gated_euclidean_knn_graph(features,proj_ids,proj_dirs,params%k_nn,params%nang_nbrs,graph, &
            &cand_min,cand_max,cand_mean,params%bandwidth_mode,params%bandwidth_tune)
        write(logfhandle,'(A,I0,A,I0,A,F8.1,A,I0)') '>>> FLEX DIFFMAP graph candidates_min=',cand_min, &
            &' max=',cand_max,' mean=',cand_mean,' directed_nnz=',graph%nnz
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP graph_seconds=',toc(t_step)
        deallocate(features,proj_dirs,proj_ids)
        t_step=tic()
        call embed_graph(graph,max_modes,coords_mode_major,eigvals,raw_mode_major,spectral_mode_major,nystrom_mode_major)
        if( size(eigvals)<1 ) THROW_HARD('diffusion embedding returned no nontrivial modes')
        call select_flex_spectral_rank(params,eigvals,nmodes,icm_converged,icm_iters,icm_score)
        nmodes=min(max(1,nmodes),min(max_modes,size(eigvals)))
        allocate(coords(nptcls,nmodes),source=transpose(coords_mode_major(1:nmodes,:)))
        allocate(raw_coords(nptcls,nmodes),source=transpose(raw_mode_major(1:nmodes,:)))
        allocate(spectral_z(nptcls,nmodes),source=transpose(spectral_mode_major(1:nmodes,:)))
        allocate(nystrom_coords(nptcls,nmodes),source=transpose(nystrom_mode_major(1:nmodes,:)))
        if( .not.all(ieee_is_finite(coords)) ) THROW_HARD('diffusion embedding produced nonfinite coordinates')
        call log_flex_coordinate_scales(spectral_z,nystrom_coords)
        write(logfhandle,'(A,I0,A,L1,A,L1,A,I0,A,ES12.4,A,F10.3)') '>>> FLEX DIFFMAP selected_modes=',nmodes, &
            &' icm_enabled=',params%l_icm,' icm_converged=',icm_converged,' icm_iters=',icm_iters,' icm_score=',icm_score, &
            &' embedding_seconds=',toc(t_step)
        call write_coordinates(pinds,coords,nptcls,nmodes)
        call write_spectrum(eigvals,nmodes,params%l_icm,icm_converged,icm_iters,icm_score)
        call write_graph_summary(graph,cand_min,cand_max,cand_mean,params)
        if( present(fit_result) ) call capture_result(fit_result,pinds,coords,eigvals,nptcls,nmodes)
        call graph%kill()
        deallocate(coords_mode_major,raw_mode_major,spectral_mode_major,nystrom_mode_major,eigvals)
    end subroutine run_flex_analysis

    subroutine read_flex_reconstruction_assignment( fname, nmodes, pinds, spectral_z )
        type(string), intent(in) :: fname
        integer, intent(in) :: nmodes
        integer, allocatable, intent(out) :: pinds(:)
        real, allocatable, intent(out) :: spectral_z(:,:)
        character(len=XLONGSTRLEN) :: line
        integer :: nvals,funit,ios,i
        if( nmodes<1 ) THROW_HARD('invalid flex reconstruction spectral rank')
        nvals=0
        open(newunit=funit,file=fname%to_char(),status='old',action='read',iostat=ios)
        call fileiochk('opening flex reconstruction assignment '//fname%to_char(),ios)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios/=0 ) exit
            if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
            nvals=nvals+1
        end do
        close(funit)
        if( nvals<1 ) THROW_HARD('empty flex reconstruction assignment file')
        allocate(pinds(nvals),spectral_z(nvals,nmodes))
        open(newunit=funit,file=fname%to_char(),status='old',action='read',iostat=ios)
        call fileiochk('opening flex reconstruction assignment '//fname%to_char(),ios)
        i=0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios/=0 ) exit
            if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
            i=i+1
            read(line,*,iostat=ios) pinds(i),spectral_z(i,:)
            if( ios/=0 ) THROW_HARD('invalid flex reconstruction assignment row')
        end do
        close(funit)
        if( any(pinds<1) .or. .not.all(ieee_is_finite(spectral_z)) ) &
            &THROW_HARD('invalid flex reconstruction assignment values')
    end subroutine read_flex_reconstruction_assignment

    subroutine collect_target_coefficients( spectral_z, nystrom_coords, medoids, target_coeffs )
        real, intent(in) :: spectral_z(:,:), nystrom_coords(:,:)
        integer, intent(in) :: medoids(:)
        real, allocatable, intent(out) :: target_coeffs(:,:)
        real(dp) :: spec_rms, nys_rms, mode_scale, mode_min_scale, mode_max_scale
        integer :: state, q
        logical :: used_fallback
        if( any(shape(spectral_z)/=shape(nystrom_coords)) ) THROW_HARD('flex spectral/Nystrom table shape mismatch')
        allocate(target_coeffs(size(medoids),size(nystrom_coords,2)))
        mode_min_scale=huge(1.d0)
        mode_max_scale=0.d0
        used_fallback=.false.
        do q=1,size(nystrom_coords,2)
            spec_rms=sqrt(sum(real(spectral_z(:,q),dp)**2)/real(max(1,size(spectral_z,1)),dp))
            nys_rms=sqrt(sum(real(nystrom_coords(:,q),dp)**2)/real(max(1,size(nystrom_coords,1)),dp))
            if( nys_rms>DTINY )then
                mode_scale=spec_rms/nys_rms
            else
                mode_scale=1.d0
                used_fallback=.true.
            endif
            mode_min_scale=min(mode_min_scale,mode_scale)
            mode_max_scale=max(mode_max_scale,mode_scale)
            do state=1,size(medoids)
                if( medoids(state)<1 .or. medoids(state)>size(nystrom_coords,1) ) &
                    &THROW_HARD('flex pre-image medoid outside spectral table')
                if( nys_rms>DTINY )then
                    target_coeffs(state,q)=real(mode_scale*nystrom_coords(medoids(state),q))
                else
                    target_coeffs(state,q)=spectral_z(medoids(state),q)
                endif
            end do
        end do
        if( .not.all(ieee_is_finite(target_coeffs)) ) THROW_HARD('invalid flex Nystrom target coefficient')
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,L1)') '>>> FLEX PRE-IMAGE target_scale min=',mode_min_scale, &
            &' max=',mode_max_scale,' fallback_to_spectral=',used_fallback
    end subroutine collect_target_coefficients

    subroutine read_int_file( fname, vals )
        type(string), intent(in) :: fname
        integer, allocatable, intent(out) :: vals(:)
        integer :: nvals,funit,ios,i
        character(len=XLONGSTRLEN) :: line
        nvals=0
        open(newunit=funit,file=fname%to_char(),status='old',action='read',iostat=ios)
        call fileiochk('opening flex particle assignment '//fname%to_char(),ios)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios/=0 ) exit
            if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
            nvals=nvals+1
        end do
        close(funit)
        if( nvals<1 ) THROW_HARD('empty flex particle assignment file')
        allocate(vals(nvals))
        open(newunit=funit,file=fname%to_char(),status='old',action='read',iostat=ios)
        call fileiochk('opening flex particle assignment '//fname%to_char(),ios)
        i=0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios/=0 ) exit
            if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
            i=i+1; read(line,*) vals(i)
        end do
        close(funit)
    end subroutine read_int_file

    subroutine write_graph_part( params, rows, nbrs, d2s, ncandidates )
        type(parameters), intent(in) :: params
        integer, intent(in) :: rows(:),nbrs(:,:),ncandidates(:)
        real, intent(in) :: d2s(:,:)
        type(string) :: fname
        integer :: u,ir,m
        if( size(nbrs,2)/=size(rows).or.any(shape(d2s)/=shape(nbrs)) ) &
            &THROW_HARD('invalid flex graph part dimensions')
        fname=string('flex_graph_neighbors_part')//int2str_pad(params%part,params%numlen)//TXT_EXT
        call del_file(fname)
        open(newunit=u,file=fname%to_char(),status='replace',action='write')
        write(u,'(A)') '# row ncandidates (neighbor squared_distance)*'
        do ir=1,size(rows)
            write(u,'(I10,1X,I10)',advance='no') rows(ir),ncandidates(ir)
            do m=1,size(nbrs,1)
                write(u,'(1X,I10,1X,ES16.8)',advance='no') nbrs(m,ir),d2s(m,ir)
            end do
            write(u,*)
        end do
        close(u); call fname%kill
    end subroutine write_graph_part

    subroutine read_graph_parts( params, nptcls, nparts, graph, cmin, cmax, cmean )
        type(parameters), intent(in) :: params
        integer, intent(in) :: nptcls,nparts
        type(diffmap_graph), intent(out) :: graph
        integer, intent(out) :: cmin,cmax
        real, intent(out) :: cmean
        integer, allocatable :: nbrs(:,:),ncandidates(:)
        real, allocatable :: d2s(:,:)
        logical, allocatable :: covered(:)
        type(string) :: fname
        character(len=XLONGSTRLEN) :: line
        integer, allocatable :: row_nbrs(:)
        real, allocatable :: row_d2s(:)
        integer :: k_used,ipart,u,ios,row,m,nread,row_ncandidates
        k_used=min(max(1,params%k_nn),nptcls-1)
        allocate(nbrs(k_used,nptcls),ncandidates(nptcls),source=0)
        allocate(d2s(k_used,nptcls),source=0.)
        allocate(covered(nptcls),source=.false.)
        allocate(row_nbrs(k_used),row_d2s(k_used))
        do ipart=1,nparts
            fname=string('flex_graph_neighbors_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            open(newunit=u,file=fname%to_char(),status='old',action='read',iostat=ios)
            call fileiochk('opening flex graph part '//fname%to_char(),ios)
            nread=0
            do
                read(u,'(A)',iostat=ios) line
                if( ios/=0 ) exit
                if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
                read(line,*) row,row_ncandidates,(row_nbrs(m),row_d2s(m),m=1,k_used)
                if( row<1 .or. row>nptcls ) THROW_HARD('flex graph part row outside table')
                if( covered(row) ) THROW_HARD('duplicate flex graph part row')
                ncandidates(row)=row_ncandidates; nbrs(:,row)=row_nbrs; d2s(:,row)=row_d2s
                covered(row)=.true.; nread=nread+1
            end do
            close(u); call fname%kill
            if( nread<1 ) THROW_HARD('empty flex graph part')
        end do
        if( any(.not.covered) ) THROW_HARD('distributed flex graph parts do not cover every particle')
        call build_gated_euclidean_graph_from_neighbors(nptcls,nbrs,d2s,ncandidates,graph, &
            &params%bandwidth_mode,params%bandwidth_tune)
        cmin=minval(ncandidates); cmax=maxval(ncandidates)
        cmean=real(sum(int(ncandidates,kind=8)),kind=sp)/real(nptcls,kind=sp)
        deallocate(nbrs,ncandidates,d2s,covered,row_nbrs,row_d2s)
    end subroutine read_graph_parts

    subroutine embed_flex_graph( params, pinds, graph, max_modes, cmin, cmax, cmean, coords, raw_coords, &
        &spectral_z, nystrom_coords, nmodes )
        type(parameters), intent(in) :: params
        integer, intent(in) :: pinds(:),max_modes,cmin,cmax
        type(diffmap_graph), intent(in) :: graph
        real, intent(in) :: cmean
        real, allocatable, intent(out) :: coords(:,:)
        real, allocatable, intent(out) :: raw_coords(:,:)
        real, allocatable, intent(out) :: spectral_z(:,:),nystrom_coords(:,:)
        integer, intent(out) :: nmodes
        real, allocatable :: coords_mode_major(:,:),raw_mode_major(:,:),eigvals(:)
        real, allocatable :: spectral_mode_major(:,:),nystrom_mode_major(:,:)
        real :: icm_score
        integer :: icm_iters
        logical :: icm_converged
        integer(timer_int_kind) :: t_step
        write(logfhandle,'(A,I0,A,I0,A,F8.1,A,I0)') '>>> FLEX DIFFMAP graph candidates_min=',cmin, &
            &' max=',cmax,' mean=',cmean,' directed_nnz=',graph%nnz
        t_step=tic()
        call embed_graph(graph,max_modes,coords_mode_major,eigvals,raw_mode_major,spectral_mode_major,nystrom_mode_major)
        if( size(eigvals)<1 ) THROW_HARD('diffusion embedding returned no nontrivial modes')
        call select_flex_spectral_rank(params,eigvals,nmodes,icm_converged,icm_iters,icm_score)
        nmodes=min(max(1,nmodes),min(max_modes,size(eigvals)))
        allocate(coords(size(pinds),nmodes),source=transpose(coords_mode_major(1:nmodes,:)))
        allocate(raw_coords(size(pinds),nmodes),source=transpose(raw_mode_major(1:nmodes,:)))
        allocate(spectral_z(size(pinds),nmodes),source=transpose(spectral_mode_major(1:nmodes,:)))
        allocate(nystrom_coords(size(pinds),nmodes),source=transpose(nystrom_mode_major(1:nmodes,:)))
        if( .not.all(ieee_is_finite(coords)) ) THROW_HARD('diffusion embedding produced nonfinite coordinates')
        call log_flex_coordinate_scales(spectral_z,nystrom_coords)
        write(logfhandle,'(A,I0,A,L1,A,L1,A,I0,A,ES12.4,A,F10.3)') '>>> FLEX DIFFMAP selected_modes=',nmodes, &
            &' icm_enabled=',params%l_icm,' icm_converged=',icm_converged,' icm_iters=',icm_iters,' icm_score=',icm_score, &
            &' embedding_seconds=',toc(t_step)
        call write_coordinates(pinds,coords,size(pinds),nmodes)
        call write_spectrum(eigvals,nmodes,params%l_icm,icm_converged,icm_iters,icm_score)
        call write_graph_summary(graph,cmin,cmax,cmean,params)
        deallocate(coords_mode_major,raw_mode_major,spectral_mode_major,nystrom_mode_major,eigvals)
    end subroutine embed_flex_graph

    subroutine select_flex_spectral_rank( params, eigvals, nmodes, converged, niters, score )
        type(parameters), intent(in) :: params
        real, intent(in) :: eigvals(:)
        integer, intent(out) :: nmodes,niters
        real, intent(out) :: score
        logical, intent(out) :: converged
        if( size(eigvals)<1 ) THROW_HARD('cannot select flex spectral rank from an empty spectrum')
        if( params%l_icm )then
            call select_spectral_rank_icm(eigvals,size(eigvals),nmodes,converged,niters,score,min_rank=1)
        else
            nmodes=size(eigvals)
            converged=.false.
            niters=0
            score=0.
        endif
    end subroutine select_flex_spectral_rank

    subroutine log_flex_coordinate_scales( spectral_z, nystrom_coords )
        real, intent(in) :: spectral_z(:,:),nystrom_coords(:,:)
        real(dp) :: spec_rms, nys_rms, ratio
        if( size(spectral_z,1)<1 .or. size(spectral_z,2)<1 ) return
        if( any(shape(spectral_z)/=shape(nystrom_coords)) )then
            write(logfhandle,'(A)') '>>> FLEX DIFFMAP coordinate_scale_warning=shape_mismatch'
            call flush(logfhandle)
            return
        endif
        if( .not.all(ieee_is_finite(spectral_z)) .or. .not.all(ieee_is_finite(nystrom_coords)) )then
            write(logfhandle,'(A)') '>>> FLEX DIFFMAP coordinate_scale_warning=nonfinite_values'
            call flush(logfhandle)
            return
        endif
        spec_rms=sqrt(sum(real(spectral_z,dp)**2)/real(size(spectral_z),dp))
        nys_rms=sqrt(sum(real(nystrom_coords,dp)**2)/real(size(nystrom_coords),dp))
        ratio=0._dp
        if( nys_rms>TINY ) ratio=spec_rms/nys_rms
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX DIFFMAP coordinate_scales spectral_rms=', &
            &spec_rms,' nystrom_rms=',nys_rms,' spectral_over_nystrom=',ratio
        call flush(logfhandle)
    end subroutine log_flex_coordinate_scales

    subroutine cleanup_distributed_analysis_parts( params, nparts )
        type(parameters), intent(in) :: params
        integer, intent(in) :: nparts
        type(string) :: fname
        integer :: ipart
        do ipart=1,nparts
            fname=string('flex_residual_features_part')//int2str_pad(ipart,params%numlen)//'.mrcs'
            call del_file(fname); call fname%kill
            fname=string('flex_registered_particle_map_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            call del_file(fname); call fname%kill
            fname=string('flex_graph_neighbors_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            call del_file(fname); call fname%kill
        end do
        call del_file('flex_mean_projections.mrcs')
    end subroutine cleanup_distributed_analysis_parts

    subroutine finish_analysis_outputs( registered_stack, registered_project )
        type(string), intent(inout) :: registered_stack,registered_project
        write(logfhandle,'(A,A)') '>>> FLEX DIFFMAP registered_stack=',registered_stack%to_char()
        write(logfhandle,'(A,A)') '>>> FLEX DIFFMAP registered_project=',registered_project%to_char()
        call flush(logfhandle)
        call registered_stack%kill; call registered_project%kill
    end subroutine finish_analysis_outputs

    subroutine validate_inputs( params, cline, max_modes )
        class(parameters), intent(inout) :: params
        class(cmdline), intent(inout) :: cline
        integer, intent(out) :: max_modes
        if( params%oritype/='ptcl3D' ) THROW_HARD('flex_analysis requires oritype=ptcl3D')
        if( .not.cline%defined('vol1') ) THROW_HARD('flex_analysis requires vol1=<mean map>')
        if( params%nstates/=1 ) THROW_HARD('flex_analysis supports one ptcl3D state')
        if( params%npreimages<2 ) THROW_HARD('flex_analysis requires npreimages>=2')
        if( params%neigs<1 ) THROW_HARD('flex_analysis requires neigs>=1')
        max_modes=params%neigs
        params%k_nn=max(1,params%k_nn)
        params%nang_nbrs=max(params%k_nn,params%nang_nbrs)
        params%bandwidth_mode=lowercase(trim(params%bandwidth_mode))
        if( params%bandwidth_mode/='median' .and. params%bandwidth_mode/='ferguson' ) &
            &THROW_HARD('flex_analysis bandwidth_mode must be median|ferguson')
        params%bandwidth_tune=max(params%bandwidth_tune,0.)
    end subroutine validate_inputs

    subroutine select_particles( params, build, pinds, nptcls )
        class(parameters), intent(in) :: params
        class(builder), intent(inout) :: build
        integer, allocatable, intent(out) :: pinds(:)
        integer, intent(out) :: nptcls
        call build%spproj_field%sample4rec([params%fromp,params%top],nptcls,pinds)
        if( nptcls<3 ) THROW_HARD('flex_analysis found fewer than three active particles')
    end subroutine select_particles

    subroutine validate_selected_particles( spproj, pinds )
        type(sp_project), intent(inout) :: spproj
        integer, intent(in) :: pinds(:)
        integer, allocatable :: stk_sizes(:),stk_offsets(:)
        logical, allocatable :: seen_slots(:)
        integer :: nptcls,nstks,nslots,istk,i,iptcl,stkind,indstk,slot,ctf_mode
        nptcls=spproj%os_ptcl3D%get_noris()
        if( spproj%os_ptcl2D%get_noris()/=nptcls ) THROW_HARD('flex_analysis requires matching ptcl2D/ptcl3D fields')
        nstks=spproj%os_stk%get_noris()
        if( nstks<1 ) THROW_HARD('flex_analysis project has no particle stacks')
        allocate(stk_sizes(nstks),stk_offsets(nstks),source=0)
        nslots=0
        do istk=1,nstks
            stk_offsets(istk)=nslots
            if( spproj%os_stk%isthere(istk,'nptcls_stk') )then
                stk_sizes(istk)=spproj%os_stk%get_int(istk,'nptcls_stk')
            else if( spproj%os_stk%isthere(istk,'fromp').and.spproj%os_stk%isthere(istk,'top') )then
                stk_sizes(istk)=spproj%os_stk%get_top(istk)-spproj%os_stk%get_fromp(istk)+1
            endif
            if( stk_sizes(istk)<1 ) THROW_HARD('flex_analysis requires valid stack particle counts')
            nslots=nslots+stk_sizes(istk)
        end do
        allocate(seen_slots(nslots),source=.false.)
        ctf_mode=spproj%get_ctfflag_type('ptcl3D',pinds(1))
        if( ctf_mode/=CTFFLAG_YES .and. ctf_mode/=CTFFLAG_FLIP ) &
            &THROW_HARD('flex_analysis requires ctf=yes or ctf=flip particle images')
        do i=1,size(pinds)
            iptcl=pinds(i)
            if( iptcl<1 .or. iptcl>nptcls ) THROW_HARD('selected flex particle index outside project')
            if( spproj%os_ptcl3D%get_state(iptcl)<=0 ) THROW_HARD('selected flex particle has state zero')
            if( .not.spproj%os_ptcl3D%isthere(iptcl,'proj') ) THROW_HARD('selected flex particle has no projection index')
            call spproj%map_ptcl_ind2stk_ind('ptcl3D',iptcl,stkind,indstk)
            if( stkind<1 .or. stkind>nstks ) THROW_HARD('selected flex particle has invalid stack index')
            if( indstk<1 .or. indstk>stk_sizes(stkind) ) THROW_HARD('selected flex particle has invalid image index')
            slot=stk_offsets(stkind)+indstk
            if( seen_slots(slot) ) THROW_HARD('selected flex particles alias the same stack image')
            seen_slots(slot)=.true.
            if( spproj%get_ctfflag_type('ptcl3D',iptcl)/=ctf_mode ) &
                &THROW_HARD('flex_analysis does not support mixed ctf=yes and ctf=flip particle stacks')
        end do
        deallocate(stk_sizes,stk_offsets,seen_slots)
    end subroutine validate_selected_particles

    subroutine capture_result( fit, pinds, coords, eigvals, nptcls, nmodes )
        type(flex_embedding_result), intent(inout) :: fit
        integer, intent(in) :: nptcls,nmodes,pinds(nptcls)
        real, intent(in) :: coords(nptcls,nmodes),eigvals(:)
        call fit%kill
        fit%method='diffusion_maps'; fit%nptcls=nptcls; fit%ncomp=nmodes
        allocate(fit%pinds(nptcls),source=pinds)
        allocate(fit%z(nptcls,nmodes),source=real(coords,dp))
        allocate(fit%eigvals(nmodes),source=real(eigvals(1:nmodes),dp))
    end subroutine capture_result

    subroutine write_coordinates( pinds, coords, nptcls, nmodes )
        integer, intent(in) :: nptcls,nmodes,pinds(nptcls)
        real, intent(in) :: coords(nptcls,nmodes)
        integer :: u,i,q
        call del_file('flex_diffmap_coordinates.txt')
        open(newunit=u,file='flex_diffmap_coordinates.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle'
        do q=1,nmodes; write(u,'(A,I0)',advance='no') ' psi',q; end do
        write(u,*)
        do i=1,nptcls
            write(u,'(I10)',advance='no') pinds(i)
            do q=1,nmodes; write(u,'(1X,ES16.8)',advance='no') coords(i,q); end do
            write(u,*)
        end do
        close(u)
    end subroutine write_coordinates

    subroutine write_preimage_particle_map( spproj, pinds, coords, labels, weights )
        type(sp_project), intent(inout) :: spproj
        integer, intent(in) :: pinds(:), labels(:)
        real, intent(in) :: coords(:,:), weights(:,:)
        integer :: u,i,q,s,iptcl,iproj
        real :: hard_weight
        if( size(pinds)<1 ) THROW_HARD('cannot write an empty flex preimage particle map')
        if( size(labels)/=size(pinds) .or. size(coords,1)/=size(pinds) .or. size(weights,1)/=size(pinds) ) &
            &THROW_HARD('flex preimage particle-map row mismatch')
        if( size(coords,2)<1 .or. size(weights,2)<1 ) THROW_HARD('invalid flex preimage particle-map columns')
        call del_file('flex_registered_particle_preimage_map.txt')
        open(newunit=u,file='flex_registered_particle_preimage_map.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# feature_row raw_particle_index projection_index preimage_index preimage_weight'
        do q=1,size(coords,2); write(u,'(A,I0)',advance='no') ' psi',q; end do
        do s=1,size(weights,2); write(u,'(A,I3.3)',advance='no') ' preimage_weight_',s; end do
        write(u,*)
        do i=1,size(pinds)
            iptcl=pinds(i)
            if( iptcl<1 .or. iptcl>spproj%os_ptcl3D%get_noris() ) THROW_HARD('flex preimage particle index outside project')
            if( labels(i)<1 .or. labels(i)>size(weights,2) ) THROW_HARD('invalid flex preimage index')
            iproj=0
            if( spproj%os_ptcl3D%isthere(iptcl,'proj') ) iproj=spproj%os_ptcl3D%get_int(iptcl,'proj')
            hard_weight=weights(i,labels(i))
            write(u,'(I10,1X,I10,1X,I8,1X,I8,1X,ES16.8)',advance='no') i,iptcl,iproj,labels(i),hard_weight
            do q=1,size(coords,2); write(u,'(1X,ES16.8)',advance='no') coords(i,q); end do
            do s=1,size(weights,2); write(u,'(1X,ES16.8)',advance='no') weights(i,s); end do
            write(u,*)
        end do
        close(u)
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE particle map: flex_registered_particle_preimage_map.txt'
        call flush(logfhandle)
    end subroutine write_preimage_particle_map

    subroutine write_spectrum( eigvals, nmodes, icm_enabled, converged, niters, score )
        real, intent(in) :: eigvals(:),score
        integer, intent(in) :: nmodes,niters
        logical, intent(in) :: icm_enabled,converged
        integer :: u,q
        call del_file('flex_diffmap_spectrum.txt')
        open(newunit=u,file='flex_diffmap_spectrum.txt',status='replace',action='write')
        write(u,'(A,I0,A,L1,A,L1,A,I0,A,ES16.8)') '# selected=',nmodes,' icm=',icm_enabled, &
            &' converged=',converged,' iterations=',niters,' score=',score
        write(u,'(A)') '# mode eigenvalue selected'
        do q=1,size(eigvals)
            write(u,'(I6,1X,ES16.8,1X,L1)') q,eigvals(q),q<=nmodes
        end do
        close(u)
    end subroutine write_spectrum

    subroutine write_graph_summary( graph, cmin, cmax, cmean, params )
        type(diffmap_graph), intent(in) :: graph
        integer, intent(in) :: cmin,cmax
        real, intent(in) :: cmean
        class(parameters), intent(in) :: params
        integer :: u
        call del_file('flex_diffmap_graph.txt')
        open(newunit=u,file='flex_diffmap_graph.txt',status='replace',action='write')
        write(u,'(A,I0)') 'particles=',graph%n
        write(u,'(A,I0)') 'directed_edges=',graph%nnz
        write(u,'(A,I0)') 'k_nn=',graph%k_nn
        write(u,'(A,I0)') 'nang_nbrs=',params%nang_nbrs
        write(u,'(A,I0)') 'candidates_min=',cmin
        write(u,'(A,I0)') 'candidates_max=',cmax
        write(u,'(A,F12.3)') 'candidates_mean=',cmean
        close(u)
    end subroutine write_graph_summary

end module simple_flex_analysis_strategy
