!@descr: sparse diffusion-map 3D variability analysis from fixed particle poses
module simple_flex_eigenvol_strategy
use simple_core_module_api
use simple_builder,               only: builder
use simple_cmdline,               only: cmdline
use simple_diff_map_denoise,      only: select_spectral_rank_icm
use simple_diff_map_graphs,       only: diffmap_graph, build_gated_euclidean_knn_graph
use simple_diffusion_maps,        only: embed_graph
use simple_flex_diffmap_features, only: prepare_flex_diffmap_features
use simple_flex_diffmap_rec3D,    only: reconstruct_flex_diffmap_modes, write_flex_diffmap_rec_part, &
    &reduce_flex_diffmap_rec_parts
use simple_parameters,            only: parameters
use simple_flex_embedding_result,  only: flex_embedding_result
use simple_qsys_env,              only: qsys_env
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_eigenvol_diffmap, run_flex_eigenvol_reconstruct_worker

integer, parameter :: FLEX_DIFFMAP_STATE_MAGIC=1180061702, FLEX_DIFFMAP_STATE_VERSION=1

contains

    subroutine run_flex_eigenvol_diffmap( params, build, cline, fit_result )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(flex_embedding_result), optional, intent(inout) :: fit_result
        type(diffmap_graph) :: graph
        type(string) :: registered_stack, registered_project
        integer, allocatable :: pinds(:), proj_ids(:)
        real, allocatable :: features(:,:), proj_dirs(:,:), coords_mode_major(:,:), eigvals(:)
        real, allocatable :: coords(:,:)
        integer :: nptcls, max_modes, nmodes, icm_iters
        integer :: cand_min, cand_max
        real :: cand_mean, icm_score
        logical :: icm_converged
        integer(timer_int_kind) :: t_total, t_step
        t_total = tic()
        call validate_inputs(params, cline, max_modes)
        call select_particles(params, build, pinds, nptcls)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,F8.3)') &
            '>>> FLEX DIFFMAP particles=',nptcls,' k_nn=',params%k_nn,' nang_nbrs=',params%nang_nbrs, &
            &' max_modes=',max_modes,' lp=',params%lp
        if( params%nparts > 1 )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX DIFFMAP parts=',params%nparts, &
                &' (registered preparation and graph are master-owned; Hermitian statistics are distributed)'
        endif
        call flush(logfhandle)

        t_step = tic()
        call prepare_flex_diffmap_features(params, build, pinds, features, proj_ids, proj_dirs, &
            &registered_stack, registered_project)
        write(logfhandle,'(A,F10.3,A,I0)') '>>> FLEX DIFFMAP registration_seconds=',toc(t_step), &
            &' feature_values=',size(features,kind=8)
        call flush(logfhandle)

        t_step = tic()
        call build_gated_euclidean_knn_graph(features, proj_ids, proj_dirs, params%k_nn, params%nang_nbrs, graph, &
            &cand_min, cand_max, cand_mean)
        write(logfhandle,'(A,I0,A,I0,A,F8.1,A,I0)') '>>> FLEX DIFFMAP graph candidates_min=',cand_min, &
            &' max=',cand_max,' mean=',cand_mean,' directed_nnz=',graph%nnz
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP graph_seconds=',toc(t_step)
        call flush(logfhandle)
        deallocate(features,proj_dirs,proj_ids)

        t_step = tic()
        call embed_graph(graph,max_modes,coords_mode_major,eigvals)
        if( size(eigvals) < 1 ) THROW_HARD('diffusion embedding returned no nontrivial modes')
        call select_spectral_rank_icm(eigvals,size(eigvals),nmodes,icm_converged,icm_iters,icm_score,min_rank=1)
        nmodes = min(max(1,nmodes),min(max_modes,size(eigvals)))
        allocate(coords(nptcls,nmodes),source=transpose(coords_mode_major(1:nmodes,:)))
        write(logfhandle,'(A,I0,A,L1,A,I0,A,ES12.4,A,F10.3)') '>>> FLEX DIFFMAP selected_modes=',nmodes, &
            &' icm_converged=',icm_converged,' icm_iters=',icm_iters,' icm_score=',icm_score, &
            &' embedding_seconds=',toc(t_step)
        call flush(logfhandle)

        call write_coordinates(pinds,coords,nptcls,nmodes)
        call write_spectrum(eigvals,nmodes,icm_converged,icm_iters,icm_score)
        call write_graph_summary(graph,cand_min,cand_max,cand_mean,params)
        if( present(fit_result) ) call capture_result(fit_result,pinds,coords,eigvals,nptcls,nmodes)

        t_step = tic()
        if( params%nparts > 1 )then
            call distributed_reconstruct(params,build,cline,pinds,coords,nmodes)
        else
            call reconstruct_flex_diffmap_modes(params,build,pinds,coords,nmodes)
        endif
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP reconstruction_seconds=',toc(t_step)
        write(logfhandle,'(A,A)') '>>> FLEX DIFFMAP registered_stack=',registered_stack%to_char()
        write(logfhandle,'(A,A)') '>>> FLEX DIFFMAP registered_project=',registered_project%to_char()
        write(logfhandle,'(A,F10.3)') '>>> FLEX DIFFMAP total_seconds=',toc(t_total)
        call flush(logfhandle)

        call graph%kill()
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(coords) ) deallocate(coords)
        if( allocated(coords_mode_major) ) deallocate(coords_mode_major)
        if( allocated(eigvals) ) deallocate(eigvals)
        call registered_stack%kill
        call registered_project%kill
    end subroutine run_flex_eigenvol_diffmap

    subroutine run_flex_eigenvol_reconstruct_worker( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, allocatable :: pinds(:)
        real, allocatable :: coords(:,:)
        integer :: nptcls,nmodes
        if( .not.cline%defined('infile') .or. .not.cline%defined('outfile') )then
            THROW_HARD('flex_eigenvol_reconstruct requires infile and outfile')
        endif
        call read_state(params%infile,pinds,coords,nptcls,nmodes)
        call write_flex_diffmap_rec_part(params,build,pinds,coords,nmodes,[params%fromp,params%top],params%outfile)
        deallocate(pinds,coords)
    end subroutine run_flex_eigenvol_reconstruct_worker

    subroutine distributed_reconstruct( params, build, cline, pinds, coords, nmodes )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        class(cmdline), intent(inout) :: cline
        integer, intent(in) :: pinds(:),nmodes
        real, intent(in) :: coords(:,:)
        type(qsys_env) :: qenv
        type(chash) :: job_descr
        type(chash), allocatable :: part_params(:)
        type(string), allocatable :: part_files(:)
        type(string) :: state_file,prg_string
        integer :: ipart,nparts_eff
        nparts_eff=min(max(1,params%nparts),size(pinds))
        state_file='flex_diffmap_reconstruct_state.dat'
        call write_state(state_file,pinds,coords,nmodes)
        allocate(part_params(nparts_eff),part_files(nparts_eff))
        do ipart=1,nparts_eff
            part_files(ipart)=string('flex_diffmap_reconstruct_part')//int2str_pad(ipart,params%numlen)//'.bin'
            call del_file(part_files(ipart))
            call part_params(ipart)%new(1)
            call part_params(ipart)%set('outfile',part_files(ipart)%to_char())
        end do
        prg_string='flex_eigenvol_reconstruct'
        call cline%gen_job_descr(job_descr,prg=prg_string)
        call job_descr%set('mkdir','no')
        call job_descr%set('infile',state_file%to_char())
        call job_descr%set('nparts',int2str(nparts_eff))
        call job_descr%set('numlen',int2str(params%numlen))
        call qenv%new(params,nparts_eff,numlen=params%numlen,nptcls=size(pinds))
        call qenv%gen_scripts_and_schedule_jobs(job_descr,part_params=part_params,array=L_USE_SLURM_ARR,extra_params=params)
        call qenv%kill
        call reduce_flex_diffmap_rec_parts(params,build,part_files,nmodes)
        do ipart=1,nparts_eff
            call del_file(part_files(ipart)); call part_files(ipart)%kill; call part_params(ipart)%kill
        end do
        call del_file(state_file)
        call state_file%kill; call prg_string%kill; call job_descr%kill
        deallocate(part_files,part_params)
    end subroutine distributed_reconstruct

    subroutine write_state( fname, pinds, coords, nmodes )
        class(string), intent(in) :: fname
        integer, intent(in) :: pinds(:),nmodes
        real, intent(in) :: coords(:,:)
        integer :: u,ios,header(4)
        header=[FLEX_DIFFMAP_STATE_MAGIC,FLEX_DIFFMAP_STATE_VERSION,size(pinds),nmodes]
        call del_file(fname)
        call fopen(u,file=fname,access='STREAM',status='REPLACE',action='WRITE',iostat=ios)
        call fileiochk('open flex diffusion state',ios)
        write(u,iostat=ios) header,pinds,coords(:,1:nmodes)
        call fileiochk('write flex diffusion state',ios); call fclose(u)
    end subroutine write_state

    subroutine read_state( fname, pinds, coords, nptcls, nmodes )
        class(string), intent(in) :: fname
        integer, allocatable, intent(out) :: pinds(:)
        real, allocatable, intent(out) :: coords(:,:)
        integer, intent(out) :: nptcls,nmodes
        integer :: u,ios,header(4)
        call fopen(u,file=fname,access='STREAM',status='OLD',action='READ',iostat=ios)
        call fileiochk('open flex diffusion state worker',ios)
        read(u,iostat=ios) header
        call fileiochk('read flex diffusion state header',ios)
        if( header(1)/=FLEX_DIFFMAP_STATE_MAGIC .or. header(2)/=FLEX_DIFFMAP_STATE_VERSION )then
            THROW_HARD('incompatible flex diffusion state')
        endif
        nptcls=header(3); nmodes=header(4)
        allocate(pinds(nptcls),coords(nptcls,nmodes))
        read(u,iostat=ios) pinds,coords
        call fileiochk('read flex diffusion state payload',ios); call fclose(u)
    end subroutine read_state

    subroutine validate_inputs( params, cline, max_modes )
        class(parameters), intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        integer,           intent(out)   :: max_modes
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_eigenvol requires oritype=ptcl3D')
        if( .not. cline%defined('vol1') ) THROW_HARD('flex_eigenvol requires vol1=<mean map>')
        if( params%nstates /= 1 ) THROW_HARD('flex_eigenvol supports one ptcl3D state')
        max_modes = min(20,max(1,params%neigs))
        params%k_nn = max(1,params%k_nn)
        params%nang_nbrs = max(params%k_nn,params%nang_nbrs)
    end subroutine validate_inputs

    subroutine select_particles( params, build, pinds, nptcls )
        class(parameters), intent(in)       :: params
        class(builder),    intent(inout)    :: build
        integer, allocatable, intent(inout) :: pinds(:)
        integer, intent(out) :: nptcls
        call build%spproj_field%sample4rec([params%fromp,params%top],nptcls,pinds)
        if( nptcls < 3 ) THROW_HARD('flex_eigenvol found fewer than three active particles')
    end subroutine select_particles

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

    subroutine write_spectrum( eigvals, nmodes, converged, niters, score )
        real, intent(in) :: eigvals(:),score
        integer, intent(in) :: nmodes,niters
        logical, intent(in) :: converged
        integer :: u,q
        call del_file('flex_diffmap_spectrum.txt')
        open(newunit=u,file='flex_diffmap_spectrum.txt',status='replace',action='write')
        write(u,'(A,I0,A,L1,A,I0,A,ES16.8)') '# selected=',nmodes,' converged=',converged,' iterations=',niters,' score=',score
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

end module simple_flex_eigenvol_strategy
