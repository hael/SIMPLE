!@descr: SIMPLE-native linear 3D variability analysis from fixed particle poses
module simple_flex_eigenvol_strategy
use simple_core_module_api
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_projected_latent_model, only: update_basis_from_latents, infer_latents_from_basis, &
    &initialize_latents, orthonormalize_latents, latent_sdev, latent_covariance, &
    &basis_fourier_energy, cleanup_planes, projected_model_kfromto, &
    &write_mstep_flexvol_part_file, update_basis_from_flexvol_part_files, &
    &write_estep_latent_part_file, reduce_estep_latent_part_files
use simple_qsys_env,        only: qsys_env
use simple_reconstructor,    only: reconstructor
implicit none

public :: run_flex_eigenvol_linear, run_flex_eigenvol_mstep_worker, run_flex_eigenvol_estep_worker
private
#include "simple_local_flags.inc"

real,     parameter :: TRAJ_SIGMAS(7) = [-3., -2., -1., 0., 1., 2., 3.]
integer,  parameter :: FLEX_STATE_MAGIC = 1180054611
integer,  parameter :: FLEX_STATE_VERSION = 1

contains

    subroutine run_flex_eigenvol_linear( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(reconstructor)              :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type),  allocatable  :: fpls(:)
        type(string)                     :: sigma2_fname, zfname
        integer, allocatable             :: pinds(:)
        real(dp), allocatable            :: z(:,:), z_postcov(:,:,:), resid_energy(:), resid_mean_energy(:), mode_vars(:)
        logical                          :: l_distributed
        integer                          :: nptcls, ncomp, niters, iter, q, kfromto_eff(2), nparts_eff
        integer(timer_int_kind)          :: t_total, t_step
        t_total = tic()
        call validate_inputs(params, cline, ncomp, niters)
        call select_particles(params, build, pinds, nptcls)
        nparts_eff = max(1, min(max(1, params%nparts), nptcls))
        if( nparts_eff /= params%nparts )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL REDUCING NPARTS FROM ', params%nparts, ' TO ', nparts_eff
            params%nparts = nparts_eff
            params%numlen = len(int2str(params%nparts))
            call cline%set('nparts', params%nparts)
            call cline%set('numlen', params%numlen)
        endif
        l_distributed = params%nparts > 1
        allocate(z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp), resid_energy(nptcls), resid_mean_energy(nptcls), &
            &mode_vars(ncomp), basis_recs(ncomp))
        call initialize_latents(z, nptcls, ncomp)
        call orthonormalize_latents(z, nptcls, ncomp)
        z_postcov = 0.d0
        mode_vars = 1.d0
        kfromto_eff = projected_model_kfromto(params)
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL PARTICLES  : ', nptcls
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL COMPONENTS : ', ncomp
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL PPCA ITERS : ', niters
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL THREADS    : ', nthr_glob
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL PARTS      : ', params%nparts
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL BOX        : ', params%box_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL SMPD       : ', params%smpd_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL MASK PX    : ', params%msk_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL BASIS LP A : ', params%lp
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL KFROMTO    : ', kfromto_eff(1), ' / ', kfromto_eff(2)
        write(logfhandle,'(A,A)')  '>>> FLEX_EIGENVOL PTCL SRC   : ', trim(params%ptcl_src)
        write(logfhandle,'(A,A)')  '>>> FLEX_EIGENVOL MEAN MAP  : ', params%vols(1)%to_char()
        call flush(logfhandle)
        t_step = tic()
        sigma2_fname = string('flex_eigenvol_sigma2.tmp')
        call build%esig%new(params, build%pftc, sigma2_fname, params%box_crop)
        call init_mean_reconstructor(params, build, mean_rec)
        do q = 1, ncomp
            call init_output_reconstructor(params, build, basis_recs(q))
        end do
        call log_seconds('>>> FLEX_EIGENVOL SETUP SECONDS', toc(t_step))
        do iter = 1, niters
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL LINEAR ITERATION ', iter, ' / ', niters
            call flush(logfhandle)
            t_step = tic()
            if( l_distributed )then
                call distributed_mstep(params, build, cline, basis_recs, z, z_postcov, &
                    &pinds, nptcls, ncomp, iter)
            else
                call update_basis_from_latents(params, build, mean_rec, basis_recs, z, z_postcov, &
                    &pinds, nptcls, ncomp, fpls, log_label='FLEX_EIGENVOL')
            endif
            call log_iter_seconds('>>> FLEX_EIGENVOL ITER M-STEP SECONDS', iter, toc(t_step))
            t_step = tic()
            if( l_distributed )then
                call distributed_estep(params, build, cline, basis_recs, z, mode_vars, &
                    &z_postcov, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, iter)
            else
                call infer_latents_from_basis(params, build, mean_rec, basis_recs, z, mode_vars, &
                    &z_postcov, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls, log_label='FLEX_EIGENVOL')
            endif
            call log_iter_seconds('>>> FLEX_EIGENVOL ITER E-STEP SECONDS', iter, toc(t_step))
            call log_residual_stats('>>> FLEX_EIGENVOL ITER MEAN-ONLY RESIDUAL ENERGY', resid_mean_energy, nptcls)
            call log_residual_stats('>>> FLEX_EIGENVOL ITER MODE RESIDUAL ENERGY', resid_energy, nptcls)
            call log_residual_reduction('>>> FLEX_EIGENVOL ITER RESIDUAL REDUCTION', resid_mean_energy, resid_energy, nptcls)
            call log_latent_sdevs('>>> FLEX_EIGENVOL ITER RAW LATENT SD', z, nptcls, ncomp)
            call log_mode_vars('>>> FLEX_EIGENVOL ITER PPCA MODE VAR', mode_vars, ncomp)
            call orthonormalize_latents(z, nptcls, ncomp)
            call log_latent_sdevs('>>> FLEX_EIGENVOL ITER NORM LATENT SD', z, nptcls, ncomp)
        end do
        write(logfhandle,'(A)') '>>> FLEX_EIGENVOL FINAL REFIT'
        call flush(logfhandle)
        t_step = tic()
        if( l_distributed )then
            call distributed_mstep(params, build, cline, basis_recs, z, z_postcov, &
                &pinds, nptcls, ncomp, niters + 1)
        else
            call update_basis_from_latents(params, build, mean_rec, basis_recs, z, z_postcov, &
                &pinds, nptcls, ncomp, fpls, log_label='FLEX_EIGENVOL')
        endif
        call log_seconds('>>> FLEX_EIGENVOL FINAL M-STEP SECONDS', toc(t_step))
        t_step = tic()
        if( l_distributed )then
            call distributed_estep(params, build, cline, basis_recs, z, mode_vars, &
                &z_postcov, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, niters + 1)
        else
            call infer_latents_from_basis(params, build, mean_rec, basis_recs, z, mode_vars, &
                &z_postcov, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls, log_label='FLEX_EIGENVOL')
        endif
        call log_seconds('>>> FLEX_EIGENVOL FINAL E-STEP SECONDS', toc(t_step))
        call log_residual_stats('>>> FLEX_EIGENVOL FINAL MEAN-ONLY RESIDUAL ENERGY', resid_mean_energy, nptcls)
        call log_residual_stats('>>> FLEX_EIGENVOL FINAL MODE RESIDUAL ENERGY', resid_energy, nptcls)
        call log_residual_reduction('>>> FLEX_EIGENVOL FINAL RESIDUAL REDUCTION', resid_mean_energy, resid_energy, nptcls)
        call log_latent_sdevs('>>> FLEX_EIGENVOL RAW FINAL LATENT SD', z, nptcls, ncomp)
        call log_latent_means('>>> FLEX_EIGENVOL RAW FINAL LATENT MEAN', z, nptcls, ncomp)
        call log_mode_vars('>>> FLEX_EIGENVOL FINAL PPCA MODE VAR', mode_vars, ncomp)
        call log_latent_covariance_eigs('>>> FLEX_EIGENVOL FINAL LATENT COV EIGVAL', z, nptcls, ncomp)
        call orthonormalize_latents(z, nptcls, ncomp)
        call log_latent_sdevs('>>> FLEX_EIGENVOL OUTPUT LATENT SD', z, nptcls, ncomp)
        call log_latent_means('>>> FLEX_EIGENVOL OUTPUT LATENT MEAN', z, nptcls, ncomp)
        call log_basis_fourier_norms('>>> FLEX_EIGENVOL OUTPUT BASIS FOURIER NORM', basis_recs, ncomp)
        zfname = output_prefix(params)//'_zcoords.txt'
        t_step = tic()
        call write_basis_outputs(params, basis_recs, ncomp)
        call write_particle_coordinates(zfname, pinds, z, resid_energy, nptcls, ncomp)
        call write_trajectory_volumes(params, basis_recs, z, nptcls, ncomp)
        call log_seconds('>>> FLEX_EIGENVOL OUTPUT SECONDS', toc(t_step))
        write(logfhandle,'(A,A)') '>>> FLEX_EIGENVOL WROTE Z TABLE : ', zfname%to_char()
        call log_seconds('>>> FLEX_EIGENVOL TOTAL SECONDS', toc(t_total))
        call flush(logfhandle)
        call cleanup_planes(fpls)
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(z) ) deallocate(z)
        if( allocated(z_postcov) ) deallocate(z_postcov)
        if( allocated(resid_energy) ) deallocate(resid_energy)
        if( allocated(resid_mean_energy) ) deallocate(resid_mean_energy)
        if( allocated(mode_vars) ) deallocate(mode_vars)
        call build%esig%kill
        call mean_rec%dealloc_rho
        call mean_rec%kill
        if( allocated(basis_recs) )then
            do q = 1, size(basis_recs)
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        endif
        call sigma2_fname%kill
        call zfname%kill
    end subroutine run_flex_eigenvol_linear

    subroutine run_flex_eigenvol_mstep_worker( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(reconstructor)              :: mean_rec
        type(fplane_type), allocatable   :: fpls(:)
        type(string)                     :: sigma2_fname
        integer, allocatable             :: pinds(:)
        real(dp), allocatable            :: z(:,:), z_postcov(:,:,:), resid_energy(:), resid_mean_energy(:), mode_vars(:)
        integer                          :: nptcls, ncomp, niters, partlims(2)
        call validate_inputs(params, cline, ncomp, niters)
        if( .not. cline%defined('infile') ) THROW_HARD('flex_eigenvol_mstep worker requires infile=<state.dat>')
        if( .not. cline%defined('outfile') ) THROW_HARD('flex_eigenvol_mstep worker requires outfile=<part.bin>')
        call read_flex_state_file(params%infile, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy, nptcls, ncomp)
        params%nptcls = nptcls
        partlims = [params%fromp, params%top]
        sigma2_fname = 'flex_eigenvol_sigma2_mstep_part'//int2str_pad(params%part, params%numlen)//'.tmp'
        call build%esig%new(params, build%pftc, sigma2_fname, params%box_crop)
        call init_mean_reconstructor(params, build, mean_rec)
        call write_mstep_flexvol_part_file(params, build, mean_rec, z, z_postcov, pinds, nptcls, ncomp, &
            &partlims, params%outfile, fpls, log_label='FLEX_EIGENVOL_WORKER')
        call cleanup_planes(fpls)
        call cleanup_worker_state(mean_rec, sigma2_fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy)
    end subroutine run_flex_eigenvol_mstep_worker

    subroutine run_flex_eigenvol_estep_worker( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(reconstructor)              :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type), allocatable   :: fpls(:)
        type(string)                     :: sigma2_fname
        integer, allocatable             :: pinds(:)
        real(dp), allocatable            :: z(:,:), z_postcov(:,:,:), resid_energy(:), resid_mean_energy(:), mode_vars(:)
        integer                          :: nptcls, ncomp, niters, partlims(2), q
        call validate_inputs(params, cline, ncomp, niters)
        if( .not. cline%defined('infile') ) THROW_HARD('flex_eigenvol_estep worker requires infile=<state.dat>')
        if( .not. cline%defined('outfile') ) THROW_HARD('flex_eigenvol_estep worker requires outfile=<part.dat>')
        if( .not. cline%defined('refs') ) THROW_HARD('flex_eigenvol_estep worker requires refs=<basis-prefix>')
        call read_flex_state_file(params%infile, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy, nptcls, ncomp)
        params%nptcls = nptcls
        allocate(basis_recs(ncomp))
        partlims = [params%fromp, params%top]
        sigma2_fname = 'flex_eigenvol_sigma2_estep_part'//int2str_pad(params%part, params%numlen)//'.tmp'
        call build%esig%new(params, build%pftc, sigma2_fname, params%box_crop)
        call init_mean_reconstructor(params, build, mean_rec)
        call init_basis_reconstructors_from_prefix(params, build, basis_recs, ncomp, params%refs)
        call write_estep_latent_part_file(params, build, mean_rec, basis_recs, mode_vars, pinds, nptcls, ncomp, &
            &partlims, params%outfile, fpls, log_label='FLEX_EIGENVOL_WORKER')
        call cleanup_planes(fpls)
        do q = 1, ncomp
            call basis_recs(q)%dealloc_rho
            call basis_recs(q)%kill
        end do
        deallocate(basis_recs)
        call cleanup_worker_state(mean_rec, sigma2_fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy)
    end subroutine run_flex_eigenvol_estep_worker

    subroutine distributed_mstep( params, build, cline, basis_recs, z, z_postcov, pinds, nptcls, ncomp, iter )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        class(cmdline),      intent(inout) :: cline
        integer,             intent(in)    :: nptcls, ncomp, iter
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        type(string), allocatable :: part_fnames(:)
        type(string) :: state_fname
        real(dp), allocatable :: dummy_mode_vars(:), dummy_resid(:), dummy_mean_resid(:)
        allocate(dummy_mode_vars(ncomp), dummy_resid(nptcls), dummy_mean_resid(nptcls), source=0.d0)
        dummy_mode_vars = 1.d0
        state_fname = flex_state_fname(iter)
        call write_flex_state_file(state_fname, pinds, z, z_postcov, &
            &dummy_mode_vars, dummy_resid, dummy_mean_resid, nptcls, ncomp)
        call make_part_fnames('flexvol_mstep', iter, params%nparts, params%numlen, '.bin', part_fnames)
        call dispatch_flex_workers(params, cline, 'flex_eigenvol_mstep', state_fname, part_fnames, nptcls)
        call update_basis_from_flexvol_part_files(params, build, basis_recs, ncomp, part_fnames, params%nparts, &
            &log_label='FLEX_EIGENVOL')
        call cleanup_string_array(part_fnames)
        deallocate(dummy_mode_vars, dummy_resid, dummy_mean_resid)
        call state_fname%kill
    end subroutine distributed_mstep

    subroutine distributed_estep( params, build, cline, basis_recs, z, mode_vars, z_postcov, resid_energy, &
        &resid_mean_energy, pinds, nptcls, ncomp, iter )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        class(cmdline),      intent(inout) :: cline
        integer,             intent(in)    :: nptcls, ncomp, iter
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(inout) :: z(nptcls,ncomp), mode_vars(ncomp)
        real(dp),            intent(inout) :: z_postcov(nptcls,ncomp,ncomp), resid_energy(nptcls), resid_mean_energy(nptcls)
        integer,             intent(in)    :: pinds(nptcls)
        type(string), allocatable :: part_fnames(:)
        type(string) :: state_fname, basis_prefix
        state_fname  = flex_state_fname(iter)
        basis_prefix = flex_basis_prefix(iter)
        call write_flex_state_file(state_fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy, nptcls, ncomp)
        call write_temp_basis_volumes(params, basis_recs, ncomp, basis_prefix)
        call make_part_fnames('flexvol_estep', iter, params%nparts, params%numlen, '.dat', part_fnames)
        call dispatch_flex_workers(params, cline, 'flex_eigenvol_estep', state_fname, part_fnames, nptcls, basis_prefix)
        call reduce_estep_latent_part_files(part_fnames, params%nparts, z, z_postcov, resid_energy, resid_mean_energy, &
            &mode_vars, nptcls, ncomp, log_label='FLEX_EIGENVOL')
        call cleanup_string_array(part_fnames)
        call state_fname%kill
        call basis_prefix%kill
    end subroutine distributed_estep

    subroutine dispatch_flex_workers( params, cline, worker_prg, state_fname, part_fnames, nptcls, basis_prefix )
        class(parameters), intent(in)    :: params
        class(cmdline),    intent(inout) :: cline
        character(len=*),  intent(in)    :: worker_prg
        class(string),     intent(in)    :: state_fname
        class(string),     intent(in)    :: part_fnames(:)
        integer,           intent(in)    :: nptcls
        class(string), optional, intent(in) :: basis_prefix
        type(qsys_env) :: qenv
        type(chash)    :: job_descr
        type(chash), allocatable :: part_params(:)
        type(string) :: prg_string
        integer :: ipart, nparts_run
        nparts_run = size(part_fnames)
        allocate(part_params(nparts_run))
        do ipart = 1, nparts_run
            call part_params(ipart)%new(3)
            call part_params(ipart)%set('outfile', part_fnames(ipart)%to_char())
            if( present(basis_prefix) ) call part_params(ipart)%set('refs', basis_prefix%to_char())
        end do
        prg_string = worker_prg
        call cline%gen_job_descr(job_descr, prg=prg_string)
        call job_descr%set('mkdir',  'no')
        call job_descr%set('infile', state_fname%to_char())
        call job_descr%set('nparts', int2str(nparts_run))
        call job_descr%set('numlen', int2str(params%numlen))
        call qenv%new(params, nparts_run, numlen=params%numlen, nptcls=nptcls)
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR, extra_params=params)
        call qenv%kill
        call job_descr%kill
        call prg_string%kill
        do ipart = 1, nparts_run
            call part_params(ipart)%kill
        end do
        deallocate(part_params)
    end subroutine dispatch_flex_workers

    subroutine write_flex_state_file( fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy, nptcls, ncomp )
        class(string), intent(in) :: fname
        integer,      intent(in) :: nptcls, ncomp, pinds(nptcls)
        real(dp),     intent(in) :: z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp)
        real(dp),     intent(in) :: mode_vars(ncomp), resid_energy(nptcls), resid_mean_energy(nptcls)
        integer :: funit, io_stat, header(4)
        header = [FLEX_STATE_MAGIC, FLEX_STATE_VERSION, nptcls, ncomp]
        call del_file(fname)
        call fopen(funit, file=fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_flex_state_file; open '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_flex_state_file; header '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) pinds, mode_vars
        call fileiochk('write_flex_state_file; particle/mode fields '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) z, z_postcov, resid_energy, resid_mean_energy
        call fileiochk('write_flex_state_file; latent fields '//fname%to_char(), io_stat)
        call fclose(funit)
    end subroutine write_flex_state_file

    subroutine read_flex_state_file( fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy, nptcls, ncomp )
        class(string), intent(in) :: fname
        integer, allocatable, intent(inout) :: pinds(:)
        real(dp), allocatable, intent(inout) :: z(:,:), z_postcov(:,:,:), mode_vars(:), resid_energy(:), resid_mean_energy(:)
        integer, intent(out) :: nptcls, ncomp
        integer :: funit, io_stat, header(4)
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(z) ) deallocate(z)
        if( allocated(z_postcov) ) deallocate(z_postcov)
        if( allocated(mode_vars) ) deallocate(mode_vars)
        if( allocated(resid_energy) ) deallocate(resid_energy)
        if( allocated(resid_mean_energy) ) deallocate(resid_mean_energy)
        if( .not. file_exists(fname) ) THROW_HARD('missing flex_eigenvol state file: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_flex_state_file; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) header
        call fileiochk('read_flex_state_file; header '//fname%to_char(), io_stat)
        if( header(1) /= FLEX_STATE_MAGIC ) THROW_HARD('bad flex_eigenvol state magic: '//fname%to_char())
        if( header(2) /= FLEX_STATE_VERSION ) THROW_HARD('bad flex_eigenvol state version: '//fname%to_char())
        nptcls = header(3)
        ncomp  = header(4)
        if( nptcls < 1 .or. ncomp < 1 ) THROW_HARD('invalid flex_eigenvol state header: '//fname%to_char())
        allocate(pinds(nptcls), z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp), mode_vars(ncomp), &
            &resid_energy(nptcls), resid_mean_energy(nptcls))
        read(funit, iostat=io_stat) pinds, mode_vars
        call fileiochk('read_flex_state_file; particle/mode fields '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) z, z_postcov, resid_energy, resid_mean_energy
        call fileiochk('read_flex_state_file; latent fields '//fname%to_char(), io_stat)
        call fclose(funit)
    end subroutine read_flex_state_file

    subroutine write_temp_basis_volumes( params, basis_recs, ncomp, prefix )
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        class(string),       intent(in)    :: prefix
        type(string) :: fname
        integer :: q
        do q = 1, ncomp
            fname = flex_basis_fname(prefix, q)
            call write_projection_ready_volume(basis_recs(q), fname)
            call fname%kill
        end do
    end subroutine write_temp_basis_volumes

    subroutine init_basis_reconstructors_from_prefix( params, build, basis_recs, ncomp, prefix )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        class(string),       intent(in)    :: prefix
        type(string) :: fname
        integer :: q
        do q = 1, ncomp
            fname = flex_basis_fname(prefix, q)
            if( .not. file_exists(fname) ) THROW_HARD('missing flex_eigenvol basis volume: '//fname%to_char())
            call basis_recs(q)%read_and_crop(fname, params%smpd_crop, params%box_crop, params%smpd_crop)
            call basis_recs(q)%fft
            call basis_recs(q)%alloc_rho(params, build%spproj, expand=.true.)
            call basis_recs(q)%expand_exp
            call fname%kill
        end do
    end subroutine init_basis_reconstructors_from_prefix

    function flex_state_fname( iter ) result( fname )
        integer, intent(in) :: iter
        type(string) :: fname
        fname = 'flex_eigenvol_state_iter'//int2str_pad(iter, 3)//'.dat'
    end function flex_state_fname

    function flex_basis_prefix( iter ) result( prefix )
        integer, intent(in) :: iter
        type(string) :: prefix
        prefix = 'flex_eigenvol_basis_iter'//int2str_pad(iter, 3)
    end function flex_basis_prefix

    function flex_basis_fname( prefix, icomp ) result( fname )
        class(string), intent(in) :: prefix
        integer,       intent(in) :: icomp
        type(string) :: fname
        fname = prefix//'_comp'//int2str_pad(icomp, 3)//'.mrc'
    end function flex_basis_fname

    subroutine make_part_fnames( fbody, iter, nparts, numlen, ext, fnames )
        character(len=*), intent(in) :: fbody, ext
        integer,          intent(in) :: iter, nparts, numlen
        type(string), allocatable, intent(inout) :: fnames(:)
        integer :: ipart
        call cleanup_string_array(fnames)
        allocate(fnames(nparts))
        do ipart = 1, nparts
            fnames(ipart) = trim(fbody)//'_iter'//int2str_pad(iter, 3)//'_part'//int2str_pad(ipart, numlen)//trim(ext)
            call del_file(fnames(ipart))
        end do
    end subroutine make_part_fnames

    subroutine cleanup_string_array( arr )
        type(string), allocatable, intent(inout) :: arr(:)
        integer :: i
        if( allocated(arr) )then
            do i = 1, size(arr)
                call arr(i)%kill
            end do
            deallocate(arr)
        endif
    end subroutine cleanup_string_array

    subroutine cleanup_worker_state( mean_rec, sigma2_fname, pinds, z, z_postcov, mode_vars, resid_energy, resid_mean_energy )
        type(reconstructor), intent(inout) :: mean_rec
        type(string),        intent(inout) :: sigma2_fname
        integer, allocatable, intent(inout) :: pinds(:)
        real(dp), allocatable, intent(inout) :: z(:,:), z_postcov(:,:,:), mode_vars(:), resid_energy(:), resid_mean_energy(:)
        call mean_rec%dealloc_rho
        call mean_rec%kill
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(z) ) deallocate(z)
        if( allocated(z_postcov) ) deallocate(z_postcov)
        if( allocated(mode_vars) ) deallocate(mode_vars)
        if( allocated(resid_energy) ) deallocate(resid_energy)
        if( allocated(resid_mean_energy) ) deallocate(resid_mean_energy)
        call sigma2_fname%kill
    end subroutine cleanup_worker_state

    subroutine validate_inputs( params, cline, ncomp, niters )
        class(parameters), intent(in)    :: params
        class(cmdline),    intent(inout) :: cline
        integer,           intent(out)   :: ncomp, niters
        if( trim(params%oritype) /= 'ptcl3D' )then
            THROW_HARD('flex_eigenvol currently requires oritype=ptcl3D')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('flex_eigenvol requires vol1=<mean map> from abinitio3D/refine3D')
        endif
        if( params%nstates /= 1 )then
            THROW_HARD('flex_eigenvol currently supports a single ptcl3D state')
        endif
        ncomp  = max(1, params%neigs)
        niters = max(1, params%maxits)
        if( ncomp > 16 )then
            THROW_WARN('flex_eigenvol capping neigs at 16 for this first PPCA implementation')
            ncomp = 16
        endif
    end subroutine validate_inputs

    subroutine select_particles( params, build, pinds, nptcls )
        class(parameters), intent(in)       :: params
        class(builder),    intent(inout)    :: build
        integer, allocatable, intent(inout) :: pinds(:)
        integer, intent(out)                :: nptcls
        call build%spproj_field%sample4rec([params%fromp, params%top], nptcls, pinds)
        if( nptcls <= 0 ) THROW_HARD('flex_eigenvol found no active particles')
    end subroutine select_particles

    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        call mean_rec%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        call mean_rec%fft
        call mean_rec%alloc_rho(params, build%spproj, expand=.true.)
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor

    subroutine init_output_reconstructor( params, build, basis_rec )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: basis_rec
        call basis_rec%new([params%box_crop, params%box_crop, params%box_crop], params%smpd_crop)
        call basis_rec%alloc_rho(params, build%spproj, expand=.true.)
        call basis_rec%reset
        call basis_rec%reset_exp
    end subroutine init_output_reconstructor

    subroutine write_basis_outputs( params, basis_recs, ncomp )
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        type(string) :: fname
        integer      :: q
        do q = 1, ncomp
            fname = basis_volume_name(params, q)
            write(logfhandle,'(A,I0,A,I0,A,A)') '>>> FLEX_EIGENVOL WRITING COMPONENT ', q, ' / ', ncomp, ': ', fname%to_char()
            call flush(logfhandle)
            call write_projection_ready_volume(basis_recs(q), fname)
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_EIGENVOL WROTE COMPONENT ', q, ': ', fname%to_char()
            call flush(logfhandle)
            call fname%kill
        end do
    end subroutine write_basis_outputs

    subroutine write_projection_ready_volume( basis_rec, fname )
        type(reconstructor), intent(inout) :: basis_rec
        class(string),       intent(in)    :: fname
        type(image) :: img
        call img%copy(basis_rec)
        if( img%is_ft() ) call img%ifft
        call log_image_stats('>>> FLEX_EIGENVOL COMPONENT VOLUME STATS', img)
        call img%write(fname, del_if_exists=.true.)
        call img%kill
    end subroutine write_projection_ready_volume

    subroutine write_particle_coordinates( zfname, pinds, z, resid_energy, nptcls, ncomp )
        class(string), intent(in) :: zfname
        integer,      intent(in) :: nptcls, ncomp
        integer,      intent(in) :: pinds(nptcls)
        real(dp),     intent(in) :: z(nptcls,ncomp), resid_energy(nptcls)
        integer :: funit, ierr, i, q
        call del_file(zfname)
        call fopen(funit, file=zfname, status='NEW', action='WRITE', iostat=ierr)
        call fileiochk('flex_eigenvol opening '//zfname%to_char(), ierr)
        write(funit,'(A)', advance='no') '# particle'
        do q = 1, ncomp
            write(funit,'(A,I0)', advance='no') ' z', q
        end do
        write(funit,'(A)') ' residual_energy'
        do i = 1, nptcls
            write(funit,'(I10)', advance='no') pinds(i)
            do q = 1, ncomp
                write(funit,'(1X,ES14.6)', advance='no') z(i,q)
            end do
            write(funit,'(1X,ES14.6)') resid_energy(i)
        end do
        call fclose(funit)
    end subroutine write_particle_coordinates

    subroutine write_trajectory_volumes( params, basis_recs, z, nptcls, ncomp )
        class(parameters),   intent(inout) :: params
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp)
        type(image)  :: mean_img, basis_img, traj_img, diff_img
        type(string) :: fname
        real         :: sig
        integer      :: q, isig
        call mean_img%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        call log_image_stats('>>> FLEX_EIGENVOL MEAN VOLUME STATS', mean_img)
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL WRITING TRAJECTORY COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            call basis_img%copy(basis_recs(q))
            if( basis_img%is_ft() ) call basis_img%ifft
            call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY BASIS STATS', basis_img)
            sig = latent_sdev(z(:,q), nptcls)
            if( sig <= TINY ) sig = 1.
            write(logfhandle,'(A,I0,A,ES12.4)') '>>> FLEX_EIGENVOL TRAJECTORY LATENT SD ', q, ': ', sig
            call flush(logfhandle)
            do isig = 1, size(TRAJ_SIGMAS)
                call traj_img%copy(mean_img)
                call traj_img%add(basis_img, TRAJ_SIGMAS(isig) * sig)
                call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY VOLUME STATS', traj_img)
                call diff_img%copy(traj_img)
                call diff_img%subtr(mean_img)
                call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY MINUS MEAN STATS', diff_img)
                fname = trajectory_diff_volume_name(params, q, TRAJ_SIGMAS(isig))
                call diff_img%write(fname, del_if_exists=.true.)
                call fname%kill
                call diff_img%kill
                fname = trajectory_volume_name(params, q, TRAJ_SIGMAS(isig))
                call traj_img%write(fname, del_if_exists=.true.)
                call traj_img%kill
                call fname%kill
            end do
            call basis_img%kill
            write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL WROTE TRAJECTORY COMPONENT ', q
            call flush(logfhandle)
        end do
        call mean_img%kill
    end subroutine write_trajectory_volumes

    subroutine log_latent_covariance_eigs( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        real(dp), allocatable :: cov(:,:), rot(:,:), eigvals(:)
        integer :: nrot, q
        if( ncomp <= 0 ) return
        allocate(cov(ncomp,ncomp), rot(ncomp,ncomp), eigvals(ncomp), source=0.d0)
        call latent_covariance(z, nptcls, ncomp, cov)
        nrot = 0
        call jacobi(cov, ncomp, ncomp, eigvals, rot, nrot)
        call eigsrt(eigvals, rot, ncomp, ncomp)
        eigvals = max(eigvals, 0.d0)
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', eigvals(q)
        end do
        call flush(logfhandle)
        deallocate(cov, rot, eigvals)
    end subroutine log_latent_covariance_eigs

    subroutine log_seconds( label, seconds )
        character(len=*),      intent(in) :: label
        real(timer_int_kind),  intent(in) :: seconds
        write(logfhandle,'(A,A,F10.3)') trim(label), ': ', seconds
        call flush(logfhandle)
    end subroutine log_seconds

    subroutine log_iter_seconds( label, iter, seconds )
        character(len=*),     intent(in) :: label
        integer,              intent(in) :: iter
        real(timer_int_kind), intent(in) :: seconds
        write(logfhandle,'(A,1X,I0,A,F10.3)') trim(label), iter, ': ', seconds
        call flush(logfhandle)
    end subroutine log_iter_seconds

    subroutine log_comp_seconds( label, icomp, seconds )
        character(len=*),     intent(in) :: label
        integer,              intent(in) :: icomp
        real(timer_int_kind), intent(in) :: seconds
        write(logfhandle,'(A,1X,I0,A,F10.3)') trim(label), icomp, ': ', seconds
        call flush(logfhandle)
    end subroutine log_comp_seconds

    subroutine log_residual_stats( label, resid_energy, nptcls )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls
        real(dp),         intent(in) :: resid_energy(nptcls)
        real(dp) :: avg, sd
        avg = sum(resid_energy) / real(max(1,nptcls), dp)
        sd  = sqrt(sum((resid_energy - avg) * (resid_energy - avg)) / real(max(1,nptcls - 1), dp))
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean=', avg, ' sd=', sd, ' min=', minval(resid_energy), ' max=', maxval(resid_energy)
        call flush(logfhandle)
    end subroutine log_residual_stats

    subroutine log_residual_reduction( label, resid_mean_energy, resid_energy, nptcls )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls
        real(dp),         intent(in) :: resid_mean_energy(nptcls), resid_energy(nptcls)
        real(dp) :: mean_only, with_modes, reduction
        mean_only  = sum(resid_mean_energy) / real(max(1,nptcls), dp)
        with_modes = sum(resid_energy)      / real(max(1,nptcls), dp)
        reduction  = 0.d0
        if( mean_only > DTINY ) reduction = (mean_only - with_modes) / mean_only
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean_only=', mean_only, ' with_modes=', with_modes, ' frac=', reduction
        call flush(logfhandle)
    end subroutine log_residual_reduction

    subroutine log_image_stats( label, img )
        character(len=*), intent(in)    :: label
        type(image),      intent(inout) :: img
        real :: ave, sdev, maxv, minv
        logical :: err
        if( img%is_ft() ) call img%ifft
        call img%stats(ave, sdev, maxv, minv, errout=err)
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean=', ave, ' sd=', sdev, ' min=', minv, ' max=', maxv, ' range=', maxv - minv
        call flush(logfhandle)
    end subroutine log_image_stats

    subroutine log_latent_sdevs( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', latent_sdev(z(:,q), nptcls)
        end do
        call flush(logfhandle)
    end subroutine log_latent_sdevs

    subroutine log_latent_means( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', &
                &sum(z(:,q)) / real(max(1,nptcls), dp)
        end do
        call flush(logfhandle)
    end subroutine log_latent_means

    subroutine log_mode_vars( label, mode_vars, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: ncomp
        real(dp),         intent(in) :: mode_vars(ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', mode_vars(q)
        end do
        call flush(logfhandle)
    end subroutine log_mode_vars

    subroutine log_basis_fourier_norms( label, basis_recs, ncomp )
        character(len=*),   intent(in) :: label
        integer,            intent(in) :: ncomp
        type(reconstructor), intent(in) :: basis_recs(ncomp)
        real(dp) :: energy
        integer  :: q
        do q = 1, ncomp
            energy = basis_fourier_energy(basis_recs(q))
            write(logfhandle,'(A,1X,I0,A,ES12.4,A,ES12.4)') trim(label), q, &
                &': norm=', sqrt(max(0.d0, energy)), ' energy=', energy
        end do
        call flush(logfhandle)
    end subroutine log_basis_fourier_norms

    function output_prefix( params ) result( prefix )
        class(parameters), intent(in) :: params
        type(string) :: prefix
        type(string) :: ext
        character(len=:), allocatable :: stem
        ext    = fname2ext(params%outvol)
        prefix = get_fbody(params%outvol, ext)
        stem   = prefix%to_char()
        if( len_trim(stem) >= 4 )then
            if( stem(len_trim(stem)-3:len_trim(stem)) == '_001' ) stem = stem(:len_trim(stem)-4)
        endif
        prefix = string(trim(stem))
        call ext%kill
    end function output_prefix

    function basis_volume_name( params, icomp ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        type(string) :: fname, prefix
        character(len=3) :: tag
        if( icomp == 1 )then
            fname = params%outvol
            return
        endif
        write(tag,'(I3.3)') icomp
        prefix = output_prefix(params)
        fname  = prefix//'_'//tag//'.mrc'
        call prefix%kill
    end function basis_volume_name

    function trajectory_volume_name( params, icomp, sigma ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        real,              intent(in) :: sigma
        type(string) :: fname, prefix
        character(len=3) :: ctag
        character(len=2) :: stag
        write(ctag,'(I3.3)') icomp
        if( sigma < 0. )then
            write(stag,'(A,I1)') 'm', nint(abs(sigma))
        else if( sigma > 0. )then
            write(stag,'(A,I1)') 'p', nint(abs(sigma))
        else
            stag = 'z0'
        endif
        prefix = output_prefix(params)
        fname  = prefix//'_mode'//ctag//'_'//stag//'.mrc'
        call prefix%kill
    end function trajectory_volume_name

    function trajectory_diff_volume_name( params, icomp, sigma ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        real,              intent(in) :: sigma
        type(string) :: fname, prefix
        character(len=3) :: ctag
        character(len=2) :: stag
        write(ctag,'(I3.3)') icomp
        if( sigma < 0. )then
            write(stag,'(A,I1)') 'm', nint(abs(sigma))
        else if( sigma > 0. )then
            write(stag,'(A,I1)') 'p', nint(abs(sigma))
        else
            stag = 'z0'
        endif
        prefix = output_prefix(params)
        fname  = prefix//'_mode'//ctag//'_'//stag//'_diff.mrc'
        call prefix%kill
    end function trajectory_diff_volume_name

end module simple_flex_eigenvol_strategy
