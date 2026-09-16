!@descr: ab initio 2D analysis
module simple_commanders_abinitio2D
use simple_commanders_api
use simple_commanders_cavgs
use simple_commanders_cluster2D
use simple_abinitio2D_controller 
use simple_gui_communicator,       only: gui_communicator
implicit none

public :: commander_abinitio2D, execute_abinitio2D_staged
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio2D
    contains
    procedure :: execute => exec_abinitio2D
end type commander_abinitio2D

! class variables
type(stage_params), allocatable :: stage_parms(:)

abstract interface
    subroutine terminal_cline_policy( cline )
        import :: cmdline
        class(cmdline), intent(inout) :: cline
    end subroutine terminal_cline_policy
end interface

contains

    subroutine exec_abinitio2D( self, cline )
        class(commander_abinitio2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call exec_abinitio2D_workflow(cline, 1, 0, 0, .false.)
    end subroutine exec_abinitio2D

    subroutine execute_abinitio2D_staged( cline, start_stage_requested, stop_stage_requested, &
        &checkpoint_last_iter, l_checkpoint, terminal_policy )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: start_stage_requested, stop_stage_requested, checkpoint_last_iter
        logical,        intent(in)    :: l_checkpoint
        procedure(terminal_cline_policy), optional :: terminal_policy
        call exec_abinitio2D_workflow(cline, start_stage_requested, stop_stage_requested, &
            &checkpoint_last_iter, l_checkpoint, terminal_policy)
    end subroutine execute_abinitio2D_staged

    subroutine exec_abinitio2D_workflow( cline, start_stage_requested, stop_stage_requested, &
        &checkpoint_last_iter, l_checkpoint, terminal_policy )
        use simple_classaverager
        use simple_timer, only: timer_int_kind, tic, toc
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: start_stage_requested, stop_stage_requested, checkpoint_last_iter
        logical,        intent(in)    :: l_checkpoint
        procedure(terminal_cline_policy), optional :: terminal_policy
        ! commanders
        type(commander_cluster2D)  :: xcluster2D
        type(commander_calc_pspec) :: xcalc_pspec
        ! command lines
        type(cmdline)              :: cline_cluster2D, cline_calc_pspec
        ! other
        type(parameters)           :: params
        type(sp_project)           :: spproj
        type(gui_communicator)     :: gui_comm
        class(oris),       pointer :: spproj_field
        integer, allocatable :: seed_parent(:), seed_pops(:)
        integer :: maxits, istage, last_iter, nptcls_eff, nstages, nsample_target_2D
        integer :: start_stage, stop_stage, it_pass
        integer(timer_int_kind) :: t_tot, t_phase
        real(timer_int_kind)    :: rt_setup, rt_calc_pspec, rt_cluster2D, rt_final_cavgs, rt_tot
        logical :: l_shmem, l_seeded
        if( L_BENCH_GLOB )then
            rt_setup       = 0.
            rt_calc_pspec  = 0.
            rt_cluster2D   = 0.
            rt_final_cavgs = 0.
            rt_tot         = 0.
            t_tot          = tic()
            t_phase        = t_tot
        endif
        call cline%set('oritype',   'ptcl2D')
        if( .not. cline%defined('sigma_est')     ) call cline%set('sigma_est',     'global')
        if( .not. cline%defined('autoscale')     ) call cline%set('autoscale',     'yes')
        if( .not. cline%defined('mkdir')         ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('center')        ) call cline%set('center',        'yes')
        if( .not. cline%defined('center_type')   ) call cline%set('center_type',   'seg')
        if( .not. cline%defined('cls_init')      ) call cline%set('cls_init',      'rand')
        if( .not. cline%defined('gauref')        ) call cline%set('gauref',        'yes')
        if( .not. cline%defined('extr_lim')      ) call cline%set('extr_lim',      EXTR_LIM_LOCAL)
        if( .not. cline%defined('nits_per_stage')) call cline%set('nits_per_stage',ITS_INCR)
        if( .not. cline%defined('eo_stage')      ) call cline%set('eo_stage',      EO_STAGE)
        if( .not. cline%defined('rank_cavgs')    ) call cline%set('rank_cavgs',    'yes')
        if( .not. cline%defined('stats')         ) call cline%set('stats',         'no')
        if( .not. cline%defined('refine')        ) call cline%set('refine',        'prob_snhc')
        if( .not. cline%defined('ml_reg')        ) call cline%set('ml_reg',        'yes')
        ! shared memory execution
        l_shmem = set_shmem_flag(cline)
        ! master parameters
        call params%new(cline)
        call gui_comm%new(params)
        if( params%l_nonuniform ) THROW_HARD('2D nonuniform filtering has been removed; exec_abinitio2D')
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        maxits = params%extr_lim
        call cline%delete('stats')
        ! check refinement flag and set stages
        call determine_abinitio2D_stages(params, nstages)
        ! override # stages
        if( cline%defined('nstages') ) nstages = min(params%nstages,NSTAGES_CLS)
        start_stage = start_stage_requested
        ! seeded restart from the previous 2D clustering: the run is entered at the
        ! first probabilistic stage, preceded by one all-particle seed pass
        l_seeded = trim(params%cls_init) == 'prev'
        if( l_seeded )then
            if( l_checkpoint ) THROW_HARD('cls_init=prev is not supported with stream checkpointing')
            if( start_stage_requested /= 1 ) THROW_HARD('cls_init=prev requires a fresh abinitio2D run')
            if( nstages < PROBREFINE_STAGE ) THROW_HARD('cls_init=prev requires nstages >= 3')
            start_stage = PROBREFINE_STAGE
        endif
        if( stop_stage_requested > 0 )then
            stop_stage = stop_stage_requested
        else
            stop_stage = nstages
        endif
        if( l_checkpoint )then
            if( start_stage < 1 .or. start_stage > nstages ) THROW_HARD('invalid stream checkpoint start stage')
            if( stop_stage < start_stage .or. stop_stage > nstages ) THROW_HARD('invalid stream checkpoint stop stage')
            if( start_stage > 1 .and. checkpoint_last_iter < 1 )&
                &THROW_HARD('stream checkpoint continuation requires the last completed iteration')
        endif
        allocate(stage_parms(nstages))
        ! read project
        call spproj%read(params%projfile)
        call set_dims                   ! set downscaling
        if( start_stage == 1 ) call inirefs ! deal with initial references only on fresh runs
        call set_lplims(nstages)        ! set resolutions limits
        call prep_command_lines(cline)  ! prepare class command lines
        if( start_stage > 1 .and. .not. l_seeded )then
            call cline_cluster2D%set('endit', checkpoint_last_iter)
            write(logfhandle,'(A,I0,A,I0)') '>>> ABINITIO2D CHECKPOINT RESUME: start_stage=', start_stage,&
                &' last_iter=', checkpoint_last_iter
        endif
        call set_sampling               ! sampling
        if( l_seeded )then
            ! the seed pass is the iteration at which the last pre-probabilistic
            ! stage would have ended; the stages from PROBREFINE_STAGE on are unchanged
            it_pass = abinitio2D_seed_pass_iter(maxits, stage_parms(1)%update_frac)
            call cline_cluster2D%set('endit', it_pass)
        endif
        if( L_BENCH_GLOB ) rt_setup = toc(t_phase)
        ! summary
        do istage = 1,nstages
            write(logfhandle,'(A,I2,A,L1,F6.1,2I8,F7.3,2L2)')'>>> STAGE ', istage,' LPSET LP MAXCLSPOP NPTCLS UFRAC ST FR: ',&
            &stage_parms(istage)%l_lpset,stage_parms(istage)%lp, stage_parms(istage)%max_cls_pop,&
            &stage_parms(istage)%nptcls, stage_parms(istage)%update_frac,&
            &stage_parms(istage)%l_sticky_sampling, stage_parms(istage)%l_frac_restore
        end do
        ! prep particles field
        if( start_stage == 1 )then
            call spproj_field%delete_2Dclustering
            call spproj_field%clean_entry('updatecnt', 'sampled')
            if( spproj_field%get_nevenodd() == 0 ) call spproj_field%partition_eo
            call spproj%write_segment_inside(params%oritype, params%projfile)
        else
            if( l_seeded ) call ensure_seed_eo ! before the sigma2 state, which is built on the e/o halves
            call ensure_resume_sigma_state
        endif
        if( l_seeded )then
            call seed_from_previous_clustering ! partition + seed references (needs the sigma2 state above)
            call execute_seed_pass
        endif

        ! Frequency marching
        do istage = start_stage,stop_stage
            write(logfhandle,'(A)')'>>>'
            if( stage_parms(istage)%l_lpset )then
                write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', stage_parms(istage)%lp
            else
                write(logfhandle,'(A,I3,A)')'>>> STAGE ', istage,' WITH GOLD STANDARD E/O'
            endif
            ! parameters update
            call set_cline_cluster2D_stage(cline_cluster2D, cline, params, stage_parms, maxits, istage)
            ! classify
            call execute_cluster2D
            ! update GUI
            call spproj%read_segment('cls2D', params%projfile)
            call spproj%read_segment('out',   params%projfile)
            call gui_comm%add_metadata(spproj, oritype='cls2D', stage2D=istage)
        enddo
        if( l_checkpoint .and. stop_stage < nstages )then
            last_iter = cline_cluster2D%get_iarg('endit')
            write(logfhandle,'(A,I0,A,I0)') '>>> ABINITIO2D CHECKPOINT READY: stage=', stop_stage,&
                &' last_iter=', last_iter
            call cline_cluster2D%kill
            call cline_calc_pspec%kill
            deallocate(stage_parms)
            call spproj%kill
            nullify(spproj_field)
            call qsys_cleanup(params)
            call simple_touch('ABINITIO2D_CHECKPOINT_STAGE'//int2str_pad(stop_stage,3))
            call simple_end('**** SIMPLE_ABINITIO2D CHECKPOINT NORMAL STOP ****')
            return
        endif
        call execute_terminal_pass( terminal_policy=terminal_policy )
        ! transfer 2D shifts to 3D field only when no prior valid 3D alignment exists
        call spproj%read_segment(params%oritype,params%projfile)
        call spproj%read_segment('ptcl3D',params%projfile)
        if( spproj%has_valid_ptcl3D_alignment() )then
            write(logfhandle,'(A)') '>>> ABINITIO2D: prior ptcl3D alignment detected; preserving ptcl3D shifts'
        else
            call spproj%os_ptcl3D%transfer_2Dshifts(spproj_field)
            call spproj%write_segment_inside('ptcl3D', params%projfile)
        endif
        ! weights & final mapping of particles
        if( trim(params%stats).eq.'yes' ) call output_stats('final')
        ! final class generation & ranking
        last_iter = cline_cluster2D%get_iarg('endit')
        if( L_BENCH_GLOB ) t_phase = tic()
        call gen_final_cavgs(last_iter)
        if( L_BENCH_GLOB )then
            rt_final_cavgs = toc(t_phase)
            call write_abinitio_benchmark(last_iter + 1, 'final_cavgs', nstages)
        endif
        if( l_seeded ) call write_seed_lineage(with_final_pops=.true.)
        ! final update GUI
        call spproj%read_segment('cls2D', params%projfile)
        call spproj%read_segment('out',   params%projfile)
        call gui_comm%add_metadata(spproj, oritype='cls2D', stage2D=0, selection=.true.) ! stage2D=0 signifies final
        ! cleanup
        call del_file('start2Drefs'//params%ext%to_char())
        call del_file('start2Drefs_even'//params%ext%to_char())
        call del_file('start2Drefs_odd'//params%ext%to_char())
        call del_files(DIST_FBODY,      params%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params%nparts,ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        call cline_cluster2D%kill
        call cline_calc_pspec%kill
        deallocate(stage_parms)
        if( allocated(seed_parent) ) deallocate(seed_parent)
        if( allocated(seed_pops)   ) deallocate(seed_pops)
        call spproj%kill
        nullify(spproj_field)
        call qsys_cleanup(params)
        call simple_touch(ABINITIO2D_FINISHED)
        call gui_comm%kill()
        call simple_end('**** SIMPLE_ABINITIO2D NORMAL STOP ****')
        
      contains

        ! Downscaling/cropping dimensions used throughout
        subroutine set_dims
            real :: smpd_target_eff, scale_factor
            if( cline%defined('box_crop') )then
                scale_factor = real(params%box_crop) / real(params%box)
                if( .not.cline%defined('smpd_crop') ) params%smpd_crop = params%smpd / scale_factor
                if( .not.cline%defined('msk_crop')  ) params%msk_crop  = round2even(params%msk * scale_factor)
                params%l_autoscale = params%box_crop < params%box
            else
                smpd_target_eff  = max(SMPD_TARGET, params%smpd)
                scale_factor     = 1.0
                params%smpd_crop = params%smpd
                params%box_crop  = params%box
                params%msk_crop  = params%msk
                if( params%l_autoscale .and. params%box >= MINBOXSZ )then
                    call autoscale(params%box, params%smpd, smpd_target_eff, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOXSZ)
                    params%l_autoscale = params%box_crop < params%box
                endif
                if( params%l_autoscale ) params%msk_crop = round2even(params%msk * scale_factor)
            endif
            if( params%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
            stage_parms(:)%box_crop  = params%box_crop
            stage_parms(:)%smpd_crop = params%smpd_crop
            stage_parms(:)%trslim    = min(5.,max(2.0, AHELIX_WIDTH/params%smpd_crop))
            if( cline%defined('trs') ) stage_parms(2:)%trslim = params%trs
        end subroutine set_dims

        ! Deals with initial references dimensions when *not* abinitio
        subroutine inirefs
            use simple_procimgstk,        only: copy_imgfile
            use simple_commanders_imgops, only: commander_scale
            use simple_commanders_volops, only: commander_noisevol
            type(commander_scale)    :: xscale
            type(cmdline)            :: cline_scalerefs
            type(string) :: refs, refs_even, refs_odd
            real         :: smpd
            integer      :: ldim(3), ncls
            logical      :: eo
            if( .not.cline%defined('refs') )then
                ! Starting from scratch: appropriate starting references
                ! will be generated by cluster2D/cluster2D_distr
                return
            endif
            refs = params%refs
            call find_ldim_nptcls(refs, ldim, ncls)
            smpd = params%smpd
            ldim(3) = 1
            if( .not.file_exists(refs) ) THROW_HARD('File does not exits: '//refs%to_char())
            if( ncls /= params%ncls )    THROW_HARD('Incompatible # of classes in: '//refs%to_char())
            refs_even        = add2fbody(refs, params%ext, '_even')
            refs_odd         = add2fbody(refs, params%ext, '_odd')
            eo               = file_exists(refs_even).and.file_exists(refs_odd)
            params%refs      = 'start2Drefs'//params%ext%to_char()    ! initial references
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
            if( ldim(1) == stage_parms(1)%box_crop )then
                call copy_imgfile(refs, params%refs, stage_parms(1)%smpd_crop, [1,params%ncls])
                if( eo )then
                    call copy_imgfile(refs_even, params%refs_even, stage_parms(1)%smpd_crop, [1,params%ncls])
                    call copy_imgfile(refs_odd,  params%refs_odd,  stage_parms(1)%smpd_crop, [1,params%ncls])
                endif
            else
                call cline_scalerefs%set('stk',    refs)
                call cline_scalerefs%set('outstk', params%refs)
                call cline_scalerefs%set('smpd',   smpd)
                call cline_scalerefs%set('newbox', stage_parms(1)%box_crop)
                call cline_scalerefs%set('nthr',   nthr_glob)
                call xscale%execute(cline_scalerefs)
                if( eo )then
                    call cline_scalerefs%set('stk',    refs_even)
                    call cline_scalerefs%set('outstk', params%refs_even)
                    call xscale%execute(cline_scalerefs)
                    call cline_scalerefs%set('stk',    refs_odd)
                    call cline_scalerefs%set('outstk', params%refs_odd)
                    call xscale%execute(cline_scalerefs)
                endif
                call cline_scalerefs%kill
            endif
            call refs%kill
            call refs_even%kill
            call refs_odd%kill
        end subroutine inirefs

        ! cls_init=prev: the previous even/odd partition is input state and is
        ! kept; one is generated only when the project has none. Note that
        ! isthere('eo') is false for eo=0 (zero-valued particle parameters count
        ! as absent), so only the field-level get_nevenodd test is meaningful.
        subroutine ensure_seed_eo
            if( spproj_field%get_nevenodd() > 0 ) return
            write(logfhandle,'(A)') '>>> ABINITIO2D SEED WARNING: no even/odd partition in the project; generating one'
            call spproj_field%partition_eo
            call spproj%write_segment_inside(params%oritype, params%projfile)
        end subroutine ensure_seed_eo

        ! cls_init=prev: seed partition and references from the previous 2D
        ! clustering, metadata only (class, state, corr; cls2D state). The only
        ! hard error is the absence of a previous clustering; everything else
        ! is repaired with a warning. Seed references are made from the labels.
        subroutine seed_from_previous_clustering
            type(commander_make_cavgs_distr) :: xmake_cavgs_distr
            type(commander_make_cavgs)       :: xmake_cavgs
            type(cmdline)                    :: cline_make_cavgs
            integer, allocatable :: cls_states(:), pops(:), clsinds(:), tmpinds(:)
            integer :: nptcls, iptcl, icls, ncls_prev, ncls_sel, nactive, nlabelled, ndropped
            integer :: nrej, nfloor, nunassigned
            nptcls    = spproj_field%get_noris()
            nactive   = 0
            nlabelled = 0
            do iptcl = 1, nptcls
                if( spproj_field%get_state(iptcl) <= 0 ) cycle
                nactive = nactive + 1
                if( spproj_field%get_class(iptcl) >= 1 ) nlabelled = nlabelled + 1
            end do
            if( nactive   == 0 ) THROW_HARD('cls_init=prev: no active particles in the project')
            if( nlabelled == 0 ) THROW_HARD('cls_init=prev requires a previous 2D clustering; no active particle carries a class label')
            if( nlabelled < nactive )then
                write(logfhandle,'(A,I0,A)') '>>> ABINITIO2D SEED WARNING: ', nactive - nlabelled,&
                    &' active particles without a class label; they enter the seed pass unassigned'
            endif
            ! accepted parents: cls2D state (when present) and the population floor
            ncls_prev = spproj_field%get_n('class')
            allocate(cls_states(ncls_prev), source=1)
            ncls_sel = spproj%os_cls2D%get_noris()
            if( ncls_sel > 0 .and. spproj%os_cls2D%isthere('state') )then
                do icls = 1, min(ncls_prev, ncls_sel)
                    cls_states(icls) = spproj%os_cls2D%get_state(icls)
                end do
                if( ncls_sel < ncls_prev )then
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> ABINITIO2D SEED WARNING: particle labels reach class ',&
                        &ncls_prev, ' but cls2D holds ', ncls_sel, ' entries; classes beyond it are accepted'
                endif
            else
                write(logfhandle,'(A)') '>>> ABINITIO2D SEED WARNING: no cls2D selection state in the project; every labelled class is accepted'
            endif
            call spproj_field%get_pops(pops, 'class', maxn=ncls_prev)
            nrej   = count(cls_states == 0 .and. pops > 0)
            nfloor = count(cls_states >  0 .and. pops > 0 .and. pops < MINCLSPOPLIM)
            tmpinds = (/(icls, icls=1,ncls_prev)/)
            clsinds = pack(tmpinds, mask=(cls_states > 0 .and. pops >= MINCLSPOPLIM))
            if( size(clsinds) == 0 ) THROW_HARD('cls_init=prev: no accepted class holds enough active particles (MINCLSPOPLIM)')
            ! the seed partition
            call spproj_field%reseed_classes(clsinds, params%ncls, seed_parent, seed_pops, ndropped)
            call spproj_field%clean_entry('updatecnt', 'sampled')
            call spproj%write_segment_inside(params%oritype, params%projfile)
            nunassigned = 0
            do iptcl = 1, nptcls
                if( spproj_field%get_state(iptcl) <= 0 ) cycle
                if( spproj_field%get_class(iptcl) < 1 ) nunassigned = nunassigned + 1
            end do
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A,I0,A,I0,A)') '>>> ABINITIO2D SEED: ', size(clsinds), ' accepted parent classes -> ',&
                &params%ncls, ' seed classes'
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> ABINITIO2D SEED: rejected classes ', nrej,&
                &', classes below the population floor ', nfloor, ', parents dropped (ncls < parents) ', ndropped
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> ABINITIO2D SEED: unassigned particles ', nunassigned,&
                &', seed populations min/median/max ', minval(seed_pops), '/', nint(median(real(seed_pops))), '/', maxval(seed_pops)
            call write_seed_lineage(with_final_pops=.false.)
            ! seed references from the labels, at the working scale, same names as inirefs
            params%refs      = 'start2Drefs'//params%ext%to_char()
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
            cline_make_cavgs = cline
            call cline_make_cavgs%delete('ptcl_src')
            call cline_make_cavgs%delete('autoscale')
            call cline_make_cavgs%delete('balance')
            call cline_make_cavgs%set('prg',       'make_cavgs')
            call cline_make_cavgs%set('refs',      params%refs)
            call cline_make_cavgs%set('box_crop',  stage_parms(1)%box_crop)
            call cline_make_cavgs%set('smpd_crop', stage_parms(1)%smpd_crop)
            call cline_make_cavgs%set('ml_reg',    'yes')
            if( l_shmem )then
                call xmake_cavgs%execute(cline_make_cavgs)
            else
                call xmake_cavgs_distr%execute(cline_make_cavgs)
            endif
            call cline_make_cavgs%kill
            deallocate(cls_states, pops, clsinds, tmpinds)
        end subroutine seed_from_previous_clustering

        ! seed_lineage.txt: seed class -> parent class, seed population and,
        ! after the final class generation, the final population
        subroutine write_seed_lineage( with_final_pops )
            logical, intent(in) :: with_final_pops
            type(sp_project)     :: lineage_proj
            integer, allocatable :: final_pops(:)
            integer :: fnr, iseed
            if( .not. allocated(seed_parent) ) return
            if( with_final_pops )then
                call lineage_proj%read_segment(params%oritype, params%projfile)
                call lineage_proj%os_ptcl2D%get_pops(final_pops, 'class', maxn=params%ncls)
                call lineage_proj%kill
            else
                allocate(final_pops(params%ncls), source=-1)
            endif
            call fopen(fnr, FILE=string('seed_lineage.txt'), STATUS='REPLACE', action='WRITE')
            write(fnr,'(A)') '# seed_class parent_class seed_pop final_pop(-1 = not yet run)'
            do iseed = 1, params%ncls
                write(fnr,'(I8,1X,I8,1X,I8,1X,I8)') iseed, seed_parent(iseed), seed_pops(iseed), final_pops(iseed)
            end do
            call fclose(fnr)
            deallocate(final_pops)
        end subroutine write_seed_lineage

        ! One dense probabilistic iteration of every active particle against the
        ! seed references at the limit of the first probabilistic stage; the
        ! stages from PROBREFINE_STAGE on then run exactly as on an unseeded run
        subroutine execute_seed_pass
            type(cmdline)    :: cline_pass
            type(sp_project) :: pass_proj
            type(oris)       :: os_before, os_after
            integer, allocatable :: nretained(:)
            real    :: sh_before(2), sh_after(2)
            integer :: iptcl, nptcls, n_active, n_class, n_shift, n_seed_ok, cls_b, cls_a
            write(logfhandle,'(A)') '>>>'
            if( stage_parms(PROBREFINE_STAGE)%l_lpset )then
                write(logfhandle,'(A,I0,A,F5.1)') '>>> ABINITIO2D SEED PASS: refine=prob, all particles, iteration ',&
                    &it_pass, ', lp = ', stage_parms(PROBREFINE_STAGE)%lp
            else
                write(logfhandle,'(A,I0)') '>>> ABINITIO2D SEED PASS: refine=prob, all particles, iteration ', it_pass
            endif
            call os_before%copy(spproj_field)
            cline_pass = cline_cluster2D
            call set_cline_cluster2D_seed_pass(cline_pass, params, stage_parms, it_pass, params%refs%to_char())
            call del_file(CLUSTER2D_FINISHED)
            if( L_BENCH_GLOB )then
                rt_calc_pspec = 0.
                rt_cluster2D  = 0.
                t_phase       = tic()
            endif
            call xcluster2D%execute(cline_pass)
            call cline_cluster2D%set('endit', it_pass)
            if( L_BENCH_GLOB )then
                rt_cluster2D = toc(t_phase)
                call write_abinitio_benchmark(it_pass, 'seed_pass', 0)
            endif
            call cline_pass%kill
            ! reassignment diagnostic against the seed partition
            call pass_proj%read_segment(params%oritype, params%projfile)
            call os_after%copy(pass_proj%os_ptcl2D)
            call pass_proj%kill
            nptcls   = os_after%get_noris()
            n_active = 0
            n_class  = 0
            n_shift  = 0
            allocate(nretained(params%ncls), source=0)
            do iptcl = 1, nptcls
                if( os_after%get_state(iptcl) < 1 ) cycle
                n_active = n_active + 1
                cls_b = os_before%get_class(iptcl)
                cls_a = os_after%get_class(iptcl)
                if( cls_a /= cls_b ) n_class = n_class + 1
                if( cls_b >= 1 .and. cls_b <= params%ncls .and. cls_a == cls_b ) nretained(cls_b) = nretained(cls_b) + 1
                sh_before = os_before%get_2Dshift(iptcl)
                sh_after  = os_after%get_2Dshift(iptcl)
                if( sqrt(sum((sh_after - sh_before)**2)) > 1.0 ) n_shift = n_shift + 1
            end do
            n_seed_ok = count(seed_pops > 0 .and. 2 * nretained >= seed_pops)
            if( n_active > 0 )then
                write(logfhandle,'(A,F6.2,A,F6.2,A,F6.2,A)') '>>> ABINITIO2D SEED PASS REASSIGNED: ',&
                    &100. * real(n_class) / real(n_active), ' % changed class, ',&
                    &100. * real(n_shift) / real(n_active), ' % moved shift > 1 px, ',&
                    &100. * real(n_seed_ok) / real(params%ncls), ' % of seed classes retained >= 50% of their members'
            endif
            deallocate(nretained)
            call os_before%kill
            call os_after%kill
        end subroutine execute_seed_pass

        ! Set resolution limits
        subroutine set_lplims( local_nstages )
            use simple_class_frcs, only: class_frcs
            integer, intent(in) :: local_nstages
            real    :: lpstart, lpstop, cenlp
            integer :: istage
            ! Resolution limits
            call mskdiam2lplimits_cluster2D(params%mskdiam, lpstart, lpstop, cenlp)
            lpstart = max(lpstart, 2.*params%smpd_crop)
            lpstop  = max(lpstop,  2.*params%smpd_crop)
            cenlp   = max(cenlp,   2.*params%smpd_crop)
            ! Stages resolution limits
            if( cline%defined('lp') )then
                ! Keep the Gaussian-reference stage coarse, then fix lp when ML regularization starts.
                stage_parms(:)%lp      = params%lp
                stage_parms(:)%l_lpset = .true.
                stage_parms(1)%lp      = lpstart
                params%lpstart = lpstart
                params%lpstop  = params%lp
            else
                ! Frequency marching
                if( .not. cline%defined('lpstart') ) params%lpstart = lpstart
                if( .not. cline%defined('lpstop')  ) params%lpstop  = lpstop
                if( trim(params%eo_stage)=='yes')then
                    if( local_nstages > 1 )then
                        stage_parms(1)%lp      = params%lpstart
                        stage_parms(1)%l_lpset = .true.
                        do istage = 2, local_nstages-1
                            stage_parms(istage)%lp      = stage_parms(istage-1)%lp - (stage_parms(istage-1)%lp - params%lpstop)/2.0
                            stage_parms(istage)%l_lpset = .true.
                        end do
                        stage_parms(local_nstages-1)%lp = params%lpstop
                    endif
                    stage_parms(local_nstages)%l_lpset = .false.
                    stage_parms(local_nstages)%lp      = params%lpstop
                else
                    stage_parms(1)%lp      = params%lpstart
                    stage_parms(:)%l_lpset = .true.
                    do istage = 2, local_nstages
                        stage_parms(istage)%lp      = stage_parms(istage-1)%lp - (stage_parms(istage-1)%lp - params%lpstop)/2.0
                    end do
                    stage_parms(local_nstages)%lp = params%lpstop
                endif
            endif
            if( .not. cline%defined('cenlp') ) params%cenlp   = cenlp
            write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', params%lpstart
            write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', params%lpstop
            write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params%cenlp
        end subroutine set_lplims

        subroutine set_sampling
            nptcls_eff = spproj%count_state_gt_zero()
            call set_abinitio2D_sampling_policy(params, stage_parms, nstages, nptcls_eff, nsample_target_2D)
        end subroutine set_sampling

        subroutine prep_command_lines( cline )
            class(cmdline), intent(in) :: cline
            cline_cluster2D  = cline
            cline_calc_pspec = cline
            call cline_cluster2D%delete('ptcl_src')
            call cline_calc_pspec%delete('ptcl_src')
            ! initial sigma2
            call cline_calc_pspec%set('prg',      'calc_pspec')
            ! cluster2D
            call cline_cluster2D%set('prg',       'cluster2D')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call cline_cluster2D%set('chunk',     'no')
            call set_automask2D_defaults( cline_cluster2D )
        end subroutine prep_command_lines

        subroutine execute_cluster2D( bench_phase, bench_stage )
            character(len=*), optional, intent(in) :: bench_phase
            integer,          optional, intent(in) :: bench_stage
            character(len=:), allocatable :: phase_label
            integer :: stage_label
            if( L_BENCH_GLOB )then
                rt_calc_pspec = 0.
                rt_cluster2D  = 0.
            endif
            phase_label = 'stage'
            stage_label = istage
            if( present(bench_phase) ) phase_label = bench_phase
            if( present(bench_stage) ) stage_label = bench_stage
            call del_file(CLUSTER2D_FINISHED)
            ! Initial sigma2
            if( istage == 1 )then
                if( L_BENCH_GLOB ) t_phase = tic()
                call xcalc_pspec%execute(cline_calc_pspec)
                if( L_BENCH_GLOB ) rt_calc_pspec = toc(t_phase)
            endif
            ! clustering
            if( L_BENCH_GLOB ) t_phase = tic()
            call xcluster2D%execute(cline_cluster2D)
            if( L_BENCH_GLOB )then
                rt_cluster2D = toc(t_phase)
                call write_abinitio_benchmark(cline_cluster2D%get_iarg('endit'), phase_label, stage_label)
            endif
            if( allocated(phase_label) ) deallocate(phase_label)
        end subroutine execute_cluster2D

        subroutine ensure_resume_sigma_state
            use, intrinsic :: iso_fortran_env, only: int64
            use simple_sigma2_state, only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
            use simple_sigma2_state_file, only: sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, &
                &SIGMA2_GROUP_STACK, SIGMA2_STATE_COMMITTED
            type(string) :: state_path
            integer(int64) :: layout_digest
            integer :: iptcl, ngroups, status
            logical :: found, rebuild
            character(len=STDLEN) :: message
            if( params%cc_objfun /= OBJFUN_EUCLID ) return
            rebuild = .true.
            call spproj%get_sigma2_state_path(state_path, found)
            if( found )then
                call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
                if( status == 0 )then
                    layout_digest = sigma2_state_project_layout_digest(spproj, spproj_field)
                    if( params%l_sigma_glob )then
                        call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, &
                            &1, fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                            &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_GLOBAL, &
                            &expected_ngroups=1)
                    else
                        ngroups = 0
                        do iptcl = 1, params%nptcls
                            if( spproj_field%get_state(iptcl) <= 0 ) cycle
                            ngroups = max(ngroups, spproj_field%get_int(iptcl, 'stkind'))
                        enddo
                        call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, &
                            &1, fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                            &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_STACK, &
                            &expected_ngroups=ngroups)
                    endif
                    rebuild = status /= 0
                endif
            endif
            if( rebuild )then
                write(logfhandle,'(A)') '>>> ABINITIO2D CHECKPOINT: rebuilding missing or stale canonical sigma2 state'
                call xcalc_pspec%execute(cline_calc_pspec)
                call spproj%read_segment('projinfo', params%projfile)
            else
                write(logfhandle,'(A)') '>>> ABINITIO2D CHECKPOINT: reusing validated canonical sigma2 state'
            endif
            call state_path%kill
        end subroutine ensure_resume_sigma_state

        subroutine execute_terminal_pass( terminal_policy )
            procedure(terminal_cline_policy), optional :: terminal_policy
            type(string) :: terminal_refs
            integer      :: terminal_start
            if( .not. any(stage_parms(:)%l_update_frac) ) return
            last_iter      = cline_cluster2D%get_iarg('endit')
            terminal_start = last_iter + 1
            terminal_refs  = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//params%ext%to_char()
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I8)')'>>> TERMINAL GREEDY ALL-PARTICLE PASS FROM ITERATION ', last_iter
            call cline_cluster2D%set('refs',          terminal_refs)
            call cline_cluster2D%set('startit',       terminal_start)
            call cline_cluster2D%set('minits',        1)
            call cline_cluster2D%set('maxits',        1)
            call cline_cluster2D%set('extr_iter',     params%extr_lim + 1)
            call cline_cluster2D%set('refine',        'greedy')
            call cline_cluster2D%set('restore_cavgs', 'yes')
            if( present(terminal_policy) ) call terminal_policy(cline_cluster2D)
            call cline_cluster2D%delete('update_frac')
            call cline_cluster2D%delete('fillin')
            call cline_cluster2D%delete('endit')
            call execute_cluster2D('terminal_greedy', nstages + 1)
            call terminal_refs%kill
        end subroutine execute_terminal_pass

        subroutine output_stats( prefix )
            character(len=*) :: prefix
            real, allocatable :: M(:,:)
            allocate(M(nptcls_eff,2))
            M(:,1) = spproj_field%get_all('class', nonzero=.true.)
            M(:,2) = spproj_field%get_all('corr',  nonzero=.true.)
            call rmat2file(M, string(trim(prefix)//'_class_scores.mat'))
            deallocate(M)
        end subroutine output_stats

        subroutine gen_final_cavgs( iter )
            type(commander_make_cavgs_distr) :: xmake_cavgs_distr
            type(commander_make_cavgs)       :: xmake_cavgs
            type(commander_rank_cavgs)       :: xrank_cavgs
            type(cmdline)                    :: cline_make_cavgs, cline_rank_cavgs
            type(string)                     :: finalcavgs, finalcavgs_ranked
            integer :: iter
            finalcavgs = CAVGS_ITER_FBODY//int2str_pad(iter,3)//params%ext%to_char()
            ! classes generation
            cline_make_cavgs = cline ! ncls is transferred here
            call cline_make_cavgs%delete('ptcl_src')
            call cline_make_cavgs%delete('autoscale')
            call cline_make_cavgs%delete('balance')
            call cline_make_cavgs%delete('smpd_crop')
            call cline_make_cavgs%delete('box_crop')
            call cline_make_cavgs%set('prg',        'make_cavgs')
            call cline_make_cavgs%set('refs',       finalcavgs)
            call cline_make_cavgs%set('which_iter', iter)
            ! Cavgs final output is regularized
            call cline_make_cavgs%set('ml_reg', 'yes')
            if( l_shmem )then
                call xmake_cavgs%execute(cline_make_cavgs)
            else
                call xmake_cavgs_distr%execute(cline_make_cavgs)
            endif
            ! adding cavgs & FRCs to project
            call spproj%read_segment('out', params%projfile)
            call spproj%add_frcs2os_out( string(FRCS_FILE), 'frc2D')
            call spproj%add_cavgs2os_out(finalcavgs, params%smpd, imgkind='cavg', mskdiam=params%mskdiam)
            call spproj%write_segment_inside('out', params%projfile)
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = CAVGS_ITER_FBODY//int2str_pad(iter,3)//'_ranked'//params%ext%to_char()
            call cline_rank_cavgs%set('projfile', params%projfile)
            call cline_rank_cavgs%set('stk',      finalcavgs)
            call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
            call xrank_cavgs%execute( cline_rank_cavgs )
            call cline_make_cavgs%kill
            call cline_rank_cavgs%kill
            call finalcavgs%kill
            call finalcavgs_ranked%kill
        end subroutine gen_final_cavgs

        subroutine write_abinitio_benchmark( iter, phase, stage )
            integer,          intent(in) :: iter, stage
            character(len=*), intent(in) :: phase
            type(string) :: benchfname, cluster_refine
            integer :: fnr
            if( .not. L_BENCH_GLOB ) return
            rt_tot = toc(t_tot)
            benchfname = string('ABINITIO2D_BENCH_ITER')//int2str_pad(iter,3)//'.txt'
            cluster_refine = cline_cluster2D%get_carg('refine')
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'abinitio2D phase                   : ', trim(phase)
            write(fnr,'(a,a)')  'abinitio2D execution mode          : ', merge('shared-memory', 'distributed  ', l_shmem)
            write(fnr,'(a,a)')  'abinitio2D requested refine mode   : ', trim(params%refine)
            write(fnr,'(a,a)')  'abinitio2D cluster2D refine mode   : ', cluster_refine%to_char()
            write(fnr,'(a,i0)') 'abinitio2D stage                   : ', stage
            write(fnr,'(a,i0)') 'abinitio2D nstages                 : ', nstages
            write(fnr,'(a,i0)') 'abinitio2D nclasses                : ', params%ncls
            write(fnr,'(a,i0)') 'abinitio2D sampled target          : ', nsample_target_2D
            write(fnr,'(a,i0)') 'abinitio2D effective particles     : ', nptcls_eff
            if( stage >= 1 .and. stage <= size(stage_parms) )then
                write(fnr,'(a,f0.3)') 'abinitio2D stage update fraction   : ', stage_parms(stage)%update_frac
                write(fnr,'(a,f0.2)') 'abinitio2D stage low-pass          : ', stage_parms(stage)%lp
            endif
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'abinitio2D setup/preparation        :', rt_setup
            write(fnr,'(a,1x,f0.2)') 'abinitio2D calc_pspec               :', rt_calc_pspec
            write(fnr,'(a,1x,f0.2)') 'abinitio2D cluster2D stage          :', rt_cluster2D
            write(fnr,'(a,1x,f0.2)') 'abinitio2D final class assembly     :', rt_final_cavgs
            write(fnr,'(a,1x,f0.2)') 'abinitio2D total time               :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'abinitio2D % accounted for          :', &
                &((rt_setup + rt_calc_pspec + rt_cluster2D + rt_final_cavgs) / rt_tot) * 100.
            call fclose(fnr)
            call cluster_refine%kill
            call benchfname%kill
        end subroutine write_abinitio_benchmark

    end subroutine exec_abinitio2D_workflow

end module simple_commanders_abinitio2D
