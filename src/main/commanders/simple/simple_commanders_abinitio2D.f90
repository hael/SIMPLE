!@descr: ab initio 2D analysis
module simple_commanders_abinitio2D
use simple_commanders_api
use simple_commanders_cavgs
use simple_commanders_cluster2D
use simple_abinitio2D_controller 
implicit none

public :: commander_abinitio2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio2D
    contains
    procedure :: execute => exec_abinitio2D
end type commander_abinitio2D

! class variables
type(stage_params), allocatable :: stage_parms(:)

contains

    subroutine exec_abinitio2D( self, cline )
        use simple_classaverager
        class(commander_abinitio2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_cluster2D)  :: xcluster2D
        type(commander_calc_pspec) :: xcalc_pspec
        ! command lines
        type(cmdline)              :: cline_cluster2D, cline_calc_pspec
        ! other
        type(parameters)           :: params
        type(sp_project)           :: spproj
        class(oris),       pointer :: spproj_field
        integer :: maxits, istage, last_iter, nptcls_eff, nstages
        logical :: l_shmem
        call cline%set('oritype',   'ptcl2D')
        call cline%set('sigma_est', 'global')
        if( .not. cline%defined('autoscale')     ) call cline%set('autoscale',     'yes')
        if( .not. cline%defined('mkdir')         ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('center')        ) call cline%set('center',        'yes')
        if( .not. cline%defined('center_type')   ) call cline%set('center_type',   'seg')
        if( .not. cline%defined('cls_init')      ) call cline%set('cls_init',      'rand')
        if( .not. cline%defined('gauref')        ) call cline%set('gauref',        'yes')
        if( .not. cline%defined('polar')         ) call cline%set('polar',         'no')
        if( .not. cline%defined('extr_lim')      ) call cline%set('extr_lim',      EXTR_LIM_LOCAL)
        if( .not. cline%defined('nits_per_stage')) call cline%set('nits_per_stage',ITS_INCR)
        if( .not. cline%defined('eo_stage')      ) call cline%set('eo_stage',      EO_STAGE)
        if( .not. cline%defined('rank_cavgs')    ) call cline%set('rank_cavgs',    'yes')
        if( .not. cline%defined('stats')         ) call cline%set('stats',         'no')
        if( .not. cline%defined('refine')        ) call cline%set('refine',        'snhc_smpl')
        if( .not. cline%defined('ref_type')      ) call cline%set('ref_type',      'polar_cavg')
        if( .not. cline%defined('ml_reg')        ) call cline%set('ml_reg',        'yes')
        ! shared memory execution
        l_shmem = set_shmem_flag(cline)
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        maxits = params%extr_lim
        call cline%delete('stats')
        ! check refinement flag and set stages
        call determine_abinitio2D_stages(params, nstages)
        ! override # stages
        if( cline%defined('nstages') ) nstages = min(params%nstages,NSTAGES_CLS)
        allocate(stage_parms(nstages))
        ! read project
        call spproj%read(params%projfile)
        call set_dims                   ! set downscaling
        call inirefs                    ! deal with initial references
        call set_lplims(nstages)        ! set resolutions limits
        call prep_command_lines(cline)  ! prepare class command lines
        call set_sampling               ! sampling
        ! summary
        do istage = 1,nstages
            write(logfhandle,'(A,I2,A,L1,F6.1,2I8)')'>>> STAGE ', istage,' LPSET LP MAXCLSPOP NPTCLS: ',&
            &stage_parms(istage)%l_lpset,stage_parms(istage)%lp, stage_parms(istage)%max_cls_pop,&
            &stage_parms(istage)%nptcls
        end do
        ! prep particles field
        call spproj_field%set_all2single('w',1.)
        call spproj_field%delete_2Dclustering
        if( spproj_field%get_nevenodd() == 0 ) call spproj_field%partition_eo
        call spproj%write_segment_inside(params%oritype, params%projfile)
        call spproj%split_stk(params%nparts, dir=string(PATH_PARENT))
        ! Frequency marching
        do istage = 1,nstages
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
        enddo
        ! transfer 2D shifts to 3D field
        call spproj%read_segment(params%oritype,params%projfile)
        call spproj%read_segment('ptcl3D',params%projfile)
        call spproj%os_ptcl3D%transfer_2Dshifts(spproj_field)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        ! weights & final mapping of particles
        if( trim(params%stats).eq.'yes' ) call output_stats('final')
        ! final class generation & ranking
        last_iter = cline_cluster2D%get_iarg('endit')
        call gen_final_cavgs(last_iter)
        ! cleanup
        call del_file('start2Drefs'//params%ext%to_char())
        call del_file('start2Drefs_even'//params%ext%to_char())
        call del_file('start2Drefs_odd'//params%ext%to_char())
        call del_files(DIST_FBODY,      params%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params%nparts,ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        deallocate(stage_parms)
        call spproj%kill
        nullify(spproj_field)
        call qsys_cleanup(params)
        call simple_touch(ABINITIO2D_FINISHED)
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
            call find_ldim_nptcls(refs, ldim, ncls, smpd=smpd)
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
        end subroutine inirefs

        ! Set resolution limits
        subroutine set_lplims( local_nstages )
            use simple_class_frcs, only: class_frcs
            integer, intent(in) :: local_nstages
            type(class_frcs) :: clsfrcs
            type(string)     :: frcs
            real    :: lpstart, lpstop, cenlp
            integer :: istage
            ! Resolution limits
            call mskdiam2lplimits_cluster2D(params%mskdiam, lpstart, lpstop, cenlp)
            lpstart = max(lpstart, 2.*params%smpd_crop)
            lpstop  = max(lpstop,  2.*params%smpd_crop)
            cenlp   = max(cenlp,   2.*params%smpd_crop)
            ! Stages resolution limits
            if( cline%defined('lp') )then
                ! Set lp throughout
                stage_parms(:)%lp      = params%lp
                stage_parms(:)%l_lpset = .true.
                params%lpstart = params%lp
                params%lpstop  = params%lp
            else
                ! Frequency marching
                if( .not. cline%defined('lpstart') ) params%lpstart = lpstart
                if( .not. cline%defined('lpstop')  ) params%lpstop  = lpstop
                if( trim(params%eo_stage)=='yes')then
                    stage_parms(1)%lp      = params%lpstart
                    stage_parms(1)%l_lpset = .true.
                    do istage = 2, local_nstages-1
                        stage_parms(istage)%lp      = stage_parms(istage-1)%lp - (stage_parms(istage-1)%lp - params%lpstop)/2.0
                        stage_parms(istage)%l_lpset = .true.
                    end do
                    stage_parms(local_nstages-1)%lp    = params%lpstop
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
            integer :: i
            nptcls_eff = spproj%count_state_gt_zero()
            stage_parms(:)%max_cls_pop = 0
            stage_parms(:)%nptcls      = nptcls_eff
        end subroutine set_sampling

        subroutine prep_command_lines( cline )
            class(cmdline), intent(in) :: cline
            cline_cluster2D  = cline
            cline_calc_pspec = cline
            ! initial sigma2
            call cline_calc_pspec%set('prg',      'calc_pspec')
            ! cluster2D
            call cline_cluster2D%set('prg',       'cluster2D')
            call cline_cluster2D%set('wiener',    'full')
            call cline_cluster2D%set('ptclw',     'no')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call cline_cluster2D%set('chunk',     'no')
            call set_automask2D_defaults( cline_cluster2D )
        end subroutine prep_command_lines

        subroutine execute_cluster2D
            call del_file(CLUSTER2D_FINISHED)
            ! Initial sigma2
            if( istage == 1 )then
                call xcalc_pspec%execute(cline_calc_pspec)
            endif
            ! clustering
            call xcluster2D%execute(cline_cluster2D)
        end subroutine execute_cluster2D

        subroutine output_stats( prefix )
            character(len=*) :: prefix
            real, allocatable :: M(:,:)
            allocate(M(nptcls_eff,2))
            M(:,1) = spproj_field%get_all('class', nonzero=.true.)
            M(:,2) = spproj_field%get_all('corr',  nonzero=.true.)
            call rmat2file(M, string(trim(prefix)//'_class_scores.mat'))
        end subroutine output_stats

        subroutine gen_final_cavgs( iter )
            type(commander_make_cavgs_distr) :: xmake_cavgs_distr
            type(commander_make_cavgs)       :: xmake_cavgs
            type(commander_rank_cavgs)       :: xrank_cavgs
            type(commander_tree_rank_cavgs)  :: xtree_rank_cavgs
            type(cmdline)                    :: cline_make_cavgs, cline_rank_cavgs, cline_tree_ranked
            type(string)                     :: finalcavgs, finalcavgs_ranked
            integer :: iter
            finalcavgs = CAVGS_ITER_FBODY//int2str_pad(iter,3)//params%ext%to_char()
            ! classes generation
            cline_make_cavgs = cline ! ncls is transferred here
            call cline_make_cavgs%delete('polar')
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
            call spproj%add_cavgs2os_out(finalcavgs, params%smpd, imgkind='cavg')
            call spproj%add_sigma22os_out(sigma2_star_from_iter(iter))
            call spproj%write_segment_inside('out', params%projfile)
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = CAVGS_ITER_FBODY//int2str_pad(iter,3)//'_ranked'//params%ext%to_char()
            call cline_rank_cavgs%set('projfile', params%projfile)
            call cline_rank_cavgs%set('stk',      finalcavgs)
            call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
            call xrank_cavgs%execute( cline_rank_cavgs )
            ! write tree-based sets for tree-based refinements
            if( trim(params%refine) .eq. 'snhc_ptree' .or. trim(params%refine) .eq. 'greedy_tree' )then
                cline_tree_ranked = cline_rank_cavgs
                call cline_tree_ranked%set('oritype', 'cls2D')
                call cline_tree_ranked%set('refine',  trim(params%refine))
                call xtree_rank_cavgs%execute( cline_tree_ranked )
            endif
        end subroutine gen_final_cavgs

    end subroutine exec_abinitio2D

end module simple_commanders_abinitio2D
