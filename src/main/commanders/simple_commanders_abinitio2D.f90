! concrete commander: high-level workflows
module simple_commanders_abinitio2D
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_commander_base,      only: commander_base
use simple_exec_helpers,        only: set_shmem_flag
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_commanders_cavgs
use simple_commanders_cluster2D
use simple_commanders_euclid
use simple_qsys_funs
implicit none

public :: commander_abinitio2D, autosample2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio2D
    contains
    procedure :: execute => exec_abinitio2D
end type commander_abinitio2D

! Dimensions
real,    parameter :: SMPD_TARGET     = 2.67
integer, parameter :: MINBOXSZ        = 88
! Stages
integer, parameter :: NSTAGES_CLS     = 6
integer, parameter :: NSTAGES_SINGLE  = 1
integer, parameter :: ITS_INCR_SINGLE = 4
integer, parameter :: ITS_INCR        = 5
integer, parameter :: PHASES(2)       = [4, 6]
! Stages variables
real,    parameter :: ICM_LAMBDA      = 1.0
integer, parameter :: EXTR_LIM_LOCAL  = 20

! convenience type
type stage_params
    real    :: lp=0., smpd_crop=0., trslim=0.
    integer :: box_crop = 0, max_cls_pop=0, nptcls=0
    logical :: l_lpset=.false.
end type stage_params

! class variables
type(stage_params), allocatable :: stage_parms(:)

contains

    subroutine exec_abinitio2D( self, cline )
        use simple_classaverager, only: cavger_kill, cavger_write, cavger_readwrite_partial_sums
        class(commander_abinitio2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_cluster2D_distr)  :: xcluster2D_distr
        type(commander_cluster2D)        :: xcluster2D
        type(commander_calc_pspec_distr) :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                    :: cline_cluster2D, cline_calc_pspec
        ! other
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        class(oris),             pointer :: spproj_field
        integer :: maxits, istage, last_iter, nptcls_eff, nstages
        logical :: l_shmem, l_inpl
        call cline%set('oritype',   'ptcl2D')
        call cline%set('sigma_est', 'global')
        if( .not. cline%defined('autoscale')  ) call cline%set('autoscale',  'yes')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('center')     ) call cline%set('center',     'yes')
        if( .not. cline%defined('center_type')) call cline%set('center_type','seg')
        if( .not. cline%defined('sh_first')   ) call cline%set('sh_first',   'no')
        if( .not. cline%defined('cls_init')   ) call cline%set('cls_init',   'rand')
        if( .not. cline%defined('icm')        ) call cline%set('icm',        'yes')
        if( .not. cline%defined('gauref')     ) call cline%set('gauref',     'yes')
        if( .not. cline%defined('polar')      ) call cline%set('polar',      'no')
        if( .not. cline%defined('lambda')     ) call cline%set('lambda',     ICM_LAMBDA)
        if( .not. cline%defined('extr_lim')   ) call cline%set('extr_lim',   EXTR_LIM_LOCAL)
        if( .not. cline%defined('rank_cavgs') ) call cline%set('rank_cavgs', 'yes')
        if( .not. cline%defined('stats')      ) call cline%set('stats',      'no')
        if( .not. cline%defined('refine')     ) call cline%set('refine',     'snhc_smpl')
        if( .not. cline%defined('ref_type')   ) call cline%set('ref_type',   'polar_cavg')
        ! shared memory execution
        l_shmem = set_shmem_flag(cline)
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        maxits = params%extr_lim
        call cline%delete('stats')
        ! check refinement flag and set stages
        l_inpl = .false.
        select case(trim(params%refine))
            case('inpl','inpl_smpl')
                ! in-plane refinement, one stage only
                nstages = NSTAGES_SINGLE
                l_inpl  = .true.
            case('snhc','snhc_smpl')
                ! usual suspects
                nstages = NSTAGES_CLS
            case('prob','prob_smpl','prob_smpl_shc','prob_greedy')
                ! prob family of algorithms
                nstages = NSTAGES_CLS
            case DEFAULT
                THROW_HARD('Unsupported REFINE argument: '//trim(params%refine))
        end select
        ! override # stages
        if( cline%defined('nstages') ) nstages = min(params%nstages,NSTAGES_CLS)
        allocate(stage_parms(nstages))
        ! read project
        call spproj%read(params%projfile)
        call set_dims                   ! set downscaling
        call inirefs                    ! deal with initial references
        call set_lplims                 ! set resolutions limits
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
        if( .not.l_inpl ) call spproj_field%delete_2Dclustering
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
            call set_cline_cluster2D(istage)
            ! classify
            call execute_cluster2D
        enddo
        ! transfer 2D shifts to 3D field
        call spproj%read_segment(params%oritype,params%projfile)
        call spproj%read_segment('ptcl3D',params%projfile)
        call spproj%os_ptcl3D%transfer_2Dshifts(spproj_field)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        ! weights & final mapping of particles
        if( nstages == 1 )then
            ! no autosampling
        else
            if( trim(params%autosample).eq.'yes' )then
                call spproj_field%set_all2single('w',1.)
                call spproj%write_segment_inside(params%oritype, params%projfile)
                ! mapping of all particles
                if( stage_parms(nstages)%nptcls < nptcls_eff )then
                    call set_final_mapping
                    call execute_cluster2D
                    if( trim(params%stats).eq.'yes' )then
                        call spproj%read_segment(params%oritype,params%projfile)
                        call output_stats('final')
                    endif
                else
                    if( trim(params%stats).eq.'yes' ) call output_stats('final')
                endif
            else
                if( trim(params%stats).eq.'yes' ) call output_stats('final')
            endif
        endif
        ! final class generation & ranking
        last_iter = cline_cluster2D%get_iarg('endit')
        if ( trim(params%rank_cavgs).eq.'yes' )then
            call gen_final_cavgs(last_iter)
        else
            if( l_shmem )then
                ! making sure e/o are written
                params%refs      = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//params%ext%to_char()
                params%refs_even = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//'_even'//params%ext%to_char()
                params%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//'_odd'//params%ext%to_char()
                call cavger_write(params%refs,     'merged')
                call cavger_write(params%refs_even,'even')
                call cavger_write(params%refs_odd, 'odd')
                ! required for remapping
                if( trim(params%chunk).eq.'yes') call cavger_readwrite_partial_sums('write')
                call cavger_kill
            endif
        endif
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
        call qsys_cleanup
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
            use simple_procimgstk,         only: copy_imgfile
            use simple_commanders_imgproc, only: commander_scale
            use simple_commanders_volops,  only: commander_noisevol
            type(commander_scale)    :: xscale
            type(commander_noisevol) :: xnoisevol
            type(cmdline)            :: cline_scalerefs, cline_noisevol
            type(string) :: refs, refs_even, refs_odd
            real         :: smpd
            integer      :: ldim(3), ncls
            logical      :: eo
            if( .not.cline%defined('refs') .and. (.not.l_inpl) )then
                ! Starting from scratch: appropriate starting references
                ! will be generated by cluster2D/cluster2D_distr
                if( trim(params%ref_type).eq.'comlin_hybrid' )then
                    call cline_noisevol%set('prg',    'noisevol')
                    call cline_noisevol%set('smpd',   stage_parms(1)%smpd_crop)
                    call cline_noisevol%set('box',    stage_parms(1)%box_crop)
                    call cline_noisevol%set('nspace', params%ncls)
                    call xnoisevol%execute_safe(cline_noisevol)
                    params%refs      = 'start2Drefs.mrc'
                    params%refs_even = 'start2Drefs_even.mrc'
                    params%refs_odd  = 'start2Drefs_odd.mrc'
                    call cline%set('refs', params%refs)
                    call cline_noisevol%kill
                endif
                return
            endif
            if( cline%defined('refs') )then
                refs = params%refs
                call find_ldim_nptcls(refs, ldim, ncls, smpd=smpd)
                ldim(3) = 1
            else
                ! l_inpl=true & not.cline%defined('refs'):
                ! refinement so getting references from previous run
                call spproj%get_cavgs_stk(refs, ncls, smpd, 'cavg')
            endif
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
                call xscale%execute_safe(cline_scalerefs)
                if( eo )then
                    call cline_scalerefs%set('stk',    refs_even)
                    call cline_scalerefs%set('outstk', params%refs_even)
                    call xscale%execute_safe(cline_scalerefs)
                    call cline_scalerefs%set('stk',    refs_odd)
                    call cline_scalerefs%set('outstk', params%refs_odd)
                    call xscale%execute_safe(cline_scalerefs)
                endif
                call cline_scalerefs%kill
            endif
        end subroutine inirefs

        ! Set resolution limits
        subroutine set_lplims
            use simple_class_frcs, only: class_frcs
            type(class_frcs)              :: clsfrcs
            type(string) :: frcs
            real         :: lpstart, lpstop, cenlp
            integer      :: istage
            ! Resolution limits
            call mskdiam2lplimits_here(params%mskdiam, lpstart, lpstop, cenlp)
            lpstart = max(lpstart, 2.*params%smpd_crop)
            lpstop  = max(lpstop,  2.*params%smpd_crop)
            cenlp   = max(cenlp,   2.*params%smpd_crop)
            ! Stages resolution limits
            if( nstages == 1 )then
                if( cline%defined('lp') )then
                    stage_parms(1)%l_lpset = .true.
                    lpstart                = params%lp
                    stage_parms(1)%lp      = params%lp
                    lpstop                 = params%lp
                else
                    call spproj%get_frcs(frcs, 'frc2D')
                    call clsfrcs%read(frcs)
                    call clsfrcs%crop(stage_parms(1)%smpd_crop, stage_parms(1)%box_crop)
                    call clsfrcs%write(string(FRCS_FILE))
                    stage_parms(1)%l_lpset = .false.
                    lpstart                = clsfrcs%estimate_lp_for_align()
                    stage_parms(1)%lp      = lpstart    ! will not be used
                    lpstop                 = 2.*params%smpd_crop
                    call clsfrcs%kill
                endif
                if( .not. cline%defined('lpstart') ) params%lpstart = lpstart
                if( .not. cline%defined('lpstop')  ) params%lpstop  = lpstop
            else
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
                    stage_parms(1)%lp      = params%lpstart
                    stage_parms(1)%l_lpset = .true.
                    do istage = 2, NSTAGES-1
                        stage_parms(istage)%lp      = stage_parms(istage-1)%lp - (stage_parms(istage-1)%lp - params%lpstop)/2.0
                        stage_parms(istage)%l_lpset = .true.
                    end do
                    stage_parms(NSTAGES-1)%lp    = params%lpstop
                    stage_parms(NSTAGES)%l_lpset = .false.
                    stage_parms(NSTAGES)%lp      = params%lpstop
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
            if( trim(params%autosample).eq.'yes' )then
                do i = 1,NSTAGES
                    call autosample2D(cline, nptcls_eff, params%ncls,&
                    &stage_parms(i)%max_cls_pop, stage_parms(i)%nptcls)
                    if( (params%stage > 0) .and. (params%stage < NSTAGES))then
                        if( i<=params%stage )then
                            call autosample2D(cline, nptcls_eff, params%ncls,&
                            &stage_parms(i)%max_cls_pop, stage_parms(i)%nptcls, popfac=0.5)
                        endif
                    endif
                enddo
            endif
        end subroutine set_sampling

        subroutine prep_command_lines( cline )
            use simple_default_clines, only: set_automask2D_defaults
            class(cmdline), intent(in) :: cline
            cline_cluster2D  = cline
            cline_calc_pspec = cline
            ! initial sigma2
            call cline_calc_pspec%set('prg',      'calc_pspec')
            ! cluster2D
            call cline_cluster2D%set('prg',       'cluster2D')
            call cline_cluster2D%set('wiener',    'full')
            call cline_cluster2D%set('ptclw',     'no')
            call cline_cluster2D%set('ml_reg',    'no')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call cline_cluster2D%set('chunk',     'no')
            call set_automask2D_defaults( cline_cluster2D )
        end subroutine prep_command_lines

        subroutine set_cline_cluster2D( istage )
            integer, intent(in) :: istage
            type(string) :: sh_first, refine, center, objfun, refs, icm, gauref
            integer      :: iphase, iter, imaxits, cc_iters, minits, extr_iter
            real         :: trs, lambda, gaufreq
            logical      :: l_gauref, l_gaufreq_input
            refine = trim(params%refine)
            ! filter for Fourier polar representation
            l_gauref        = .false.
            l_gaufreq_input = .false.
            if( trim(params%polar).eq.'yes' )then
                l_gauref        = trim(params%gauref).eq.'yes'
                l_gaufreq_input = cline%defined('gaufreq')
                ! ICM used for cartesian filtering of random refs
                params%l_icm    = istage==1
            endif
            ! objective function
            if( params%cc_objfun == OBJFUN_CC )then
                objfun   = 'cc'
                cc_iters = 999
            else
                objfun   = 'euclid'
                cc_iters = 0
            endif
            ! iteration number book-keeping
            iter = 0
            if( cline_cluster2D%defined('endit') ) iter = cline_cluster2D%get_iarg('endit')
            iter = iter + 1
            call cline_cluster2D%delete('which_iter')
            ! Phase & stages
            if( nstages == 1 )then
                ! Single stage refinement
                minits    = ITS_INCR_SINGLE
                imaxits   = ITS_INCR_SINGLE+2
                extr_iter = 999
                trs       = stage_parms(1)%trslim
                sh_first  = trim(params%sh_first)
                center    = trim(params%center)
                cc_iters  = 0
                objfun    = 'euclid'
                ! Filters deactivated
                icm       = 'no'
                gauref    = 'no'
                ! continue from previous references
                refs      = params%refs
                ! resolution limit
                call cline_cluster2D%set('lpstart', params%lpstart)
                call cline_cluster2D%set('lpstop',  params%lpstop)
                if( stage_parms(1)%l_lpset ) call cline_cluster2D%set('lp', params%lp)
            else
                ! phase logics
                if(      istage <= PHASES(1) )then
                    iphase = 1
                else if( istage <= PHASES(2) )then
                    iphase = 2
                else
                    iphase = 0
                    THROW_HARD('Invalid istage index')
                endif
                ! stages
                select case(iphase)
                case(1)
                    ! phase constants
                    extr_iter = 0
                    ! phase variables
                    imaxits   = nint(real(istage)*real(maxits)/real(PHASES(1)))
                    minits    = imaxits
                    select case(istage)
                    case(1)
                        trs      = 0.
                        sh_first = 'no'
                        center   = 'no'
                        if( cline%defined('refs') )then
                            refs     = params%refs
                        else
                            refs     = NIL
                        endif
                        if( params%l_icm )then
                            icm      = 'yes'
                            lambda   = params%lambda
                        else
                            icm      = 'no'
                        endif
                        if( l_gauref )then
                            gauref   = 'yes'
                            if( l_gaufreq_input )then
                                gaufreq  = params%gaufreq
                            else
                                gaufreq  = stage_parms(istage)%lp
                            endif
                        else
                            gauref   = 'no'
                        endif
                    case(2)
                        trs          = stage_parms(istage)%trslim
                        sh_first     = trim(params%sh_first)
                        center       = trim(params%center)
                        refs         = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                        if( params%l_icm )then
                            icm      = 'yes'
                            lambda   = params%lambda/2.
                        else
                            icm      = 'no'
                        endif
                        if( l_gauref )then
                            gauref   = 'yes'
                            if( l_gaufreq_input )then
                                gaufreq  = params%gaufreq * stage_parms(2)%lp / stage_parms(1)%lp
                            else
                                gaufreq  = stage_parms(istage)%lp
                            endif
                        else
                            gauref   = 'no'
                        endif
                    case(3)
                        trs          = stage_parms(istage)%trslim
                        sh_first     = trim(params%sh_first)
                        center       = trim(params%center)
                        refs         = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                        if( params%l_icm )then
                            icm      = 'yes'
                            lambda   = params%lambda/4.
                        else
                            icm      = 'no'
                        endif
                        if( l_gauref )then
                            gauref   = 'yes'
                            if( l_gaufreq_input )then
                                gaufreq  = params%gaufreq * stage_parms(3)%lp / stage_parms(2)%lp
                            else
                                gaufreq  = stage_parms(istage)%lp
                            endif
                        else
                            gauref   = 'no'
                        endif
                    case(4)
                        trs          = stage_parms(istage)%trslim
                        sh_first     = trim(params%sh_first)
                        center       = trim(params%center)
                        refs         = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                        icm          = 'no'
                        gauref       = 'no'
                        if( l_gauref )then
                            gauref   = 'yes'
                            gaufreq  = 2.0 * params%smpd
                        else
                            gauref   = 'no'
                        endif
                    end select
                case(2)
                    ! phase constants
                    imaxits   = iter+ITS_INCR-1
                    sh_first  = trim(params%sh_first)
                    trs       = stage_parms(istage)%trslim
                    center    = trim(params%center)
                    extr_iter = params%extr_lim+1
                    refs      = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                    icm       = 'no'
                    gauref    = 'no'
                    minits    = iter+1
                end select
                if( stage_parms(istage)%max_cls_pop > 0 )then
                    call cline_cluster2D%set('maxpop', stage_parms(istage)%max_cls_pop)
                    if( cline%defined('update_frac') )then
                        call cline_cluster2D%set('update_frac', params%update_frac)
                    endif
                endif
            endif
            ! command line update
            if( stage_parms(istage)%l_lpset )then
                call cline_cluster2D%set('lp',    stage_parms(istage)%lp)
            else
                call cline_cluster2D%delete('lp')
            endif
            if( refs .ne. NIL ) call cline_cluster2D%set('refs', refs)
            if( extr_iter > 0 )then
                call cline_cluster2D%set('extr_iter', extr_iter)
            else
                call cline_cluster2D%delete('extr_iter')
            endif
            call cline_cluster2D%set('icm', icm)
            if( icm.eq.'yes' )then
                call cline_cluster2D%set('lambda', lambda)
            else
                call cline_cluster2D%delete('lambda')
            endif
            call cline_cluster2D%set('gauref', gauref)
            if( gauref.eq.'yes' )then
                call cline_cluster2D%set('gaufreq', gaufreq)
            else
                call cline_cluster2D%delete('gaufreq')
            endif
            if( trim(params%polar).eq.'yes')then
                if( trim(params%ref_type)=='comlin_hybrid') center = 'no' ! because the references form a volume
            endif
            call cline_cluster2D%set('minits',    minits)
            call cline_cluster2D%set('maxits',    imaxits)
            call cline_cluster2D%set('startit',   iter)
            call cline_cluster2D%set('refine',    refine)
            call cline_cluster2D%set('objfun',    objfun)
            call cline_cluster2D%set('cc_iters',  cc_iters)
            call cline_cluster2D%set('trs',       trs)
            call cline_cluster2D%set('sh_first',  sh_first)
            call cline_cluster2D%set('center',    center)
            call cline_cluster2D%set('box_crop',  stage_parms(istage)%box_crop)
            call cline_cluster2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
            call cline_cluster2D%delete('maxpop')
            call cline_cluster2D%delete('nsample_max')
            call cline_cluster2D%delete('nsample')
            call cline_cluster2D%delete('autosample')
            call cline_cluster2D%delete('update_frac')
            call cline_cluster2D%delete('endit')
        end subroutine set_cline_cluster2D

        subroutine execute_cluster2D
            call del_file(CLUSTER2D_FINISHED)
            ! Initial sigma2
            if( istage == 1 )then
                call xcalc_pspec_distr%execute_safe(cline_calc_pspec)
            endif
            ! clustering
            if( l_shmem )then
                call xcluster2D%execute_safe(cline_cluster2D)
            else
                call xcluster2D_distr%execute_safe(cline_cluster2D)
            endif
        end subroutine execute_cluster2D

        subroutine output_stats( prefix )
            character(len=*) :: prefix
            real, allocatable :: M(:,:)
            allocate(M(nptcls_eff,2))
            M(:,1) = spproj_field%get_all('class', nonzero=.true.)
            M(:,2) = spproj_field%get_all('corr',  nonzero=.true.)
            call rmat2file(M, string(trim(prefix)//'_class_scores.mat'))
        end subroutine output_stats

        subroutine set_final_mapping
            type(string) :: refs
            integer :: iter, minits
            iter      = cline_cluster2D%get_iarg('endit') + 1
            refs      = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
            minits    = iter
            if( params%l_autoscale ) call cline_cluster2D%set('restore_cavgs', 'no')
            call cline_cluster2D%set('startit',   iter)
            call cline_cluster2D%set('minits',    minits)
            call cline_cluster2D%set('maxits',    minits)
            call cline_cluster2D%set('refine',    'greedy_smpl')
            call cline_cluster2D%set('objfun',    'euclid')
            call cline_cluster2D%set('cc_iters',  0)
            call cline_cluster2D%set('trs',       stage_parms(NSTAGES)%trslim)
            call cline_cluster2D%set('sh_first',  params%sh_first)
            call cline_cluster2D%set('center',    'no')
            call cline_cluster2D%set('box_crop',  stage_parms(NSTAGES)%box_crop)
            call cline_cluster2D%set('smpd_crop', stage_parms(NSTAGES)%smpd_crop)
            call cline_cluster2D%set('refs',      refs)
            call cline_cluster2D%set('extr_iter', params%extr_lim+1)
            call cline_cluster2D%set('icm',       'no')
            call cline_cluster2D%delete('lp')
            call cline_cluster2D%delete('lambda')
            call cline_cluster2D%delete('nsample_max')
            call cline_cluster2D%delete('nsample')
            call cline_cluster2D%delete('autosample')
            call cline_cluster2D%delete('update_frac')
            call cline_cluster2D%delete('maxpop')
            call cline_cluster2D%delete('endit')
            write(logfhandle,'(A,/,A)')'>>>','>>> FINAL MAPPING'
        end subroutine set_final_mapping

        subroutine gen_final_cavgs( iter )
            type(commander_make_cavgs_distr) :: xmake_cavgs_distr
            type(commander_make_cavgs)       :: xmake_cavgs
            type(commander_rank_cavgs)       :: xrank_cavgs
            type(cmdline)                    :: cline_make_cavgs, cline_rank_cavgs
            type(string)                     :: finalcavgs, finalcavgs_ranked
            integer :: iter
            finalcavgs = CAVGS_ITER_FBODY//int2str_pad(iter,3)//params%ext%to_char()
            ! classes generation
            if( params%l_autoscale )then
                cline_make_cavgs = cline ! ncls is transferred here
                call cline_make_cavgs%delete('polar')
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%delete('smpd_crop')
                call cline_make_cavgs%delete('box_crop')
                call cline_make_cavgs%set('prg',        'make_cavgs')
                call cline_make_cavgs%set('refs',       finalcavgs)
                call cline_make_cavgs%set('which_iter', iter)
                if( l_shmem )then
                    call xmake_cavgs%execute_safe(cline_make_cavgs)
                else
                    call xmake_cavgs_distr%execute_safe(cline_make_cavgs)
                endif
            endif
            ! adding cavgs & FRCs to project
            call spproj%read_segment('out', params%projfile)
            call spproj%add_frcs2os_out( string(FRCS_FILE), 'frc2D')
            call spproj%add_cavgs2os_out(finalcavgs, params%smpd, imgkind='cavg')
            call spproj%write_segment_inside('out', params%projfile)
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = CAVGS_ITER_FBODY//int2str_pad(iter,3)//'_ranked'//params%ext%to_char()
            call cline_rank_cavgs%set('projfile', params%projfile)
            call cline_rank_cavgs%set('stk',      finalcavgs)
            call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
            call xrank_cavgs%execute( cline_rank_cavgs )
        end subroutine gen_final_cavgs

        ! attempt at defining resolution limits
        subroutine mskdiam2lplimits_here( mskdiam, lpstart,lpstop, lpcen )
            real, intent(in)    :: mskdiam
            real, intent(inout) :: lpstart,lpstop, lpcen
            lpstart = min(max(mskdiam/12., 15.), 20.)
            lpstop  = min(max(mskdiam/22.,  6.),  8.)
            lpcen   = min(max(mskdiam/6.,  20.), 30.)
        end subroutine mskdiam2lplimits_here

    end subroutine exec_abinitio2D

    subroutine autosample2D( cline, nptcls, ncls, max_pop, max_nptcls, popfac )
        use simple_parameters, only: params_glob
        class(cmdline), intent(in)  :: cline
        integer,        intent(in)  :: nptcls, ncls
        integer,        intent(out) :: max_pop, max_nptcls
        real, optional, intent(in)  :: popfac ! testing purpose
        max_pop    = 0
        max_nptcls = nptcls
        if( trim(params_glob%autosample).ne.'yes' ) return
        ! Class population limit
        max_pop = MAXPOP_CLS
        if( cline%defined('nsample') ) max_pop = params_glob%nsample
        ! Total population limit
        max_nptcls = min(nptcls, 2*max_pop*ncls)
        if( cline%defined('nsample_max') )then
            max_nptcls = min(params_glob%nsample_max, max_nptcls)
        else
            if( present(popfac) )then
                max_nptcls = min(nint(popfac*real(MAXPOP_PTCLS)), nptcls)
            else
                max_nptcls = min(MAXPOP_PTCLS, nptcls)
            endif
        endif
    end subroutine autosample2D

end module simple_commanders_abinitio2D
