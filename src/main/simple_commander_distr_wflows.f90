! concrete commander: distributed workflows
module simple_commander_distr_wflows
include 'simple_lib.f08'
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup, qsys_watcher
use simple_commander_base, only: commander_base
use simple_sp_project,     only: sp_project
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_builder,        only: builder
implicit none

public :: preprocess_distr_commander
public :: extract_distr_commander
public :: reextract_distr_commander
public :: motion_correct_distr_commander
public :: gen_pspecs_and_thumbs_distr_commander
public :: motion_correct_tomo_distr_commander
public :: ctf_estimate_distr_commander
public :: pick_distr_commander
public :: make_cavgs_distr_commander
public :: cluster2D_distr_commander
public :: refine3D_distr_commander
public :: reconstruct3D_distr_commander
public :: tseries_track_distr_commander
public :: scale_project_distr_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_distr_commander
  contains
    procedure :: execute      => exec_preprocess_distr
end type preprocess_distr_commander
type, extends(commander_base) :: extract_distr_commander
  contains
    procedure :: execute      => exec_extract_distr
end type extract_distr_commander
type, extends(commander_base) :: reextract_distr_commander
  contains
    procedure :: execute      => exec_reextract_distr
end type reextract_distr_commander
type, extends(commander_base) :: motion_correct_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_distr
end type motion_correct_distr_commander
type, extends(commander_base) :: gen_pspecs_and_thumbs_distr_commander
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs_distr
end type gen_pspecs_and_thumbs_distr_commander
type, extends(commander_base) :: motion_correct_tomo_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_tomo_distr
end type motion_correct_tomo_distr_commander
type, extends(commander_base) :: ctf_estimate_distr_commander
  contains
    procedure :: execute      => exec_ctf_estimate_distr
end type ctf_estimate_distr_commander
type, extends(commander_base) :: pick_distr_commander
  contains
    procedure :: execute      => exec_pick_distr
end type pick_distr_commander
type, extends(commander_base) :: make_cavgs_distr_commander
  contains
    procedure :: execute      => exec_make_cavgs_distr
end type make_cavgs_distr_commander
type, extends(commander_base) :: cluster2D_distr_commander
  contains
    procedure :: execute      => exec_cluster2D_distr
end type cluster2D_distr_commander
type, extends(commander_base) :: refine3D_distr_commander
  contains
    procedure :: execute      => exec_refine3D_distr
end type refine3D_distr_commander
type, extends(commander_base) :: reconstruct3D_distr_commander
  contains
    procedure :: execute      => exec_reconstruct3D_distr
end type reconstruct3D_distr_commander
type, extends(commander_base) :: tseries_track_distr_commander
  contains
    procedure :: execute      => exec_tseries_track_distr
end type tseries_track_distr_commander
type, extends(commander_base) :: scale_project_distr_commander
  contains
    procedure :: execute      => exec_scale_project_distr
end type scale_project_distr_commander

contains

    subroutine exec_preprocess_distr( self, cline )
        use simple_commander_preprocess, only: preprocess_commander
        class(preprocess_distr_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(cmdline)                 :: cline_make_pickrefs
        type(chash)                   :: job_descr
        type(sp_project)              :: spproj
        logical                       :: l_pick
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          20.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            6.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read in movies
        call spproj%read(params%projfile)
        ! DISTRIBUTED EXECUTION
        params%nptcls = spproj%get_nmovies()
        if( params%nptcls == 0 )then
            THROW_HARD('no movie to process! exec_preprocess_distr')
        endif
        if( params%nparts > params%nptcls ) THROW_HARD('# partitions (nparts) must be < number of entries in filetable')
        ! deal with numlen so that length matches JOB_FINISHED indicator files
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        l_pick = .false.
        if( cline%defined('refs') .or. cline%defined('vol1') )then
            l_pick = .true.
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', trim(PICKREFS)//params%ext)
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess_distr

    !> for extracting particle images from integrated DDD movies
    subroutine exec_extract_distr( self, cline )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        class(extract_distr_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj, spproj_part
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o_tmp
        type(oris)                              :: os_stk
        character(len=LONGSTRLEN),  allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        character(len=:),           allocatable :: mic_name, imgkind, boxfile_name
        integer :: boxcoords(2), lfoo(3)
        integer :: nframes,imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart
        if( .not. cline%defined('outside')   ) call cline%set('outside',   'no')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('stream')    ) call cline%set('stream',    'no')
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        call cline%set('nthr', 1.)
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read(params%projfile)
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot < params%nparts ) params%nparts = nmics_tot
        ! wipes previous stacks & particles
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        call spproj%os_cls2D%kill
        call spproj%os_cls3D%kill
        call spproj%os_out%kill
        call spproj%write
        ! input directory
        if( cline%defined('dir_box') )then
            if( params%mkdir.eq.'yes' .and. params%dir_box(1:1).ne.'/')then
                params%dir_box = trim(filepath(PATH_PARENT,params%dir_box))
            endif
            params%dir_box = simple_abspath(params%dir_box)
            if( file_exists(params%dir_box) )then
                call simple_list_files(trim(params%dir_box)//'/*.box', boxfiles)
                if(.not.allocated(boxfiles))then
                    write(logfhandle,*)'No box file found in ', trim(params%dir_box), '; simple_commander_preprocess::exec_extract 1'
                    THROW_HARD('No box file found; exec_extract, 1')
                endif
                if(size(boxfiles)==0)then
                    write(logfhandle,*)'No box file found in ', trim(params%dir_box), '; simple_commander_preprocess::exec_extract 2'
                    THROW_HARD('No box file found; exec_extract 2')
                endif
            else
                write(logfhandle,*)'Directory does not exist: ', trim(params%dir_box), 'simple_commander_preprocess::exec_extract'
                THROW_HARD('box directory does not exist; exec_extract')
            endif
            call cline%set('dir_box', params%dir_box)
        endif
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! box input
            if( cline%defined('dir_box') )then
                boxfile_name = boxfile_from_mic(mic_name)
                if(trim(boxfile_name).eq.NIL)cycle
            else
                call o_mic%getter('boxfile', boxfile_name)
                if( .not.file_exists(boxfile_name) )cycle
            endif
            ! get number of frames from stack
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract')
            ! update counter
            nmics = nmics + 1
        enddo
        if( nmics == 0 ) THROW_HARD('No particles to extract! exec_extract')
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! ASSEMBLY
        allocate(parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1,params%nparts
            parts_fname(ipart) = trim(ALGN_FBODY)//int2str_pad(ipart,numlen)//trim(METADATA_EXT)
        enddo
        ! copy updated micrographs
        cnt   = 0
        nstks = 0
        do ipart = 1,params%nparts
            call spproj_part%read_segment('mic',parts_fname(ipart))
            do imic = 1,spproj_part%os_mic%get_noris()
                cnt = cnt + 1
                call spproj_part%os_mic%get_ori(imic, o_mic)
                call spproj%os_mic%set_ori(cnt,o_mic)
                if( nint(o_mic%get('nptcls')) > 0 ) nstks = nstks + 1
            enddo
            call spproj_part%kill
        enddo
        if( cnt /= nmics_tot ) THROW_HARD('Inconstistent number of micrographs in individual projects')
        ! fetch stacks table
        if( nstks > 0 )then
            call os_stk%new(nstks)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_part%read_segment('stk',parts_fname(ipart))
                do istk = 1,spproj_part%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj_part%os_stk%get_ori(istk, o_tmp)
                    call os_stk%set_ori(cnt,o_tmp)
                    stktab(cnt) = os_stk%get_static(cnt,'stk')
                enddo
                call spproj_part%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            ! transfer particles locations to ptcl2D
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_part%read_segment('ptcl2D',parts_fname(ipart))
                do i = 1,spproj_part%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    call spproj_part%get_boxcoords(i, boxcoords)
                    call spproj%set_boxcoords(cnt, boxcoords)
                enddo
                call spproj_part%kill
            enddo
            call os_stk%kill
        endif
        ! final write
        call spproj%write
        ! clean
        call o_mic%kill
        call o_tmp%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_EXTRACT_DISTR NORMAL STOP ****')

        contains

            character(len=LONGSTRLEN) function boxfile_from_mic(mic)
                character(len=*), intent(in) :: mic
                character(len=LONGSTRLEN)    :: box_from_mic
                integer :: ibox
                box_from_mic     = fname_new_ext(basename(mic),'box')
                boxfile_from_mic = NIL
                do ibox=1,size(boxfiles)
                    if(trim(basename(boxfiles(ibox))).eq.trim(box_from_mic))then
                        boxfile_from_mic = trim(boxfiles(ibox))
                        return
                    endif
                enddo
            end function boxfile_from_mic

    end subroutine exec_extract_distr

    !> for extracting particle images from integrated DDD movies
    subroutine exec_reextract_distr( self, cline )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        class(reextract_distr_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj
        type(sp_project),           allocatable :: spproj_parts(:)
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o, o_tmp
        type(oris)                              :: os_stk
        character(len=LONGSTRLEN),  allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        character(len=:),           allocatable :: mic_name, imgkind
        integer :: boxcoords(2)
        integer :: imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        call cline%set('nthr',1.)
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read(params%projfile)
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        if( spproj%get_nstks() == 0 ) THROW_HARD('This project file does not contain stacks!')
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot < params%nparts )then
            params%nparts = nmics_tot
        endif
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! update counter
            nmics = nmics + 1
            ! removes boxfile from micrographs
            call spproj%os_mic%delete_entry(imic,'boxfile')
        enddo
        if( nmics == 0 ) THROW_HARD('No particles to re-extract! exec_extract')
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! ASSEMBLY
        call spproj%os_mic%kill
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        allocate(spproj_parts(params%nparts),parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1,params%nparts
            parts_fname(ipart) = trim(ALGN_FBODY)//int2str_pad(ipart,numlen)//trim(METADATA_EXT)
        enddo
        ! copy updated micrographs
        cnt   = 0
        nmics = 0
        do ipart = 1,params%nparts
            call spproj_parts(ipart)%read_segment('mic',parts_fname(ipart))
            nmics = nmics + spproj_parts(ipart)%os_mic%get_noris()
        enddo
        if( nmics > 0 )then
            call spproj%os_mic%new(nmics)
            ! transfer stacks
            cnt   = 0
            nstks = 0
            do ipart = 1,params%nparts
                do imic = 1,spproj_parts(ipart)%os_mic%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_mic%get_ori(imic,o_tmp)
                    call spproj%os_mic%set_ori(cnt,o_tmp)
                enddo
                call spproj_parts(ipart)%kill
                call spproj_parts(ipart)%read_segment('stk',parts_fname(ipart))
                nstks = nstks + spproj_parts(ipart)%os_stk%get_noris()
            enddo
            if( nstks /= nmics ) THROW_HARD('Inconstistent number of stacks in individual projects')
            ! generates stacks table
            call os_stk%new(nstks)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                do istk = 1,spproj_parts(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_stk%get_ori(istk, o_tmp)
                    call os_stk%set_ori(cnt,o_tmp)
                    stktab(cnt) = os_stk%get_static(cnt,'stk')
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            ! transfer 2D parameters
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl2D',parts_fname(ipart))
                do i = 1,spproj_parts(ipart)%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    ! particles coordinates
                    call spproj_parts(ipart)%get_boxcoords(i, boxcoords)
                    call spproj%set_boxcoords(cnt, boxcoords)
                    ! search history & parameters
                    call spproj_parts(ipart)%os_ptcl2D%get_ori(i, o)
                    call spproj%os_ptcl2D%transfer_2Dparams(cnt, o)
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! transfer 3D parameters
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl3D',parts_fname(ipart))
                do i = 1,spproj_parts(ipart)%os_ptcl3D%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_ptcl3D%get_ori(i, o)
                    call spproj%os_ptcl3D%transfer_3Dparams(cnt, o)
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! clean-up
            call os_stk%kill
        endif
        ! final write
        call spproj%write
        ! clean-up
        call qsys_cleanup
        call o_mic%kill
        call o%kill
        call o_tmp%kill
        ! end gracefully
        call simple_end('**** SIMPLE_REEXTRACT_DISTR NORMAL STOP ****')

        contains

            character(len=LONGSTRLEN) function boxfile_from_mic(mic)
                character(len=*), intent(in) :: mic
                character(len=LONGSTRLEN)    :: box_from_mic
                integer :: ibox
                box_from_mic     = fname_new_ext(basename(mic),'box')
                boxfile_from_mic = NIL
                do ibox=1,size(boxfiles)
                    if(trim(basename(boxfiles(ibox))).eq.trim(box_from_mic))then
                        boxfile_from_mic = trim(boxfiles(ibox))
                        return
                    endif
                enddo
            end function boxfile_from_mic

    end subroutine exec_reextract_distr

    subroutine exec_motion_correct_distr( self, cline )
        class(motion_correct_distr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   20.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     6.)
        call cline%set('oritype', 'mic')
        call params%new(cline)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nmovies() ==0 )then
            THROW_HARD('no movie to process! exec_motion_correct_distr')
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct_distr

    subroutine exec_gen_pspecs_and_thumbs_distr( self, cline )
        class(gen_pspecs_and_thumbs_distr_commander), intent(inout) :: self
        class(cmdline),                               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: nintgs
        call cline%set('oritype', 'mic')
        call params%new(cline)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nintgs = spproj%get_nintgs()
        if( nintgs ==0 )then
            THROW_HARD('no integrated movies to process! exec_gen_pspecs_and_thumbs_distr')
        endif
        if( params%nparts > nintgs )then
            call cline%set('nparts', real(nintgs))
            params%nparts = nintgs
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs_distr

    subroutine exec_motion_correct_tomo_distr( self, cline )
        use simple_oris, only: oris
        class(motion_correct_tomo_distr_commander), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: tomonames(:)
        type(parameters)         :: params
        type(oris)               :: exp_doc
        integer                  :: nseries, ipart
        type(qsys_env)           :: qenv
        character(len=KEYLEN)    :: str
        type(chash)              :: job_descr
        type(chash), allocatable :: part_params(:)
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   20.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     6.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',    'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        call cline%set('prg', 'motion_correct')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('tomoseries') )then
            call read_filetable(params%tomoseries, tomonames)
        else
            THROW_HARD('need tomoseries (filetable of filetables) to be part of the command line when tomo=yes')
        endif
        nseries = size(tomonames)
        call exp_doc%new(nseries)
        if( cline%defined('exp_doc') )then
            if( file_exists(params%exp_doc) )then
                call exp_doc%read(params%exp_doc)
            else
                THROW_HARD('the required parameter file (flag exp_doc): '//trim(params%exp_doc)//' not in cwd')
            endif
        else
            THROW_HARD('need exp_doc (line: exp_time=X dose_rate=Y) to be part of the command line when tomo=yes')
        endif
        params%nparts = nseries
        params%nptcls = nseries
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts), stat=alloc_stat) ! -1. is default excluded value
        if(alloc_stat.ne.0)call allocchk("simple_commander_distr_wflows::motion_correct_tomo_moview_distr ", alloc_stat)
        do ipart=1,params%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('filetab', trim(tomonames(ipart)))
            call part_params(ipart)%set('fbody', 'tomo'//int2str_pad(ipart,params%numlen_tomo))
            str = real2str(exp_doc%get(ipart,'exp_time'))
            call part_params(ipart)%set('exp_time', trim(str))
            str = real2str(exp_doc%get(ipart,'dose_rate'))
            call part_params(ipart)%set('dose_rate', trim(str))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT_TOMO NORMAL STOP ****')
    end subroutine exec_motion_correct_tomo_distr

    subroutine exec_ctf_estimate_distr( self, cline )
        class(ctf_estimate_distr_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(chash)                   :: job_descr
        type(qsys_env)                :: qenv
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       10.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',       2.5)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('no micrograph to process! exec_ctf_estimate_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call qsys_cleanup
        ! graceful ending
        call simple_end('**** SIMPLE_DISTR_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate_distr

    subroutine exec_pick_distr( self, cline )
        class(pick_distr_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(cmdline)    :: cline_make_pickrefs
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        if( .not.cline%defined('refs') .and. .not.cline%defined('vol1') )then
            THROW_HARD('one of REFS and VOL1 must be provided!')
        endif
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'mic')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('No micrograph to process! exec_pick_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        cline_make_pickrefs = cline
        call cline_make_pickrefs%set('prg','make_pickrefs')
        call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
        call cline%set('refs', trim(PICKREFS)//params%ext)
        call cline%delete('vol1')
        write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call qsys_cleanup
        ! graceful exit
        call simple_end('**** SIMPLE_DISTR_PICK NORMAL STOP ****')
    end subroutine exec_pick_distr

    subroutine exec_make_cavgs_distr( self, cline )
        class(make_cavgs_distr_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(cmdline)    :: cline_cavgassemble
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure the use of all resources in assembly
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble class averages
        call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_cluster2D_distr( self, cline )
        use simple_procimgfile
        use simple_commander_cluster2D, only: check_2Dconv_commander
        class(cluster2D_distr_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(check_2Dconv_commander)     :: xcheck_2Dconv
        type(make_cavgs_distr_commander) :: xmake_cavgs
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        character(len=LONGSTRLEN) :: refs, refs_even, refs_odd, str, str_iter
        integer                   :: iter
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        logical                   :: l_stream
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        ! streaming gets its own logics because it is an exception to rules in parameters
        l_stream = .false.
        if( cline%defined('stream') ) l_stream = trim(cline%get_carg('stream')).eq.'yes'
        call cline%set('stream','no')
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        if( build%spproj%get_nptcls() == 0 )then
            THROW_HARD('no particles found! exec_cluster2D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        if( l_stream )then
            call job_descr%set('stream','yes')
            ! need to be explicity defined for streaming
            call job_descr%set('box',   int2str(params%box))
            call job_descr%set('smpd',  real2str(params%smpd))
            call job_descr%set('nptcls',int2str(params%nptcls))
        endif
        ! splitting
        call build%spproj%split_stk(params%nparts)
        ! prepare command lines from prototype master
        cline_check_2Dconv = cline
        cline_cavgassemble = cline
        cline_make_cavgs   = cline ! ncls is transferred here
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure use of all resources in assembly
        call cline_make_cavgs%set('prg',   'make_cavgs')
        if( l_stream )call cline_check_2Dconv%set('stream','yes')
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs' // params%ext
            params%refs      = trim(refs)
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
            if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                if( params%tseries .eq. 'yes' )then
                    call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                else
                    call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                endif
                call copy_imgfile(trim(params%refs), trim(params%refs_even), params%smpd, [1,params%ncls])
                call copy_imgfile(trim(params%refs), trim(params%refs_odd),  params%smpd, [1,params%ncls])
            else
                call cline_make_cavgs%set('refs', params%refs)
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            refs = trim(params%refs)
        endif
        ! variable neighbourhood size
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! deal with eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! main loop
        iter = params%startit - 1
        do
            iter = iter + 1
            str_iter = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            ! updates
            call job_descr%set('refs', trim(refs))
            call job_descr%set('startit', int2str(iter))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', trim(FRCS_FILE))
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
            ! merge orientation documents
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // params%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // params%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // params%ext
            call cline_cavgassemble%set('refs', trim(refs))
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
            ! check convergence
            call xcheck_2Dconv%execute(cline_check_2Dconv)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check_2Dconv%get_rarg('frac')
            ! the below activates shifting & automasking
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check_2Dconv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check_2Dconv%get_rarg('trs'))
                    call job_descr%set('trs', trim(str) )
                endif
            endif
            if( cline_check_2Dconv%get_carg('converged').eq.'yes' .or. iter==params%maxits )then
                if( cline_check_2Dconv%get_carg('converged').eq.'yes' )call cline%set('converged','yes')
                exit
            endif
        end do
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_distr

    subroutine exec_refine3D_distr( self, cline )
        use simple_commander_refine3D, only: check_3Dconv_commander
        use simple_commander_volops,   only: postprocess_commander
        class(refine3D_distr_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(reconstruct3D_distr_commander) :: xreconstruct3D_distr
        type(check_3Dconv_commander)        :: xcheck_3Dconv
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline)    :: cline_reconstruct3D_distr
        type(cmdline)    :: cline_check_3Dconv
        type(cmdline)    :: cline_volassemble
        type(cmdline)    :: cline_postprocess
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        character(len=:),          allocatable :: vol_fname, prev_refine_path, target_name
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        integer,                   allocatable :: state_pops(:)
        character(len=STDLEN)     :: vol, vol_even, vol_odd, vol_iter, vol_iter_even
        character(len=STDLEN)     :: vol_iter_odd, str, str_iter, optlp_file
        character(len=STDLEN)     :: str_state, fsc_file
        character(len=LONGSTRLEN) :: volassemble_output
        real    :: corr, corr_prev, smpd
        integer :: i, state, iter, iostat, box, nfiles, niters, iter_switch2euclid
        logical :: err, vol_defined, have_oris, do_abinitio, converged, fall_over
        logical :: l_projection_matching, l_switch2euclid, l_continue, l_matchfilt_ini
        if( .not. cline%defined('refine') )then
            call cline%set('refine',  'single')
        else
            if( cline%get_carg('refine').eq.'multi' .and. .not. cline%defined('nstates') )then
                THROW_HARD('refine=MULTI requires specification of NSTATES')
            endif
        endif
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp', 30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! objfun=euclid logics, part 1
        l_switch2euclid  = .false.
        if( cline%defined('objfun') )then
            l_continue = .false.
            if( cline%defined('continue') ) l_continue = trim(cline%get_carg('continue')).eq.'yes'
            if( (trim(cline%get_carg('objfun')).eq.'euclid') .and. .not.l_continue )then
                l_switch2euclid = .true.
                call cline%set('objfun','cc')
            endif
        endif
        ! init
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            write(logfhandle,*)'Unsupported ORITYPE; simple_commander_distr_wflows::exec_refine3D_distr'
        end select
        if( fall_over )then
            THROW_HARD('no particles found! :exec_refine3D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' )then
            call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        endif
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr = cline
        cline_check_3Dconv        = cline
        cline_volassemble         = cline
        cline_postprocess         = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' ) ! required for distributed call
        call cline_postprocess%set('prg', 'postprocess' )            ! required for local call
        if( trim(params%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set( 'pgrp', 'c1' )
        call cline_postprocess%set('mirr',    'no')
        call cline_postprocess%set('mkdir',   'no')
        call cline_postprocess%set('imgkind', 'vol')
        if( trim(params%oritype).eq.'cls3D' ) call cline_postprocess%set('imgkind', 'vol_cavg')
        call cline_volassemble%set('nthr', 0.)  ! to ensure use of all resources in assembly
        ! for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates) , stat=alloc_stat)
        if(alloc_stat /= 0)call allocchk("simple_commander_distr_wflows::exec_refine3D_distr state_assemble ",alloc_stat)
        ! removes unnecessary volume keys and generates volassemble finished names
        do state = 1,params%nstates
            vol = 'vol'//int2str( state )
            call cline_check_3Dconv%delete( vol )
            call cline_volassemble%delete( vol )
            call cline_postprocess%delete( vol )
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        DebugPrint ' In exec_refine3D_distr; begin starting models'
        ! E/O PARTITIONING
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! GENERATE STARTING MODELS & ORIENTATIONS
        if( params%continue .eq. 'yes' )then
            ! we are continuing from a previous refinement round,
            ! i.e. projfile is fetched from a X_refine3D dir
            ! set starting volume(s), iteration number & previous refinement path
            do state=1,params%nstates
                ! volume(s)
                vol = 'vol' // int2str(state)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%get_vol('vol_cavg', state, vol_fname, smpd, box)
                else
                    call build%spproj%get_vol('vol', state, vol_fname, smpd, box)
                endif
                call cline%set(trim(vol), vol_fname)
                params%vols(state) = vol_fname
                if( state == 1 )then
                    ! get the iteration number
                    iter = fname2iter(basename(vol_fname))
                    ! startit becomes the next
                    params%startit = iter + 1
                    call cline%set('startit', real(params%startit))
                endif
            end do
            prev_refine_path = get_fpath(vol_fname)
            ! carry over FRCs/FSCs
            ! one FSC file per state
            do state=1,params%nstates
                str_state = int2str_pad(state,2)
                fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                call simple_copy_file(trim(prev_refine_path)//trim(fsc_file), fsc_file)
            end do
            ! one FRC file for all states
            call simple_copy_file(trim(prev_refine_path)//trim(FRCS_FILE), trim(FRCS_FILE))
            ! carry over the oridistributions_part* files
            call simple_list_files(prev_refine_path//'oridistributions_part*', list)
            nfiles = size(list)
            err    = params%nparts /= nfiles
            if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
            do i=1,nfiles
                target_name = PATH_HERE//basename(trim(list(i)))
                call simple_copy_file(trim(list(i)), target_name)
            end do
            deallocate(list)
            ! if we are doing fractional volume update, partial reconstructions need to be carried over
            if( params%l_frac_update )then
                call simple_list_files(prev_refine_path//'*recvol_state*part*', list)
                nfiles = size(list)
                err = params%nparts * 4 /= nfiles
                if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
                deallocate(list)
            endif
            ! if we are doing objfun=euclid the sigm estimates need to be carried over
            if( trim(params%objfun) .eq. 'euclid' )then
                call simple_list_files(prev_refine_path//trim(SIGMA2_FBODY)//'*', list)
                nfiles = size(list)
                if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
                deallocate(list)
            endif
        endif
        vol_defined = .false.
        do state = 1,params%nstates
            vol = 'vol' // int2str(state)
            if( cline%defined(trim(vol)) ) vol_defined = .true.
        enddo
        have_oris   = .not. build%spproj%is_virgin_field(params%oritype)
        do_abinitio = .not. have_oris .and. .not. vol_defined
        if( do_abinitio )then
            call build%spproj_field%rnd_oris
            call build%spproj_field%zero_shifts
            have_oris = .true.
            call build%spproj%write_segment_inside(params%oritype)
        endif
        l_projection_matching = .false.
        if( have_oris .and. .not. vol_defined )then
            ! reconstructions needed
            call xreconstruct3D_distr%execute( cline_reconstruct3D_distr )
            do state = 1,params%nstates
                ! rename volumes and update cline
                str_state = int2str_pad(state,2)
                vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//params%ext
                iostat = simple_rename( trim(vol), trim(str) )
                vol = 'vol'//trim(int2str(state))
                call cline%set( trim(vol), trim(str) )
                vol_even = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params%ext
                iostat= simple_rename( trim(vol_even), trim(str) )
                vol_odd  = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params%ext
                iostat =  simple_rename( trim(vol_odd), trim(str) )
            enddo
        else if( vol_defined .and. params%continue .ne. 'yes' )then
            ! projection matching
            l_projection_matching = .true.
            if( .not. have_oris )then
                select case( params%neigh )
                    case( 'yes' )
                        THROW_HARD('refinement method requires input orientations')
                    case DEFAULT
                        ! all good
                end select
            endif
            if( .not.cline%defined('lp') ) THROW_HARD('LP needs be defined for the first step of projection matching!')
            call cline%delete('update_frac')
            if( params%neigh .ne. 'yes' )then
                ! this forces the first round of alignment on the starting model(s)
                ! to be greedy and the subseqent ones to be whatever the refinement flag is set to
                call build%spproj%os_ptcl3D%delete_3Dalignment(keepshifts=.true.)
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! EXTREMAL DYNAMICS
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! objfun=euclid logics, part 2
        iter_switch2euclid = -1
        if( l_switch2euclid )then
            iter_switch2euclid = 1
            if( cline%defined('update_frac') ) iter_switch2euclid = ceiling(1./(params%update_frac+0.001))
            call cline%set('needs_sigma','yes')
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        niters = 0
        iter   = params%startit - 1
        corr   = -1.
        do
            niters            = niters + 1
            iter              = iter + 1
            params%which_iter = iter
            str_iter          = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
            if( have_oris .or. iter > params%startit )then
                call build%spproj%read(params%projfile)
                if( params%refine .eq. 'snhc' )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = build%spproj_field%get_avg('corr')
                    if( iter > 1 .and. corr <= corr_prev )then
                        params%szsn = min(SZSN_MAX,params%szsn + SZSN_STEP)
                    endif
                    call job_descr%set('szsn', int2str(params%szsn))
                    call cline%set('szsn', real(params%szsn))
                endif
            endif
            ! exponential cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            call job_descr%set('which_iter', trim(int2str(params%which_iter)))
            call cline%set('which_iter', real(params%which_iter))
            call job_descr%set( 'startit', trim(int2str(iter)))
            call cline%set('startit', real(iter))
            ! FRCs
            if( cline%defined('frcs') )then
                ! all good
            else
                call job_descr%set('frcs', trim(FRCS_FILE))
            endif
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
            ! ASSEMBLE ALIGNMENT DOCS
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY)
            ! ASSEMBLE VOLUMES
            select case(trim(params%refine))
            case('eval')
                ! nothing to do
            case DEFAULT
                call cline_volassemble%set( 'prg', 'volassemble' ) ! required for cmdline exec
                do state = 1,params%nstates
                    str_state = int2str_pad(state,2)
                    volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
                    call cline_volassemble%set( 'state', real(state) )
                    if( params%nstates>1 )call cline_volassemble%set('part', real(state))
                    call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                    &'simple_script_state'//trim(str_state), volassemble_output)
                end do
                call qsys_watcher(state_assemble_finished)
                ! rename & add volumes to project & update job_descr
                call build%spproj_field%get_pops(state_pops, 'state')
                do state = 1,params%nstates
                    str_state = int2str_pad(state,2)
                    if( state_pops(state) == 0 )then
                        ! cleanup for empty state
                        vol = 'vol'//trim(int2str(state))
                        call cline%delete( vol )
                        call job_descr%delete( vol )
                    else
                        ! rename state volume
                        vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                        if( params%refine .eq. 'snhc' )then
                            vol_iter = trim(SNHCVOL)//trim(str_state)//params%ext
                        else
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//params%ext
                        endif
                        iostat = simple_rename( vol, vol_iter )
                        vol_even      = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                        vol_odd       = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                        vol_iter_even = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_even'//params%ext
                        vol_iter_odd  = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_odd' //params%ext
                        iostat        = simple_rename( vol_even, vol_iter_even )
                        iostat        = simple_rename( vol_odd,  vol_iter_odd  )
                        fsc_file      = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                        optlp_file    = ANISOLP_FBODY//trim(str_state)//params%ext
                        ! add filters to os_out
                        call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                        call build%spproj%add_vol2os_out(optlp_file, params%smpd, state, 'vol_filt', box=params%box)
                        ! add state volume to os_out
                        if( trim(params%oritype).eq.'cls3D' )then
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol_cavg')
                        else
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol')
                        endif
                        ! updates cmdlines & job description
                        vol = 'vol'//trim(int2str(state))
                        call job_descr%set( vol, vol_iter )
                        call cline%set(vol, vol_iter )
                    endif
                enddo
                ! volume mask, one for all states
                if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                ! writes os_out
                call build%spproj%write_segment_inside('out')
                ! per state post-process
                do state = 1,params%nstates
                    if( state_pops(state) == 0 )cycle
                    call cline_postprocess%set('state', real(state))
                    if( cline%defined('lp') ) call cline_postprocess%set('lp', params%lp)
                    call xpostprocess%execute(cline_postprocess)
                enddo
            end select
            ! CONVERGENCE
            converged = .false.
            select case(trim(params%refine))
                case('eval')
                    ! nothing to do
                case DEFAULT
                    if( str_has_substr(params%refine,'cluster')) call cline_check_3Dconv%delete('update_res')
                    call xcheck_3Dconv%execute(cline_check_3Dconv)
                    if( iter >= params%startit + 2 )then
                        ! after a minimum of 2 iterations
                        if( cline_check_3Dconv%get_carg('converged') .eq. 'yes' ) converged = .true.
                    endif
            end select
            if( iter >= params%maxits ) converged = .true.
            if( converged )then
                ! safest to write the whole thing here as multiple fields updated
                call build%spproj%write
                exit ! main loop
            endif
            ! ITERATION DEPENDENT UPDATES
            if( cline_check_3Dconv%defined('trs') .and. .not.job_descr%isthere('trs') )then
                ! activates shift search if frac >= 90
                str = real2str(cline_check_3Dconv%get_rarg('trs'))
                call job_descr%set( 'trs', trim(str) )
                call cline%set( 'trs', cline_check_3Dconv%get_rarg('trs') )
            endif
            if( l_projection_matching .and. cline%defined('lp_iters') .and. (niters == params%lp_iters ) )then
                ! e/o projection matching
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO EVEN/ODD RESOLUTION LIMIT'
                write(logfhandle,'(A)')'>>>'
                l_projection_matching = .false.
                call cline%delete('lp')
                call job_descr%delete('lp')
                call cline_postprocess%delete('lp')
                if( params%l_frac_update )then
                    call job_descr%set('update_frac', real2str(params%update_frac))
                    call cline%set('update_frac', params%update_frac)
                    call cline_check_3Dconv%set('update_frac', params%update_frac)
                    call cline_volassemble%set('update_frac', params%update_frac)
                endif
            endif
            ! objfun=euclid, part 3: actual switch
            if( l_switch2euclid .and. niters.eq.iter_switch2euclid )then
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO OBJFUN=EUCLID'
                call cline%set('objfun','euclid')
                call cline%set('match_filt','no')
                call job_descr%set('objfun','euclid')
                call job_descr%set('match_filt','no')
                call cline_volassemble%set('objfun','euclid')
                params%objfun    = 'euclid'
                params%cc_objfun = OBJFUN_EUCLID
                l_switch2euclid  = .false.
            endif
        end do
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D_distr

    subroutine exec_reconstruct3D_distr( self, cline )
        class(reconstruct3D_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=:),          allocatable :: target_name
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        character(len=LONGSTRLEN) :: refine_path
        character(len=STDLEN)     :: volassemble_output, str_state, fsc_file, optlp_file
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        integer          :: state, ipart
        logical          :: fall_over
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        ! soft reconstruction from o_peaks in dir_refine?
        if( params%l_rec_soft )then
            call make_relativepath(CWD_GLOB,params%dir_refine,refine_path)
            call simple_list_files(trim(refine_path)//'/oridistributions_part*', list)
            if( size(list) == 0 )then
                THROW_HARD('No oridistributions can be found in '//trim(params%dir_refine))
            elseif( size(list) /= params%nparts )then
                THROW_HARD('# partitions not consistent with that in '//trim(params%dir_refine))
            endif
            ! copy the orientation peak distributions
            if( trim(params%dir_refine).eq.trim(CWD_GLOB) )then
                ! already here
            else
                do ipart=1,params%nparts
                    target_name = PATH_HERE//basename(trim(list(ipart)))
                    call simple_copy_file(trim(list(ipart)), target_name)
                end do
            endif
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' )then
            call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        endif
        ! eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble volumes
        ! this is for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates) )
        do state = 1, params%nstates
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        cline_volassemble = cline
        call cline_volassemble%set('prg', 'volassemble')
        call cline_volassemble%set('nthr', 0.) ! to ensure the use of all resources in assembly
        ! parallel assembly
        do state = 1,params%nstates
            str_state = int2str_pad(state,2)
            volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
            call cline_volassemble%set( 'state', real(state) )
            if( params%nstates>1 )call cline_volassemble%set('part', real(state))
            call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
            'simple_script_state'//trim(str_state), trim(volassemble_output))
        end do
        call qsys_watcher(state_assemble_finished)
        ! updates project file only if called from another workflow
        if( params%mkdir.eq.'yes' )then
            do state = 1,params%nstates
                fsc_file      = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                optlp_file    = ANISOLP_FBODY//trim(str_state)//params%ext
                call build%spproj%add_fsc2os_out(trim(fsc_file), state, params%box)
                call build%spproj%add_vol2os_out(trim(optlp_file), params%smpd, state, 'vol_filt', box=params%box)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd, state, 'vol_cavg')
                else
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd, state, 'vol')
                endif
            enddo
            call build%spproj%write_segment_inside('out',params%projfile)
            if( params%l_rec_soft )then
                ! ptcl3D segment may have been updated and needs to be written
                call build%spproj%write_segment_inside('ptcl3D',params%projfile)
            endif
        endif
        ! termination
        call qsys_cleanup
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D_distr

    subroutine exec_tseries_track_distr( self, cline )
        use simple_nrtxtfile,         only: nrtxtfile
        class(tseries_track_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, alloc_stat, j, orig_box, ipart
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')      ) call cline%set('neg',     'yes')
        if( .not. cline%defined('lp')       ) call cline%set('lp',       2.0)
        if( .not. cline%defined('lp_backgr')) call cline%set('lp_backgr',1.1)
        if( .not. cline%defined('hp')       ) call cline%set('hp',       5.0)
        if( .not. cline%defined('width')    ) call cline%set('width',    1.1)
        if( .not. cline%defined('cenlp')    ) call cline%set('cenlp',    5.0)
        if( .not. cline%defined('ctf')      ) call cline%set('ctf',      'no')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        !if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
        if( nlines(params%boxfile) > 0 )then
            call boxfile%new(params%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track', alloc_stat)
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    THROW_HARD('Only square windows allowed!')
                endif
            end do
        else
            THROW_HARD('inputted boxfile is empty; exec_tseries_track')
        endif
        call boxfile%kill
        call cline%delete('boxfile')
        params%nptcls = ndatlines
        params%nparts = params%nptcls
        !if( params%ncunits > params%nparts )&
        !&THROW_HARD('# computational units (ncunits) mjust be <= number of entries in boxfiles')
        ! box and numlen need to be part of command line
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(3)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('ind',    int2str(ipart))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track_distr

    recursive subroutine exec_scale_project_distr( self, cline )
        class(scale_project_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(scale_project_distr_commander) :: xscale_distr
        type(qsys_env)                     :: qenv
        type(chash)                        :: job_descr
        type(cmdline)                      :: cline_scale
        type(chash),               allocatable :: part_params(:)
        character(len=LONGSTRLEN), allocatable :: part_stks(:)
        type(parameters)              :: params
        type(builder)                 :: build
        character(len=:), allocatable :: projfile_sc
        character(len=STDLEN) :: filetab
        integer, allocatable  :: parts(:,:)
        real                  :: smpd, smpd_target
        integer               :: istk, ipart, nparts, nstks, cnt, partsz, box, newbox
        logical               :: gen_sc_project
        ! mkdir=yes: a new *_sc project + stacks are generated
        ! mkdir=no : only stacks are scaled
        gen_sc_project = cline%get_carg('mkdir').eq.'yes'
        ! make parameters and project
        call build%init_params_and_build_spproj(cline, params)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! copy command line
        cline_scale = cline
        ! prepare part-dependent parameters
        nstks = build%spproj%os_stk%get_noris()
        if( nstks == 0 ) THROW_HARD('os_stk field of spproj empty; exec_scale_distr')
        if( cline%defined('nparts') )then
            nparts = min(params%nparts, nstks)
            call cline_scale%set('nparts', real(nparts))
        else
            nparts = 1
        endif
        smpd = build%spproj%get_smpd()
        box  = build%spproj%get_box()
        call cline_scale%set('smpd', smpd)
        call cline_scale%set('box',  real(box))
        if( gen_sc_project )then
            ! make new project & scales
            smpd_target = max(smpd, smpd * real(box)/real(params%newbox))
            call simple_mkdir(filepath(PATH_PARENT,'stack_parts_sc'), errmsg="commander_distr_wflows::exec_scale_project_distr ")
            call build%spproj%scale_projfile(smpd_target, projfile_sc, cline, cline_scale,&
                dir=filepath(PATH_PARENT,'stack_parts_sc'))
            newbox = nint(cline_scale%get_rarg('newbox'))
            if( newbox == box )then
                write(logfhandle,*)'Inconsistent input dimensions: from ',box,' to ',newbox
                THROW_HARD('inconsistent input dimensions; exec_scale_project_distr')
            endif
            call cline_scale%set('newbox', real(newbox))
            call xscale_distr%execute( cline_scale )
            ! delete copy in working directory
            call del_file(params%projfile)
        else
            ! scaling only
            params%nparts = nparts
            parts = split_nobjs_even(nstks, nparts)
            allocate(part_params(nparts))
            cnt = 0
            do ipart=1,nparts
                call part_params(ipart)%new(1)
                partsz = parts(ipart,2) - parts(ipart,1) + 1
                allocate(part_stks(partsz))
                ! creates part filetab
                filetab = 'scale_stktab_part'//int2str(ipart)//trim(TXT_EXT)
                do istk=1,partsz
                    cnt = cnt + 1
                    part_stks(istk) = build%spproj%get_stkname(cnt)
                enddo
                ! write part filetab & update part parameters
                call write_filetable( filetab, part_stks )
                call part_params(ipart)%set('filetab', filetab)
                deallocate(part_stks)
            end do
            deallocate(parts)
            ! setup the environment for distributed execution
            call qenv%new(nparts)
            ! prepare job description
            call cline_scale%gen_job_descr(job_descr)
            call job_descr%set('prg', 'scale')
            call job_descr%set('autoscale', 'no')
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
            ! clean
            call qsys_cleanup
            ! removes temporary split stktab lists
            do ipart=1,nparts
                filetab = 'scale_stktab_part'//int2str(ipart)//trim(TXT_EXT)
                call del_file( filetab )
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_SCALE NORMAL STOP ****')
    end subroutine exec_scale_project_distr

end module simple_commander_distr_wflows
