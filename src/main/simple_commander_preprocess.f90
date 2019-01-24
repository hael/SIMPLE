! concrete commander: pre-processing routines
module simple_commander_preprocess
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_qsys_funs,      only: qsys_job_finished
use simple_sp_project,     only: sp_project
use simple_commander_base, only: commander_base
implicit none

public :: preprocess_commander
public :: motion_correct_commander
public :: gen_pspecs_and_thumbs_commander
public :: ctf_estimate_commander
public :: map_cavgs_selection_commander
public :: pick_commander
public :: pick_commander_chiara
public :: extract_commander
public :: reextract_commander
public :: pick_extract_commander
public :: make_pickrefs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_commander
  contains
    procedure :: execute      => exec_preprocess
end type preprocess_commander
type, extends(commander_base) :: motion_correct_commander
  contains
    procedure :: execute      => exec_motion_correct
end type motion_correct_commander
type, extends(commander_base) :: gen_pspecs_and_thumbs_commander
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs
end type gen_pspecs_and_thumbs_commander
type, extends(commander_base) :: ctf_estimate_commander
  contains
    procedure :: execute      => exec_ctf_estimate
end type ctf_estimate_commander
type, extends(commander_base) :: map_cavgs_selection_commander
  contains
    procedure :: execute      => exec_map_cavgs_selection
end type map_cavgs_selection_commander
type, extends(commander_base) :: pick_commander
  contains
    procedure :: execute      => exec_pick
end type pick_commander
type, extends(commander_base) :: pick_commander_chiara
  contains
    procedure :: execute      => exec_pick_chiara
end type pick_commander_chiara
type, extends(commander_base) :: extract_commander
  contains
    procedure :: execute      => exec_extract
end type extract_commander
type, extends(commander_base) :: reextract_commander
  contains
    procedure :: execute      => exec_reextract
end type reextract_commander
type, extends(commander_base) :: pick_extract_commander
  contains
    procedure :: execute      => exec_pick_extract
end type pick_extract_commander
type, extends(commander_base) :: make_pickrefs_commander
  contains
    procedure :: execute      => exec_make_pickrefs
end type make_pickrefs_commander

contains

    subroutine exec_preprocess( self, cline )
        use simple_ori,                 only: ori
        use simple_sp_project,          only: sp_project
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        use simple_picker_iter,         only: picker_iter
        use simple_binoris_io,          only: binwrite_oritab
        class(preprocess_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)              :: params
        type(ori)                     :: o_mov
        type(ctf_estimate_iter)       :: cfiter
        type(motion_correct_iter)     :: mciter
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        type(ctfparams)               :: ctfvars
        character(len=:), allocatable :: imgkind, moviename, output_dir_picker, fbody
        character(len=:), allocatable :: moviename_forctf, moviename_intg, output_dir_motion_correct
        character(len=:), allocatable :: output_dir_ctf_estimate, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, nptcls_out
        logical :: l_pick
        call cline%set('oritype', 'mic')
        call params%new(cline)
        if( params%scale > 1.05 )then
            THROW_HARD('scale cannot be > 1; exec_preprocess')
        endif
        if( params%tomo .eq. 'yes' )then
            THROW_HARD('tomography mode (tomo=yes) not yet supported!')
        endif
        l_pick = cline%defined('refs')
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nmovies()==0 .and. spproj%get_nintgs()==0 ) THROW_HARD('No movie/micrograph to process!')
        ! output directories & naming
        output_dir_ctf_estimate   = PATH_HERE
        output_dir_motion_correct = PATH_HERE
        if( l_pick )then
            output_dir_picker  = PATH_HERE
        endif
        if( params%stream.eq.'yes' )then
            output_dir_ctf_estimate   = trim(DIR_CTF_ESTIMATE)
            output_dir_motion_correct = trim(DIR_MOTION_CORRECT)
            call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_preprocess :: preprocess; ")
            call simple_mkdir(output_dir_motion_correct, errmsg="commander_preprocess :: preprocess;")
            if( l_pick )then
                output_dir_picker  = trim(DIR_PICKER)
                output_dir_extract = trim(DIR_EXTRACT)
                call simple_mkdir(output_dir_picker, errmsg="commander_preprocess :: preprocess; ")
                call simple_mkdir(output_dir_extract, errmsg="commander_preprocess :: preprocess;")
            endif
        endif
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! range
        if( params%stream.eq.'yes' )then
            ! STREAMING MODE
            fromto(:) = 1
        else
            ! DISTRIBUTED MODE
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_preprocess')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! numlen
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            params%numlen = len(int2str(nmovies))
        endif
        frame_counter = 0
        ! loop over exposures (movies)
        do imovie = fromto(1),fromto(2)
            ! fetch movie orientation
            o_mov = spproj%os_mic%get_ori(imovie)
            ! sanity check
            if(.not.o_mov%isthere('imgkind') )cycle
            if(.not.o_mov%isthere('movie') .and. .not.o_mov%isthere('intg'))cycle
            call o_mov%getter('imgkind', imgkind)
            select case(trim(imgkind))
                case('movie')
                    ! motion_correct
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('movie', moviename)
                    if( .not.file_exists(moviename)) cycle
                    if( cline%defined('gainref') )then
                        call mciter%iterate(cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct, gainref_fname=params%gainref)
                    else
                        call mciter%iterate(cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct)
                    endif
                    moviename_forctf = mciter%get_moviename('forctf')
                case('mic')
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('intg', moviename_forctf)
                case DEFAULT
                    cycle
            end select
            ! ctf_estimate
            params_glob%hp = params%hp_ctf_estimate
            params_glob%lp = max(params%fny, params%lp_ctf_estimate)
            call cfiter%iterate(ctfvars, moviename_forctf, o_mov, output_dir_ctf_estimate, .false.)
            ! update project
            call spproj%os_mic%set_ori(imovie, o_mov)
            ! picker
            if( l_pick )then
                params_glob%lp = max(params%fny, params%lp_pick)
                moviename_intg = mciter%get_moviename('intg')
                call piter%iterate(cline,  moviename_intg, boxfile, nptcls_out, output_dir_picker)
                call o_mov%set('boxfile', trim(boxfile)   )
                call o_mov%set('nptcls',  real(nptcls_out))
                ! update project
                call spproj%os_mic%set_ori(imovie, o_mov)
                ! extract particles
                if( trim(params%stream) .eq. 'yes' )then
                    ! needs to write and re-read project at the end as extract overwrites it
                    call spproj%write_segment_inside(params%oritype)
                    cline_extract = cline
                    call cline_extract%set('dir', trim(output_dir_extract))
                    call cline_extract%set('pcontrast', params%pcontrast)
                    if( cline%defined('box_extract') )call cline_extract%set('box', real(params%box_extract))
                    call xextract%execute(cline_extract)
                    call spproj%kill
                endif
            endif
        end do
        if( params%stream .eq. 'yes' )then
            if( .not.l_pick )then
                ! because extract performs the writing otherwise
                call spproj%write_segment_inside(params%oritype)
            endif
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_preprocess' )
        call simple_end('**** SIMPLE_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess

    subroutine exec_motion_correct( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_binoris_io,          only: binwrite_oritab
        use simple_ori,                 only: ori
        use simple_motion_correct_iter, only: motion_correct_iter
        class(motion_correct_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(motion_correct_iter)     :: mciter
        type(ctfparams)               :: ctfvars
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename, imgkind, fbody
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, lfoo(3), nframes, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('No movie to process!')
        endif
        if( params%scale > 1.05 )then
            THROW_HARD('scale cannot be > 1; exec_motion_correct')
        endif
        if( params%tomo .eq. 'yes' )then
            if( .not. params%l_dose_weight )then
                write(logfhandle,*) 'tomo=yes only supported with dose weighting!'
                THROW_HARD('give total exposure time: exp_time (in seconds) and dose_rate (in e/A2/s)')
            endif
        endif
        if( cline%defined('gainref') )then
            if(.not.file_exists(params%gainref) )then
                THROW_HARD('gain reference: '//trim(params%gainref)//' not found; motion_correct')
            endif
        endif
        ! output directory & names
        output_dir = PATH_HERE
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! determine loop range & fetch movies oris object
        if( params%tomo .eq. 'no' )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; motion_correct')
            endif
        else
            ! all movies
            fromto(1) = 1
            fromto(2) = spproj%os_mic%get_noris()
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! for series of tomographic movies we need to calculate the time_per_frame
        if( params%tomo .eq. 'yes' )then
            ! get number of frames & dim from stack
            call spproj%os_mic%getter(1, 'movie', moviename)
            call find_ldim_nptcls(moviename, lfoo, nframes)
            ! calculate time_per_frame
            params%time_per_frame = params%exp_time/real(nframes*nmovies)
        endif
        ! align
        frame_counter = 0
        cnt = 0
        do imovie=fromto(1),fromto(2)
            o = spproj%os_mic%get_ori(imovie)
            if( o%isthere('imgkind').and.o%isthere('movie') )then
                cnt = cnt + 1
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'movie' )cycle
                call o%getter('movie', moviename)
                ctfvars = spproj%get_micparams(imovie)
                if( cline%defined('gainref') )then
                    call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, trim(output_dir), gainref_fname=params%gainref)
                else
                    call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, trim(output_dir))
                endif
                call spproj%os_mic%set_ori(imovie, o)
                write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_motion_correct' )
        call simple_end('**** SIMPLE_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct

    subroutine exec_gen_pspecs_and_thumbs( self, cline )
        use simple_sp_project,       only: sp_project
        use simple_binoris_io,       only: binwrite_oritab
        use simple_ori,              only: ori
        use simple_pspec_thumb_iter, only: pspec_thumb_iter
        class(gen_pspecs_and_thumbs_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(pspec_thumb_iter)        :: ptiter
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename_intg, imgkind
        integer :: nintgs, fromto(2), iintg, ntot, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nintgs = spproj%get_nintgs()
        if( nintgs == 0 )then
            THROW_HARD('No integrated movies to process!')
        endif
        ! output directory
        output_dir = PATH_HERE
        ! determine loop range & fetch movies oris object
        if( params%l_distr_exec )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; gen_pspecs_and_thumbs')
            endif
        else
            fromto = [1,nintgs]
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! align
        cnt = 0
        do iintg=fromto(1),fromto(2)
            o = spproj%os_mic%get_ori(iintg)
            if( o%isthere('imgkind').and.o%isthere('intg') )then
                cnt = cnt + 1
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', moviename_intg)
                call ptiter%iterate(o, moviename_intg, trim(output_dir))
                call spproj%os_mic%set_ori(iintg, o)
                write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the integrated movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished('simple_commander_preprocess :: exec_gen_pspecs_and_thumbs')
        call simple_end('**** SIMPLE_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs

    subroutine exec_ctf_estimate( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_binoris_io,          only: binwrite_oritab
        use simple_ori,                 only: ori
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        class(ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline  !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(ctf_estimate_iter)       :: cfiter
        type(ctfparams)               :: ctfvars
        type(ori)                     :: o
        character(len=:), allocatable :: intg_forctf, output_dir, imgkind
        integer                       :: fromto(2), imic, ntot, cnt, state
        logical                       :: l_gen_thumb
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! read in integrated movies
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        ! output directory
        output_dir = PATH_HERE
        ! parameters & loop range
        if( params%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_ctf_estimate')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! loop over exposures (movies)
        cnt = 0
        do imic = fromto(1),fromto(2)
            cnt   = cnt + 1
            o     = spproj%os_mic%get_ori(imic)
            state = 1
            if( o%isthere('state') ) state = nint(o%get('state'))
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                if( o%isthere('forctf') )then
                    call o%getter('forctf', intg_forctf)
                else if( o%isthere('intg') )then
                    call o%getter('intg', intg_forctf)
                else
                    THROW_HARD('no image available (forctf|intg) for CTF fittings :: exec_ctf_estimate')
                endif
                l_gen_thumb = .not. o%isthere('thumb')
                ctfvars     = o%get_ctfvars()
                call cfiter%iterate( ctfvars, intg_forctf, o, trim(output_dir), l_gen_thumb)
                call spproj%os_mic%set_ori(imic, o)
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_ctf_estimate' )
        call simple_end('**** SIMPLE_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate

    subroutine exec_map_cavgs_selection( self, cline )
        use simple_image,               only: image
        use simple_corrmat,             only: calc_cartesian_corrmat
        class(map_cavgs_selection_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline !< command line input
        type(parameters)         :: params
        type(builder)            :: build
        type(image), allocatable :: imgs_sel(:), imgs_all(:)
        integer,     allocatable :: states(:)
        real,        allocatable :: correlations(:,:)
        integer :: iimg, isel, nall, nsel, loc(1), lfoo(3)
        call build%init_params_and_build_spproj(cline,params)
        ! find number of selected cavgs
        call find_ldim_nptcls(params%stk2, lfoo, nsel)
        ! find number of original cavgs
        call find_ldim_nptcls(params%stk, lfoo, nall)
        ! read images
        allocate(imgs_sel(nsel), imgs_all(nall))
        do isel=1,nsel
            call imgs_sel(isel)%new([params%box,params%box,1], params%smpd)
            call imgs_sel(isel)%read(params%stk2, isel)
        end do
        do iimg=1,nall
            call imgs_all(iimg)%new([params%box,params%box,1], params%smpd)
            call imgs_all(iimg)%read(params%stk, iimg)
        end do
        write(logfhandle,'(a)') '>>> CALCULATING CORRELATIONS'
        call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
        ! create the states array for mapping the selection
        allocate(states(nall))
        states = 0
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            states(loc(1)) = 1
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_selection

    subroutine exec_pick( self, cline )
        use simple_binoris_io,     only: binwrite_oritab
        use simple_ori,            only: ori
        use simple_picker_iter,    only: picker_iter
        use simple_image,          only: image
        class(pick_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(picker_iter)             :: piter
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, intg_name, imgkind
        character(len=LONGSTRLEN)     :: boxfile
        integer :: fromto(2), imic, ntot, nptcls_out, cnt, state
        call cline%set('oritype', 'mic')
        call params%new(cline)
        ! output directory
        output_dir = PATH_HERE
        ! parameters & loop range
        if( params%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_pick')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! read project file
        call spproj%read(params%projfile)
        ! look for movies
        if( spproj%get_nintgs() == 0 )then
            THROW_HARD('No integrated micrograph to process!')
        endif
        ! perform picking
        cnt = 0
        do imic=fromto(1),fromto(2)
            cnt   = cnt + 1
            o     = spproj%os_mic%get_ori(imic)
            state = 1
            if( o%isthere('state') ) state = nint(o%get('state'))
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', intg_name)
                call piter%iterate(cline, intg_name, boxfile, nptcls_out, output_dir)
                call spproj%os_mic%set(imic, 'boxfile', trim(boxfile))
                call spproj%os_mic%set(imic, 'nptcls', real(nptcls_out))
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_pick' )
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick

    subroutine exec_pick_chiara( self, cline)
        use simple_picker_chiara
        use simple_micops
        use simple_image
        use simple_stackops
        use simple_segmentation, only: sobel, canny, automatic_thresh_sobel
        class(pick_commander_chiara), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline !< command line input
        character(len=:), allocatable :: output_dir
        type(parameters)   :: params
        type(image)        :: mic_shrunken, mic_bin, mic_copy
        type(image)        :: imgcc, imgwi
        integer            :: ldim_shrunken(3), box_shrunken, winsz, min_sz, max_sz
        real               :: part_radius
        real,    parameter :: SHRINK = 4.
        real               :: smpd_shrunken, lp, ave, sdev, maxv, minv, thresh(1)
        if( .not. cline%defined('fname') )then
            THROW_HARD('ERROR! fname needs to be present; exec_pick_chiara')
        endif
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_pick_chiara')
        endif
        if( .not. cline%defined('part_radius') )then
            THROW_HARD('ERROR! part_radius to be present; exec_pick_chiara')
        endif
        call params%new(cline)
        if( cline%defined('thres') .and. params%detector .ne. 'sobel')then
            THROW_HARD('ERROR! thres is compatible only with sobel detector; exec_pick_chiara')
        endif
        if( cline%defined('part_concentration') .and. params%detector .ne. 'sobel')then
            THROW_HARD('ERROR! part_concentration is compatible only with sobel detector; exec_pick_chiara')
        endif
        if( cline%defined('thres') .and. cline%defined('part_concentration'))then
            THROW_HARD('ERROR! thres and part_concentration are not compatible; exec_pick_chiara')
        endif
        ! output directory
        output_dir = PATH_HERE
        ! 0) Reading and saving original micrograph
        call read_micrograph(micfname = params%fname, smpd = params%smpd)
        ! 1) Shrink and high pass filtering
        part_radius = params%part_radius
        call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
        call set_box(int(SHRINK*( 4*part_radius+10 )), box_shrunken)
        call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
        call mic_shrunken%read('shrunken_hpassfiltered.mrc')
        call mic_copy%new(ldim_shrunken, smpd_shrunken) !work on a copy not to modify the original mic
        call mic_copy%read('shrunken_hpassfiltered.mrc')
        ! 2) Low pass filtering
        if( cline%defined('lp')) then
            lp = params%lp
        else
            lp = 35.
        endif
        call mic_copy%bp(0.,lp)
        call mic_copy%ifft()
        ! 3) Edge Detection
        call mic_copy%stats( ave, sdev, maxv, minv )
        if(params%detector .eq. 'sobel') then
            if(cline%defined('thres')) then
                thresh(1) = params%thres
                write(logfhandle,*) 'threshold selected = ', thresh
                call sobel(mic_copy,thresh)
            elseif(cline%defined('part_concentration')) then
                call automatic_thresh_sobel(mic_copy, params%part_concentration)
            else
                thresh(1) = ave+.8*sdev
                call sobel(mic_copy,thresh)
            endif
        else if (params%detector .eq. 'bin') then
            call mic_copy%bin(ave+.8*sdev)
        else if (params%detector .eq. 'canny') then
            call canny(mic_copy)
        endif
        call mic_copy%write('Bin1.mrc')
        if( cline%defined('winsz')) then
            winsz = nint(params%winsz)
        else
            winsz = 4
        endif
        call mic_copy%real_space_filter(winsz,'median') !median filtering allows easy calculation of cc
        call mic_copy%write('Bin1Median.mrc')
        ! 5) Connected components (cc) identification
        call imgcc%new(ldim_shrunken, smpd_shrunken)
        call mic_copy%find_connected_comps(imgcc)
        call imgcc%write('ConnectedComponents.mrc')
        ! 6) cc filtering
        min_sz = 10*int(part_radius)  !CHOSE A PERCENTAGE OF THE of the size of the particle
        max_sz = 70*int(part_radius)
        call imgcc%elim_cc([min_sz,max_sz])
        call imgcc%write('ConnectedComponentsElimin.mrc')
        call extract_particles(mic_shrunken, imgcc, int(part_radius))
        open(unit = 17, file = "PickerInfo.txt")
        write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING>>>>>>>>>>>>>>>>>>'
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = '(a,a)') 'Working Directory ', output_dir
        write(unit = 17, fmt = "(a,f0.0)")  'Mic Shrunken by factor ', SHRINK
        write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Dim after shrink ', ldim_shrunken
        write(unit = 17, fmt = "(a,f0.0)")  'Smpd after shrink ', smpd_shrunken
        write(unit = 17, fmt = "(a,i0)")  'Hp box ', int(SHRINK*( 4*part_radius+10 ))
        write(unit = 17, fmt = "(a,f0.0)")  'Lp  ', lp
        write(unit = 17, fmt = "(a,i0,tr1, i0)")  'Connected Components size filtering  ', min_sz, max_sz
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  'SELECTED PARAMETERS '
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,tr1,f0.0)")  'smpd ', params%smpd
        write(unit = 17, fmt = "(a,tr1,f0.0)")  'part_radius  ', params%part_radius
        write(unit = 17, fmt = "(a,a)")  'detector  ', params%detector
        close(17, status = "keep")
        ! end gracefully
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick_chiara

    !> for extracting particle images from integrated DDD movies
    subroutine exec_extract( self, cline )
        use simple_image, only: image
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        use simple_ctf,   only: ctf
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(builder)                           :: build
        type(parameters)                        :: params
        type(sp_project)                        :: spproj_in, spproj
        type(nrtxtfile)                         :: boxfile
        type(image)                             :: micrograph
        type(ori)                               :: o_mic
        type(ctf)                               :: tfun
        type(ctfparams)                         :: ctfparms
        character(len=:),           allocatable :: output_dir, mic_name, imgkind
        real,                       allocatable :: boxdata(:,:)
        logical,                    allocatable :: oris_mask(:), mics_mask(:)
        character(len=LONGSTRLEN) :: stack, boxfile_name, box_fname
        integer                   :: nframes, imic, iptcl, ldim(3), nptcls,nmics,nmics_here,box, box_first, fromto(2)
        integer                   :: cnt, nmics_tot, lfoo(3), ifoo, noutside, nptcls_eff, state, iptcl_glob
        real                      :: particle_position(2)
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'no')
        call params%new(cline)
        ! init
        output_dir = PATH_HERE
        fromto(:)  = [params%fromp, params%top]
        nmics_here = fromto(2)-fromto(1)+1
        if( params%stream.eq.'yes' )then
            output_dir = DIR_EXTRACT
            fromto(:)  = [1,1]
            nmics_here = 1
            ! read in integrated movies, output project = input project
            call spproj%read(params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj%get_nintgs() /= 1 ) THROW_HARD('Incompatible # of integrated micrograph to process!')
        else
            ! read in integrated movies
            call spproj_in%read_segment(params%oritype, params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
            ! init output project
            call spproj%read_non_data_segments(params%projfile)
            call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile', params%outfile)
            params%projfile = trim(params%outfile) ! for builder later
            call spproj%os_mic%new(nmics_here)
            cnt = 0
            do imic = fromto(1),fromto(2)
                cnt = cnt + 1
                call spproj%os_mic%set_ori(cnt, spproj_in%os_mic%get_ori(imic))
            enddo
            call spproj_in%kill
        endif
        ! input boxes
        if( cline%defined('dir_box') )then
            if( .not.file_exists(params%dir_box) )then
                write(logfhandle,*)'Directory does not exist: ', trim(params%dir_box), 'simple_commander_preprocess::exec_extract'
                THROW_HARD('box directory does not exist; exec_extract')
            endif
        endif
        ! sanity checks
        allocate(mics_mask(1:nmics_here), source=.false.)
        nmics  = 0
        nptcls = 0
        do imic = 1,nmics_here
            o_mic = spproj%os_mic%get_ori(imic)
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
                box_fname = trim(params%dir_box)//'/'//fname_new_ext(basename(mic_name),'box')
                call make_relativepath(CWD_GLOB,trim(box_fname),boxfile_name)
                if( .not.file_exists(boxfile_name) )cycle
                call spproj%os_mic%set(imic, 'boxfile', trim(boxfile_name))
            else
                boxfile_name = trim(o_mic%get_static('boxfile'))
                if( .not.file_exists(boxfile_name) )cycle
            endif
            ! get number of frames from stack
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract')
            ! update mask
            mics_mask(imic) = .true.
            nmics = nmics + 1
            ! image & box dimensions
            if( nmics == 1 )call find_ldim_nptcls(mic_name, ldim, ifoo)
            if( nptcls == 0 .and. .not.cline%defined('box') )then
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 )cycle
                allocate( boxdata(nptcls,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                call boxfile%readNextDataLine(boxdata(1,:))
                call boxfile%kill
                params%box = nint(boxdata(1,3))
            endif
        enddo
        call spproj%write
        call spproj%kill
        ! actual extraction
        if( nmics == 0 )then
            ! nothing to do!
        else
            if( params%box == 0 )THROW_HARD('box cannot be zero; exec_extract')
            call build%build_general_tbox(params, cline, do3d=.false.)
            call micrograph%new([ldim(1),ldim(2),1], params%smpd)
            noutside  = 0
            box_first = 0
            ! main loop
            iptcl_glob = 0 ! extracted particle index among ALL stacks
            do imic = 1,nmics_here
                if( .not.mics_mask(imic) )cycle
                ! fetch micrograph
                o_mic = build%spproj_field%get_ori(imic)
                call o_mic%getter('imgkind', imgkind)
                call o_mic%getter('intg', mic_name)
                boxfile_name = trim(o_mic%get_static('boxfile'))
                ! box file
                nptcls = 0
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 ) cycle
                call progress(imic,nmics_tot)
                ! box checks
                if(allocated(oris_mask))deallocate(oris_mask)
                allocate(oris_mask(nptcls), source=.false., stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk("In exec_extract oris_mask")
                ! read box data & update mask
                if(allocated(boxdata))deallocate(boxdata)
                allocate( boxdata(nptcls,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk('In: exec_extract; boxdata etc., 2',alloc_stat)
                do iptcl=1,nptcls
                    call boxfile%readNextDataLine(boxdata(iptcl,:))
                    box = nint(boxdata(iptcl,3))
                    if( nint(boxdata(iptcl,3)) /= nint(boxdata(iptcl,4)) )then
                        THROW_HARD('only square windows allowed; exec_extract')
                    endif
                    ! modify coordinates if change in box (shift by half the difference)
                    if( box /= params%box ) boxdata(iptcl,1:2) = boxdata(iptcl,1:2) - real(params%box-box)/2.
                    if( .not.cline%defined('box') .and. nint(boxdata(iptcl,3)) /= params%box )then
                        write(logfhandle,*) 'box_current: ', nint(boxdata(iptcl,3)), 'box in params: ', params%box
                        THROW_HARD('inconsistent box sizes in box files; exec_extract')
                    endif
                    ! update particle mask & movie index
                    if( box_inside(ldim, nint(boxdata(iptcl,1:2)), params%box) )oris_mask(iptcl) = .true.
                end do
                nptcls_eff = count(oris_mask)
                ! extract ctf info
                ctfparms      = o_mic%get_ctfvars()
                ctfparms%smpd = params%smpd
                if( o_mic%isthere('dfx') )then
                    if( .not.o_mic%isthere('cs') .or. .not.o_mic%isthere('kv') .or. .not.o_mic%isthere('fraca') )then
                        THROW_HARD('input lacks at least cs, kv or fraca; exec_extract')
                    endif
                    ! clean micrograph stats before transfer to particles
                    call o_mic%delete_entry('xdim')
                    call o_mic%delete_entry('ydim')
                    call o_mic%delete_entry('nframes')
                    call o_mic%delete_entry('nptcls')
                    call o_mic%delete_entry('ctf_estimatecc')
                    call o_mic%delete_entry('ctfscore')
                    call o_mic%delete_entry('dferr')
                endif
                ! output stack
                stack = trim(output_dir)//trim(EXTRACT_STK_FBODY)//trim(basename(mic_name))
                ! extract windows from integrated movie
                call micrograph%read(mic_name, 1)
                ! phase-flip micrograph
                if( cline%defined('ctf') )then
                    if( trim(params%ctf).eq.'flip' .and. o_mic%isthere('dfx') )then
                        ! phase flip micrograph
                        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                        call micrograph%zero_edgeavg
                        call micrograph%fft
                        call tfun%apply_serial(micrograph, 'flip', ctfparms)
                        ! update stack ctf flag, mic flag unchanged
                        ctfparms%ctfflag = CTFFLAG_FLIP
                    endif
                endif
                ! filter out frequencies lower than the box can express to avoid aliasing
                call micrograph%bp(real(params%box) * params%smpd, 0.)
                call micrograph%ifft ! need to be here in case it was flipped
                ! write new stack
                cnt = 0
                do iptcl=1,nptcls ! loop over boxes
                    if( oris_mask(iptcl) )then
                        cnt = cnt + 1
                        ! extract the window
                        particle_position = boxdata(iptcl,1:2)
                        call micrograph%window(nint(particle_position), params%box, build%img, noutside)
                        if( params%pcontrast .eq. 'black' ) call build%img%neg()
                        call build%img%noise_norm(build%lmsk)
                        call build%img%write(trim(adjustl(stack)), cnt)
                    endif
                end do
                ! IMPORT INTO PROJECT
                call build%spproj%add_stk(trim(adjustl(stack)), ctfparms)
                ! add box coordinates to ptcl2D field only
                do iptcl=1,nptcls
                    if( .not.oris_mask(iptcl) )cycle
                    iptcl_glob = iptcl_glob + 1
                    particle_position = boxdata(iptcl,1:2)
                    call build%spproj%os_ptcl2D%set(iptcl_glob,'xpos',boxdata(iptcl,1))
                    call build%spproj%os_ptcl2D%set(iptcl_glob,'ypos',boxdata(iptcl,2))
                end do
                ! clean
                call boxfile%kill()
            enddo
            ! write
            call build%spproj%write
        endif
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_extract' )
        call simple_end('**** SIMPLE_EXTRACT NORMAL STOP ****')

        contains

            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                if( params%outside .eq. 'yes' )then
                    inside = .true.
                    return
                endif
                fromc  = coord+1       ! compensate for the c-range that starts at 0
                toc    = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.        ! box is inside
                if( any(fromc < 1) .or. toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function box_inside

    end subroutine exec_extract

    !> for extracting particle images from integrated DDD movies
    subroutine exec_reextract( self, cline )
        use simple_image, only: image
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        use simple_ctf,   only: ctf
        class(reextract_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(image)                   :: micrograph, img, mskimg
        type(ori)                     :: o_mic, o_stk
        type(ctf)                     :: tfun
        type(ctfparams)               :: ctfparms
        character(len=:), allocatable :: mic_name, imgkind
        logical,          allocatable :: ptcl_msk(:,:,:)
        integer,          allocatable :: micstk_inds(:)
        character(len=LONGSTRLEN)     :: stack, rel_stack
        integer  :: nframes,imic,iptcl,istk,nstks,nmics,prev_box,box_foo,cnt,ntot,nptcls
        integer :: prev_pos(2),new_pos(2),center(2),ishift(2),ldim(3),ldim_foo(3),noutside,fromp,top,noutside_tot
        real    :: prev_shift(2), shift2d(2), shift3d(2)
        logical :: l_2d
        call params%new(cline)
        if(cline%defined('oritype'))then
            ! whether to use shifts from 2D or 3D
            l_2d = params%oritype=='ptcl2D'
        else
            l_2d= .true.
        endif
        ! read in integrated movies
        call spproj%read(params_glob%projfile)
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        ntot = spproj%os_mic%get_noris()
        allocate(micstk_inds(ntot), source=0)
        ! sanity checks & dimensions
        box_foo  = 0
        prev_box = 0
        ldim_foo = 0
        ldim     = 0
        cnt      = 0
        do imic = 1,ntot
            ! sanity checks
            o_mic = spproj%os_mic%get_ori(imic)
            if( o_mic%isthere('state') )then
                if( o_mic%get_state() == 0 )cycle
            endif
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            call find_ldim_nptcls(mic_name, ldim_foo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_reextract')
            if( ldim_foo(1).ne.nint(o_mic%get('xdim')) )THROW_HARD('Inconsistent x dimensions; exec_reextract')
            if( ldim_foo(2).ne.nint(o_mic%get('ydim')) )THROW_HARD('Inconsistent y dimensions; exec_reextract')
            if( any(ldim == 0) ) ldim = ldim_foo
            if( .not.o_mic%isthere('nptcls') ) THROW_HARD('micrographs not imported correctly 1; exec_reextract')
            o_stk  = spproj%os_stk%get_ori(cnt+1)
            fromp = nint(o_stk%get('fromp'))
            top   = nint(o_stk%get('top'))
            do iptcl=fromp,top
                if(.not.spproj%os_ptcl2D%isthere(iptcl,'xpos') ) THROW_HARD('missing particle coordinates 1; exec_reextract')
                if(.not.spproj%os_ptcl2D%isthere(iptcl,'ypos') ) THROW_HARD('missing particle coordinates 2; exec_reextract')
            enddo
            box_foo = nint(o_stk%get('box'))
            if( prev_box == 0 ) prev_box = box_foo
            if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_reextract')
            ! update mask
            cnt = cnt+1
            micstk_inds(cnt) = imic ! index to os_mic (including movies)
        enddo
        nmics = count(micstk_inds>0)
        if( nmics == 0 )  THROW_HARD('No particles to extract! exec_extract')
        if( .not.cline%defined('box') ) params%box = prev_box
        if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_extract')
        ! extraction
        call mskimg%disc([params%box,params%box,1], params%smpd, params%msk, ptcl_msk)
        call mskimg%kill
        call micrograph%new([ldim(1),ldim(2),1], params%smpd)
        noutside_tot = 0
        nptcls   = 0
        nstks    = spproj%os_stk%get_noris()
        do istk = 1,nstks
            call progress(istk,nstks)
            imic = micstk_inds(istk)
            if( imic == 0 )cycle
            ! fetch & read micrograph
            o_mic = spproj%os_mic%get_ori(imic)
            o_stk = spproj%os_stk%get_ori(istk)
            call o_mic%getter('intg', mic_name)
            ctfparms = o_mic%get_ctfvars()
            fromp = nint(o_stk%get('fromp'))
            top   = nint(o_stk%get('top'))
            call micrograph%read(mic_name)
            ! phase-flip micrograph. To-add as an option to not flip?
            if( ctfparms%ctfflag == CTFFLAG_FLIP )then
                if( o_mic%isthere('dfx') )then
                    ! phase flip micrograph
                    tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                    call micrograph%zero_edgeavg
                    call micrograph%fft
                    call tfun%apply_serial(micrograph, 'flip', ctfparms)
                endif
            endif
            call micrograph%bp(real(params%box) * params%smpd, 0.)
            call micrograph%ifft ! need to be here in case it was flipped
            stack  = trim(EXTRACT_STK_FBODY)//trim(basename(mic_name))
            cnt    = 0
            do iptcl=fromp,top
                cnt = cnt + 1
                prev_pos(1) = nint(spproj%os_ptcl2D%get(iptcl,'xpos'))
                prev_pos(2) = nint(spproj%os_ptcl2D%get(iptcl,'ypos'))
                if( l_2d )then
                    prev_shift = spproj%os_ptcl2D%get_2Dshift(iptcl)
                else
                    prev_shift = spproj%os_ptcl3D%get_2Dshift(iptcl)
                endif
                ! calc new position & shift
                center  = prev_pos + prev_box/2 - 1
                ishift  = nint(prev_shift)
                center  = center - ishift
                new_pos = center - params%box/2 + 1
                if( .not.box_inside(ldim, new_pos, params%box) )then
                    if( params%outside.eq.'yes')then
                        ! all done, border guaranteed
                    else
                        ! if out of mic keep previous position with new box
                        ! at the risk of including a border
                        center  = prev_pos + prev_box/2 - 1
                        new_pos = center - params%box/2 + 1
                    endif
                else
                    ! all done
                endif
                ! updates doc
                call spproj%os_ptcl2D%set(iptcl,'xpos',real(new_pos(1)))
                call spproj%os_ptcl2D%set(iptcl,'ypos',real(new_pos(2)))
                call spproj%os_stk%set(istk,'box',     real(params%box))
                if( l_2d )then
                    shift2d = prev_shift - real(ishift)
                    shift3d = spproj%os_ptcl3D%get_2Dshift(iptcl) - real(ishift)
                else
                    shift2d = spproj%os_ptcl2D%get_2Dshift(iptcl) - real(ishift)
                    shift3d = prev_shift - real(ishift)
                endif
                call spproj%os_ptcl2D%set_shift(iptcl,shift2d)
                call spproj%os_ptcl3D%set_shift(iptcl,shift3d)
                ! image extraction
                noutside = 0
                call micrograph%window(new_pos, params%box, img, noutside)
                if( noutside > 0 )then
                    noutside_tot = noutside_tot+1
                    if( params%outside.ne.'yes')then
                        call spproj%os_ptcl2D%set(iptcl,'state',0.)
                        call spproj%os_ptcl3D%set(iptcl,'state',0.)
                    endif
                endif
                if( params%pcontrast .eq. 'black' ) call img%neg()
                call img%noise_norm(ptcl_msk)
                call img%write(trim(adjustl(stack)), cnt)
                nptcls = nptcls+1
            enddo
            call make_relativepath(CWD_GLOB, stack, rel_stack)
            call spproj%os_stk%set(istk,'stk',rel_stack)
        enddo
        call progress(1,1)
        write(logfhandle,'(A,I8)')'>>> RE-EXTRACTED  PARTICLES: ', nptcls
        write(logfhandle,'(A,I8)')'>>> OUT OF LIMITS PARTICLES: ', noutside_tot
        call spproj%write
        call simple_end('**** SIMPLE_REEXTRACT NORMAL STOP ****')

        contains

            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                if( params%outside .eq. 'yes' )then
                    inside = .true.
                    return
                endif
                fromc  = coord+1       ! compensate for the c-range that starts at 0
                toc    = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.        ! box is inside
                if( any(fromc < 1) .or. toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function box_inside

    end subroutine exec_reextract

    subroutine exec_pick_extract( self, cline )
        use simple_ori,                 only: ori
        use simple_sp_project,          only: sp_project
        use simple_picker_iter,         only: picker_iter
        use simple_binoris_io,          only: binwrite_oritab
        class(pick_extract_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)              :: params
        type(ori)                     :: o_mic
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        character(len=:), allocatable :: micname, output_dir_picker, fbody, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        integer :: fromto(2), imic, ntot, nptcls_out, state
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        if( params%scale > 1.01 )then
            THROW_HARD('scale cannot be > 1; exec_preprocess')
        endif
        if( params%stream.ne.'yes' ) THROW_HARD('streaming only application')
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('no micrograph to process!')
        ! output directories
        if( params%stream.eq.'yes' )then
            output_dir_picker  = trim(DIR_PICKER)
            output_dir_extract = trim(DIR_EXTRACT)
            call simple_mkdir(output_dir_picker,errmsg="commander_preprocess :: preprocess; ")
            call simple_mkdir(output_dir_extract,errmsg="commander_preprocess :: preprocess; ")
        else
            output_dir_picker  = PATH_HERE
            output_dir_extract = PATH_HERE
        endif
        ! command lines
        cline_extract = cline
        call cline_extract%set('dir', trim(output_dir_extract))
        call cline_extract%set('pcontrast', params%pcontrast)
        if( cline%defined('box_extract') )call cline_extract%set('box', real(params%box_extract))
        call cline%delete('box')
        call cline_extract%delete('box_extract')
         ! file name
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! range
        if( params%stream.eq.'yes' )then
            fromto(:) = 1
        else
            fromto(:) = [params%fromp, params%top]
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! loop over exposures (movies)
        do imic = fromto(1),fromto(2)
            ! fetch movie orientation
            o_mic = spproj%os_mic%get_ori(imic)
            ! sanity check
            state = 1
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
            if( state == 0 ) cycle
            if( .not.o_mic%isthere('intg')   )cycle
            call o_mic%getter('intg', micname)
            if( .not.file_exists(micname)) cycle
            ! picker
            params_glob%lp = max(params%fny, params%lp_pick)
            call piter%iterate(cline, micname, boxfile, nptcls_out, output_dir_picker)
            call o_mic%set('boxfile', trim(boxfile)   )
            call o_mic%set('nptcls',  real(nptcls_out))
            ! update project
            call spproj%os_mic%set_ori(imic, o_mic)
            call spproj%write_segment_inside(params%oritype)
            ! extract particles
            call xextract%execute(cline_extract)
            call spproj%kill
        end do
        if( params%stream .eq. 'yes' )then
            ! nothing to do, extract did it
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_pick_extract' )
        call simple_end('**** SIMPLE_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_pick_extract

    subroutine exec_make_pickrefs( self, cline )
        use simple_oris,                only: oris
        use simple_sp_project,          only: sp_project
        use simple_sym,                 only: sym
        use simple_image,               only: image
        use simple_projector_hlev,      only: reproject
        class(make_pickrefs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(oris)                    :: os
        type(sym)                     :: pgrpsyms
        type(image)                   :: ref3D, ref2D
        type(image),      allocatable :: projs(:)
        character(len=:), allocatable :: imgkind
        character(len=:), allocatable :: volname, volcavgname, cavgname
        real,             allocatable :: states(:)
        integer, parameter :: NREFS=100, NPROJS=20
        real    :: ang, rot, smpd_here
        integer :: nrots, iref, irot, ldim(3), ldim_here(3), ifoo, ncavgs, icavg
        integer :: ncls2D, cnt, n_os_out, i, nsel
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        if( params%stream.eq.'yes' ) THROW_HARD('not a streaming application')
        if( .not. cline%defined('pgrp') ) params%pgrp = 'd1' ! only northern hemisphere
        ! point-group object
        call pgrpsyms%new(trim(params%pgrp))
        ! user input overrides volume/cavgs in project file
        if( cline%defined('vol1') )then
            volname  = trim(params%vols(1))
        else if( cline%defined('refs') )then
            cavgname = trim(params%refs)
        else
            ! read project file
            call spproj%read(params%projfile)
            n_os_out = spproj%os_out%get_noris()
            if( n_os_out == 0 )then
                THROW_HARD('Nothing in os_out. Need vol / cavgs for creating picking references!')
            endif
            ! interrogate project for vol / cavgs
            do i=1,n_os_out
                if( spproj%os_out%isthere(i,'imgkind') )then
                    call spproj%os_out%getter(i,'imgkind',imgkind)
                    select case(imgkind)
                        case( 'vol' )
                            call spproj%os_out%getter(i, 'vol',      volname)
                        case( 'vol_cavg' )
                            call spproj%os_out%getter(i, 'vol_cavg', volcavgname)
                        case( 'cavg' )
                            call spproj%os_out%getter(i, 'stk',     cavgname)
                    end select
                endif
            enddo
        endif
        if( allocated(volname) )then
            ! have 3D reference
        else if( allocated(volcavgname) )then
            ! have 3D reference
            volname = volcavgname
        else if( allocated(cavgname) )then
            ! have 2D references
        else
            THROW_HARD('Could not identify vol / cavgs in project for creating picking references!')
        endif
        if( allocated(volname) )then
            ! find logical dimension & read reference volume
            call find_ldim_nptcls(volname, ldim_here, ifoo, smpd=smpd_here)
            call ref3D%new(ldim_here, smpd_here)
            call ref3D%read(volname)
            call scale_ref(ref3D, params%smpd)
            ! make projection directions
            call os%new(NPROJS)
            call pgrpsyms%build_refspiral(os)
            ! generate reprojections
            projs = reproject(ref3D, os)
            nrots = NREFS / NPROJS
            nsel  = NPROJS
        else
            ! find logical dimension & read class averages
            print *,trim(cavgname)
            call find_ldim_nptcls(cavgname, ldim_here, ncavgs, smpd=smpd_here)
            ldim_here(3) = 1
            if( cline%defined('refs') )then
                nsel = ncavgs
                allocate(states(nsel), source=1.)
            else
                ! check consistency with cls2D field
                ncls2D = spproj%os_cls2D%get_noris()
                if( ncavgs /= ncls2D )then
                    print *, '# cavgs in file:          ', ncls2D
                    print *, '# entries in cls2D field: ', ncavgs
                    THROW_HARD('inconsistent # cavgs in file vs. project field')
                endif
                ! manage selection
                states = spproj%os_cls2D%get_all('state')
                nsel   = count(states > 0.5)
                if( nsel > NREFS / 4 )then
                    THROW_HARD('too many class averages selected as references for picking, select max 25')
                endif
            endif
            ! read selected cavgs
            allocate( projs(nsel) )
            cnt = 0
            do icavg=1,ncavgs
                if( states(icavg) > 0.5 )then
                    cnt = cnt + 1
                    call projs(cnt)%new(ldim_here, smpd_here)
                    call projs(cnt)%read(cavgname, icavg)
                    call scale_ref(projs(cnt), params%smpd)
                endif
            end do
            nrots  = NREFS / nsel
        endif
        ! expand in in-plane rotation and write to file
        if( nrots > 1 )then
            print *,ldim
            call ref2D%new([ldim(1),ldim(2),1], params%smpd)
            ang = 360./real(nrots)
            rot = 0.
            cnt = 0
            do iref=1,nsel
                do irot=1,nrots
                    cnt = cnt + 1
                    call projs(iref)%rtsq(rot, 0., 0., ref2D)
                    if(params%pcontrast .eq. 'black') call ref2D%neg
                    call ref2D%write(trim(PICKREFS)//params%ext, cnt)
                    rot = rot + ang
                end do
            end do
        else
            ! should never happen
            do iref=1,nsel
                if(params%pcontrast .eq. 'black') call projs(iref)%neg
                call projs(iref)%write(trim(PICKREFS)//params%ext, iref)
            end do
        endif
        ! end gracefully
        call simple_touch('MAKE_PICKREFS_FINISHED', errmsg='In: commander_preprocess::exec_make_pickrefs')
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
        contains

            subroutine scale_ref(refimg, smpd_target)
                class(image), intent(inout) :: refimg
                real,         intent(in)    :: smpd_target
                type(image) :: targetimg
                integer     :: ldim_ref(3), ldim_target(3)
                real        :: smpd_ref, scale
                smpd_ref = refimg%get_smpd()
                scale    = smpd_target / smpd_ref
                if( is_equal(scale,1.) )then
                    ldim = ldim_here
                    return
                endif
                ldim_ref       = refimg%get_ldim()
                ldim_target(1) = round2even(real(ldim_ref(1))/scale)
                ldim_target(2) = ldim_target(1)
                ldim_target(3) = 1
                if( refimg%is_3d() )ldim_target(3) = ldim_target(1)
                call refimg%norm
                if( scale > 1. )then
                    ! downscaling
                    call refimg%fft
                    call refimg%clip_inplace(ldim_target)
                    call refimg%ifft
                else
                    call targetimg%new(ldim_target, smpd_target)
                    call refimg%fft
                    call refimg%pad(targetimg, backgr=0.)
                    call targetimg%ifft
                    refimg = targetimg
                    call targetimg%kill
                endif
                ! updates dimensions
                ldim = ldim_target
            end subroutine

    end subroutine exec_make_pickrefs

end module simple_commander_preprocess
