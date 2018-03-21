! concrete commander: pre-processing routines
module simple_commander_preprocess
#include "simple_lib.f08"
use simple_cmdline,             only: cmdline
use simple_params,              only: params
use simple_build,               only: build
use simple_commander_base,      only: commander_base
use simple_image,               only: image
use simple_procimgfile,         only: neg_imgfile
use simple_nrtxtfile,           only: nrtxtfile
use simple_imgfile,             only: imgfile
use simple_oris,                only: oris
use simple_ori,                 only: ori
use simple_imghead,             only: find_ldim_nptcls
use simple_corrmat,             only: calc_cartesian_corrmat
use simple_motion_correct_iter, only: motion_correct_iter
use simple_ctf_estimate_iter,   only: ctf_estimate_iter
use simple_picker_iter,         only: picker_iter
use simple_qsys_funs,           only: qsys_job_finished
use simple_binoris_io           ! use all in there
implicit none

public :: preprocess_commander
public :: select_frames_commander
public :: boxconvs_commander
public :: powerspecs_commander
public :: motion_correct_commander
public :: ctf_estimate_commander
public :: select_commander
public :: make_pickrefs_commander
public :: pick_commander
public :: extract_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_commander
  contains
    procedure :: execute      => exec_preprocess
end type preprocess_commander
type, extends(commander_base) :: select_frames_commander
  contains
    procedure :: execute      => exec_select_frames
end type select_frames_commander
type, extends(commander_base) :: boxconvs_commander
  contains
    procedure :: execute      => exec_boxconvs
end type boxconvs_commander
type, extends(commander_base) :: powerspecs_commander
 contains
   procedure :: execute       => exec_powerspecs
end type powerspecs_commander
type, extends(commander_base) :: motion_correct_commander
  contains
    procedure :: execute      => exec_motion_correct
end type motion_correct_commander
type, extends(commander_base) :: ctf_estimate_commander
  contains
    procedure :: execute      => exec_ctf_estimate
end type ctf_estimate_commander
type, extends(commander_base) :: select_commander
  contains
    procedure :: execute      => exec_select
end type select_commander
type, extends(commander_base) :: make_pickrefs_commander
  contains
    procedure :: execute      => exec_make_pickrefs
end type make_pickrefs_commander
type, extends(commander_base) :: pick_commander
  contains
    procedure :: execute      => exec_pick
end type pick_commander
type, extends(commander_base) :: extract_commander
  contains
    procedure :: execute      => exec_extract
end type extract_commander

contains

    subroutine exec_preprocess( self, cline )
        class(preprocess_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(ctf_estimate_iter)    :: cfiter
        type(motion_correct_iter)  :: mciter
        type(picker_iter)          :: piter
        type(extract_commander)    :: xextract
        type(cmdline)              :: cline_extract
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:),      allocatable :: fname_unidoc_output
        character(len=:),      allocatable :: moviename_forctf, moviename_intg
        character(len=:),      allocatable :: fname_stk_extract, fname_ctf_extract
        character(len=:),      allocatable :: output_dir, output_dir_ctf_estimate
        character(len=:),      allocatable :: output_dir_motion_correct, output_dir_extract
        character(len=:),      allocatable :: output_dir_picker, output_dir_unidoc
        character(len=STDLEN) :: boxfile, dir_ptcls, movie_fbody, movie_ext, movie_fname
        type(params) :: p
        type(oris)   :: os_uni
        type(ori)    :: orientation
        real         :: smpd_original, smpd_scaled
        integer      :: nmovies, fromto(2), imovie, ntot, movie_counter
        integer      :: frame_counter, movie_ind, nptcls_out
        logical      :: l_pick
        p = params(cline, spproj_a_seg=STK_SEG) ! constants & derived constants produced
        if( p%scale > 1.05 )then
            stop 'scale cannot be > 1; simple_commander_preprocess :: exec_preprocess'
        endif
        if( p%tomo .eq. 'yes' )then
            stop 'tomography mode (tomo=yes) not yet supported!'
        endif
        if( cline%defined('refs') )then
            l_pick = .true.
        else
            l_pick = .false.
        endif
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            p%numlen = len(int2str(nmovies))
        endif
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            if( p%stream.eq.'yes' )then
                output_dir = trim(DIR_PREPROC_STREAM)
            else
                output_dir = trim(DIR_PREPROC)
            endif
        endif
        output_dir_ctf_estimate   = trim(output_dir)//trim(DIR_CTF_ESTIMATE)
        output_dir_motion_correct = trim(output_dir)//trim(DIR_MOTION_CORRECT)
        output_dir_unidoc         = trim(output_dir)//trim(DIR_UNIDOC)
        call mkdir(output_dir)
        call mkdir(output_dir_ctf_estimate)
        call mkdir(output_dir_motion_correct)
        call mkdir(output_dir_unidoc)
        if( l_pick )then
            output_dir_picker  = trim(output_dir)//trim(DIR_PICKER)
            output_dir_extract = trim(output_dir)//trim(DIR_EXTRACT)
            call mkdir(output_dir_picker)
            call mkdir(output_dir_extract)
        endif
        ! filenaming and range
        if( p%stream.eq.'yes' )then
            ! STREAMING MODE
            movie_fname = remove_abspath(trim(movienames(1)))
            movie_ext   = fname2ext(trim(movie_fname))
            movie_fbody = get_fbody(trim(movie_fname), trim(movie_ext))
            if( cline%defined('fbody') )movie_fbody = trim(p%fbody)//trim(movie_fbody)
            p%fbody = trim(movie_fbody)
            ! ctf: on call and output set below
            allocate(fname_unidoc_output, source=trim(DIR_UNIDOC)//trim(UNIDOC_OUTPUT)//&
                &trim(movie_fbody)//trim(METADATA_EXT))
            ! picker: on call
            allocate(fname_stk_extract, source=trim(EXTRACT_STK_FBODY)//trim(movie_fbody)//'.'//trim(movie_ext))
            allocate(fname_ctf_extract, source=trim(EXTRACT_PARAMS_FBODY)//trim(movie_fbody)//trim(METADATA_EXT))
            call cline%set('fbody', trim(p%fbody))
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
        else
            if( p%l_distr_exec )then
                ! DISTRIBUTED MODE
                if( cline%defined('outfile') )then
                    allocate(fname_unidoc_output, source=trim(output_dir_unidoc)//trim(p%outfile))
                else
                    allocate(fname_unidoc_output, source=trim(output_dir_unidoc)//&
                        &trim(UNIDOC_OUTPUT)//'_part'//int2str_pad(p%part,p%numlen)//trim(METADATA_EXT))
                endif
                if( cline%defined('fromp') .and. cline%defined('top') )then
                    fromto(1) = p%fromp
                    fromto(2) = p%top
                else
                    stop 'fromp & top args need to be defined in parallel execution; exec_preprocess'
                endif
            else
                ! PRIVATE
                allocate(fname_unidoc_output,source=trim(output_dir_unidoc)//trim(UNIDOC_OUTPUT)//trim(METADATA_EXT))
                fromto(1) = 1
                if( cline%defined('startit') ) fromto(1) = p%startit
                fromto(2) = nmovies
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        frame_counter = 0
        movie_counter = 0
        call orientation%new
        call os_uni%new(ntot)
        smpd_original = p%smpd
        ! loop over exposures (movies)
        do imovie=fromto(1),fromto(2)
            p%smpd    = smpd_original
            if( ntot == 1 )then
                movie_ind = p%part ! streaming mode
            else
                movie_ind = imovie ! standard mode
            endif
            ! motion_correct
            call mciter%iterate(cline, p, orientation, movie_ind, movie_counter,&
                &frame_counter, movienames(imovie), smpd_scaled, output_dir_motion_correct)
            call os_uni%set_ori(movie_counter, orientation)
            p%smpd           = smpd_scaled
            movie_counter    = movie_counter - 1
            ! ctf_estimate
            moviename_forctf = mciter%get_moviename('forctf')    !! realloc warning
            moviename_intg   = mciter%get_moviename('intg')
            p%hp             = p%hp_ctfestimate
            p%lp             = p%lp_ctf_estimate
            call cfiter%iterate(p, movie_ind, movie_counter, moviename_forctf, os_uni, output_dir_ctf_estimate)
            if( p%stream .eq. 'yes' ) call os_uni%write(fname_unidoc_output)
            ! picker
            if( l_pick )then
                movie_counter = movie_counter - 1
                p%lp          = p%lp_pick
                call piter%iterate(cline, p, movie_counter, moviename_intg, boxfile, nptcls_out, output_dir_picker)
                call os_uni%set(movie_counter, 'boxfile', trim(boxfile)   )
                call os_uni%set(movie_counter, 'nptcls',  real(nptcls_out))
                if( p%stream .eq. 'yes' ) call os_uni%write(fname_unidoc_output)
            endif
            if( p%stream .eq. 'yes' )then
                ! extract particles & params
                if( l_pick )then
                    cline_extract = cline
                    call cline_extract%set('dir',       trim(output_dir_extract))
                    call cline_extract%set('smpd',      p%smpd)
                    call cline_extract%set('unidoc',    fname_unidoc_output)
                    call cline_extract%set('outfile',   fname_ctf_extract)
                    call cline_extract%set('outstk',    fname_stk_extract)
                    call cline_extract%set('pcontrast', p%pcontrast)
                    if( cline%defined('box_extract') )then
                        call cline_extract%set('box', real(p%box_extract))
                    endif
                    call xextract%execute(cline_extract)
                endif
            endif
        end do
        if( p%stream .eq. 'no' )then
            ! write unidoc
            call os_uni%write(fname_unidoc_output)
        endif
        ! destruct
        call os_uni%kill
        deallocate(fname_unidoc_output)
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_preprocess' )
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess

    subroutine exec_select_frames( self, cline )
        class(select_frames_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nframes, frame, numlen, ldim(3)
        integer                            :: fromto(2), movie, cnt, cnt2, ntot, lfoo(3), ifoo
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:), allocatable      :: new_name
        type(image)                        :: img_frame
        p = params(cline) ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general objects built
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        ! find ldim and numlen (length of number string)
        if( cline%defined('startit') )then
            call find_ldim_nptcls(movienames(p%startit), ldim, ifoo)
        endif
        DebugPrint  'logical dimension: ', ldim
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        numlen = len(int2str(nmovies))
        DebugPrint  'length of number string: ', numlen
        ! determine loop range
        if( cline%defined('part') )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
                ntot = fromto(2)-fromto(1)+1
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_select_frames'
            endif
        else
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
            ntot      = nmovies
        endif
        DebugPrint  'fromto: ', fromto(1), fromto(2)
        call img_frame%new([ldim(1),ldim(2),1], p%smpd)
        ! loop over exposures (movies)
        cnt2 = 0
        do movie=fromto(1),fromto(2)
            if( .not. file_exists(movienames(movie)) )then
                write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
            endif
            cnt2 = cnt2+1
            ! get number of frames from stack
            call find_ldim_nptcls(movienames(movie), lfoo, nframes)
            DebugPrint  'number of frames: ', nframes
            ! create new name
            allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
            cnt = 0
            do frame=p%fromp,p%top
                cnt = cnt+1
                call img_frame%read(movienames(movie),frame)
                call img_frame%write(new_name,cnt)
            end do
            deallocate(new_name)
            write(*,'(f4.0,1x,a)') 100.*(real(cnt2)/real(ntot)), 'percent of the movies processed'
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT_FRAMES NORMAL STOP ****')
    end subroutine exec_select_frames

    subroutine exec_boxconvs( self, cline )
        class(boxconvs_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: tmp
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline) ! parameters generated
        call b%build_general_tbox( p, cline, do3d=.false. ) ! general objects built
        ! do the work
        if( cline%defined('stk') )then
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg)
                tmp = b%img%boxconv(p%boxconvsz)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            DebugPrint  'read the img filenames'
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
            DebugPrint  'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            iimg_start = 1
            if( cline%defined('startit') ) iimg_start = p%startit
            iimg_stop  = nimgs
            DebugPrint  'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1)
                tmp = b%img%boxconv(p%boxconvsz)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, nimgs)
            end do
            call tmp%kill
            deallocate(imgnames)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_BOXCONVS NORMAL STOP ****')
    end subroutine exec_boxconvs

    subroutine exec_powerspecs( self, cline )
        class(powerspecs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: powspec, tmp, mask
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline, spproj_a_seg=STK_SEG)         ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general toolbox built
        ! create mask
        call tmp%new([p%clip,p%clip,1], p%smpd)
        tmp = cmplx(1.,0.)
        call tmp%bp(0.,p%lp,0.)
        call tmp%ft2img('real', mask)
        if( p%l_distr_exec )then
            if( p%part == 1 ) call mask%write('resolution_mask.mrc', 1)
        else
            call mask%write('resolution_mask.mrc', 1)
        endif
        ! do the work
        if( cline%defined('stk') )then
            if( p%l_distr_exec )then
                stop 'stk input incompatible with distributed exection; commander_preprocess :: powerspecs'
            endif
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg)
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            DebugPrint  'read the img filenames'
            if( cline%defined('numlen') )then
                ! nothing to do
            else
                p%numlen = len(int2str(nimgs))
            endif
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
            DebugPrint  'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            if( p%l_distr_exec )then
                iimg_start = p%fromp
                iimg_stop  = p%top
            else
                iimg_start = 1
                if( cline%defined('startit') ) iimg_start = p%startit
                iimg_stop  = nimgs
            endif
            DebugPrint  'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1)
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                if( p%l_distr_exec )then
                    call tmp%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(iimg,p%numlen)//p%ext)
                    call progress(iimg, nimgs)
                else
                    call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                    call progress(iimg, nimgs)
                endif
            end do
        endif
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_powerspecs' )
        ! end gracefully
        call simple_end('**** SIMPLE_POWERSPECS NORMAL STOP ****')
    end subroutine exec_powerspecs

    subroutine exec_motion_correct( self, cline )
        use simple_sp_project, only: sp_project
        class(motion_correct_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(params)                  :: p
        type(motion_correct_iter)     :: mciter
        type(sp_project)              :: spproj
        type(ori)                     :: orientation, o_mov
        character(len=:), allocatable :: output_dir, moviename
        real    :: smpd_scaled
        integer :: nmovies, fromto(2), imovie, ntot, movie_counter, frame_counter, lfoo(3), nframes
        p = params(cline, spproj_a_seg=MIC_SEG) ! constants & derived constants produced
        if( p%scale > 1.05 )then
            stop 'scale cannot be > 1; simple_commander_preprocess :: exec_motion_correct'
        endif
        if( p%tomo .eq. 'yes' )then
            if( .not. p%l_dose_weight )then
                write(*,*) 'tomo=yes only supported with dose weighting!'
                stop 'give total exposure time: exp_time (in seconds) and dose_rate (in e/A2/s)'
            endif
        endif
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_MOTION_CORRECT)
        endif
        call mkdir(output_dir)
        ! determine loop range & fetch movies oris object
        if( p%tomo .eq. 'no' )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_motion_correct'
            endif
            call spproj%read_segment( 'mic', p%projfile, fromto )
        else
            ! all movies
            call spproj%read_segment( 'mic', p%projfile)
            fromto(1) = 1
            fromto(2) = spproj%os_mic%get_noris()
            !!!!!!!!!!!!!!!!!!!!!!!!!!!! should check for 'movie' field here? !!!!!!!!!!!!!!!!!!1
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! for series of tomographic movies we need to calculate the time_per_frame
        if( p%tomo .eq. 'yes' )then
            ! get number of frames & dim from stack
            !!!!!!!!!!!!!!!!!!!!!!!!!!!! should check for 'movie' field here? !!!!!!!!!!!!!!!!!!1
            call spproj%os_mic%getter(1, 'movie', moviename)
            call find_ldim_nptcls(moviename, lfoo, nframes)
            ! calculate time_per_frame
            p%time_per_frame = p%exp_time/real(nframes*nmovies)
        endif
        ! align
        frame_counter = 0
        movie_counter = 0
        do imovie=fromto(1),fromto(2)
            call orientation%new_ori_clean
            ! movie name
            call spproj%os_mic%getter(imovie, 'movie', moviename)
            ! parameters transfer
            if( spproj%os_mic%isthere(imovie, 'cs') )    call orientation%set('cs', spproj%os_mic%get(imovie, 'cs'))
            if( spproj%os_mic%isthere(imovie, 'kv') )    call orientation%set('kv', spproj%os_mic%get(imovie, 'kv'))
            if( spproj%os_mic%isthere(imovie, 'fraca') ) call orientation%set('fraca', spproj%os_mic%get(imovie, 'fraca'))
            ! actual alignment
            call mciter%iterate(cline, p, orientation, imovie, movie_counter,&
            &frame_counter, moviename, smpd_scaled, trim(output_dir))
            call spproj%os_mic%set_ori(imovie, orientation)
            call orientation%print_ori()
            write(*,'(f4.0,1x,a)') 100.*(real(movie_counter)/real(ntot)), 'percent of the movies processed'
        end do
        ! output
        call binwrite_oritab(p%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_motion_correct' )
        call simple_end('**** SIMPLE_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct

    subroutine exec_ctf_estimate( self, cline )
        use simple_sp_project, only: sp_project
        class(ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline  !< command line input
        type(params)                       :: p
        type(ctf_estimate_iter)            :: cfiter
        type(sp_project)                   :: spproj
        character(len=:),      allocatable :: intg_forctf, output_dir
        character(len=STDLEN) :: intg_fbody, intg_ext, intg_name
        type(oris)            :: os
        integer               :: nmovies, fromto(2), imic, ntot, intg_counter
        p = params(cline, spproj_a_seg=MIC_SEG) ! constants & derived constants produced
        !call read_filetable(p%filetab, movienames_forctf)
        !nmovies = size(movienames_forctf)
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_CTF_ESTIMATE)
        endif
        call mkdir(output_dir)
        ! parameters & loop range
        if( p%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:)   = 1 !!!!!!!!!!!!!!!!!!!!!!! NOW SHOULD POINT TO MICROGRAPH OF ORIGIN ????
            call spproj%os_mic%getter(1, 'intg', intg_forctf)
            intg_name  = remove_abspath(trim(intg_forctf))
            intg_ext   = fname2ext(trim(intg_name))
            intg_fbody = get_fbody(trim(intg_name), trim(intg_ext))
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_ctf_estimate'
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        call os%new(ntot)
        ! read in integrated movies
        call spproj%read_segment( 'mic', p%projfile, fromto ) !!!!!!!!!!!!!!!!!!!!!!! SO THIS MAKES SENSE
        ! loop over exposures (movies)
        intg_counter = 0
        do imic = fromto(1),fromto(2)
            call spproj%os_mic%getter(imic, 'intg', intg_forctf)
            call cfiter%iterate(p, imic, intg_counter, intg_forctf, os, trim(output_dir))
            write(*,'(f4.0,1x,a)') 100.*(real(intg_counter)/real(ntot)), 'percent of the micrographs processed'
        end do
        !call os%write(fname_ctf_estimate_output)
        call os%kill
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_ctf_estimate' )
        ! end gracefully
        call simple_end('**** SIMPLE_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate

    subroutine exec_select( self, cline )
        class(select_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline !< command line input
        type(params)                       :: p
        type(image)                        :: stk3_img
        type(image),           allocatable :: imgs_sel(:), imgs_all(:)
        character(len=STDLEN), allocatable :: imgnames(:)
        integer,               allocatable :: selected(:)
        real,                  allocatable :: correlations(:,:)
        logical,               allocatable :: lselected(:)
        character(len=STDLEN)              :: cmd_str
        integer                            :: iimg, isel, nall, nsel, loc(1), ios
        integer                            :: funit, ldim(3), ifoo, lfoo(3), io_stat
        ! error check
        if( cline%defined('stk3') .or. cline%defined('filetab') )then
            ! all good
        else
            stop 'Need either stk3 or filetab are part of the command line!'
        endif
        p = params(cline) ! parameters generated
        ! find number of selected cavgs
        call find_ldim_nptcls(p%stk2, lfoo, nsel)
        ! find number of original cavgs
        call find_ldim_nptcls(p%stk, lfoo, nall)
        ! read images
        allocate(imgs_sel(nsel), imgs_all(nall))
        do isel=1,nsel
            call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
            call imgs_sel(isel)%read(p%stk2, isel)
        end do
        do iimg=1,nall
            call imgs_all(iimg)%new([p%box,p%box,1], p%smpd)
            call imgs_all(iimg)%read(p%stk, iimg)
        end do
        if( file_exists('corrmat_select.bin') )then
            allocate(correlations(nsel,nall), stat=alloc_stat)
            allocchk('In: exec_select; simple_commander_preprocess')
            ! read matrix
            call fopen(funit, status='OLD', action='READ', file='corrmat_select.bin', access='STREAM', iostat=io_stat)
            call fileio_errmsg('simple_commander_preprocess ; fopen error when opening corrmat_select.bin  ', io_stat)
            read(unit=funit,pos=1,iostat=io_stat) correlations
            ! Check if the read was successful
            if(io_stat/=0) then
                call fileio_errmsg('**ERROR(simple_commander_preprocess): I/O error reading corrmat_select.bin. Remove the file to override the memoization.', io_stat)
            endif

            call fclose(funit,errmsg='simple_commander_preprocess ; error when closing corrmat_select.bin  ')
        else
            write(*,'(a)') '>>> CALCULATING CORRELATIONS'
            call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
            ! write matrix
            call fopen(funit, status='REPLACE', action='WRITE', file='corrmat_select.bin', access='STREAM', iostat=io_stat)
            call fileio_errmsg('simple_commander_preprocess ; error when opening corrmat_select.bin  ', io_stat)
            write(unit=funit,pos=1,iostat=io_stat) correlations
            ! Check if the write was successful
            if(io_stat/=0) then
                call fileio_errmsg('**ERROR(simple_commander_preprocess): I/O error writing corrmat_select.bin. Remove the file to override the memoization.', io_stat)
            endif
            call fclose(funit,errmsg='simple_commander_preprocess ; error when closing corrmat_select.bin  ')
        endif
        ! find selected
        ! in addition to the index array, also make a logical array encoding the selection (to be able to reject)
        allocate(selected(nsel), lselected(nall),stat=alloc_stat)
        allocchk("In commander_preprocess::select selected lselected ")
        lselected = .false.
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            selected(isel) = loc(1)
            lselected(selected(isel)) = .true.
            DebugPrint 'selected: ', loc(1), ' with corr: ', correlations(isel,loc(1))
        end do
        if( cline%defined('filetab') )then
            ! read filetable
            call read_filetable(p%filetab, imgnames)
            if( size(imgnames) /= nall ) stop 'nr of entries in filetab and stk not consistent'
            call fopen(funit, file=p%outfile,status="replace", action="write", access="sequential", iostat=io_stat)
            call fileio_errmsg('simple_commander_preprocess ; fopen error when opening '//trim(p%outfile), io_stat)
            call mkdir(p%dir_select)
            call mkdir(p%dir_reject)
            ! write outoput & move files
            do iimg=1,nall
                if( lselected(iimg) )then
                    write(funit,'(a)') trim(adjustl(imgnames(iimg)))
                    cmd_str = 'mv '//trim(adjustl(imgnames(iimg)))//' '//trim(adjustl(p%dir_select))
                    call exec_cmdline(cmd_str)
                else
                    cmd_str = 'mv '//trim(adjustl(imgnames(iimg)))//' '//trim(adjustl(p%dir_reject))
                    call exec_cmdline(cmd_str)
                endif
            end do
            call fclose(funit,errmsg='simple_commander_preprocess ; fopen error when closing '//trim(p%outfile))
            deallocate(imgnames)
        endif
        if( cline%defined('stk3') )then
            call find_ldim_nptcls(p%stk3, ldim, ifoo)
            ldim(3) = 1
            call stk3_img%new(ldim,1.)
            do isel=1,nsel
                call stk3_img%read(p%stk3, selected(isel))
                call stk3_img%write(p%outstk, isel)
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT NORMAL STOP ****')
    end subroutine exec_select

    !> for making picker references
    subroutine exec_make_pickrefs( self, cline )
        use simple_commander_volops,  only: project_commander
        class(make_pickrefs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline !< command line input
        integer,          parameter :: NREFS=100, NPROJS=20
        character(len=*), parameter :: ORIFILE='pickrefs_oris'//trim(TXT_EXT)
        type(params)                :: p
        type(build)                 :: b
        type(cmdline)               :: cline_project
        type(project_commander)     :: xproject
        integer                     :: nrots, cnt, iref, irot
        real                        :: ang, rot
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            p = params(cline) ! constants & derived constants produced
            if( cline%defined('vol1') )then
                p%nptcls = NPROJS
                call b%a%new(NPROJS)
                call b%a%spiral( p%nsym, p%eullims )
                call b%a%write(trim(ORIFILE), [1,NPROJS])
                cline_project = cline
                call cline_project%set('nspace', real(NPROJS))
                p%stk = 'even_projs'//p%ext
                call cline_project%set('outstk', trim(p%stk)  )
                call cline_project%set('oritab', trim(ORIFILE))
                call cline_project%set('smpd',   PICKER_SHRINK)
                call cline_project%set('neg',    'no'         )
                call cline_project%set('msk',    real(p%box/2-5))
                call xproject%execute(cline_project)
            endif
            ! expand in in-plane rotation
            nrots = NREFS/p%nptcls
            if( nrots > 1 )then
                ang = 360./real(nrots)
                rot = 0.
                cnt  = 0
                do iref=1,p%nptcls
                    call b%img%read(p%stk, iref)
                    do irot=1,nrots
                        cnt = cnt + 1
                        call b%img%rtsq(rot, 0., 0., b%img_copy)
                        call b%img_copy%write('rotated_from_make_pickrefs'//p%ext, cnt)
                        rot = rot + ang
                    end do
                end do
                call cline%set('stk', 'rotated_from_make_pickrefs'//p%ext)
            endif
            if( p%pcontrast .eq. 'black' )then
                call neg_imgfile('rotated_from_make_pickrefs'//p%ext, 'pickrefs'//p%ext, p%smpd)
            else
                call simple_rename('rotated_from_make_pickrefs'//p%ext, 'pickrefs'//p%ext)
            endif
        else
            stop 'need input volume (vol1) or class averages (stk) to generate picking references'
        endif
        call del_file('rotated_from_make_pickrefs'//p%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
    end subroutine exec_make_pickrefs

    subroutine exec_pick( self, cline)
        class(pick_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline !< command line input
        type(params)                         :: p
        type(picker_iter)                    :: piter
        character(len=STDLEN),   allocatable :: movienames_intg(:)
        character(len=:),        allocatable :: output_dir
        character(len=STDLEN)                :: boxfile
        integer :: nmovies, fromto(2), imovie, ntot, movie_counter, nptcls_out
        p = params(cline, spproj_a_seg=STK_SEG) ! constants & derived constants produced
        ! check filetab existence
        call read_filetable(p%filetab, movienames_intg)
        nmovies = size(movienames_intg)
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_PICKER)
        endif
        call mkdir(output_dir)
        ! determine loop range
        if( cline%defined('part') )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_pick'
            endif
        else
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
        endif
        ntot          = fromto(2) - fromto(1) + 1
        movie_counter = 0
        do imovie=fromto(1),fromto(2)
            call piter%iterate(cline, p, movie_counter, movienames_intg(imovie), boxfile, nptcls_out, output_dir)
            write(*,'(f4.0,1x,a)') 100.*(real(movie_counter)/real(ntot)), 'percent of the micrographs processed'
        end do
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_pick' )
        ! end gracefully
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick

    !> for extracting particle images from integrated DDD movies
    subroutine exec_extract( self, cline )
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nboxfiles, nframes, pind, noris
        integer                            :: i, j, ldim(3), box_current, movie, ndatlines, nptcls
        integer                            :: cnt, niter, fromto(2), orig_box
        integer                            :: movie_ind, ntot, lfoo(3), ifoo, noutside
        type(nrtxtfile)                    :: boxfile
        character(len=STDLEN)              :: stack, outfile, moviename, stack_frames, fbody
        character(len=STDLEN), allocatable :: movienames(:), boxfilenames(:), movienames_frames(:)
        character(len=:),      allocatable :: output_dir
        real,                  allocatable :: boxdata(:,:)
        logical,               allocatable :: oris_mask(:)
        real                               :: kv, cs, fraca, dfx, dfy, angast, ctfres
        real                               :: particle_position(2), ctf_estimatecc, phshift
        type(image)                        :: micrograph
        type(oris)                         :: outoris, os_uni
        logical                            :: params_present(3), ctf_estimatecc_is_there, phshift_is_there
        logical                            :: ctfres_is_there
        noutside = 0
        p = params(cline, spproj_a_seg=STK_SEG) ! constants & derived constants produced
        if( p%stream .eq. 'yes' )then
            if( cline%defined('outstk') .and. cline%defined('outfile') )then
                ! all ok
            else
                stop 'need outstk & outfile to be part of the command line in stream=yes mode'
            endif
        endif
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_EXTRACT)
        endif
        call mkdir(output_dir)
        ! check file inout existence and read filetables
        if( .not. cline%defined('filetab') ) stop 'need filetab input to extract'
        if( .not. cline%defined('boxtab')  ) stop 'need boxtab input to extract'
        if( .not. file_exists(p%filetab)   ) stop 'inputted filetab does not exist in cwd'
        if( .not. file_exists(p%boxtab)    ) stop 'inputted boxtab does not exist in cwd'
        nmovies   = nlines(p%filetab)
        nboxfiles = nlines(p%boxtab)
        DebugPrint  'nboxfiles: ', nboxfiles
        if( nmovies /= nboxfiles ) stop 'number of entries in inputted files do not match!'
        call read_filetable(p%filetab, movienames)
        call read_filetable(p%boxtab,  boxfilenames)
        DebugPrint  'nmovies: ', nmovies
        ! determine loop range
        fromto(1) = 1
        fromto(2) = nmovies
        ntot = fromto(2) - fromto(1) + 1
        ! find ldim
        call find_ldim_nptcls(movienames(1), ldim, ifoo)
        ! initialize
        call micrograph%new([ldim(1),ldim(2),1], p%smpd)
        niter    = 0
        noutside = 0
        ! loop over exposures (movies)
        do movie=fromto(1),fromto(2)
            moviename = trim(adjustl(movienames(movie)))
            ! get movie index (parsing the number string from the filename)
            call fname2ind(moviename, movie_ind)
            ! process boxfile
            nptcls = 0
            if( file_exists(boxfilenames(movie)) )then
                if( nlines(boxfilenames(movie)) > 0 )then
                    call boxfile%new(boxfilenames(movie), 1)
                    nptcls = boxfile%get_ndatalines()
                endif
            endif
            if( nptcls == 0 ) cycle
            ! update iteration counter
            niter = niter + 1
            ! show progress
            if( niter > 1 ) call progress(niter,ntot)
            ! oris object for params
            call outoris%new(nptcls)
            if(allocated(oris_mask))deallocate(oris_mask)
            allocate(oris_mask(nptcls), source=.false., stat=alloc_stat)
            allocchk("In simple_extract oris_mask")
            ! read box data & update mask
            if(allocated(boxdata))deallocate(boxdata)
            allocate( boxdata(nptcls,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            allocchk('In: simple_extract; boxdata etc., 2')
            do j=1,nptcls
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    stop 'Only square windows are currently allowed!'
                endif
                if( j == 1 .and. .not. cline%defined('box') ) p%box = nint(boxdata(1,3)) ! use the box size from the box file
                ! modify coordinates if change in box (shift by half the difference)
                if( orig_box /= p%box ) boxdata(j,1:2) = boxdata(j,1:2)-real(p%box-orig_box)/2.
                ! update particle mask & movie index
                if( box_inside(ldim, nint(boxdata(j,1:2)), p%box) )then
                    call outoris%set(j, 'movie', real(movie_ind))
                    oris_mask(j) = .true.
                endif
            end do
            ! check box parsing
            if( .not. cline%defined('box') )then
                if( niter == 1 )then
                    p%box = nint(boxdata(1,3))   ! use the box size from the box file
                else
                    box_current = nint(boxdata(1,3))
                    if( box_current /= p%box )then
                        write(*,*) 'box_current: ', box_current, 'box in params: ', p%box
                        stop 'inconsistent box sizes in box files'
                    endif
                endif
                if( p%box == 0 )then
                    write(*,*) 'ERROR, box cannot be zero!'
                    stop
                endif
            endif
            DebugPrint  'did check box parsing'
            ! get number of frames from stack
            call find_ldim_nptcls(moviename, lfoo, nframes )
            if( nframes > 1 ) stop 'multi-frame extraction no longer supported; simple_extract'
            ! build general objects
            if( niter == 1 )then
                call b%build_general_tbox(p,cline,do3d=.false.)
            endif
            ! extract ctf info
            if( b%a%isthere('dfx') )then
                params_present(1) = b%a%isthere('kv')
                params_present(2) = b%a%isthere('cs')
                params_present(3) = b%a%isthere('fraca')
                if( all(params_present) )then
                    ! alles ok
                else
                    if( .not. params_present(1) ) write(*,*) 'ERROR! input doc lacks kv'
                    if( .not. params_present(2) ) write(*,*) 'ERROR! input doc lacks cs'
                    if( .not. params_present(3) ) write(*,*) 'ERROR! input doc lacks fraca'
                    stop
                endif
                kv    = b%a%get(movie,'kv')
                cs    = b%a%get(movie,'cs')
                fraca = b%a%get(movie,'fraca')
                dfx   = b%a%get(movie,'dfx')
                ctf_estimatecc_is_there = b%a%isthere('ctf_estimatecc')
                phshift_is_there  = b%a%isthere('phshift')
                ctfres_is_there   = b%a%isthere('ctfres')
                if( ctf_estimatecc_is_there ) ctf_estimatecc = b%a%get(movie,'ctf_estimatecc')
                if( phshift_is_there  ) phshift  = b%a%get(movie,'phshift')
                if( ctfres_is_there   ) ctfres   = b%a%get(movie,'ctfres')
                angast = 0.
                if( b%a%isthere('dfy') )then ! astigmatic CTF
                    if( .not. b%a%isthere('angast') ) stop 'need angle of astigmatism for CTF correction'
                    dfy    = b%a%get(movie,'dfy')
                    angast = b%a%get(movie,'angast')
                endif
                do j=1,nptcls
                    if( oris_mask(j) )then
                        call outoris%set(j, 'kv',         kv)
                        call outoris%set(j, 'cs',         cs)
                        call outoris%set(j, 'fraca',   fraca)
                        call outoris%set(j, 'smpd',   p%smpd)
                        call outoris%set(j, 'dfx',       dfx)
                        if( b%a%isthere('dfy') )then
                            call outoris%set(j, 'angast', angast)
                            call outoris%set(j, 'dfy',       dfy)
                        endif
                        if( ctf_estimatecc_is_there ) call outoris%set(j, 'ctf_estimatecc', ctf_estimatecc)
                        if( phshift_is_there  ) call outoris%set(j, 'phshift',  phshift)
                        if( ctfres_is_there   ) call outoris%set(j, 'ctfres',   ctfres)
                    endif
                end do
                DebugPrint  'did set CTF parameters dfx/dfy/angast/ctfres: ', dfx, dfy, angast, ctfres
                ! compress & write params (remove entries for boxes outside micrograph)
                call outoris%compress(oris_mask)
                call outoris%kill_chash() ! remove chash part
                ! output stk & params file names
                if( p%stream.eq.'yes' )then
                    outfile = trim(output_dir)//p%outfile
                else
                    fbody   = get_fbody(trim(remove_abspath(moviename)), p%ext, separator=.false.)
                    outfile = trim(output_dir)//trim(EXTRACT_PARAMS_FBODY)//trim(fbody)//trim(METADATA_EXT)
                endif
                call del_file(outfile)
                call binwrite_oritab(outfile, b%spproj, outoris, [1,count(oris_mask)])
            endif
            ! output stack
            if( p%stream.eq.'yes' )then
                stack   = trim(output_dir)//p%outstk
            else
                stack   = trim(output_dir)//trim(EXTRACT_STK_FBODY) // trim(remove_abspath(moviename))
            endif
            ! extract windows from integrated movie
            call micrograph%read(moviename, 1)
            ! filter out frequencies lower than the box can express to avoid aliasing
            call micrograph%bp(real(p%box) * p%smpd, 0.)
            cnt = 0
            do j=1,nptcls ! loop over boxes
                if( oris_mask(j) )then
                    cnt = cnt + 1
                    ! extract the window
                    particle_position = boxdata(j,1:2)
                    call micrograph%window(nint(particle_position), p%box, b%img, noutside)
                    if( p%pcontrast .eq. 'black' ) call b%img%neg()
                    call b%img%norm()
                    call b%img%edges_norm()
                    call b%img%write(trim(adjustl(stack)), cnt)
                endif
            end do
            ! extract windows from integrated movie (subset of frames)
            if( allocated(movienames_frames) )then
                fbody        = get_fbody(trim(stack), p%ext, separator=.false.)
                stack_frames = add2fbody(fbody, p%ext, '_frames_subset')
                call micrograph%read(movienames_frames(movie),1)
                cnt = 0
                do j=1,nptcls ! loop over boxes
                    if( oris_mask(j) )then
                        cnt = cnt + 1
                        ! extract the window
                        particle_position = boxdata(j,1:2)
                        call micrograph%window(nint(particle_position), p%box, b%img)
                        if( p%pcontrast .eq. 'black' ) call b%img%neg()
                        call b%img%norm()
                        call b%img%edges_norm
                        call b%img%write(trim(adjustl(stack_frames)), cnt)
                    endif
                end do
            endif
            ! destruct
            call boxfile%kill
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_EXTRACT NORMAL STOP ****')

        contains

            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                if( p%outside .eq. 'yes' )then
                    inside = .true.
                    return
                endif
                fromc  = coord+1       ! compensate for the c-range that starts at 0
                toc    = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.        ! box is inside
                if( any(fromc < 1) .or. toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function box_inside

    end subroutine exec_extract

end module simple_commander_preprocess
