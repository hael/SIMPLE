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
use simple_sp_project,          only: sp_project
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
        class(cmdline),              intent(inout) :: cline
        type(params)                  :: p
        type(ori)                     :: o_mov
        type(ctf_estimate_iter)       :: cfiter
        type(motion_correct_iter)     :: mciter
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        character(len=:), allocatable :: imgkind, moviename, output_dir_picker, fbody
        character(len=:), allocatable :: moviename_forctf, moviename_intg, output_dir_motion_correct
        character(len=:), allocatable :: output_dir, output_dir_ctf_estimate, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        real    :: smpd_original, smpd_scaled
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, nptcls_out
        logical :: l_pick
        ! constants & derived constants
        p = params(cline, spproj_a_seg=MIC_SEG)
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
        ! output directories & name
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
        call mkdir(output_dir)
        call mkdir(output_dir_ctf_estimate)
        call mkdir(output_dir_motion_correct)
        if( l_pick )then
            output_dir_picker  = trim(output_dir)//trim(DIR_PICKER)
            call mkdir(output_dir_picker)
            if( p%stream.eq.'yes' )then
                output_dir_extract = trim(output_dir)//trim(DIR_EXTRACT)
                call mkdir(output_dir_extract)
            endif
        endif
        if( cline%defined('fbody') )then
            fbody = trim(p%fbody)
        else
            fbody = ''
        endif
        ! read in movies
        call spproj%read( p%projfile )
        ! range
        if( p%stream.eq.'yes' )then
            ! STREAMING MODE
            fromto(:) = 1
        else
            ! DISTRIBUTED MODE
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; exec_preprocess'
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! numlen
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            p%numlen = len(int2str(nmovies))
        endif
        !
        frame_counter = 0
        smpd_original = p%smpd
        ! loop over exposures (movies)
        do imovie = fromto(1),fromto(2)
            ! fetch movie orientation
            o_mov = spproj%os_mic%get_ori(imovie)
            ! sanity check
            if( .not.o_mov%isthere('imgkind') )cycle
            if( .not.o_mov%isthere('movie')   )cycle
            call o_mov%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'movie' )cycle
            call o_mov%getter('movie', moviename)
            if( .not.file_exists(moviename)) cycle
            !
            p%smpd    = smpd_original
            ! motion_correct
            call mciter%iterate(cline, p, o_mov, fbody, frame_counter, moviename,&
                &smpd_scaled, output_dir_motion_correct)
            p%smpd = smpd_scaled
            ! ctf_estimate
            moviename_forctf = mciter%get_moviename('forctf')
            p%hp             = p%hp_ctf_estimate
            p%lp             = max(p%fny, p%lp_ctf_estimate) ! should be in params?
            call cfiter%iterate(p, moviename_forctf, o_mov, output_dir_ctf_estimate)
            ! update project
            call spproj%os_mic%set_ori(imovie, o_mov)
            ! picker
            if( l_pick )then
                p%lp = max(p%fny, p%lp_pick) ! should be in params?
                moviename_intg = mciter%get_moviename('intg')
                call piter%iterate(cline, p, moviename_intg, boxfile, nptcls_out, output_dir_picker)
                call o_mov%set('boxfile', trim(boxfile)   )
                call o_mov%set('nptcls',  real(nptcls_out))
                ! update project
                call spproj%os_mic%set_ori(imovie, o_mov)
                ! extract particles
                if( p%stream .eq. 'yes' )then
                    ! needs to write and re-read project at the end as extract overwrites it
                    call spproj%write()
                    cline_extract = cline
                    call cline_extract%set('dir', trim(output_dir_extract))
                    call cline_extract%set('pcontrast', p%pcontrast)
                    if( cline%defined('box_extract') )call cline_extract%set('box', real(p%box_extract))
                    call xextract%execute(cline_extract)
                    call spproj%kill
                endif
            endif
        end do
        if( p%stream .eq. 'yes' )then
            call spproj%print_info
            if( .not.l_pick )call spproj%write()
        else
            call binwrite_oritab(p%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        ! end gracefully
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_preprocess' )
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
        class(motion_correct_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(params)                  :: p
        type(motion_correct_iter)     :: mciter
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename, imgkind, fbody
        real    :: smpd_scaled
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, lfoo(3), nframes, cnt
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
        ! output directory & names
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_MOTION_CORRECT)
        endif
        call mkdir(output_dir)
        if( cline%defined('fbody') )then
            fbody = trim(p%fbody)
        else
            fbody = ''
        endif
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
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! for series of tomographic movies we need to calculate the time_per_frame
        if( p%tomo .eq. 'yes' )then
            ! get number of frames & dim from stack
            call spproj%os_mic%getter(1, 'movie', moviename)
            call find_ldim_nptcls(moviename, lfoo, nframes)
            ! calculate time_per_frame
            p%time_per_frame = p%exp_time/real(nframes*nmovies)
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
                call mciter%iterate(cline, p, o, fbody, frame_counter, moviename, smpd_scaled, trim(output_dir))
                call spproj%os_mic%set_ori(imovie, o)
                write(*,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(p%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_motion_correct' )
        call simple_end('**** SIMPLE_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct

    subroutine exec_ctf_estimate( self, cline )
        class(ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline  !< command line input
        type(params)                  :: p
        type(ctf_estimate_iter)       :: cfiter
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: intg_forctf, output_dir, imgkind
        integer                       :: fromto(2), imic, ntot, cnt
        p = params(cline, spproj_a_seg=MIC_SEG) ! constants & derived constants produced
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
            fromto(:) = 1
            call spproj%os_mic%getter(1, 'intg', intg_forctf)
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_ctf_estimate'
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! read in integrated movies
        call spproj%read_segment( 'mic', p%projfile, fromto )
        ! loop over exposures (movies)
        cnt= 0
        do imic = fromto(1),fromto(2)
            o = spproj%os_mic%get_ori(imic)
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                cnt = cnt + 1
                call o%getter('intg', intg_forctf)
                call cfiter%iterate(p, intg_forctf, o, trim(output_dir))
                call spproj%os_mic%set_ori(imic, o)
                write(*,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
            endif
        end do
        ! output
        call binwrite_oritab(p%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_ctf_estimate' )
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
        integer                            :: iimg, isel, nall, nsel, loc(1)
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
        type(sp_project)                     :: spproj
        type(picker_iter)                    :: piter
        type(ori)                            :: o
        character(len=:),        allocatable :: output_dir, intg_name, imgkind
        character(len=LONGSTRLEN)            :: boxfile
        integer :: fromto(2), imic, ntot, nptcls_out, cnt
        p = params(cline, spproj_a_seg=MIC_SEG) ! constants & derived constants produced
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_PICKER)
        endif
        call mkdir(output_dir)
        ! parameters & loop range
        if( p%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:)   = 1 !!!!!!!!!!!!!!!!!!!!!!! NOW SHOULD POINT TO MICROGRAPH OF ORIGIN ????
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_pick'
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! read in integrated movies
        call spproj%read_segment('mic', p%projfile, fromto)
        ! main loop
        cnt = 0
        do imic=fromto(1),fromto(2)
            o = spproj%os_mic%get_ori(imic)
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                cnt = cnt + 1
                call o%getter('intg', intg_name)
                call piter%iterate(cline, p, intg_name, boxfile, nptcls_out, output_dir)
                call spproj%os_mic%set(imic, 'boxfile', trim(boxfile))
                call spproj%os_mic%set(imic, 'nptcls', real(nptcls_out))
                write(*,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
            endif
        end do
        ! output
        call binwrite_oritab(p%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! end gracefully
        call qsys_job_finished( p, 'simple_commander_preprocess :: exec_pick' )
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick

    !> for extracting particle images from integrated DDD movies
    subroutine exec_extract( self, cline )
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(params)                  :: p
        type(build)                   :: b
        type(sp_project)              :: spproj
        type(nrtxtfile)               :: boxfile
        type(image)                   :: micrograph
        type(oris)                    :: o_ptcls
        type(ori)                     :: o_mic
        type(ctfparams)               :: ctfparms
        character(len=:), allocatable :: output_dir, mic_name, boxfile_name, imgkind, ctfstr, phplate
        real,             allocatable :: boxdata(:,:)
        logical,          allocatable :: oris_mask(:), mics_mask(:)
        character(len=LONGSTRLEN)     :: stack
        integer                       :: nframes, imic, iptcl, ldim(3), nptcls, nmics, box
        integer                       :: cnt, niter, ntot, lfoo(3), ifoo, noutside, nptcls_eff
        real                          :: particle_position(2)
        p = params(cline, spproj_a_seg=MIC_SEG) ! constants & derived constants produced
        call cline%printline
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p%dir)//'/'
        else
            output_dir = trim(DIR_EXTRACT)
        endif
        call mkdir(output_dir)
        ! read in integrated movies
        call spproj%read_segment('mic', p%projfile)
        ntot  = spproj%os_mic%get_noris()
        ! sanity checks
        allocate(mics_mask(ntot), source=.false.)
        nmics  = 0
        nptcls = 0
        do imic = 1, ntot
            o_mic = spproj%os_mic%get_ori(imic)
            if( .not.o_mic%isthere('imgkind') )cycle
            if( .not.o_mic%isthere('intg')    )cycle
            if( .not.o_mic%isthere('boxfile') )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            call o_mic%getter('boxfile', boxfile_name)
            if( .not.file_exists(boxfile_name) )cycle
            ! get number of frames from stack
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) stop 'multi-frame extraction no longer supported; commander_preproc :: exec_extract'
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
                call cline%set('box', boxdata(1,3))
            endif
        enddo
        call spproj%kill
        if( nmics == 0 ) stop 'No particles to extract! commander_preproc :: exec_extract'
        p%box = nint(cline%get_rarg('box'))
        if( p%box == 0 )stop 'ERROR! box cannot be zero; commander_preproc :: exec_extract'
        ! init
        call b%build_general_tbox(p, cline, do3d=.false.)
        call micrograph%new([ldim(1),ldim(2),1], p%smpd)
        niter    = 0
        noutside = 0
        ! main loop
        do imic = 1, ntot
            if( .not.mics_mask(imic) )cycle
            ! fetch micrograph
            o_mic = b%a%get_ori(imic)
            call o_mic%getter('imgkind', imgkind)
            call o_mic%getter('intg', mic_name)
            call o_mic%getter('boxfile', boxfile_name)
            ! box file
            nptcls = 0
            if( nlines(boxfile_name) > 0 )then
                call boxfile%new(boxfile_name, 1)
                nptcls = boxfile%get_ndatalines()
            endif
            if( nptcls == 0 ) cycle
            ! ...
            niter = niter + 1
            call progress(imic,ntot)
            ! box checks
            call o_ptcls%new(nptcls)
            if(allocated(oris_mask))deallocate(oris_mask)
            allocate(oris_mask(nptcls), source=.false., stat=alloc_stat)
            ! read box data & update mask
            if( allocated(boxdata) )deallocate(boxdata)
            allocate( boxdata(nptcls, boxfile%get_nrecs_per_line()), stat=alloc_stat)
            do iptcl=1,nptcls
                call boxfile%readNextDataLine(boxdata(iptcl,:))
                box = nint(boxdata(iptcl,3))
                if( nint(boxdata(iptcl,3)) /= nint(boxdata(iptcl,4)) )then
                    stop 'ERROR! Only square windows are currently allowed; commander_preproc :: exec_extract'
                endif
                ! modify coordinates if change in box (shift by half the difference)
                if( box /= p%box ) boxdata(iptcl,1:2) = boxdata(iptcl,1:2) - real(p%box-box)/2.
                if( nint(boxdata(iptcl,3)) /= p%box )then
                    write(*,*) 'box_current: ', nint(boxdata(iptcl,3)), 'box in params: ', p%box
                    stop 'ERROR! inconsistent box sizes in box files; commander_preproc :: exec_extract'
                endif
                ! update particle mask & movie index
                if( box_inside(ldim, nint(boxdata(iptcl,1:2)), p%box) )oris_mask(iptcl) = .true.
            end do
            nptcls_eff = count(oris_mask)
            call o_ptcls%new_clean(nptcls)
            ! extract ctf info
            if( o_mic%isthere('dfx') )then
                if( .not.o_mic%isthere('cs') .or. .not.o_mic%isthere('kv') .or. .not.o_mic%isthere('fraca') )then
                    stop 'ERROR! Input lacks at least cs, kv or fraca field; commander_preproc :: exec_extract'
                endif
                ctfparms%smpd = p%smpd
                ! prepare CTF vars
                call o_mic%getter('ctf', ctfstr)
                ctfparms%ctfflag = 1
                select case( trim(ctfstr) )
                    case('no')
                        ctfparms%ctfflag = 0
                    case('yes')
                        ctfparms%ctfflag = 1
                    case('flip')
                        ctfparms%ctfflag = 2
                end select
                ctfparms%kv           = o_mic%get('kv')
                ctfparms%cs           = o_mic%get('cs')
                ctfparms%fraca        = o_mic%get('fraca')
                if( o_mic%isthere('phaseplate'))then
                    call o_mic%getter('phaseplate', phplate)
                    ctfparms%l_phaseplate = trim(phplate) .eq. 'yes'
                else
                    ctfparms%l_phaseplate = .false.
                endif
                if( ctfparms%ctfflag > 0 )then
                    if( .not.o_mic%isthere('dfx') )then
                        stop 'ERROR! ctf .ne. no and input lacks dfx; commander_preproc :: exec_extract'
                    endif
                endif
                ! transfer to particles
                do iptcl=1,nptcls
                    if( .not.oris_mask(iptcl) )cycle
                    call o_ptcls%set_ori(iptcl, o_mic)
                end do
                call o_ptcls%compress(oris_mask)
            endif
            ! output stack
            stack = trim(output_dir)//trim(EXTRACT_STK_FBODY)//trim(remove_abspath(mic_name))
            ! extract windows from integrated movie
            call micrograph%read(mic_name, 1)
            ! filter out frequencies lower than the box can express to avoid aliasing
            call micrograph%bp(real(p%box) * p%smpd, 0.)
            ! write new stack
            cnt = 0
            do iptcl=1,nptcls ! loop over boxes
                if( oris_mask(iptcl) )then
                    cnt = cnt + 1
                    ! extract the window
                    particle_position = boxdata(iptcl,1:2)
                    call micrograph%window(nint(particle_position), p%box, b%img, noutside)
                    if( p%pcontrast .eq. 'black' ) call b%img%neg()
                    call b%img%norm()
                    call b%img%edges_norm()
                    call b%img%write(trim(adjustl(stack)), cnt)
                endif
            end do
            ! IMPORT INTO PROJECT
            call b%spproj%add_stk(trim(adjustl(stack)), ctfparms, o_ptcls)
            ! clean
            call boxfile%kill()
        enddo
        call b%spproj%write()
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
