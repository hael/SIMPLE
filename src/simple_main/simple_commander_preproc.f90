!==Class simple_commander_preproc
!
! This class contains the set of concrete preprocessing commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_preproc
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_syscalls,       only: exec_cmdline
use simple_jiffys          ! use all in there
use simple_filehandling    ! use all in there
implicit none

public :: preproc_commander
public :: select_frames_commander
public :: boxconvs_commander
public :: powerspecs_commander
public :: unblur_commander
public :: ctffind_commander
public :: select_commander
public :: makepickrefs_commander
public :: pick_commander
public :: extract_commander
private

type, extends(commander_base) :: preproc_commander
  contains
    procedure :: execute      => exec_preproc
end type preproc_commander
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
type, extends(commander_base) :: unblur_commander
  contains
    procedure :: execute      => exec_unblur
end type unblur_commander
type, extends(commander_base) :: ctffind_commander
  contains
    procedure :: execute      => exec_ctffind
end type ctffind_commander
type, extends(commander_base) :: select_commander
  contains
    procedure :: execute      => exec_select
end type select_commander
type, extends(commander_base) :: makepickrefs_commander
  contains
    procedure :: execute      => exec_makepickrefs
end type makepickrefs_commander
type, extends(commander_base) :: pick_commander
  contains
    procedure :: execute      => exec_pick
end type pick_commander
type, extends(commander_base) :: extract_commander
  contains
    procedure :: execute      => exec_extract
end type extract_commander

contains

    ! UNBLUR + CTFFIND + PICK + EXTRACT IN SEQUENCE

    subroutine exec_preproc( self, cline )
        use simple_unblur_iter,  only: unblur_iter
        use simple_ctffind_iter, only: ctffind_iter
        use simple_pick_iter,    only: pick_iter
        use simple_oris,         only: oris
        use simple_ori,          only: ori
        use simple_math,         only: round2even
        use simple_qsys_funs,    only: qsys_job_finished
        class(preproc_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(ctffind_iter)      :: cfiter
        type(unblur_iter)       :: ubiter
        type(pick_iter)         :: piter
        type(extract_commander) :: xextract
        type(cmdline)           :: cline_extract
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:),      allocatable :: fname_ctffind_ctrl, fname_unidoc_output
        character(len=:),      allocatable :: moviename_forctf, moviename_intg, outfile
        logical, parameter    :: DEBUG = .true.
        character(len=STDLEN) :: boxfile, dir_ptcls
        type(params) :: p
        type(oris)   :: os_uni
        type(ori)    :: orientation
        real         :: smpd_native, smpd_scaled
        integer      :: nmovies, fromto(2), imovie, ntot, movie_counter
        integer      :: frame_counter, lfoo(3), nframes, i, movie_ind, nptcls_out
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced
        if( p%scale > 1.05 )then
            stop 'scale cannot be > 1; simple_commander_preproc :: exec_preproc'
        endif
        if( p%tomo .eq. 'yes' )then
            stop 'tomography mode (tomo=yes) not yet supported!'
        endif
        if( p%l_pick )then
            if( .not. cline%defined('refs') )then
                stop 'need references for picker or turn off picking with dopick=no'
            endif
        endif
        dir_ptcls = trim(p%dir_target)//'/particles'
        call exec_cmdline('mkdir -p '//trim(adjustl(dir_ptcls)))
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            p%numlen = len(int2str(nmovies))
        endif
        if( p%l_distr_exec )then
            if( cline%defined('dir_target') )then
                allocate(fname_ctffind_ctrl,  source=trim(p%dir_target)//'/'//&
                &'ctffind_ctrl_file_part'//int2str_pad(p%part,p%numlen)//'.txt')
                if( cline%defined('outfile') )then
                    allocate(fname_unidoc_output, source=trim(p%dir_target)//'/'//&
                    &trim(p%outfile))
                else
                    allocate(fname_unidoc_output, source=trim(p%dir_target)//'/'//&
                    &'unidoc_output_part'//int2str_pad(p%part,p%numlen)//'.txt')
                endif
            else
                allocate(fname_ctffind_ctrl,  source='ctffind_ctrl_file_part'//&
                &int2str_pad(p%part,p%numlen)//'.txt')
                if( cline%defined('outfile') )then
                    allocate(fname_unidoc_output, source=trim(p%outfile))
                else
                    allocate(fname_unidoc_output, source='unidoc_output_part'//&
                    &int2str_pad(p%part,p%numlen)//'.txt')
                endif
            endif
            ! determine loop range
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; exec_preproc'
            endif
        else
            allocate(fname_ctffind_ctrl,  source='ctffind_ctrl_file.txt')
            allocate(fname_unidoc_output, source='unidoc_output.txt')
            ! determine loop range
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
        endif
        ntot = fromto(2) - fromto(1) + 1
        frame_counter = 0
        movie_counter = 0
        call orientation%new
        call os_uni%new(ntot)
        smpd_native = p%smpd
        ! loop over exposures (movies)
        do imovie=fromto(1),fromto(2)
            p%smpd    = smpd_native
            p%pspecsz = p%pspecsz_unblur
            if( ntot == 1 )then
                movie_ind = p%part ! streaming mode
            else
                movie_ind = imovie ! standard mode
            endif 
            call ubiter%iterate(cline, p, orientation, movie_ind, movie_counter,&
            &frame_counter, movienames(imovie), smpd_scaled)
            call os_uni%set_ori(movie_counter, orientation)
            p%smpd           = smpd_scaled     
            movie_counter    = movie_counter - 1
            moviename_forctf = ubiter%get_moviename('forctf')
            moviename_intg   = ubiter%get_moviename('intg')
            p%pspecsz        = p%pspecsz_ctffind
            p%hp             = p%hp_ctffind
            p%lp             = p%lp_ctffind 
            call cfiter%iterate(p, movie_ind, movie_counter, moviename_forctf,&
            &fname_ctffind_ctrl, fname_unidoc_output, os_uni)
            if( p%l_pick )then
                movie_counter = movie_counter - 1
                p%lp          = p%lp_pick
                call piter%iterate(cline, p, movie_counter, moviename_intg, boxfile, nptcls_out)
                call os_uni%set(movie_counter, 'boxfile', trim(boxfile)   )
                call os_uni%set(movie_counter, 'nptcls',  real(nptcls_out))
            endif
            if( p%stream .eq. 'yes' )then
                ! write unidoc
                call os_uni%write(fname_unidoc_output)
                cline_extract = cline
                call cline_extract%set('dir_ptcls', trim(dir_ptcls))
                call cline_extract%set('smpd',      p%smpd)
                call cline_extract%set('unidoc',    fname_unidoc_output)
                call cline_extract%set('outfile',   'extract_params_movie'//int2str_pad(movie_ind,p%numlen)//'.txt')
                call cline_extract%set('outstk',    'ptcls_from_movie'//int2str_pad(movie_ind,p%numlen)//p%ext)
                call xextract%execute(cline_extract)
            endif
        end do
        if( p%stream .eq. 'no' )then
            ! write unidoc
            call os_uni%write(fname_unidoc_output)
        endif
        ! destruct
        call os_uni%kill
        deallocate(fname_ctffind_ctrl,fname_unidoc_output)
        call qsys_job_finished( p, 'simple_commander_preproc :: exec_preproc' )
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROC NORMAL STOP ****')
    end subroutine exec_preproc

    subroutine exec_select_frames( self, cline )
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(select_frames_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nframes, frame, numlen, ldim(3)
        integer                            :: fromto(2), movie, cnt, cnt2, ntot, lfoo(3), ifoo
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:), allocatable      :: new_name
        type(image)                        :: img_frame
        logical, parameter                 :: debug = .false.
        p = params(cline, checkdistr=.false.)            ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general objects built
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        ! find ldim and numlen (length of number string)
        if( cline%defined('startit') )then
            call find_ldim_nptcls(movienames(p%startit), ldim, ifoo)
        endif
        if( debug ) write(*,*) 'logical dimension: ', ldim
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        numlen = len(int2str(nmovies))
        if( debug ) write(*,*) 'length of number string: ', numlen
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
        if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
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
            if( debug ) write(*,*) 'number of frames: ', nframes
            ! create new name
            allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
            cnt = 0
            do frame=p%fromp,p%top
                cnt = cnt+1
                call img_frame%read(movienames(movie),frame,rwaction='READ')
                call img_frame%write(new_name,cnt)
            end do
            deallocate(new_name)
            write(*,'(f4.0,1x,a)') 100.*(real(cnt2)/real(ntot)), 'percent of the movies processed'
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT_FRAMES NORMAL STOP ****')
    end subroutine exec_select_frames
    
    subroutine exec_boxconvs( self, cline )
        use simple_image, only: image
        class(boxconvs_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: tmp
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        logical                            :: debug=.false.
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline, checkdistr=.false.)               ! parameters generated
        call b%build_general_tbox( p, cline, do3d=.false. ) ! general objects built
        ! do the work
        if( cline%defined('stk') )then
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg, rwaction='READ')
                tmp = b%img%boxconv(p%boxconvsz)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            if( debug ) write(*,*) 'read the img filenames'
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( debug ) write(*,*) 'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            iimg_start = 1
            if( cline%defined('startit') ) iimg_start = p%startit
            iimg_stop  = nimgs
            if( debug ) write(*,*) 'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1, rwaction='READ')
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
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(powerspecs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: powspec, tmp, mask
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        logical                            :: debug=.false.
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline, checkdistr=.false.)           ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general toolbox built
        ! create mask
        call tmp%new([p%clip,p%clip,1], p%smpd)
        tmp = cmplx(1.,0.)
        call tmp%bp(0.,p%lp,0.)
        call tmp%ft2img('real', mask)
        call mask%write('resolution_mask.mrc', 1)
        ! do the work
        if( cline%defined('stk') )then
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg, rwaction='READ')
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            if( debug ) write(*,*) 'read the img filenames'
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( debug ) write(*,*) 'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            iimg_start = 1
            if( cline%defined('startit') ) iimg_start = p%startit
            iimg_stop  = nimgs
            if( debug ) write(*,*) 'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1, rwaction='READ')
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, nimgs)
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_POWERSPECS NORMAL STOP ****')
    end subroutine exec_powerspecs

    subroutine exec_unblur( self, cline )
        use simple_unblur_iter, only: unblur_iter
        use simple_oris,        only: oris
        use simple_ori,         only: ori
        use simple_math,        only: round2even
        use simple_qsys_funs,   only: qsys_job_finished
        class(unblur_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)      :: p
        type(unblur_iter) :: ubiter
        type(oris)        :: os_uni
        type(ori)         :: orientation
        character(len=STDLEN), allocatable :: movienames(:)
        real    :: smpd_scaled
        integer :: nmovies, fromto(2), imovie, ntot, movie_counter
        integer :: frame_counter, lfoo(3), nframes
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced
        if( p%scale > 1.05 )then
            stop 'scale cannot be > 1; simple_commander_preproc :: exec_unblur'
        endif
        if( p%tomo .eq. 'yes' )then
            if( .not. p%l_dose_weight )then
                write(*,*) 'tomo=yes only supported with dose weighting!'
                stop 'give total exposure time: exp_time (in seconds) and dose_rate (in e/A2/s)'
            endif
        endif
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            p%numlen = len(int2str(nmovies))
        endif
        ! determine loop range
        if( p%l_distr_exec )then
            if( p%tomo .eq. 'no' )then
                if( cline%defined('fromp') .and. cline%defined('top') )then
                    fromto(1) = p%fromp
                    fromto(2) = p%top
                else
                    stop 'fromp & top args need to be defined in parallel execution; simple_unblur'
                endif
            else
                fromto(1) = 1
                fromto(2) = nmovies
            endif
        else
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2)  = nmovies
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! for series of tomographic movies we need to calculate the time_per_frame
        if( p%tomo .eq. 'yes' )then
            ! get number of frames & dim from stack
            call find_ldim_nptcls(movienames(1), lfoo, nframes, endconv=endconv)
            ! calculate time_per_frame
            p%time_per_frame = p%exp_time/real(nframes*nmovies)
        endif
        ! align
        frame_counter = 0
        movie_counter = 0
        call orientation%new
        call os_uni%new(ntot)
        do imovie=fromto(1),fromto(2)
            call ubiter%iterate(cline, p, orientation, imovie, movie_counter,&
            &frame_counter, movienames(imovie), smpd_scaled)
            call os_uni%set_ori(movie_counter, orientation)
            write(*,'(f4.0,1x,a)') 100.*(real(movie_counter)/real(ntot)), 'percent of the movies processed'
        end do
        call os_uni%write(p%outfile)
        call qsys_job_finished( p, 'simple_commander_preproc :: exec_unblur' )
        ! end gracefully
        call simple_end('**** SIMPLE_UNBLUR NORMAL STOP ****')
    end subroutine exec_unblur

    subroutine exec_ctffind( self, cline )
        use simple_ctffind_iter, only: ctffind_iter
        use simple_oris,         only: oris
        use simple_qsys_funs,    only: qsys_job_finished
        class(ctffind_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                       :: p
        type(ctffind_iter)                 :: cfiter
        character(len=STDLEN), allocatable :: movienames_forctf(:)
        character(len=:),      allocatable :: fname_ctffind_ctrl, fname_ctffind_output      
        type(oris) :: os       
        integer    :: nmovies, fromto(2), imovie, ntot, movie_counter
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced
        call read_filetable(p%filetab, movienames_forctf)
        nmovies = size(movienames_forctf)
        if( p%l_distr_exec )then
            allocate(fname_ctffind_ctrl,  source='ctffind_ctrl_file_part'//int2str_pad(p%part,p%numlen)//'.txt')
            allocate(fname_ctffind_output, source='ctffind_output_part'//int2str_pad(p%part,p%numlen)//'.txt')
            ! determine loop range
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_ctffind'
            endif
        else
            allocate(fname_ctffind_ctrl,   source='ctffind_ctrl_file.txt')
            allocate(fname_ctffind_output, source='ctffind_output.txt')
            ! determine loop range
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
        endif
        ntot = fromto(2) - fromto(1) + 1
        call os%new(ntot)
        ! loop over exposures (movies)
        movie_counter = 0
        do imovie=fromto(1),fromto(2)
            call cfiter%iterate(p, imovie, movie_counter, movienames_forctf(imovie),&
            &fname_ctffind_ctrl, fname_ctffind_output, os)
            write(*,'(f4.0,1x,a)') 100.*(real(movie_counter)/real(ntot)), 'percent of the micrographs processed'
        end do
        call os%write(fname_ctffind_output)
        call os%kill
        deallocate(fname_ctffind_ctrl,fname_ctffind_output)
        call qsys_job_finished( p, 'simple_commander_preproc :: exec_ctffind' )
        ! end gracefully
        call simple_end('**** SIMPLE_CTFFIND NORMAL STOP ****')
    end subroutine exec_ctffind

    subroutine exec_select( self, cline )
        use simple_image,    only: image
        use simple_corrmat   ! use all in there
        class(select_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: stk3_img
        type(image),           allocatable :: imgs_sel(:), imgs_all(:)
        character(len=STDLEN), allocatable :: imgnames(:)
        integer,               allocatable :: selected(:)
        real,                  allocatable :: correlations(:,:)
        logical,               allocatable :: lselected(:)
        character(len=STDLEN)              :: cmd_str
        integer                            :: iimg, isel, nall, nsel, loc(1), ios, alloc_stat
        integer                            :: funit, ldim(3), ifoo, lfoo(3), io_stat
        logical, parameter                 :: debug=.false.
        ! error check
        if( cline%defined('stk3') .or. cline%defined('filetab') )then
            ! all good
        else
            stop 'Need either stk3 or filetab are part of the command line!'
        endif
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
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
            call alloc_err('In: exec_select; simple_commander_preproc', alloc_stat)
            ! read matrix
            funit = get_fileunit()
            open(unit=funit, status='OLD', action='READ', file='corrmat_select.bin', access='STREAM')
            read(unit=funit,pos=1,iostat=io_stat) correlations
            ! Check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(exec_select): I/O error ',&
                io_stat, ' when reading corrmat_select.bin. Remove the file to override the memoization.'
                stop 'I/O error; exec_select; simple_commander_preproc'
            endif
            close(funit)
        else
            write(*,'(a)') '>>> CALCULATING CORRELATIONS'
            call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
            ! write matrix
            funit = get_fileunit()
            open(unit=funit, status='REPLACE', action='WRITE', file='corrmat_select.bin', access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) correlations
            ! Check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(exec_select): I/O error ',&
                io_stat, ' when writing to corrmat_select.bin'
                stop 'I/O error; exec_select; simple_commander_preproc'
            endif
            close(funit)
        endif
        ! find selected
        ! in addition to the index array, also make a logical array encoding the selection (to be able to reject)
        allocate(selected(nsel), lselected(nall))
        lselected = .false.
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            selected(isel) = loc(1)
            lselected(selected(isel)) = .true.
            if( debug ) print *, 'selected: ', loc(1), ' with corr: ', correlations(isel,loc(1))
        end do
        if( cline%defined('filetab') )then
            ! read filetable
            call read_filetable(p%filetab, imgnames)
            if( size(imgnames) /= nall ) stop 'nr of entries in filetab and stk not consistent'
            funit = get_fileunit()
            open(unit=funit, file=p%outfile, iostat=ios, status="replace", action="write", access="sequential")
            if( ios /= 0 )then
                write(*,*) "Error opening file name", trim(adjustl(p%outfile))
                stop
            endif
            call exec_cmdline('mkdir -p '//trim(adjustl(p%dir_select)))
            call exec_cmdline('mkdir -p '//trim(adjustl(p%dir_reject)))
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
            close(funit)
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

    subroutine exec_makepickrefs( self, cline )
        use simple_commander_volops,  only: projvol_commander
        use simple_commander_imgproc, only: stackops_commander, scale_commander
        use simple_procimgfile,       only: neg_imgfile
        class(makepickrefs_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        character(STDLEN), parameter :: ORIFILE='pickrefs_oris.txt'
        integer, parameter           :: NREFS=100
        type(params)                 :: p
        type(build)                  :: b
        type(cmdline)                :: cline_projvol, cline_stackops, cline_scale
        type(projvol_commander)      :: xprojvol
        type(stackops_commander)     :: xstackops
        integer                      :: nrots, cnt, iref, irot
        real                         :: ang, rot, rclip
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            p = params(cline, checkdistr=.false.) ! constants & derived constants produced
            if( cline%defined('vol1') )then
                call b%a%new(NREFS)
                call b%a%gen_diverse
                call b%a%write(trim(ORIFILE))
                cline_projvol = cline
                call cline_projvol%set('nspace', real(NREFS))
                call cline_projvol%set('outstk', 'pickrefs'//p%ext)
                call cline_projvol%set('oritab', trim(ORIFILE))
                call cline_projvol%set('neg', trim(p%neg))
                call cline_projvol%set('smpd', PICKER_SHRINK)
                call xprojvol%execute(cline_projvol)
            else
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
                            call b%img_copy%write('rotated_from_makepickrefs'//p%ext, cnt)
                            rot = rot + ang
                        end do
                    end do
                    call cline%set('stk', 'rotated_from_makepickrefs'//p%ext)
                endif
                if( p%neg .eq. 'yes' )then
                    call neg_imgfile('rotated_from_makepickrefs'//p%ext, 'pickrefs'//p%ext, p%smpd)
                else
                     call rename('rotated_from_makepickrefs'//p%ext, 'pickrefs'//p%ext)
                endif
            endif
        else
            stop 'need input volume (vol1) or class averages (stk) to generate picking references'
        endif
        call del_file('rotated_from_makepickrefs'//p%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_MAKEPICKREFS NORMAL STOP ****')
    end subroutine exec_makepickrefs

    subroutine exec_pick( self, cline)
        use simple_pick_iter, only: pick_iter
        class(pick_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(params)    :: p
        type(pick_iter) :: piter
        character(len=STDLEN), allocatable :: movienames_intg(:)
        character(len=STDLEN) :: boxfile
        integer :: nmovies, fromto(2), imovie, ntot, movie_counter, nptcls_out
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced
        ! check filetab existence
        call read_filetable(p%filetab, movienames_intg)
        nmovies = size(movienames_intg)
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
            call piter%iterate(cline, p, movie_counter, movienames_intg(imovie), boxfile, nptcls_out)
            write(*,'(f4.0,1x,a)') 100.*(real(movie_counter)/real(ntot)), 'percent of the micrographs processed'
        end do
    end subroutine exec_pick
    
    subroutine exec_extract( self, cline )
        use simple_nrtxtfile, only: nrtxtfile
        use simple_imgfile,   only: imgfile
        use simple_image,     only: image
        use simple_math,      only: euclid, hpsort, median
        use simple_oris,      only: oris
        use simple_stat,      only: moment
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nboxfiles, nframes, pind
        integer                            :: i, j, k, alloc_stat, ldim(3), box_current, movie, ndatlines, nptcls
        integer                            :: cnt, niter, fromto(2), orig_box
        integer                            :: movie_ind, ntot, lfoo(3), ifoo, noutside
        type(nrtxtfile)                    :: boxfile
        character(len=STDLEN)              :: mode, sumstack, outfile, sumstack_frames
        character(len=STDLEN), allocatable :: movienames(:), boxfilenames(:), movienames_frames(:)
        real, allocatable                  :: boxdata(:,:)
        integer, allocatable               :: pinds(:)
        real                               :: x, y, kv, cs, fraca, dfx, dfy, angast, ctfres
        real                               :: med, ave, sdev, var, particle_position(2)
        type(image)                        :: micrograph
        type(oris)                         :: outoris, os_uni
        logical                            :: err, params_present(3)
        logical, parameter                 :: debug = .false.
        noutside = 0
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced
        if( p%stream .eq. 'yes' )then
            if( cline%defined('outstk') )then
                ! all ok
            else
                stop 'need output stack (outstk) to be part of the command line in stream=yes mode'
            endif
            if( cline%defined('outfile') )then
                ! all ok
            else
                stop 'need output file (outfile) to be part of the command line in stream=yes mode'
            endif
        endif

        ! set output files
        if( cline%defined('outstk') )then
            if( cline%defined('dir_ptcls') )then
                sumstack = trim(p%dir_ptcls)//trim(p%outstk)
            else
                sumstack = trim(p%outstk)
            endif
        else
            if( cline%defined('dir_ptcls') )then
                sumstack = trim(p%dir_ptcls)//'sumstack'//p%ext
            else
                sumstack = 'sumstack'//p%ext
            endif
        endif
        sumstack_frames = add2fbody(sumstack, p%ext, '_frames_subset')
        if( cline%defined('outfile') )then
            if( cline%defined('dir_ptcls') )then
                outfile = trim(p%dir_ptcls)//trim(p%outfile)
            else
                outfile = trim(p%outfile)
            endif
        else
            if( cline%defined('dir_ptcls') )then
                outfile = trim(p%dir_ptcls)//'extract_params.txt'
            else
                outfile = 'extract_params.txt'
            endif
        endif

        ! check file inout existence and read filetables
        if( cline%defined('unidoc') )then
            if( .not. file_exists(p%unidoc)  ) stop 'inputted unidoc does not exist in cwd'
            nmovies = nlines(p%unidoc)
            call os_uni%new(nmovies)
            call os_uni%read(p%unidoc)
            movienames   = os_uni%extract_table('intg')
            if( os_uni%isthere('boxfile') )then
                boxfilenames = os_uni%extract_table('boxfile')
            else
                if( .not. cline%defined('boxtab')  ) stop 'need boxtab input to extract'
                nboxfiles = nlines(p%boxtab)
                if( nmovies /= nboxfiles ) stop 'number of entries in inputted files do not match!'
                call read_filetable(p%boxtab,  boxfilenames)
            endif
            if( os_uni%isthere('intg_frames') )then
                movienames_frames = os_uni%extract_table('intg_frames')
            endif
        else
            if( .not. cline%defined('filetab') ) stop 'need filetab input to extract'
            if( .not. cline%defined('boxtab')  ) stop 'need boxtab input to extract'
            if( .not. file_exists(p%filetab)   ) stop 'inputted filetab does not exist in cwd'
            if( .not. file_exists(p%boxtab)    ) stop 'inputted boxtab does not exist in cwd'
            nmovies = nlines(p%filetab)
            nboxfiles = nlines(p%boxtab)
            if( debug ) write(*,*) 'nboxfiles: ', nboxfiles
            if( nmovies /= nboxfiles ) stop 'number of entries in inputted files do not match!'
            call read_filetable(p%filetab, movienames)
            call read_filetable(p%boxtab,  boxfilenames)
        endif
        if( debug ) write(*,*) 'nmovies: ', nmovies

        ! remove possibly pre-existing output file
        call del_file(outfile)

        ! determine loop range
        fromto(1) = 1
        fromto(2) = nmovies
        ntot = fromto(2) - fromto(1) + 1

        ! count the number of particles & find ldim
        nptcls = 0
        do movie=fromto(1),fromto(2)
            if( file_exists(boxfilenames(movie)) )then
                if( nlines(boxfilenames(movie)) > 0 )then
                    call boxfile%new(boxfilenames(movie), 1)
                    ndatlines = boxfile%get_ndatalines()
                    nptcls    = nptcls + ndatlines
                    call boxfile%kill
                endif
            else
                write(*,*) 'WARNING! The inputted boxfile (below) does not exist'
                write(*,*) trim(boxfilenames(movie))
            endif
        end do
        
        ! initialize
        call find_ldim_nptcls(movienames(1), ldim, ifoo)
        call micrograph%new([ldim(1),ldim(2),1], p%smpd)
        call outoris%new(nptcls)
        pind = 0
        if( debug ) write(*,*) 'made outoris'

        ! loop over exposures (movies)
        niter    = 0
        noutside = 0
        do movie=fromto(1),fromto(2)
    
            ! get movie index (parsing the number string from the filename)
            call fname2ind(trim(adjustl(movienames(movie))), movie_ind)

            ! process boxfile
            ndatlines = 0
            if( file_exists(boxfilenames(movie)) )then
                if( nlines(boxfilenames(movie)) > 0 )then
                    call boxfile%new(boxfilenames(movie), 1)
                    ndatlines = boxfile%get_ndatalines()
                endif       
            endif
            if( ndatlines == 0 ) cycle
    
            ! update iteration counter (has to be after the cycle statements or suffer bug!!!)
            niter = niter + 1

            ! show progress
            if( niter > 1 ) call progress(niter,ntot)
    
            ! read box data
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), pinds(ndatlines), stat=alloc_stat)
            call alloc_err('In: simple_extract; boxdata etc., 2', alloc_stat)
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    stop 'Only square windows are currently allowed!'
                endif
                if( j == 1 .and. .not. cline%defined('box') ) p%box = nint(boxdata(1,3)) ! use the box size from the box file
                ! modify coordinates if change in box (shift by half the difference)
                if( orig_box /= p%box ) boxdata(j,1:2) = boxdata(j,1:2)-real(p%box-orig_box)/2.
            end do
    
            ! create particle index list and set movie index
            do j=1,ndatlines
                if( box_inside(ldim, nint(boxdata(j,1:2)), p%box) )then
                    pind     = pind + 1
                    pinds(j) = pind
                    call outoris%set(pinds(j), 'movie', real(movie_ind))
                else
                    pinds(j) = 0
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
            if( debug ) write(*,*) 'did check box parsing'
    
            ! get number of frames from stack
            call find_ldim_nptcls(movienames(movie), lfoo, nframes )
            if( nframes > 1 ) stop 'multi-frame extraction no longer supported; simple_extract'
    
            ! build general objects
            if( niter == 1 ) call b%build_general_tbox(p,cline,do3d=.false.)

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
                kv     = b%a%get(movie,'kv')
                cs     = b%a%get(movie,'cs')
                fraca  = b%a%get(movie,'fraca')
                dfx    = b%a%get(movie,'dfx')
                ctfres = b%a%get(movie,'ctfres')
                angast = 0.
                if( b%a%isthere('dfy') )then ! astigmatic CTF
                    if( .not. b%a%isthere('angast') ) stop 'need angle of astigmatism for CTF correction'
                    dfy    = b%a%get(movie,'dfy')
                    angast = b%a%get(movie,'angast')
                endif
                ! set CTF info in outoris
                do i=1,ndatlines
                    if( pinds(i) > 0 )then
                        call outoris%set(pinds(i), 'kv',         kv)
                        call outoris%set(pinds(i), 'cs',         cs)
                        call outoris%set(pinds(i), 'fraca',   fraca)
                        call outoris%set(pinds(i), 'smpd',   p%smpd)
                        call outoris%set(pinds(i), 'dfx',       dfx)
                        call outoris%set(pinds(i), 'ctfres', ctfres)
                        if( b%a%isthere('dfy') )then
                            call outoris%set(pinds(i), 'angast', angast)
                            call outoris%set(pinds(i), 'dfy',       dfy)
                        endif
                    endif
                end do
                if( debug ) write(*,*) 'did set CTF parameters dfx/dfy/angast/ctfres: ', dfx, dfy, angast, ctfres
            endif
    
            ! extract windows
            ! integrated movie
            call micrograph%read(movienames(movie),1,rwaction='READ')
            cnt = 0
            do j=1,ndatlines ! loop over boxes
                if( pinds(j) > 0 )then
                    ! extract the window
                    particle_position = boxdata(j,1:2)
                    call micrograph%window(nint(particle_position), p%box, b%img, noutside)
                    if( p%neg .eq. 'yes' ) call b%img%neg
                    call b%img%norm
                    call b%img%write(trim(adjustl(sumstack)), pinds(j))
                endif
            end do
            ! integrated movie (subset of frames)
            if( allocated(movienames_frames) )then
                call micrograph%read(movienames_frames(movie),1,rwaction='READ')
                cnt = 0
                do j=1,ndatlines ! loop over boxes
                    if( pinds(j) > 0 )then
                        ! extract the window
                        particle_position = boxdata(j,1:2)
                        call micrograph%window(nint(particle_position), p%box, b%img)
                        if( p%neg .eq. 'yes' ) call b%img%neg
                        call b%img%norm
                        call b%img%write(trim(adjustl(sumstack_frames)), pinds(j))
                    endif
                end do
            endif
            
            ! write output
            do j=1,ndatlines ! loop over boxes
                if( pinds(j) > 0 )then
                    call outoris%write(pinds(j), outfile)
                endif
            end do
    
            ! destruct
            call boxfile%kill
            deallocate(boxdata, pinds)
        end do
        if( p%outside .eq. 'yes' .and. noutside > 0 )then
            write(*,'(a,1x,i5,1x,a)') 'WARNING!', noutside, 'boxes extend outside micrograph'
        endif

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

end module simple_commander_preproc
