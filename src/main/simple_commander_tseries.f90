! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
implicit none

public :: tseries_import_commander
public :: tseries_track_commander
public :: tseries_ctf_estimate_commander
public :: tseries_split_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: tseries_import_commander
  contains
    procedure :: execute      => exec_tseries_import
end type tseries_import_commander
type, extends(commander_base) :: tseries_track_commander
  contains
    procedure :: execute      => exec_tseries_track
end type tseries_track_commander
type, extends(commander_base) :: tseries_backgr_subtr_commander
  contains
    procedure :: execute      => exec_tseries_backgr_subtr
end type tseries_backgr_subtr_commander
type, extends(commander_base) :: tseries_ctf_estimate_commander
  contains
    procedure :: execute      => exec_tseries_ctf_estimate
end type tseries_ctf_estimate_commander
type, extends(commander_base) :: tseries_split_commander
  contains
    procedure :: execute      => exec_tseries_split
end type tseries_split_commander

contains

    !> for creating a new project
    subroutine exec_tseries_import( self, cline )
        use simple_sp_project,        only: sp_project
        use simple_image,             only: image
        class(tseries_import_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)            :: params
        type(sp_project)            :: spproj
        type(image)                 :: frame_img
        type(ctfparams)             :: ctfvars
        character(len=LONGSTRLEN), allocatable :: filenames(:), movfnames(:)
        character(len=LONGSTRLEN)              :: outfname
        integer          :: ldim(3), nframes, frame_from, frame_to, numlen, cnt
        integer          :: iframe, jframe, nfiles, endit
        call cline%set('oritype','mic')
        ! check input
        if( cline%defined('nframesgrp') )then
            if( nint(cline%get_rarg('nframesgrp')) < 3 )then
                THROW_HARD('nframesgrp integer (nr of frames to average) needs to be >= 3; exec_tseries_extract')
            endif
        else
            THROW_HARD('need nframesgrp integer input = nr of frames to average; exec_tseries_extract')
        endif
        ! build parameters
        call params%new(cline)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        numlen = len(int2str(nfiles))
        if( params%mkdir.eq.'yes' )then
            do iframe = 1,nfiles
                if(filenames(iframe)(1:1).ne.'/') filenames(iframe) = '../'//trim(filenames(iframe))
            enddo
        endif
        call find_ldim_nptcls(filenames(1),ldim,nframes)
        if( nframes == 1 .and. ldim(3) == 1 )then
            ! all ok
            call frame_img%new(ldim, params%smpd)
        else
            write(logfhandle,*) 'ldim(3): ', ldim(3)
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('exec_tseries_import assumes one frame per file')
        endif
        ! generates 'movies' to import
        call spproj%read(params%projfile)
        endit = nfiles - params%nframesgrp + 1
        allocate(movfnames(endit))
        do iframe=1,endit
            call progress(iframe, endit)
            if( cline%defined('fbody') )then
                outfname = 'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            else
                outfname = trim(params%fbody)//'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            endif
            movfnames(iframe) = trim(outfname)
            frame_from = iframe
            frame_to   = iframe + params%nframesgrp - 1
            cnt = 0
            do jframe=frame_from,frame_to
                cnt = cnt + 1
                call frame_img%read(filenames(jframe),1)
                call frame_img%write(movfnames(iframe),cnt)
            end do
        end do
        call frame_img%kill
        ! import
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = params%smpd
        ctfvars%cs      = params%cs
        ctfvars%kv      = params%kv
        ctfvars%fraca   = params%fraca
        ctfvars%l_phaseplate = .false.
        ctfvars%dfx     = 0.
        ctfvars%dfy     = 0.
        ctfvars%angast  = 0.
        ctfvars%phshift = 0.
        call spproj%add_movies( movfnames, ctfvars )
        !final write
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT NORMAL STOP ****')
    end subroutine exec_tseries_import

    subroutine exec_tseries_track( self, cline )
        use simple_tseries_tracker
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_sp_project, only: sp_project
        class(tseries_track_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(tseries_backgr_subtr_commander)   :: backgr_subtracter
        type(cmdline)                          :: cline_backgrsubtr
        type(sp_project)                       :: spproj
        type(parameters)                       :: params
        type(ctfparams)                        :: ctfvars
        type(nrtxtfile)                        :: boxfile
        type(oris)                             :: deftab
        character(len=:),          allocatable :: dir, stk, outstk, ext, stkneigh
        character(len=LONGSTRLEN), allocatable :: framenames(:)
        real,                      allocatable :: boxdata(:,:)
        integer :: i, iframe, ndatlines, j, orig_box, numlen, nmovs, ineigh
        call params%new(cline)
        numlen   = 5 ! default value
        orig_box = params%box
        ! boxfile input
        if( cline%defined('boxfile') )then
            if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
            if( nlines(params%boxfile) > 0 )then
                call boxfile%new(params%boxfile, 1)
                ndatlines = boxfile%get_ndatalines()
                numlen    = len(int2str(ndatlines))
                allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track boxdata')
                do j=1,ndatlines
                    call boxfile%readNextDataLine(boxdata(j,:))
                    orig_box = nint(boxdata(j,3))
                    if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )THROW_HARD('Only square windows allowed!')
                end do
            else
                THROW_HARD('inputted boxfile is empty; exec_tseries_track')
            endif
        else if( cline%defined('xcoord') .and. cline%defined('ycoord') )then
            if( .not. cline%defined('box') ) THROW_HARD('need box to be part of command line for this mode of execution; exec_tseries_track')
            allocate( boxdata(1,2) )
            boxdata(1,1) = real(params%xcoord)
            boxdata(1,2) = real(params%ycoord)
            ndatlines = 1
        else
            THROW_HARD('need either boxfile or xcoord/ycoord to be part of command line; exec_tseries_track')
        endif
        ! frames input
        call spproj%read(params%projfile)
        nmovs = spproj%get_nintgs()
        allocate(framenames(nmovs))
        iframe = 0
        do i = 1,spproj%os_mic%get_noris()
            if( spproj%os_mic%isthere(i,'intg') )then
                iframe = iframe + 1
                framenames(iframe) = trim(spproj%os_mic%get_static(i,'intg'))
            endif
        enddo
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        cline_backgrsubtr = cline
        call cline_backgrsubtr%set('projfile',params%projfile)
        if( cline%defined('lp_backgr') )then
            call cline_backgrsubtr%set('lp', params%lp_backgr)
        else
            call cline_backgrsubtr%set('lp', 1.1)
        endif
        if( .not.cline%defined('hp')    ) call cline_backgrsubtr%set('hp',5.0)
        if( .not.cline%defined('width') ) call cline_backgrsubtr%set('width',5.0)
        if( .not.cline%defined('ctf')   ) call cline_backgrsubtr%set('ctf', 'no')
        ! actual tracking
        do j=1,ndatlines
            ! extract stack
            call init_tracker( nint(boxdata(j,1:2)), framenames)
            call track_particle
            if( cline%defined('ind') )then
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
                if( .not. cline%defined('numlen') ) THROW_HARD('need numlen to be part of command line if ind is; exec_tseries_track')
                call write_tracked_series(dir,trim(params%fbody)//int2str_pad(params%ind,params%numlen))
            else
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
                call write_tracked_series(dir,trim(params%fbody)//int2str_pad(j,numlen))
            endif
            ! perform background substraction
            stk = trim(tracker_get_stkname())
            call cline_backgrsubtr%set('stk', stk)
            do ineigh=1,tracker_get_nnn()
                stkneigh = trim(tracker_get_nnfname(ineigh))
                call cline_backgrsubtr%set('stk_backgr', stkneigh)
                ext    = fname2ext(stkneigh)
                outstk = trim(get_fbody(stkneigh,ext,.true.))//'_stk.'//trim(ext)
                call cline_backgrsubtr%set('outstk', outstk)
                call backgr_subtracter%execute( cline_backgrsubtr )
            enddo
            ! clean
            call kill_tracker
        end do
        call qsys_job_finished(  'simple_commander_tseries :: exec_tseries_track')
        ! updates project
        if( spproj%os_mic%isthere('dfx') )then
            if( params%ctf.eq.'flip' )then
                ! the stacks have been phase-flipped in backgr_subtracter
                do i = 1,nmovs
                    call spproj%os_mic%set(i,'ctf','flip')
                enddo
                call spproj%write
            endif
            ! write defocus for convenience
            call deftab%new(nmovs)
            do i = 1,nmovs
                ctfvars = spproj%os_mic%get_ctfvars(i)
                call deftab%set_ctfvars(i,ctfvars)
            enddo
            call deftab%write('deftab.txt')
            call deftab%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track

    subroutine exec_tseries_backgr_subtr( self, cline )
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        use simple_image, only: image
        use simple_ctf,   only: ctf
        class(tseries_backgr_subtr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: img_backgr, img_backgr_wctf, ave_img
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfvars
        character(len=:), allocatable :: ext,imgname
        real             :: ave,sdev,minv,maxv
        integer          :: iptcl, nptcls, ldim(3)
        logical          :: do_flip, err
        call cline%set('oritype','mic')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! dimensions
        call find_ldim_nptcls(params%stk,ldim,nptcls)
        ! CTF logics
        do_flip = .false.
        if( (build%spproj_field%isthere('ctf')) .and. (cline%defined('ctf').and.params%ctf.eq.'flip') )then
            if( build%spproj%get_nintgs() /= nptcls ) THROW_HARD('Incompatible # of images and micrographs!')
            if( .not.build%spproj_field%isthere('dfx') ) THROW_HARD('Missing CTF parameters')
            do_flip = .true.
        endif
        ! get background image
        call img_backgr%new([params%box,params%box,1], params%smpd)
        call img_backgr_wctf%new([params%box,params%box,1], params%smpd)
        call img_backgr%read(params%stk_backgr, 1)
        ! background image is skipped if sdev is small because this neighbour was
        ! outside the micrograph
        call img_backgr%stats(ave, sdev, maxv, minv, errout=err)
        if( sdev>TINY .and. .not.err )then
            call ave_img%new([params%box,params%box,1], params%smpd)
            ave_img = 0.
            ! main loop
            do iptcl=1,nptcls
                call progress(iptcl,nptcls)
                ! read particle image
                call build%img%read(params%stk, iptcl)
                ! copy background image & CTF
                img_backgr_wctf = img_backgr
                ! fwd ft
                call build%img%fft()
                call img_backgr_wctf%fft()
                if( do_flip )then
                    ctfvars = build%spproj_field%get_ctfvars(iptcl)
                    tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
                    call tfun%apply_serial(build%img, 'flip', ctfvars)
                    call tfun%apply_serial(img_backgr_wctf, 'flip', ctfvars)
                endif
                ! filter background image
                call img_backgr_wctf%bp(params%hp,params%lp,width=params%width)
                ! subtract background
                call build%img%subtr(img_backgr_wctf)
                ! bwd ft
                call build%img%ifft()
                ! normalise
                call build%img%norm()
                call ave_img%add(build%img)
                ! output corrected image
                call build%img%write(params%outstk, iptcl)
            end do
            ! generates average
            call ave_img%div(real(nptcls))
            ext     = fname2ext(params%outstk)
            imgname = trim(get_fbody(params%outstk,ext,.true.))//'_ave.'//trim(ext)
            call ave_img%write(imgname)
        endif
        call ave_img%kill
        call img_backgr%kill
        call img_backgr_wctf%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_BACKGR_SUBTR NORMAL STOP ****')
    end subroutine exec_tseries_backgr_subtr

    subroutine exec_tseries_ctf_estimate( self, cline )
        use simple_image,             only: image
        use simple_ctf_estimate
        class(tseries_ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        character(len=LONGSTRLEN), parameter :: pspec_fname  = 'tseries_ctf_estimate_pspec.mrc'
        character(len=LONGSTRLEN), parameter :: diag_fname   = 'tseries_ctf_estimate_diag'//JPG_EXT
        integer,                   parameter :: nmics4ctf    = 10
        type(parameters)              :: params
        type(builder)                 :: build
        type(image)                   :: mic_img, pspec_avg, pspec_up_avg, pspec_low_avg
        type(image)                   :: pspec, pspec_up, pspec_low
        type(ctfparams)               :: ctfvars
        character(len=:), allocatable :: micname
        real    :: cc90, dfx, dfy, angast, phshift, dferr, cc, ctfscore
        integer :: nmics, imic, dims(3), nselmics
        call cline%set('oritype','mic')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! generates ave
        dims(1) = nint(build%spproj_field%get(1,'xdim'))
        dims(2) = nint(build%spproj_field%get(1,'ydim'))
        dims(3) = 1
        call mic_img%new(dims, params%smpd)
        call pspec_low%new([params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec_up%new( [params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec%new(    [params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec_low_avg%new([params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec_up_avg%new( [params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec_avg%new(    [params%pspecsz,params%pspecsz,1], params%smpd)
        nmics = build%spproj%get_nintgs()
        if( nmics /= build%spproj%os_mic%get_noris() ) THROW_HARD('Invalid number of micrographs')
        write(logfhandle,'(A)')'>>> BUILDING AVERAGE POWER SPECTRUM'
        ! cumulative sums of power
        nselmics = 0
        do imic = 1,nmics,10
            call progress(imic,nmics)
            nselmics = nselmics + 1
            call build%spproj_field%getter(imic,'forctf',micname)
            call mic_img%read(micname)
            pspec_low = 0.
            pspec_up  = 0.
            pspec     = 0.
            call mic_img%mic2eospecs(params%pspecsz, 'power', params%hp, pspec_low, pspec_up, pspec, postproc=.false.)
            call pspec_low_avg%add(pspec_low,w=0.01)
            call pspec_up_avg%add(pspec_up,  w=0.01)
            call pspec_avg%add(pspec,        w=0.01)
        enddo
        call progress(1,1)
        call mic_img%kill
        call pspec%kill
        call pspec_up%kill
        call pspec_low%kill
        ! averages powers
        call pspec_low_avg%div(0.01*real(nselmics))
        call pspec_up_avg%div( 0.01*real(nselmics))
        call pspec_avg%div(    0.01*real(nselmics))
        ! square root
        call pspec_low_avg%sq_rt
        call pspec_up_avg%sq_rt
        call pspec_avg%sq_rt
        ! post-process
        call pspec_low_avg%dampen_pspec_central_cross
        call pspec_up_avg%dampen_pspec_central_cross
        call pspec_avg%dampen_pspec_central_cross
        call pspec_low_avg%subtr_backgr(params%hp)
        call pspec_up_avg%subtr_backgr(params%hp)
        call pspec_avg%subtr_backgr(params%hp)
        ! ctf estimation
        write(logfhandle,'(A)')'>>> CTF ESTIMATION'
        ctfvars = build%spproj_field%get_ctfvars(1)
        call ctf_estimate_init(pspec_avg, pspec_low_avg, pspec_up_avg, params%smpd, ctfvars%kv,&
            &ctfvars%cs, ctfvars%fraca, [params%dfmin,params%dfmax], [params%hp,params%lp],&
            &params%astigtol, ctfvars%l_phaseplate, cc90)
        call ctf_estimate_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, ctfscore, diag_fname)
        call ctf_estimate_kill
        ctfvars%dfx     = dfx
        ctfvars%dfy     = dfy
        ctfvars%angast  = angast
        ctfvars%phshift = phshift
        do imic = 1,nmics
            call build%spproj_field%set_ctfvars(imic, ctfvars)
        enddo
        ! final write
        call build%spproj%write(params%projfile)
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_tseries_ctf_estimate

    subroutine exec_tseries_split( self, cline )
        class(tseries_split_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(oris)       :: os
        character(len=:), allocatable :: stkname, oriname
        integer :: i, iptcl, numlen, istart, istop, cnt, cnt2, ntot
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        numlen = len(int2str(params%nptcls))
        ! count the number of chunks
        ntot = 0
        do i=1,params%nptcls,params%stepsz
            istart = i
            istop  = i + params%chunksz - 1
            if( istop > params%nptcls ) exit
            ntot = ntot + 1
        end do
        ! split
        do i=1,params%nptcls,params%stepsz
            istart = i
            istop  = i + params%chunksz - 1
            if( istop > params%nptcls ) exit
            cnt    = cnt + 1
            call progress(cnt,ntot)
            call os%new(params%chunksz)
            call simple_mkdir('tseries_chunk'//int2str_pad(cnt,numlen),errmsg="commander_tseries::exec_tseries_split")
            stkname = filepath(PATH_HERE,'tseries_chunk'//int2str_pad(cnt,numlen),'imgs'//params%ext)
            oriname = filepath('tseries_chunk'//int2str_pad(cnt,numlen),'oris'//trim(TXT_EXT))
            call del_file( stkname )
            call del_file( oriname )
            cnt2 = 0
            do iptcl=istart,istop
                cnt2 = cnt2 + 1
                call os%set_ori(cnt2, build%spproj_field%get_ori(iptcl) )
                call build%img%read(params%stk, iptcl)
                call build%img%write(stkname,cnt2)
            end do
            call os%write(oriname)
            call os%kill
            deallocate(stkname, oriname)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_SPLIT NORMAL STOP ****')
    end subroutine exec_tseries_split

end module simple_commander_tseries
