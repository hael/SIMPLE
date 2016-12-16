!==Class simple_commander_distr
!
! This class contains the set of concrete distr commanders of the SIMPLE library used to provide pre/post processing routines
! for SIMPLE when executed in distributed mode. This class provides the glue between the reciver (main reciever is simple_exec 
! program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Frederic Bonnet, Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_distr
use simple_defs            ! singleton
use simple_jiffys          ! singleton
use simple_timing          ! singleton
use simple_cuda            ! singleton
use simple_cuda_defs       ! singleton
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: merge_algndocs_commander
public :: merge_crefine_out_commander
public :: merge_nnmat_commander
public :: merge_shellweights_commander
public :: merge_similarities_commander
public :: split_pairs_commander
public :: split_commander
private

type, extends(commander_base) :: merge_algndocs_commander
  contains
    procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander
type, extends(commander_base) :: merge_crefine_out_commander
  contains
    procedure :: execute      => exec_merge_crefine_out
end type merge_crefine_out_commander
type, extends(commander_base) :: merge_nnmat_commander
  contains
    procedure :: execute      => exec_merge_nnmat
end type merge_nnmat_commander
type, extends(commander_base) :: merge_shellweights_commander
  contains
    procedure :: execute      => exec_merge_shellweights
end type merge_shellweights_commander
type, extends(commander_base) :: merge_similarities_commander
  contains
    procedure :: execute      => exec_merge_similarities
end type merge_similarities_commander
type, extends(commander_base) :: split_pairs_commander
  contains
    procedure :: execute      => exec_split_pairs
end type split_pairs_commander
type, extends(commander_base) :: split_commander
  contains
    procedure :: execute      => exec_split
end type split_commander

contains
    
    subroutine exec_merge_algndocs( self, cline )
        use simple_oris, only: oris
        use simple_map_reduce  ! singleton
        class(merge_algndocs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)          :: p
        type(oris)            :: o, o_read
        integer               :: i, j, nentries, cnt, nentries_all, numlen
        integer, allocatable  :: parts(:,:)
        character(len=STDLEN) :: fname
        logical               :: here, useoritab
        p = params(cline) ! parameters generated
        useoritab = .false.
        if( cline%defined('oritab') )then
            if( file_exists(p%oritab) ) useoritab = .true.
        endif
        if( useoritab )then
            if( nlines(p%oritab) /= p%nptcls )then
                stop 'the inputted nptcls is not consistent with the nptcls in oritab!'
            endif
            ! create object for orientations
            o = oris(p%nptcls)
            ! read previous orientations
            call o%read(p%oritab)
        endif
        select case(p%split_mode)
            case('chunk')
                parts = split_nobjs_in_chunks(p%nptcls, p%chunksz)
            case('even')
                parts = split_nobjs_even(p%nptcls, p%ndocs)
            case DEFAULT
                write(*,*) 'split mode: ', trim(p%split_mode)
                stop 'unsupported split_mode; simple_commander_distr :: exec_merge_algndocs'
        end select
        numlen = len(int2str(p%ndocs))
        do i=1,p%ndocs
            fname = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
            ! calculate the number of all entries
            nentries_all = parts(i,2)-parts(i,1)+1 
            ! calculate the actual number of entries
            inquire(FILE=fname, EXIST=here)
            if( here )then
                nentries = nlines(fname)
            else
                nentries = 0
            endif
            ! check if oritab is there to fill-in blanks
            if( nentries < nentries_all )then
                write(*,*) 'nentries: ',     nentries
                write(*,*) 'nentries_all: ', nentries_all
                if( .not. useoritab )then
                    stop 'need previous oritab to fill-in blanks; simple_merge_algndocs'
                endif
            endif
            ! print partition info
            write(*,'(a,1x,i3,1x,a,1x,i6,1x,i6)') 'partition:', i, 'from/to:', parts(i,1), parts(i,2)
            if( nentries > 0 )then
                o_read = oris(nentries)
                call o_read%read(fname)
            endif
            ! read
            if( useoritab )then ! read and fill-in from oritab
                cnt = 0
                do j= parts(i,1),parts(i,2)
                    cnt = cnt+1
                    if( cnt <= nentries )then
                        call o%set_ori(j,o_read%get_ori(cnt))
                    else
                        exit
                    endif
                end do
            else ! just merge (all ptcls is there)
                call o%merge(o_read)
            endif
        end do
        call o%write(p%outfile)
        deallocate(parts)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****')
    end subroutine exec_merge_algndocs

    subroutine exec_merge_crefine_out( self, cline )
        use simple_oris,       only: oris
        use simple_ori,        only: ori
        use simple_math,       only: calc_lowpass_lim
        use simple_online_var, only: online_var
        class(merge_crefine_out_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(params)                  :: p
        type(oris)                    :: o, o_read
        type(ori)                     :: osingle
        type(online_var), allocatable :: resvars(:)
        integer, allocatable          :: pinds(:)
        real,    allocatable          :: shweights(:,:)
        integer                       :: icls, nentries, ncls, ipind, find, filtsz, iimg
        character(len=STDLEN)         :: fname_doc, fname_shweights
        integer, parameter            :: NUMLEN=5
        real                          :: mv(2), sdev
        p = params(cline) ! parameters generated
        if( nlines(p%oritab) /= p%nptcls )then
            stop 'the inputted nptcls is not consistent with the nptcls in oritab!'
        endif
        ! create object for orientations
        o = oris(p%nptcls)
        ! read previous orientations
        call o%read(p%oritab)
        ! get number of classes 
        ncls = o%get_ncls()
        ! loop over classes
        do icls=1,ncls
            fname_doc  = 'classrefine_doc_class'//int2str_pad(icls,NUMLEN)//'.txt'
            if( file_exists(fname_doc) )then
                nentries = nlines(fname_doc)
                ! extract the particle indices of icls
                pinds = o%get_cls_pinds(icls)
                if( size(pinds) .eq. nentries )then
                    ! replace the corresponding oris in the mother doc
                    call o_read%new(nentries)
                    call o_read%read(fname_doc)
                    do ipind=1,nentries
                        osingle = o_read%get_ori(ipind)
                        call o%set_ori(pinds(ipind),osingle)
                    end do
                    call del_txtfile(fname_doc)
                else
                    write(*,*) 'doc: ', trim(fname_doc)
                    stop 'nr of entries in classdoc not equal to class population in mother doc, aborting'
                endif
            endif
            fname_shweights = 'classrefine_shweights_class'//int2str_pad(icls,NUMLEN)//'.bin'
            if( file_exists(fname_shweights) )then
                shweights = file2arr2D(fname_shweights)
                filtsz    = size(shweights,2)
                if( .not. allocated(resvars) ) allocate(resvars(filtsz))
                do find=1,filtsz
                    do iimg=1,size(shweights,1)
                        call resvars(find)%add(shweights(iimg,find))
                    end do
                end do
                call del_binfile(fname_shweights)
            endif
        end do
        ! finalize and print the weight variances (per-resolution-shell)
        do find=1,filtsz
            call resvars(find)%finalize
            sdev = 0.
            mv = resvars(find)%get_mean_var()
            if( mv(2) > TINY ) sdev = sqrt(mv(2))
            write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:',&
            calc_lowpass_lim(find, p%box, p%smpd), '>>> SDEV(WEIGHT):', sdev
        end do
        call o%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_CLASSDOCS NORMAL STOP ****')
    end subroutine exec_merge_crefine_out

    subroutine exec_merge_nnmat( self, cline )
        use simple_map_reduce, only: merge_nnmat_from_parts
        class(merge_nnmat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)         :: p
        integer, allocatable :: nnmat(:,:)
        integer :: filnum, io_stat
        p      = params(cline) ! parameters generated
        nnmat  = merge_nnmat_from_parts(p%nptcls, p%nparts, p%nnn)
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='nnmat.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) nnmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to nnmat.bin'
            stop 'I/O error; simple_merge_nnmat'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_NNMAT NORMAL STOP ****')
    end subroutine exec_merge_nnmat

    subroutine exec_merge_shellweights( self, cline )
        use simple_map_reduce, only: merge_rmat_from_parts
        use simple_stat,       only: corrs2weights, moment, normalize_sigm
        class(merge_shellweights_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: wmat(:,:),  weights_tmp(:), wsums(:)
        integer :: filnum, io_stat, filtsz, ishell, nptcls, alloc_stat, iptcl
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built (assumes stk input)
        filtsz = b%img%get_filtsz()         ! nr of resolution elements
        wmat   = merge_rmat_from_parts(p%nptcls, p%nparts, filtsz, 'shellweights_part')
        write(*,'(A)') '>>> CALCULATING SHELL-BY-SHELL WEIGHTS'
        ! create weights, shell-by-shell
        nptcls = size(wmat,1)
        do ishell=1,filtsz
            weights_tmp = corrs2weights(wmat(:,ishell))*real(nptcls)
            wmat(:,ishell) = weights_tmp
            deallocate(weights_tmp)
        end do
        ! set and write the per-particle wsums
        allocate(wsums(nptcls), stat=alloc_stat)
        call alloc_err('In: simple_commander_distr :: exec_merge_shellweights', alloc_stat)
        do iptcl=1,nptcls
            wsums(iptcl) = sum(wmat(iptcl,:))               
        end do
        call normalize_sigm(wsums)
        do iptcl=1,nptcls
            wmat(iptcl,:) = wmat(iptcl,:)*wsums(iptcl)
        end do
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='shellweights.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) wmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to shellweights.bin'
            stop 'I/O error; merge_shellweights'
        endif
        close(filnum)
        deallocate(wsums,wmat)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SHELLWEIGHTS NORMAL STOP ****')
    end subroutine exec_merge_shellweights
    
    subroutine exec_merge_similarities( self, cline )
        use simple_map_reduce, only: merge_similarities_from_parts
        class(merge_similarities_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        real, allocatable :: simmat(:,:)
        integer           :: filnum, io_stat
        p      = params(cline) ! parameters generated
        simmat = merge_similarities_from_parts(p%nptcls, p%nparts)
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            stop 'I/O error; simple_merge_similarities'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****')
    end subroutine exec_merge_similarities

    subroutine exec_split_pairs( self, cline )
        use simple_map_reduce, only: split_pairs_in_parts
        class(split_pairs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call split_pairs_in_parts(p%nptcls, p%nparts)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****')
    end subroutine exec_split_pairs

    subroutine exec_split( self, cline )
        use simple_map_reduce   ! singleton
        use simple_image, only: image
        class(split_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params)         :: p
        type(image)          :: img
        integer              :: iptcl, ipart, ldim(3), cnt, nimgs
        character(len=4)     :: ext
        integer, allocatable :: parts(:,:)
        logical              :: either_defined
        p = params(cline) ! parameters generated
        ext = '.'//fname2ext(p%stk)
        call find_ldim_nptcls(p%stk, ldim, nimgs)
        ldim(3) = 1
        call img%new(ldim,1.)
        select case(p%split_mode)
            case('chunk')
                parts = split_nobjs_in_chunks(nimgs, p%chunksz)
            case('even')
                parts = split_nobjs_even(nimgs, p%nparts)
            case DEFAULT
                write(*,*) 'split mode: ', trim(p%split_mode)
                stop 'unsupported split_mode; simple_commander_distr :: exec_split'
        end select
        if( size(parts,1) /= p%nparts ) stop 'ERROR! generated number of parts not same as inputted nparts'
        do ipart=1,p%nparts
            call progress(ipart,p%nparts)
            cnt = 0
            do iptcl=parts(ipart,1),parts(ipart,2)
                cnt = cnt+1
                call img%read(p%stk, iptcl)
                if( p%neg .eq. 'yes' ) call img%neg
                call img%write('stack_part'//int2str_pad(ipart,p%numlen)//ext, cnt)
            end do
        end do
        deallocate(parts)
        call img%kill
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****')
    end subroutine exec_split

end module simple_commander_distr
