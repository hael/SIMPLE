!==Class simple_commander_comlin
!
! This class contains the set of concrete common-lines commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Frederic Bonnet, Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_comlin
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

public :: comlin_smat_commander
public :: symsrch_commander
private

type, extends(commander_base) :: comlin_smat_commander
  contains
    procedure :: execute      => exec_comlin_smat
end type comlin_smat_commander
type, extends(commander_base) :: symsrch_commander
  contains
    procedure :: execute      => exec_symsrch
end type symsrch_commander

contains
    
    subroutine exec_comlin_smat( self, cline )
        use simple_comlin_sym  ! singleton
        use simple_comlin_corr ! singleton
        use simple_ori,        only: ori
        use simple_imgfile,    only: imgfile
        use simple_comlin,     only: comlin
        class(comlin_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target          :: p
        type(build), target           :: b
        type(ori)                     :: orientation_best
        integer                       :: iptcl, jptcl, alloc_stat, funit, io_stat
        integer                       :: cnt, ntot, npairs, ipair, fnr
        real, allocatable             :: corrmat(:,:), corrs(:)
        integer, allocatable          :: pairs(:,:)
        logical                       :: debug=.false.
        character(len=:), allocatable :: fname
        p = params(cline, .false.)                           ! constants & derived constants produced
        call b%build_general_tbox(p, cline, .false., .true.) ! general objects built (no oritab reading)
        allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
        call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
        if( debug ) print *, 'analysing this number of objects: ', p%nptcls
        do iptcl=1,p%nptcls
            call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd, p%imgkind)
            call b%imgs_sym(iptcl)%read(p%stk, iptcl)
            ! apply a soft-edged mask
            call b%imgs_sym(iptcl)%mask(p%msk, 'soft')
            ! Fourier transform
            call b%imgs_sym(iptcl)%fwd_ft
        end do
        b%clins = comlin(b%a, b%imgs_sym)
        if( cline%defined('part') )then
            npairs = p%top-p%fromp+1
            if( debug ) print *, 'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
            ! read the pairs
            funit = get_fileunit()
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
            if( .not. file_exists(fname) )then
                write(*,*) 'file: ', fname, ' does not exist!'
                write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
                stop 'I/O error; simple_comlin_smat'
            endif
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            if( debug ) print *, 'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_comlin_smat): I/O error ', io_stat, ' when reading file: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname)
            ! calculate the similarities
            call comlin_sym_init(b, p)
            cnt = 0
            do ipair=p%fromp,p%top
                cnt = cnt+1
                call progress(cnt, npairs)
                p%iptcl = pairs(ipair,1)
                p%jptcl = pairs(ipair,2)
                call comlin_sym_axis(orientation_best, 'pair', .false.)
                corrs(ipair) = orientation_best%get('corr')
            end do
            if( debug ) print *, 'did set this number of similarities: ', cnt
            ! write the similarities
            funit = get_fileunit()
            allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin')
            open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
            ! Check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_comlin_smat): I/O error ', io_stat, ' when writing to: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname, corrs, pairs)
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen), STATUS='REPLACE', action='WRITE', iostat=io_stat)
            call fopen_err( 'In: simple_comlin_smat', io_stat )
            write(fnr,'(A)') '**** SIMPLE_COMLIN_SMAT NORMAL STOP ****'
            close(fnr)
        else
            allocate(corrmat(p%nptcls,p%nptcls), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 3', alloc_stat)
            corrmat = 1.
            ntot = (p%nptcls*(p%nptcls-1))/2
            call comlin_sym_init(b, p)
            cnt = 0
            do iptcl=1,p%nptcls-1
                do jptcl=iptcl+1,p%nptcls
                    cnt = cnt+1
                    call progress(cnt, ntot)
                    p%iptcl = iptcl
                    p%jptcl = jptcl
                    call comlin_sym_axis(orientation_best, 'pair', .false.)
                    corrmat(iptcl,jptcl) = orientation_best%get('corr')
                    corrmat(jptcl,iptcl) = corrmat(iptcl,jptcl)
                end do
            end do
            call progress(ntot, ntot)
            funit = get_fileunit()
            open(unit=funit, status='REPLACE', action='WRITE', file='clin_smat.bin', access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to clin_smat.bin'
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(corrmat)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_COMLIN_SMAT NORMAL STOP ****')
    end subroutine exec_comlin_smat
    
    subroutine exec_symsrch( self, cline )
        use simple_oris,      only: oris
        use simple_symsrcher, only: symsrch_master
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(oris)   :: o
        p = params(cline) ! parameters generated
        if( cline%defined('stk') )then
            p%nptcls = 1
            p%nspace = 1
        endif
        call b%build_general_tbox(p, cline, .true., .true.) ! general objects built (no oritab reading)
        call b%build_comlin_tbox(p)                         ! objects for common lines based alignment built
        call symsrch_master( cline, p, b, o )
        call o%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
    end subroutine exec_symsrch

end module simple_commander_comlin
