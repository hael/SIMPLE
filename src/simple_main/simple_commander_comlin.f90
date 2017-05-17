!==Class simple_commander_comlin
!
! This class contains the set of concrete common-lines commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_comlin
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
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
        use simple_comlin_srch   ! use all in there
        use simple_ori,          only: ori
        use simple_imgfile,      only: imgfile
        use simple_comlin,       only: comlin
        use simple_qsys_funs,    only: qsys_job_finished
        use simple_strings,      only: int2str_pad
        class(comlin_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target          :: p
        type(build),  target          :: b
        type(ori)                     :: orientation_best
        integer                       :: iptcl, jptcl, alloc_stat, funit, io_stat
        integer                       :: ntot, npairs, ipair, fnr
        real,    allocatable          :: corrmat(:,:), corrs(:)
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
        b%clins = comlin(b%a, b%imgs_sym, p%lp)
        if( cline%defined('part') )then
            npairs = p%top - p%fromp + 1
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
            call comlin_srch_init( b, p, 'simplex', 'pair' )
            do ipair=p%fromp,p%top
                p%iptcl = pairs(ipair,1)
                p%jptcl = pairs(ipair,2)
                corrs(ipair) = comlin_srch_pair()
            end do
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
            call qsys_job_finished(p,'simple_commander_comlin :: exec_comlin_smat')
        else
            allocate(corrmat(p%nptcls,p%nptcls), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 3', alloc_stat)
            corrmat = 1.
            ntot = (p%nptcls*(p%nptcls-1))/2
            ! calculate the similarities
            call comlin_srch_init( b, p, 'simplex', 'pair' )
            do iptcl=1,p%nptcls-1
                do jptcl=iptcl+1,p%nptcls
                    p%iptcl = iptcl
                    p%jptcl = jptcl
                    corrmat(iptcl,jptcl) = comlin_srch_pair()
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
        use simple_strings,        only: int2str_pad
        use simple_oris,           only: oris
        use simple_ori,            only: ori
        use simple_projector_hlev, only: projvol
        use simple_comlin_srch     ! use all in there
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        type(ori)                     :: symaxis_ori
        type(oris)                    :: os, oshift
        integer                       :: fnr, file_stat, comlin_srch_nbest, cnt, i, j
        real                          :: shvec(3)
        character(len=STDLEN)         :: fname_finished
        character(len=32), parameter  :: SYMSHTAB = 'sym_3dshift.txt'
        p = params(cline)                                   ! parameters generated
        call b%build_general_tbox(p, cline, .true., .true.) ! general objects built (no oritab reading)
        call b%build_comlin_tbox(p)                         ! objects for common lines based alignment built
        ! center volume
        call b%vol%read(p%vols(1))
        shvec = b%vol%center(p%cenlp,'no',p%msk)
        if( p%l_distr_exec .and. p%part.eq.1 )then
            ! writes shifts for distributed execution
            call oshift%new(1)
            call oshift%set(1,'x',shvec(1))
            call oshift%set(1,'y',shvec(2))
            call oshift%set(1,'z',shvec(3))
            call oshift%write(trim(SYMSHTAB))
            call oshift%kill
        endif
        ! generate projections
        call b%vol%mask(p%msk, 'soft')
        b%ref_imgs(1,:) = projvol(b%vol, b%e, p)
        ! expand over symmetry group
        cnt = 0
        do i=1,p%nptcls
            do j=1,p%nsym
                cnt = cnt+1
                b%imgs_sym(cnt) = b%ref_imgs(1,i)
                call b%imgs_sym(cnt)%fwd_ft
            end do
        end do
        ! search for the axes
        call comlin_srch_init( b, p, 'simplex', 'sym')
        if(  p%l_distr_exec )then
            call comlin_srch_symaxis( symaxis_ori, [p%fromp,p%top])
            call comlin_srch_write_resoris('symaxes_part'//int2str_pad(p%part,p%numlen)//'.txt', [p%fromp,p%top])
        else
            comlin_srch_nbest = comlin_srch_get_nbest()
            call comlin_srch_symaxis( symaxis_ori )
            call comlin_srch_write_resoris('symaxes.txt', [1,comlin_srch_nbest])
            if( cline%defined('oritab') )then
                ! rotate the orientations & transfer the 3d shifts to 2d
                shvec = -1.*shvec
                os    = oris(nlines(p%oritab))
                call os%read(p%oritab)
                if( cline%defined('state') )then
                    call b%se%apply_sym_with_shift(os, symaxis_ori, shvec, p%state )
                else
                    call b%se%apply_sym_with_shift(os, symaxis_ori, shvec )
                endif
                call os%write(p%outfile)
            endif
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
        ! indicate completion (when run in a qsys env)
        if(p%l_distr_exec)then
            fname_finished = 'JOB_FINISHED_'//int2str_pad(p%part,p%numlen)
        else
            fname_finished = 'SYMSRCH_FINISHED'
        endif
        fnr = get_fileunit()
        open(unit=fnr, FILE=trim(fname_finished), STATUS='REPLACE', action='WRITE', iostat=file_stat)
        call fopen_err('In: commander_comlin :: symsrch', file_stat )
        close( unit=fnr )
    end subroutine exec_symsrch

end module simple_commander_comlin
