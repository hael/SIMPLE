! concrete commander: common-lines based clustering and search

module simple_commander_comlin
#include "simple_lib.f08"

use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: comlin_smat_commander
public :: symsrch_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: comlin_smat_commander
  contains
    procedure :: execute      => exec_comlin_smat
end type comlin_smat_commander
type, extends(commander_base) :: symsrch_commander
  contains
    procedure :: execute      => exec_symsrch
end type symsrch_commander

contains

    ! calculates a similarity matrix based on common line correlations to measure 3D similarity of 2D images
    subroutine exec_comlin_smat( self, cline )
        use simple_comlin_srch   ! use all in there
        use simple_ori,          only: ori
        use simple_imgfile,      only: imgfile
        use simple_comlin,       only: comlin
        use simple_qsys_funs,      only: qsys_job_finished
        class(comlin_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target          :: p
        type(build),  target          :: b
        integer                       :: iptcl, jptcl, funit, io_stat
        integer                       :: ntot, npairs, ipair
        real,    allocatable          :: corrmat(:,:), corrs(:)
        integer, allocatable          :: pairs(:,:)
        character(len=:), allocatable :: fname
        p = params(cline, .false.)                           ! constants & derived constants produced
        call b%build_general_tbox(p, cline, .false., .true.) ! general objects built (no oritab reading)
        allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
        allocchk('In: simple_comlin_smat, 1')
        DebugPrint  'analysing this number of objects: ', p%nptcls
        do iptcl=1,p%nptcls
            call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd)
            call b%imgs_sym(iptcl)%read(p%stk, iptcl)
            ! apply a soft-edged mask
            call b%imgs_sym(iptcl)%mask(p%msk, 'soft')
            ! Fourier transform
            call b%imgs_sym(iptcl)%fwd_ft
        end do
        b%clins = comlin(b%a, b%imgs_sym, p%lp)
        if( cline%defined('part') )then
            npairs = p%top - p%fromp + 1
            DebugPrint  'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            allocchk('In: simple_comlin_smat, 2')
            ! read the pairs
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin', stat=alloc_stat)
            allocchk("In:  simple_comlin_smat, 3")
            if( .not. file_exists(fname) )then
                write(*,*) 'file: ', fname, ' does not exist!'
                write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
                stop 'I/O error; simple_comlin_smat'
            endif
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM', iostat=io_stat)
            call fileio_errmsg('simple_comlin_smat opening  '//trim(fname), io_stat)
            DebugPrint  'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! check if the read was successful
            if( io_stat .ne. 0 ) call fileio_errmsg('simple_comlin_smat reading  '//trim(fname), io_stat)
            call fclose(funit,errmsg='simple_comlin_smat closing  '//trim(fname))
            deallocate(fname)
            ! calculate the similarities
            call comlin_srch_init( b, p, 'simplex', 'pair' )
            do ipair=p%fromp,p%top
                p%iptcl = pairs(ipair,1)
                p%jptcl = pairs(ipair,2)
                corrs(ipair) = comlin_srch_pair()
            end do
            ! write the similarities
            allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin')
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, access='STREAM', iostat=io_stat)
            call fileio_errmsg('simple_comlin_smat opening  '//trim(fname), io_stat)
            write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
            ! Check if the write was successful
            if( io_stat .ne. 0 )call fileio_errmsg('simple_comlin_smat writing  '//trim(fname), io_stat)

            call fclose(funit, errmsg='simple_comlin_smat opening  '//trim(fname))
            deallocate(fname, corrs, pairs)
            call qsys_job_finished(p,'simple_commander_comlin :: exec_comlin_smat')
        else
            allocate(corrmat(p%nptcls,p%nptcls), stat=alloc_stat)
            allocchk('In: simple_comlin_smat, 3')
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

            call fopen(funit, status='REPLACE', action='WRITE', file='clin_smat.bin', access='STREAM', iostat=io_stat)
            call fileio_errmsg('simple_comlin_smat opening  clin_smat.bin', io_stat)
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if( io_stat .ne. 0 ) call fileio_errmsg('simple_comlin_smat writing  clin_smat.bin', io_stat)

            call fclose(funit, errmsg='simple_comlin_smat closing clin_smat.bin ')
            deallocate(corrmat)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_COMLIN_SMAT NORMAL STOP ****')
    end subroutine exec_comlin_smat

    !> for identification of the principal symmetry axis
    subroutine exec_symsrch( self, cline )
        use simple_strings,        only: int2str_pad
        use simple_oris,           only: oris
        use simple_ori,            only: ori
        use simple_projector_hlev, only: projvol
        use simple_comlin_srch     ! use all in there
        use simple_binoris_io,     only: binwrite_oritab, binread_oritab, binread_nlines
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                 :: p
        type(build)                  :: b
        type(ori)                    :: symaxis, orientation
        type(oris)                   :: os, oshift, symaxes, orientation_best, tmp_os
        real,            allocatable :: corrs(:)
        integer,         allocatable :: order(:)
        integer                      :: fnr, file_stat, comlin_srch_nbest, cnt, i, j, nl
        integer                      :: bestloc(1), nbest_here, noris
        real                         :: shvec(3)
        character(len=STDLEN)        :: fname_finished
        character(len=32), parameter :: SYMSHTAB   = 'sym_3dshift'//METADATEXT
        character(len=32), parameter :: SYMPROJSTK = 'sym_projs.mrc'
        character(len=32), parameter :: SYMPROJTAB = 'sym_projs'//METADATEXT
        integer,           parameter :: NBEST = 30
        p = params(cline)                                   ! parameters generated
        call b%build_general_tbox(p, cline, .true., nooritab=.true.) ! general objects built (no oritab reading)
        call b%build_comlin_tbox(p)  ! objects for common lines based alignment built
        ! SETUP
        if( (p%l_distr_exec .and. p%refine.eq.'no') .or. .not.p%l_distr_exec )then
            ! center volume
            call b%vol%read(p%vols(1))
            shvec = b%vol%center(p%cenlp,'no',p%msk)
            if( p%l_distr_exec .and. p%part.eq.1 )then
                ! writes shifts for distributed execution
                call oshift%new(1)
                call oshift%set(1,'x',shvec(1))
                call oshift%set(1,'y',shvec(2))
                call oshift%set(1,'z',shvec(3))
                call binwrite_oritab(trim(SYMSHTAB), oshift, [1,1])
                call oshift%kill
            endif
            ! generate projections
            call b%vol%mask(p%msk, 'soft')
            call b%vol%fwd_ft
            call b%vol%expand_cmat
            b%ref_imgs(1,:) = projvol(b%vol, b%e, p)
            if( p%l_distr_exec .and. p%part > 1 )then
                ! do nothing
            else
                ! writes projections images and orientations for subsequent reconstruction
                ! only in local and distributed (part=1) modes
                noris = b%e%get_noris()
                do i=1,noris 
                    call b%ref_imgs(1,i)%write(SYMPROJSTK, i)
                enddo
                call binwrite_oritab(SYMPROJTAB, b%e, [1,noris])
            endif
            ! expand over symmetry group
            cnt = 0
            do i=1,p%nptcls
                do j=1,p%nsym
                    cnt = cnt+1
                    b%imgs_sym(cnt) = b%ref_imgs(1,i)
                    call b%imgs_sym(cnt)%fwd_ft
                end do
            end do
        endif
        ! COARSE SEARCH
        if( (p%l_distr_exec .and. p%refine.eq.'no') .or. .not.p%l_distr_exec )then
            call comlin_srch_init( b, p, 'simplex', 'sym')
            call comlin_coarsesrch_symaxis( [p%fromp,p%top], symaxes)
            if( p%l_distr_exec )then
                call binwrite_oritab(trim(p%fbody)//int2str_pad(p%part,p%numlen)//METADATEXT, symaxes, [p%fromp,p%top])
            else
                noris      = symaxes%get_noris()
                nbest_here = min(NBEST, noris)
                order      = symaxes%order_corr()
                call tmp_os%new(nbest_here)
                cnt = 0
                do i = noris, noris-nbest_here+1, -1
                    cnt = cnt + 1
                    call tmp_os%set_ori(cnt, symaxes%get_ori(order(i)))
                enddo
                symaxes = tmp_os
                call binwrite_oritab('sympeaks'//METADATEXT, symaxes, [1,nbest_here])
                deallocate(order)
                call tmp_os%kill
            endif
        endif
        ! FINE SEARCH
        if( p%refine.eq.'yes' .or. .not.p%l_distr_exec)then
            call orientation_best%new(1)
            if( p%l_distr_exec )then
                ! fetch orientation to refine
                nl = binread_nlines(p%oritab)
                call symaxes%new(nl)
                call binread_oritab(p%oritab, symaxes, [1,nl])
                orientation = symaxes%get_ori(p%part)
                ! fetch refernce orientations
                nl = binread_nlines(SYMPROJTAB)
                call b%e%new(nl)
                call binread_oritab(SYMPROJTAB, b%e, [1,nl])
                do i=1,p%nptcls
                    call b%ref_imgs(1,i)%new([p%box, p%box, 1], p%smpd)
                    call b%ref_imgs(1,i)%read(SYMPROJSTK, i)
                enddo
                ! expand over symmetry group
                cnt = 0
                do i=1,p%nptcls
                    do j=1,p%nsym
                        cnt = cnt+1
                        b%imgs_sym(cnt) = b%ref_imgs(1,i)
                        call b%imgs_sym(cnt)%fwd_ft
                    end do
                end do
                ! search
                call comlin_srch_init( b, p, 'simplex', 'sym')
                call comlin_singlesrch_symaxis(orientation)
                call orientation_best%set_ori(1, orientation)
                call binwrite_oritab(trim(p%fbody)//int2str_pad(p%part, p%numlen)//METADATEXT, orientation_best, [1,1])
            else
                ! search selected peaks in non-distributed modes
                write(*,'(A)') '>>> CONTINOUS SYMMETRY AXIS REFINEMENT'
                do  i = 1, nbest_here
                    call progress(i, nbest_here)
                    orientation = symaxes%get_ori(i)
                    call comlin_singlesrch_symaxis(orientation)
                    call symaxes%set_ori(i, orientation)
                enddo
                corrs   = symaxes%get_all('corr')
                bestloc = maxloc(corrs)
                symaxis = symaxes%get_ori(bestloc(1))
                write(*,'(A)') '>>> FOUND SYMMETRY AXIS ORIENTATION:'
                call symaxis%print_ori()
                if( cline%defined('oritab') )then
                    ! rotate the orientations & transfer the 3d shifts to 2d
                    shvec = -1.*shvec
                    if( cline%defined('state') )then
                        call b%se%apply_sym_with_shift(b%a, symaxis, shvec, p%state )
                    else
                        call b%se%apply_sym_with_shift(b%a, symaxis, shvec )
                    endif
                    call binwrite_oritab(p%outfile, b%a, [1,b%a%get_noris()])
                endif
            endif
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
        ! indicate completion (when run in a qsys env)
        if(p%l_distr_exec)then
            fname_finished = 'JOB_FINISHED_'//int2str_pad(p%part,p%numlen)
            call simple_touch(trim(fname_finished), errmsg='In: commander_comlin :: exec_symsrch finished' )
         endif
    end subroutine exec_symsrch

end module simple_commander_comlin
