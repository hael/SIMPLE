! concrete commander: common-lines based clustering and search
module simple_commander_comlin
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_projector_hlev, only: project
use simple_sp_project,     only: sp_project
use simple_comlin_srch     ! use all in there
use simple_binoris_io,     only: binwrite_oritab, binread_oritab, binread_nlines
!use simple_binoris_io      ! use all in there
implicit none

public :: symsrch_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: symsrch_commander
  contains
    procedure :: execute      => exec_symsrch
end type symsrch_commander

contains

    !> for identification of the principal symmetry axis
    subroutine exec_symsrch( self, cline )
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                 :: p
        type(build)                  :: b
        type(ori)                    :: symaxis, orientation
        class(oris), pointer         :: symaxes => null()
        type(sp_project)             :: spproj
        type(oris)                   :: oshift, orientation_best, tmp_os
        real,            allocatable :: corrs(:)
        integer,         allocatable :: order(:)
        !integer                      :: fnr, file_stat, comlin_srch_nbest
        integer                      :: bestloc(1), nbest_here, noris, cnt, i, j, nl
        real                         :: shvec(3)
        character(len=STDLEN)        :: fname_finished
        character(len=32), parameter :: SYMSHTAB   = 'sym_3dshift'//trim(TXT_EXT)
        character(len=32), parameter :: SYMPROJSTK = 'sym_projs.mrc'
        character(len=32), parameter :: SYMPROJTAB = 'sym_projs'//trim(TXT_EXT)
        integer,           parameter :: NBEST = 30
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        p = params(cline)                                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.true., nooritab=.true.) ! general objects built (no oritab reading)
        call b%build_comlin_tbox(p) ! objects for common lines based alignment built
        ! SETUP
        if( (p%l_distr_exec .and. p%refine.eq.'no') .or. .not.p%l_distr_exec )then
            ! center volume
            call b%vol%read(p%vols(1))
            shvec = 0.
            if( p%center.eq.'yes' ) shvec = b%vol%center(p%cenlp,p%msk)
            if( p%l_distr_exec .and. p%part.eq.1 )then
                ! writes shifts for distributed execution
                call oshift%new(1)
                call oshift%set(1,'x',shvec(1))
                call oshift%set(1,'y',shvec(2))
                call oshift%set(1,'z',shvec(3))
                call oshift%write(trim(SYMSHTAB), [1,1])
                call oshift%kill
            endif
            ! generate projections
            call b%vol%mask(p%msk, 'soft')
            call b%vol%fft()
            call b%vol%expand_cmat(p%alpha)
            b%ref_imgs(1,:) = project(b%vol, b%e, p)
            if( p%l_distr_exec .and. p%part > 1 )then
                ! do nothing
            else
                ! writes projections images and orientations for subsequent reconstruction
                ! only in local and distributed (part=1) modes
                noris = b%e%get_noris()
                do i=1,noris
                    call b%ref_imgs(1,i)%write(SYMPROJSTK, i)
                enddo
                call b%e%write(SYMPROJTAB, [1,noris])
            endif
            ! expand over symmetry group
            cnt = 0
            do i=1,p%nptcls
                do j=1,p%nsym
                    cnt = cnt+1
                    b%imgs_sym(cnt) = b%ref_imgs(1,i)
                    call b%imgs_sym(cnt)%fft()
                end do
            end do
        endif
        ! COARSE SEARCH
        if( (p%l_distr_exec .and. p%refine.eq.'no') .or. .not.p%l_distr_exec )then
            call comlin_srch_init( b, p, 'simplex', 'sym')
            call comlin_coarsesrch_symaxis( [p%fromp,p%top], symaxes)
            if( p%l_distr_exec )then
                call symaxes%write(trim(p%fbody)//int2str_pad(p%part,p%numlen)//trim(TXT_EXT), [p%fromp,p%top])
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
                call symaxes%write('sympeaks'//trim(TXT_EXT), [1,nbest_here])
                deallocate(order)
                call tmp_os%kill
            endif
        endif
        ! FINE SEARCH
        if( p%refine.eq.'yes' .or. .not.p%l_distr_exec)then
            call orientation_best%new(1)
            if( p%l_distr_exec )then
                ! fetch orientation to refine
                nl = binread_nlines(p, p%oritab)
                call spproj%new_seg_with_ptr(nl, p%oritype, symaxes)
                call binread_oritab(p%oritab, spproj, symaxes, [1,nl])
                orientation = symaxes%get_ori(p%part)
                ! fetch refernce orientations
                nl = nlines(SYMPROJTAB)
                call b%e%new(nl)
                call b%e%read(SYMPROJTAB, [1,nl])
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
                        call b%imgs_sym(cnt)%fft()
                    end do
                end do
                ! search
                call comlin_srch_init( b, p, 'simplex', 'sym')
                call comlin_singlesrch_symaxis(orientation)
                call orientation_best%set_ori(1, orientation)
                call orientation_best%write(trim(p%fbody)//int2str_pad(p%part, p%numlen)//trim(TXT_EXT), [1,1])
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
                    call binwrite_oritab(p%outfile, b%spproj, b%a, [1,b%a%get_noris()])
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
