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
        type(params)                   :: p
        type(build)                    :: b
        type(ori)                      :: symaxis, orientation
        class(oris), pointer           :: psymaxes
        type(sp_project)               :: spproj
        type(oris)                     :: oshift, orientation_best, tmp_os, symaxes
        real,              allocatable :: corrs(:)
        integer,           allocatable :: order(:)
        integer                        :: bestloc(1), nbest_here, noris, cnt, i, j
        real                           :: shvec(3)
        character(len=:),  allocatable :: symaxes_str
        character(len=STDLEN)          :: fname_finished
        character(len=32), parameter   :: SYMSHTAB   = 'sym_3dshift'//trim(TXT_EXT)
        character(len=32), parameter   :: SYMPROJSTK = 'sym_projs.mrc'
        character(len=32), parameter   :: SYMPROJTAB = 'sym_projs'//trim(TXT_EXT)
        integer,           parameter   :: NBEST = 30
        nullify(psymaxes)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        p = params(cline)                                ! parameters generated
        ! SETUP
        call b%build_general_tbox(p, cline, do3d=.true.) ! general objects built (no oritab reading)
        if( p%refine.eq.'no' )then
            call b%build_comlin_tbox(p) ! objects for common lines based alignment built
            ! center volume
            call b%vol%read(p%vols(1))
            shvec = 0.
            if( p%center.eq.'yes' ) shvec = b%vol%center(p%cenlp,p%msk)
            if( p%part.eq.1 )then
                ! writes shifts for distributed execution
                call oshift%new_clean(1)
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
            if( p%part == 1 )then
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
        if( p%refine.eq.'no' )then
            call comlin_srch_init( b, p, 'simplex', 'sym')
            call comlin_coarsesrch_symaxis( [p%fromp,p%top], symaxes)
            b%spproj%os_cls3D = symaxes
            symaxes_str       = trim(p%fbody)//int2str_pad(p%part,p%numlen)//trim(METADATA_EXT)
            call binwrite_oritab(symaxes_str, b%spproj, symaxes, [p%fromp,p%top], isegment=CLS3D_SEG)
        endif
        ! FINE SEARCH
        if( p%refine.eq.'yes' )then
            ! fetch orientation to refine
            orientation = b%spproj%os_cls3D%get_ori(p%part)
            ! b%a contained the coarse solutions, needs to be rebuilt for common line toolbox
            call b%a%new_clean(p%nptcls)
            call b%a%spiral
            call b%build_comlin_tbox(p) ! objects for common lines based alignment built
            ! fetch reference orientations
            call b%e%new_clean(p%nspace)
            call b%e%read(SYMPROJTAB, [1,p%nspace])
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
            call orientation_best%new_clean(1)
            call orientation_best%set_ori(1, orientation)
            symaxes_str       = trim(p%fbody)//int2str_pad(p%part,p%numlen)//trim(METADATA_EXT)
            b%spproj%os_cls3D = orientation_best
            call binwrite_oritab(symaxes_str, b%spproj, orientation_best, [1,1], isegment=CLS3D_SEG)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
        ! indicate completion (when run in a qsys env)
        fname_finished = 'JOB_FINISHED_'//int2str_pad(p%part,p%numlen)
        call simple_touch(trim(fname_finished), errmsg='In: commander_comlin :: exec_symsrch finished' )
    end subroutine exec_symsrch

end module simple_commander_comlin
