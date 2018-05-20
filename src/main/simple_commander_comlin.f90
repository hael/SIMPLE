! concrete commander: common-lines based clustering and search
module simple_commander_comlin
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_commander_base, only: commander_base
use simple_binoris_io,     only: binwrite_oritab, binread_oritab, binread_nlines
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
        use simple_comlin_srch     ! use all in there
        use simple_projector_hlev, only: reproject
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters)               :: params
        type(builder)                  :: build
        type(ori)                      :: orientation
        class(oris), pointer           :: psymaxes
        type(oris)                     :: oshift, orientation_best,symaxes
        integer                        :: noris, cnt, i, j
        real                           :: shvec(3)
        character(len=:),  allocatable :: symaxes_str
        character(len=STDLEN)          :: fname_finished
        character(len=32), parameter   :: SYMSHTAB   = 'sym_3dshift'//trim(TXT_EXT)
        character(len=32), parameter   :: SYMPROJSTK = 'sym_projs.mrc'
        character(len=32), parameter   :: SYMPROJTAB = 'sym_projs'//trim(TXT_EXT)
        nullify(psymaxes)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( params%refine.eq.'no' )then
            call build%build_comlin_tbox(params) ! objects for common lines based alignment built
            ! center volume
            call build%vol%read(params%vols(1))
            shvec = 0.
            if( params%center.eq.'yes' ) shvec = build%vol%center(params%cenlp,params%msk)
            if( params%part.eq.1 )then
                ! writes shifts for distributed execution
                call oshift%new(1)
                call oshift%set(1,'x',shvec(1))
                call oshift%set(1,'y',shvec(2))
                call oshift%set(1,'z',shvec(3))
                call oshift%write(trim(SYMSHTAB), [1,1])
                call oshift%kill
            endif
            ! generate projections
            call build%vol%mask(params%msk, 'soft')
            call build%vol%fft()
            call build%vol%expand_cmat(params%alpha)
            build%ref_imgs(1,:) = reproject(build%vol, build%e)
            if( params%part == 1 )then
                ! writes projections images and orientations for subsequent reconstruction
                ! only in local and distributed (part=1) modes
                noris = build%e%get_noris()
                do i=1,noris
                    call build%ref_imgs(1,i)%write(SYMPROJSTK, i)
                enddo
                call build%e%write(SYMPROJTAB, [1,noris])
            endif
            ! expand over symmetry group
            cnt = 0
            do i=1,params%nptcls
                do j=1,params%nsym
                    cnt = cnt+1
                    build%imgs_sym(cnt) = build%ref_imgs(1,i)
                    call build%imgs_sym(cnt)%fft()
                end do
            end do
        endif
        ! COARSE SEARCH
        if( params%refine.eq.'no' )then
            call comlin_srch_init(  'simplex', 'sym')
            call comlin_coarsesrch_symaxis( [params%fromp,params%top], symaxes)
            build%spproj%os_cls3D = symaxes
            symaxes_str       = trim(params%fbody)//int2str_pad(params%part,params%numlen)//trim(METADATA_EXT)
            call binwrite_oritab(symaxes_str, build%spproj, symaxes, [params%fromp,params%top], isegment=CLS3D_SEG)
        endif
        ! FINE SEARCH
        if( params%refine.eq.'yes' )then
            ! fetch orientation to refine
            orientation = build%spproj%os_cls3D%get_ori(params%part)
            ! build%a contained the coarse solutions, needs to be rebuilt for common line toolbox
            call build%a%new(params%nptcls)
            call build%a%spiral
            call build%build_comlin_tbox(params) ! objects for common lines based alignment built
            ! fetch reference orientations
            call build%e%new(params%nspace)
            call build%e%read(SYMPROJTAB, [1,params%nspace])
            do i=1,params%nptcls
                call build%ref_imgs(1,i)%new([params%box, params%box, 1], params%smpd)
                call build%ref_imgs(1,i)%read(SYMPROJSTK, i)
            enddo
            ! expand over symmetry group
            cnt = 0
            do i=1,params%nptcls
                do j=1,params%nsym
                    cnt = cnt+1
                    build%imgs_sym(cnt) = build%ref_imgs(1,i)
                    call build%imgs_sym(cnt)%fft()
                end do
            end do
            ! search
            call comlin_srch_init( 'simplex', 'sym')
            call comlin_singlesrch_symaxis(orientation)
            call orientation_best%new(1)
            call orientation_best%set_ori(1, orientation)
            symaxes_str       = trim(params%fbody)//int2str_pad(params%part,params%numlen)//trim(METADATA_EXT)
            build%spproj%os_cls3D = orientation_best
            call binwrite_oritab(symaxes_str, build%spproj, orientation_best, [1,1], isegment=CLS3D_SEG)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
        ! indicate completion (when run in a qsys env)
        fname_finished = 'JOB_FINISHED_'//int2str_pad(params%part,params%numlen)
        call simple_touch(trim(fname_finished), errmsg='In: commander_comlin :: exec_symsrch finished' )
    end subroutine exec_symsrch

end module simple_commander_comlin
