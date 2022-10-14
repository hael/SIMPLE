! concrete commander: masking routines
module simple_commander_mask
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
use simple_default_clines
implicit none

public :: mask_commander
public :: automask2D_commander
public :: automask_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: mask_commander
 contains
   procedure :: execute      => exec_mask
end type mask_commander

type, extends(commander_base) :: automask2D_commander
  contains
    procedure :: execute      => exec_automask2D
end type automask2D_commander

type, extends(commander_base) :: automask_commander
 contains
   procedure :: execute      => exec_automask
end type automask_commander

contains

    !> for masking images and volumes
    subroutine exec_mask( self, cline )
        use simple_procimgstk, only: mask_imgfile, taper_edges_imgfile
        use simple_atoms,      only: atoms
        class(mask_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(parameters)           :: params
        type(builder)              :: build
        type(image)                :: mskvol
        type(atoms)                :: pdb
        type(masker)               :: msker
        character(len=STDLEN)      :: pdbout_fname
        integer                    :: ldim(3)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        if( cline%defined('stk') .and. cline%defined('vol1') ) THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
            if( cline%defined('mskdiam') )then
                ! spherical
                call mask_imgfile(params%stk, params%outstk, params%msk, params%smpd, which=params%msktype)
            else if( params%taper_edges.eq.'yes' )then
                call taper_edges_imgfile(params%stk, params%outstk, params%smpd)
            else
                THROW_HARD('Nothing to do!')
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
            if( .not. file_exists(params%vols(1)) ) THROW_HARD('Cannot find input volume')
            call build%vol%read(params%vols(1))
            if( params%l_filemsk )then
                ! from file
                ldim = build%vol%get_ldim()
                call mskvol%new(ldim, params%smpd)
                call mskvol%read(params%mskfile)
                if( cline%defined('lp_backgr') )then
                    call build%vol%lp_background(mskvol,params%lp_backgr)
                else
                    call build%vol%zero_background
                    call build%vol%mul(mskvol)
                endif
                call mskvol%kill
                if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
            else if( cline%defined('mskdiam') )then
                ! spherical
                call build%vol%mask(params%msk, params%msktype)
                if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
            else if( cline%defined('pdbfile') )then
                ! focus masking
                call pdb%new(params%pdbfile)
                pdbout_fname = trim(get_fbody(params%pdbfile, 'pdb')) // '_centered'
                if( params%center.eq.'yes' )then
                    call msker%mask_from_pdb( pdb, build%vol, os=build%spproj_field, pdbout=pdbout_fname)
                else
                    call msker%mask_from_pdb( pdb, build%vol)
                endif
                call build%spproj%write_segment_inside(params%oritype,params%projfile)
                call build%spproj_field%write(params%outfile)
                call build%vol%write(params%outvol)
                call msker%write('maskfile'//params%ext)
            else
                THROW_HARD('Nothing to do!')
            endif
        else
            THROW_HARD('No input images(s) or volume provided')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASK NORMAL STOP ****')
    end subroutine exec_mask

    !> for automasking of class averages
    subroutine exec_automask2D( self, cline )
        use simple_masker, only: automask2D
        class(automask2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)         :: params
        type(image), allocatable :: imgs(:), masks(:)
        real,        allocatable :: diams(:)
        integer :: ldim(3), n, i
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call set_automask2D_defaults(cline)
        call params%new(cline)
        call find_ldim_nptcls(params%stk, ldim, n)
        ldim(3) = 1 ! 2D
        allocate(imgs(n), masks(n), diams(n))
        diams = 0.
        do i = 1,n
            call imgs(i)%new(ldim, params%smpd)
            call imgs(i)%read(params%stk, i)
            call masks(i)%copy(imgs(i))
        end do
        call automask2D(masks, params%ngrow, nint(params%winsz), params%edge, diams)
        do i = 1,n
            call imgs(i)%mul(masks(i))
            call imgs(i)%write('automasked.mrc', i)
            call masks(i)%kill
            call imgs(i)%kill
        end do
        deallocate(imgs, masks, diams)
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK2D NORMAL STOP ****')
    end subroutine exec_automask2D

    !> for automasking of volume
    subroutine exec_automask( self, cline )
        class(automask_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(builder)    :: build
        type(parameters) :: params
        type(masker)     :: mskvol
        character(len=:), allocatable :: fname_out
        call params%new(cline)
        call build%build_spproj(params, cline)
        call build%build_general_tbox(params, cline)
        call build%vol%read(params%vols(1))
        if( cline%defined('mw') .and. cline%defined('thres') )then
            call mskvol%automask3D(build%vol)
        else
            call mskvol%automask3D_otsu(build%vol)
        endif
        call mskvol%write('automask'//params%ext)
        fname_out = basename(add2fbody(trim(params%vols(1)), params%ext, '_automsk'))
        call build%vol%write(fname_out)
        write(logfhandle,'(A)') '>>> WROTE OUTPUT '//'automask'//params%ext//' & '//fname_out
        call simple_end('**** SIMPLE_AUTOMASK NORMAL STOP ****')
    end subroutine exec_automask

end module simple_commander_mask
