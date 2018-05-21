! concrete commander: masking routines
module simple_commander_mask
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
implicit none

public :: mask_commander
public :: automask2D_commander
public :: resmask_commander
private

type, extends(commander_base) :: mask_commander
 contains
   procedure :: execute      => exec_mask
end type mask_commander
type, extends(commander_base) :: automask2D_commander
  contains
    procedure :: execute      => exec_automask2D
end type automask2D_commander
type, extends(commander_base) :: resmask_commander
  contains
    procedure :: execute      => exec_resmask
end type resmask_commander

contains

    !> for masking images and volumes
    subroutine exec_mask( self, cline )
        use simple_image,       only: image
        use simple_procimgfile, only: mask_imgfile, taper_edges_imgfile
        use simple_atoms,       only: atoms
        use simple_masker,      only: masker
        class(mask_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(parameters)           :: params
        type(builder)              :: build
        type(automask2D_commander) :: automask2D
        type(image)                :: mskvol
        type(atoms)                :: pdb
        type(masker)               :: msker
        character(len=STDLEN)      :: pdbout_fname
        integer                    :: ldim(3)
        if( cline%defined('stk') .and. cline%defined('vol1')   ) stop 'Cannot operate on images AND volume at once'
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline,params,do3d=.false.,boxmatch_off=.true.)
            if( params%automsk.eq.'yes' )then
                ! auto
                if( .not. cline%defined('amsklp') )call cline%set('amsklp', 25.)
                if( .not. cline%defined('edge')   )call cline%set('edge', 10.)
                call exec_automask2D( automask2D, cline )
            else if( cline%defined('msk') .or. cline%defined('inner') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call mask_imgfile(params%stk, params%outstk, params%msk, params%smpd, inner=params%inner, width=params%width, which=params%msktype)
                    else
                        call mask_imgfile(params%stk, params%outstk, params%msk, params%smpd, inner=params%inner, which=params%msktype)
                    endif
                else
                    call mask_imgfile(params%stk, params%outstk, params%msk, params%smpd, which=params%msktype)
                endif
            else if( params%taper_edges.eq.'yes' )then
                call taper_edges_imgfile(params%stk, params%outstk, params%smpd)
            else
                stop 'Nothing to do!'
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call build%init_params_and_build_general_tbox(cline,params,do3d=.true.,boxmatch_off=.true.)
            if( .not. file_exists(params%vols(1)) ) stop 'Cannot find input volume'
            call build%vol%read(params%vols(1))
            if( cline%defined('mskfile') )then
                ! from file
                if( .not. file_exists(params%mskfile) ) stop 'Cannot find input mskfile'
                ldim = build%vol%get_ldim()
                call mskvol%new(ldim, params%smpd)
                call mskvol%read(params%mskfile)
                call build%vol%mul(mskvol)
                call mskvol%kill
                if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
            else if( cline%defined('msk') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call build%vol%mask(params%msk, params%msktype, inner=params%inner, width=params%width)
                    else
                        call build%vol%mask(params%msk, params%msktype, inner=params%inner)
                    endif
                else
                    call build%vol%mask(params%msk, params%msktype)
                endif
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
                call build%spproj_field%write(params%outfile)
                call build%vol%write(params%outvol)
                call msker%write('maskfile'//params%ext)
            else
                stop 'Nothing to do!'
            endif
        else
            stop 'No input images(s) or volume provided'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASK NORMAL STOP ****')
    end subroutine exec_mask

    !> for solvent flattening of class averages
    subroutine exec_automask2D( self, cline )
        class(automask2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: iptcl
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.,boxmatch_off=.true.)
        write(*,'(A,F8.2,A)') '>>> AUTOMASK LOW-PASS:',        params%amsklp, ' ANGSTROMS'
        write(*,'(A,I3,A)')   '>>> AUTOMASK SOFT EDGE WIDTH:', params%edge,   ' PIXELS'
        params%outstk = add2fbody(params%stk, params%ext, 'msk')
        do iptcl=1,params%nptcls
            call build%img%read(params%stk, iptcl)
            call build%mskimg%apply_2Denvmask22Dref(build%img)
            call build%img%write(params%outstk, iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK2D NORMAL STOP ****')
    end subroutine exec_automask2D

    !> for generating an envelope mask for resolution estimation
    subroutine exec_resmask( self, cline )
        use simple_masker, only: masker
        class(resmask_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(masker)     :: mskvol
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.,boxmatch_off=.true.)
        call mskvol%new([params%box,params%box,params%box], params%smpd)
        if( file_exists(params%mskfile) )then
            call mskvol%read(params%mskfile)
            call mskvol%resmask()
            call mskvol%write('resmask'//params%ext)
            call mskvol%kill
        else
            write(*,*) 'the inputted mskfile: ', trim(params%mskfile)
            stop 'does not exists in cwd; commander_mask :: exec_resmask'
        endif
         ! end gracefully
        call simple_end('**** SIMPLE_RESMASK NORMAL STOP ****')
    end subroutine exec_resmask

end module simple_commander_mask
