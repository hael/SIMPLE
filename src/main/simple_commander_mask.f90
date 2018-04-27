! concrete commander: masking routines
module simple_commander_mask
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
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
        type(build)                :: b
        type(params)               :: p
        type(automask2D_commander) :: automask2D
        type(image)                :: mskvol
        type(atoms)                :: pdb
        type(masker)               :: msker
        character(len=STDLEN)      :: pdbout_fname
        integer                    :: ldim(3)
        p = params(cline)  ! parameters generated
        p%boxmatch = p%box ! turns off boxmatch logics
        if( cline%defined('stk') .and. cline%defined('vol1')   ) stop 'Cannot operate on images AND volume at once'
        if( p%automsk.eq.'yes'   .and..not.cline%defined('mw') ) stop 'Missing mw argument for automasking'
        if( p%msktype.ne.'soft'  .and. p%msktype.ne.'hard' .and. p%msktype.ne.'cavg' ) stop 'Invalid mask type'
        if( cline%defined('stk') )then
            ! 2D
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            if( p%automsk.eq.'yes' )then
                ! auto
                if( .not. cline%defined('amsklp') )call cline%set('amsklp', 25.)
                if( .not. cline%defined('edge')   )call cline%set('edge', 10.)
                call exec_automask2D( automask2D, cline )
            else if( cline%defined('msk') .or. cline%defined('inner') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call mask_imgfile(p%stk, p%outstk, p%msk, p%smpd, inner=p%inner, width=p%width, which=p%msktype)
                    else
                        call mask_imgfile(p%stk, p%outstk, p%msk, p%smpd, inner=p%inner, which=p%msktype)
                    endif
                else
                    call mask_imgfile(p%stk, p%outstk, p%msk, p%smpd, which=p%msktype)
                endif
            else if( p%taper_edges.eq.'yes' )then
                call taper_edges_imgfile(p%stk, p%outstk, p%smpd)
            else
                stop 'Nothing to do!'
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call b%build_general_tbox(p, cline) ! general objects built
            ! reallocate vol (boxmatch issue)
            call b%vol%new([p%box, p%box, p%box], p%smpd)
            if( .not. file_exists(p%vols(1)) ) stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( cline%defined('mskfile') )then
                ! from file
                if( .not. file_exists(p%mskfile) ) stop 'Cannot find input mskfile'
                ldim = b%vol%get_ldim()
                call mskvol%new(ldim, p%smpd)
                call mskvol%read(p%mskfile)
                call b%vol%mul(mskvol)
                call mskvol%kill
                if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
            else if( p%automsk.eq.'yes' )then
                stop '3D automasking now deferred to program: postprocess'
            else if( cline%defined('msk') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner, width=p%width)
                    else
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner)
                    endif
                else
                    call b%vol%mask(p%msk, p%msktype)
                endif
                if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
            else if( cline%defined('pdbfile') )then
                ! focus masking
                call pdb%new(p%pdbfile)
                pdbout_fname = trim(get_fbody(p%pdbfile, 'pdb')) // '_centered'
                if( p%center.eq.'yes' )then
                    call msker%mask_from_pdb(p, pdb, b%vol, os=b%a, pdbout=pdbout_fname)
                else
                    call msker%mask_from_pdb(p, pdb, b%vol)
                endif
                call b%a%write(p%outfile)
                call b%vol%write(p%outvol)
                call msker%write('maskfile'//p%ext)
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
        type(params) :: p
        type(build)  :: b
        integer      :: iptcl
        p = params(cline)                                 ! parameters generated
        p%boxmatch = p%box                                ! turns off boxmatch logics
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        write(*,'(A,F8.2,A)') '>>> AUTOMASK LOW-PASS:',        p%amsklp, ' ANGSTROMS'
        write(*,'(A,I3,A)')   '>>> AUTOMASK SOFT EDGE WIDTH:', p%edge,   ' PIXELS'
        p%outstk = add2fbody(p%stk, p%ext, 'msk')
        do iptcl=1,p%nptcls
            call b%img%read(p%stk, iptcl)
            call b%mskimg%apply_2Denvmask22Dref(b%img)
            call b%img%write(p%outstk, iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK2D NORMAL STOP ****')
    end subroutine exec_automask2D

    !> for generating an envelope mask for resolution estimation
    subroutine exec_resmask( self, cline )
        use simple_masker, only: masker
        class(resmask_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)  :: p
        !type(build)   :: b
        type(masker)  :: mskvol
        p = params(cline)
        p%boxmatch = p%box                  ! turns off boxmatch logics
        !call b%build_general_tbox(p, cline) ! general objects built
        call mskvol%new([p%box,p%box,p%box], p%smpd)
        if( file_exists(p%mskfile) )then
            call mskvol%read(p%mskfile)
            call mskvol%resmask(p)
            call mskvol%write('resmask'//p%ext)
            call mskvol%kill
        else
            write(*,*) 'the inputted mskfile: ', trim(p%mskfile)
            stop 'does not exists in cwd; commander_mask :: exec_resmask'
        endif
         ! end gracefully
        call simple_end('**** SIMPLE_RESMASK NORMAL STOP ****')
    end subroutine exec_resmask

end module simple_commander_mask
