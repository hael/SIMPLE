module simple_commander_quant
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_image,          only: image
use simple_binimage,       only: binimage
use simple_nanoparticle
implicit none

public :: detect_atoms_commander
public :: atoms_rmsd_commander
public :: radial_dependent_stats_commander
public :: atom_cluster_analysis_commander
public :: nano_softmask_commander
public :: geometry_analysis_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: detect_atoms_commander
  contains
    procedure :: execute      => exec_detect_atoms
end type detect_atoms_commander
type, extends(commander_base) :: atoms_rmsd_commander
  contains
    procedure :: execute      => exec_atoms_rmsd
end type atoms_rmsd_commander
type, extends(commander_base) :: radial_dependent_stats_commander
  contains
    procedure :: execute      => exec_radial_dependent_stats
end type radial_dependent_stats_commander
type, extends(commander_base) :: atom_cluster_analysis_commander
  contains
    procedure :: execute      => exec_atom_cluster_analysis
end type atom_cluster_analysis_commander
type, extends(commander_base) :: nano_softmask_commander
  contains
    procedure :: execute      => exec_nano_softmask
end type nano_softmask_commander
type, extends(commander_base) :: geometry_analysis_commander
  contains
    procedure :: execute      => exec_geometry_analysis
end type geometry_analysis_commander

contains

    ! Performs preprocessing on the nanoparticle map, identifies atomic positions,
    ! validates them and write them to disk
    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_detect_atoms')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_detect_atoms')
        endif
        call nano%new(params%vols(1), params%smpd, params%element)
        ! execute
        if(cline%defined('thres')) then
          call nano%identify_atomic_pos(nint(params%thres))
        else
          call nano%identify_atomic_pos()
        endif
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

    subroutine exec_atoms_rmsd( self, cline )
        use simple_commander_volops, only: dock_volpair_commander
        use simple_ori,              only: ori
        use simple_atoms,            only: atoms
        class(atoms_rmsd_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline !< command line input
        type(dock_volpair_commander) :: xdock_volpair
        character(len=STDLEN)  :: fname1, fname2
        type(parameters)       :: params
        type(nanoparticle) :: nano1, nano2
        type(cmdline)          :: cline_dock
        type(ori)   :: orientation
        type(oris)  :: ori2read
        type(atoms) :: atom_coord
        real        :: cxyz(3)
        real        :: smpd, mat(3,3), shift(3)
        integer     :: ldim(3), nptcls
        integer     :: i
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_atoms_rmsd')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_atoms_rmsd')
        endif
        if( .not. cline%defined('vol2') )then
            THROW_HARD('ERROR! vol2 needs to be present; exec_atoms_rmsd')
        endif
        if( .not. cline%defined('dock') )then
            THROW_HARD('ERROR! dock needs to be present; exec_atoms_rmsd')
        endif
        if(cline%defined('element')) then
            call nano1%new(params%vols(1), params%smpd,params%element)
            call nano2%new(params%vols(2), params%smpd,params%element)
        else
            call nano1%new(params%vols(1), params%smpd)
            call nano2%new(params%vols(2), params%smpd)
        endif
        fname1 = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        fname2 = get_fbody(trim(basename(params%vols(2))), trim(fname2ext(params%vols(2))))
        call nano1%set_atomic_coords(trim(fname1)//'_atom_centers.pdb')
        ! execute
        if(params%dock .eq. 'yes') then
            cline_dock = cline
            call cline_dock%set('lpstart', 1.)
            call cline_dock%set('lpstop',  3.)
            call cline_dock%set('msk',  60.)
            call cline_dock%set('mkdir', 'no')
            call cline_dock%set('nthr',0.)
            call cline_dock%set('outfile', 'algndoc.txt')
            call cline_dock%set('outvol', './'//trim(fname2)//'docked.mrc')
            call xdock_volpair%execute(cline_dock)
            !1) Center the coords
            call find_ldim_nptcls(params%vols(2),ldim, nptcls,smpd)
            cxyz = (real(ldim)/2.)*smpd
            call atom_coord%new(trim(fname2)//'_atom_centers.pdb')
            call atom_coord%translate(-cxyz)
            !2) Rotate the coords
            call ori2read%new(1)
            call ori2read%read('algndoc.txt')
            call ori2read%get_ori(1,orientation)
            call ori2read%kill
            mat = orientation%get_mat()
            call atom_coord%rotate(mat)
            !3) Translate back
            call atom_coord%translate(cxyz)
            !4) Shift
            shift = orientation%get_3Dshift()
            call orientation%kill
            call atom_coord%translate(shift)
            !5) Set new coords
            call atom_coord%writePDB(    trim(fname2)//'_ROTATEDatom_centers')
            call atom_coord%kill
            call nano2%set_atomic_coords(trim(fname2)//'_ROTATEDatom_centers.pdb')
        else ! no DOCKING
            call nano2%set_atomic_coords(trim(fname2)//'_atom_centers.pdb')
        endif
        ! RMSD calculation
        call nano1%atoms_rmsd(nano2)
        ! kill
        call nano1%kill
        call nano2%kill
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_RMSD NORMAL STOP ****')
    end subroutine exec_atoms_rmsd

    ! Calculates distances distribution across the whole nanoparticle
    ! and radial dependent statistics.
    subroutine exec_radial_dependent_stats( self, cline )
        class(radial_dependent_stats_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle) :: nano
        integer :: ldim(3), nptcls
        real    :: smpd
        real    :: min_rad, max_rad, step
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_radial_dependent_stats')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_radial_dependent_stats')
        endif
        min_rad = params%min_rad
        max_rad = params%max_rad
        step    = params%stepsz
        if(min_rad > max_rad) THROW_HARD('Minimum radius has to be smaller then maximum radius! exec_radial_sym_test')
        if(step > max_rad-min_rad) THROW_HARD('Inputted too big stepsz! exec_radial_sym_test')
        if(cline%defined('element')) then
            call nano%new(params%vols(1), params%smpd,params%element)
        else
            call nano%new(params%vols(1), params%smpd)
        endif
        ! execute
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords(trim(fname)//'_atom_centers.pdb')
        call nano%set_img(trim(fname)//'CC.mrc', 'img_cc')
        call nano%radial_dependent_stats(min_rad,max_rad,step)
        ! fetch again, after killing
        if(cline%defined('element')) then
            call nano%new(params%vols(1), params%smpd,params%element)
        else
            call nano%new(params%vols(1), params%smpd)
        endif
        call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
        call nano%set_img('../'//trim(fname)//'CC.mrc', 'img_cc')
        ! calculate intensity statistics
        call nano%atom_intensity_stats()
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_RADIAL_DEPENDENT_STATS NORMAL STOP ****')
    end subroutine exec_radial_dependent_stats

    subroutine exec_atom_cluster_analysis( self, cline )
        class(atom_cluster_analysis_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle)     :: nano
        real, allocatable      :: max_intensity(:)
        real    :: smpd
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_atom_cluster_analysis')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_atom_cluster_analysis')
        endif
        if( .not. cline%defined('clustermode') )then
            THROW_HARD('ERROR! clustermode needs to be present; exec_atom_cluster_analysis')
        endif
        if( .not. cline%defined('thres') )then
            THROW_HARD('ERROR! thres needs to be present; exec_atom_cluster_analysis')
        endif
        call nano%new(params%vols(1), params%smpd)
        ! execute
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords(trim(fname)//'_atom_centers.pdb')
        call nano%set_img(trim(fname)//'CC.mrc', 'img_cc')
        select case(trim(params%clustermode))
            case('ar')
              call nano%cluster_ar(params%thres)
            case('dist')
              call nano%cluster_interdist(params%thres)
            case('ang')
              call nano%cluster_ang(params%thres)
            case('maxint')
              call nano%atom_intensity_stats(max_intensity)
              call nano%cluster_atom_intensity(max_intensity)
            case DEFAULT
                write(logfhandle,*) 'clustermode: ', trim(params%clustermode)
                THROW_HARD('unsupported clustermode; exec_atom_cluster_analysis')
        end select
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_ATOM_CLUSTER_ANALYSIS NORMAL STOP ****')
    end subroutine exec_atom_cluster_analysis

    subroutine exec_nano_softmask( self, cline )
        class(nano_softmask_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle) :: nano
        real  :: smpd
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_nano_softmask')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_nano_softmask')
        endif
        if(cline%defined('element')) then
            call nano%new(params%vols(1), params%smpd,params%element)
        else
            call nano%new(params%vols(1), params%smpd)
        endif
        ! fetch img_bin
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_img(trim(fname)//'BIN.mrc','img_bin')
        ! execute
        call nano%make_soft_mask()
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NANO_SOFTMASK NORMAL STOP ****')
    end subroutine exec_nano_softmask

    subroutine exec_geometry_analysis( self, cline )
        class(geometry_analysis_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle) :: nano
        real    :: smpd
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_geometry_analysis')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_geometry_analysis')
        endif
        if( .not. cline%defined('pdbfile') )then
            THROW_HARD('ERROR! pdbfile needs to be present; exec_geometry_analysis')
        endif
        if(cline%defined('element')) then
            call nano%new(params%vols(1), params%smpd,params%element)
        else
            call nano%new(params%vols(1), params%smpd)
        endif
        ! fetch img_bin, img_cc and atomic positions
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_img('../'//trim(fname)//'BIN.mrc','img_bin')
        call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
        call nano%set_img('../'//trim(fname)//'CC.mrc', 'img_cc')
        ! execute
          call nano%geometry_analysis(trim(params%pdbfile))
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_GEOMETRY_ANALYSIS NORMAL STOP ****')
    end subroutine exec_geometry_analysis

end module simple_commander_quant
