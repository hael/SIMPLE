module simple_commander_quant
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_image,          only: image
use simple_binimage,       only: binimage
use simple_nanoparticle,   only: nanoparticle
use simple_dock_coords,    only: dock_coords_init, dock_coords_minimize
implicit none

public :: detect_atoms_commander
public :: atoms_stats_commander
public :: atom_cluster_analysis_commander
public :: nano_softmask_commander
public :: geometry_analysis_commander
public :: plot_atom_commander
public :: dock_coords_commander
public :: atoms_mask_commander
public :: strain_analysis_commander

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: detect_atoms_commander
  contains
    procedure :: execute      => exec_detect_atoms
end type detect_atoms_commander
type, extends(commander_base) :: dock_coords_commander
  contains
    procedure :: execute      => exec_dock_coords
end type dock_coords_commander
type, extends(commander_base) :: atoms_stats_commander
  contains
    procedure :: execute      => exec_atoms_stats
end type atoms_stats_commander
type, extends(commander_base) :: atom_cluster_analysis_commander
  contains
    procedure :: execute      => exec_atom_cluster_analysis
end type atom_cluster_analysis_commander
type, extends(commander_base) :: nano_softmask_commander
  contains
    procedure :: execute      => exec_nano_softmask
end type nano_softmask_commander
type, extends(commander_base) :: atoms_mask_commander
  contains
    procedure :: execute      => exec_atoms_mask
end type atoms_mask_commander
type, extends(commander_base) :: geometry_analysis_commander
  contains
    procedure :: execute      => exec_geometry_analysis
end type geometry_analysis_commander
type, extends(commander_base) :: plot_atom_commander
  contains
    procedure :: execute      => exec_plot_atom
end type plot_atom_commander
type, extends(commander_base) :: strain_analysis_commander
  contains
    procedure :: execute      => exec_strain_analysis
end type strain_analysis_commander


contains

    ! Performs preprocessing on the nanoparticle map, identifies atomic positions,
    ! validates them and write them to disk
    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano
        integer            :: cn_thresh
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_detect_atoms')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_detect_atoms')
        endif
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        call nano%new(params%vols(1), params%smpd, params%element)
        ! volume soft-edge masking
        call nano%mask(params%msk)
        ! execute
        if(cline%defined('thres')) then      ! threshold for binarization
          if(cline%defined('cn_thres')) then ! contact score threshold for outliers removal
            if(.not. cline%defined('cn_type')) then
              THROW_HARD('ERROR! cn_type needs to be specified when cn_thres is used; exec_detect_atoms')
            else if (params%cn_type .ne. 'cn_gen' .and. params%cn_type .ne. 'cn_std') then
              THROW_WARN('Unvalid cn_type, proceeding with generalised coordination number (cn_gen)')
              call cline%set('cn_type', 'cn_gen')
            endif
              call nano%identify_atomic_pos_thresh(params%thres, nint(params%cn_thres), params%cn_type)
          else
              call cline%set('cn_type', 'cn_nan')
              cn_thresh = 0
              call nano%identify_atomic_pos_thresh(params%thres, cn_thresh, params%cn_type)
          endif
        else
          if(cline%defined('cn_thres')) then
            if(.not. cline%defined('cn_type')) then
              THROW_HARD('ERROR! cn_type needs to be specified when cn_thres is used; exec_detect_atoms')
            elseif (params%cn_type .ne. 'cn_gen' .and. params%cn_type .ne. 'cn_std') then
              THROW_WARN('Unvalid cn_type, proceeding with generalised coordination number (cn_gen)')
              call cline%set('cn_type', 'cn_gen')
            endif
            call nano%identify_atomic_pos(nint(params%cn_thres), params%cn_type)
          else
            call cline%set('cn_type', 'cn_nan')
              cn_thresh = 0
              call nano%identify_atomic_pos(cn_thresh,params%cn_type)
          endif
        endif
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

    subroutine exec_dock_coords(self, cline)
      use simple_ori ! for generation of the rotation matrix
      use simple_atoms, only : atoms
      class(dock_coords_commander), intent(inout) :: self
      class(cmdline),               intent(inout) :: cline !< command line input
      type(parameters)       :: params
      type(atoms)            :: a_ref, a_targ
      real :: rot_trans(7) ! cost, rotation angles, shift
      real :: rmat(3,3)
      call params%new(cline)
      if(.not. cline%defined('pdbfile')) then
          THROW_HARD('ERROR! pdbfile needs to be present; exec_dock_coords')
      endif
      if(.not. cline%defined('pdbfile2')) then
          THROW_HARD('ERROR! pdbfile2 needs to be present; exec_dock_coords')
      endif
      if(.not. cline%defined('thres')) then
          THROW_HARD('ERROR! thres needs to be present; exec_dock_coords')
      endif
      call a_ref%new(basename(params%pdbfile))
      call a_targ%new(basename(params%pdbfile2))
      call dock_coords_init( a_ref, a_targ, params%thres )
      rot_trans = dock_coords_minimize()
      write(logfhandle, *) 'Docked coords have rmsd', rot_trans(1)
      write(logfhandle, *) 'Rotation angles        ', rot_trans(2:4)
      write(logfhandle, *) 'Translation vector     ', rot_trans(5:7)
      ! now perform rotation and translation of the coords
      ! and output rotshited coords
      rmat = euler2m(rot_trans(2:4))
      call a_targ%translate(-rot_trans(5:7))
      call a_targ%writePDB('TranslatedCoords')
      call a_targ%rotate(rmat)
      call a_targ%writePDB('DockedCoords')
      call a_ref%kill
      call a_targ%kill
    end subroutine exec_dock_coords

    ! Calculates distances distribution across the whole nanoparticle
    ! and radial dependent statistics.
    subroutine exec_atoms_stats( self, cline )
        class(atoms_stats_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle) :: nano
        real    :: min_rad, max_rad, step
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_atoms_stats')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_atoms_stats')
        endif
        if(cline%defined('cn_min') .and. .not. cline%defined('cn_max')) then
          THROW_HARD('ERROR! define cn_max too!; exec_atoms_stats')
        endif
        if(cline%defined('cn_max') .and. .not. cline%defined('cn_min')) then
          THROW_HARD('ERROR! define cn_min too!; exec_atoms_stats')
        endif
        min_rad = params%min_rad
        max_rad = params%max_rad
        step    = params%stepsz
        if(min_rad > max_rad) THROW_HARD('Minimum radius has to be smaller then maximum radius! exec_atoms_stats')
        if(abs(max_rad-min_rad)>TINY .and. step > max_rad-min_rad) THROW_HARD('Inputted too big stepsz! exec_atoms_stats')
        call nano%new(params%vols(1), params%smpd,params%element)
        ! execute
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
        call nano%set_img('../'//trim(fname)//'CC.mrc', 'img_cc')
        call nano%update_self_ncc()
        if(cline%defined('cn_min')) then
          if(cline%defined('cn')) then
            call nano%radial_dependent_stats(min_rad,max_rad,step,params%cn_min,params%cn_max,params%cn)
          else
            call nano%radial_dependent_stats(min_rad,max_rad,step,params%cn_min,params%cn_max)
          endif
        else
          if(cline%defined('cn')) then
            call nano%radial_dependent_stats(min_rad,max_rad,step,cn=params%cn)
          else
            call nano%radial_dependent_stats(min_rad,max_rad,step)
          endif
        endif
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_STATS NORMAL STOP ****')
    end subroutine exec_atoms_stats

    subroutine exec_atom_cluster_analysis( self, cline )
        class(atom_cluster_analysis_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle)     :: nano
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
        call nano%new(params%vols(1), params%smpd, params%element)
        ! execute
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords(trim(fname)//'_atom_centers.pdb')
        call nano%set_img(trim(fname)//'CC.mrc', 'img_cc')
        call nano%update_self_ncc()
        select case(trim(params%clustermode))
            case('ar')
              call nano%cluster_ar(params%thres)
            case('dist')
              call nano%cluster_interdist(params%thres)
            case('ang')
              call nano%cluster_ang(params%thres)
            case('maxint')
              call nano%atoms_stats(.false.)
              call nano%cluster_atom_maxint()
            case('intint')
              call nano%atoms_stats(.false.)
              call nano%cluster_atom_intint()
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
        class(cmdline),                 intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle) :: nano
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_nano_softmask')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_nano_softmask')
        endif
        call nano%new(params%vols(1), params%smpd,params%element)
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

    subroutine exec_atoms_mask( self, cline )
        use simple_nano_utils, only : atoms_mask
        class(atoms_mask_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline !< command line input
        type(parameters) :: params
        integer          :: nremoved
        call params%new(cline)
        if( .not. cline%defined('pdbfile') )then
            THROW_HARD('ERROR! pdbfile needs to be present; exec_atoms_mask')
        endif
        if( .not. cline%defined('pdbfile2') )then
            THROW_HARD('ERROR! pdbfile2 needs to be present; exec_atoms_mask')
        endif
        if( .not. cline%defined('max_rad') )then
            THROW_HARD('ERROR! max_rad needs to be present; exec_atoms_mask')
        endif
        ! execute
        call atoms_mask(params%pdbfile,params%max_rad,params%pdbfile2,nremoved)
        write(logfhandle,*) 'REMOVED ', nremoved, 'ATOMS FROM THE PDBFILE'
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_MASK NORMAL STOP ****')
    end subroutine exec_atoms_mask

    subroutine exec_geometry_analysis( self, cline )
        use simple_atoms, only : atoms
        class(geometry_analysis_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle)     :: nano
        type(atoms)            :: a
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_geometry_analysis')
        endif
        if( .not. cline%defined('pdbfile2') )then
            THROW_HARD('ERROR! pdbfile2 needs to be present; exec_geometry_analysis')
        endif
        if( .not. cline%defined('pdbfile2') )then
            THROW_HARD('ERROR! pdbfile2 needs to be present; exec_geometry_analysis')
        endif
        if( .not. cline%defined('pdbfile') )then
          if(.not. cline%defined('vol1')) THROW_HARD('ERROR! pdbfile or vol1 need to be present; exec_geometry_analysis')
        endif
        if(cline%defined('vol1')) then
          call nano%new(params%vols(1), params%smpd,params%element)
          ! fetch img_bin, img_cc and atomic positions
          fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
          call nano%set_img('../../'//trim(fname)//'BIN.mrc','img_bin')
          call nano%set_atomic_coords('../../'//trim(fname)//'_atom_centers.pdb')
          call nano%set_img('../../'//trim(fname)//'CC.mrc', 'img_cc')
          call nano%update_self_ncc()
          ! execute
          if(cline%defined('thres')) then
              call nano%geometry_analysis(trim(params%pdbfile2), params%thres)
          else
              call nano%geometry_analysis(trim(params%pdbfile2))
          endif
          ! kill
          call nano%kill
        elseif(cline%defined('pdbfile')) then
          call a%new(params%pdbfile)
          if(cline%defined('thres')) then
              call a%geometry_analysis_pdb(trim(params%pdbfile2), params%thres)
          else
              call a%geometry_analysis_pdb(trim(params%pdbfile2))
          endif
          call a%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_GEOMETRY_ANALYSIS NORMAL STOP ****')
    end subroutine exec_geometry_analysis

    subroutine exec_plot_atom( self, cline )
        class(plot_atom_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle)     :: nano
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_plot_atom')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_plot_atom')
        endif
        call nano%new(params%vols(1), params%smpd,params%element)
        ! fetch img_bin, img_cc and atomic positions
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
        call nano%set_img('../'//trim(fname)//'CC.mrc', 'img_cc')
        call nano%update_self_ncc()
        ! execute
        call nano%print_atoms()
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PLOT_ATOM NORMAL STOP ****')
      end subroutine exec_plot_atom


      subroutine exec_strain_analysis( self, cline )
          use simple_lattice_fitting, only : run_lattice_fit
          use simple_strain_mapping,  only : strain_analysis
          class(strain_analysis_commander), intent(inout) :: self
          class(cmdline),                   intent(inout) :: cline !< command line input
          type(parameters)       :: params
          real                   :: a(3) ! fitted lattice parameters
          real, allocatable      :: model(:,:)
          call cline%set('mkdir', 'yes')
          call params%new(cline)
          if( .not. cline%defined('pdbfile') )then
              THROW_HARD('ERROR! pdbfile needs to be present; exec_strain_analysis')
          endif
          call run_lattice_fit(params%pdbfile,model,a)
          call strain_analysis(model, a)
          ! end gracefully
          call simple_end('**** SIMPLE_STRAIN_ANALYSIS NORMAL STOP ****')
        end subroutine exec_strain_analysis
end module simple_commander_quant
