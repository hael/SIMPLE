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
public :: atoms_rmsd_commander
public :: radial_dependent_stats_commander
public :: atom_cluster_analysis_commander
public :: nano_softmask_commander
public :: geometry_analysis_commander
public :: radial_sym_test_commander
public :: plot_atom_commander
public :: dock_coords_commander

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
type, extends(commander_base) :: radial_sym_test_commander
  contains
    procedure :: execute      => exec_radial_sym_test
end type radial_sym_test_commander
type, extends(commander_base) :: plot_atom_commander
  contains
    procedure :: execute      => exec_plot_atom
end type plot_atom_commander

contains

    ! Performs preprocessing on the nanoparticle map, identifies atomic positions,
    ! validates them and write them to disk
    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_detect_atoms')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_detect_atoms')
        endif
        call params%new(cline)
        call nano%new(params%vols(1), params%smpd, params%element)
        ! volume soft-edge masking
        if(cline%defined('msk')) call nano%mask(params%msk)
        ! execute
        if(cline%defined('thres')) then
          if(cline%defined('cs_thres')) then ! contact score threshold for outliers removal
              call nano%identify_atomic_pos_thresh(params%thres, nint(params%cs_thres))
          else
              call nano%identify_atomic_pos_thresh(params%thres)
          endif
        else
          if(cline%defined('cs_thres')) then ! contact score threshold for outliers removal
              call nano%identify_atomic_pos(nint(params%cs_thres))
          else
              call nano%identify_atomic_pos()
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

    ! Example of cline:
    ! simple_exec prg=dock_volpair lpstart=3 lpstop=1 smpd=0.358 vol1=vol_Nov26.mrc
    ! vol2=vol_Nov28_2.mrc  msk=40 nthr=8
    subroutine exec_atoms_rmsd( self, cline )
        use simple_ori,              only: ori
        use simple_atoms,            only: atoms
        use simple_nano_utils,       only: kabsch, read_pdb2matrix, write_matrix2pdb, find_couples
        class(atoms_rmsd_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname1, fname2
        type(parameters)       :: params
        type(nanoparticle)     :: nano1, nano2
        integer, parameter     :: D = 3 ! dimension of the space (3D)
        real(dp), allocatable  :: points_P(:,:), points_Q(:,:), P(:,:), Q(:,:)
        type(ori)   :: orientation
        type(oris)  :: ori2read
        real        :: smpd, cxyz(D), mat(D,D), shift(D)
        real(dp)    :: U(D,D), r(D), lrms
        integer     :: ldim(3), nptcls, i
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
        if(.not. cline%defined('element')) then
          THROW_HARD('ERROR! element needs to be present; exec_atoms_rmsd')
        endif
        call nano1%new(params%vols(1), params%smpd,params%element)
        call nano2%new(params%vols(2), params%smpd,params%element)
        fname1 = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        fname2 = get_fbody(trim(basename(params%vols(2))), trim(fname2ext(params%vols(2))))
        call nano1%set_atomic_coords(trim(fname1)//'_atom_centers.pdb')
        ! execute
        call nano2%set_atomic_coords(trim(fname2)//'_atom_centers.pdb')
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
        if(min_rad > max_rad) THROW_HARD('Minimum radius has to be smaller then maximum radius! exec_radial_dependent_stats')
        if(step > max_rad-min_rad) THROW_HARD('Inputted too big stepsz! exec_radial_dependent_stats')
        call nano%new(params%vols(1), params%smpd,params%element)
        ! execute
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords(trim(fname)//'_atom_centers.pdb')
        call nano%set_img(trim(fname)//'CC.mrc', 'img_cc')
        call nano%update_self_ncc()
        call nano%radial_dependent_stats(min_rad,max_rad,step)
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
              call nano%atom_intensity_stats(.false.)
              call nano%cluster_atom_maxint()
            case('intint')
              call nano%atom_intensity_stats(.false.)
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
        real  :: smpd
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

    subroutine exec_geometry_analysis( self, cline )
        use simple_atoms, only : atoms
        class(geometry_analysis_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline !< command line input
        character(len=STDLEN)  :: fname
        type(parameters)       :: params
        type(nanoparticle)     :: nano
        type(atoms)            :: a
        real    :: smpd
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
          call nano%set_img('../'//trim(fname)//'BIN.mrc','img_bin')
          call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
          call nano%set_img('../'//trim(fname)//'CC.mrc', 'img_cc')
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

    subroutine exec_radial_sym_test( self, cline )
        use simple_nanoparticle,    only : nanoparticle
        use simple_atoms,           only : atoms
        use simple_commander_volops, only: symmetry_test_commander
        class(radial_sym_test_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(symmetry_test_commander) :: symtstcmd
        type(parameters)       :: params
        character(len=STDLEN)  :: fname
        character(len=3)       :: pgrp
        type(nanoparticle) :: nano
        type(atoms) :: atom
        type(image) :: simulated_distrib
        real        :: cutoff
        real        :: min_rad, max_rad, step
        integer     :: i, ldim(3)
        real        :: radius
        character(len=100) :: fname_conv
        character(len=100) :: out_pdbfile
        if( .not. cline%defined('lp'))    call cline%set('lp', 1.)
        if( .not. cline%defined('cenlp')) call cline%set('cenlp', 5.)
        if( .not. cline%defined('mkdir')) call cline%set('mkdir','yes')
        if( .not. cline%defined('center'))call cline%set('center','yes')
        call params%new(cline)
        call nano%new(params%vols(1), params%smpd,params%element)
        !identifiy atomic positions
        fname = get_fbody(trim(basename(params%vols(1))), trim(fname2ext(params%vols(1))))
        call nano%set_atomic_coords('../'//trim(fname)//'_atom_centers.pdb')
        !Translate so that the center of mass coincide with its closest atom, WE DECIDED NOT TO DO IT
        ! out_pdbfile = '../'//trim(fname)//'_atom_centeredOnAtom'
        ! call nano%center_on_atom('../'//trim(fname)//'_atom_centers.pdb', out_pdbfile )
        ! call nano%set_atomic_coords('../'//trim(fname)//'_atom_centeredOnAtom.pdb')
        call nano%get_ldim(ldim)
        min_rad = params%min_rad
        max_rad = params%max_rad
        step    = params%stepsz
        if(min_rad > max_rad) THROW_HARD('Minimum radius has to be smaller then maximum radius! exec_radial_sym_test')
        if(step > max_rad-min_rad) THROW_HARD('Inputted too big stepsz! exec_radial_sym_test')
        cutoff = 8.*params%smpd
        call simulated_distrib%new(ldim,params%smpd)
        do i =1, nint((max_rad-min_rad)/step)
            radius = min_rad+real(i-1)*step
            if(radius <= max_rad) then
              call cline%set('msk', radius/params%smpd+3.) !+3 to be sure
              if(i > 1) call cline%set('mkdir','no')
              fname_conv = 'atomic_coords_'//trim(int2str(nint(radius))//'A')
              call nano%keep_atomic_pos_at_radius(radius, params%element, fname_conv)
              ! Generate distribution based on atomic position
              call atom%new('atomic_coords_'//trim(int2str(nint(radius))//'A.pdb'))
              call atom%convolve(simulated_distrib, cutoff)
              call simulated_distrib%write('density_'//trim(int2str(nint(radius))//'A.mrc'))
              ! Prepare for calling exec_symmetry_test
              call cline%set('vol1','density_'//trim(int2str(nint(radius))//'A.mrc'))
              call cline%set('fname','sym_test_'//int2str(nint(radius))//'A.txt')
              call symtstcmd%execute(cline)
              if(i == 1) then
                  call del_file('../'//'atomic_coords_'//trim(int2str(nint(radius))//'A.pdb'))
                  call del_file('../'//'density_'//trim(int2str(nint(radius))//'A.mrc'))
              else
                  call del_file('atomic_coords_'//trim(int2str(nint(radius))//'A.pdb'))
                  call del_file('density_'//trim(int2str(nint(radius))//'A.mrc'))
              endif
              call atom%kill
            endif
        enddo
        ! end gracefully
        call simple_end('**** SIMPLE_RADIAL_SYM_TEST NORMAL STOP ****')
    end subroutine exec_radial_sym_test

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
end module simple_commander_quant
