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
public :: tseries_atoms_analysis_commander
public :: dock_coords_commander
public :: atoms_mask_commander

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: detect_atoms_commander
  contains
    procedure :: execute      => exec_detect_atoms
end type detect_atoms_commander
type, extends(commander_base) :: atoms_stats_commander
  contains
    procedure :: execute      => exec_atoms_stats
end type atoms_stats_commander
type, extends(commander_base) :: tseries_atoms_analysis_commander
  contains
    procedure :: execute      => exec_tseries_atoms_analysis
end type tseries_atoms_analysis_commander
type, extends(commander_base) :: dock_coords_commander
  contains
    procedure :: execute      => exec_dock_coords
end type dock_coords_commander
type, extends(commander_base) :: atoms_mask_commander
  contains
    procedure :: execute      => exec_atoms_mask
end type atoms_mask_commander

integer, parameter :: CNMIN             = 5
integer, parameter :: CNMAX             = 12
integer, parameter :: CN_THRESH_DEFAULT = 5

contains

    ! Performs preprocessing on the nanoparticle map, identifies atomic positions,
    ! validates them and write them to disk
    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano
        real               :: a(3) ! lattice parameters
        logical            :: prefit_lattice, use_cs_thres, use_auto_corr_thres
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_detect_atoms')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_detect_atoms')
        endif
        prefit_lattice = cline%defined('vol2')
        call params%new(cline)
        use_cs_thres        = trim(params%use_thres) .eq. 'yes'
        use_auto_corr_thres = .not.cline%defined('corr_thres')
        if( prefit_lattice )then
            call nano%new(params%vols(2), params%smpd, params%element, params%msk)
            ! execute
            call nano%identify_lattice_params(a, use_auto_corr_thres=use_auto_corr_thres)
            ! kill
            call nano%kill
            call nano%new(params%vols(1), params%smpd, params%element, params%msk)
            ! execute
            if( cline%defined('cs_thres') )then
                call nano%identify_atomic_pos(a, l_fit_lattice=.false., use_cs_thres=use_cs_thres, use_auto_corr_thres=use_auto_corr_thres, cs_thres=params%cs_thres)
            else
                call nano%identify_atomic_pos(a, l_fit_lattice=.false., use_cs_thres=use_cs_thres, use_auto_corr_thres=use_auto_corr_thres)
            endif
            ! kill
            call nano%kill
        else
            call nano%new(params%vols(1), params%smpd, params%element, params%msk)
            ! execute
            if( cline%defined('cs_thres') )then
                call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres, use_auto_corr_thres=use_auto_corr_thres, cs_thres=params%cs_thres)
            else
                call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres, use_auto_corr_thres=use_auto_corr_thres)
            endif
            ! kill
            call nano%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

    subroutine exec_atoms_stats( self, cline )
        class(atoms_stats_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline !< command line input
        character(len=STDLEN) :: fname
        type(parameters)      :: params
        type(nanoparticle)    :: nano
        real                  :: a(3) ! lattice parameters
        logical               :: prefit_lattice, use_subset_coords, use_auto_corr_thres
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_atoms_stats')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 (raw map) needs to be present; exec_atoms_stats')
        endif
        if( .not. cline%defined('vol2') )then
            THROW_HARD('ERROR! vol2 (connected components map *CC.mrc) needs to be present; exec_atoms_stats')
        endif
        if( .not. cline%defined('pdbfile') )then
            THROW_HARD('ERROR! pdbfile needs to be present; exec_atoms_stats')
        endif
        prefit_lattice      = cline%defined('vol3')
        use_subset_coords   = cline%defined('pdbfile2')
        use_auto_corr_thres = .not.cline%defined('corr_thres')
        call params%new(cline)
        if( prefit_lattice )then
            ! fit lattice using vol3
            call nano%new(params%vols(3), params%smpd, params%element, params%msk)
            call nano%identify_lattice_params(a, use_auto_corr_thres=use_auto_corr_thres)
            call nano%kill
            ! calc stats
            call nano%new(params%vols(1), params%smpd, params%element, params%msk)
            call nano%set_atomic_coords(params%pdbfile)
            if( use_subset_coords ) call nano%set_coords4stats(params%pdbfile2)
            call nano%set_img(params%vols(2), 'img_cc')
            call nano%update_ncc()
            call nano%fillin_atominfo( a )
            call nano%write_csv_files
            call nano%kill
        else
            ! calc stats
            call nano%new(params%vols(1), params%smpd, params%element, params%msk)
            call nano%set_atomic_coords(params%pdbfile)
            if( use_subset_coords ) call nano%set_coords4stats(params%pdbfile2)
            call nano%set_img(params%vols(2), 'img_cc')
            call nano%update_ncc()
            call nano%fillin_atominfo
            call nano%write_csv_files
            call nano%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_STATS NORMAL STOP ****')
    end subroutine exec_atoms_stats

    subroutine exec_tseries_atoms_analysis( self, cline )
        use simple_nanoparticle_utils
        type :: common_atoms
            integer           :: ind1, ind2, ncommon
            real, allocatable :: coords1(:,:), coords2(:,:)
            real, allocatable :: common1(:,:), common2(:,:)
            real, allocatable :: different1(:,:), different2(:,:)
            real, allocatable :: displacements(:,:), dists(:)
        end type common_atoms
        class(tseries_atoms_analysis_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline !< command line input
        character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
        type(common_atoms),        allocatable :: atms_common(:,:)
        character(len=:),          allocatable :: fname1, fname2
        real, allocatable  :: pdbmat(:,:), dists_all(:)
        type(parameters)   :: params
        integer            :: npdbs, i, j, k, ndists, cnt
        character(len=2)   :: el
        type(stats_struct) :: dist_stats
        call params%new(cline)
        call read_filetable(params%pdbfiles, pdbfnames)
        npdbs = size(pdbfnames)
        allocate( atms_common(npdbs,npdbs) )
        el = trim(adjustl(params%element))
        ! identify common atoms across pairs
        do i = 1, npdbs - 1
            do j = i + 1, npdbs
                atms_common(i,j)%ind1 = i
                atms_common(i,j)%ind2 = j
                call read_pdb2matrix(trim(pdbfnames(i)), pdbmat)
                allocate(atms_common(i,j)%coords1(3,size(pdbmat,dim=2)), source=pdbmat)
                deallocate(pdbmat)
                call read_pdb2matrix(trim(pdbfnames(j)), pdbmat)
                allocate(atms_common(i,j)%coords2(3,size(pdbmat,dim=2)), source=pdbmat)
                deallocate(pdbmat)
                ! identify couples
                call find_couples( atms_common(i,j)%coords1, atms_common(i,j)%coords2, el,&
                                  &atms_common(i,j)%common1, atms_common(i,j)%common2)
                atms_common(i,j)%ncommon = size(atms_common(i,j)%common1, dim=2)
                ! calculate displacements and distances
                allocate( atms_common(i,j)%displacements(3,atms_common(i,j)%ncommon), atms_common(i,j)%dists(atms_common(i,j)%ncommon) )
                do k = 1, atms_common(i,j)%ncommon
                    atms_common(i,j)%displacements(:,k) = atms_common(i,j)%common2(:,k) - atms_common(i,j)%common1(:,k)
                    atms_common(i,j)%dists(k) = sqrt(sum((atms_common(i,j)%displacements(:,k))**2.))
                end do
                ! write PDB files
                if( j == i + 1 )then
                    fname1 = 'common_atoms_'//int2str_pad(i,2)//'in'//int2str_pad(j,2)//'.pdb'
                    fname2 = 'common_atoms_'//int2str_pad(j,2)//'in'//int2str_pad(i,2)//'.pdb'
                    call write_matrix2pdb(el, atms_common(i,j)%common1, fname1)
                    call write_matrix2pdb(el, atms_common(i,j)%common2, fname2)
                endif
            end do
        end do
        ! calculate distance statistics
        ndists = 0
        do i = 1, npdbs - 1
            do j = i + 1, npdbs
                ndists = ndists + atms_common(i,j)%ncommon
            end do
        end do
        allocate( dists_all(ndists), source=0. )
        cnt = 0
        do i = 1, npdbs - 1
            do j = i + 1, npdbs
                do k = 1, atms_common(i,j)%ncommon
                    cnt = cnt + 1
                    dists_all(cnt) = atms_common(i,j)%dists(k)
                end do
            end do
        end do
        call calc_stats(dists_all, dist_stats)
        write(logfhandle,'(A)') '>>> DISTANCE STATS FOR COMMON ATOMS BELOW'
        write(logfhandle,'(A,F8.4)') 'Average: ', dist_stats%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', dist_stats%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', dist_stats%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', dist_stats%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', dist_stats%minv
        ! identify different atoms across pairs
        do i = 1, npdbs - 1
            do j = i + 1, npdbs
                call remove_atoms( atms_common(i,j)%common1, atms_common(i,j)%coords1, atms_common(i,j)%different1 )
                call remove_atoms( atms_common(i,j)%common2, atms_common(i,j)%coords2, atms_common(i,j)%different2 )
                ! write PDB files
                if( j == i + 1 )then
                    fname1 = 'different_atoms_'//int2str_pad(i,2)//'not_in'//int2str_pad(j,2)//'.pdb'
                    fname2 = 'different_atoms_'//int2str_pad(j,2)//'not_in'//int2str_pad(i,2)//'.pdb'
                    call write_matrix2pdb(el, atms_common(i,j)%different1, fname1)
                    call write_matrix2pdb(el, atms_common(i,j)%different2, fname2)
                endif
            end do
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_ATOMS_ANALYSIS NORMAL STOP ****')
    end subroutine exec_tseries_atoms_analysis

    subroutine exec_dock_coords( self, cline )
        use simple_ori ! for generation of the rotation matrix
        use simple_atoms, only: atoms
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

    subroutine exec_atoms_mask( self, cline )
        use simple_nanoparticle_utils, only: atoms_mask
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

end module simple_commander_quant
