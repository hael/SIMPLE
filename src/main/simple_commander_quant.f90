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
use simple_nanoparticle_utils
implicit none

public :: detect_atoms_commander
public :: detect_atoms_eo_commander
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

type, extends(commander_base) :: detect_atoms_eo_commander
  contains
    procedure :: execute      => exec_detect_atoms_eo
end type detect_atoms_eo_commander

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

type :: common_atoms
    integer           :: ind1, ind2, ncommon
    real, allocatable :: coords1(:,:), coords2(:,:)
    real, allocatable :: common1(:,:), common2(:,:)
    real, allocatable :: different1(:,:), different2(:,:)
    real, allocatable :: displacements(:,:), dists(:)
end type common_atoms

contains

    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)   :: params
        type(nanoparticle) :: nano
        real               :: a(3) ! lattice parameters
        logical            :: prefit_lattice, use_cs_thres, use_auto_corr_thres
        prefit_lattice = cline%defined('vol2')
        call params%new(cline)
        use_cs_thres        = trim(params%use_thres) .eq. 'yes'
        use_auto_corr_thres = .not.cline%defined('corr_thres')
        if( prefit_lattice )then
            call nano%new(params%vols(2), params%msk)
            ! execute
            call nano%identify_lattice_params(a, use_auto_corr_thres=use_auto_corr_thres)
            ! kill
            call nano%kill
            call nano%new(params%vols(1), params%msk)
            ! execute
            if( cline%defined('cs_thres') )then
                call nano%identify_atomic_pos(a, l_fit_lattice=.false., use_cs_thres=use_cs_thres,&
                &use_auto_corr_thres=use_auto_corr_thres, cs_thres=params%cs_thres)
            else
                call nano%identify_atomic_pos(a, l_fit_lattice=.false., use_cs_thres=use_cs_thres,&
                &use_auto_corr_thres=use_auto_corr_thres)
            endif
            ! kill
            call nano%kill
        else
            call nano%new(params%vols(1), params%msk)
            ! execute
            if( cline%defined('cs_thres') )then
                call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
                &use_auto_corr_thres=use_auto_corr_thres, cs_thres=params%cs_thres)
            else
                call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
                &use_auto_corr_thres=use_auto_corr_thres)
            endif
            ! kill
            call nano%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

    subroutine exec_detect_atoms_eo( self, cline )
        use simple_opt_filter, only: opt_filter
        class(detect_atoms_eo_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters)   :: params
        type(nanoparticle) :: nano
        type(image)        :: even, odd, mskvol, sim_density
        type(common_atoms) :: atms_common
        real, allocatable  :: pdbmat(:,:)
        character(len=2)   :: el
        integer            :: k
        type(stats_struct) :: dist_stats
        character(len=:), allocatable :: add2fn, ename_filt, oname_filt, eatms, oatms
        character(len=:), allocatable :: fname_avg, map_avg_filt, tmp, eatms_common
        character(len=:), allocatable :: oatms_common, oatms_sim, eatms_sim, atms_avg, atms_avg_sim
        call cline%set('use_thres', 'no')
        call cline%set('corr_thres', 0.)
        call cline%set('lp_lb',      3.)
        call params%new(cline)
        ! read e/o:s
        call odd %new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols_odd(1))
        call even%read(params%vols_even(1))
        ! spherical masking
        call even%mask(params%msk, 'soft')
        call odd%mask(params%msk, 'soft')
        call mskvol%disc([params%box,params%box,params%box], params%smpd,&
        &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        ! nonuniform filter
        ! ... nonuniform=.true., smooth_ext=1, filter_type='butterworth', lp_lb=3. ...
        call opt_filter(odd, even, mskvol)
        call even%mask(params%msk, 'soft')
        call odd%mask(params%msk, 'soft')
        add2fn       = '_filt'
        oname_filt   = add2fbody(params%vols_odd(1),     trim(params%ext), add2fn)
        ename_filt   = add2fbody(params%vols_even(1),    trim(params%ext), add2fn)
        tmp          = rm_from_fbody(params%vols_odd(1), trim(params%ext), '_odd')
        map_avg_filt = add2fbody(tmp,                    trim(params%ext), '_filt_AVG')
        atms_avg_sim = add2fbody(tmp,                    trim(params%ext), '_ATMS_AVG_SIM')
        fname_avg    = swap_suffix(tmp, '.pdb',          trim(params%ext) )
        atms_avg     = add2fbody(fname_avg ,             '.pdb',           '_ATMS_AVG')
        call odd%write(oname_filt)
        call even%write(ename_filt)
        ! detect atoms in odd
        call nano%new(oname_filt)
        call nano%identify_atomic_pos_slim(.false.)
        ! detect atoms in even
        call nano%new(ename_filt)
        call nano%identify_atomic_pos_slim(.true.)
        ! compare independent atomic models
        el               = trim(adjustl(params%element))
        tmp              = add2fbody(oname_filt,    trim(params%ext), '_ATMS')
        oatms            = swap_suffix(tmp, '.pdb', trim(params%ext) )
        tmp              = add2fbody(ename_filt,    trim(params%ext), '_ATMS')
        eatms            = swap_suffix(tmp, '.pdb', trim(params%ext) )
        tmp              = add2fbody(oname_filt,    trim(params%ext), '_ATMS_COMMON')
        oatms_sim        = add2fbody(oname_filt,    trim(params%ext), '_ATMS_COMMON_SIM')
        oatms_common     = swap_suffix(tmp, '.pdb', trim(params%ext) )
        tmp              = add2fbody(ename_filt,    trim(params%ext), '_ATMS_COMMON')
        eatms_sim        = add2fbody(ename_filt,    trim(params%ext), '_ATMS_COMMON_SIM')
        eatms_common     = swap_suffix(tmp, '.pdb', trim(params%ext) )
        atms_common%ind1 = 1
        atms_common%ind2 = 2
        call read_pdb2matrix(oatms, pdbmat)
        allocate(atms_common%coords1(3,size(pdbmat,dim=2)), source=pdbmat)
        deallocate(pdbmat)
        call read_pdb2matrix(eatms, pdbmat)
        allocate(atms_common%coords2(3,size(pdbmat,dim=2)), source=pdbmat)
        deallocate(pdbmat)
        ! identify couples
        call find_couples( atms_common%coords1, atms_common%coords2, el,&
        &atms_common%common1, atms_common%common2)
        atms_common%ncommon = size(atms_common%common1, dim=2)
        ! calculate displacements and distances
        allocate( atms_common%displacements(3,atms_common%ncommon), atms_common%dists(atms_common%ncommon) )
        do k = 1, atms_common%ncommon
            atms_common%displacements(:,k) = atms_common%common2(:,k) - atms_common%common1(:,k)
            atms_common%dists(k) = sqrt(sum((atms_common%displacements(:,k))**2.))
        end do
        ! write pdb files for the e/o:s
        call write_matrix2pdb(el, atms_common%common1, oatms_common)
        call write_matrix2pdb(el, atms_common%common2, eatms_common)
        ! write the average atomic positions
        atms_common%common1 = (atms_common%common1 + atms_common%common2) / 2.
        call write_matrix2pdb(el, atms_common%common1, atms_avg)
        ! RMSD reporting
        call calc_stats(atms_common%dists(:), dist_stats)
        write(logfhandle,'(A)') '>>> DISTANCE STATS (IN A) FOR COMMON ATOMS BELOW'
        write(logfhandle,'(A,F8.4)') 'Average: ', dist_stats%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', dist_stats%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', dist_stats%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', dist_stats%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', dist_stats%minv
        ! generate simulated densities for the e/o pairs
        call sim_density%new([params%box,params%box,params%box], params%smpd)
        call nano%set_atomic_coords(oatms_common)
        call nano%simulate_atoms(sim_density)
        call sim_density%write(oatms_sim)
        call nano%set_atomic_coords(eatms_common)
        call nano%simulate_atoms(sim_density)
        call sim_density%write(eatms_sim)
        call nano%set_atomic_coords(oatms_common)
        ! simulate density for the average atomic positions
        call nano%set_atomic_coords(atms_avg)
        call nano%simulate_atoms(sim_density)
        call sim_density%write(atms_avg_sim)
        ! write average filtered map
        call even%add(odd)
        call even%mul(0.5)
        call even%write(map_avg_filt)
        ! kill
        call nano%kill
        call even%kill
        call odd%kill
        call mskvol%kill
        call sim_density%kill
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS_EO NORMAL STOP ****')
    end subroutine exec_detect_atoms_eo

    subroutine exec_atoms_stats( self, cline )
        class(atoms_stats_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline !< command line input
        character(len=STDLEN) :: fname
        type(parameters)      :: params
        type(nanoparticle)    :: nano
        real                  :: a(3) ! lattice parameters
        logical               :: prefit_lattice, use_subset_coords, use_auto_corr_thres
        prefit_lattice      = cline%defined('vol3')
        use_subset_coords   = cline%defined('pdbfile2')
        use_auto_corr_thres = .not.cline%defined('corr_thres')
        call params%new(cline)
        if( prefit_lattice )then
            ! fit lattice using vol3
            call nano%new(params%vols(3), params%msk)
            call nano%identify_lattice_params(a, use_auto_corr_thres=use_auto_corr_thres)
            call nano%kill
            ! calc stats
            call nano%new(params%vols(1), params%msk)
            call nano%set_atomic_coords(params%pdbfile)
            if( use_subset_coords ) call nano%set_coords4stats(params%pdbfile2)
            call nano%set_img(params%vols(2), 'img_cc')
            call nano%update_ncc()
            call nano%fillin_atominfo( a )
            call nano%write_csv_files
            call nano%kill
        else
            ! calc stats
            call nano%new(params%vols(1), params%msk)
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
        class(tseries_atoms_analysis_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline !< command line input
        character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
        type(common_atoms),        allocatable :: atms_common(:,:)
        character(len=:),          allocatable :: fname1, fname2
        real, allocatable  :: pdbmat(:,:), dists_all(:)
        type(parameters)   :: params
        integer            :: npdbs, i, j, k, ndists, cnt, ipdb
        character(len=2)   :: el
        type(stats_struct) :: dist_stats
        call params%new(cline)
        call read_filetable(params%pdbfiles, pdbfnames)
        npdbs = size(pdbfnames)
        if( params%mkdir.eq.'yes' )then
            do ipdb = 1,npdbs
                if(pdbfnames(ipdb)(1:1).ne.'/') pdbfnames(ipdb) = '../'//trim(pdbfnames(ipdb))
            enddo
        endif
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
        class(atoms_mask_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline !< command line input
        type(parameters) :: params
        integer          :: nremoved
        call params%new(cline)
        ! execute
        call atoms_mask(params%pdbfile,params%max_rad,params%pdbfile2,nremoved)
        write(logfhandle,*) 'REMOVED ', nremoved, 'ATOMS FROM THE PDBFILE'
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_MASK NORMAL STOP ****')
    end subroutine exec_atoms_mask

end module simple_commander_quant
