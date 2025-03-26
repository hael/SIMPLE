module simple_commander_atoms
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_image,          only: image
use simple_binimage,       only: binimage
use simple_nanoparticle,   only: nanoparticle
use simple_nanoparticle_utils
use simple_atoms,          only: atoms
implicit none

public :: atoms_stats_commander
public :: conv_atom_denoise_commander
public :: detect_atoms_commander
public :: map_validation_commander
public :: model_validation_commander
public :: model_validation_eo_commander
public :: pdb2mrc_commander
public :: tseries_atoms_rmsd_commander
public :: tseries_core_atoms_analysis_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: atoms_stats_commander
  contains
    procedure :: execute      => exec_atoms_stats
end type atoms_stats_commander

type, extends(commander_base) :: conv_atom_denoise_commander
  contains
    procedure :: execute      => exec_conv_atom_denoise
end type conv_atom_denoise_commander

type, extends(commander_base) :: detect_atoms_commander
  contains
    procedure :: execute      => exec_detect_atoms
end type detect_atoms_commander

type, extends(commander_base) :: map_validation_commander
  contains
    procedure :: execute      => exec_map_validation
end type map_validation_commander

type, extends(commander_base) :: model_validation_commander
  contains
    procedure :: execute      => exec_model_validation
end type model_validation_commander

type, extends(commander_base) :: model_validation_eo_commander
  contains
    procedure :: execute      => exec_model_validation_eo
end type model_validation_eo_commander

type, extends(commander_base) :: pdb2mrc_commander
  contains
    procedure :: execute      => exec_pdb2mrc
end type pdb2mrc_commander

type, extends(commander_base) :: tseries_atoms_rmsd_commander
  contains
    procedure :: execute      => exec_tseries_atoms_rmsd
end type tseries_atoms_rmsd_commander

type, extends(commander_base) :: tseries_core_atoms_analysis_commander
  contains
    procedure :: execute      => exec_tseries_core_atoms_analysis
end type tseries_core_atoms_analysis_commander

type :: common_atoms
    integer           :: ind1, ind2, ncommon, ndiff
    real, allocatable :: coords1(:,:), coords2(:,:)
    real, allocatable :: common1(:,:), common2(:,:)
    real, allocatable :: different1(:,:), different2(:,:)
    real, allocatable :: displacements(:,:), dists(:)
end type common_atoms

contains

    subroutine exec_atoms_stats( self, cline )
        class(atoms_stats_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline !< command line input
        type(parameters)      :: params
        type(nanoparticle)    :: nano
        real                  :: a(3) ! lattice parameters
        logical               :: prefit_lattice, use_subset_coords
        prefit_lattice    = cline%defined('vol3')
        use_subset_coords = cline%defined('pdbfile2')
        call params%new(cline)
        if( prefit_lattice )then
            ! fit lattice using vol3
            call nano%new(params%vols(3))
            call nano%identify_lattice_params(a)
            call nano%kill
            ! calc stats
            call nano%new(params%vols(1))
            call nano%set_atomic_coords(params%pdbfile)
            if( use_subset_coords ) call nano%set_coords4stats(params%pdbfile2)
            call nano%set_img(params%vols(2), 'img_cc')
            call nano%update_ncc()
            call nano%fillin_atominfo(a)
            call nano%write_csv_files
            call nano%kill
        else
            ! calc stats
            call nano%new(params%vols(1))
            call nano%set_atomic_coords(params%pdbfile)
            if( use_subset_coords ) call nano%set_coords4stats(params%pdbfile2)
            call nano%set_img(params%vols(2), 'img_cc')
            call nano%update_ncc()
            call nano%fillin_atominfo()
            call nano%write_csv_files
            call nano%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_STATS NORMAL STOP ****')
    end subroutine exec_atoms_stats

    subroutine exec_conv_atom_denoise( self, cline )
        class(conv_atom_denoise_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)   :: params
        type(nanoparticle) :: nano
        if( .not. cline%defined('outvol') ) call cline%set('outvol', 'denoised.mrc')
        call params%new(cline)
        call nano%new(params%vols(1), params%msk)
        call nano%conv_denoise(params%outvol)
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_CONV_ATOM_DENOISE NORMAL STOP ****')
    end subroutine exec_conv_atom_denoise

    subroutine exec_detect_atoms( self, cline )
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)   :: params
        type(nanoparticle) :: nano
        real               :: a(3) ! lattice parameters
        logical            :: prefit_lattice
        prefit_lattice = cline%defined('vol2')
        call params%new(cline)
        if( prefit_lattice )then
            call nano%new(params%vols(2))
            ! execute
            call nano%identify_lattice_params(a)
            ! kill
            call nano%kill
            call nano%new(params%vols(1))
            ! execute
            call nano%identify_atomic_pos(a, l_fit_lattice=.false., l_atom_thres=trim(params%atom_thres).eq.'yes')
            ! kill
            call nano%kill
        else
            call nano%new(params%vols(1))
            ! execute
            call nano%identify_atomic_pos(a, l_fit_lattice=.true., l_atom_thres=trim(params%atom_thres).eq.'yes')
            ! kill
            call nano%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

    subroutine exec_map_validation( self, cline )
        class(map_validation_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        type(image)      :: exp_vol, sim_vol
        integer          :: ldim(3), ldim_new(3), ifoo, box, box_new
        real             :: upscaling_factor, smpd_new
        character(len=STDLEN), allocatable :: sim_vol_file, pdbout, upscale_vol_file
        call params%new(cline)
        upscale_vol_file = trim(get_fbody(params%vols(1),'mrc'))//'_upscale.mrc'
        call find_ldim_nptcls(params%vols(1), ldim, ifoo)
        write(logfhandle,'(a,3i6,a,f8.3,a)') 'Original dimensions (', ldim,' ) voxels, smpd: ', params%smpd, ' Angstrom'
        box              = ldim(1)
        upscaling_factor = params%smpd / params%smpd_target
        box_new          = round2even(real(ldim(1)) * upscaling_factor)
        ldim_new(:)      = box_new
        upscaling_factor = real(box_new) / real(box)
        smpd_new         = params%smpd / upscaling_factor
        write(logfhandle,'(a,3i6,a,f8.3,a)') 'Scaled dimensions   (', ldim_new,' ) voxels, smpd: ', smpd_new, ' Angstrom'
        call molecule%new(params%pdbfile)
        sim_vol_file     = trim(get_fbody(params%pdbfile,'pdb'))//'_sim_vol.mrc'
        pdbout           = trim(get_fbody(params%pdbfile,'pdb'))//'_centered.pdb'
        call molecule%pdb2mrc( params%pdbfile, sim_vol_file, smpd_new, pdb_out=pdbout, vol_dim=ldim_new)
        call sim_vol%new([box_new, box_new, box_new], smpd_new)
        call sim_vol%read(sim_vol_file)
        call exp_vol%read_and_crop(params%vols(1), params%smpd, box_new, smpd_new)
        call exp_vol%write(upscale_vol_file)
        call molecule%map_validation(exp_vol, sim_vol, filename='map_val_coor')
        call molecule%kill()
        call exp_vol%kill()
        call sim_vol%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_VALIDATION NORMAL STOP ****')
    end subroutine exec_map_validation

    subroutine exec_model_validation( self, cline )
        class(model_validation_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        call params%new(cline)
        call molecule%new(params%pdbfile)
        call molecule%model_validation(params%pdbfile, params%vols(1), params%smpd, params%smpd_target)
        call molecule%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_MODEL_VALIDATION NORMAL STOP ****')
    end subroutine exec_model_validation

    subroutine exec_model_validation_eo( self, cline )
        class(model_validation_eo_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        call params%new(cline)
        call molecule%new(params%pdbfile)
        if( .not.cline%defined('vol2') ) params%vols(2) = trim(get_fbody(params%vols(1),'mrc'))//'_even.mrc'
        if( .not.cline%defined('vol3') ) params%vols(3) = trim(get_fbody(params%vols(1),'mrc'))//'_odd.mrc'
        call molecule%model_validation_eo(params%pdbfile, params%vols(1), params%vols(2), params%vols(3), params%smpd, params%smpd_target)
        call molecule%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_MODEL_VALIDATION_EO NORMAL STOP ****')
    end subroutine exec_model_validation_eo

    subroutine exec_pdb2mrc( self, cline )
        class(pdb2mrc_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        integer          :: vol_dims(3)
        call params%new(cline)
        call molecule%new(params%pdbfile)
        if( .not.cline%defined('pdbout') )then
            params%pdbout = trim(get_fbody(params%pdbfile,'pdb'))//'_centered.pdb'
        endif
        if( .not.cline%defined('outvol') ) params%outvol = swap_suffix(params%pdbfile,'mrc','pdb')
        if( cline%defined('vol_dim') ) vol_dims(:) = params%vol_dim
        if( cline%defined('vol_dim') )then  
            if( params%center_pdb .eq. 'yes' )then
                call molecule%pdb2mrc(params%pdbfile, params%outvol, params%smpd, center_pdb=.true., pdb_out=params%pdbout, vol_dim=vol_dims)
            else
                call molecule%pdb2mrc(params%pdbfile, params%outvol, params%smpd, pdb_out=params%pdbout, vol_dim=vol_dims)
            endif
        else
            if( params%center_pdb .eq. 'yes' )then
                call molecule%pdb2mrc(params%pdbfile, params%outvol, params%smpd, center_pdb=.true., pdb_out=params%pdbout)
            else
                call molecule%pdb2mrc(params%pdbfile, params%outvol, params%smpd, pdb_out=params%pdbout)
            endif
        endif
        call molecule%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PDB2MRC NORMAL STOP ****')
    end subroutine exec_pdb2mrc

    subroutine exec_tseries_atoms_rmsd( self, cline )
        class(tseries_atoms_rmsd_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline !< command line input
        character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
        real,               allocatable :: pdbmat1(:,:), pdbmat2(:,:), dists(:), distmat(:,:), dists_neigh(:)
        integer,            allocatable :: inds(:), pairs(:,:), inds_neigh(:)
        type(stats_struct), allocatable :: dist_stats(:)
        type(parameters) :: params
        integer          :: npdbs, i, ii, j, ind, ipdb, npairs, cnt
        character(len=2) :: el
        call cline%set('mkdir', 'no')
        call params%new(cline)
        call read_filetable(params%pdbfiles, pdbfnames)
        npdbs = size(pdbfnames)
        if( params%mkdir.eq.'yes' )then
            do ipdb = 1,npdbs
                if(pdbfnames(ipdb)(1:1).ne.'/') pdbfnames(ipdb) = '../'//trim(pdbfnames(ipdb))
            enddo
        endif
        el     = trim(adjustl(params%element))
        npairs = (npdbs * (npdbs - 1))/2
        allocate(dists(npairs), distmat(npdbs,npdbs), dists_neigh(npdbs), source=0.)
        allocate(inds(npairs),  pairs(npairs,2), inds_neigh(npdbs), source=0 )
        allocate(dist_stats(npairs))
        ! identify common atoms across pairs
        cnt = 0
        do i = 1, npdbs - 1
            do j = i + 1, npdbs
                cnt = cnt + 1
                call read_pdb2matrix(trim(pdbfnames(i)), pdbmat1)
                call read_pdb2matrix(trim(pdbfnames(j)), pdbmat2)
                dist_stats(cnt) = atm_rmsd_stats(pdbmat1, pdbmat2)
                dists(cnt)      = dist_stats(cnt)%med
                distmat(i,j)    = dist_stats(cnt)%med
                distmat(j,i)    = dist_stats(cnt)%med
                inds(cnt)       = cnt
                pairs(cnt,:)    = [i,j]      
                deallocate(pdbmat1, pdbmat2)
            end do
        end do
        ! identify nearest neighbors
        do i = 1, npdbs
            do j = 1, npdbs
                inds_neigh(j)  = j
                dists_neigh(j) = distmat(i,j) 
            end do
            call hpsort(dists_neigh, inds_neigh)
            write(logfhandle,"(A,*(G0,:,','))") 'NEAREST NEIGBORS: ', inds_neigh
            write(logfhandle,"(A,*(G0,:,','))") 'NEAREST DISTANCE: ', dists_neigh
        end do
        deallocate(dist_stats, dists, dists_neigh, inds, inds_neigh, pairs)
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_ATOMS_RMSD NORMAL STOP ****')
    end subroutine exec_tseries_atoms_rmsd

    subroutine exec_tseries_core_atoms_analysis( self, cline )
        class(tseries_core_atoms_analysis_commander), intent(inout) :: self
        class(cmdline),                               intent(inout) :: cline !< command line input
        character(len=LONGSTRLEN), allocatable :: pdbfnames(:), pdbfnames_core(:), pdbfnames_bfac(:), pdbfnames_fringe(:)
        type(common_atoms),        allocatable :: atms_common(:)
        character(len=:),          allocatable :: fname1, fname2
        real,                      allocatable :: betas(:), pdbmat(:,:)
        logical,                   allocatable :: mask_core(:)
        type(parameters)   :: params
        integer            :: npdbs, i, j, k, cnt, ipdb, natoms
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
        allocate(pdbfnames_core(npdbs), pdbfnames_bfac(npdbs), pdbfnames_fringe(npdbs))
        do ipdb = 1,npdbs
            pdbfnames_core(ipdb)   = add2fbody(basename(pdbfnames(ipdb)), '.pdb', '_core')
            pdbfnames_bfac(ipdb)   = add2fbody(basename(pdbfnames(ipdb)), '.pdb', '_core_b1_fringe_b0')
            pdbfnames_fringe(ipdb) = add2fbody(basename(pdbfnames(ipdb)), '.pdb', '_fringe')
        end do
        el = trim(adjustl(params%element))
        allocate( atms_common(npdbs) )
        do i = 1, npdbs
            call read_pdb2matrix(trim(pdbfnames(i)), pdbmat)
            allocate(atms_common(i)%coords1(3,size(pdbmat,dim=2)), source=pdbmat)
            allocate(atms_common(i)%coords2(3,size(pdbmat,dim=2)), source=pdbmat)
            deallocate(pdbmat)
        end do
        do i = 2, npdbs
            ! identify couples
            call find_couples( atms_common(i)%coords2, atms_common(i-1)%coords2, el,&
                &atms_common(i)%common1, atms_common(i-1)%common1, frac_diam=params%frac_diam)
            atms_common(i)%ncommon = size(atms_common(i)%common1, dim=2)
            deallocate(atms_common(i)%coords2)
            allocate(atms_common(i)%coords2(3,atms_common(i)%ncommon), source=(atms_common(i)%common1 + atms_common(i-1)%common1)/2.)
        end do
        ! write PDB file
        call write_matrix2pdb(el, atms_common(npdbs)%coords2, 'core.pdb')
        ! identify common atoms
        do i = 1, npdbs
            ! identify couples
            call find_couples( atms_common(npdbs)%coords2, atms_common(i)%coords1, el,&
                &atms_common(i)%common2, atms_common(i)%common1, frac_diam=params%frac_diam)
            atms_common(i)%ncommon = size(atms_common(i)%common1, dim=2)
            ! write PDB file
            call write_matrix2pdb(el, atms_common(i)%common1, pdbfnames_core(i))
        end do
        ! calculate displacements and distances
        do i = 1, npdbs
            if( i == 1 )then
                allocate(atms_common(i)%dists(atms_common(i)%ncommon), source=0.)
            else
                atms_common(i)%dists = dists_btw_common(atms_common(i)%common1, atms_common(i-1)%common1)
            endif
            call calc_stats(atms_common(i)%dists, dist_stats)
            write(logfhandle,'(A)') '>>> DISTANCE STATS, COMMON ATOMS '//basename(pdbfnames(i))
            write(logfhandle,'(A,F8.4)') 'Average: ', dist_stats%avg
            write(logfhandle,'(A,F8.4)') 'Median : ', dist_stats%med
            write(logfhandle,'(A,F8.4)') 'Sigma  : ', dist_stats%sdev
            write(logfhandle,'(A,F8.4)') 'Max    : ', dist_stats%maxv
            write(logfhandle,'(A,F8.4)') 'Min    : ', dist_stats%minv
        end do
        ! identify different atoms
        do i = 1, npdbs
            natoms = size(atms_common(i)%coords1,dim=2)
            allocate(mask_core(natoms), source=.true.)
            allocate(betas(natoms),     source=0.)
            call find_atoms_subset(atms_common(i)%common1, atms_common(i)%coords1, mask_core)
            where( mask_core )
                betas = 1.
            elsewhere
                betas = 0.
            endwhere
            ! write PDB file
            call write_matrix2pdb(el, atms_common(i)%coords1, pdbfnames_bfac(i), betas)
            deallocate(mask_core, betas)
            call remove_atoms(atms_common(i)%common1, atms_common(i)%coords1, atms_common(i)%different1)
            call write_matrix2pdb(el, atms_common(i)%different1, pdbfnames_fringe(i))
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_ATOMS_ANALYSIS NORMAL STOP ****')
    end subroutine exec_tseries_core_atoms_analysis

end module simple_commander_atoms
