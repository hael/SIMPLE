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
public :: tseries_atoms_analysis_commander
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

type, extends(commander_base) :: tseries_atoms_analysis_commander
  contains
    procedure :: execute      => exec_tseries_atoms_analysis
end type tseries_atoms_analysis_commander

type :: common_atoms
    integer           :: ind1, ind2, ncommon
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
            call nano%identify_atomic_pos(a, l_fit_lattice=.false., l_discard=trim(params%discard_atoms).eq.'yes')
            ! kill
            call nano%kill
        else
            call nano%new(params%vols(1))
            ! execute
            call nano%identify_atomic_pos(a, l_fit_lattice=.true., l_discard=trim(params%discard_atoms).eq.'yes')
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
                                  &atms_common(i,j)%common1, atms_common(i,j)%common2, frac_diam=params%frac_diam)
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

end module simple_commander_atoms
