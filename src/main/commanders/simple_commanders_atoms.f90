module simple_commanders_atoms
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_image,          only: image
use simple_image_bin,      only: image_bin
use simple_nanoparticle,   only: nanoparticle
use simple_atoms,          only: atoms
use simple_nanoparticle_utils
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_atoms_stats
  contains
    procedure :: execute      => exec_atoms_stats
end type commander_atoms_stats

type, extends(commander_base) :: commander_atoms_register
  contains
    procedure :: execute      => exec_atoms_register
end type commander_atoms_register

type, extends(commander_base) :: commander_crys_score
  contains
    procedure :: execute      => exec_crys_score
end type commander_crys_score

type, extends(commander_base) :: commander_conv_atom_denoise
  contains
    procedure :: execute      => exec_conv_atom_denoise
end type commander_conv_atom_denoise

type, extends(commander_base) :: commander_detect_atoms
  contains
    procedure :: execute      => exec_detect_atoms
end type commander_detect_atoms

type, extends(commander_base) :: commander_map2model_fsc
  contains
    procedure :: execute      => exec_map2model_fsc
end type commander_map2model_fsc

type, extends(commander_base) :: commander_map_validation
  contains
    procedure :: execute      => exec_map_validation
end type commander_map_validation

type, extends(commander_base) :: commander_model_validation
  contains
    procedure :: execute      => exec_model_validation
end type commander_model_validation

type, extends(commander_base) :: commander_model_validation_eo
  contains
    procedure :: execute      => exec_model_validation_eo
end type commander_model_validation_eo

type, extends(commander_base) :: commander_pdb2mrc
  contains
    procedure :: execute      => exec_pdb2mrc
end type commander_pdb2mrc

type, extends(commander_base) :: commander_tseries_atoms_rmsd
  contains
    procedure :: execute      => exec_tseries_atoms_rmsd
end type commander_tseries_atoms_rmsd

type, extends(commander_base) :: commander_tseries_core_atoms_analysis
  contains
    procedure :: execute      => exec_tseries_core_atoms_analysis
end type commander_tseries_core_atoms_analysis

type :: common_atoms
    integer           :: ind1, ind2, ncommon, ndiff
    real, allocatable :: coords1(:,:), coords2(:,:)
    real, allocatable :: common1(:,:), common2(:,:)
    real, allocatable :: different1(:,:), different2(:,:)
    real, allocatable :: displacements(:,:), dists(:)
end type common_atoms

contains

    subroutine exec_atoms_stats( self, cline )
        class(commander_atoms_stats), intent(inout) :: self
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

    subroutine exec_atoms_register( self, cline )
        class(commander_atoms_register), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(parameters)   :: params
        real, allocatable  :: matrix1(:,:), matrix2(:,:), matrix_rot(:,:)
        character(len=100) :: io_message
        type(string)       :: line, ref_pdb, cur_pdb, params_file
        integer            :: i, j, file_stat, nl, fnr, io_stat, addr, funit
        real               :: rot_mat(3,3), trans_vec(3), scale
        call params%new(cline)
        if( .not. cline%defined('maxits') )params%maxits = 5
        addr = 1 + sizeof(scale) * 9    ! rot_mat is of size 3x3
        nl   = nlines(params%fname)
        if( nl == 0 ) return
        call fopen(fnr, FILE=params%fname, STATUS='OLD', action='READ', iostat=file_stat, iomsg=io_message)
        call fileiochk("exec_atoms_register: error when opening file for reading: "//params%fname%to_char()//':'//trim(io_message), file_stat)
        do i = 1, nl
            call line%readline(fnr, io_stat)
            ! identity rotation matrix
            rot_mat=0.; do j=1,3 ; rot_mat(j,j)=1. ; end do
            ! zero translation
            scale       = 1.
            trans_vec   = 0.
            cur_pdb     = line
            params_file = 'PARAMS_'//cur_pdb%to_char()
            if( i == 1 )then
                ! the first one is the reference pdb
                ref_pdb = cur_pdb
                call read_pdb2matrix( ref_pdb, matrix1 )
                call write_matrix2pdb( 'Pt', matrix1, string('DOCKED_')//ref_pdb )
            else
                print *, 'Docking nanoparticle ', cur_pdb%to_char(), ' to the reference nanoparticle ', ref_pdb%to_char()
                if( allocated(matrix_rot) )deallocate(matrix_rot)
                call read_pdb2matrix( line, matrix2 )
                allocate(matrix_rot(3,size(matrix2,2)))
                if( file_exists(params_file) )then
                    print *, 'Using cached rotation matrix and translation vector from ', params_file%to_char()
                    call fopen(funit,params_file,access='STREAM',action='READ',status='OLD', iostat=io_stat)
                    read(unit=funit,pos=1)    rot_mat
                    read(unit=funit,pos=addr) trans_vec
                    call fclose(funit)
                    do j = 1, size(matrix2,2)
                        matrix_rot(:,j) = scale * matmul(rot_mat, matrix2(:,j)) + trans_vec
                    enddo
                else
                    call atoms_register(matrix2, matrix1, matrix_rot, maxits=params%maxits,&
                        &out_mat=rot_mat, out_trans=trans_vec, out_scale=scale)
                endif
                call write_matrix2pdb( 'Pt', matrix_rot, string('DOCKED_')//line )
            endif
            if( .not. file_exists(params_file) )then
                ! write out the rotation matrix and translation vector
                call fopen(funit,params_file,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
                write(funit,pos=1)    rot_mat
                write(funit,pos=addr) trans_vec
                call fclose(funit)
            endif
            print *, '--------------------------'
        enddo
        call fclose(fnr)
        ! end gracefully
        call simple_end('**** SIMPLE_ATOMS_REGISTER NORMAL STOP ****')
    end subroutine exec_atoms_register

    subroutine exec_crys_score( self, cline )
        use simple_commanders_sim, only: commander_simulate_atoms
        class(commander_crys_score), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(string),                allocatable   :: pdbfnames(:), corenames(:)
        real,                        allocatable   :: dists(:), mat(:,:), dists_cen(:), vars(:), ref_stats(:), cur_mat(:,:), cur_stats(:)
        integer,                     allocatable   :: cnts(:)
        logical,                     allocatable   :: l_dists(:)
        type(commander_simulate_atoms) :: xsim_atoms
        type(nanoparticle)             :: nano
        type(parameters)               :: p
        type(cmdline)                  :: cline_sim_atoms
        type(string)                   :: t_name, cmd
        integer :: npdbs, ipdb, icore, ncores, Natoms, cnt, i, j, nclus
        real    :: min_dist, dist, prob, d, moldiam, a(3)
        call p%new(cline)
        call read_filetable(p%fname, pdbfnames)
        npdbs = size(pdbfnames)
        do ipdb = 1, npdbs
            ! finding max_dist for moldiam
            call read_pdb2matrix( p%pdbfile, mat )
            Natoms  = size(mat,2)
            moldiam = 0.
            do i = 1, Natoms-1
                do j = i+1, Natoms
                    dist = sqrt(sum((mat(:,i) - mat(:,j))**2))
                    if( dist > moldiam ) moldiam = dist
                enddo
            enddo
            ! simulating nanoparticle with the found moldiam
            call cline_sim_atoms%set('mkdir',   'no')
            call cline_sim_atoms%set('smpd',    p%smpd)
            call cline_sim_atoms%set('element', p%element)
            call cline_sim_atoms%set('outvol',  'tmp.mrc')
            call cline_sim_atoms%set('moldiam', moldiam)
            call cline_sim_atoms%set('box',     p%box)
            call cline_sim_atoms%set('nthr',    real(p%nthr))
            call xsim_atoms%execute(cline_sim_atoms)
            ! detecting atom positions of the simulated nanoparticle
            call nano%new(string('tmp.mrc'))
            call nano%identify_atomic_pos(a, l_fit_lattice=.true., l_atom_thres=trim(p%atom_thres).eq.'yes', l_print=.false.)
            call nano%kill
            ! generating stats for the reference lattice
            call read_pdb2matrix(string('tmp_ATMS.pdb'), mat )
            Natoms = size(mat,2)
            if( allocated(l_dists) )deallocate(l_dists)
            if( allocated(  dists) )deallocate(  dists)
            allocate(l_dists(Natoms**2), source=.false.)
            allocate(  dists(Natoms**2), source=0.)
            ! finding all dists and min_dist
            cnt      = 0
            min_dist = huge(min_dist)
            do i = 1, Natoms-1
                do j = i+1, Natoms
                    cnt          = cnt + 1
                    dists(cnt)   = sqrt(sum((mat(:,i) - mat(:,j))**2))
                    l_dists(cnt) = .true.
                    if( dists(cnt) < min_dist .and. dists(cnt) > TINY ) min_dist = dists(cnt)
                enddo
            enddo
            dists = pack(dists, mask=l_dists)
            call sort_cluster(dists, nclus, dists_cen, cnts, vars, delta=min_dist / 10.)
            call hpsort(dists_cen)
            if( allocated(ref_stats) )deallocate(ref_stats)
            if( allocated(cur_stats) )deallocate(cur_stats)
            allocate(ref_stats(nclus), cur_stats(nclus), source=0.)
            call compute_stats(mat, dists_cen(1:nclus), ref_stats)
            ! generating stats for the core
            t_name = pdbfnames(ipdb)
            print *, '--------------------------'
            print *, 'Computing crystal score of cores in folder : ', t_name%to_char()
            cmd = 'ls '// t_name%to_char() // '/*_core.pdb | xargs -n 1 basename > output_'//int2str(ipdb)//'.log'
            call execute_command_line(cmd%to_char())
            call read_filetable(string('output_'//int2str(ipdb)//'.log'), corenames)
            ncores = size(corenames)
            do icore = 1, ncores
                call read_pdb2matrix( t_name // '/' // corenames(icore), cur_mat )
                call compute_stats(cur_mat, dists_cen(1:nclus), cur_stats)
                call kstwo( ref_stats, nclus, cur_stats, nclus, d, prob )
                print *, 'Core : ', corenames(icore)%to_char(), ' with score : ', prob
            enddo
            call execute_command_line('rm largest_cc.mrc split_ccs.mrc binarized.mrc tmp* output_'//int2str(ipdb)//'.log')
            print *, '--------------------------'
        enddo
        ! end gracefully
        call simple_end('**** SIMPLE_CRYS_SCORE NORMAL STOP ****')

      contains

        ! compute the distribution of the dists around the reference distance points
        subroutine compute_stats(pos_mat, ref_dists, out_stats, msk)
            real,    allocatable, intent(in)    :: pos_mat(:,:)
            real,                 intent(in)    :: ref_dists(:)
            real,    allocatable, intent(inout) :: out_stats(:)
            logical, optional,    intent(in)    :: msk(:)
            integer :: Na, Nr, i, j, k
            real    :: eps_g, eps_l
            Na        = size(pos_mat,2)
            Nr        = size(ref_dists)
            out_stats = 0.
            do i = 1, Nr
                if( i > 1 .and. i < Nr )then
                    eps_g = (ref_dists(i+1) - ref_dists(i))/2.
                    eps_l = (ref_dists(i)   - ref_dists(i-1))/2.
                elseif( i == 1 )then
                    eps_g = (ref_dists(i+1) - ref_dists(i))/2.
                    eps_l =  huge(eps_l)
                elseif( i == Nr )then
                    eps_g =  huge(eps_g)
                    eps_l = (ref_dists(i)   - ref_dists(i-1))/2.
                endif
                do j = 1, Na-1
                    do k = j+1, Na
                        if( present(msk) )then
                            if( .not.msk(j) .or. .not.msk(k) )cycle
                        endif
                        dist = sqrt(sum((pos_mat(:,j) - pos_mat(:,k))**2))
                        if( (dist <= ref_dists(i) .and. abs(dist - ref_dists(i)) < eps_l ) .or. &
                            (dist >= ref_dists(i) .and. abs(dist - ref_dists(i)) < eps_g ) )then
                            out_stats(i) = out_stats(i) + 1.
                        endif
                    enddo
                enddo
            enddo
            out_stats = out_stats * 100. / sum(out_stats)
        end subroutine compute_stats

        subroutine sort_cluster( points_in, cur_clusN, out_points, out_cnts, out_vars, delta )
            real,                 intent(in)    :: points_in(:)
            integer,              intent(inout) :: cur_clusN
            real,    allocatable, intent(inout) :: out_points(:)
            integer, allocatable, intent(inout) :: out_cnts(:)
            real,    allocatable, intent(inout) :: out_vars(:)
            real,                 intent(in)    :: delta
            real,    allocatable :: points(:), avg_points(:), vars(:)
            integer, allocatable :: cnts(:)
            integer :: cur_ind, n_points, i, j
            real    :: cur_avg
            n_points = size(points_in)
            allocate(points(n_points), source=points_in)
            allocate(avg_points(n_points),cnts(n_points),vars(n_points))
            call hpsort(points)
            cur_clusN = 0
            cur_ind   = 1
            cur_avg   = 0.
            do i = 1, n_points - 1
                cur_avg = cur_avg+points(i)
                if( abs(points(i+1) - cur_avg/real(i-cur_ind+1)) > delta )then
                    cur_clusN = cur_clusN + 1
                    avg_points(cur_clusN) = cur_avg/real(i-cur_ind+1)
                    cnts(cur_clusN)       = i-cur_ind+1
                    vars(cur_clusN)       = 0.
                    do j = cur_ind, i
                        vars(cur_clusN) = vars(cur_clusN) + (points(i) - avg_points(cur_clusN))**2
                    enddo
                    vars(cur_clusN) = vars(cur_clusN)/real(i-cur_ind+1)
                    cur_ind = i + 1
                    cur_avg = 0.
                endif
            enddo
            ! last cluster
            cur_clusN = cur_clusN + 1
            avg_points(cur_clusN) = sum(points(cur_ind:i))/real(i-cur_ind+1)
            cnts(cur_clusN)       = i-cur_ind+1
            vars(cur_clusN)       = 0.
            do j = cur_ind, i
                vars(cur_clusN) = vars(cur_clusN) + (points(i) - avg_points(cur_clusN))**2
            enddo
            vars(cur_clusN) = vars(cur_clusN)/real(i-cur_ind+1)
            if( allocated(out_points) )deallocate(out_points)
            if( allocated(out_cnts)   )deallocate(out_cnts)
            if( allocated(out_vars)   )deallocate(out_vars)
            allocate(out_points(1:cur_clusN), source=avg_points(1:cur_clusN))
            allocate(out_cnts(1:cur_clusN),   source=cnts(1:cur_clusN))
            allocate(out_vars(1:cur_clusN),   source=vars(1:cur_clusN))
        end subroutine sort_cluster
    end subroutine exec_crys_score

    subroutine exec_conv_atom_denoise( self, cline )
        class(commander_conv_atom_denoise), intent(inout) :: self
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
        class(commander_detect_atoms), intent(inout) :: self
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

    subroutine exec_map2model_fsc( self, cline )
        use simple_commanders_resolest, only: commander_fsc
        class(commander_map2model_fsc), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(commander_pdb2mrc) :: xpdb2mrc
        type(commander_fsc)     :: xfsc
        type(cmdline)           :: cline_pdb2mrc, cline_fsc
        type(parameters)        :: params
        ! parse parameters
        call params%new(cline)
        call cline_pdb2mrc%set('prg', 'pdb2mrc')
        call cline_pdb2mrc%set('pdbfile', params%pdbfile)
        call cline_pdb2mrc%set('smpd', params%smpd)
        call cline_pdb2mrc%set('vol_dim', params%vol_dim)
        params%outvol = get_fbody(params%pdbfile,'pdb')//'_sim.mrc'
        call cline_pdb2mrc%set('outvol', params%outvol)
        call xpdb2mrc%execute(cline_pdb2mrc)
        call cline_fsc%set('prg', 'fsc')
        call cline_fsc%set('smpd', params%smpd)
        call cline_fsc%set('vol1', params%vols(1))
        call cline_fsc%set('vol2', params%outvol)
        call cline_fsc%set('mskdiam', params%mskdiam)
        call cline_fsc%set('nthr', params%nthr)
        call xfsc%execute(cline_fsc)
        ! end gracefully
        call simple_end('**** SIMPLE_MAP2MODEL_FSC NORMAL STOP ****')
    end subroutine exec_map2model_fsc
   
    subroutine exec_map_validation( self, cline )
        class(commander_map_validation), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        type(image)      :: exp_vol, sim_vol
        integer          :: ldim(3), ldim_new(3), ifoo, box, box_new
        real             :: upscaling_factor, smpd_new
        type(string)     :: sim_vol_file, pdbout, upscale_vol_file
        call params%new(cline)
        upscale_vol_file = get_fbody(params%vols(1),'mrc')//'_upscale.mrc'
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
        sim_vol_file     = get_fbody(params%pdbfile,'pdb')//'_sim_vol.mrc'
        pdbout           = get_fbody(params%pdbfile,'pdb')//'_centered.pdb'
        call molecule%pdb2mrc( params%pdbfile, sim_vol_file, smpd_new, pdb_out=pdbout, vol_dim=ldim_new)
        call sim_vol%new([box_new, box_new, box_new], smpd_new)
        call sim_vol%read(sim_vol_file)
        call exp_vol%read_and_crop(params%vols(1), params%smpd, box_new, smpd_new)
        call exp_vol%write(upscale_vol_file)
        call molecule%map_validation(exp_vol, sim_vol)!, filename='map_val_coor')
        call molecule%kill()
        call exp_vol%kill()
        call sim_vol%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_VALIDATION NORMAL STOP ****')
    end subroutine exec_map_validation

    subroutine exec_model_validation( self, cline )
        class(commander_model_validation), intent(inout) :: self
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
        use simple_commanders_volops, only: commander_sharpvol
        class(commander_model_validation_eo), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(commander_sharpvol) :: xsharpvol
        type(cmdline)    :: cline_sharpvol1, cline_sharpvol2
        type(parameters) :: params
        type(atoms)      :: molecule
        call params%new(cline)
        ! sharpening even/odd volume pairs prior to comparison
        call cline_sharpvol1%set('vol1', params%vols(2))
        call cline_sharpvol1%set('pdbfile', params%pdbfile)
        call cline_sharpvol1%set('smpd', params%pdbfile)
        call cline_sharpvol1%set('nthr', params%nthr)
        call xsharpvol%execute(cline_sharpvol1)
        call cline_sharpvol2%set('vol1', params%vols(3))
        call cline_sharpvol2%set('pdbfile', params%pdbfile)
        call cline_sharpvol2%set('smpd', params%pdbfile)
        call cline_sharpvol2%set('nthr', params%nthr)
        call xsharpvol%execute(cline_sharpvol2)
        params%vols(2) =  basename(add2fbody(params%vols(2), params%ext, PPROC_SUFFIX))
        params%vols(3) =  basename(add2fbody(params%vols(3), params%ext, PPROC_SUFFIX))
        call molecule%new(params%pdbfile)
        if( .not.cline%defined('vol2') ) params%vols(2) = get_fbody(params%vols(1),'mrc')//'_even.mrc'
        if( .not.cline%defined('vol3') ) params%vols(3) = get_fbody(params%vols(1),'mrc')//'_odd.mrc'
        call molecule%model_validation_eo(params%pdbfile, params%vols(1), params%vols(2), params%vols(3), params%smpd, params%smpd_target)
        call molecule%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_MODEL_VALIDATION_EO NORMAL STOP ****')
    end subroutine exec_model_validation_eo

    subroutine exec_pdb2mrc( self, cline )
        class(commander_pdb2mrc), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(atoms)      :: molecule
        integer          :: vol_dims(3)
        call params%new(cline)
        call molecule%new(params%pdbfile)
        if( .not.cline%defined('pdbout') )then
            params%pdbout = get_fbody(params%pdbfile,'pdb')//'_centered.pdb'
        endif
        if( .not.cline%defined('outvol') ) params%outvol = get_fbody(params%pdbfile,'pdb')//'_sim.mrc'
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
        class(commander_tseries_atoms_rmsd), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline !< command line input
        type(string),       allocatable :: pdbfnames(:)
        real,               allocatable :: pdbmat1(:,:), pdbmat2(:,:), dists(:), distmat(:,:), dists_neigh(:)
        integer,            allocatable :: inds(:), pairs(:,:), inds_neigh(:)
        type(stats_struct), allocatable :: dist_stats(:)
        type(parameters) :: params
        integer          :: npdbs, i, j, ipdb, npairs, cnt
        character(len=2) :: el
        call cline%set('mkdir', 'no')
        call params%new(cline)
        call read_filetable(params%pdbfiles, pdbfnames)
        npdbs = size(pdbfnames)
        if( params%mkdir.eq.'yes' )then
            do ipdb = 1,npdbs
                if(pdbfnames(ipdb)%to_char([1,1]).ne.'/') pdbfnames(ipdb) = string('../')//pdbfnames(ipdb)
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
                call read_pdb2matrix(pdbfnames(i), pdbmat1)
                call read_pdb2matrix(pdbfnames(j), pdbmat2)
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
        class(commander_tseries_core_atoms_analysis), intent(inout) :: self
        class(cmdline),                               intent(inout) :: cline !< command line input
        type(string),       allocatable :: pdbfnames(:), pdbfnames_core(:), pdbfnames_bfac(:), pdbfnames_fringe(:)
        type(common_atoms), allocatable :: atms_common(:)
        real,               allocatable :: betas(:), pdbmat(:,:)
        logical,            allocatable :: mask_core(:)
        type(parameters)   :: params
        integer            :: npdbs, i, ipdb, natoms
        character(len=2)   :: el
        type(stats_struct) :: dist_stats
        call params%new(cline)
        call read_filetable(params%pdbfiles, pdbfnames)
        npdbs = size(pdbfnames)
        if( params%mkdir.eq.'yes' )then
            do ipdb = 1,npdbs
                if(pdbfnames(ipdb)%to_char([1,1]).ne.'/') pdbfnames(ipdb) = string('../')//pdbfnames(ipdb)
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
            call read_pdb2matrix(pdbfnames(i), pdbmat)
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
        call write_matrix2pdb(el, atms_common(npdbs)%coords2, string('core.pdb'))
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
            write(logfhandle,'(A)') '>>> DISTANCE STATS, COMMON ATOMS '//pdbfnames(i)%to_char()
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

end module simple_commanders_atoms
