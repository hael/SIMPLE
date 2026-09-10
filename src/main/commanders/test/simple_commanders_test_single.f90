!@descr: tests for single 
module simple_commanders_test_single
use simple_commanders_api
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_atoms_stats
  contains
    procedure :: execute      => exec_test_atoms_stats
end type commander_test_atoms_stats

type, extends(commander_base) :: commander_test_detect_atoms
  contains
    procedure :: execute      => exec_test_detect_atoms
end type commander_test_detect_atoms

type, extends(commander_base) :: commander_test_detect_calpha
  contains
    procedure :: execute      => exec_test_detect_calpha
end type commander_test_detect_calpha

type, extends(commander_base) :: commander_test_detect_calpha_molecules
  contains
    procedure :: execute      => exec_test_detect_calpha_molecules
end type commander_test_detect_calpha_molecules

type, extends(commander_base) :: commander_test_simulate_nanoparticle
  contains
    procedure :: execute      => exec_test_simulate_nanoparticle
end type commander_test_simulate_nanoparticle

type, extends(commander_base) :: commander_test_single_workflow
  contains
    procedure :: execute      => exec_test_single_workflow
end type commander_test_single_workflow

integer, parameter :: BOX          = 160
integer, parameter :: MOLDIAM      = 20
integer, parameter :: NTHR         = 40

contains

subroutine exec_test_atoms_stats( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    use simple_commanders_atoms, only: commander_atoms_stats
    class(commander_test_atoms_stats), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_detat, cline_atstats
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_detect_atoms)          :: xdetat
    type(commander_atoms_stats)           :: xatstats
    write(logfhandle,'(a)') '>>> TEST_ATOMS_STATS:'
    call params%new(cline)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call cline_detat%set('prg',             'detect_atoms')
    call cline_detat%set('vol1',              'outvol.mrc')
    call cline_detat%set('smpd',               params%smpd)
    call cline_detat%set('element',         params%element)
    call cline_detat%set('nthr',                      NTHR)
    call xdetat%execute(cline_detat)
    call cline_atstats%set('prg',             'atoms_stats')
    call cline_atstats%set('vol1',             'outvol.mrc')
    call cline_atstats%set('vol2',          'outvol_CC.mrc')
    call cline_atstats%set('pdbfile',     'outvol_ATMS.pdb')
    call cline_atstats%set('smpd',              params%smpd)
    call cline_atstats%set('element',        params%element)
    call cline_atstats%set('nthr',                     NTHR)
    call xatstats%execute(cline_atstats)
    call simple_end('**** SIMPLE_TEST_ATOMS_STATS NORMAL STOP ****')
end subroutine exec_test_atoms_stats

subroutine exec_test_detect_atoms( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    class(commander_test_detect_atoms), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_detat
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_detect_atoms)          :: xdetat
    write(logfhandle,'(a)') '>>> TEST_DETECT_ATOMS:'
    call params%new(cline)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call cline_detat%set('prg',             'detect_atoms')
    call cline_detat%set('vol1',              'outvol.mrc')
    call cline_detat%set('smpd',               params%smpd)
    call cline_detat%set('element',         params%element)
    call cline_detat%set('nthr',                      NTHR)
    call xdetat%execute(cline_detat)
    call simple_end('**** SIMPLE_TEST_DETECT_ATOMS NORMAL STOP ****')
end subroutine exec_test_detect_atoms

subroutine exec_test_detect_calpha( self, cline )
    use simple_atoms,         only: atoms
    use simple_calpha_finder, only: calpha_finder
    class(commander_test_detect_calpha), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    character(len=*), parameter :: PDB_FILE = 'test_calpha_candidates.pdb'
    character(len=*), parameter :: CSV_FILE = 'test_calpha_candidates.csv'
    character(len=*), parameter :: MRC_FILE = 'test_calpha_scores.mrc'
    integer,          parameter :: TEST_BOX = 32, NRES = 3
    type(image)         :: workvol
    type(atoms)         :: candidates
    type(calpha_finder) :: finder
    real(kind=c_float), pointer :: density(:,:,:)
    real    :: centers(3,NRES), atom_sites(3,3), amplitudes(3), rotation(3,3)
    real    :: xyz(3), site(3), delta(3), distance, closest
    integer :: ldim(3), ires, iatom, ix, iy, iz

    write(logfhandle,'(A)') '>>> TEST_DETECT_CALPHA:'
    ldim            = [TEST_BOX,TEST_BOX,TEST_BOX]
    centers(:,1)    = [8.,8.,8.]
    centers(:,2)    = [16.,16.,16.]
    centers(:,3)    = [24.,24.,24.]
    atom_sites(:,1) = [0.,0.,0.]
    atom_sites(:,2) = [1.458*cos(111.2*PI/180.), 1.458*sin(111.2*PI/180.), 0.]
    atom_sites(:,3) = [1.525,0.,0.]
    amplitudes      = [1.25,1.0,1.0]
    rotation(:,1)   = [0.,1.,0.]
    rotation(:,2)   = [-0.5,0.,sqrt(0.75)]
    rotation(:,3)   = [sqrt(0.75),0.,0.5]
    call workvol%new(ldim, 1.0)
    call workvol%get_rmat_ptr(density)
    density         = 0.
    do ires         = 1, NRES
        do iatom = 1, size(atom_sites,2)
            site = centers(:,ires) + matmul(rotation, atom_sites(:,iatom))
            do iz = 1, TEST_BOX
                do iy = 1, TEST_BOX
                    do ix = 1, TEST_BOX
                        xyz = real([ix,iy,iz] - 1)
                        delta = xyz - site
                        density(ix,iy,iz) = density(ix,iy,iz) + amplitudes(iatom) * &
                            exp(-0.5 * sum(delta * delta) / 0.85**2)
                    enddo
                enddo
            enddo
        enddo
    enddo

    call finder%new(1.0, 4.0)
    call finder%search(workvol, 180.0, 10, 0.5, string(PDB_FILE), string(MRC_FILE))
    if(nlines(string(PDB_FILE)) == 0) THROW_HARD('TEST_DETECT_CALPHA FAILED: no candidates')
    call candidates%new(string(PDB_FILE))
    closest = huge(1.)
    do iatom = 1, candidates%get_n()
        do ires = 1, NRES
            distance = sqrt(sum((candidates%get_coord(iatom) - centers(:,ires))**2))
            closest  = min(closest, distance)
        enddo
    enddo
    if( closest > 1.5 ) THROW_HARD('TEST_DETECT_CALPHA FAILED: peak is displaced')

    call candidates%kill()
    call finder%kill()
    call workvol%kill()
    if(file_exists(PDB_FILE)) call del_file(PDB_FILE)
    if(file_exists(CSV_FILE)) call del_file(CSV_FILE)
    if(file_exists(MRC_FILE)) call del_file(MRC_FILE)
    write(logfhandle,'(A,F7.3,A)') '>>> TEST_DETECT_CALPHA: PASS (closest peak ', closest, ' A)'
    call simple_end('**** SIMPLE_TEST_DETECT_CALPHA NORMAL STOP ****')

end subroutine exec_test_detect_calpha

subroutine exec_test_detect_calpha_molecules( self, cline )
    use simple_atoms,         only: atoms
    use simple_calpha_finder, only: calpha_finder
    use simple_molecule_data, only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
    class(commander_test_detect_calpha_molecules), intent(inout) :: self
    class(cmdline),                                intent(inout) :: cline
    type(parameters)    :: params
    type(molecule_data) :: mol

    if( .not.cline%defined('smpd') )    call cline%set('smpd', 1.3)
    if( .not.cline%defined('angstep') ) call cline%set('angstep', 45)
    if( .not.cline%defined('thres') )   call cline%set('thres', 0.25)
    call params%new(cline)

    write(logfhandle,'(A)') '>>> C-ALPHA MOLECULE BENCHMARK:'
    write(logfhandle,'(A,F6.2,A,I0,A,F6.3)') '    smpd=', params%smpd, &
        ' A, angstep=', params%angstep, ' degrees, threshold=', params%thres
    mol = sars_cov2_spkgp_6vxx()
    call evaluate_molecule('6VXX', mol, 2916, params%smpd, params%angstep, params%thres)
    mol = betagal_1jyx()
    call evaluate_molecule('1JYX', mol, 4044, params%smpd, params%angstep, params%thres)
    call simple_end('**** SIMPLE_TEST_DETECT_CALPHA_MOLECULES NORMAL STOP ****')

contains

    subroutine evaluate_molecule( label, molecule_data_in, expected_truth, smpd, angstep, threshold )
        character(len=*),    intent(in) :: label
        type(molecule_data), intent(in) :: molecule_data_in
        integer,             intent(in) :: expected_truth, angstep
        real,                intent(in) :: smpd, threshold
        real, parameter      :: MAP_PADDING = 12.0, MATCH_RADIUS = 2.0
        type(atoms)          :: molecule, candidates
        type(calpha_finder)  :: finder
        type(image)          :: workvol
        type(string)         :: source_file, truth_file, vol_file, candidate_file, score_file
        real, allocatable    :: truth_xyz(:,:)
        logical, allocatable :: truth_matched(:)
        real    :: span(3), delta(3), best_distance_sq, recall, precision
        real    :: recall_top_n, precision_top_n
        integer :: ldim(3), iatom, itruth, ipred, ntruth, npred, nmatched, best_truth
        integer :: top_n_count, nmatched_top_n

        source_file    = trim(label)//'.pdb'
        truth_file     = trim(label)//'_calpha_truth.pdb'
        vol_file       = trim(label)//'_calpha_input.mrc'
        candidate_file = trim(label)//'_calpha_candidates.pdb'
        score_file     = trim(label)//'_calpha_scores.mrc'
        span           = maxval(molecule_data_in%xyz, dim=1) - minval(molecule_data_in%xyz, dim=1)
        ldim           = max(round2even((span + 2. * MAP_PADDING) / smpd), 16)
        call molecule%pdb2mrc(pdbfile=source_file, volfile=vol_file, smpd=smpd, &
            center_pdb=.true., pdb_out=truth_file, vol_dim=ldim, mol=molecule_data_in)

        ntruth = 0
        do iatom = 1, molecule%get_n()
            if(molecule%get_name(iatom) == ' CA ' .and. molecule%get_element(iatom) == 'C ') &
                ntruth = ntruth + 1
        enddo
        if(ntruth /= expected_truth) THROW_HARD('Unexpected built-in C-alpha count')
        allocate(truth_xyz(3,ntruth), source=0.)
        allocate(truth_matched(ntruth), source=.false.)
        itruth = 0
        do iatom = 1, molecule%get_n()
            if(molecule%get_name(iatom) /= ' CA ' .or. molecule%get_element(iatom) /= 'C ') cycle
            itruth = itruth + 1
            truth_xyz(:,itruth) = molecule%get_coord(iatom)
        enddo

        call workvol%new(ldim, smpd)
        call workvol%read(vol_file)
        call finder%new(smpd, 4.0)
        call finder%search(workvol, real(angstep), 2 * ntruth, threshold, candidate_file, score_file)

        npred = 0
        if(nlines(candidate_file) > 0)then
            call candidates%new(candidate_file)
            npred = candidates%get_n()
        endif
        nmatched       = 0
        nmatched_top_n = 0
        top_n_count = min(ntruth, npred)
        do ipred = 1, npred
            best_truth       = 0
            best_distance_sq = huge(1.)
            do itruth = 1, ntruth
                if(truth_matched(itruth)) cycle
                delta = candidates%get_coord(ipred) - truth_xyz(:,itruth)
                if(sum(delta * delta) < best_distance_sq)then
                    best_distance_sq = sum(delta * delta)
                    best_truth       = itruth
                endif
            enddo
            if(best_truth > 0 .and. best_distance_sq <= MATCH_RADIUS**2)then
                truth_matched(best_truth) = .true.
                nmatched                  = nmatched + 1
            endif
            if(ipred == top_n_count) nmatched_top_n = nmatched
        enddo
        recall_top_n    = real(nmatched_top_n) / real(ntruth)
        precision_top_n = 0.
        if(top_n_count > 0) precision_top_n = real(nmatched_top_n) / real(top_n_count)
        recall    = real(nmatched) / real(ntruth)
        precision = 0.
        if(npred > 0) precision = real(nmatched) / real(npred)

        write(logfhandle,'(A,A)') '>>> ', trim(label)
        write(logfhandle,'(A,I0,A,I0)') '    truth=', ntruth, ', candidate cap=', 2 * ntruth
        write(logfhandle,'(A,I0,A,I0,A,F7.3,A,F7.3)') '    top-N: predicted=', top_n_count, &
            ', matched=', nmatched_top_n, ', recall=', recall_top_n, ', precision=', precision_top_n
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.3,A,F7.3)') '    top-2N: predicted=', npred, &
            ', matched=', nmatched, ', missed=', ntruth - nmatched, ', recall=', recall, &
            ', precision=', precision
        write(logfhandle,'(A,3(I0,1X))') '    map dimensions=', ldim
        write(logfhandle,'(A,A)') '    candidates: ', candidate_file%to_char()
        write(logfhandle,'(A,A)') '    score volume: ', score_file%to_char()

        if(npred > 0) call candidates%kill()
        call finder%kill()
        call workvol%kill()
        call molecule%kill()
        deallocate(truth_xyz, truth_matched)
    end subroutine evaluate_molecule

end subroutine exec_test_detect_calpha_molecules

subroutine exec_test_simulate_nanoparticle( self, cline )
    use simple_commanders_sim, only: commander_simulate_nanoparticle
    class(commander_test_simulate_nanoparticle), intent(inout) :: self
    class(cmdline),                              intent(inout) :: cline
    type(cmdline)                         :: cline_sim
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_NANOPARTICLE:'
    call params%new(cline)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call simple_end('**** SIMPLE_TEST_SIMULATE_NANOPARTICLE NORMAL STOP ****')
end subroutine exec_test_simulate_nanoparticle

subroutine exec_test_single_workflow( self, cline )
    use single_commanders_nano2D,       only: commander_analysis2D_nano
    use simple_commanders_sim,          only: commander_simulate_nanoparticle
    use simple_commanders_reproject,    only: commander_reproject
    use simple_commanders_stkops,       only: commander_stackops
    use single_commanders_trajectory,   only: commander_trajectory_denoise
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_commanders_project_core, only: commander_new_project
    use single_commanders_nano3D,       only: commander_autorefine3D_nano
    class(commander_test_single_workflow), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_reproject, cline_trajectory, cline_denoise
    type(cmdline)                         :: cline_nproj, cline_imptcls, cline_an2Dnano, cline_aref3Dnano
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_reproject)             :: xreproject
    type(commander_stackops)              :: xtrajectory
    type(commander_trajectory_denoise)    :: xdenoise
    type(commander_new_project)           :: xnproj
    type(commander_import_particles)      :: ximptcls
    type(commander_analysis2D_nano)       :: xan2Dnano
    type(commander_autorefine3D_nano)     :: xaref3Dnano
    type(string)                          :: projname, projfile, project_dir, startvol
    type(string)                          :: simulated_vol, reprojections, trajectory, denoised_trajectory
    character(len=*), parameter           :: REPROJECTIONS_STK = 'reprojections.mrc'
    character(len=*), parameter           :: TRAJECTORY_STK    = 'simulated_trajectory.mrc'
    character(len=*), parameter           :: DENOISED_STK      = 'denoised_trajectory.mrc'
    character(len=*), parameter           :: TRAJECTORY_ORITAB = 'glc_trajectory_oris.txt'
    character(len=*), parameter           :: SIMULATION_DIR    = '1_simulate_nanoparticle'
    character(len=*), parameter           :: REPROJECTION_DIR  = '2_generate_reprojections'
    character(len=*), parameter           :: TRAJECTORY_DIR    = '3_generate_trajectory'
    character(len=*), parameter           :: DENOISE_DIR       = '4_trajectory_denoise'
    character(len=*), parameter           :: IMPORT_DIR        = '5_import_particles'
    character(len=*), parameter           :: ANALYSIS2D_DIR    = '6_analysis2D_nano'
    character(len=*), parameter           :: AUTOREFINE3D_DIR  = '7_autorefine3D_nano'
    integer,          parameter           :: NREPROJS = 5000, MASKDIAM = 40
    integer,          parameter           :: NFRAMES_PER_GROUP = 50
    integer                               :: chdir_status
    real,             parameter           :: TRAJECTORY_SNR    = 0.2
    write(logfhandle,'(a)') '>>> TEST_SINGLE_WORKFLOW:'
    projname = 'test_single_workflow'
    call params%new(cline)
    projfile = projname%to_char()//'.simple'
    call cline_nproj%set('prg',                       'new_project')
    call cline_nproj%set('projname',             projname%to_char())
    call xnproj%execute(cline_nproj)
    call simple_getcwd(project_dir)
    projfile = filepath(project_dir, projfile)
    simulated_vol       = filepath(filepath(project_dir, SIMULATION_DIR), 'outvol.mrc')
    reprojections       = filepath(filepath(project_dir, REPROJECTION_DIR), REPROJECTIONS_STK)
    trajectory          = filepath(filepath(project_dir, TRAJECTORY_DIR), TRAJECTORY_STK)
    denoised_trajectory = filepath(filepath(project_dir, DENOISE_DIR), DENOISED_STK)
    startvol            = filepath(filepath(project_dir, ANALYSIS2D_DIR), 'startvol.mrc')
    
    call enter_workflow_stage(SIMULATION_DIR)
    call cline_sim%set('prg',               'simulate_nanoparticle')
    call cline_sim%set('box',                                   BOX)
    call cline_sim%set('smpd',                          params%smpd)
    call cline_sim%set('moldiam',                           MOLDIAM)
    call cline_sim%set('element',                    params%element)
    call cline_sim%set('nthr',                                 NTHR)
    call xsim_nptcl%execute(cline_sim)
    call return_to_project_dir

    call enter_workflow_stage(REPROJECTION_DIR)
    call make_glc_trajectory_oris(TRAJECTORY_ORITAB, NREPROJS, NFRAMES_PER_GROUP)
    call cline_reproject%set('prg',                     'reproject')
    call cline_reproject%set('pgrp',                           'c1')
    call cline_reproject%set('vol1',        simulated_vol%to_char())
    call cline_reproject%set('smpd',                    params%smpd)
    call cline_reproject%set('oritab',            TRAJECTORY_ORITAB)
    call cline_reproject%set('mskdiam',                          20)
    call cline_reproject%set('outstk',            REPROJECTIONS_STK)
    call cline_reproject%set('nthr',                           NTHR)
    call xreproject%execute(cline_reproject)
    call return_to_project_dir

    call enter_workflow_stage(TRAJECTORY_DIR)
    call cline_trajectory%set('prg',                     'stackops')
    call cline_trajectory%set('mkdir',                         'no')
    call cline_trajectory%set('stk',         reprojections%to_char())
    call cline_trajectory%set('outstk',              TRAJECTORY_STK)
    call cline_trajectory%set('smpd',                   params%smpd)
    call cline_trajectory%set('snr',                 TRAJECTORY_SNR)
    call cline_trajectory%set('nthr',                          NTHR)
    call xtrajectory%execute(cline_trajectory)
    call return_to_project_dir

    call enter_workflow_stage(DENOISE_DIR)
    call cline_denoise%set('prg',              'trajectory_denoise')
    call cline_denoise%set('mkdir',                            'no')
    call cline_denoise%set('stk',              trajectory%to_char())
    call cline_denoise%set('outstk',                   DENOISED_STK)
    call cline_denoise%set('smpd',                      params%smpd)
    call cline_denoise%set('nthr',                             NTHR)
    call xdenoise%execute(cline_denoise)
    call return_to_project_dir

    call enter_workflow_stage(IMPORT_DIR)
    call cline_imptcls%set('prg',                'import_particles')
    call cline_imptcls%set('mkdir',                            'no')
    call cline_imptcls%set('projfile',           projfile%to_char())
    call cline_imptcls%set('stk',     denoised_trajectory%to_char())
    call cline_imptcls%set('smpd',                      params%smpd)
    call cline_imptcls%set('ctf',                              'no')
    call ximptcls%execute(cline_imptcls)
    call return_to_project_dir

    call enter_workflow_stage(ANALYSIS2D_DIR)
    call cline_an2Dnano%set('prg',                'analysis2D_nano')
    call cline_an2Dnano%set('mkdir',                           'no')
    call cline_an2Dnano%set('projfile',          projfile%to_char())
    call cline_an2Dnano%set('element',               params%element)
    call cline_an2Dnano%set('nthr',                            NTHR)
    call xan2Dnano%execute(cline_an2Dnano)
    if( .not. file_exists(startvol) ) THROW_HARD('analysis2D_nano did not generate '//startvol%to_char())
    call return_to_project_dir

    call enter_workflow_stage(AUTOREFINE3D_DIR)
    call cline_aref3Dnano%set('prg',            'autorefine3D_nano')
    call cline_aref3Dnano%set('projfile',        projfile%to_char())
    call cline_aref3Dnano%set('vol1',            startvol%to_char())
    call cline_aref3Dnano%set('smpd',                   params%smpd)
    call cline_aref3Dnano%set('element',             params%element)
    call cline_aref3Dnano%set('nthr',                          NTHR)
    call cline_aref3Dnano%set('pgrp',                          'c1')
    call cline_aref3Dnano%set('lp',                             1.5)  
    call cline_aref3Dnano%set('mskdiam',                   MASKDIAM)
    call xaref3Dnano%execute(cline_aref3Dnano)
    call return_to_project_dir
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')
contains
    subroutine enter_workflow_stage( stage )
        character(len=*), intent(in) :: stage
        call simple_mkdir(filepath(project_dir, stage))
        call simple_chdir(filepath(project_dir, stage), chdir_status)
        if( chdir_status /= 0 ) THROW_HARD('Could not enter single_workflow stage')
    end subroutine enter_workflow_stage

    subroutine return_to_project_dir
        call simple_chdir(project_dir, chdir_status)
        if( chdir_status /= 0 ) THROW_HARD('Could not return to the single_workflow project directory')
    end subroutine return_to_project_dir
end subroutine exec_test_single_workflow

subroutine make_glc_trajectory_oris( oritab, nreprojs, nframes_per_group )
    character(len=*), intent(in) :: oritab
    integer,          intent(in) :: nreprojs, nframes_per_group
    real, parameter              :: ANGULAR_SPAN = 2.0
    type(oris)                   :: group_oris, trajectory_oris
    real                         :: base_euls(3), euls(3), frame_frac
    integer                      :: iframe, igroup, iproj, ngroups
    if( nframes_per_group < 2 ) THROW_HARD('GLC orientation groups require at least two frames')
    if( mod(nreprojs, nframes_per_group) /= 0 ) THROW_HARD('GLC reprojection count must divide into equal frame groups')
    ngroups = nreprojs / nframes_per_group
    call group_oris%new(ngroups, is_ptcl=.false.)
    call group_oris%spiral()
    call trajectory_oris%new(nreprojs, is_ptcl=.false.)
    do igroup = 1, ngroups
        base_euls = group_oris%get_euler(igroup)
        do iframe = 1, nframes_per_group
            iproj      = (igroup - 1) * nframes_per_group + iframe
            frame_frac = real(iframe - 1) / real(nframes_per_group - 1)
            euls       = base_euls + ANGULAR_SPAN * (frame_frac - 0.5) * [1.0, 0.5, 1.0]
            call trajectory_oris%set_euler(iproj, euls)
        enddo
    enddo
    call trajectory_oris%set_all2single('state', 1.0)
    call trajectory_oris%write(string(oritab), [1, nreprojs])
    call trajectory_oris%kill
    call group_oris%kill
end subroutine make_glc_trajectory_oris

end module simple_commanders_test_single
