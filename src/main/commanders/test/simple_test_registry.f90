!@descr: the registry of high-level test entries: command line, threads, workers, CPU reservation, deadline, kind, identity outputs
! One row per high-level CTest entry. test=nightly and test=profile read these rows, and
! scripts/check_test_registry.py holds production/CMakeLists.txt to them (name, command line,
! TIMEOUT, OMP_NUM_THREADS). ncpu is the peak CPU reservation, master and workers together: the
! largest CPU use test=profile sampled, less 0.1 for sampling jitter, rounded up (the profile is in
! the Phase 1 record of the test-suite plan); outputs are the patterns test=identity compares; mem_gb
! and disk_gb the peak resident memory and what the run directory grows to, which the dispatcher
! packs against what the machine has free; scales whether the master's threads follow a profile's
! thread count.
module simple_test_registry
implicit none

public :: test_entry, get_test_registry, find_test_entry, ENTRY_GATE, ENTRY_PROFILE
private

integer, parameter :: ENTRY_GATE    = 1 !< a correctness gate, run by test=nightly
integer, parameter :: ENTRY_PROFILE = 2 !< a robustness profile only, run by test=profile

type :: test_entry
    character(len=32)   :: name     = '' !< the CTest name
    character(len=32)   :: prg      = '' !< simple_test_exec test=<prg>
    character(len=64)   :: args     = '' !< the other keys of its command line
    integer             :: nthr     = 1  !< OMP_NUM_THREADS of the entry
    integer             :: nworkers = 0  !< worker processes it starts at its busiest
    integer             :: wthr     = 0  !< OpenMP threads of each worker
    integer             :: ncpu     = 1  !< peak CPU reservation
    integer             :: deadline = 0  !< seconds, the CTest TIMEOUT
    integer             :: kind     = ENTRY_GATE
    character(len=1024) :: outputs  = '' !< blank-separated patterns relative to the run directory
    real                :: mem_gb   = 0. !< peak resident memory (GiB), master and workers
    real                :: disk_gb  = 0. !< what its run directory grows to (GiB)
    logical             :: scales   = .true. !< its master threads follow nthr (a profile's thread count)
end type test_entry

! the declared outputs of the one-thread workflows, below their workflow directory: maps, stacks, picked
! coordinates, FSC and sigma2 files, solve3D manifests, projects and orientation tables (the CTF text
! files hold absolute paths and the BENCH files timings, so they are not outputs)
character(len=*), parameter :: SIMWF = 'test_simulated_workflow_'
character(len=48), parameter :: SIMWF_OUT(7) = [character(len=48) :: &
    &'/simulated_workflow/*/*.mrc', '/simulated_workflow/*/*.mrcs', '/simulated_workflow/*/*.box', &
    &'/simulated_workflow/*/*.bin', '/simulated_workflow/*/solve3D_manifest.txt', &
    &'/simulated_workflow/*/*.simple', '/simulated_workflow/*.simple']
! the entries that work in a directory named by their process id (test_<name>_<pid>) declare their
! outputs relative to it, and test=identity is given that directory as dir1 and dir2; the star files
! and motion models of the stream hold absolute movie paths, so they are not outputs; state_fixture
! and species_discovery keep no files but their gate table
character(len=48), parameter :: MINI_OUT(13) = [character(len=48) :: '*.mrc', '*.mrcs', '*.pdb', '*.txt', &
    &'1_simulate_micrographs/*.mrc', '1_simulate_micrographs/*.txt', '2_mini_stream/*.mrc', '2_mini_stream/*.mrcs', &
    &'2_mini_stream/*.box', '2_mini_stream/*.bin', '2_mini_stream/*.star', '2_mini_stream/*.txt', '2_mini_stream/*.simple']
character(len=48), parameter :: SIMPTCLS_OUT(6) = [character(len=48) :: '*.mrc', '*.mrcs', '*.pdb', '*.txt', &
    &'*/*.mrc', '*/*.txt']
character(len=48), parameter :: ATOMS_OUT(3) = [character(len=48) :: '*.mrc', '*.pdb', '*.csv']
character(len=48), parameter :: PREPROC_OUT(8) = [character(len=48) :: '*.mrcs', '*.simple', '1_simulate_movies/*.mrc', &
    &'simulation_truth/*.mrc', 'simulation_truth/*.txt', 'stream_preproc/*.simple', 'stream_preproc/motion_correct/*.mrc', &
    &'stream_preproc/spprojs_completed/*.simple']
character(len=*), parameter :: SINGLEWF = 'test_single_workflow_'
character(len=48), parameter :: SINGLEWF_OUT(9) = [character(len=48) :: &
    &'/*/*.mrc', '/final_results/*.mrc', '/*.mrc', '/*/*.bin', '/*/*.pdb', '/*/*_oris.txt', '/ptcls_oris.txt', &
    &'/*/*.simple', '/*.simple']

contains

    !> every high-level entry, in the order of production/CMakeLists.txt
    function get_test_registry() result( entries )
        type(test_entry), allocatable :: entries(:)
        entries = [ &
        ! mem_gb: the profiles' peak resident memory plus 20 %, in steps of 0.5 GiB; disk_gb: the run
        ! directory sampled during the test=nightly of 2026-10-18 plus 25 %, in steps of 0.5 GiB
        ! every commander runs with nthr=1; the 2D search and its jobs overlap
        &test_entry('mini_stream_6vxx', 'mini_stream', 'suite=6vxx', 8, 0, 0, 2, 7200, ENTRY_GATE, outputs_of('', MINI_OUT), &
            &mem_gb=1.5, disk_gb=0.5, scales=.false.), &
        &test_entry('mini_stream_1jxy', 'mini_stream', 'suite=1jxy', 8, 0, 0, 3, 7200, ENTRY_GATE, outputs_of('', MINI_OUT), &
            &mem_gb=1.0, disk_gb=0.5, scales=.false.), &
        ! at master's sizes (160 particles) its run directory reached 4.1 GiB in the profile of 2026-10-19
        &test_entry('simulated_workflow_6vxx', 'simulated_workflow', 'suite=6vxx', 1, 0, 0, 2, 7200, ENTRY_GATE, &
            &outputs_of(SIMWF//'6vxx', SIMWF_OUT), mem_gb=2.5, disk_gb=5.5, scales=.false.), &
        ! its distributed solve3D runs two one-thread workers while the master waits
        &test_entry('simulated_workflow_1jxy', 'simulated_workflow', 'suite=1jxy', 1, 2, 1, 3, 7200, ENTRY_GATE, &
            &outputs_of(SIMWF//'1jxy', SIMWF_OUT), mem_gb=5.0, disk_gb=5.5, scales=.false.), &
        &test_entry('single_workflow_fcc', 'single_workflow', 'suite=fcc smpd=0.358', 1, 0, 0, 2, 7200, ENTRY_GATE, &
            &outputs_of(SINGLEWF//'fcc', SINGLEWF_OUT), mem_gb=2.0, disk_gb=1.0, scales=.false.), &
        &test_entry('single_workflow_wurtzite', 'single_workflow', 'suite=wurtzite smpd=0.358', 1, 0, 0, 2, 7200, &
            &ENTRY_GATE, outputs_of(SINGLEWF//'wurtzite', SINGLEWF_OUT), mem_gb=3.0, disk_gb=3.0, scales=.false.), &
        &test_entry('simulate_particles', 'simulate_particles', 'nthr=8', 8, 0, 0, 8, 3600, ENTRY_GATE, &
            &outputs_of('', SIMPTCLS_OUT), mem_gb=2.0, disk_gb=0.5), &
        &test_entry('solve3D_addon', 'solve3D_addon', 'nthr=8', 8, 0, 0, 11, 7200, ENTRY_GATE, '', mem_gb=1.5, disk_gb=2.0), &
        ! its distributed step runs two workers of four threads
        &test_entry('cont_refine3D_1jxy', 'cont_refine3D_1jxy', 'nthr=8', 8, 2, 4, 10, 7200, ENTRY_GATE, '', mem_gb=11.0, disk_gb=8.0), &
        ! its reconstruction overlaps the eight search threads (peak 10.15 cores)
        &test_entry('multistate_assignment', 'multistate_assignment', 'nthr=8', 8, 0, 0, 11, 3600, ENTRY_GATE, '', mem_gb=1.5, disk_gb=1.0), &
        ! as above: reconstruction overlaps the eight search threads (peak 10.8 cores); its four routes take
        ! 32 minutes alone, an hour beside other entries (2026-10-19)
        &test_entry('multistate_workflows', 'multistate_workflows', 'nthr=8', 8, 0, 0, 11, 7200, ENTRY_GATE, '', mem_gb=4.0, disk_gb=1.0), &
        &test_entry('state_fixture', 'state_fixture', 'nthr=8', 8, 0, 0, 8, 1800, ENTRY_GATE, 'metrics.tsv', mem_gb=2.0, &
            &disk_gb=1.0), &
        &test_entry('single_atoms_stats', 'single_atoms_stats', 'smpd=0.358 element=Pt', 8, 0, 0, 3, 3600, &
            &ENTRY_GATE, outputs_of('', ATOMS_OUT), mem_gb=0.5, disk_gb=1.0), &
        &test_entry('single_atoms_stats_cdsew', 'single_atoms_stats', 'smpd=0.358 element=CdSeW', 8, 0, 0, 4, 3600, &
            &ENTRY_GATE, outputs_of('', ATOMS_OUT), mem_gb=0.5, disk_gb=1.0), &
        &test_entry('single_atoms_stats_cdsez', 'single_atoms_stats', 'smpd=0.358 element=CdSeZ', 8, 0, 0, 3, 3600, &
            &ENTRY_GATE, outputs_of('', ATOMS_OUT), mem_gb=0.5, disk_gb=1.0), &
        &test_entry('single_atoms_stats_cdser', 'single_atoms_stats', 'smpd=0.358 element=CdSeR', 8, 0, 0, 3, 3600, &
            &ENTRY_GATE, outputs_of('', ATOMS_OUT), mem_gb=0.5, disk_gb=1.0), &
        &test_entry('species_discovery', 'species_discovery', 'nthr=8', 8, 0, 0, 8, 3600, ENTRY_GATE, 'metrics.tsv', &
            &mem_gb=1.0, disk_gb=0.5), &
        ! the preprocessing stage submits one-thread jobs to the local queue, one at a time
        &test_entry('stream_preproc', 'stream_preproc', '', 8, 1, 1, 2, 3600, ENTRY_GATE, outputs_of('', PREPROC_OUT), &
            &mem_gb=0.5, disk_gb=0.5, scales=.false.), &
        ! two workers of four threads in its distributed FLEX run; the state refinements take eight
        &test_entry('flex_pca_blobs', 'flex_pca_blobs', 'nthr=4', 8, 2, 4, 9, 7200, ENTRY_GATE, '', mem_gb=3.5, disk_gb=2.5), &
        &test_entry('solve3D_pyramid_blobs', 'solve3D_pyramid_blobs', 'nthr=8', 8, 0, 0, 9, 7200, ENTRY_GATE, '', mem_gb=2.0, disk_gb=0.5) ]
    end function get_test_registry

    !> the registry index of the entry called name, 0 if there is none
    integer function find_test_entry( entries, name )
        type(test_entry), intent(in) :: entries(:)
        character(len=*), intent(in) :: name
        integer :: i
        find_test_entry = 0
        do i = 1, size(entries)
            if( trim(entries(i)%name) == trim(name) )then
                find_test_entry = i
                return
            endif
        enddo
    end function find_test_entry

    !> the blank-separated patterns of a workflow directory root
    function outputs_of( root, patterns ) result( outputs )
        character(len=*), intent(in) :: root, patterns(:)
        character(len=1024) :: outputs
        character(len=:), allocatable :: joined
        integer :: i
        joined = ''
        do i = 1, size(patterns)
            joined = joined//' '//root//trim(patterns(i))
        enddo
        if( len(joined) > len(outputs) ) error stop 'simple_test_registry: declared outputs too long'
        outputs = adjustl(joined)
    end function outputs_of

end module simple_test_registry
