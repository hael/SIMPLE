module simple_user_interface
use simple_defs
implicit none

public :: cluster2D, new_cluster2D
private

type simple_input_param
    character(len=:), allocatable :: key
    character(len=:), allocatable :: keytype ! (binary|multi|num|str|file)
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: descr_placeholder
    logical :: required = .true.
    integer :: tie      = 0
    integer :: group    = 0
end type simple_input_param

type :: simple_program
    character(len=:),         allocatable :: name
    character(len=:),         allocatable :: descr_short
    character(len=:),         allocatable :: descr_long
    character(len=:),         allocatable :: executable
    type(simple_input_param), allocatable :: inputs(:)
    integer :: n_ties   = 0
    integer :: n_groups = 0
contains
    procedure, private :: new
    procedure, private :: set_input
    procedure          :: print
end type simple_program

type(simple_program), protected :: cluster2D

contains

    subroutine new( self, name, descr_short, descr_long, executable, n_input_params )
        class(simple_program), intent(inout) :: self
        character(len=*),      intent(in)    :: name, descr_short, descr_long, executable
        integer,               intent(in)    :: n_input_params
        if( allocated(self%name) )        deallocate(self%name)
        if( allocated(self%descr_short) ) deallocate(self%descr_short)
        if( allocated(self%descr_long) )  deallocate(self%descr_long)
        if( allocated(self%executable) )  deallocate(self%executable)
        if( allocated(self%inputs) )      deallocate(self%inputs)
        allocate(self%name,               source=trim(name))
        allocate(self%descr_short,        source=trim(descr_short))
        allocate(self%descr_long,         source=trim(descr_long))
        allocate(self%executable,         source=trim(executable))
        allocate(self%inputs(n_input_params))
    end subroutine new

    subroutine set_input( self, i, key, keytype, descr_short, descr_long, descr_placeholder, required, tie, group )
        class(simple_program), intent(inout) :: self
        integer,               intent(in)    :: i, tie, group
        character(len=*),      intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,               intent(in)    :: required
        integer :: sz
        sz = size(self%inputs)
        if( i < 1 .or. i > sz )then
            write(*,*) 'input index i: ', i, ' is out of bound, ubound is: ', sz
            stop 'ERROR, simple_user_interface :: simple_program :: set_input'
        endif
        if( allocated(self%inputs(i)%key) )               deallocate(self%inputs(i)%key)
        if( allocated(self%inputs(i)%keytype) )           deallocate(self%inputs(i)%keytype)
        if( allocated(self%inputs(i)%descr_short) )       deallocate(self%inputs(i)%descr_short)
        if( allocated(self%inputs(i)%descr_long) )        deallocate(self%inputs(i)%descr_long)
        if( allocated(self%inputs(i)%descr_placeholder) ) deallocate(self%inputs(i)%descr_placeholder)
        allocate(self%inputs(i)%key,                      source=trim(key))
        allocate(self%inputs(i)%keytype,                  source=trim(keytype))
        allocate(self%inputs(i)%descr_short,              source=trim(descr_short))
        allocate(self%inputs(i)%descr_long,               source=trim(descr_long))
        allocate(self%inputs(i)%descr_placeholder,        source=trim(descr_placeholder))
        self%inputs(i)%required = required
        self%inputs(i)%tie      = tie
        self%inputs(i)%group    = group
    end subroutine set_input

    subroutine print( self )
        use simple_chash,   only: chash
        use simple_strings, only: int2str
        class(simple_program), intent(in) :: self
        logical, allocatable :: required(:)
        type(chash)          :: ch
        integer              :: sz, n_required, n_optional, cnt, i
        write(*,'(a)') ''
        write(*,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name)
        call ch%push('descr_short', self%descr_short)
        call ch%push('descr_long',  self%descr_long)
        call ch%push('executable',  self%executable)
        call ch%print_key_val_pairs
        call ch%kill
        sz = size(self%inputs)
        allocate(required(sz))
        do i=1,sz
            required(i) = self%inputs(i)%required
        end do
        n_required = count(required)
        n_optional = sz - n_required
        write(*,'(a)') ''
        cnt = 0
        do i=1,sz
            if( self%inputs(i)%required )then
                cnt = cnt + 1
                write(*,'(a,1x,i3)') '>>> REQUIRED PARAMETER #', cnt
                call ch%new(8)
                call ch%push('key',               self%inputs(i)%key)
                call ch%push('keytype',           self%inputs(i)%keytype)
                call ch%push('descr_short',       self%inputs(i)%descr_short)
                call ch%push('descr_long',        self%inputs(i)%descr_long)
                call ch%push('descr_placeholder', self%inputs(i)%descr_placeholder)
                call ch%push('required',          '.true.')
                call ch%push('tie',               int2str(self%inputs(i)%tie))
                call ch%push('group',             int2str(self%inputs(i)%group))
                call ch%print_key_val_pairs
                call ch%kill
            endif
        end do
        write(*,'(a)') ''
        cnt = 0
        do i=1,sz
            if( .not. self%inputs(i)%required )then
                cnt = cnt + 1
                write(*,'(a,1x,i3)') '>>> OPTIONAL PARAMETER #', cnt
                call ch%new(8)
                call ch%push('key',               self%inputs(i)%key)
                call ch%push('keytype',           self%inputs(i)%keytype)
                call ch%push('descr_short',       self%inputs(i)%descr_short)
                call ch%push('descr_long',        self%inputs(i)%descr_long)
                call ch%push('descr_placeholder', self%inputs(i)%descr_placeholder)
                call ch%push('required',          '.false.')
                call ch%push('tie',               int2str(self%inputs(i)%tie))
                call ch%push('group',             int2str(self%inputs(i)%group))
                call ch%print_key_val_pairs
                call ch%kill
            endif
        end do
    end subroutine print

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',& ! descr_long
        &'simple_distr_exec', 26) ! executable, n_input_params

        ! INPUT PARAMETER SPECIFICATION
        ! required parameters
        call cluster2D%set_input(1, 'ctf', 'multi', 'CTF correction', 'Contrast Transfer Function correction; &
        &flip indicates that images have been phase-flipped prior(yes|no|flip){no}', '(yes|no|flip){no}', .true., 0, 0)

        call cluster2D%set_input(2, 'msk', 'num', 'Mask radius', 'Mask radius in pixels for application of a soft-edged &
        &circular mask to remove background noise', 'mask radius in pixels', .true., 0, 0)

        call cluster2D%set_input(3, 'ncls', 'num', 'Number of 2D clusters', 'Number of groups to sort the particles &
        &into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .true., 0, 0)

        call cluster2D%set_input(4, 'smpd', 'num', 'Sampling distance', 'Distance between neighbouring pixels in Angstroms',&
        &'pixel size in Angstroms', .true., 0, 0)

        call cluster2D%set_input(5, 'nparts', 'num', 'Number of parts', 'Number of partitions for distrbuted memory execution.&
        & One part typically corresponds to one CPU socket in the distributed system. On a single-socket machine there may&
        & be speed benfits to dividing the jobs into a few (2-4) partitions, depending on memory capacity',&
        &'divide job into # parts', .true., 0, 0)

        ! optional parameters
        call cluster2D%set_input(6, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false., 0, 0)

        call cluster2D%set_input(7, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 0, 1)

        call cluster2D%set_input(8, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 0, 1)

        call cluster2D%set_input(9, 'deftab', 'file', 'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE &
        &project (*.simple) format with dfx, dfy and angast values', '.simple|.txt parameter file', .false., 0, 0)

        call cluster2D%set_input(10, 'dyncls', 'binary', 'Dynamic reallocation of clusters', 'Dynamic reallocation of clusters &
        &that fall below a minimum population by randomization(yes|no){yes}', '(yes|no){yes}', .false., 0, 0)

        call cluster2D%set_input(11, 'hp', 'num', 'High-pass limit', 'High-pass resolution limit for omitting inelastic &
        &scattering contributions at low resolution', 'high-pass limit in Angstroms', .false., 0, 0)

        call cluster2D%set_input(12, 'inner', 'num', 'Inner mask radius', 'Inner mask radius for omitting unordered cores of &
        &particles with high radial symmetry, typically icosahedral viruses', 'inner mask radius in pixles', .false., 0, 0)

        call cluster2D%set_input(13, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false., 0, 0)

        call cluster2D%set_input(14, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &few iterations of search, before the automatic scheme kicks in. Also controls the degree of downsampling in the first &
        phase', 'initial low-pass limit in Angstroms', .false., 0, 0)

        call cluster2D%set_input(15, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 0, 0)

        call cluster2D%set_input(16, 'match_filt', 'binary', 'Matched filter', 'Filter to maximize the signal-to-noise &
        &ratio (SNR) in the presence of additive stochastic noise. Sometimes causes over-fitting and needs to be turned off(yes|no){yes}',&
        '(yes|no){yes}', .false., 0, 0)

        call cluster2D%set_input(17, 'maxits', 'num', 'Max iterations', 'Maximum number of iterations', 'Max # iterations{50}',&
        .false., 0, 0)

        call cluster2D%set_input(18, 'nthr', 'num', 'Number of threads per part', 'Number of shared-memory OpenMP threads &
        &with close affinity per partition. Typically the same as the number of logical threads in a socket.', '# shared-&
        &memory CPU threads', .false., 0, 0)

        call cluster2D%set_input(19, 'oritab', 'file', 'Orientation and CTF parameter file', 'Orientation and CTF parameter &
        &file in plain text (.txt) or SIMPLE project (*.simple) format', '.simple|.txt parameter file', .false., 0, 0)

        call cluster2D%set_input(20, 'phaseplate', 'binary', 'Phase-plate images', 'Images obtained with Volta &
        &phase-plate(yes|no){no}', '(yes|no){no}', .false., 0, 0)

        call cluster2D%set_input(21, 'refs', 'file', 'Initial references', 'Initial 2D references used to bootstrap the 2D search',&
        &'xxx.mrc file containing references', .false., 0, 0)

        call cluster2D%set_input(22, 'startit', 'num', 'First iteration', 'Index of first iteration when starting from a &
        &previous solution', 'start iterations from here', .false., 0, 0)

        call cluster2D%set_input(23, 'stk', 'file', 'Particle image stack', 'Particle image stack', 'xxx.mrc file with &
        &particles', .false., 1, 0)

        call cluster2D%set_input(24, 'stktab', 'file', 'List of per-micrograph particle stacks', 'List of per-micrograph &
        &particle stacks', 'stktab.txt file containing file names', .false., 1, 0)

        call cluster2D%set_input(25, 'trs', 'num', 'Maximum translational shift', 'Maximum half-width for bund-constrained &
        &search of rotational origin shifts', 'max shift per iteration in pixels', .false., 0, 0)

        call cluster2D%set_input(26, 'weights2D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false., 0, 0)

        ! SET # GROUPS
        cluster2D%n_groups = 1

        ! SET # TIES
        cluster2D%n_ties   = 1
    end subroutine new_cluster2D


end module simple_user_interface
