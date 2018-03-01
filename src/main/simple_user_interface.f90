module simple_user_interface
use simple_defs
implicit none

public :: cluster2D, new_cluster2D
private

logical, parameter :: DEBUG = .false.

type simple_input_param
    character(len=:), allocatable :: key
    character(len=:), allocatable :: keytype ! (binary|multi|num|str|file)
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: descr_placeholder
    logical :: required = .true.
end type simple_input_param

type :: simple_program
    character(len=:), allocatable :: name
    character(len=:), allocatable :: descr_short
    character(len=:), allocatable :: descr_long
    character(len=:), allocatable :: executable
    ! image input/output
    type(simple_input_param), allocatable :: img_ios(:)
    ! parameter input/output
    type(simple_input_param), allocatable :: parm_ios(:)
    ! alternative inputs
    type(simple_input_param), allocatable :: alt_ios(:)
    ! search controls
    type(simple_input_param), allocatable :: srch_ctrls(:)
    ! filter controls
    type(simple_input_param), allocatable :: filt_ctrls(:)
    ! mask controls
    type(simple_input_param), allocatable :: mask_ctrls(:)
    ! computer controls
    type(simple_input_param), allocatable :: comp_ctrls(:)
    ! existence flag
    logical :: exists = .false.
contains
    procedure, private :: new
    procedure, private :: set_input
    procedure          :: print_ui
    procedure, private :: kill
    procedure          :: write2json

end type simple_program

type(simple_program), protected :: cluster2D

contains

    subroutine new( self, name, descr_short, descr_long, executable, n_img_ios, n_parm_ios,&
        &n_alt_ios, n_srch_ctrls, n_filt_ctrls, n_mask_ctrls, n_comp_ctrls )
        class(simple_program), intent(inout) :: self
        character(len=*),      intent(in)    :: name, descr_short, descr_long, executable
        integer,               intent(in)    :: n_img_ios, n_parm_ios, n_alt_ios, n_srch_ctrls
        integer,               intent(in)    :: n_filt_ctrls, n_mask_ctrls, n_comp_ctrls
        call self%kill
        allocate(self%name,        source=trim(name)       )
        allocate(self%descr_short, source=trim(descr_short))
        allocate(self%descr_long,  source=trim(descr_long) )
        allocate(self%executable,  source=trim(executable) )
        if( n_img_ios    > 0 ) allocate(self%img_ios(n_img_ios)      )
        if( n_parm_ios   > 0 ) allocate(self%parm_ios(n_parm_ios)    )
        if( n_alt_ios    > 0 ) allocate(self%alt_ios(n_alt_ios)      )
        if( n_srch_ctrls > 0 ) allocate(self%srch_ctrls(n_srch_ctrls))
        if( n_filt_ctrls > 0 ) allocate(self%filt_ctrls(n_filt_ctrls))
        if( n_mask_ctrls > 0 ) allocate(self%mask_ctrls(n_mask_ctrls))
        if( n_comp_ctrls > 0 ) allocate(self%comp_ctrls(n_comp_ctrls))
        self%exists = .true.
    end subroutine new

    subroutine set_input( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder, required )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                       intent(in)    :: required
        select case(trim(which))
            case('img_ios')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field img_ios'
                call set(self%img_ios, i)
            case('parm_ios')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field parm_ios'
                call set(self%parm_ios, i)
            case('alt_ios')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field alt_ios'
                call set(self%alt_ios, i)
            case('srch_ctrls')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field srch_ctrls'
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field filt_ctrls'
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field mask_ctrls'
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input, updating field comp_ctrls'
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                if( DEBUG ) print *, 'DEBUG user_interface :: set_input :: set updating entry: ', i
                allocate(arr(i)%key,               source=trim(key))
                allocate(arr(i)%keytype,           source=trim(keytype))
                allocate(arr(i)%descr_short,       source=trim(descr_short))
                allocate(arr(i)%descr_long,        source=trim(descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(descr_placeholder))
                arr(i)%required = required
                if( DEBUG ) print *, arr(i)%key, arr(i)%keytype, arr(i)%descr_short, arr(i)%descr_long, arr(i)%descr_placeholder, arr(i)%required
            end subroutine set

    end subroutine set_input

    subroutine print_ui( self )
        use simple_chash,   only: chash
        use simple_strings, only: int2str
        class(simple_program), intent(in) :: self
        type(chash) :: ch
        integer     :: i
        write(*,'(a)') ''
        write(*,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name)
        call ch%push('descr_short', self%descr_short)
        call ch%push('descr_long',  self%descr_long)
        call ch%push('executable',  self%executable)
        call ch%print_key_val_pairs
        call ch%kill
        write(*,'(a)') ''
        write(*,'(a)') '>>> IMAGE INPUT/OUTPUT'
        call print_param_hash(self%img_ios)
        write(*,'(a)') ''
        write(*,'(a)') '>>> PARAMETER INPUT/OUTPUT'
        call print_param_hash(self%parm_ios)
        write(*,'(a)') ''
        write(*,'(a)') '>>> ALTERNATIVE INPUTS'
        call print_param_hash(self%alt_ios)
        write(*,'(a)') ''
        write(*,'(a)') '>>> SEARCH CONTROLS'
        call print_param_hash(self%srch_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') '>>> FILTER CONTROLS'
        call print_param_hash(self%filt_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') '>>> MASK CONTROLS'
        call print_param_hash(self%mask_ctrls)
        write(*,'(a)') ''
        write(*,'(a)') '>>> COMPUTER CONTROLS'
        call print_param_hash(self%comp_ctrls)

        contains

            subroutine print_param_hash( arr )
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        write(*,'(a,1x,i3)') '>>> PARAMETER #', i
                        if( DEBUG ) print *, arr(i)%key, arr(i)%keytype, arr(i)%descr_short, arr(i)%descr_long, arr(i)%descr_placeholder, arr(i)%required
                        call ch%new(6)
                        call ch%push('key',               arr(i)%key)
                        call ch%push('keytype',           arr(i)%keytype)
                        call ch%push('descr_short',       arr(i)%descr_short)
                        call ch%push('descr_long',        arr(i)%descr_long)
                        call ch%push('descr_placeholder', arr(i)%descr_placeholder)
                        if( arr(i)%required )then
                            call ch%push('required', 'T')
                        else
                            call ch%push('required', 'F')
                        endif
                        call ch%print_key_val_pairs
                        call ch%kill
                    end do
                endif
            end subroutine print_param_hash

    end subroutine print_ui

    subroutine kill( self )
        class(simple_program), intent(inout) :: self
        integer :: i, sz
        if( self%exists )then
            deallocate(self%name, self%descr_short, self%descr_long, self%executable)
            call dealloc_field(self%img_ios)
            call dealloc_field(self%parm_ios)
            call dealloc_field(self%alt_ios)
            call dealloc_field(self%srch_ctrls)
            call dealloc_field(self%filt_ctrls)
            call dealloc_field(self%mask_ctrls)
            call dealloc_field(self%comp_ctrls)
            self%exists = .false.
        endif

        contains

            subroutine dealloc_field( arr )
                type(simple_input_param), allocatable, intent(inout) :: arr(:)
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        if( allocated(arr(i)%key)               ) deallocate(arr(i)%key              )
                        if( allocated(arr(i)%keytype)           ) deallocate(arr(i)%keytype          )
                        if( allocated(arr(i)%descr_short)       ) deallocate(arr(i)%descr_short      )
                        if( allocated(arr(i)%descr_long)        ) deallocate(arr(i)%descr_long       )
                        if( allocated(arr(i)%descr_placeholder) ) deallocate(arr(i)%descr_placeholder)
                    end do
                    deallocate(arr)
                endif
            end subroutine dealloc_field

    end subroutine kill

    subroutine write2json( self )
        ! use, intrinsic :: iso_fortran_env, only: wp => real64
        ! use json_module
        ! use simple_strings, only: int2str
        class(simple_program), intent(in) :: self
        ! logical, allocatable     :: required(:)
        ! type(json_core)          :: json
        ! type(json_value),pointer :: pjson, json_header, json_required, json_optional, entry
        ! integer                  :: sz, n_required, n_optional, i
        ! ! JSON init
        ! call json%initialize()
        ! call json%create_object(pjson,'')
        ! call json%create_object(json_header,'header')
        ! call json%create_object(json_required,'required')
        ! call json%create_object(json_optional,'optional')
        ! call json%create_object(json_optional,'optional')
        ! call json%add(pjson, json_header)
        ! call json%add(pjson, json_required)
        ! call json%add(pjson, json_optional)
        ! ! init
        ! sz = size(self%inputs)
        ! allocate(required(sz))
        ! do i=1,sz
        !     required(i) = self%inputs(i)%required
        ! end do
        ! n_required = count(required)
        ! n_optional = sz - n_required
        ! ! header
        ! call json%add(json_header, 'name',        self%name)
        ! call json%add(json_header, 'descr_short', self%descr_short)
        ! call json%add(json_header, 'descr_short', self%descr_long)
        ! call json%add(json_header, 'executable',  self%executable)
        ! ! required
        ! do i=1,sz
        !     call json%create_object(entry, self%inputs(i)%key)
        !     call json%add(entry,'keytype',           self%inputs(i)%keytype)
        !     call json%add(entry,'descr_short',       self%inputs(i)%descr_short)
        !     call json%add(entry,'descr_long',        self%inputs(i)%descr_long)
        !     call json%add(entry,'descr_placeholder', self%inputs(i)%descr_placeholder)
        !     if( self%inputs(i)%required )then
        !         call json%add(entry,'required',      '.true.')
        !         call json%add(json_required, entry)
        !     else
        !         call json%add(entry, 'required',     '.false.')
        !         call json%add(json_optional, entry)
        !     endif
        ! end do
        ! ! write & clean
        ! call json%print(pjson, trim(adjustl(self%name))//'.json')
        ! call json%destroy(pjson)
    end subroutine write2json

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',& ! descr_long
        &'simple_distr_exec',& ! executable
        &3, 5, 0, 7, 7, 2, 2)  ! # entries in each group

        ! INPUT PARAMETER SPECIFICATIONS

        ! image input/output
        call cluster2D%set_input('img_ios', 1, 'stk', 'file', 'Particle image stack', 'Particle image stack',&
        &'xxx.mrc file with particles', .false.)

        call cluster2D%set_input('img_ios', 2, 'stktab', 'file', 'List of per-micrograph particle stacks',&
        &'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false.)

        call cluster2D%set_input('img_ios', 3, 'refs', 'file', 'Initial references',&
        &'Initial 2D references used to bootstrap the 2D search', 'xxx.mrc file containing references', .false.)

        ! parameter input/output
        call cluster2D%set_input('parm_ios', 1, 'ctf', 'multi', 'CTF correction', 'Contrast Transfer Function correction; &
        &flip indicates that images have been phase-flipped prior(yes|no|flip){no}', '(yes|no|flip){no}', .true.)

        call cluster2D%set_input('parm_ios', 2, 'smpd', 'num', 'Sampling distance', 'Distance between neighbouring pixels in Angstroms',&
        &'pixel size in Angstroms', .true.)

        call cluster2D%set_input('parm_ios', 3, 'phaseplate', 'binary', 'Phase-plate images', 'Images obtained with Volta &
        &phase-plate(yes|no){no}', '(yes|no){no}', .false.)

        call cluster2D%set_input('parm_ios', 4, 'deftab', 'file', 'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE &
        &project (*.simple) format with dfx, dfy and angast values', '.simple|.txt parameter file', .false.)

        call cluster2D%set_input('parm_ios', 5, 'oritab', 'file', 'Orientation and CTF parameter file', 'Orientation and CTF parameter &
        &file in plain text (.txt) or SIMPLE project (*.simple) format', '.simple|.txt parameter file', .false.)

        ! alternative inputs
        !<empty>

        ! search controls
        call cluster2D%set_input('srch_ctrls', 1, 'ncls', 'num', 'Number of 2D clusters', 'Number of groups to sort the particles &
        &into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .true.)

        call cluster2D%set_input('srch_ctrls', 2, 'startit', 'num', 'First iteration', 'Index of first iteration when starting from a &
        &previous solution', 'start iterations from here', .false.)

        call cluster2D%set_input('srch_ctrls', 3, 'trs', 'num', 'Maximum translational shift', 'Maximum half-width for bund-constrained &
        &search of rotational origin shifts', 'max shift per iteration in pixels', .false.)

        call cluster2D%set_input('srch_ctrls', 4, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false.)

        call cluster2D%set_input('srch_ctrls', 5, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false.)

        call cluster2D%set_input('srch_ctrls', 6, 'dyncls', 'binary', 'Dynamic reallocation of clusters', 'Dynamic reallocation of clusters &
        &that fall below a minimum population by randomization(yes|no){yes}', '(yes|no){yes}', .false.)

        call cluster2D%set_input('srch_ctrls', 7, 'maxits', 'num', 'Max iterations', 'Maximum number of iterations',&
        &'Max # iterations{50}', .false.)

        ! filter controls
        call cluster2D%set_input('filt_ctrls', 1, 'hp', 'num', 'High-pass limit', 'High-pass resolution limit for omitting inelastic &
        &scattering contributions at low resolution', 'high-pass limit in Angstroms', .false.)

        call cluster2D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false.)

        call cluster2D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit to apply to diagnose possible &
        &issues with the dynamic update scheme used by default', 'low-pass limit in Angstroms', .false.)

        call cluster2D%set_input('filt_ctrls', 4, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &few iterations of search, before the automatic scheme kicks in. Also controls the degree of downsampling in the first &
        phase', 'initial low-pass limit in Angstroms', .false.)

        call cluster2D%set_input('filt_ctrls', 5, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false.)

        call cluster2D%set_input('filt_ctrls', 6, 'match_filt', 'binary', 'Matched filter', 'Filter to maximize the signal-to-noise &
        &ratio (SNR) in the presence of additive stochastic noise. Sometimes causes over-fitting and needs to be turned off(yes|no){yes}',&
        '(yes|no){yes}', .false.)

        call cluster2D%set_input('filt_ctrls', 7, 'weights2D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false.)

        ! mask controls
        call cluster2D%set_input('mask_ctrls', 1, 'msk', 'num', 'Mask radius', 'Mask radius in pixels for application of a soft-edged &
        &circular mask to remove background noise', 'mask radius in pixels', .true.)

        call cluster2D%set_input('mask_ctrls', 2, 'inner', 'num', 'Inner mask radius', 'Inner mask radius for omitting unordered cores of &
        &particles with high radial symmetry, typically icosahedral viruses', 'inner mask radius in pixles', .false.)

        ! computer controls
        call cluster2D%set_input('comp_ctrls', 1, 'nparts', 'num', 'Number of parts', 'Number of partitions for distrbuted memory execution.&
        & One part typically corresponds to one CPU socket in the distributed system. On a single-socket machine there may&
        & be speed benfits to dividing the jobs into a few (2-4) partitions, depending on memory capacity',&
        &'divide job into # parts', .true.)

        call cluster2D%set_input('comp_ctrls', 2, 'nthr', 'num', 'Number of threads per part', 'Number of shared-memory OpenMP threads &
        &with close affinity per partition. Typically the same as the number of logical threads in a socket.', '# shared-&
        &memory CPU threads', .false.)

    end subroutine new_cluster2D

end module simple_user_interface
