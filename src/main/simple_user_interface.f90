module simple_user_interface
use simple_defs
implicit none

public :: make_user_interface
public :: cluster2D, refine3D, initial_3Dmodel
public :: postprocess, extract, scale
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
    procedure, private :: set_input_1
    procedure, private :: set_input_2
    generic,   private :: set_input => set_input_1, set_input_2
    procedure          :: print_ui
    procedure, private :: kill
    procedure          :: write2json

end type simple_program

! declare protected program specifications here
type(simple_program), protected :: cluster2D
type(simple_program), protected :: refine3D
type(simple_program), protected :: initial_3Dmodel
type(simple_program), protected :: postprocess
type(simple_program), protected :: extract
type(simple_program), protected :: scale

! declare common params here, with name same as flag
type(simple_input_param) :: ctf
type(simple_input_param) :: deftab
type(simple_input_param) :: frac
type(simple_input_param) :: hp
type(simple_input_param) :: inner
type(simple_input_param) :: maxits
type(simple_input_param) :: msk
type(simple_input_param) :: mskfile
type(simple_input_param) :: nspace
type(simple_input_param) :: nparts
type(simple_input_param) :: nthr
type(simple_input_param) :: objfun
type(simple_input_param) :: oritab
type(simple_input_param) :: phaseplate
type(simple_input_param) :: pgrp
type(simple_input_param) :: smpd
type(simple_input_param) :: startit
type(simple_input_param) :: stk
type(simple_input_param) :: stktab
type(simple_input_param) :: trs
type(simple_input_param) :: update_frac

contains

    ! class methods

    subroutine make_user_interface
        call set_common_params
        call new_cluster2D
        call new_refine3D
        call new_initial_3Dmodel
        call new_postprocess
        call new_extract
        call new_scale
        ! ...
    end subroutine make_user_interface

    subroutine set_common_params
        call set_param(stk,         'stk',         'file',   'Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false.)
        call set_param(stktab,      'stktab',      'file',   'List of per-micrograph particle stacks', 'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false.)
        call set_param(ctf,         'ctf',         'multi',  'CTF correction', 'Contrast Transfer Function correction; flip indicates that images have been phase-flipped prior(yes|no|flip){no}',&
        &'(yes|no|flip){no}', .true.)
        call set_param(smpd,        'smpd',        'num',    'Sampling distance', 'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms', .true.)
        call set_param(phaseplate,  'phaseplate',  'binary', 'Phase-plate images', 'Images obtained with Volta phase-plate(yes|no){no}', '(yes|no){no}', .false.)
        call set_param(deftab,      'deftab',      'file',   'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format with dfx, dfy and angast values',&
        &'.simple|.txt parameter file', .false.)
        call set_param(oritab,      'oritab',      'file',   'Orientation and CTF parameter file', 'Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format',&
        &'.simple|.txt parameter file', .false.)
        call set_param(startit,     'startit',     'num',    'First iteration', 'Index of first iteration when starting from a previous solution', 'start iterations from here', .false.)
        call set_param(trs,         'trs',         'num',    'Maximum translational shift', 'Maximum half-width for bund-constrained search of rotational origin shifts',&
        &'max shift per iteration in pixels', .false.)
        call set_param(maxits,      'maxits',      'num',    'Max iterations', 'Maximum number of iterations', 'Max # iterations', .false.)
        call set_param(hp,          'hp',          'num',    'High-pass limit', 'High-pass resolution limit', 'high-pass limit in Angstroms', .false.)
        call set_param(msk,         'msk',         'num',    'Mask radius', 'Mask radius in pixels for application of a soft-edged circular mask to remove background noise', 'mask radius in pixels', .true.)
        call set_param(inner,       'inner',       'num',    'Inner mask radius', 'Inner mask radius for omitting unordered cores of particles with high radial symmetry, typically icosahedral viruses',&
        &'inner mask radius in pixles', .false.)
        call set_param(nparts,      'nparts',      'num',    'Number of parts', 'Number of partitions for distrbuted memory execution. One part typically corresponds to one CPU socket in the distributed &
        &system. On a single-socket machine there may be speed benfits to dividing the jobs into a few (2-4) partitions, depending on memory capacity', 'divide job into # parts', .true.)
        call set_param(nthr,        'nthr',        'num',    'Number of threads per part', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .false.)
        call set_param(update_frac, 'update_frac', 'num',    'Fractional update per iteration', 'Fraction of particles to update per iteration in incremental learning scheme for accelerated convergence &
        &rate(0.1-0.5){1.}', 'update this fraction per iter{1.0}', .false.)
        call set_param(frac,        'frac',        'num',    'Fraction of particles to include', 'Fraction of particles to include based on spectral score (median of FRC between reference and particle)',&
        'fraction of particles used{1.0}', .false.)
        call set_param(mskfile,     'mskfile',     'file',   'Input mask file', 'Input mask file to apply to reference volume(s) before projection', 'e.g. automask.mrc from postprocess', .false.)
        call set_param(pgrp,        'pgrp',        'str',    'Point-group symmetry', 'Point-group symmetry of particle(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .true.)
        call set_param(nspace,      'nspace',      'num',    'Number of projection directions', 'Number of projection directions &
        &used in discrete 3D orientation search', '# projections used', .false.)
        call set_param(objfun,      'objfun',      'binary', 'Objective function', 'Objective function(cc|ccres){cc}', '(cc|ccres){cc}', .false.)

        contains

            subroutine set_param( self, key, keytype, descr_short, descr_long, descr_placeholder, required )
                type(simple_input_param), intent(inout) :: self
                character(len=*),         intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
                logical,                  intent(in)    :: required
                allocate(self%key,               source=trim(key))
                allocate(self%keytype,           source=trim(keytype))
                allocate(self%descr_short,       source=trim(descr_short))
                allocate(self%descr_long,        source=trim(descr_long))
                allocate(self%descr_placeholder, source=trim(descr_placeholder))
                self%required = required
            end subroutine set_param

    end subroutine set_common_params

    ! instance methods

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

    subroutine set_input_1( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder, required )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                       intent(in)    :: required
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input_1'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                allocate(arr(i)%key,               source=trim(key))
                allocate(arr(i)%keytype,           source=trim(keytype))
                allocate(arr(i)%descr_short,       source=trim(descr_short))
                allocate(arr(i)%descr_long,        source=trim(descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(descr_placeholder))
                arr(i)%required = required
            end subroutine set

    end subroutine set_input_1

    subroutine set_input_2( self, which, i, param )
        class(simple_program), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        integer,                       intent(in)    :: i
        type(simple_input_param),      intent(in)    :: param
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                write(*,*) 'which field selector: ', trim(which)
                stop 'unsupported parameter field; simple_user_interface :: simple_program :: set_input_2'
        end select

        contains

            subroutine set( arr, i )
                type(simple_input_param), intent(inout) :: arr(i)
                integer,                  intent(in)  :: i
                allocate(arr(i)%key,               source=trim(param%key))
                allocate(arr(i)%keytype,           source=trim(param%keytype))
                allocate(arr(i)%descr_short,       source=trim(param%descr_short))
                allocate(arr(i)%descr_long,        source=trim(param%descr_long))
                allocate(arr(i)%descr_placeholder, source=trim(param%descr_placeholder))
                arr(i)%required = param%required
            end subroutine set

    end subroutine set_input_2

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
        use json_module
        use simple_strings, only: int2str
        class(simple_program), intent(in) :: self
        type(json_core)           :: json
        type(json_value), pointer :: pjson, program
        ! JSON init
        call json%initialize()
        call json%create_object(pjson,'')
        call json%create_object(program, trim(self%name))
        call json%add(pjson, program)
        ! program section
        call json%add(program, 'name',        self%name)
        call json%add(program, 'descr_short', self%descr_short)
        call json%add(program, 'descr_long',  self%descr_long)
        call json%add(program, 'executable',  self%executable)
        ! all sections
        call create_section( 'image input/output',     self%img_ios )
        call create_section( 'parameter input/output', self%parm_ios )
        call create_section( 'alternative inputs',     self%alt_ios )
        call create_section( 'search controls',        self%srch_ctrls )
        call create_section( 'filter controls',        self%filt_ctrls )
        call create_section( 'mask controls',          self%mask_ctrls )
        call create_section( 'computer controls',      self%comp_ctrls )
        ! write & clean
        call json%print(pjson, trim(adjustl(self%name))//'.json')
        if( json%failed() )then
            print *, 'json input/output error for program: ', trim(self%name)
            stop
        endif
        call json%destroy(pjson)

        contains

            subroutine create_section( name, arr )
                use simple_strings, only: split, parsestr
                character(len=*),          intent(in) :: name
                type(simple_input_param), allocatable, intent(in) :: arr(:)
                type(json_value), pointer :: entry, section, options
                character(len=STDLEN)     :: options_str, before
                character(len=KEYLEN)     :: args(8)
                integer                   :: i, j, sz, nargs
                logical :: found, param_is_multi, param_is_binary, exception
                call json%create_array(section, trim(name))
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        call json%create_object(entry, trim(arr(i)%key))
                        call json%add(entry, 'key', trim(arr(i)%key))
                        call json%add(entry, 'keytype', trim(arr(i)%keytype))
                        call json%add(entry, 'descr_short', trim(arr(i)%descr_short))
                        call json%add(entry, 'descr_long', trim(arr(i)%descr_long))
                        call json%add(entry, 'descr_placeholder', trim(arr(i)%descr_placeholder))
                        call json%add(entry, 'required', arr(i)%required)
                        param_is_multi  = trim(arr(i)%keytype).eq.'multi'
                        param_is_binary = trim(arr(i)%keytype).eq.'binary'
                        if( param_is_multi .or. param_is_binary )then
                            options_str = trim(arr(i)%descr_placeholder)
                            call split( options_str, '(', before )
                            call split( options_str, ')', before )
                            call parsestr(before, '|', args, nargs)
                            exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 3)
                            if( exception )then
                                write(*,*)'Poorly formatted options string for entry ', trim(arr(i)%key)
                                write(*,*)trim(arr(i)%descr_placeholder)
                                stop
                            endif
                            call json%add(entry, 'options', args(1:nargs))
                            do j = 1, nargs
                                call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                            enddo
                        endif
                        call json%add(section, entry)
                    enddo
                endif
                call json%add(pjson, section)
            end subroutine create_section

    end subroutine write2json

    subroutine new_refine3D
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',& ! name
        &'3D volume refinement',&                                                                          ! descr_short
        &'is a distributed workflow for 3D volume refinement based on probabilistic projection matching',& ! descr_long
        &'simple_distr_exec',&                                                                             ! executable
        &3, 5, 0, 12, 7, 5, 2)                                                                              ! # entries in each group

        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%set_input('img_ios', 1, stk)
        call refine3D%set_input('img_ios', 2, stktab)
        call refine3D%set_input('img_ios', 3, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. recvol.mrc', .false.)
        ! parameter input/output
        call refine3D%set_input('parm_ios', 1, ctf)
        call refine3D%set_input('parm_ios', 2, smpd)
        call refine3D%set_input('parm_ios', 3, phaseplate)
        call refine3D%set_input('parm_ios', 4, deftab)
        call refine3D%set_input('parm_ios', 5, oritab)
        ! alternative inputs
        !<empty>
        ! search controls
        call refine3D%set_input('srch_ctrls', 1, nspace)
        call refine3D%set_input('srch_ctrls', 2, startit)
        call refine3D%set_input('srch_ctrls', 3, trs)
        call refine3D%set_input('srch_ctrls', 4, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false.)
        call refine3D%set_input('srch_ctrls', 5, maxits)
        call refine3D%set_input('srch_ctrls', 6, update_frac)
        call refine3D%set_input('srch_ctrls', 7, frac)
        call refine3D%set_input('srch_ctrls', 8, pgrp)
        call refine3D%set_input('srch_ctrls', 9, 'nnn', 'num', 'Number of nearest neighbours', 'Number of nearest projection direction &
        &neighbours in neigh=yes refinement', '# projection neighbours{10% of search space}', .false.)
        call refine3D%set_input('srch_ctrls', 10, 'nstates', 'num', 'Number of states', 'Number of conformational/compositional states to reconstruct',&
        '# states to reconstruct', .false.)
        call refine3D%set_input('srch_ctrls', 11, objfun)
        call refine3D%set_input('srch_ctrls', 12, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc|single|multi|greedy_single|greedy_multi|cluster|&
        &clustersym){no}', '(snhc|single|multi|greedy_single|greedy_multi|cluster|clustersym){no}', .false.)
        ! filter controls
        call refine3D%set_input('filt_ctrls', 1, hp)
        call refine3D%set_input('filt_ctrls', 2, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false.)
        call refine3D%set_input('filt_ctrls', 3, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false.)
        call refine3D%set_input('filt_ctrls', 4, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false.)
        call refine3D%set_input('filt_ctrls', 5, 'lplim_crit', 'num', 'Low-pass limit FSC criterion', 'FSC criterion for determining the low-pass limit(0.143-0.5){0.3}',&
        &'low-pass FSC criterion(0.143-0.5){0.3}', .false.)
        call refine3D%set_input('filt_ctrls', 6, 'eo', 'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for &
        &filtering and resolution estimation(yes|no){yes}', '(yes|no){yes}', .false.)
        call refine3D%set_input('filt_ctrls', 7, 'weights3D', 'binary', 'Spectral weighting', 'Weighted particle contributions based on &
        &the median FRC between the particle and its corresponding reference(yes|no){no}', '(yes|no){no}', .false.)
        ! mask controls
        call refine3D%set_input('mask_ctrls', 1, msk)
        call refine3D%set_input('mask_ctrls', 2, inner)
        call refine3D%set_input('mask_ctrls', 3, mskfile)
        call refine3D%set_input('mask_ctrls', 4, 'focusmsk', 'num', 'Mask radius in focused refinement', 'Mask radius in pixels for application of a soft-edged circular &
        &mask to remove background noise in focused refinement', 'focused mask radius in pixels', .false.)
        call refine3D%set_input('mask_ctrls', 5, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false. )
        ! computer controls
        call refine3D%set_input('comp_ctrls', 1, nparts)
        call refine3D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_refine3D

    subroutine new_cluster2D
        ! PROGRAM SPECIFICATION
        call cluster2D%new(&
        &'cluster2D',& ! name
        &'Simultaneous 2D alignment and clustering of single-particle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm adopted from the prime3D &
        &probabilistic ab initio 3D reconstruction algorithm',&                 ! descr_long
        &'simple_distr_exec',&                                                  ! executable
        &3, 5, 0, 9, 7, 2, 2)                                                   ! # entries in each group

        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster2D%set_input('img_ios', 1, stk)
        call cluster2D%set_input('img_ios', 2, stktab)
        call cluster2D%set_input('img_ios', 3, 'refs', 'file', 'Initial references',&
        &'Initial 2D references used to bootstrap the search', 'xxx.mrc file containing references', .false.)
        ! parameter input/output
        call cluster2D%set_input('parm_ios', 1, ctf)
        call cluster2D%set_input('parm_ios', 2, smpd)
        call cluster2D%set_input('parm_ios', 3, phaseplate)
        call cluster2D%set_input('parm_ios', 4, deftab)
        call cluster2D%set_input('parm_ios', 5, oritab)
        ! alternative inputs
        !<empty>
        ! search controls
        call cluster2D%set_input('srch_ctrls', 1, 'ncls', 'num', 'Number of 2D clusters', 'Number of groups to sort the particles &
        &into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .true.)
        call cluster2D%set_input('srch_ctrls', 2, startit)
        call cluster2D%set_input('srch_ctrls', 3, trs)
        call cluster2D%set_input('srch_ctrls', 4, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Initial/Final low-pass limits control the degree of down-scaling(yes|no){yes}',&
        &'(yes|no){yes}', .false.)
        call cluster2D%set_input('srch_ctrls', 5, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false.)
        call cluster2D%set_input('srch_ctrls', 6, 'dyncls', 'binary', 'Dynamic reallocation of clusters', 'Dynamic reallocation of clusters &
        &that fall below a minimum population by randomization(yes|no){yes}', '(yes|no){yes}', .false.)
        call cluster2D%set_input('srch_ctrls', 7, maxits)
        call cluster2D%set_input('srch_ctrls', 8, update_frac)
        call cluster2D%set_input('srch_ctrls', 9, frac)
        ! filter controls
        call cluster2D%set_input('filt_ctrls', 1, hp)
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
        call cluster2D%set_input('mask_ctrls', 1, msk)
        call cluster2D%set_input('mask_ctrls', 2, inner)
        ! computer controls
        call cluster2D%set_input('comp_ctrls', 1, nparts)
        call cluster2D%set_input('comp_ctrls', 2, nthr)
    end subroutine new_cluster2D

    subroutine new_initial_3Dmodel
        ! PROGRAM SPECIFICATION
        call initial_3Dmodel%new(&
        &'initial_3Dmodel',& ! name
        &'3D abinitio model generation from class averages',&                           ! descr_short
        &'is a distributed workflow for generating an initial 3D model from class'&
        &' averages obtained with cluster2D',&                                          ! descr_long
        &'simple_distr_exec',&                                                          ! executable
        &1, 1, 0, 9, 3, 3, 2)                                                           ! # entries in each group

        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call initial_3Dmodel%set_input('img_ios', 1, stk)
        ! parameter input/output
        call initial_3Dmodel%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        !<empty>
        ! search controls
        call initial_3Dmodel%set_input('srch_ctrls', 1, nspace)
        call initial_3Dmodel%set_input('srch_ctrls', 2, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false.)
        call initial_3Dmodel%set_input('srch_ctrls', 3, maxits)
        call initial_3Dmodel%set_input('srch_ctrls', 4, update_frac)
        call initial_3Dmodel%set_input('srch_ctrls', 5, frac)
        call initial_3Dmodel%set_input('srch_ctrls', 6, pgrp)
        call initial_3Dmodel%set_input('srch_ctrls', 7, 'pgrp_known', 'binary', 'Point-group applied directly', 'Point-group applied direclty rather than first doing a reconstruction &
        &in c1 and searching for the symmerty axis(yes|no){no}', '(yes|no){no}', .false.)
        call initial_3Dmodel%set_input('srch_ctrls', 8, objfun)
        call initial_3Dmodel%set_input('srch_ctrls', 9, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated convergence rate. Final low-pass limit controls the degree of down-scaling(yes|no){yes}','(yes|no){yes}', .false.)
        ! filter controls
        call initial_3Dmodel%set_input('filt_ctrls', 1, hp)
        call initial_3Dmodel%set_input('filt_ctrls', 2, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms', .false.)
        call initial_3Dmodel%set_input('filt_ctrls', 3, 'lpstop',  'num', 'Final low-pass limit',   'Final low-pass limit',   'low-pass limit in Angstroms', .false.)
        ! mask controls
        call initial_3Dmodel%set_input('mask_ctrls', 1, msk)
        call initial_3Dmodel%set_input('mask_ctrls', 2, inner)
        call initial_3Dmodel%set_input('mask_ctrls', 3, 'width', 'num', 'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false. )
        ! computer controls
        call initial_3Dmodel%set_input('comp_ctrls', 1, nparts)
        call initial_3Dmodel%set_input('comp_ctrls', 2, nthr)
    end subroutine new_initial_3Dmodel

    subroutine new_postprocess
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',& ! name
        &'is a program for post-processing of volumes',&                            ! descr_short
        &'Use program volops to estimate the B-factor with the Guinier plot',&      ! descr_long
        &'simple_exec',&                                                            ! executable
        &1, 1, 0, 0, 7, 9, 1)                                                       ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess%set_input('img_ios', 1, 'vol1', 'file', 'Volume', 'Volume to post-process &
        & sections for particle image matching', 'input volume e.g. recvol.mrc', .true.)
        ! parameter input/output
        call postprocess%set_input('parm_ios', 1, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%set_input('filt_ctrls', 1, hp)
        call postprocess%set_input('filt_ctrls', 2, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false.)
        call postprocess%set_input('filt_ctrls', 3, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false.)
        call postprocess%set_input('filt_ctrls', 4, 'vol_filt', 'file', 'Input filter volume', 'Input filter volume',&
        & 'input filter volume e.g. aniso_optlp_state01.mrc', .false.)
        call postprocess%set_input('filt_ctrls', 5, 'fsc', 'file', 'Binary file with FSC info', 'Binary file with FSC info&
        & for filtering', 'input binary file e.g. fsc_state01.bin', .false.)
        call postprocess%set_input('filt_ctrls', 6, 'bfac', 'num', 'B-factor for sharpening',&
        &'B-factor for sharpening in Angstroms^2', 'B-factor in Angstroms^2', .false.)
        call postprocess%set_input('filt_ctrls', 7, 'mirr', 'multi', 'Perform mirroring',&
        &'Whether to mirror and along which axis(no|x|y){no}', '(no|x|y){no}', .false.)
        ! mask controls
        call postprocess%set_input('mask_ctrls', 1, msk)
        call postprocess%set_input('mask_ctrls', 2, inner)
        call postprocess%set_input('mask_ctrls', 3, mskfile)
        call postprocess%set_input('mask_ctrls', 4, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false.)
        call postprocess%set_input('mask_ctrls', 5, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for enevloppe mask generation', 'Volume threshold', .false.)
        call postprocess%set_input('mask_ctrls', 6, 'automsk', 'binary', 'Perform envelope masking',&
        &'Whether to generate an envelope mask(yes|no){no}', '(yes|no){no}', .false.)
        call postprocess%set_input('mask_ctrls', 7, 'mw', 'num', 'Molecular weight','Molecular weight in kDa', 'in kDa', .false.)
        call postprocess%set_input('mask_ctrls', 8, 'width', 'num', 'Inner mask falloff',&
        &'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge', .false. )
        call postprocess%set_input('mask_ctrls', 9, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false. )
        ! computer controls
        call postprocess%set_input('comp_ctrls', 1, nthr)
    end subroutine new_postprocess

    subroutine new_extract
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', & ! name
        &'is a program that extracts particle images from integrated movies',&                  ! descr_short
        &'Boxfiles are assumed to be in EMAN format but we provide a conversion script'//&      ! descr long
        &' (relion2emanbox.pl) for *.star files containing particle coordinates obtained&
        & with Relion. In addition to one single-particle image stack per micrograph the&
        & program produces a parameter files that should be concatenated for use in&
        & conjunction with other SIMPLE programs.',&
        &'simple_exec',&                                                                        ! executable
        &0, 6, 3, 0, 0, 0, 0)                                                                   ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract%set_input('parm_ios', 1, smpd)
        call extract%set_input('parm_ios', 2, 'unidoc', 'file', 'Unified resources doc', 'Unified resources doc&
        & mapping micrographs, box files & CTF parameters', 'input unidoc e.g. unidoc_001.txt', .false.)
        call extract%set_input('parm_ios', 3, 'ctffind_doc', 'file', 'CTFFIND CTF parameters', 'list of per-micrograph &
        & CTFFIND CTF parameters to transfer', 'list input *.txt', .false.)
        call extract%set_input('parm_ios', 4, 'box', 'num', 'Box size', 'Square box size in pixels', 'in pixels', .false.)
        call extract%set_input('parm_ios', 5, 'pcontrast', 'binary', 'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false.)
        call extract%set_input('parm_ios', 6, 'outside', 'binary', 'Extract outside boundaries', 'Extract boxes outside the micrograph boundaries(yes|no){no}', '(yes|no){no}', .false.)
        ! alternative inputs
        call extract%set_input('alt_ios', 1, 'filetab', 'file', 'Micrographs list',&
        &'List of integrated micrographs', 'list input e.g. mics.txt', .false.)
        call extract%set_input('alt_ios', 2, 'boxtab', 'file', 'List of boxes',&
        &'List of single-particles boxes (EMAN format)', 'list input e.g. boxes.txt', .false.)
        call extract%set_input('alt_ios', 3, 'dir', 'string', 'Ouput directory',&
        &'Ouput directory for single-particle images & CTF parameters', '...', .false.)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_extract

    subroutine new_scale
        ! PROGRAM SPECIFICATION
        call scale%new(&
        &'scale', & ! name
        &'is a program for re-scaling MRC & SPIDER stacks and volumes',&                        ! descr_short
        &'is a program for re-scaling, clipping and padding MRC & SPIDER stacks and volumes',&  ! descr_long
        &'simple_exec',&                                                                        ! executable
        &3, 5, 3, 0, 0, 1, 1)                                                                   ! # entries in each group
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call scale%set_input('img_ios', 1, stk)
        scale%img_ios(1)%required = .false.
        call scale%set_input('img_ios', 2, 'vol1', 'file', 'Input volume', 'Input volume to re-scale',&
        &'input volume e.g. recvol.mrc', .false.)
        call scale%set_input('img_ios', 3, 'filetab', 'file', 'Stacks list',&
        &'List of stacks of images to rescale', 'list input e.g. stktab.txt', .false.)
        ! parameter input/output
        call scale%set_input('parm_ios', 1, smpd)
        call scale%set_input('parm_ios', 2, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'in pixels', .false.)
        call scale%set_input('parm_ios', 3, 'scale', 'num', 'Scaling ratio', 'Target box ratio for scaling(0-1)', '(0-1)', .false.)
        call scale%set_input('parm_ios', 4, 'scale2', 'num', 'Second scaling ratio', 'Second target box ratio for scaling(0-1)', '(0-1)', .false.)
        call scale%set_input('parm_ios', 5, 'clip', 'num', 'Clipped box size', 'Target box size for clipping in pixels', 'in pixels', .false.)
        ! alternative inputs
        call scale%set_input('alt_ios', 1, 'outvol', 'file', 'Output volume name', 'Output volume name', 'e.g. outvol.mrc', .false.)
        call scale%set_input('alt_ios', 2, 'outstk', 'file', 'Output stack name', 'Output images stack name', 'e.g. outstk.mrc', .false.)
        call scale%set_input('alt_ios', 3, 'outstk2', 'file', 'Second output stack name', 'Second output images stack name', 'e.g. outstk2.mrc', .false.)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call scale%set_input('mask_ctrls', 1, msk)
        scale%mask_ctrls(1)%required = .false.
        ! computer controls
        call scale%set_input('comp_ctrls', 1, nthr)
    end subroutine new_scale

end module simple_user_interface
