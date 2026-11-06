!@descr: private program interface defintions (those executed by simple_private_exec)
module simple_private_prgs
use simple_core_module_api
implicit none

public :: make_private_ui, print_private_cmdline, get_private_keys_required, get_n_private_keys_required
private
#include "simple_local_flags.inc"

! private program instance
type simple_private_prg
    private
    type(string)          :: name
    character(len=KEYLEN) :: keys_required(MAXNKEYS)
    character(len=KEYLEN) :: keys_optional(MAXNKEYS)
    integer               :: nreq=0, nopt=0
  contains
    procedure :: set_name
    procedure :: push_req_key
    procedure :: push_opt_key
end type simple_private_prg

! array of simple_private_exec program specifications
integer, parameter       :: NMAX_PRIVATE_PRGS = 100
integer                  :: n_private_prgs    = 0
type(simple_private_prg) :: private_prgs(NMAX_PRIVATE_PRGS)

! command-line dictionary describing the keys of the private programs (usage text)
integer, parameter :: NMAX_CMD_DICT = 300
type(chash)        :: cmd_dict

contains

    ! instance methods

    subroutine set_name( self, name )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: name
        self%name = trim(name)
    end subroutine set_name

    subroutine push_req_key( self, key )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: key
        self%nreq = self%nreq + 1
        self%keys_required(self%nreq) = trim(key)
    end subroutine push_req_key

    subroutine push_opt_key( self, key )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: key
        self%nopt = self%nopt + 1
        self%keys_optional(self%nopt) = trim(key)
    end subroutine push_opt_key

    ! class methods

    subroutine make_private_ui
        call init_cmd_dict
        call new_private_prgs
    end subroutine make_private_ui

    function get_n_private_keys_required( prg ) result( nreq )
        character(len=*), intent(in)  :: prg
        integer :: nreq, iprg, i
        nreq = 0
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name%to_char() )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nreq = private_prgs(iprg)%nreq
    end function get_n_private_keys_required

    function get_private_keys_required( prg ) result( keys_required )
        character(len=*), intent(in)  :: prg
        type(string), allocatable     :: keys_required(:)
        integer :: iprg, i, nreq
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name%to_char() )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nreq = private_prgs(iprg)%nreq
        if( nreq == 0 ) return
        allocate(keys_required(nreq))
        do i=1,nreq
            keys_required(i) =trim(private_prgs(iprg)%keys_required(i))
        end do
    end function get_private_keys_required

    subroutine print_private_cmdline( prg )
        character(len=*), intent(in) :: prg
        character(len=KEYLEN), allocatable :: sorted_keys(:)
        integer :: iprg, i, nreq, nopt
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name%to_char() )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 )then
            THROW_WARN(trim(prg)//' lacks description in the private_prgs class')
            return
        endif
        write(logfhandle,'(a)') 'USAGE:'
        write(logfhandle,'(a)') 'bash-3.2$ simple_private_exec prg=simple_program key1=val1 key2=val2 ...'
        ! print required
        nreq = private_prgs(iprg)%nreq
        if( nreq > 0 )then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'REQUIRED'
            allocate(sorted_keys(nreq), source=private_prgs(iprg)%keys_required(:nreq))
            call lex_sort(sorted_keys)
            call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
            deallocate(sorted_keys)
        endif
        ! print optionals
        nopt = private_prgs(iprg)%nopt
        if( nopt > 0 )then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'OPTIONAL'
            allocate(sorted_keys(nopt), source=private_prgs(iprg)%keys_optional(:nopt))
            call lex_sort(sorted_keys)
            call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
            deallocate(sorted_keys)
        endif
        write(logfhandle,'(a)') ''
    end subroutine print_private_cmdline

    subroutine init_cmd_dict
        call cmd_dict%new(NMAX_CMD_DICT)
        call cmd_dict%push('infile',        'file with inputs(.txt)')
        call cmd_dict%push('infile2',       'file with inputs(.txt)')
        call cmd_dict%push('keys',          'keys of values to print')
        call cmd_dict%push('lp',            'low-pass limit(in A)')
        call cmd_dict%push('maxits',        'maximum # iterations')
        call cmd_dict%push('mirr',          'mirror(no|x|y){no}')
        call cmd_dict%push('mskdiam',       'mask diameter(in A)')
        call cmd_dict%push('ncls',          '# clusters')
        call cmd_dict%push('neg',           'invert contrast of images(yes|no)')
        call cmd_dict%push('nparts',        '# partitions in distributed exection')
        call cmd_dict%push('nptcls',        '# images in stk/# orientations in oritab')
        call cmd_dict%push('nrots',         '# rotations in 2D analysis{0}')
        call cmd_dict%push('nspace',        '# projection directions')
        call cmd_dict%push('nstates',       '# states to reconstruct')
        call cmd_dict%push('nthr',          '# OpenMP threads{1}')
        call cmd_dict%push('oritype',       'SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D|projinfo|jobproc|compenv)')
        call cmd_dict%push('outfile',       'output document')
        call cmd_dict%push('outstk',        'output image stack')
        call cmd_dict%push('projfile_target', 'SIMPLE *.simple project file to merge into')
        call cmd_dict%push('scale',         'image scale factor')
        call cmd_dict%push('smpd',          'sampling distance(in A)')
        call cmd_dict%push('dir_target',    'output directory')
        call cmd_dict%push('pickrefs',       'MRC stack of picking references')
        call cmd_dict%push('projfile',      'SIMPLE *.simple project file')
        call cmd_dict%push('refs',          'initial2Dreferences.ext')
        call cmd_dict%push('stk',           'particle stack with all images(ptcls.ext)')
        call cmd_dict%push('trs',           'maximum halfwidth shift(in pixels)')
        call cmd_dict%push('trust_header',  'trust the smpd value in the image header(yes|no){no}')
        call cmd_dict%push('vol1',          'input volume no1(invol1.ext)')
        call cmd_dict%push('vol2',          'input volume no2(invol2.ext)')
        call cmd_dict%push('which_iter',    'iteration nr')
    end subroutine init_cmd_dict

    !> the next entry of the table; push_req and push_opt add to it (n_private_prgs counts them)
    subroutine add_private_prg( name )
        character(len=*), intent(in) :: name
        type(string) :: errmsg
        integer      :: i
        do i=1,n_private_prgs
            if( private_prgs(i)%name%to_char() == trim(name) )then
                errmsg = 'private program '//trim(name)//' has two entries; add_private_prg'
                THROW_HARD(errmsg%to_char())
            endif
        end do
        if( n_private_prgs == NMAX_PRIVATE_PRGS ) THROW_HARD('raise NMAX_PRIVATE_PRGS; add_private_prg')
        n_private_prgs = n_private_prgs + 1
        call private_prgs(n_private_prgs)%set_name(name)
    end subroutine add_private_prg

    subroutine push_req( key )
        character(len=*), intent(in) :: key
        call private_prgs(n_private_prgs)%push_req_key(key)
    end subroutine push_req

    subroutine push_opt( key )
        character(len=*), intent(in) :: key
        call private_prgs(n_private_prgs)%push_opt_key(key)
    end subroutine push_opt

    !> every program of simple_private_exec without a public descriptor has an entry here
    !! (scripts/check_cli_chain.py); a program with one is checked by its descriptor
    subroutine new_private_prgs
        private_prgs(:)%nreq = 0
        private_prgs(:)%nopt = 0
        n_private_prgs       = 0

        ! CAVGASSEMBLE, for assembling class averages
        call add_private_prg('cavgassemble')
        ! required keys
        call push_req('projfile')
        call push_req('nparts')
        call push_req('ncls')
        ! optional keys
        call push_opt('nthr')
        call push_opt('refs')

        ! CHECK_BOX
        call add_private_prg('check_box')
        ! optional keys
        call push_opt('stk')
        call push_opt('vol1')

        ! CHECK_NPTCLS
        call add_private_prg('check_nptcls')
        ! required keys
        call push_req('stk')

        ! EXPORT_CAVGS
        call add_private_prg('export_cavgs')
        ! required keys
        call push_req('projfile')
        ! optional keys
        call push_opt('outstk')

        ! KSTEST, Kolmogorov-Smirnov test to deduce equivalence or
        ! non-equivalence between two distributions in a non-parametric manner
        call add_private_prg('kstest')
        ! required keys
        call push_req('infile')
        call push_req('infile2')

        ! MAKE_PICKREFS, for preparing templates for particle picking
        call add_private_prg('make_pickrefs')
        ! optional keys
        call push_opt('nthr')
        call push_opt('pickrefs')
        call push_opt('neg')
        call push_opt('nrots')
        call push_opt('mirr')
        call push_opt('trust_header')

        ! PRINT_PROJECT_VALS, for printing specific values in a project file field
        call add_private_prg('print_project_vals')
        ! required keys
        call push_req('projfile')
        call push_req('keys')
        call push_req('oritype')

        ! RANK_CAVGS, for ranking class averages
        call add_private_prg('rank_cavgs')
        ! required keys
        call push_req('projfile')
        call push_req('stk')
        ! set optional keys
        call push_opt('outstk')

        ! ROTMATS2ORIS, for converting a text file (9 records per line) describing rotation matrices into a SIMPLE oritab
        call add_private_prg('rotmats2oris')
        ! required keys
        call push_req('infile')
        ! optional keys
        call push_opt('outfile')
        call push_opt('oritype')

        ! VOLASSEMBLE, for asssembling subvolumes generated in distributed execution
        call add_private_prg('volassemble')
        ! required keys
        call push_req('nparts')
        call push_req('projfile')
        call push_req('mskdiam')
        ! optional keys
        call push_opt('nthr')
        call push_opt('nstates')
        call push_opt('which_iter')

        ! CALC_PSPEC, for asssembling power spectra for refine3D
        call add_private_prg('calc_pspec')
        ! required keys
        call push_req('nparts')
        call push_req('projfile')
        call push_req('nthr')

        ! COMMIT_SIGMA2, the commit of the canonical sigma2 state after a matcher pass
        call add_private_prg('commit_sigma2')
        ! required keys
        call push_req('nparts')
        call push_req('projfile')
        call push_req('nthr')
        call push_req('which_iter')

        ! Pearson's correlation coefficient
        call add_private_prg('pearsn')
        ! required keys
        call push_req('infile')
        call push_req('infile2')

        ! check stochastic update scheme
        call add_private_prg('check_stoch_update')
        ! required keys
        call push_req('maxits')
        call push_req('nptcls')

        ! check fractional update scheme
        call add_private_prg('check_update_frac')
        ! required keys
        call push_req('nptcls')

        ! suggest picking references
        call add_private_prg('shape_rank_cavgs')
        ! required keys
        call push_req('projfile')
        ! optional keys
        call push_opt('nthr')

        ! the project's class averages matched against state volumes vol1..volN (stream multistate 3D)
        call add_private_prg('assign_cavgs_states')
        ! required keys
        call push_req('projfile')
        call push_req('nstates')
        call push_req('vol1')
        call push_req('mskdiam')
        ! optional keys
        call push_opt('vol2')
        call push_opt('nthr')
        call push_opt('lp')
        call push_opt('trs')
        call push_opt('nspace')

        ! REFINE2D_DISTR, the distributed refine2D master of the stream's 2D pool
        call add_private_prg('refine2D_distr')
        ! required keys
        call push_req('projfile')
        ! optional keys
        call push_opt('nparts')
        call push_opt('nthr')

        ! CALC_PSPEC_ASSEMBLE, the assembly of the partitions' power spectra (sigma2 bootstrap)
        call add_private_prg('calc_pspec_assemble')
        ! required keys
        call push_req('projfile')
        call push_req('nparts')
        ! optional keys
        call push_opt('nthr')

        ! AGGREGATE_CHUNKS, the merge of two chunk projects and their class averages
        call add_private_prg('aggregate_chunks')
        ! required keys
        call push_req('projfile')
        call push_req('projfile_target')
        ! optional keys
        call push_opt('mskdiam')
        call push_opt('ncls')

        ! SCALE_PROJECT, a project's particles scaled to a new box (written into a new project)
        call add_private_prg('scale_project')
        ! required keys
        call push_req('projfile')
        ! optional keys
        call push_opt('scale')
        call push_opt('smpd')
        call push_opt('dir_target')
    end subroutine new_private_prgs

end module simple_private_prgs
