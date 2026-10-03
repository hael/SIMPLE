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
        call cmd_dict%push('maxits',        'maximum # iterations')
        call cmd_dict%push('mirr',          'mirror(no|x|y){no}')
        call cmd_dict%push('mskdiam',       'mask diameter(in A)')
        call cmd_dict%push('ncls',          '# clusters')
        call cmd_dict%push('neg',           'invert contrast of images(yes|no)')
        call cmd_dict%push('nparts',        '# partitions in distributed exection')
        call cmd_dict%push('nptcls',        '# images in stk/# orientations in oritab')
        call cmd_dict%push('nrots',         '# rotations in 2D analysis{0}')
        call cmd_dict%push('nstates',       '# states to reconstruct')
        call cmd_dict%push('nthr',          '# OpenMP threads{1}')
        call cmd_dict%push('oritype',       'SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D|projinfo|jobproc|compenv)')
        call cmd_dict%push('outfile',       'output document')
        call cmd_dict%push('outstk',        'output image stack')
        call cmd_dict%push('pickrefs',       'MRC stack of picking references')
        call cmd_dict%push('projfile',      'SIMPLE *.simple project file')
        call cmd_dict%push('refs',          'initial2Dreferences.ext')
        call cmd_dict%push('stk',           'particle stack with all images(ptcls.ext)')
        call cmd_dict%push('trust_header',  'trust the smpd value in the image header(yes|no){no}')
        call cmd_dict%push('vol1',          'input volume no1(invol1.ext)')
        call cmd_dict%push('which_iter',    'iteration nr')
    end subroutine init_cmd_dict

    subroutine new_private_prgs
        private_prgs(:)%nreq = 0
        private_prgs(:)%nopt = 0

        ! CAVGASSEMBLE, for assembling class averages
        call private_prgs(1)%set_name('cavgassemble')
        ! required keys
        call private_prgs(1)%push_req_key('projfile')
        call private_prgs(1)%push_req_key('nparts')
        call private_prgs(1)%push_req_key('ncls')
        ! optional keys
        call private_prgs(1)%push_opt_key('nthr')
        call private_prgs(1)%push_opt_key('refs')

        ! CHECK_BOX
        call private_prgs(2)%set_name('check_box')
        ! optional keys
        call private_prgs(2)%push_opt_key('stk')
        call private_prgs(2)%push_opt_key('vol1')

        ! CHECK_NPTCLS
        call private_prgs(3)%set_name('check_nptcls')
        ! required keys
        call private_prgs(3)%push_req_key('stk')

        ! EXPORT_CAVGS
        call private_prgs(4)%set_name('export_cavgs')
        ! required keys
        call private_prgs(4)%push_req_key('projfile')
        ! optional keys
        call private_prgs(4)%push_opt_key('outstk')

        ! KSTEST, Kolmogorov-Smirnov test to deduce equivalence or
        ! non-equivalence between two distributions in a non-parametric manner
        call private_prgs(5)%set_name('kstest')
        ! required keys
        call private_prgs(5)%push_req_key('infile')
        call private_prgs(5)%push_req_key('infile2')

        ! MAKE_PICKREFS, for preparing templates for particle picking
        call private_prgs(6)%set_name('make_pickrefs')
        ! optional keys
        call private_prgs(6)%push_opt_key('nthr')
        call private_prgs(6)%push_opt_key('pickrefs')
        call private_prgs(6)%push_opt_key('neg')
        call private_prgs(6)%push_opt_key('nrots')
        call private_prgs(6)%push_opt_key('mirr')
        call private_prgs(6)%push_opt_key('trust_header')

        ! PRINT_PROJECT_VALS, for printing specific values in a project file field
        call private_prgs(7)%set_name('print_project_vals')
        ! required keys
        call private_prgs(7)%push_req_key('projfile')
        call private_prgs(7)%push_req_key('keys')
        call private_prgs(7)%push_req_key('oritype')

        ! RANK_CAVGS, for ranking class averages
        call private_prgs(8)%set_name('rank_cavgs')
        ! required keys
        call private_prgs(8)%push_req_key('projfile')
        call private_prgs(8)%push_req_key('stk')
        ! set optional keys
        call private_prgs(8)%push_opt_key('outstk')

        ! ROTMATS2ORIS, for converting a text file (9 records per line) describing rotation matrices into a SIMPLE oritab
        call private_prgs(9)%set_name('rotmats2oris')
        ! required keys
        call private_prgs(9)%push_req_key('infile')
        ! optional keys
        call private_prgs(9)%push_opt_key('outfile')
        call private_prgs(9)%push_opt_key('oritype')

        ! VOLASSEMBLE, for asssembling subvolumes generated in distributed execution
        call private_prgs(10)%set_name('volassemble')
        ! required keys
        call private_prgs(10)%push_req_key('nparts')
        call private_prgs(10)%push_req_key('projfile')
        call private_prgs(10)%push_req_key('mskdiam')
        ! optional keys
        call private_prgs(10)%push_opt_key('nthr')
        call private_prgs(10)%push_opt_key('nstates')
        call private_prgs(10)%push_opt_key('which_iter')

        ! CALC_PSPEC, for asssembling power spectra for refine3D
        call private_prgs(11)%set_name('calc_pspec')
        ! required keys
        call private_prgs(11)%push_req_key('nparts')
        call private_prgs(11)%push_req_key('projfile')
        call private_prgs(11)%push_req_key('nthr')

        ! CALC_GROUP_SIGMAS, for asssembling sigmas for refine3D
        call private_prgs(12)%set_name('calc_group_sigmas')
        ! required keys
        call private_prgs(12)%push_req_key('nparts')
        call private_prgs(12)%push_req_key('projfile')
        call private_prgs(12)%push_req_key('nthr')
        call private_prgs(12)%push_req_key('which_iter')

        ! Pearson's correlation coefficient
        call private_prgs(13)%set_name('pearsn')
        ! required keys
        call private_prgs(13)%push_req_key('infile')
        call private_prgs(13)%push_req_key('infile2')

        ! check stochastic update scheme
        call private_prgs(14)%set_name('check_stoch_update')
        ! required keys
        call private_prgs(14)%push_req_key('maxits')
        call private_prgs(14)%push_req_key('nptcls')

        ! check fractional update scheme
        call private_prgs(15)%set_name('check_update_frac')
        ! required keys
        call private_prgs(15)%push_req_key('nptcls')

        ! suggest picking references
        call private_prgs(16)%set_name('shape_rank_cavgs')
        ! required keys
        call private_prgs(16)%push_req_key('projfile')
        ! optional keys
        call private_prgs(16)%push_opt_key('nthr')

        n_private_prgs = 16
    end subroutine new_private_prgs

end module simple_private_prgs
