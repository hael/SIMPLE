module simple_sp_project
#include "simple_lib.f08"
use simple_oris,    only: oris
use simple_binoris, only: binoris
implicit none

public :: sp_project
private

integer, parameter :: MAXN_OS_SEG = 13
character(len=4)   :: NULL = 'null'

type sp_project
    ! ORIS REPRESENTATIONS OF BINARY FILE SEGMENTS
    ! segments 1-10 reserved for simple program outputs, orientations and files
    type(oris)        :: os_stk    ! per-micrograph stack os, segment 1
    type(oris)        :: os_ptcl2D ! per-particle 2D os,      segment 2
    type(oris)        :: os_cls2D  ! per-cluster 2D os,       segment 3
    type(oris)        :: os_cls3D  ! per-cluster 3D os,       segment 4
    type(oris)        :: os_ptcl3D ! per-particle 3D os,      segment 5

    !...

    ! ARRAY REPRESENTATIONS OF BINARY FILE SEGMENTS FOR FRCS & FSCS
    real, allocatable :: frcs(:,:) ! Fourier Ring  Corrs      segment 9
    real, allocatable :: fscs(:,:) ! Fourier Shell Corrs      segment 10

    ! ORIS REPRESENTATIONS OF PROJECT DATA / DISTRIBUTED SYSTEM INFO / SYSTEM MANAGEMENT STUFF
    ! segments 11-20 reserved for project info, job management etc.
    type(oris)        :: projinfo  ! project information      segment 11
    type(oris)        :: jobproc   ! jobid + PID + etc.       segment 12
    type(oris)        :: compenv   ! computing environment    segment 13

    ! binary file-handler
    type(binoris) :: bos

    ! globals
    real    :: smpd, kv, cs, fraca
    integer :: nptcls, box
    logical :: l_phaseplate
contains
    ! field constructor
    procedure          :: new_seg_with_ptr
    ! field updaters
    procedure          :: update_projinfo
    procedure          :: update_compenv
    ! modifiers
    procedure          :: new_sp_oris
    procedure          :: set_sp_oris
    procedure          :: add_stk
    procedure, private :: projinfo2globals
    ! printers
    procedure          :: print_info
    ! I/O
    ! readers
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure          :: read_segment
    procedure, private :: segreader
    procedure, private :: read_2Darray_segment
    ! writers
    procedure          :: write
    procedure          :: write_segment
    procedure, private :: segwriter
    ! destructor
    procedure          :: kill
end type sp_project

contains

    ! private supporting subroutines / functions

    function which_flag2isgement( which ) result( isegment )
        character(len=*),  intent(in) :: which
        integer :: isegment
        select case(trim(which))
            case('stk')
                isegment = STK_SEG
            case('ptcl2D')
                isegment = PTCL2D_SEG
            case('cls2D')
                isegment = CLS2D_SEG
            case('cls3D')
                isegment = CLS3D_SEG
            case('ptcl3D')
                isegment = PTCL3D_SEG
            case('frcs')
                isegment = FRCS_SEG
            case('fscs')
                isegment = FSCS_SEG
            case('projinfo')
                isegment = PROJINFO_SEG
            case('jobproc')
                isegment = JOBPROC_SEG
            case('compenv')
                isegment = COMPENV_SEG
            case DEFAULT
                stop 'unsupported which flag; sp_project :: which_flag2isgement'
        end select
    end function which_flag2isgement

    ! field constructor

    subroutine new_seg_with_ptr( self, n, oritype, os_ptr )
        class(sp_project), target, intent(inout) :: self
        integer,                   intent(in)    :: n
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer,      intent(inout) :: os_ptr
        select case(trim(oritype))
            case('stk')
                call self%os_stk%new_clean(n)
                os_ptr => self%os_stk
            case('ptcl2D')
                call self%os_ptcl2D%new_clean(n)
                os_ptr => self%os_ptcl2D
            case('cls2D')
                call self%os_cls2D%new_clean(n)
                os_ptr => self%os_cls2D
            case('cls3D')
                call self%os_cls2D%new_clean(n)
                os_ptr => self%os_cls2D
            case('ptcl3D')
                call self%os_ptcl3D%new_clean(n)
                os_ptr => self%os_ptcl3D
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype)
                stop 'unsupported oritype; sp_project :: new_with_os_segptr'
        end select
    end subroutine new_seg_with_ptr

    ! field updaters

    subroutine update_projinfo( self, cline )
        use simple_cmdline, only: cmdline
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        character(len=:), allocatable :: projname_old
        character(len=STDLEN)         :: projname_new, projfile, projname, cwd
        character(len=3)              :: phaseplate
        if( self%projinfo%get_noris() > 1 )then
            ! no need to construct field
        else
            call self%projinfo%new_clean(1)
        endif
        ! projname & profile
        if( self%projinfo%isthere('projname') )then
            if( cline%defined('projname') )then
                projname_new = cline%get_carg('projname')
                call self%projinfo%getter(1, 'projname', projname_old)
                write(*,*) 'Changing project name from ', trim(projname_old), ' to ', trim(projname_new)
                call self%projinfo%set(1, 'projname', trim(projname_new))
                call self%projinfo%set(1, 'projfile', trim(projname_new)//'.simple')
            endif
        else
            if( .not. cline%defined('projname') .and. .not. cline%defined('projfile') )then
                stop 'ERROR, the project needs a name, inputted via projname or projfile!'
            endif
            if( cline%defined('projfile') )then
                projfile = cline%get_carg('projfile')
                select case(fname2format(projfile))
                    case('O')
                        call self%projinfo%set(1, 'projfile', trim(projfile) )
                    case DEFAULT
                        write(*,*) 'Inputted projfile: ', trim(projfile)
                        stop 'has unsupported format'
                end select
                projname = get_fbody(projfile, '.simple')
                call self%projinfo%set(1, 'projname', trim(projname))
            endif
            if( cline%defined('projname') )then
                projname = cline%get_carg('projname')
                call self%projinfo%set(1, 'projname', trim(projname))
                call self%projinfo%set(1, 'projfile', trim(projname)//'.simple')
            endif
        endif
        ! hard requirements
        if( .not. cline%defined('smpd')     ) stop 'smpd (sampling distance in A) input required to create new project; sp_project :: new'
        if( .not. cline%defined('kv')       ) stop 'kv (acceleration voltage in kV{300}) input required to create new project; sp_project :: new'
        if( .not. cline%defined('cs')       ) stop 'cs (spherical aberration constant in mm{2.7}) input required to create new project; sp_project :: new'
        if( .not. cline%defined('fraca')    ) stop 'fraca (fraction of amplitude contrast{0.1}) input required to create new project; sp_project :: new'
        self%smpd  = cline%get_rarg('smpd')
        call self%projinfo%set(1, 'smpd',  self%smpd )
        self%kv    = cline%get_rarg('kv')
        call self%projinfo%set(1, 'kv',    self%kv   )
        self%cs    = cline%get_rarg('cs')
        call self%projinfo%set(1, 'cs',    self%cs   )
        self%fraca = cline%get_rarg('fraca')
        call self%projinfo%set(1, 'fraca', self%fraca )
        ! phaseplate flag is optional
        if( cline%defined('phaseplate') )then
            phaseplate        = cline%get_carg('phaseplate')
            self%l_phaseplate = trim(phaseplate) .eq. 'yes'
            call self%projinfo%set(1, 'phaseplate', trim(phaseplate))
        else
            self%l_phaseplate = .false.
            call self%projinfo%set(1, 'phaseplate', 'no')
        endif
        ! it is assumed that the project is created in the "project directory", i.e. stash cwd
        call simple_getcwd(cwd)
        call self%projinfo%set(1, 'cwd', trim(cwd))
    end subroutine update_projinfo

    subroutine update_compenv( self, cline )
        use simple_cmdline, only: cmdline
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        character(len=STDLEN)         :: env_var
        character(len=:), allocatable :: projname
        if( self%compenv%get_noris() > 1 )then
            ! no need to construct field
        else
            call self%compenv%new_clean(1)
        endif
        ! compenv has to be filled as strings as it is used as a string only dictionnary
        ! get from environment
        env_var = trim(simple_getenv('SIMPLE_PATH'))
        if( env_var.eq.'' )then
            write(*,*) 'ERROR! SIMPLE_PATH is not defined in your shell environment!'
            write(*,*) 'Please refer to installation documentation for correct system configuration'
            stop
        else
            call self%compenv%set(1, 'simple_path', trim(env_var))
        endif
        env_var = trim(simple_getenv('SIMPLE_QSYS'))
        if( env_var.eq.'' )then
            stop 'SIMPLE_QSYS is not defined in your environment.'
        else
            call self%compenv%set(1, 'qsys_name', trim(env_var))
        endif
        env_var = trim(simple_getenv('SIMPLE_EMAIL'))
        if( env_var.eq.'' ) env_var = 'my.name@uni.edu'
        call self%compenv%set(1, 'user_email', trim(env_var))
        ! get from command line
        if( cline%defined('time_per_image') )then
            call self%compenv%set(1, 'time_per_image', real2str(cline%get_rarg('time_per_image')))
        else
            if( .not. self%compenv%isthere('time_per_image') )then
                call self%compenv%set(1, 'time_per_image', int2str(TIME_PER_IMAGE_DEFAULT))
            endif
        endif
        if( cline%defined('user_account') )then
            call self%compenv%set(1, 'user_account', cline%get_carg('user_account'))
        else
            if( .not. self%compenv%isthere('user_account') )then
                call self%compenv%set(1, 'user_account', NULL)
            endif
        endif
        if( cline%defined('user_project') )then
            call self%compenv%set(1, 'user_project', cline%get_carg('user_project'))
        else
            if( .not. self%compenv%isthere('user_project') )then
                call self%compenv%set(1, 'user_project', NULL)
            endif
        endif
        if( cline%defined('qsys_partition') )then
            call self%compenv%set(1, 'qsys_partition', cline%get_carg('qsys_partition'))
        else
            if( .not. self%compenv%isthere('qsys_partition') )then
                call self%compenv%set(1, 'qsys_partition', NULL)
            endif
        endif
        if( cline%defined('qsys_qos') )then
            call self%compenv%set(1, 'qsys_qos', cline%get_carg('qsys_qos'))
        else
            if( .not. self%compenv%isthere('qsys_qos') )then
                call self%compenv%set(1, 'qsys_qos', NULL)
            endif
        endif
        if( cline%defined('qsys_reservation') )then
            call self%compenv%set(1, 'qsys_reservation', cline%get_carg('qsys_reservation'))
        else
            if( .not. self%compenv%isthere('qsys_reservation') )then
                call self%compenv%set(1, 'qsys_reservation', NULL)
            endif
        endif
        if( .not. self%compenv%isthere('job_name') )then
            call self%projinfo%getter(1, 'projname', projname)
            call self%compenv%set(1, 'job_name', 'simple_'//trim(projname) )
        endif
        if( cline%defined('job_memory_per_task') )then
            call self%compenv%set(1, 'job_memory_per_task', real2str(cline%get_rarg('job_memory_per_task')) )
        else
            if( .not. self%compenv%isthere('job_memory_per_task') )then
                call self%compenv%set(1, 'job_memory_per_task', int2str(JOB_MEMORY_PER_TASK_DEFAULT) )
            endif
        endif
    end subroutine update_compenv

    ! modifiers

    subroutine new_sp_oris( self, which, n )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer,           intent(in)    :: n
        select case(trim(which))
            case('stk')
                call self%os_stk%new_clean(n)
            case('ptcl2D')
                call self%os_ptcl2D%new_clean(n)
            case('cls2D')
                call self%os_cls2D%new_clean(n)
            case('cls3D')
                call self%os_cls3D%new_clean(n)
            case('ptcl3D')
                call self%os_ptcl3D%new_clean(n)
            case('projinfo')
                call self%projinfo%new_clean(n)
            case('jobproc')
                call self%jobproc%new_clean(n)
            case('compenv')
                call self%compenv%new_clean(n)
            case DEFAULT
                stop 'unsupported which flag; sp_project :: new_sp_oris'
        end select
    end subroutine new_sp_oris

    subroutine set_sp_oris( self, which, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        class(oris),       intent(inout) :: os
        select case(trim(which))
            case('stk')
                self%os_stk    = os
            case('ptcl2D')
                self%os_ptcl2D = os
            case('cls2D')
                self%os_cls2D  = os
            case('cls3D')
                self%os_cls3D  = os
            case('ptcl3D')
                self%os_ptcl3D = os
            case('projinfo')
                self%projinfo  = os
            case('jobproc')
                self%jobproc   = os
            case('compenv')
                self%compenv   = os
            case DEFAULT
                stop 'unsupported which flag; sp_project :: set_sp_oris'
        end select
    end subroutine set_sp_oris

    subroutine add_stk( self, stk, os )
        use simple_imghead, only: find_ldim_nptcls
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: stk
        class(oris),       intent(in)    :: os ! parameters associated with stk
        integer :: ldim(3), nptcls, box, n_os, n_os_stk, n_os_ptcl2D, n_os_ptcl3D
        ! file exists?
        if( .not. file_exists(stk) )then
            write(*,*) 'Inputted stack (stk): ', trim(stk)
            stop 'does not exist in cwd; sp_project :: add_stk'
        endif
        ! check if field is already populated
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk > 1 )then
            write(*,*) 'stack field (self%os_stk) already populated with # entries: ', n_os_stk
            stop 'ABORTING, this mode is not intended for data appendation!'
        endif
        ! find dimension of inputted stack and compare with os
        call find_ldim_nptcls(stk, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(*,*) 'xdim: ', ldim(1)
            write(*,*) 'ydim: ', ldim(2)
            stop 'nonsquare particle images not supported; sp_project :: add_stk'
        endif
        ! check that inputs are of conforming sizes
        n_os = os%get_noris()
        if( n_os /= nptcls )then
            write(*,*) '# input oris: ', n_os
            write(*,*) '# ptcl imgs in stk: ', nptcls
            stop 'ERROR, nonconforming sizes of inputs; sp_project :: add_stk'
        endif
        ! make stk field
        call self%os_stk%new_clean(1)
        call self%os_stk%set(1, 'stk',    trim(stk))
        call self%os_stk%set(1, 'box',    real(ldim(1)))
        call self%os_stk%set(1, 'nptcls', real(nptcls))
        ! update globals
        self%nptcls = nptcls
        self%box    = box
        ! update particle fields
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        if( n_os_ptcl2D > 1 )then
            write(*,*) 'ptcl2D field (self%os_ptcl2D) already populated with # entries: ', n_os_ptcl2D
            stop 'ABORTING, this mode assumes empty particle fields in project file'
        endif
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_ptcl3D > 1 )then
            write(*,*) 'ptcl3D field (self%os_ptcl3D) already populated with # entries: ', n_os_ptcl3D
            stop 'ABORTING, this mode assumes empty particle fields in project file'
        endif
        self%os_ptcl2D = os
        self%os_ptcl3D = os
    end subroutine add_stk

    subroutine projinfo2globals( self )
        class(sp_project), intent(inout) :: self
        character(len=:), allocatable    :: phaseplate
        if( self%projinfo%get_noris() > 0 )then
            if( self%projinfo%isthere('smpd')       ) self%smpd  = self%projinfo%get(1, 'smpd'  )
            if( self%projinfo%isthere('kv')         ) self%kv    = self%projinfo%get(1, 'kv'    )
            if( self%projinfo%isthere('cs')         ) self%cs    = self%projinfo%get(1, 'cs'    )
            if( self%projinfo%isthere('fraca')      ) self%fraca = self%projinfo%get(1, 'fraca' )
            if( self%projinfo%isthere('phaseplate') )then
                call self%projinfo%getter(1, 'phaseplate', phaseplate)
                self%l_phaseplate = trim(phaseplate) .eq. 'yes'
            endif
        endif
    end subroutine projinfo2globals

    ! printers

    subroutine print_info( self )
        class(sp_project), intent(in) :: self
        integer :: n
        n = self%os_stk%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in per-micrograph stack segment (1) :', n
        n = self%os_ptcl2D%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in per-particle 2D      segment (2) :', n
        n = self%os_cls2D%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in per-cluster  2D      segment (3) :', n
        n = self%os_cls3D%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in per-cluster  3D      segment (4) :', n
        n = self%os_ptcl3D%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in per-particle 3D      segment (5) :', n
        n = size(self%frcs,1)
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in FRCs                 segment (9) :', n
        n = size(self%fscs,1)
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in FSCs                 segment (9) :', n
        n = self%projinfo%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in project info         segment (11):', n
        n = self%jobproc%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in jobproc              segment (12):', n
        n = self%compenv%get_noris()
        if( n > 1 ) write(*,'(a,1x,i10)') '# entries in compenv              segment (12):', n
    end subroutine print_info

    ! readers

    subroutine read( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment, n
        if( .not. file_exists(trim(fname)) )then
            write(*,*) 'fname: ', trim(fname)
            stop 'inputted file does not exist; sp_project :: read'
        endif
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: read'
        endif
        call self%bos%open(fname)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment)
        end do
        call self%bos%close
        call self%projinfo2globals
    end subroutine read

    subroutine read_ctfparams_state_eo( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment, n
        if( .not. file_exists(trim(fname)) )then
            write(*,*) 'fname: ', trim(fname)
            stop 'inputted file does not exist; sp_project :: read'
        endif
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: read'
        endif
        call self%bos%open(fname)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment, only_ctfparams_state_eo=.true.)
        end do
        call self%bos%close
        call self%projinfo2globals
    end subroutine read_ctfparams_state_eo

    subroutine read_segment( self, which, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        integer :: isegment
        if( .not. file_exists(trim(fname)) )then
            write(*,*) 'fname: ', trim(fname)
            stop 'inputted file does not exist; sp_project :: read_segment'
        endif
        select case(fname2format(fname))
            case('O')
                ! *.simple project file
                isegment = which_flag2isgement(which)
                call self%bos%open(fname)
                call self%segreader(isegment)
                call self%bos%close
            case('T')
                ! *.txt plain text ori file
                select case(trim(which))
                    case('stk')
                        call self%os_stk%read(fname)
                    case('ptcl2D')
                        call self%os_ptcl2D%read(fname, fromto)
                    case('cls2D')
                        call self%os_cls2D%read(fname)
                    case('cls3D')
                        call self%os_cls3D%read(fname,  fromto)
                    case('ptcl3D')
                        call self%os_ptcl3D%read(fname, fromto)
                    case('projinfo')
                        call self%projinfo%read(fname)
                        call self%projinfo2globals
                    case('jobproc')
                        call self%jobproc%read(fname)
                    case('compenv')
                        call self%compenv%read(fname)
                    case DEFAULT
                        stop 'unsupported which flag; sp_project :: read_segment'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: read_segment'
        end select
    end subroutine read_segment

    subroutine segreader( self, isegment, only_ctfparams_state_eo )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: isegment
        logical, optional, intent(in)    :: only_ctfparams_state_eo
        integer :: n
        n = self%bos%get_n_records(isegment)
        select case(isegment)
            case(STK_SEG)
                call self%os_stk%new_clean(n)
                call self%bos%read_segment(isegment, self%os_stk,    only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(PTCL2D_SEG)
                call self%os_ptcl2D%new_clean(n)
                call self%bos%read_segment(isegment, self%os_ptcl2D, only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(CLS2D_SEG)
                call self%os_cls2D%new_clean(n)
                call self%bos%read_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%os_cls3D%new_clean(n)
                call self%bos%read_segment(isegment, self%os_cls3D)
            case(PTCL3D_SEG)
                call self%os_ptcl3D%new_clean(n)
                call self%bos%read_segment(isegment, self%os_ptcl3D, only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(FRCS_SEG)
                call self%read_2Darray_segment(FRCS_SEG, self%frcs)
            case(FSCS_SEG)
                call self%read_2Darray_segment(FSCS_SEG, self%fscs)
            case(PROJINFO_SEG)
                call self%projinfo%new_clean(n)
                call self%bos%read_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%jobproc%new_clean(n)
                call self%bos%read_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%compenv%new_clean(n)
                call self%bos%read_segment(isegment, self%compenv)
        end select
    end subroutine segreader

    subroutine read_2Darray_segment( self, isegment, array )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: isegment
        real, allocatable, intent(out)   :: array(:,:)
        real    :: rval
        integer :: ndim1, ndim2
        ndim1 = self%bos%get_n_records(isegment)
        ndim2 = self%bos%get_n_bytes_per_record(isegment) / sizeof(rval)
        if( allocated(array) ) deallocate(array)
        allocate( array(ndim1,ndim2), source=0. )
        call self%bos%read_segment(isegment, array)
    end subroutine read_2Darray_segment

    ! writers

    subroutine write( self, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*), optional, intent(in) :: fname
        integer,          optional, intent(in) :: fromto(2)
        character(len=:), allocatable :: projfile
        integer :: isegment
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: write'
            endif
            projfile = trim(fname)
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        call self%bos%open(projfile, del_if_exists=.true.)
        do isegment=1,MAXN_OS_SEG
            call self%segwriter(isegment, fromto)
        end do
        ! update header
        call self%bos%write_header
        call self%bos%close
    end subroutine write

    subroutine write_segment( self, which, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        integer :: isegment
        select case(fname2format(fname))
            case('O')
                stop 'write_segment is not supported for *.simple project files; sp_project :: write_segment'
            case('T')
                ! *.txt plain text ori file
                select case(trim(which))
                    case('stk')
                        if( self%os_stk%get_noris() > 0 )then
                            call self%os_stk%write(fname)
                        else
                            write(*,*) 'WARNING, no stk-type oris available to write; sp_project :: write_segment'
                        endif
                    case('ptcl2D')
                        if( self%os_ptcl2D%get_noris() > 0 )then
                            call self%os_ptcl2D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl2D-type oris available to write; sp_project :: write_segment'
                        endif
                    case('cls2D')
                        if( self%os_cls2D%get_noris() > 0 )then
                            call self%os_cls2D%write(fname)
                        else
                            write(*,*) 'WARNING, no cls2D-type oris available to write; sp_project :: write_segment'
                        endif
                    case('cls3D')
                        if( self%os_cls3D%get_noris() > 0 )then
                            call self%os_cls3D%write(fname,  fromto)
                        else
                            write(*,*) 'WARNING, no cls3D-type oris available to write; sp_project :: write_segment'
                        endif
                    case('ptcl3D')
                        if( self%os_ptcl3D%get_noris() > 0 )then
                            call self%os_ptcl3D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl3D-type oris available to write; sp_project :: write_segment'
                        endif
                    case('projinfo')
                        if( self%projinfo%get_noris() > 0 )then
                            call self%projinfo%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no projinfo-type oris available to write; sp_project :: write_segment'
                        endif
                    case('jobproc')
                        if( self%jobproc%get_noris() > 0 )then
                            call self%jobproc%write(fname)
                        else
                            write(*,*) 'WARNING, no jobproc-type oris available to write; sp_project :: write_segment'
                        endif
                    case('compenv')
                        if( self%compenv%get_noris() > 0 )then
                            call self%compenv%write(fname)
                        else
                            write(*,*) 'WARNING, no compenv-type oris available to write; sp_project :: write_segment'
                        endif
                    case DEFAULT
                        stop 'unsupported which flag; sp_project :: write_segment'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: write_segment'
        end select
    end subroutine write_segment

    subroutine segwriter( self, isegment, fromto )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: isegment
        integer, optional, intent(in)    :: fromto(2)
        logical :: fromto_present
        fromto_present = present(fromto)
        select case(isegment)
            case(STK_SEG)
                call self%bos%write_segment(isegment, self%os_stk)
            case(PTCL2D_SEG)
                if( fromto_present )then
                    call self%bos%write_segment(isegment, self%os_ptcl2D, fromto)
                else
                    call self%bos%write_segment(isegment, self%os_ptcl2D)
                endif
            case(CLS2D_SEG)
                call self%bos%write_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                if( fromto_present )then
                    call self%bos%write_segment(isegment, self%os_cls3D, fromto)
                else
                    call self%bos%write_segment(isegment, self%os_cls3D)
                endif
            case(PTCL3D_SEG)
                if( fromto_present )then
                    call self%bos%write_segment(isegment, self%os_ptcl3D, fromto)
                else
                    call self%bos%write_segment(isegment, self%os_ptcl3D)
                endif
            case(FRCS_SEG)
                call self%bos%write_segment(FRCS_SEG, self%frcs)
            case(FSCS_SEG)
                call self%bos%write_segment(FSCS_SEG, self%fscs)
            case(PROJINFO_SEG)
                call self%bos%write_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%bos%write_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%bos%write_segment(isegment, self%compenv)
        end select
    end subroutine segwriter

    ! destructor

    subroutine kill( self )
        class(sp_project), intent(inout) :: self
        call self%os_stk%kill
        call self%os_ptcl2D%kill
        call self%os_cls2D%kill
        call self%os_cls3D%kill
        call self%os_ptcl3D%kill
        call self%projinfo%kill
        call self%jobproc%kill
        call self%compenv%kill
    end subroutine kill

end module simple_sp_project
