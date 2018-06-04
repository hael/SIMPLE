module simple_sp_project
include 'simple_lib.f08'
use simple_ori,     only: ori
use simple_oris,    only: oris
use simple_binoris, only: binoris
implicit none

public :: sp_project, oritype2segment
private

integer, parameter :: MAXN_OS_SEG = 13

type sp_project
    ! ORIS REPRESENTATIONS OF BINARY FILE SEGMENTS
    ! segments 1-10 reserved for simple program outputs, orientations and files
    ! In segment 7 we stash class averages, ranked class averages, final volumes etc.
    type(oris)        :: os_mic    ! micrographs,              segment 1
    type(oris)        :: os_stk    ! per-micrograph stack os,  segment 2
    type(oris)        :: os_ptcl2D ! per-particle 2D os,       segment 3
    type(oris)        :: os_cls2D  ! per-cluster 2D os,        segment 4
    type(oris)        :: os_cls3D  ! per-cluster 3D os,        segment 5
    type(oris)        :: os_ptcl3D ! per-particle 3D os,       segment 6
    type(oris)        :: os_out    ! critical project outputs, segment 7

    ! ORIS REPRESENTATIONS OF PROJECT DATA / DISTRIBUTED SYSTEM INFO / SYSTEM MANAGEMENT STUFF
    ! segments 11-20 reserved for project info, job management etc.
    type(oris)        :: projinfo  ! project information      segment 11
    type(oris)        :: jobproc   ! jobid + PID + etc.       segment 12
    type(oris)        :: compenv   ! computing environment    segment 13

    ! binary file-handler
    type(binoris) :: bos
contains
    ! field constructor
    procedure          :: new_seg_with_ptr
    ! field updaters
    procedure          :: update_projinfo
    procedure          :: update_compenv
    procedure          :: append_project
    procedure          :: append_job_descr2jobproc
    ! index management
    procedure, private :: map_ptcl_ind2stk_ind
    procedure          :: map_cavgs_selection
    ! os_mic related methods
    procedure          :: add_single_movie
    procedure          :: add_movies
    procedure          :: get_movies_table
    procedure          :: get_micparams
    ! os_stk related methods
    procedure          :: add_stk
    procedure          :: add_stktab
    procedure          :: add_single_stk
    procedure          :: get_stkname
    procedure          :: get_stkname_and_ind
    procedure, private :: add_scale_tag
    ! os_out related methods
    procedure          :: add_cavgs2os_out
    procedure          :: add_fsc2os_out
    procedure          :: add_frcs2os_out
    procedure          :: add_vol2os_out
    procedure, private :: add_entry2os_out
    procedure          :: get_cavgs_stk
    procedure          :: get_vol
    procedure          :: get_fsc
    procedure          :: get_frcs
    procedure          :: get_nptcls_from_osout
    ! getters
    procedure          :: get_nptcls
    procedure          :: get_box
    procedure          :: get_smpd
    procedure          :: get_nmics
    procedure          :: get_nmovies
    procedure          :: get_nintgs
    procedure          :: get_ctfflag
    procedure          :: get_ctfflag_type
    procedure          :: has_phaseplate
    procedure          :: get_ctfparams
    procedure          :: is_virgin_field
    ! modifiers
    procedure          :: split_stk
    procedure          :: set_sp_oris
    procedure          :: scale_projfile
    procedure          :: merge_algndocs
    procedure          :: map2ptcls
    ! I/O
    ! printers
    procedure          :: print_info
    procedure          :: print_segment
    ! readers
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure          :: read_segment
    procedure, private :: segreader
    ! writers
    procedure          :: write
    procedure          :: write_segment2txt
    procedure          :: write_segment_inside
    procedure, private :: segwriter
    procedure          :: segwriter_inside
    ! destructor
    procedure          :: kill
end type sp_project

contains

    ! field constructor

    subroutine new_seg_with_ptr( self, n, oritype, os_ptr )
        class(sp_project), target, intent(inout) :: self
        integer,                   intent(in)    :: n
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer,      intent(inout) :: os_ptr
        select case(trim(oritype))
            case('mic')
                call self%os_mic%new(n)
                os_ptr => self%os_mic
            case('stk')
                call self%os_stk%new(n)
                os_ptr => self%os_stk
            case('ptcl2D')
                call self%os_ptcl2D%new(n)
                os_ptr => self%os_ptcl2D
            case('cls2D')
                call self%os_cls2D%new(n)
                os_ptr => self%os_cls2D
            case('cls3D')
                call self%os_cls3D%new(n)
                os_ptr => self%os_cls3D
            case('ptcl3D')
                call self%os_ptcl3D%new(n)
                os_ptr => self%os_ptcl3D
            case('out')
                call self%os_out%new(n)
                os_ptr => self%os_out
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype)
                stop 'unsupported oritype; sp_project :: new_seg_with_ptr'
        end select
    end subroutine new_seg_with_ptr

    ! field updaters

    subroutine update_projinfo( self, cline )
        use simple_cmdline, only: cmdline
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        character(len=:), allocatable :: projname
        character(len=STDLEN)         :: projfile, cwd
        if( self%projinfo%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%projinfo%new(1)
        endif
        ! projname & profile
        if( self%projinfo%isthere('projname') )then
            if( cline%defined('projname') )then
                projname = cline%get_carg('projname')
                call self%projinfo%set(1, 'projname', trim(projname))
                call self%projinfo%set(1, 'projfile', trim(projname)//'.simple')
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
                projname = get_fbody(projfile, 'simple')
                call self%projinfo%set(1, 'projname', trim(projname))
            endif
            if( cline%defined('projname') )then
                projname = cline%get_carg('projname')
                call self%projinfo%set(1, 'projname', trim(projname))
                call self%projinfo%set(1, 'projfile', trim(projname)//'.simple')
            endif
        endif
        ! it is assumed that the project is created in the root "project directory", i.e. stash cwd
        call simple_getcwd(cwd)
        call self%projinfo%set(1, 'cwd', trim(cwd))
    end subroutine update_projinfo

    subroutine update_compenv( self, cline )
        use simple_cmdline, only: cmdline
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        character(len=STDLEN)            :: env_var
        character(len=:), allocatable    :: projname
        integer :: iostat
        if( self%compenv%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%compenv%new(1)
        endif
        ! compenv has to be filled as strings as it is used as a string only dictionnary
        ! get from environment
        iostat  = simple_getenv('SIMPLE_PATH', env_var)
        if( iostat /= 0 )then
            write(*,*) 'ERROR! SIMPLE_PATH is not defined in your shell environment!'
            write(*,*) 'Please refer to installation documentation for correct system configuration'
            stop
        else
            call self%compenv%set(1, 'simple_path', trim(env_var))
        endif
        iostat  = simple_getenv('SIMPLE_QSYS', env_var)
        if( iostat /= 0 )then
            stop 'SIMPLE_QSYS is not defined in your environment.'
        else
            iostat  = simple_getenv('SIMPLE_QSYS', env_var)
            call self%compenv%set(1, 'qsys_name', trim(env_var))
        endif
        iostat = simple_getenv('SIMPLE_EMAIL', env_var)
        if( iostat/=0 ) env_var = 'my.name@uni.edu'
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
        endif
        if( cline%defined('user_project') )then
            call self%compenv%set(1, 'user_project', cline%get_carg('user_project'))
        endif
        if( cline%defined('qsys_partition') )then
            call self%compenv%set(1, 'qsys_partition', cline%get_carg('qsys_partition'))
        endif
        if( cline%defined('qsys_qos') )then
            call self%compenv%set(1, 'qsys_qos', cline%get_carg('qsys_qos'))
        endif
        if( cline%defined('qsys_reservation') )then
            call self%compenv%set(1, 'qsys_reservation', cline%get_carg('qsys_reservation'))
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

    !> append segment to current project. BOTH projects must be read in first!
    subroutine append_project( self, proj, oritype )
        class(sp_project), target, intent(inout) :: self, proj
        character(len=*),          intent(in)    :: oritype
        class(oris),          pointer :: os_ptr, os_append_ptr
        type(oris)                    :: os
        type(ctfparams)               :: ctfvar
        character(len=:), allocatable :: stk
        real                          :: smpd, smpd_self
        integer                       :: i, cnt, n, n2append
        select case(trim(oritype))
            case('mic')
                os_ptr => self%os_mic
                os_append_ptr => proj%os_mic
            case('stk')
                os_ptr => self%os_stk
                os_append_ptr => proj%os_stk
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype)
                stop 'unsupported oritype for this purpose; sp_project :: append_project'
        end select
        n2append = os_append_ptr%get_noris()
        if( n2append == 0 )return
        smpd = os_append_ptr%get(1, 'smpd')
        n    = os_ptr%get_noris()
        if( n == 0 )then
            ! first entry
        else
            smpd_self = os_ptr%get(1, 'smpd')
            if( abs(smpd-smpd_self) > 0.001 )then
                write(*,*) 'smpd self', smpd_self
                write(*,*) 'smpd 2 append', smpd
                stop ' Only a project with the same smpd can be appended to the project; simple_sp_project :: append_project'
            endif
        endif
        select case(trim(oritype))
            case('mic')
                if( n == 0 )then
                    os_ptr = os_append_ptr
                else
                    ! append
                    call os%new(n + n2append)
                    do i=1,n
                        call os%set_ori(i, os_ptr%get_ori(i))
                    enddo
                    cnt = n
                    do i=1,n2append
                        cnt = cnt + 1
                        call os%set_ori(cnt, os_append_ptr%get_ori(i))
                    enddo
                    os_ptr = os
                endif
            case('stk')
                ! this assumes there's only one stack in the project to append
                call os_append_ptr%getter(1, 'stk', stk)
                ctfvar = proj%get_ctfparams('ptcl2D', 1)
                call self%add_stk(stk, ctfvar)
        end select
    end subroutine append_project

    subroutine append_job_descr2jobproc( self, exec_dir, job_descr, did_update )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: exec_dir
        class(chash),      intent(inout) :: job_descr
        logical,           intent(out)   :: did_update
        character(len=:), allocatable :: edir
        type(ori)     :: o
        integer       :: njobs, ijob, ind
        character(8)  :: date
        character(10) :: time
        did_update = .true.
        njobs = self%jobproc%get_noris()
        if( njobs > 0 )then
            if( trim(exec_dir) .ne. './')then
                do ijob=1,njobs
                    if( self%jobproc%isthere(ijob, 'exec_dir') )then
                        call self%jobproc%getter(ijob, 'exec_dir', edir)
                        if( str_has_substr(exec_dir,edir) )then
                            ! we already have a job description for this exec dir stored
                            did_update = .false.
                            return
                        endif
                    endif
                end do
            endif
        endif
        ! store job description along with exec_dir & execution time
        call o%chash2ori(job_descr)
        call o%set('exec_dir', trim(exec_dir))
        call date_and_time(date=date, time=time)
        call o%set('date', date)
        call o%set('time', time)
        ! update jobproc field
        if( njobs > 0 )then
            ind = njobs + 1
            call self%jobproc%reallocate(ind)
        else
            call self%jobproc%new(1)
            ind = 1
        endif
        call self%jobproc%set_ori(ind,o)
    end subroutine append_job_descr2jobproc

    ! index management

    subroutine map_ptcl_ind2stk_ind( self, oritype, iptcl, stkind, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        integer,                   intent(out)   :: stkind
        integer,                   intent(out)   :: ind_in_stk
        class(oris), pointer                     :: ptcl_field
        integer :: nptcls, fromp, top
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype), ' is not supported by this method'
                stop 'sp_project :: map_ptcl_ind2stk_ind'
        end select
        nptcls = ptcl_field%get_noris()
        ! first sanity check, range
        if( iptcl < 1 .or. iptcl > nptcls )then
            print *, 'iptcl : ', iptcl
            print *, 'nptcls: ', nptcls
            stop 'iptcl index out of range; sp_project :: map_ptcl_ind2stk_ind'
        endif
        ! second sanity check, stack index present in ptcl_field
        if( .not. ptcl_field%isthere(iptcl, 'stkind') )then
            print *, 'iptcl: ', iptcl
            print *, 'ERROR, stkind not present in field: ', trim(oritype)
            stop 'sp_project :: map_ptcl_ind2stk_ind'
        endif
        stkind = nint(ptcl_field%get(iptcl, 'stkind'))
        ! third sanity check, particle index in range
        fromp = nint(self%os_stk%get(stkind, 'fromp'))
        top   = nint(self%os_stk%get(stkind, 'top'))
        if( iptcl < fromp .or. iptcl > top )then
            print *, 'iptcl            : ', iptcl
            print *, 'prange for micstk: ', fromp, top
            stop 'iptcl index out of micstk range; sp_project :: map_ptcl_ind2stk_ind'
        endif
        ! output index in stack
        ind_in_stk = iptcl - fromp + 1
    end subroutine map_ptcl_ind2stk_ind

    subroutine map_cavgs_selection( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer, allocatable :: pinds(:)
        integer :: icls, iptcl, sz_cls2D, sz_states, sz_cls3D, ncls
        real    :: rstate
        sz_cls2D  = self%os_cls2D%get_noris()
        sz_states = size(states)
        if( sz_cls2D /= sz_states )then
            write(*,*) 'size(cls2D): ', sz_cls2D
            write(*,*) 'sz_states  : ', sz_states
            write(*,*) 'ERROR! size(cls2D) not consistent with size(states), aborting...'
            stop 'simple_sp_project :: map_cavgs_selection'
        endif
        ! map selection to self%os_cls2D
        do icls=1,sz_cls2D
            call self%os_cls2D%set(icls, 'state', real(states(icls)))
        end do
        ! map selection to self%os_cls3D
        sz_cls3D = self%os_cls3D%get_noris()
        if( sz_cls3D /= sz_cls2D ) call self%os_cls3D%new(sz_cls2D)
        sz_cls3D = sz_cls2D
        do icls=1,sz_cls3D
            call self%os_cls3D%set(icls, 'state', real(states(icls)))
        end do
        ! map selection to self%os_ptcl2D
        ncls = sz_states
        if( self%os_ptcl2D%get_noris() > 0 )then
            do icls=1,ncls
                call self%os_ptcl2D%get_pinds(icls, 'class', pinds)
                if( allocated(pinds) )then
                    rstate = real(states(icls))
                    do iptcl=1,size(pinds)
                        call self%os_ptcl2D%set(iptcl, 'state', rstate)
                    end do
                    deallocate(pinds)
                endif
            end do
        endif
    end subroutine map_cavgs_selection

    ! os_mic related methods

    subroutine add_single_movie( self, moviename, ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: moviename
        type(ctfparams),           intent(in)    :: ctfvars
        class(oris),      pointer     :: os_ptr
        character(len=:), allocatable :: fname
        integer :: n_os_mic, ldim(3), nframes
        ! oris object pointer
        os_ptr => self%os_mic
        ! check that stk field is empty
        n_os_mic = os_ptr%get_noris()
        if( n_os_mic > 0 )then
            write(*,*) 'stack field (self%os_stk) already populated with # entries: ', n_os_mic
            stop 'ABORTING! sp_project :: add_single_movie'
        endif
        ! update ori
        call os_ptr%new(1)
        call simple_full_path(moviename, fname, 'simple_sp_project::add_single_movie')
        call find_ldim_nptcls(trim(fname), ldim, nframes)
        if( nframes <= 0 )then
            write(*,*) 'WARNING! # frames in movie ', trim(fname), ' <= zero, ommitting'
        else if( nframes > 1 )then
            call os_ptr%set(1, 'movie', trim(fname))
            call os_ptr%set(1, 'imgkind', 'movie')
            call os_ptr%set(1, 'nframes',    real(nframes))
        else
            call os_ptr%set(1, 'intg',  trim(fname))
            call os_ptr%set(1, 'imgkind', 'mic')
        endif
        ! updates segment
        call os_ptr%set(1, 'xdim',       real(ldim(1)))
        call os_ptr%set(1, 'ydim',       real(ldim(2)))
        call os_ptr%set(1, 'smpd',       ctfvars%smpd)
        call os_ptr%set(1, 'kv',         ctfvars%kv)
        call os_ptr%set(1, 'cs',         ctfvars%cs)
        call os_ptr%set(1, 'fraca',      ctfvars%fraca)
        if( ctfvars%l_phaseplate )then
            call os_ptr%set(1, 'phaseplate', 'yes')
        else
            call os_ptr%set(1, 'phaseplate', 'no')
        endif
        select case(ctfvars%ctfflag)
            case(0)
                call os_ptr%set(1, 'ctf', 'no')
            case(1)
                call os_ptr%set(1, 'ctf', 'yes')
            case(2)
                call os_ptr%set(1, 'ctf', 'flip')
            case DEFAULT
                write(*,*) 'ctfvars%ctfflag: ', ctfvars%ctfflag
                stop 'ERROR, unsupported ctfflag; sp_project :: add_single_movie'
        end select
    end subroutine add_single_movie

    ! subroutine add_movies( self, filetab, smpd, kv, cs, fraca, phaseplate )
    subroutine add_movies( self, filetab, ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: filetab
        type(ctfparams),           intent(in)    :: ctfvars
        class(oris),               pointer     :: os_ptr
        character(len=LONGSTRLEN), allocatable :: movienames(:)
        character(len=:),          allocatable :: name, moviename
        integer               :: imic, ldim(3), nframes, nmics, nprev_mics, cnt, ntot
        logical               :: is_movie
        ! file exists?
        if( .not. file_exists(filetab) )then
            write(*,*) 'Inputted movie list (filetab): ', trim(filetab)
            stop 'does not exist in cwd; sp_project :: add_movies'
        endif
        ! oris object pointer
        os_ptr => self%os_mic
        ! read movie names
        call read_filetable(filetab, movienames)
        nmics = size(movienames)
        ! update oris
        nprev_mics = os_ptr%get_noris()
        ntot       = nmics + nprev_mics
        if( nprev_mics == 0 )then
            call os_ptr%new(ntot)
        else
            call os_ptr%reallocate(ntot)
        endif
        cnt = 0
        do imic=nprev_mics + 1,ntot
            cnt = cnt + 1
            call simple_full_path(movienames(cnt), moviename, 'simple_sp_project::add_movies')
            call find_ldim_nptcls(trim(moviename), ldim, nframes)
            if( nframes <= 0 )then
                write(*,*) 'WARNING! # frames in movie ', trim(moviename), ' <= zero, ommitting'
                cycle
            else if( nframes > 1 )then
                call os_ptr%set(imic, 'movie', trim(moviename))
                call os_ptr%set(imic, 'imgkind', 'movie')
                is_movie = .true.
            else
                call os_ptr%set(imic, 'intg',  trim(moviename))
                call os_ptr%set(imic, 'imgkind', 'mic')
                is_movie = .false.
            endif
            ! updates segment
            call os_ptr%set(imic, 'xdim',    real(ldim(1)))
            call os_ptr%set(imic, 'ydim',    real(ldim(2)))
            call os_ptr%set(imic, 'nframes', real(nframes))
            call os_ptr%set(imic, 'smpd',    ctfvars%smpd)
            call os_ptr%set(imic, 'kv',      ctfvars%kv)
            call os_ptr%set(imic, 'cs',      ctfvars%cs)
            call os_ptr%set(imic, 'fraca',   ctfvars%fraca)
            call os_ptr%set(imic, 'state',   1.0) ! default on import
            if( ctfvars%l_phaseplate )then
                call os_ptr%set(imic, 'phaseplate', 'yes')
            else
                call os_ptr%set(imic, 'phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(0)
                    call os_ptr%set(imic, 'ctf', 'no')
                case(1)
                    call os_ptr%set(imic, 'ctf', 'yes')
                case(2)
                    call os_ptr%set(imic, 'ctf', 'flip')
                case DEFAULT
                    write(*,*) 'ctfvars%ctfflag: ', ctfvars%ctfflag
                    stop 'ERROR, unsupported ctfflag; sp_project :: add_movies'
            end select
            deallocate(moviename)
        enddo
        if( is_movie )then
            name = 'MOVIE(S)'
        else
            name = 'MICROGRAPH(S)'
        endif
        write(*,'(A13,I6,A1,A)')'>>> IMPORTED ', nmics,' ', trim(name)
        write(*,'(A20,A,A1,I6)')'>>> TOTAL NUMBER OF ', trim(name),':',ntot
    end subroutine add_movies

    subroutine get_movies_table( self, moviestab )
        class(sp_project),                      intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(out)   :: moviestab(:)
        character(len=:), allocatable :: imgkind, mov
        integer :: i,n,cnt
        if(allocated(moviestab))deallocate(moviestab)
        n = self%get_nmovies()
        if( n==0 )return
        allocate(moviestab(n))
        cnt = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere('imgkind'))then
                call self%os_mic%getter(i,'imgkind',imgkind)
                if( trim(imgkind).eq.'movie' )then
                    cnt = cnt + 1
                    call self%os_mic%getter(i,'movie',mov)
                    moviestab(cnt) = trim(mov)
                endif
            endif
        enddo
    end subroutine get_movies_table

    function get_micparams( self, imic ) result( ctfvars )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(ctfparams) :: ctfvars
        ctfvars = self%os_mic%get_ctfvars(imic)
    end function get_micparams

    ! os_stk related methods

    subroutine add_stk( self, stk, ctfvars )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars ! All CTF parameters associated with stk
        type(ori)                     :: o
        character(len=:), allocatable :: stk_abspath
        integer :: ldim(3), nptcls, n_os_stk, n_os_ptcl2D, n_os_ptcl3D
        integer :: i, fromp, top
        ! full path and existence check
        call simple_full_path(stk, stk_abspath, 'sp_project :: add_stk')
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_abspath, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(*,*) 'xdim: ', ldim(1)
            write(*,*) 'ydim: ', ldim(2)
            stop 'ERROR! nonsquare particle images not supported; sp_project :: add_stk'
        endif
        ! updates_fields
        n_os_stk    = self%os_stk%get_noris() + 1
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 1 )then
            call self%os_stk%new(1)
            call self%os_ptcl2D%new(nptcls)
            call self%os_ptcl3D%new(nptcls)
            fromp = 1
            top   = nptcls
        else
            ! stk
            if( .not.self%os_stk%isthere(n_os_stk-1,'top') )then
                stop 'FROMP/TOP keys should always be informed; simple_sp_project :: add_stk'
            endif
            call self%os_stk%reallocate(n_os_stk)
            ! 2d
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls)
            ! 3d
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls)
            fromp = nint(self%os_stk%get(n_os_stk-1,'top')) + 1
            top   = fromp + nptcls - 1
        endif
        ! updates oris_objects
        call self%os_stk%set(n_os_stk, 'stk',     trim(stk_abspath))
        call self%os_stk%set(n_os_stk, 'box',     real(ldim(1)))
        call self%os_stk%set(n_os_stk, 'nptcls',  real(nptcls))
        call self%os_stk%set(n_os_stk, 'fromp',   real(fromp))
        call self%os_stk%set(n_os_stk, 'top',     real(top))
        call self%os_stk%set(n_os_stk, 'stkkind', 'split')
        call self%os_stk%set(n_os_stk, 'imgkind', 'ptcl')
        call self%os_stk%set(n_os_stk, 'state',   1.0) ! default on import
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO)
                call self%os_stk%set(n_os_stk, 'ctf', 'no')
            case(CTFFLAG_YES)
                call self%os_stk%set(n_os_stk, 'ctf', 'yes')
            case(CTFFLAG_FLIP)
                call self%os_stk%set(n_os_stk, 'ctf', 'flip')
            case DEFAULT
                write(*,*) 'ctfvars%ctfflag: ', ctfvars%ctfflag
                stop 'ERROR, unsupported ctfflag; sp_project :: add_stk'
        end select
        call self%os_stk%set(n_os_stk, 'smpd',    ctfvars%smpd)
        call self%os_stk%set(n_os_stk, 'kv',      ctfvars%kv)
        call self%os_stk%set(n_os_stk, 'cs',      ctfvars%cs)
        call self%os_stk%set(n_os_stk, 'fraca',   ctfvars%fraca)
        if( ctfvars%l_phaseplate )then
            call self%os_stk%set(n_os_stk, 'phaseplate', 'yes')
        else
            call self%os_stk%set(n_os_stk, 'phaseplate', 'no')
        endif
        ! preprocessign / streaming adds pairs: one micrograph and one stack
        ! so this keeps track of the index in this setting
        if( self%os_mic%get_noris() == n_os_stk )then
            call self%os_stk%set(n_os_stk, 'micind',  real(n_os_stk))
        endif
        ! update particle oris objects
        do i = 1, nptcls
            call o%new
            call o%set('dfx',    ctfvars%dfx)
            call o%set('dfy',    ctfvars%dfy)
            call o%set('angast', ctfvars%angast)
            if( ctfvars%l_phaseplate ) call o%set('phshift', ctfvars%phshift)
            call o%set('stkind', real(n_os_stk))
            call o%set('state',1.) ! default on import
            call self%os_ptcl2D%set_ori(n_os_ptcl2D+i, o)
            call self%os_ptcl3D%set_ori(n_os_ptcl3D+i, o)
        enddo
    end subroutine add_stk

    subroutine add_single_stk( self, stk, ctfvars, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars ! CTF parameters associated with stk (smpd,kv,cs,fraca,phaseplate)
        class(oris),       intent(inout) :: os      ! parameters associated with stk (dfx,dfy,angast,phshift)
        integer :: n_os_stk, n_os_ptcl2D, n_os_ptcl3D, ldim(3), nptcls
        character(len=:), allocatable :: stk_abspath
        ! check that stk field is empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk > 0 )then
            write(*,*) 'stack field (self%os_stk) already populated with # entries: ', n_os_stk
            stop 'ABORTING! sp_project :: add_single_stk'
        endif
        ! check that particle fields are empty
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        if( n_os_ptcl2D > 0 )then
            write(*,*) 'ptcl2D field (self%os_ptcl2D) already populated with # entries: ', n_os_ptcl2D
            stop 'ABORTING! empty particle fields in project file assumed; sp_project :: add_single_stk'
        endif
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_ptcl3D > 0 )then
            write(*,*) 'ptcl3D field (self%os_ptcl3D) already populated with # entries: ', n_os_ptcl3D
            stop 'ABORTING! empty particle fields in project file assumed; sp_project :: add_single_stk'
        endif
        self%os_ptcl2D = os
        self%os_ptcl3D = os
        call self%os_ptcl2D%set_all2single('stkind', 1.)
        call self%os_ptcl3D%set_all2single('stkind', 1.)
        call self%os_ptcl2D%set_all2single('state',  1.) ! default on import
        call self%os_ptcl3D%set_all2single('state',  1.) ! default on import
        ! full path and existence check
        call simple_full_path(stk, stk_abspath, 'sp_project :: add_single_stk')
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_abspath, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(*,*) 'xdim: ', ldim(1)
            write(*,*) 'ydim: ', ldim(2)
            stop 'ERROR! nonsquare particle images not supported; sp_project :: add_single_stk'
        endif
        call self%os_stk%new(1)
        call self%os_stk%set(1, 'stk',     trim(stk_abspath))
        call self%os_stk%set(1, 'box',     real(ldim(1)))
        call self%os_stk%set(1, 'nptcls',  real(nptcls))
        call self%os_stk%set(1, 'fromp',   1.0)
        call self%os_stk%set(1, 'top',     real(nptcls))
        call self%os_stk%set(1, 'stkkind', 'single')
        call self%os_stk%set(1, 'imgkind', 'ptcl')
        call self%os_stk%set(1, 'smpd',    ctfvars%smpd)
        call self%os_stk%set(1, 'kv',      ctfvars%kv)
        call self%os_stk%set(1, 'cs',      ctfvars%cs)
        call self%os_stk%set(1, 'fraca',   ctfvars%fraca)
        call self%os_stk%set(1, 'state',   1.0) ! default on import
        if( ctfvars%l_phaseplate )then
            if( .not. os%isthere('phshift') ) stop 'ERROR! phaseplate=yes & input oris lack phshift; sp_project :: add_single_stk'
            call self%os_stk%set(1, 'phaseplate', 'yes')
        else
            call self%os_stk%set(1, 'phaseplate', 'no')
        endif
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO)
                call self%os_stk%set(1, 'ctf', 'no')
            case(CTFFLAG_YES)
                call self%os_stk%set(1, 'ctf', 'yes')
            case(CTFFLAG_FLIP)
                call self%os_stk%set(1, 'ctf', 'flip')
            case DEFAULT
                write(*,*) 'ctfvars%ctfflag: ', ctfvars%ctfflag
                stop 'ERROR, unsupported ctfflag; sp_project :: add_single_stk'
        end select
    end subroutine add_single_stk

    subroutine add_stktab( self, stktab, os )
        class(sp_project),   intent(inout) :: self
        character(len=*),    intent(in)    :: stktab
        class(oris),         intent(inout) :: os ! parameters associated with stktab
        type(ctfparams) :: ctfvars
        type(ori)       :: o_stk
        character(len=LONGSTRLEN), allocatable :: stknames(:)
        integer :: istk, ldim(3), ldim_here(3), nptcls, n_os, iptcl, nstks
        ! file exists?
        if( .not. file_exists(stktab) )then
            write(*,*) 'Inputted stack list (stktab): ', trim(stktab)
            stop 'does not exist in cwd; sp_project :: add_stktab'
        endif
        ! read micrograph stack names
        call read_filetable(stktab, stknames)
        nstks = size(stknames)
        ! check that inputs are of conforming sizes
        n_os = os%get_noris()
        if( n_os /= nstks )then
            write(*,*) '# input oris      : ', n_os
            write(*,*) '# stacks in stktab: ', nstks
            stop 'ERROR! nonconforming sizes of inputs; sp_project :: add_stktab'
        endif
        do istk=1,nstks
            if( .not.file_exists(stknames(istk)) )then
                write(*,*) 'Inputted stack: ', trim(stknames(istk))
                stop 'does not exist in cwd; sp_project :: add_stktab'
            endif
            o_stk = os%get_ori(istk)
            ! logical dimension management
            call find_ldim_nptcls(trim(stknames(istk)), ldim, nptcls)
            ldim(3) = 1
            if( istk == 1 )then
                ldim_here = ldim
            else
                if( .not. all(ldim_here == ldim) )then
                    write(*,*) 'micrograph stack #  : ', istk
                    write(*,*) 'stk name            : ', trim(stknames(istk))
                    write(*,*) 'ldim in object      : ', ldim_here
                    write(*,*) 'ldim read from stack: ', ldim
                    stop 'inconsistent logical dimensions; sp_project :: add_stktab'
                endif
            endif
            if( ldim(1) /= ldim(2) )then
                write(*,*) 'stk name: ', trim(stknames(istk))
                write(*,*) 'xdim:     ', ldim(1)
                write(*,*) 'ydim:     ', ldim(2)
                stop 'ERROR! nonsquare particle images not supported; sp_project :: add_stktab'
            endif
            ! check variable presence
            if( .not. o_stk%isthere('ctf') )    stop 'ERROR! ctf flag missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('smpd') )   stop 'ERROR! smpd missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('kv') )     stop 'ERROR! kv missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('cs') )     stop 'ERROR! cs missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('fraca') )  stop 'ERROR! fraca missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('dfx') )    stop 'ERROR! dfx missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('dfy') )    stop 'ERROR! dfy missing in os input; sp_project :: add_stktab'
            if( .not. o_stk%isthere('angast') ) stop 'ERROR! angast missing in os input; sp_project :: add_stktab'
            ctfvars = o_stk%get_ctfvars()
            call self%add_stk(stknames(istk), ctfvars)
        enddo
    end subroutine add_stktab

    !>  Only commits to disk when a change to the project is made
    subroutine split_stk( self, nparts )
        use simple_map_reduce, only: split_nobjs_even
        use simple_image,      only: image
        class(sp_project),     intent(inout) :: self
        integer,               intent(in)    :: nparts
        character(len=4),   parameter :: EXT = '.mrc'
        type(image)                   :: img
        type(ori)                     :: orig_stk
        character(len=:), allocatable :: stk, tmp_dir, imgkind, stkpart, dest_stkpart
        character(len=:), allocatable :: ctfstr
        character(len=STDLEN) :: cwd
        integer    :: parts(nparts,2), ind_in_stk, iptcl, cnt, istk, box, n_os_stk
        integer    :: nptcls, nptcls_part, numlen, status
        real       :: smpd, cs, kv, fraca
        if( nparts < 2 )return
        ! check that stk field is not empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk==0 )then
            stop 'No stack to split! sp_project :: split_single_stk'
        else if( n_os_stk >= nparts )then
            return
        endif
        smpd    = self%os_stk%get(1,'smpd')
        box     = nint(self%os_stk%get(1,'box'))
        ! copy prep
        call self%os_stk%getter(1,'imgkind', imgkind)
        nptcls  = self%get_nptcls()
        parts   = split_nobjs_even( nptcls, nparts )
        numlen  = len_trim(int2str(nparts))
        ! images copy
        call img%new([box,box,1], smpd)
        call simple_getcwd(cwd)
        tmp_dir = trim(cwd) // '/tmp_stacks/'
        call simple_mkdir(trim(tmp_dir))
        do istk = 1,nparts
            allocate(stkpart, source=tmp_dir//'stack_part'//int2str_pad(istk,numlen)//EXT)
            cnt = 0
            do iptcl = parts(istk,1), parts(istk,2)
                cnt = cnt + 1
                call self%get_stkname_and_ind( 'ptcl2D', iptcl, stk, ind_in_stk )
                call img%read(stk, ind_in_stk)
                call img%write(stkpart, cnt)
            enddo
            deallocate(stkpart)
        enddo
        call img%kill
        if( n_os_stk > 1 )then
            ! wipe previous stack parts
            do istk = 1,n_os_stk
                call self%os_stk%getter(istk,'stk', stkpart)
                call del_file(stkpart)
            enddo
        endif
        ! updates new stack parts
        orig_stk = self%os_stk%get_ori(1)
        call self%os_stk%getter(1,'ctf', ctfstr)
        cs    = self%os_stk%get(1,'cs')
        kv    = self%os_stk%get(1,'kv')
        fraca = self%os_stk%get(1,'fraca')
        call self%os_stk%new(nparts)
        call simple_mkdir(trim(STKPARTSDIR), status=status)
        do istk = 1,nparts
            allocate(stkpart, source=tmp_dir//'stack_part'//int2str_pad(istk,numlen)//EXT)
            allocate(dest_stkpart, source=trim(STKPARTFBODY)//int2str_pad(istk,numlen)//EXT)
            status = simple_rename(trim(stkpart), trim(dest_stkpart))
            deallocate(stkpart)
            call simple_full_path(dest_stkpart, stkpart, 'sp_project :: split_stk')
            nptcls_part = parts(istk,2)-parts(istk,1)+1
            call self%os_stk%set(istk, 'ctf',   ctfstr)
            call self%os_stk%set(istk, 'cs',    cs)
            call self%os_stk%set(istk, 'kv',    kv)
            call self%os_stk%set(istk, 'fraca', fraca)
            call self%os_stk%set(istk, 'stk',     trim(stkpart))
            call self%os_stk%set(istk, 'box',     real(box))
            call self%os_stk%set(istk, 'smpd',    smpd)
            call self%os_stk%set(istk, 'nptcls',  real(nptcls_part))
            call self%os_stk%set(istk, 'fromp',   real(parts(istk,1)))
            call self%os_stk%set(istk, 'top',     real(parts(istk,2)))
            call self%os_stk%set(istk, 'imgkind', trim(imgkind))
            call self%os_stk%set(istk, 'stkkind', 'split')
            do iptcl=parts(istk,1),parts(istk,2)
                call self%os_ptcl2D%set(iptcl,'stkind',real(istk))
                call self%os_ptcl3D%set(iptcl,'stkind',real(istk))
            enddo
            deallocate(stkpart, dest_stkpart)
        enddo
        call self%write
        call simple_rmdir(tmp_dir)
    end subroutine split_stk

    function get_stkname( self, imic ) result( stkname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        character(len=:), allocatable    :: stkname
        integer :: nmics
        nmics = self%os_stk%get_noris()
        if( imic < 1 .or. imic > nmics )then
            print *, 'imic : ', imic
            print *, 'nmics: ', nmics
            stop 'imic index out of range; sp_project :: get_stkname'
        endif
        call self%os_stk%getter(imic, 'stk', stkname)
    end function get_stkname

    subroutine get_stkname_and_ind( self, oritype, iptcl, stkname, ind_in_stk )
        class(sp_project), target,     intent(inout) :: self
        character(len=*),              intent(in)    :: oritype
        integer,                       intent(in)    :: iptcl
        character(len=:), allocatable, intent(out)   :: stkname
        integer,                       intent(out)   :: ind_in_stk
        integer :: stkind
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk )
        ! output name
        if( allocated(stkname) ) deallocate(stkname)
        call self%os_stk%getter(stkind, 'stk', stkname)
    end subroutine get_stkname_and_ind

    subroutine add_scale_tag( self, dir )
        class(sp_project),          intent(inout) :: self
        character(len=*), optional, intent(in)    :: dir
        character(len=:), allocatable :: ext, newname, stkname
        character(len=4) :: ext_out
        integer :: imic, nmics
        nmics = self%os_stk%get_noris()
        do imic=1,nmics
            call self%os_stk%getter(imic, 'stk', stkname)
            ext = fname2ext(trim(stkname))
            select case(fname2format(stkname))
                case('M','D','B')
                    ext_out = '.mrc'
                case('S')
                    ext_out = '.spi'
                case DEFAULT
                    write(*,*)'format: ', trim(ext)
                    call simple_stop('This file format is not supported by SIMPLE; simple_sp_project::add_scale_tag')
            end select
            if(present(dir))then
                newname = trim(dir)//basename(add2fbody(stkname, '.'//trim(ext), trim(SCALE_SUFFIX)))
            else
                newname = add2fbody(stkname, '.'//trim(ext), trim(SCALE_SUFFIX))
            endif
            newname = fname_new_ext(newname, ext_out(2:4))
            call self%os_stk%set(imic, 'stk', newname)
        end do
    end subroutine add_scale_tag

    ! os_out related methods

    subroutine add_cavgs2os_out( self, stk, smpd)
        class(sp_project),     intent(inout) :: self
        character(len=*),      intent(in)    :: stk
        real,                  intent(in)    :: smpd ! sampling distance of images in stk
        character(len=:), allocatable :: stk_abspath
        integer :: ldim(3), nptcls, ind
        ! full path and existence check
        call simple_full_path(stk, stk_abspath, 'sp_project :: add_cavgs2os_out')
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_abspath, ldim, nptcls)
        ! add os_out entry
        call self%add_entry2os_out('cavg', ind)
        ! fill-in field
        call self%os_out%set(ind, 'stk',     trim(stk_abspath))
        call self%os_out%set(ind, 'box',     real(ldim(1)))
        call self%os_out%set(ind, 'nptcls',  real(nptcls))
        call self%os_out%set(ind, 'fromp',   1.0)
        call self%os_out%set(ind, 'top',     real(nptcls))
        call self%os_out%set(ind, 'smpd',    real(smpd))
        call self%os_out%set(ind, 'stkkind', 'single')
        call self%os_out%set(ind, 'imgkind', 'cavg')
        call self%os_out%set(ind, 'ctf',     'no')
    end subroutine add_cavgs2os_out

    subroutine add_frcs2os_out( self, frc, which_imgkind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: frc, which_imgkind
        character(len=:), allocatable :: frc_abspath
        integer :: ind
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                write(*,*)'Invalid FRC kind: ', trim(which_imgkind)
                stop 'sp_project :: add_frcs2os_out'
        end select
        ! full path and existence check
        call simple_full_path(frc, frc_abspath, 'sp_project :: add_frcs2os_out')
        ! add os_out entry
        call self%add_entry2os_out(which_imgkind, ind)
        ! fill-in field
        call self%os_out%set(ind, 'frcs', trim(frc_abspath))
        call self%os_out%set(ind, 'imgkind', trim(which_imgkind))
    end subroutine add_frcs2os_out

    subroutine add_fsc2os_out( self, fsc, state, box)
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fsc
        integer,           intent(in)    :: state, box
        character(len=:), allocatable :: fsc_abspath, imgkind
        integer :: i, ind, n_os_out
        ! full path and existence check
        call simple_full_path(fsc, fsc_abspath, 'sp_project :: add_fsc2os_out')
        ! add os_out entry
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out)
        else
            ind = 0
            do i=1,n_os_out
                if( self%os_out%isthere(i,'imgkind') )then
                    call self%os_out%getter(i,'imgkind',imgkind)
                    if(trim(imgkind).eq.'fsc')then
                        if( self%os_out%get_state(i) == state )then
                            ind = i
                            exit
                        endif
                    endif
                endif
            end do
            if( ind == 0 )then
                n_os_out = n_os_out + 1
                call self%os_out%reallocate(n_os_out)
                ind = n_os_out
            endif
        endif
        ! fill-in field
        call self%os_out%set(ind, 'fsc',     trim(fsc_abspath))
        call self%os_out%set(ind, 'imgkind', 'fsc')
        call self%os_out%set(ind, 'state',   real(state))
        call self%os_out%set(ind, 'box',   real(box))
    end subroutine add_fsc2os_out

    subroutine add_vol2os_out( self, vol, smpd, state, which_imgkind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: vol, which_imgkind
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: state
        character(len=:), allocatable :: vol_abspath, imgkind
        integer :: n_os_out, ind, i, ldim(3), ifoo
        select case(trim(which_imgkind))
            case('vol_cavg','vol','vol_filt','vol_msk')
                ! all good
            case DEFAULT
                write(*,*)'Invalid VOL kind: ', trim(which_imgkind), '; sp_project :: add_vol2os_out'
                stop 'sp_project :: add_vol2os_out'
        end select
        ! full path and existence check
        call simple_full_path(vol, vol_abspath, 'sp_project :: add_vol2os_out')
        ! find_dimension of inputted volume
        call find_ldim_nptcls(vol_abspath, ldim, ifoo)
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out)
        else
            select case(trim(which_imgkind))
                case('vol_msk','vol_cavg')
                    ! one volume type for all states
                    ind = 0
                    do i=1,n_os_out
                        if( self%os_out%isthere(i,'imgkind') )then
                            call self%os_out%getter(i,'imgkind',imgkind)
                            if(trim(imgkind).eq.trim(which_imgkind))then
                                ind = i
                                exit
                            endif
                        endif
                    end do
                    if( ind == 0 )then
                        n_os_out = n_os_out + 1
                        call self%os_out%reallocate(n_os_out)
                        ind = n_os_out
                    endif
                case DEFAULT
                    ! one volume per state
                    ind = 0
                    do i=1,n_os_out
                        if( self%os_out%isthere(i,'imgkind') )then
                            call self%os_out%getter(i,'imgkind',imgkind)
                            if(trim(imgkind).eq.trim(which_imgkind))then
                                if( self%os_out%get_state(i) == state )then
                                    ind = i
                                    exit
                                endif
                            endif
                        endif
                    end do
                    if( ind == 0 )then
                        n_os_out = n_os_out + 1
                        call self%os_out%reallocate(n_os_out)
                        ind = n_os_out
                    endif
            end select
        endif
        ! fill-in field
        call self%os_out%set(ind, 'vol',     trim(vol_abspath))
        call self%os_out%set(ind, 'box',     real(ldim(1)))
        call self%os_out%set(ind, 'smpd',    smpd)
        call self%os_out%set(ind, 'imgkind', trim(which_imgkind))
        call self%os_out%set(ind, 'state',   real(state))
    end subroutine add_vol2os_out

    subroutine add_entry2os_out( self, which_imgkind, ind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        integer,           intent(out)   :: ind
        character(len=:), allocatable :: imgkind
        integer :: n_os_out, i
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out)
        else
            ind = 0
            do i=1,n_os_out
                if( self%os_out%isthere(i,'imgkind') )then
                    call self%os_out%getter(i,'imgkind',imgkind)
                    if(trim(imgkind).eq.trim(which_imgkind))then
                        ind = i
                        exit
                    endif
                endif
            end do
            if( ind == 0 )then
                n_os_out = n_os_out + 1
                call self%os_out%reallocate(n_os_out)
                ind = n_os_out
            endif
        endif
    end subroutine add_entry2os_out

    subroutine get_cavgs_stk( self, stkname, ncls, smpd, fail )
        class(sp_project),             intent(inout) :: self
        character(len=:), allocatable, intent(inout) :: stkname
        integer,                       intent(out)   :: ncls
        real,                          intent(out)   :: smpd
        logical,          optional,    intent(in)    :: fail
        character(len=:), allocatable :: imgkind
        integer :: n_os_out, ind, i, cnt
        logical :: fail_here
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) stop 'ERROR! trying to fetch from empty os_out field; sp_project :: get_cavgs_stk'
        ! look for cavgs
        ind = 0
        cnt = 0
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind)
                if(trim(imgkind).eq.'cavg')then
                    ind = i
                    cnt = cnt + 1
                endif
            endif
        end do
        if( fail_here )then
            if( cnt > 1 )  stop 'ERROR! multiple os_out entries with imgkind=cavg, aborting...; sp_project :: get_cavgs_stk'
            if( cnt == 0 ) stop 'ERROR! no os_out entry with imgkind=cavg identified, aborting...; sp_project :: get_cavgs_stk'
            ! set return values
            if( allocated(stkname) ) deallocate(stkname)
            call self%os_out%getter(ind,'stk',stkname)
            ncls = nint(self%os_out%get(ind, 'nptcls'))
            smpd = self%os_out%get(ind, 'smpd')
        else
            stkname = NIL
            ncls    = 0
            smpd    = 0.
        endif
    end subroutine get_cavgs_stk

    subroutine get_vol( self, imgkind, state, vol_fname, smpd, box)
        class(sp_project),             intent(inout) :: self
        character(len=*),              intent(in)    :: imgkind
        integer,                       intent(in)    :: state
        character(len=:), allocatable, intent(inout) :: vol_fname
        real,                          intent(out)   :: smpd
        integer,                       intent(out)   :: box
        character(len=:), allocatable :: imgkind_here
        integer :: i, ind, cnt
        select case(trim(imgkind))
            case('vol_cavg','vol','vol_filt','vol_msk')
                ! all good
            case DEFAULT
                write(*,*)'Invalid VOL kind: ', trim(imgkind), '; sp_project :: get_vol'
                stop 'sp_project :: get_vol'
        end select
        ! defaults
        if( allocated(vol_fname) ) deallocate(vol_fname)
        allocate(vol_fname, source='')
        box  = 0
        smpd = 0.
        ! fetch index
        ind = 0
        cnt = 0
        select case(trim(imgkind))
            case('vol_cavg', 'vol_msk')
                do i=1,self%os_out%get_noris()
                    if( self%os_out%isthere(i,'imgkind') )then
                        call self%os_out%getter(i,'imgkind',imgkind_here)
                        if(trim(imgkind).eq.trim(imgkind_here))then
                            ind = i
                            cnt = cnt + 1
                        endif
                    endif
                enddo
            case DEFAULT
                do i=1,self%os_out%get_noris()
                    if( self%os_out%isthere(i,'imgkind') )then
                        call self%os_out%getter(i,'imgkind',imgkind_here)
                        if(trim(imgkind).eq.trim(imgkind_here))then
                            if( self%os_out%get_state(i).eq.state )then
                                ind = i
                                cnt = cnt + 1
                            endif
                        endif
                    endif
                enddo
        end select
        if( cnt == 0 )then
            if( trim(imgkind).eq.'vol_msk')then
                ! we do not fall over if the volume mask is absent
                return
            else
                stop 'ERROR! no os_out entry with imgkind=volXXX identified, aborting...; sp_project :: get_vol'
            endif
        endif
        if( cnt > 1 )  stop 'ERROR! multiple os_out entries with imgkind=volXXX, aborting...; sp_project :: get_vol'
        ! set output
        deallocate(vol_fname)
        call self%os_out%getter(ind, 'vol', vol_fname)
        smpd = self%os_out%get(ind,  'smpd')
        box  = self%os_out%get(ind,  'box')
    end subroutine get_vol

    subroutine get_fsc( self, state, fsc_fname, box )
        class(sp_project),             intent(inout) :: self
        integer,                       intent(in)    :: state
        character(len=:), allocatable, intent(inout) :: fsc_fname
        integer,                       intent(out)   :: box
        character(len=:), allocatable :: imgkind_here
        integer :: i, ind, cnt
        ! defaults
        if( allocated(fsc_fname) ) deallocate(fsc_fname)
        allocate(fsc_fname, source=NIL)
        ! fetch index
        ind   = 0
        cnt   = 0
        do i=1,self%os_out%get_noris()
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind_here)
                if(trim(imgkind_here).eq.'fsc')then
                    if( self%os_out%get_state(i).eq.state )then
                        ind = i
                        cnt = cnt + 1
                    endif
                endif
            endif
        enddo
        if( cnt == 0 )stop 'ERROR! no os_out entry with imgkind=fsc identified, aborting...; sp_project :: get_fsc'
        if( cnt > 1 ) stop 'ERROR! multiple os_out entries with imgkind=fsc, aborting...; sp_project :: get_fsc'
        ! set output
        deallocate(fsc_fname)
        call self%os_out%getter(ind, 'fsc', fsc_fname)
        box = nint(self%os_out%get(ind, 'box'))
    end subroutine get_fsc

    subroutine get_frcs( self, frcs, which_imgkind, fail )
        class(sp_project),             intent(inout) :: self
        character(len=:), allocatable, intent(inout) :: frcs
        character(len=*),              intent(in)    :: which_imgkind
        logical,          optional,    intent(in)    :: fail
        character(len=:), allocatable :: imgkind
        integer :: n_os_out, ind, i, cnt
        logical :: fail_here
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                write(*,*)'Invalid FRC kind: ', trim(which_imgkind)
                stop 'sp_project :: get_frcs'
        end select
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) stop 'ERROR! trying to fetch from empty os_out field; sp_project :: get_frcs'
        ! look for cavgs
        ind = 0
        cnt = 0
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind)
                if(trim(imgkind).eq.trim(which_imgkind))then
                    ind = i
                    cnt = cnt + 1
                endif
            endif
        end do
        if( allocated(frcs) ) deallocate(frcs)
        if( fail_here )then
            if( cnt > 1 )  stop 'ERROR! multiple os_out entries with imgkind=frcXD, aborting...; sp_project :: get_frcs'
            if( cnt == 0 ) stop 'ERROR! no os_out entry with imgkind=frcsXD identified, aborting...; sp_project :: get_frcs'
            ! set return values
            call self%os_out%getter(ind,'frcs',frcs)
        else
            frcs = NIL
        endif
    end subroutine get_frcs

    integer function get_nptcls_from_osout( self )
        class(sp_project), target, intent(inout) :: self
        character(len=LONGSTRLEN) :: imgkind
        integer :: i, nos
        get_nptcls_from_osout = 0
        nos = self%os_out%get_noris()
        if( nos == 0 )return
        ! look for cavgs, defaults to zero for other entries (frcs, volumes)
        do i=1,nos
            if( self%os_out%isthere(i,'imgkind') )then
                imgkind = self%os_out%get_static(i,'imgkind')
                if(trim(imgkind).eq.'cavg')then
                    if( self%os_out%isthere(i,'fromp').and.self%os_out%isthere(i,'top') )then
                        get_nptcls_from_osout = get_nptcls_from_osout +&
                            &nint(self%os_out%get(i,'top')) - nint(self%os_out%get(i,'fromp')) + 1
                    else
                        write(*,*) 'Missing fromp and top entries in cavg ', i
                        stop 'Missing fromp and top entries in cavg ; sp_project :: get_nptcls_from_osout'
                    endif
                endif
            endif
        end do
    end function get_nptcls_from_osout

    ! getters

    integer function get_nptcls( self )
        class(sp_project), target, intent(inout) :: self
        integer :: i, nos
        get_nptcls = 0
        nos        = self%os_stk%get_noris()
        if( nos == 0 )return
        do i=1,nos
            get_nptcls = get_nptcls + nint(self%os_stk%get(i,'nptcls'))
        enddo
        ! sanity check
        if( self%os_stk%isthere(nos,'top') )then
            if( nint(self%os_stk%get(nos,'top')) /=  get_nptcls )then
                write(*,*) 'nptcls from ptcls', get_nptcls
                write(*,*) 'nptcls from top', nint(self%os_stk%get(nos,'top'))
                stop 'ERROR! total # particles .ne. last top index; sp_project :: get_nptcls'
            endif
        endif
    end function get_nptcls

    integer function get_box( self )
        class(sp_project), target, intent(inout) :: self
        integer :: n_os_stk
        get_box  = 0
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            stop 'ERROR! empty os_stk field! sp_project :: get_box'
        endif
        get_box = nint( self%os_stk%get(1,'box') )
    end function get_box

    real function get_smpd( self )
        class(sp_project), target, intent(inout) :: self
        integer :: n_os_stk
        get_smpd  = 0.
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            stop 'ERROR! empty os_stk field! sp_project :: get_smpd'
        endif
        get_smpd = self%os_stk%get(1,'smpd')
    end function get_smpd

    integer function get_nmics( self )
        class(sp_project), target, intent(inout) :: self
        character(len=:), allocatable :: imgkind
        integer :: i
        get_nmics = 0
        do i=1,self%os_mic%get_noris()
            call self%os_mic%getter(i,'imgkind',imgkind)
            if( trim(imgkind).eq.'mic' ) get_nmics = get_nmics + 1
        enddo
    end function get_nmics

    integer function get_nmovies( self )
        class(sp_project), target, intent(inout) :: self
        character(len=:), allocatable :: imgkind
        integer :: i
        get_nmovies = 0
        do i=1,self%os_mic%get_noris()
            call self%os_mic%getter(i,'imgkind',imgkind)
            if( trim(imgkind).eq.'movie' ) get_nmovies = get_nmovies + 1
        enddo
    end function get_nmovies

    integer function get_nintgs( self )
        class(sp_project), target, intent(inout) :: self
        integer :: i
        get_nintgs = 0
        do i=1,self%os_mic%get_noris()
            if( self%os_mic%isthere(i,'intg') )get_nintgs = get_nintgs + 1
        enddo
    end function get_nintgs

    character(len=STDLEN) function get_ctfflag( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer          :: ptcl_field
        character(len=:), allocatable :: ctfflag
        integer              :: stkind, ind_in_stk
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case('cls2D', 'cls3D')
                get_ctfflag = 'no'
                return
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype), ' is not supported by this method'
                stop 'sp_project :: get_ctfflag'
        end select
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, 1, stkind, ind_in_stk)
        ! CTF flag
        if( self%os_stk%isthere(stkind, 'ctf') )then
            call self%os_stk%getter(stkind, 'ctf', ctfflag)
        else if( ptcl_field%isthere(1, 'ctf') )then
            call ptcl_field%getter(1, 'ctf', ctfflag)
        else
            ctfflag = 'no'
        endif
        get_ctfflag = trim(ctfflag)
    end function get_ctfflag

    integer function get_ctfflag_type( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        character(len=:), allocatable :: ctfflag
        ctfflag = self%get_ctfflag(oritype)
        select case(trim(ctfflag))
            case('no')
                get_ctfflag_type = CTFFLAG_NO
            case('yes')
                get_ctfflag_type = CTFFLAG_YES
            case('mul')
                stop 'ERROR ctf=mul deprecated; simple_sp_project :: get_ctfflag_type'
            case('flip')
                get_ctfflag_type = CTFFLAG_FLIP
            case DEFAULT
                print *, 'ctf flag:', trim(ctfflag)
                stop 'Unsupported ctf flag; simple_sp_project :: get_ctfflag_type'
        end select
    end function get_ctfflag_type

    logical function has_phaseplate( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer          :: ptcl_field
        character(len=:), allocatable :: phaseplate
        integer              :: stkind, ind_in_stk
        nullify(ptcl_field)
        select case(trim(oritype))
            case('ptcl2D', 'ptcl3D')
                ! all good
            case('cls3D')
                has_phaseplate = .false.
                return
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype), ' is not supported by this method'
                stop 'sp_project :: has_phaseplate'
        end select
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, 1, stkind, ind_in_stk)
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
        end select
        ! get info
        if( self%os_stk%isthere(stkind, 'phaseplate') )then
            call self%os_stk%getter(stkind, 'phaseplate', phaseplate)
        else if( ptcl_field%isthere(1, 'phaseplate') )then
            call ptcl_field%getter(1, 'phaseplate', phaseplate)
        else
            phaseplate = 'no'
        endif
        has_phaseplate = trim(phaseplate).eq.'yes'
    end function has_phaseplate

    function get_ctfparams( self, oritype, iptcl ) result( ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        class(oris), pointer          :: ptcl_field
        character(len=:), allocatable :: ctfflag
        type(ctfparams)      :: ctfvars
        integer              :: stkind, ind_in_stk
        logical              :: dfy_was_there, l_noctf
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype), ' is not supported by this method'
                stop 'sp_project :: get_ctfparams'
        end select
        ! extract the CTF parameters
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk)
        ! sampling distance
        if( self%os_stk%isthere(stkind, 'smpd') )then
            ctfvars%smpd = self%os_stk%get(stkind, 'smpd')
        else
            write(*,*) 'ERROR! smpd (sampling distance) lacking in os_stk and ptcl fields'
            stop 'sp_project :: get_ctfparams'
        endif
        ! CTF flag
        if( self%os_stk%isthere(stkind, 'ctf') )then
            ctfflag = trim(self%os_stk%get_static(stkind, 'ctf'))
        else
            ctfflag = NIL
        endif
        l_noctf = .false.
        select case(trim(ctfflag))
            case(NIL)
                write(*,*) 'ERROR! ctf key lacking in os_stk_field & ptcl fields'
                stop 'sp_project :: get_ctfparams'
            case('no')
                ctfvars%ctfflag = CTFFLAG_NO
                l_noctf = .true.
            case('yes')
                ctfvars%ctfflag = CTFFLAG_YES
            case('mul')
                stop 'ERROR ctf=mul deprecated; simple_sp_project :: get_ctfparams'
            case('flip')
                ctfvars%ctfflag = CTFFLAG_FLIP
            case DEFAULT
                write(*,*)'unsupported ctf flag:', trim(ctfflag), stkind, iptcl
                stop 'unsupported ctf flag; simple_sp_project :: get_ctfparams'
        end select
        ! acceleration voltage
        if( self%os_stk%isthere(stkind, 'kv') )then
            ctfvars%kv = self%os_stk%get(stkind, 'kv')
        else
            write(*,*) 'ERROR! kv (acceleration voltage) lacking in os_stk_field'
            stop 'sp_project :: get_ctfparams'
        endif
        ! spherical aberration constant
        if( self%os_stk%isthere(stkind, 'cs') )then
            ctfvars%cs = self%os_stk%get(stkind, 'cs')
        else
            write(*,*) 'ERROR! cs (spherical aberration constant) lacking in os_stk_field'
            stop 'sp_project :: get_ctfparams'
        endif
        ! fraction of amplitude contrast
        if( self%os_stk%isthere(stkind, 'fraca') )then
            ctfvars%fraca = self%os_stk%get(stkind, 'fraca')
        else
            write(*,*) 'ERROR! fraca (fraction of amplitude contrast) lacking in os_stk_field'
            stop 'sp_project :: get_ctfparams'
        endif
        if( l_noctf )then
            ctfvars%dfx     = 0.
            ctfvars%dfy     = 0.
            ctfvars%angast  = 0.
            ctfvars%phshift = 0.
            return
        endif
        ! defocus in x
        if( ptcl_field%isthere(iptcl, 'dfx') )then
            ctfvars%dfx = ptcl_field%get(iptcl, 'dfx')
        else
            write(*,*) 'ERROR! dfx (defocus in x) lacking in ptcl_field'
            call ptcl_field%print_(iptcl)
            stop 'sp_project :: get_ctfparams'
        endif
        ! defocus in y
        dfy_was_there = .false.
        if( ptcl_field%isthere(iptcl, 'dfy') )then
            ctfvars%dfy = ptcl_field%get(iptcl, 'dfy')
            dfy_was_there = .true.
        else
            ctfvars%dfy = ctfvars%dfx
        endif
        ! angle of astigmatism
        if( ptcl_field%isthere(iptcl, 'angast') )then
            ctfvars%angast = ptcl_field%get(iptcl, 'angast')
        else
            if( dfy_was_there )then
                write(*,*) 'ERROR! astigmatic CTF model requires angast (angle of astigmatism) lacking in os_stk field'
                stop 'sp_project :: get_ctfparams'
            else
                ctfvars%angast = 0.
            endif
        endif
        ! additional phase shift
        if( ptcl_field%isthere(iptcl, 'phshift') )then
            ctfvars%phshift = ptcl_field%get(iptcl, 'phshift')
        else
            ctfvars%phshift = 0.
        endif
    end function get_ctfparams

    logical function is_virgin_field( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer :: os
        integer :: i, n
        nullify(os)
        is_virgin_field = .false.
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                os => self%os_ptcl2D
            case('ptcl3D')
                os => self%os_ptcl3D
            case('cls3D')
                os => self%os_cls3D
            case DEFAULT
                write(*,*) 'oritype: ', trim(oritype), ' is not supported by this method'
                stop 'sp_project :: is_virgin_field'
        end select
        n = os%get_noris()
        if( n == 0 )then
            write(*,*) 'WARNING! cannot check virginity of non-existent field (touched for the very first time???)'
            print *,'WARNING! cannot check virginity of non-existent field (touched for the very first time???)'
            return
        endif
        do i=1,n
            if( os%has_been_searched(i) )return
        enddo
        is_virgin_field = .true.
    end function is_virgin_field

    ! modifiers

    subroutine set_sp_oris( self, which_imgkind, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        class(oris),       intent(inout) :: os
        select case(trim(which_imgkind))
            case('mic')
                self%os_mic    = os
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
            case('out')
                self%os_out    = os
            case('projinfo')
                self%projinfo  = os
            case('jobproc')
                self%jobproc   = os
            case('compenv')
                self%compenv   = os
            case DEFAULT
                stop 'unsupported which_imgkind flag; sp_project :: set_sp_oris'
        end select
    end subroutine set_sp_oris

    subroutine scale_projfile( self, smpd_target, new_projfile, cline, cline_scale, dir )
        ! this probably needs an oritype input for dealing with scale class averages
        use simple_cmdline, only: cmdline
        class(sp_project),             intent(inout) :: self
        real,                          intent(inout) :: smpd_target
        character(len=:), allocatable, intent(out)   :: new_projfile
        class(cmdline),                intent(inout) :: cline
        class(cmdline),                intent(out)   :: cline_scale
        character(len=*), optional,    intent(in)    :: dir
        character(len=:), allocatable :: projfile, projname, new_projname
        real    :: scale_factor, smpd_sc, msk_sc, smpd, msk
        integer :: box, box_sc, istk, n_os_stk
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            stop 'Empty stack object! simple_sp_project :: scale_projfile'
        endif
        call self%projinfo%getter(1, 'projfile', projfile)
        call self%projinfo%getter(1, 'projname', projname)
        if( trim(projname).eq.'' )then
            projname = get_fbody(projfile, METADATA_EXT, separator=.false.)
        endif
        ! dimensions
        smpd = self%get_smpd()
        box  = self%get_box()
        call autoscale(box, smpd, smpd_target, box_sc, smpd_sc, scale_factor)
        call cline_scale%set('prg',      'scale_project')
        call cline_scale%set('scale',    scale_factor)
        call cline_scale%set('projfile', projfile)
        call cline_scale%set('smpd',     smpd_sc)
        if(present(dir))call cline_scale%set('dir_target',trim(dir)//'/')
        if( box == box_sc )then
            ! no scaling
            new_projfile = trim(projfile)
            return
        endif
        ! parameter updates
        if( cline%defined('msk') )then
            msk = cline%get_rarg('msk')
            msk_sc = msk * scale_factor
            call cline%set('msk', msk_sc)
        endif
        do istk = 1,n_os_stk
            call self%os_stk%set(istk, 'smpd', real(smpd_sc))
            call self%os_stk%set(istk, 'box', real(box_sc))
        enddo
        call self%os_ptcl2D%mul_shifts(scale_factor)
        call self%os_ptcl3D%mul_shifts(scale_factor)
        ! name changes and list for scaling job
        new_projname = trim(projname)//SCALE_SUFFIX
        new_projfile = trim(new_projname)//'.simple'
        call cline%set('projname', trim(new_projname))
        call cline%delete('projfile')
        call self%update_projinfo( cline )
        if(present(dir))then
            call self%add_scale_tag(dir=trim(dir)//'/')
        else
            call self%add_scale_tag
        endif
        ! save
        call self%write()
        ! command line for scaling
        call cline_scale%set('newbox', real(box_sc))
        if( cline%defined('nthr') )   call cline_scale%set('nthr', cline%get_rarg('nthr'))
        if( cline%defined('nparts') ) call cline_scale%set('nparts', cline%get_rarg('nparts'))
    end subroutine scale_projfile

    !> for merging alignment documents from SIMPLE runs in distributed mode
    subroutine merge_algndocs( self, nptcls, ndocs, oritype, fbody, numlen_in )
        use simple_map_reduce, only: split_nobjs_even
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: nptcls, ndocs
        character(len=*),  intent(in)    :: oritype, fbody
        integer, optional, intent(in)    :: numlen_in
        class(oris),          pointer :: os_ptr
        integer,          allocatable :: parts(:,:)
        character(len=:), allocatable :: fname
        type(binoris) :: bos_doc
        type(oris)    :: os_part
        integer       :: i, iptcl, cnt, numlen, n_records, partsz, isegment
        numlen = len(int2str(ndocs))
        if( present(numlen_in) ) numlen = numlen_in
        parts  = split_nobjs_even(nptcls, ndocs)
        ! convert from flag to enumerator to integer
        isegment = oritype2segment(oritype)
        ! allocate merged oris
        call self%new_seg_with_ptr( nptcls, oritype, os_ptr )
        ! read & transfer
        do i=1,ndocs
            ! read part
            fname     = trim(adjustl(fbody))//int2str_pad(i,numlen)//'.simple'
            call bos_doc%open(trim(fname))
            n_records = bos_doc%get_n_records(isegment)
            partsz    = parts(i,2) - parts(i,1) + 1
            if( n_records /= partsz )then
                write(*,*) 'ERROR, # records does not match expectation'
                write(*,*) 'EXTRACTED FROM file: ', trim(fname)
                write(*,*) 'n_records: ', n_records
                write(*,*) 'CALCULATED FROM input p%nptcls/p%ndocs'
                write(*,*) 'fromto: ', parts(i,1), parts(i,2)
                write(*,*) 'partsz: ', partsz
                stop
            endif
            call os_part%new(n_records)
            call bos_doc%read_segment(isegment, os_part)
            call bos_doc%close()
            ! transfer to self
            cnt = 0
            do iptcl = parts(i,1), parts(i,2)
                cnt = cnt + 1
                call os_ptr%set_ori(iptcl, os_part%get_ori(cnt))
            enddo
        end do
        call self%write_segment_inside(oritype)
    end subroutine merge_algndocs

    ! this map2ptcls routine assumes that any selection of class averages is done
    ! exclusively by state=0 flagging without any physical deletion
    subroutine map2ptcls( self )
        class(sp_project), intent(inout) :: self
        integer, allocatable :: particles(:)
        type(ori) :: ori2d, ori_comp, o
        integer   :: ncls, icls, iptcl, pind, noris_ptcl3D, noris_ptcl2D
        real      :: corr, rproj, rstate
        if( self%is_virgin_field('cls3D') )then
            write(*,*) 'ERROR! os_cls3D is virgin field; nothing to map back'
            stop 'sp_project :: map2ptcls'
        endif
        if( self%is_virgin_field('ptcl2D') )then
            ! 2D was not done with SIMPLE but class averages imported from elsewhere
            return
        endif
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D ) call self%os_ptcl3D%new(noris_ptcl2D)
        ! do the mapping
        ncls = self%os_cls3D%get_noris()
        do icls=1,ncls
            ! get particle indices
            call self%os_ptcl2D%get_pinds(icls, 'class', particles)
            ! get 3d ori info
            o      = self%os_cls3D%get_ori(icls)
            rproj  = o%get('proj')
            rstate = o%get('state')
            corr   = o%get('corr')
            if( allocated(particles) )then
                do iptcl=1,size(particles)
                    ! get particle index
                    pind  = particles(iptcl)
                    ! get 2d ori
                    ori2d = self%os_ptcl2D%get_ori(pind)
                    ! transfer original parameters in self%os_ptcl2D
                    ori_comp = self%os_ptcl2D%get_ori(pind)
                    ! compose ori3d and ori2d
                    call o%compose3d2d(ori2d, ori_comp)
                    ! update state in self%os_ptcl2D
                    call self%os_ptcl2D%set(pind, 'state', rstate)
                    ! set 3D orientation in self%os_ptcl3D
                    call self%os_ptcl3D%set_ori(pind, ori_comp)
                    ! set proj/state/corr
                    call self%os_ptcl3D%set(pind, 'proj',  rproj)
                    call self%os_ptcl3D%set(pind, 'state', rstate)
                    call self%os_ptcl3D%set(pind, 'corr',  corr)
                end do
                deallocate(particles)
            endif
        end do
        ! state = 0 all entries that don't have a state label
        do iptcl=1,noris_ptcl2D
            if( .not. self%os_ptcl2D%isthere(iptcl, 'state') ) call self%os_ptcl2D%set(iptcl, 'state', 0.)
            if( .not. self%os_ptcl3D%isthere(iptcl, 'state') ) call self%os_ptcl3D%set(iptcl, 'state', 0.)
        end do
    end subroutine map2ptcls

    ! printers

    subroutine print_info( self )
        class(sp_project), intent(in) :: self
        integer :: n
        n = self%os_mic%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in micrographs          segment (1) :', n
        n = self%os_stk%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in per-micrograph stack segment (2) :', n
        n = self%os_ptcl2D%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in per-particle 2D      segment (3) :', n
        n = self%os_cls2D%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in per-cluster  2D      segment (4) :', n
        n = self%os_cls3D%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in per-cluster  3D      segment (5) :', n
        n = self%os_ptcl3D%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in per-particle 3D      segment (6) :', n
        n = self%os_out%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in out                  segment (7) :', n
        n = self%projinfo%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in project info         segment (11):', n
        n = self%jobproc%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in jobproc              segment (12):', n
        n = self%compenv%get_noris()
        if( n > 0 ) write(*,'(a,1x,i10)') '# entries in compenv              segment (13):', n
    end subroutine print_info

    ! readers

    subroutine read( self, fname )
        class(sp_project),          intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        character(len=:), allocatable :: projfile
        integer :: isegment
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: read'
            endif
            projfile = trim(fname)
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        if( .not. file_exists(trim(projfile)) )then
            write(*,*) 'fname: ', trim(projfile)
            stop 'inputted file does not exist; sp_project :: read'
        endif
        call self%bos%open(projfile)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment)
        end do
        call self%bos%close
    end subroutine read

    subroutine read_ctfparams_state_eo( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment
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
    end subroutine read_ctfparams_state_eo

    subroutine read_segment( self, oritype, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
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
                isegment = oritype2segment(oritype)
                call self%bos%open(fname)
                call self%segreader(isegment)
                call self%bos%close
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        call self%os_mic%read(fname)
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
                    case('out')
                        call self%os_out%read(fname)
                    case('projinfo')
                        call self%projinfo%read(fname)
                    case('jobproc')
                        call self%jobproc%read(fname)
                    case('compenv')
                        call self%compenv%read(fname)
                    case DEFAULT
                        stop 'unsupported oritype flag; sp_project :: read_segment'
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
            case(MIC_SEG)
                call self%os_mic%new(n)
                call self%bos%read_segment(isegment, self%os_mic)
            case(STK_SEG)
                call self%os_stk%new(n)
                call self%bos%read_segment(isegment, self%os_stk,    only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(PTCL2D_SEG)
                call self%os_ptcl2D%new(n)
                call self%bos%read_segment(isegment, self%os_ptcl2D, only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(CLS2D_SEG)
                call self%os_cls2D%new(n)
                call self%bos%read_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%os_cls3D%new(n)
                call self%bos%read_segment(isegment, self%os_cls3D)
            case(PTCL3D_SEG)
                call self%os_ptcl3D%new(n)
                call self%bos%read_segment(isegment, self%os_ptcl3D, only_ctfparams_state_eo=only_ctfparams_state_eo)
            case(OUT_SEG)
                call self%os_out%new(n)
                call self%bos%read_segment(isegment, self%os_out)
            case(PROJINFO_SEG)
                call self%projinfo%new(n)
                call self%bos%read_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%jobproc%new(n)
                call self%bos%read_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%compenv%new(n)
                call self%bos%read_segment(isegment, self%compenv)
        end select
    end subroutine segreader

    ! writers

    subroutine write( self, fname, fromto, isegment )
        class(sp_project), intent(inout) :: self
        character(len=*), optional, intent(in) :: fname
        integer,          optional, intent(in) :: fromto(2)
        integer,          optional, intent(in) :: isegment
        character(len=:), allocatable :: projfile
        integer :: iseg
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
        if( present(isegment) )then
            call self%segwriter(isegment, fromto)
        else
            do iseg=1,MAXN_OS_SEG
                call self%segwriter(iseg, fromto)
            end do
        endif
        ! update header
        call self%bos%write_header
        call self%bos%close
    end subroutine write

    subroutine write_segment_inside( self, oritype, fname, fromto )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: oritype
        character(len=*), optional, intent(in)    :: fname
        integer,          optional, intent(in)    :: fromto(2)
        character(len=:), allocatable :: projfile
        integer :: iseg
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: write'
            endif
            projfile = trim(fname)
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        if( file_exists(projfile) )then
            iseg = oritype2segment(oritype)
            call self%bos%open(projfile, del_if_exists=.false.)
            call self%segwriter_inside(iseg, fromto)
        else
            call self%write(fname, fromto)
        endif
        ! no need to update header (taken care of in binoris object)
        call self%bos%close
    end subroutine write_segment_inside

    subroutine write_segment2txt( self, oritype, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        select case(fname2format(fname))
            case('O')
                stop 'write_segment2txt is not supported for *.simple project files; sp_project :: write_segment2txt'
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        if( self%os_mic%get_noris() > 0 )then
                            call self%os_mic%write(fname)
                        else
                            write(*,*) 'WARNING, no mic-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('stk')
                        if( self%os_stk%get_noris() > 0 )then
                            call self%os_stk%write(fname)
                        else
                            write(*,*) 'WARNING, no stk-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('ptcl2D')
                        if( self%os_ptcl2D%get_noris() > 0 )then
                            call self%os_ptcl2D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl2D-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('cls2D')
                        if( self%os_cls2D%get_noris() > 0 )then
                            call self%os_cls2D%write(fname)
                        else
                            write(*,*) 'WARNING, no cls2D-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('cls3D')
                        if( self%os_cls3D%get_noris() > 0 )then
                            call self%os_cls3D%write(fname,  fromto)
                        else
                            write(*,*) 'WARNING, no cls3D-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('ptcl3D')
                        if( self%os_ptcl3D%get_noris() > 0 )then
                            call self%os_ptcl3D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl3D-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('out')
                        if( self%os_out%get_noris() > 0 )then
                            call self%os_out%write(fname)
                        else
                            write(*,*) 'WARNING, no out-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('projinfo')
                        if( self%projinfo%get_noris() > 0 )then
                            call self%projinfo%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no projinfo-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('jobproc')
                        if( self%jobproc%get_noris() > 0 )then
                            call self%jobproc%write(fname)
                        else
                            write(*,*) 'WARNING, no jobproc-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case('compenv')
                        if( self%compenv%get_noris() > 0 )then
                            call self%compenv%write(fname)
                        else
                            write(*,*) 'WARNING, no compenv-type oris available to write; sp_project :: write_segment2txt'
                        endif
                    case DEFAULT
                        stop 'unsupported oritype flag; sp_project :: write_segment2txt'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: write_segment2txt'
        end select
    end subroutine write_segment2txt

    subroutine print_segment( self, oritype, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        integer, optional, intent(in)    :: fromto(2)
        integer :: ffromto(2), iori, noris
        logical :: fromto_present
        fromto_present = present(fromto)
        if( fromto_present ) ffromto = fromto
        select case(trim(oritype))
            case('mic')
                noris = self%os_mic%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_mic%ori2str(iori)
                    end do
                else
                    write(*,*) 'No mic-type oris available to print; sp_project :: print_segment'
                endif
            case('stk')
                noris = self%os_stk%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_stk%ori2str(iori)
                    end do
                else
                    write(*,*) 'No stk-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl2D')
                noris = self%os_ptcl2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_ptcl2D%ori2str(iori)
                    end do
                else
                    write(*,*) 'No ptcl2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls2D')
                noris = self%os_cls2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_cls2D%ori2str(iori)
                    end do
                else
                    write(*,*) 'No cls2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls3D')
                noris = self%os_cls3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_cls3D%ori2str(iori)
                    end do
                else
                    write(*,*) 'No cls3D-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl3D')
                noris = self%os_ptcl3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_ptcl3D%ori2str(iori)
                    end do
                else
                    write(*,*) 'No ptcl3D-type oris available to print; sp_project :: print_segment'
                endif
            case('out')
                noris = self%os_out%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%os_out%ori2str(iori)
                    end do
                else
                    write(*,*) 'No out-type oris available to print; sp_project :: print_segment'
                endif
            case('projinfo')
                noris = self%projinfo%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%projinfo%ori2str(iori)
                    end do
                else
                    write(*,*) 'No projinfo-type oris available to print; sp_project :: print_segment'
                endif
            case('jobproc')
                noris = self%jobproc%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%jobproc%ori2str(iori)
                    end do
                else
                    write(*,*) 'No jobproc-type oris available to print; sp_project :: print_segment'
                endif
            case('compenv')
                noris = self%compenv%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(*,'(a)') self%compenv%ori2str(iori)
                    end do
                else
                    write(*,*) 'No compenv-type oris available to print; sp_project :: print_segment'
                endif
            case DEFAULT
                stop 'unsupported oritype flag; sp_project :: print_segment'
        end select
    end subroutine print_segment

    subroutine segwriter( self, isegment, fromto )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: isegment
        integer, optional, intent(in)    :: fromto(2)
        select case(isegment)
            case(MIC_SEG)
                call self%bos%write_segment(isegment, self%os_mic, fromto)
            case(STK_SEG)
                call self%bos%write_segment(isegment, self%os_stk)
            case(PTCL2D_SEG)
                call self%bos%write_segment(isegment, self%os_ptcl2D, fromto)
            case(CLS2D_SEG)
                call self%bos%write_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%bos%write_segment(isegment, self%os_cls3D, fromto)
            case(PTCL3D_SEG)
                call self%bos%write_segment(isegment, self%os_ptcl3D, fromto)
            case(OUT_SEG)
                call self%bos%write_segment(isegment, self%os_out)
            case(PROJINFO_SEG)
                call self%bos%write_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%bos%write_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%bos%write_segment(isegment, self%compenv)
        end select
    end subroutine segwriter

    subroutine segwriter_inside( self, isegment, fromto )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: isegment
        integer, optional, intent(in)    :: fromto(2)
        select case(isegment)
            case(MIC_SEG)
                call self%bos%write_segment_inside(isegment, self%os_mic, fromto)
            case(STK_SEG)
                call self%bos%write_segment_inside(isegment, self%os_stk)
            case(PTCL2D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_ptcl2D, fromto)
            case(CLS2D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_cls3D, fromto)
            case(PTCL3D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_ptcl3D, fromto)
            case(OUT_SEG)
                call self%bos%write_segment_inside(isegment, self%os_out)
            case(PROJINFO_SEG)
                call self%bos%write_segment_inside(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%bos%write_segment_inside(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%bos%write_segment_inside(isegment, self%compenv)
        end select
    end subroutine segwriter_inside

    ! destructor

    subroutine kill( self )
        class(sp_project), intent(inout) :: self
        call self%os_stk%kill
        call self%os_ptcl2D%kill
        call self%os_cls2D%kill
        call self%os_cls3D%kill
        call self%os_ptcl3D%kill
        call self%os_out%kill
        call self%projinfo%kill
        call self%jobproc%kill
        call self%compenv%kill
    end subroutine kill

    ! private supporting subroutines / functions

    integer function oritype2segment( oritype )
        character(len=*),  intent(in) :: oritype
        select case(trim(oritype))
            case('mic')
                oritype2segment = MIC_SEG
            case('stk')
                oritype2segment = STK_SEG
            case('ptcl2D')
                oritype2segment = PTCL2D_SEG
            case('cls2D')
                oritype2segment = CLS2D_SEG
            case('cls3D')
                oritype2segment = CLS3D_SEG
            case('ptcl3D')
                oritype2segment = PTCL3D_SEG
            case('out')
                oritype2segment = OUT_SEG
            case('projinfo')
                oritype2segment = PROJINFO_SEG
            case('jobproc')
                oritype2segment = JOBPROC_SEG
            case('compenv')
                oritype2segment = COMPENV_SEG
            case DEFAULT
                stop 'unsupported oritype flag; sp_project :: oritype_flag2isgement'
        end select
    end function oritype2segment

end module simple_sp_project
